use std::fs::{self, File, Metadata, OpenOptions};
use std::io::{self, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use aiips::fastx::{summarize, Summary};

use crate::{
    cli::{is_sketch, Command},
    report,
};

fn error(path: &Path, error: impl std::fmt::Display) -> String {
    format!("{}: {error}", path.display())
}

fn open_input(path: &Path) -> Result<(File, Metadata), String> {
    let metadata = fs::metadata(path).map_err(|e| error(path, e))?;
    if !metadata.is_file() {
        return Err(error(path, "input must be a regular file"));
    }
    let file = File::open(path).map_err(|e| error(path, e))?;
    Ok((file, metadata))
}

fn sequence_summary(path: &Path, k: NonZeroUsize) -> Result<(Summary, Metadata), String> {
    let (file, metadata) = open_input(path)?;
    let summary = summarize(file, k).map_err(|e| error(path, e))?;
    Ok((summary, metadata))
}

fn compare_input(path: &Path, k: NonZeroUsize) -> Result<Option<Summary>, String> {
    let (file, metadata) = open_input(path)?;
    if !is_sketch(path) {
        return summarize(file, k).map(Some).map_err(|e| error(path, e));
    }
    if metadata.len() != 0 {
        return Err(error(
            path,
            "only empty .aiips placeholder sketches are supported",
        ));
    }
    Ok(None)
}

#[cfg(unix)]
fn same_file(a: &Metadata, b: &Metadata, _: &Path, _: &Path) -> bool {
    use std::os::unix::fs::MetadataExt;
    a.dev() == b.dev() && a.ino() == b.ino()
}

#[cfg(not(unix))]
fn same_file(_: &Metadata, _: &Metadata, a: &Path, b: &Path) -> bool {
    // No portable file-ID API exists in std; conservatively disallow forced
    // replacement on these platforms below rather than risk an undetected alias.
    fs::canonicalize(a)
        .ok()
        .zip(fs::canonicalize(b).ok())
        .is_some_and(|(a, b)| a == b)
}

// Owns only paths created by this invocation. Cleanup must never panic.
struct StagedOutput(Option<PathBuf>);

impl StagedOutput {
    fn create(parent: &Path, destination: &Path) -> Result<Self, String> {
        static NEXT_ID: AtomicU64 = AtomicU64::new(0);
        loop {
            let id = NEXT_ID.fetch_add(1, Ordering::Relaxed);
            let path = parent.join(format!(".aiips-{}-{id}.tmp", std::process::id()));
            match OpenOptions::new().write(true).create_new(true).open(&path) {
                Ok(_) => return Ok(Self(Some(path))),
                Err(e) if e.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(error(destination, e)),
            }
        }
    }

    fn path(&self) -> &Path {
        self.0.as_deref().unwrap()
    }
}

impl Drop for StagedOutput {
    fn drop(&mut self) {
        if let Some(path) = &self.0 {
            let _ = fs::remove_file(path);
        }
    }
}

struct OutputTransaction {
    destination: PathBuf,
    backup: Option<StagedOutput>,
    installed: bool,
}

impl OutputTransaction {
    fn install(
        input: &Path,
        source: &Metadata,
        output: &Path,
        force: bool,
    ) -> Result<Self, String> {
        let existing = match fs::symlink_metadata(output) {
            Ok(metadata) => Some(metadata),
            Err(e) if e.kind() == io::ErrorKind::NotFound => None,
            Err(e) => return Err(error(output, e)),
        };
        if let Some(metadata) = &existing {
            if metadata.is_dir() {
                return Err(error(output, "output is a directory"));
            }
            if !metadata.is_file() && !metadata.file_type().is_symlink() {
                return Err(error(output, "output must be a regular file or symlink"));
            }
            let target = match fs::metadata(output) {
                Ok(target) => Some(target),
                Err(e)
                    if e.kind() == io::ErrorKind::NotFound && metadata.file_type().is_symlink() =>
                {
                    None
                }
                Err(e) => return Err(error(output, e)),
            };
            if target
                .as_ref()
                .is_some_and(|target| same_file(source, target, input, output))
            {
                return Err(error(output, "output resolves to the input file"));
            }
            if !force {
                return Err(error(
                    output,
                    "output already exists; use -f or --force to overwrite",
                ));
            }
            #[cfg(not(unix))]
            return Err(error(
                output,
                "forced replacement requires Unix file-identity safeguards",
            ));
        }
        let parent = output
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let staged = StagedOutput::create(parent, output)?;
        if let Some(metadata) = existing.as_ref().filter(|m| m.is_file()) {
            fs::set_permissions(staged.path(), metadata.permissions())
                .map_err(|e| error(output, e))?;
        }
        let mut transaction = Self {
            destination: output.into(),
            backup: None,
            installed: false,
        };
        if force && existing.is_some() {
            let backup = StagedOutput::create(parent, output)?;
            fs::rename(output, backup.path()).map_err(|e| error(output, e))?;
            transaction.backup = Some(backup);
        }
        let result = if force {
            fs::rename(staged.path(), output)
        } else {
            // Atomic no-clobber publication: a competing creator cannot be overwritten.
            fs::hard_link(staged.path(), output)
        };
        if let Err(e) = result {
            let message = if e.kind() == io::ErrorKind::AlreadyExists {
                error(
                    output,
                    "output already exists; use -f or --force to overwrite",
                )
            } else {
                error(output, e)
            };
            return Err(transaction.with_rollback(message));
        }
        transaction.installed = true;
        Ok(transaction)
    }

    fn rollback(&mut self) -> Result<(), String> {
        if let Some(mut backup) = self.backup.take() {
            if let Err(e) = fs::rename(backup.path(), &self.destination) {
                // Retain the backup on disk so a rollback failure cannot destroy it.
                let path = backup.0.take().unwrap();
                self.installed = false;
                return Err(format!(
                    "{}: rollback failed: {e}; original preserved at {}",
                    self.destination.display(),
                    path.display()
                ));
            }
            backup.0 = None;
        } else if self.installed {
            fs::remove_file(&self.destination).map_err(|e| error(&self.destination, e))?;
        }
        self.installed = false;
        Ok(())
    }

    fn with_rollback(&mut self, message: String) -> String {
        match self.rollback() {
            Ok(()) => message,
            Err(e) => format!("{message}; {e}"),
        }
    }

    fn finish(mut self) -> Result<(), String> {
        if let Some(backup) = &mut self.backup {
            if let Err(e) = fs::remove_file(backup.path()) {
                let message = error(backup.path(), e);
                return Err(self.with_rollback(message));
            }
            backup.0 = None;
        }
        self.backup = None;
        self.installed = false;
        Ok(())
    }
}

impl Drop for OutputTransaction {
    fn drop(&mut self) {
        let _ = self.rollback();
    }
}

pub fn execute(command: Command, output: &mut impl Write) -> Result<(), String> {
    match command {
        Command::Help(help) => writeln!(output, "{help}")
            .and_then(|()| output.flush())
            .map_err(|e| format!("stdout: {e}")),
        Command::Sketch {
            input,
            output: destination,
            k,
            force,
        } => {
            let (summary, metadata) = sequence_summary(&input, k)?;
            let mut transaction =
                OutputTransaction::install(&input, &metadata, &destination, force)?;
            let report = (|| -> io::Result<()> {
                report::statistics(output, &[(&input, Some(&summary))], k)?;
                writeln!(output, "output\t{}", report::label(&destination))?;
                output.flush()
            })();
            if let Err(e) = report {
                return Err(transaction.with_rollback(format!("stdout: {e}")));
            }
            transaction.finish()
        }
        Command::Compare { inputs, k } => {
            let first = compare_input(&inputs[0], k)?;
            let second = compare_input(&inputs[1], k)?;
            report::statistics(
                output,
                &[(&inputs[0], first.as_ref()), (&inputs[1], second.as_ref())],
                k,
            )
            .map_err(|e| e.to_string())?;
            report::comparison(output, &inputs)
                .and_then(|()| output.flush())
                .map_err(|e| format!("stdout: {e}"))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct FailingWriter {
        flush_only: bool,
    }

    impl Write for FailingWriter {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            if self.flush_only {
                Ok(bytes.len())
            } else {
                Err(io::Error::other("injected write failure"))
            }
        }
        fn flush(&mut self) -> io::Result<()> {
            Err(io::Error::other("injected flush failure"))
        }
    }

    fn stdout_failure(flush_only: bool, force: bool) {
        let parent = std::env::temp_dir();
        let input = StagedOutput::create(&parent, &parent).unwrap();
        fs::write(input.path(), b">a\nACGT\n").unwrap();
        let destination = StagedOutput::create(&parent, &parent).unwrap();
        if force {
            fs::write(destination.path(), b"original result").unwrap();
        } else {
            fs::remove_file(destination.path()).unwrap();
        }
        let result = execute(
            Command::Sketch {
                input: input.path().into(),
                output: destination.path().into(),
                k: NonZeroUsize::new(3).unwrap(),
                force,
            },
            &mut FailingWriter { flush_only },
        );
        assert!(result.unwrap_err().contains(if flush_only {
            "injected flush failure"
        } else {
            "injected write failure"
        }));
        if force {
            assert_eq!(fs::read(destination.path()).unwrap(), b"original result");
        } else {
            assert!(!destination.path().exists());
        }
    }

    #[test]
    fn stdout_write_failure_removes_new_output() {
        stdout_failure(false, false);
    }
    #[cfg(unix)]
    #[test]
    fn stdout_write_failure_restores_replaced_output() {
        stdout_failure(false, true);
    }
    #[test]
    fn stdout_flush_failure_removes_new_output() {
        stdout_failure(true, false);
    }
    #[cfg(unix)]
    #[test]
    fn stdout_flush_failure_restores_replaced_output() {
        stdout_failure(true, true);
    }
    #[cfg(unix)]
    #[test]
    fn rollback_failure_keeps_original_backup_for_recovery() {
        let parent = std::env::temp_dir();
        let input = StagedOutput::create(&parent, &parent).unwrap();
        fs::write(input.path(), b">a\nACGT\n").unwrap();
        let destination = StagedOutput::create(&parent, &parent).unwrap();
        fs::write(destination.path(), b"original result").unwrap();
        let source = fs::metadata(input.path()).unwrap();
        let mut transaction =
            OutputTransaction::install(input.path(), &source, destination.path(), true).unwrap();
        let backup = transaction.backup.as_ref().unwrap().path().to_path_buf();
        fs::remove_file(destination.path()).unwrap();
        fs::create_dir(destination.path()).unwrap();
        let message = transaction.rollback().unwrap_err();
        assert!(message.contains("original preserved at"));
        assert!(message.contains(&backup.display().to_string()));
        drop(transaction);
        assert_eq!(fs::read(&backup).unwrap(), b"original result");
        fs::remove_dir(destination.path()).unwrap();
        fs::rename(backup, destination.path()).unwrap();
    }
}
