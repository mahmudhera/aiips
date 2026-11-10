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

fn preflight_output(
    input: &Path,
    source: &Metadata,
    output: &Path,
    force: bool,
) -> Result<Option<Metadata>, String> {
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
            Err(e) if e.kind() == io::ErrorKind::NotFound && metadata.file_type().is_symlink() => {
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
    if !fs::metadata(parent).map_err(|e| error(output, e))?.is_dir() {
        return Err(error(output, "destination parent must be a directory"));
    }
    Ok(existing)
}

impl OutputTransaction {
    fn install(
        input: &Path,
        source: &Metadata,
        output: &Path,
        force: bool,
    ) -> Result<Self, String> {
        let existing = preflight_output(input, source, output, force)?;
        if !force {
            // The placeholder is empty, so create it exclusively without requiring
            // hard-link support. create_new is atomic even after the preflight.
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(output)
                .map_err(|e| {
                    if e.kind() == io::ErrorKind::AlreadyExists {
                        error(
                            output,
                            "output already exists; use -f or --force to overwrite",
                        )
                    } else {
                        error(output, e)
                    }
                })?;
            return Ok(Self {
                destination: output.into(),
                backup: None,
                installed: true,
            });
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
        if existing.is_some() {
            let backup = StagedOutput::create(parent, output)?;
            fs::rename(output, backup.path()).map_err(|e| error(output, e))?;
            transaction.backup = Some(backup);
        }
        let result = fs::rename(staged.path(), output);
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

    fn finish(self) -> Option<String> {
        self.finish_with_cleanup(|path| fs::remove_file(path))
    }

    fn finish_with_cleanup(
        mut self,
        remove: impl FnOnce(&Path) -> io::Result<()>,
    ) -> Option<String> {
        // Successful stdout flush commits the new output. Cleanup errors must not
        // undo an already-reported success or let Drop remove the retained backup.
        self.installed = false;
        let mut backup = self.backup.take()?;
        let path = backup.0.take().unwrap();
        remove(&path).err().map(|e| format!(
            "{}: could not remove backup: {e}; new output retained at {}; original backup retained for cleanup",
            path.display(), self.destination.display()
        ))
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
            let (file, metadata) = open_input(&input)?;
            preflight_output(&input, &metadata, &destination, force)?;
            let summary = summarize(file, k).map_err(|e| error(&input, e))?;
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
            if let Some(warning) = transaction.finish() {
                eprintln!("aiips: warning: {warning}");
            }
            Ok(())
        }
        Command::Compare { inputs, k } => {
            let first = compare_input(&inputs[0], k)?;
            let second = compare_input(&inputs[1], k)?;
            report::statistics(
                output,
                &[(&inputs[0], first.as_ref()), (&inputs[1], second.as_ref())],
                k,
            )
            .map_err(|e| format!("stdout: {e}"))?;
            report::comparison(output, &inputs)
                .and_then(|()| output.flush())
                .map_err(|e| format!("stdout: {e}"))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::TestDirectory;

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
        let parent = TestDirectory::new();
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
        let parent = TestDirectory::new();
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
    #[cfg(unix)]
    #[test]
    fn backup_cleanup_failure_keeps_committed_output_and_warns() {
        let directory = TestDirectory::new();
        let input = directory.join("in.fa");
        let destination = directory.join("out.aiips");
        fs::write(&input, b">a\nACGT\n").unwrap();
        fs::write(&destination, b"original result").unwrap();
        let transaction =
            OutputTransaction::install(&input, &fs::metadata(&input).unwrap(), &destination, true)
                .unwrap();
        let backup = transaction.backup.as_ref().unwrap().path().to_path_buf();
        let warning = transaction
            .finish_with_cleanup(|path| {
                assert_eq!(path, backup);
                Err(io::ErrorKind::PermissionDenied.into())
            })
            .unwrap();
        assert!(warning.contains("could not remove backup"));
        assert!(warning.contains(&backup.display().to_string()));
        assert!(warning.contains(&destination.display().to_string()));
        assert_eq!(fs::read(destination).unwrap(), b"");
        assert_eq!(fs::read(backup).unwrap(), b"original result");
    }

    #[test]
    fn new_output_does_not_clobber_a_file_created_after_preflight() {
        let directory = TestDirectory::new();
        let input = directory.join("in.fa");
        let destination = directory.join("out.aiips");
        fs::write(&input, b">a\nACGT\n").unwrap();
        let metadata = fs::metadata(&input).unwrap();
        preflight_output(&input, &metadata, &destination, false).unwrap();
        fs::write(&destination, b"competing result").unwrap();
        let result = OutputTransaction::install(&input, &metadata, &destination, false);
        assert!(result.err().unwrap().contains("output already exists"));
        assert_eq!(fs::read(destination).unwrap(), b"competing result");
    }

    #[test]
    fn compare_statistics_and_flush_failures_identify_stdout() {
        let directory = TestDirectory::new();
        let input = directory.join("in.fa");
        fs::write(&input, b">a\nACGT\n").unwrap();
        for flush_only in [false, true] {
            let result = execute(
                Command::Compare {
                    inputs: [input.clone(), input.clone()],
                    k: NonZeroUsize::new(3).unwrap(),
                },
                &mut FailingWriter { flush_only },
            );
            assert!(result.unwrap_err().starts_with("stdout: "));
        }
    }

    #[test]
    fn compare_comparison_table_failure_identifies_stdout() {
        struct LimitedWriter(usize);
        impl Write for LimitedWriter {
            fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
                if self.0 == 0 {
                    return Err(io::Error::other("injected comparison write failure"));
                }
                let count = bytes.len().min(self.0);
                self.0 -= count;
                Ok(count)
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let directory = TestDirectory::new();
        let input = directory.join("in.fa");
        fs::write(&input, b">a\nACGT\n").unwrap();
        let k = NonZeroUsize::new(3).unwrap();
        let summary = summarize(File::open(&input).unwrap(), k).unwrap();
        let mut statistics = Vec::new();
        report::statistics(
            &mut statistics,
            &[(&input, Some(&summary)), (&input, Some(&summary))],
            k,
        )
        .unwrap();
        let result = execute(
            Command::Compare {
                inputs: [input.clone(), input],
                k,
            },
            &mut LimitedWriter(statistics.len()),
        );
        assert_eq!(
            result.unwrap_err(),
            "stdout: injected comparison write failure"
        );
    }
}
