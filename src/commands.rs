use std::fs::{self, File, OpenOptions};
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

fn open_input(path: &Path) -> Result<File, String> {
    let metadata = fs::metadata(path).map_err(|e| error(path, e))?;
    if !metadata.is_file() {
        return Err(error(path, "input must be a regular file"));
    }
    File::open(path).map_err(|e| error(path, e))
}

fn sequence_summary(path: &Path, k: NonZeroUsize) -> Result<Summary, String> {
    summarize(open_input(path)?, k).map_err(|e| error(path, e))
}

fn compare_input(path: &Path, k: NonZeroUsize) -> Result<Option<Summary>, String> {
    if !is_sketch(path) {
        return sequence_summary(path, k).map(Some);
    }
    let file = open_input(path)?;
    if file.metadata().map_err(|e| error(path, e))?.len() != 0 {
        return Err(error(
            path,
            "only empty .aiips placeholder sketches are supported",
        ));
    }
    Ok(None)
}

// Removes only a temporary file created by this invocation, including on errors.
struct StagedOutput(PathBuf);

impl Drop for StagedOutput {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}

fn create_placeholder(input: &Path, output: &Path, force: bool) -> Result<(), String> {
    let source = fs::canonicalize(input).map_err(|e| error(input, e))?;
    if fs::canonicalize(output).is_ok_and(|destination| destination == source) {
        return Err(error(output, "output resolves to the input file"));
    }
    if !force {
        return OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output)
            .map(|_| ())
            .map_err(|e| {
                if e.kind() == io::ErrorKind::AlreadyExists {
                    error(
                        output,
                        "output already exists; use -f or --force to overwrite",
                    )
                } else {
                    error(output, e)
                }
            });
    }

    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    static NEXT_ID: AtomicU64 = AtomicU64::new(0);
    let staged = loop {
        let id = NEXT_ID.fetch_add(1, Ordering::Relaxed);
        let path = parent.join(format!(".aiips-{}-{id}.tmp", std::process::id()));
        match OpenOptions::new().write(true).create_new(true).open(&path) {
            Ok(file) => {
                drop(file);
                break StagedOutput(path);
            }
            Err(e) if e.kind() == io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(error(output, e)),
        }
    };
    fs::rename(&staged.0, output).map_err(|e| error(output, e))
}

pub fn execute(command: Command, output: &mut impl Write) -> Result<(), String> {
    match command {
        Command::Help(help) => writeln!(output, "{help}").map_err(|e| e.to_string()),
        Command::Sketch {
            input,
            output: destination,
            k,
            force,
        } => {
            let summary = sequence_summary(&input, k)?;
            create_placeholder(&input, &destination, force)?;
            report::statistics(output, &[(&input, Some(&summary))], k)
                .map_err(|e| e.to_string())?;
            writeln!(output, "output\t{}", report::label(&destination)).map_err(|e| e.to_string())
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
            report::comparison(output, &inputs).map_err(|e| e.to_string())
        }
    }
}
