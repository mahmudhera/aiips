use std::fs::{self, File, Metadata, OpenOptions};
use std::io::{self, Read, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use aiips::fastx::Summary;
use aiips::sketch::fracminhash::{self, FracMinHashSketch};
use aiips::sketch::minhash::{self, MinHashSketch};
use aiips::sketch::{self as sketch_format, Estimates, Sketch};

use crate::{
    cli::{is_sketch, AlgorithmOptions, Command, SketchParams},
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

/// Reads a sketch file after bounding its size by its own header (see
/// `sketch::max_file_bytes`): a FracMinHash file by its declared size, every
/// other file by the largest MinHash sketch.
fn read_sketch(path: &Path) -> Result<Sketch, String> {
    let (mut file, metadata) = open_input(path)?;
    let mut bytes = Vec::new();
    (&mut file)
        .take(sketch_format::HEADER_BYTES)
        .read_to_end(&mut bytes)
        .map_err(|e| error(path, e))?;
    let limit = sketch_format::max_file_bytes(&bytes);
    if metadata.len() > limit {
        return Err(error(path, "too large to be an aiips sketch"));
    }
    bytes.reserve(metadata.len().saturating_sub(bytes.len() as u64) as usize);
    file.take(limit.saturating_add(1).saturating_sub(bytes.len() as u64))
        .read_to_end(&mut bytes)
        .map_err(|e| error(path, e))?;
    if bytes.len() as u64 > limit {
        return Err(error(path, "too large to be an aiips sketch"));
    }
    sketch_format::parse(&bytes).map_err(|e| error(path, e))
}

/// A compare input: a sketch file, or a FASTA/FASTQ file sketched in memory.
enum CompareInput {
    Sketch(Sketch),
    Sequence,
}

/// `.aiips` names must be sketches. Other files are sketches when they start
/// with the sketch magic and are otherwise read as FASTA/FASTQ.
fn classify(path: &Path) -> Result<CompareInput, String> {
    if is_sketch(path) {
        return read_sketch(path).map(CompareInput::Sketch);
    }
    let (file, _) = open_input(path)?;
    let mut prefix = Vec::with_capacity(sketch_format::MAGIC.len());
    file.take(sketch_format::MAGIC.len() as u64)
        .read_to_end(&mut prefix)
        .map_err(|e| error(path, e))?;
    if prefix == sketch_format::MAGIC.as_bytes() {
        read_sketch(path).map(CompareInput::Sketch)
    } else {
        Ok(CompareInput::Sequence)
    }
}

fn unsupported_pair(paths: &[PathBuf; 2], kinds: [&str; 2]) -> String {
    format!(
        "comparison requires two sketches of one algorithm, minhash or fracminhash ({}: {}, {}: {})",
        paths[0].display(),
        kinds[0],
        paths[1].display(),
        kinds[1]
    )
}

/// Builds both comparison sketches. Sequence inputs take the algorithm and
/// parameters of a sketch input (given options must match it), or the options
/// and defaults when both inputs are sequences.
fn compare_sketches(
    paths: &[PathBuf; 2],
    options: &AlgorithmOptions,
) -> Result<[Sketch; 2], String> {
    let inputs = [classify(&paths[0])?, classify(&paths[1])?];
    let kind = |input: &CompareInput| match input {
        CompareInput::Sketch(sketch) => sketch.algorithm().name(),
        CompareInput::Sequence => "fasta/fastq",
    };
    if inputs.iter().any(
        |input| matches!(input, CompareInput::Sketch(sketch) if !sketch.algorithm().is_implemented()),
    ) {
        return Err(unsupported_pair(paths, [kind(&inputs[0]), kind(&inputs[1])]));
    }
    let template = paths
        .iter()
        .zip(&inputs)
        .find_map(|(path, input)| match input {
            CompareInput::Sketch(sketch) => SketchParams::of(sketch).map(|p| (path, p)),
            CompareInput::Sequence => None,
        });
    let params = match template {
        Some((path, template)) => options.matching(&template, path)?,
        None => options.params()?,
    };
    let mut sketches = Vec::with_capacity(2);
    for (path, input) in paths.iter().zip(inputs) {
        sketches.push(match input {
            CompareInput::Sketch(sketch) => sketch,
            CompareInput::Sequence => sketch_sequence(path, params.clone())?,
        });
    }
    Ok(sketches.try_into().unwrap())
}

fn sketch_sequence(path: &Path, params: SketchParams) -> Result<Sketch, String> {
    let (file, _) = open_input(path)?;
    sketch_file(path, file, params)
}

/// Sketches an already-opened FASTA/FASTQ file; refuses inputs without k-mers.
fn sketch_file(path: &Path, file: File, params: SketchParams) -> Result<Sketch, String> {
    let k = params.k();
    let sketch = match params {
        SketchParams::MinHash(p) => MinHashSketch::build(file, p).map(Sketch::MinHash),
        SketchParams::FracMinHash(p) => FracMinHashSketch::build(file, p).map(Sketch::FracMinHash),
    }
    .map_err(|e| error(path, e))?;
    if sketch.summary().is_some_and(|summary| summary.kmers == 0) {
        return Err(error(
            path,
            format!("no valid k-mers of length {k}; nothing to sketch"),
        ));
    }
    Ok(sketch)
}

/// Estimates for two sketches of one algorithm with identical parameters;
/// halts otherwise, naming each input's type or every differing field.
fn estimate_pair(paths: &[PathBuf; 2], sketches: &[Sketch; 2]) -> Result<Estimates, String> {
    let differ = |differences: Vec<String>| {
        format!(
            "{} and {}: sketch parameters differ: {}; re-sketch with identical parameters",
            paths[0].display(),
            paths[1].display(),
            differences.join(", ")
        )
    };
    match sketches {
        [Sketch::MinHash(a), Sketch::MinHash(b)] => {
            let differences = a.params.mismatches(&b.params);
            if !differences.is_empty() {
                return Err(differ(differences));
            }
            minhash::estimate(a, b)
        }
        [Sketch::FracMinHash(a), Sketch::FracMinHash(b)] => {
            let differences = a.params.mismatches(&b.params);
            if !differences.is_empty() {
                return Err(differ(differences));
            }
            if let Some(path) = paths
                .iter()
                .zip([a, b])
                .find_map(|(path, s)| s.hashes.is_empty().then_some(path))
            {
                return Err(error(path, fracminhash::EMPTY));
            }
            fracminhash::estimate(a, b)
        }
        [a, b] => Err(unsupported_pair(
            paths,
            [a.algorithm().name(), b.algorithm().name()],
        )),
    }
}

fn summary_row<'a>(path: &'a Path, sketch: &'a Sketch) -> (&'a Path, &'a Summary, NonZeroUsize) {
    // Only implemented sketches reach reporting.
    (path, sketch.summary().unwrap(), sketch.k().unwrap())
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

const MAX_UNIQUE_ATTEMPTS: u32 = 1000;

/// Retries `attempt` while it reports AlreadyExists, with a fixed bound so a
/// filesystem that rejects every new name cannot hang the command.
fn retry_unique<T>(mut attempt: impl FnMut() -> io::Result<T>) -> io::Result<T> {
    for _ in 0..MAX_UNIQUE_ATTEMPTS {
        match attempt() {
            Err(e) if e.kind() == io::ErrorKind::AlreadyExists => {}
            result => return result,
        }
    }
    Err(io::Error::new(
        io::ErrorKind::AlreadyExists,
        format!("no unused temporary name after {MAX_UNIQUE_ATTEMPTS} attempts"),
    ))
}

fn parent_dir(output: &Path) -> &Path {
    output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."))
}

fn exists_error(output: &Path) -> String {
    error(
        output,
        "output already exists; use -f or --force to overwrite",
    )
}

/// Diagnoses a failed create/rename at the destination. Suggests force only
/// when it was not already given and the occupant is not a directory.
fn install_error(output: &Path, e: io::Error, force: bool) -> String {
    if fs::symlink_metadata(output).is_ok_and(|m| m.is_dir()) {
        error(output, "output is a directory")
    } else if !force && e.kind() == io::ErrorKind::AlreadyExists {
        exists_error(output)
    } else {
        error(output, e)
    }
}

// Owns only paths created by this invocation. Cleanup must never panic.
struct StagedOutput(Option<PathBuf>);

impl StagedOutput {
    fn create(parent: &Path, destination: &Path) -> Result<Self, String> {
        static NEXT_ID: AtomicU64 = AtomicU64::new(0);
        let path = retry_unique(|| {
            let id = NEXT_ID.fetch_add(1, Ordering::Relaxed);
            let path = parent.join(format!(".aiips-{}-{id}.tmp", std::process::id()));
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&path)
                .map(|_| path)
        })
        .map_err(|e| error(destination, e))?;
        Ok(Self(Some(path)))
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

/// Full preflight before parsing input: existing-output safeguards plus a
/// destination-parent check.
fn preflight_output(
    input: &Path,
    source: &Metadata,
    output: &Path,
    force: bool,
) -> Result<Option<Metadata>, String> {
    let existing = check_existing(input, source, output, force)?;
    if !fs::metadata(parent_dir(output))
        .map_err(|e| error(output, e))?
        .is_dir()
    {
        return Err(error(output, "destination parent must be a directory"));
    }
    Ok(existing)
}

/// Existence, type, and alias safeguards. Rechecked just before install
/// because the output may change while the input is being parsed.
fn check_existing(
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
            return Err(exists_error(output));
        }
        #[cfg(not(unix))]
        return Err(error(
            output,
            "forced replacement requires Unix file-identity safeguards",
        ));
    }
    Ok(existing)
}

impl OutputTransaction {
    fn install(
        input: &Path,
        source: &Metadata,
        output: &Path,
        force: bool,
        contents: &[u8],
    ) -> Result<Self, String> {
        // A missing parent is reported by the create/rename below, so only the
        // existing-output safeguards need rechecking here.
        let existing = check_existing(input, source, output, force)?;
        Self::install_checked(output, force, existing, contents)
    }

    /// Installs `contents` after safeguards passed; `existing` is the
    /// recheck's result.
    fn install_checked(
        output: &Path,
        force: bool,
        existing: Option<Metadata>,
        contents: &[u8],
    ) -> Result<Self, String> {
        if !force {
            // Exclusive creation needs no hard-link support and stays atomic
            // even after the preflight. A failed write removes the new file.
            let mut file = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(output)
                .map_err(|e| install_error(output, e, false))?;
            let mut transaction = Self {
                destination: output.into(),
                backup: None,
                installed: true,
            };
            if let Err(e) = file.write_all(contents).and_then(|()| file.flush()) {
                drop(file);
                return Err(transaction.with_rollback(error(output, e)));
            }
            return Ok(transaction);
        }
        let parent = parent_dir(output);
        let staged = StagedOutput::create(parent, output)?;
        fs::write(staged.path(), contents).map_err(|e| error(output, e))?;
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
            if let Err(e) = fs::rename(output, backup.path()) {
                // Nothing was moved; remove both staging files before reporting.
                drop(backup);
                drop(staged);
                return Err(format!(
                    "{}: could not move existing output aside for replacement: {e}; existing output left unchanged",
                    output.display()
                ));
            }
            transaction.backup = Some(backup);
        }
        if let Err(e) = fs::rename(staged.path(), output) {
            let message = install_error(output, e, true);
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
            params,
            force,
        } => {
            let (file, metadata) = open_input(&input)?;
            preflight_output(&input, &metadata, &destination, force)?;
            let sketch = sketch_file(&input, file, params)?;
            let bytes = match &sketch {
                Sketch::MinHash(s) => s.to_bytes(),
                Sketch::FracMinHash(s) => s.to_bytes(),
                Sketch::Unsupported(_) => unreachable!("only implemented sketches are built"),
            };
            let mut transaction =
                OutputTransaction::install(&input, &metadata, &destination, force, &bytes)?;
            let report = (|| -> io::Result<()> {
                report::statistics(output, &[summary_row(&input, &sketch)])?;
                writeln!(output, "output\t{}", report::label(&destination))?;
                output.flush()
            })();
            if let Err(e) = report {
                return Err(transaction.with_rollback(format!("stdout: {e}")));
            }
            if let Some(warning) = transaction.finish() {
                eprintln!("aiips: warning: {warning}");
            }
            if let Sketch::FracMinHash(s) = &sketch {
                if s.hashes.is_empty() {
                    eprintln!(
                        "aiips: warning: {}",
                        error(
                            &destination,
                            format!(
                                "no k-mer hash is at or below the scaled={} threshold; the sketch is empty and cannot be compared",
                                s.params.scaled
                            )
                        )
                    );
                }
            }
            Ok(())
        }
        Command::Compare { inputs, options } => {
            let sketches = compare_sketches(&inputs, &options)?;
            let estimates = estimate_pair(&inputs, &sketches)?;
            report::statistics(
                output,
                &[
                    summary_row(&inputs[0], &sketches[0]),
                    summary_row(&inputs[1], &sketches[1]),
                ],
            )
            .map_err(|e| format!("stdout: {e}"))?;
            report::comparison(output, &inputs, &estimates)
                .and_then(|()| output.flush())
                .map_err(|e| format!("stdout: {e}"))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::TestDirectory;
    use aiips::sketch::fracminhash::FracMinHashParams;
    use aiips::sketch::minhash::MinHashParams;

    fn params(k: usize) -> MinHashParams {
        MinHashParams::new(NonZeroUsize::new(k).unwrap(), 16, 42).unwrap()
    }

    /// Writes a minhash sketch of `fasta` to `path` and returns it.
    fn write_sketch(path: &Path, fasta: &[u8], k: usize) -> MinHashSketch {
        let sketch = MinHashSketch::build(fasta, params(k)).unwrap();
        fs::write(path, sketch.to_bytes()).unwrap();
        sketch
    }

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
                params: SketchParams::MinHash(params(3)),
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
            OutputTransaction::install(input.path(), &source, destination.path(), true, b"new")
                .unwrap();
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
        let transaction = OutputTransaction::install(
            &input,
            &fs::metadata(&input).unwrap(),
            &destination,
            true,
            b"new sketch",
        )
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
        assert_eq!(fs::read(destination).unwrap(), b"new sketch");
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
        let result = OutputTransaction::install(&input, &metadata, &destination, false, b"new");
        assert!(result.err().unwrap().contains("output already exists"));
        assert_eq!(fs::read(destination).unwrap(), b"competing result");
    }

    #[test]
    fn compare_statistics_and_flush_failures_identify_stdout() {
        let directory = TestDirectory::new();
        let input = directory.join("in.aiips");
        write_sketch(&input, b">a\nACGT\n", 3);
        for flush_only in [false, true] {
            let result = execute(
                Command::Compare {
                    inputs: [input.clone(), input.clone()],
                    options: AlgorithmOptions::default(),
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
        let input = directory.join("in.aiips");
        let sketch = write_sketch(&input, b">a\nACGT\n", 3);
        let k = sketch.params.k;
        let mut statistics = Vec::new();
        report::statistics(
            &mut statistics,
            &[(&input, &sketch.summary, k), (&input, &sketch.summary, k)],
        )
        .unwrap();
        let result = execute(
            Command::Compare {
                inputs: [input.clone(), input],
                options: AlgorithmOptions::default(),
            },
            &mut LimitedWriter(statistics.len()),
        );
        assert_eq!(
            result.unwrap_err(),
            "stdout: injected comparison write failure"
        );
    }

    fn assert_no_staging_files(directory: &Path) {
        for entry in fs::read_dir(directory).unwrap() {
            let name = entry.unwrap().file_name();
            assert!(
                !name.to_string_lossy().starts_with(".aiips-"),
                "leftover staging file {name:?}"
            );
        }
    }

    #[test]
    fn directory_appearing_after_recheck_never_suggests_force() {
        for force in [false, true] {
            let directory = TestDirectory::new();
            let destination = directory.join("out.aiips");
            // Simulates a directory created between the recheck and install.
            fs::create_dir(&destination).unwrap();
            fs::write(destination.join("keep"), b"keep").unwrap();
            let message = OutputTransaction::install_checked(&destination, force, None, b"new")
                .err()
                .unwrap();
            assert!(message.contains("output is a directory"), "{message}");
            assert!(!message.contains("use -f"), "{message}");
            assert_eq!(fs::read(destination.join("keep")).unwrap(), b"keep");
            assert_no_staging_files(&directory);
        }
    }

    #[test]
    fn already_exists_suggests_force_only_when_not_forced() {
        let directory = TestDirectory::new();
        let destination = directory.join("out.aiips");
        fs::write(&destination, b"").unwrap();
        let unforced = install_error(&destination, io::ErrorKind::AlreadyExists.into(), false);
        assert!(unforced.contains("use -f or --force"), "{unforced}");
        let forced = install_error(&destination, io::ErrorKind::AlreadyExists.into(), true);
        assert!(!forced.contains("-f"), "{forced}");
        assert!(forced.starts_with(&destination.display().to_string()));
    }

    #[cfg(unix)]
    #[test]
    fn backup_move_failure_reports_unchanged_output_and_cleans_staging() {
        let directory = TestDirectory::new();
        let destination = directory.join("out.aiips");
        fs::write(&destination, b"original result").unwrap();
        let existing = fs::metadata(&destination).unwrap();
        // Make rename(output, backup) fail after the recheck observed a file.
        fs::remove_file(&destination).unwrap();
        let message =
            OutputTransaction::install_checked(&destination, true, Some(existing), b"new")
                .err()
                .unwrap();
        assert!(
            message.contains("could not move existing output aside"),
            "{message}"
        );
        assert!(
            message.contains("existing output left unchanged"),
            "{message}"
        );
        assert!(!destination.exists());
        assert_no_staging_files(&directory);
    }

    #[test]
    fn installs_exact_contents_with_and_without_force() {
        let directory = TestDirectory::new();
        let destination = directory.join("out.aiips");
        OutputTransaction::install_checked(&destination, false, None, b"first")
            .unwrap()
            .finish();
        assert_eq!(fs::read(&destination).unwrap(), b"first");
        let existing = fs::metadata(&destination).unwrap();
        OutputTransaction::install_checked(&destination, true, Some(existing), b"second")
            .unwrap()
            .finish();
        assert_eq!(fs::read(&destination).unwrap(), b"second");
        assert_no_staging_files(&directory);
    }

    #[test]
    fn compare_halts_unless_both_sketches_share_an_implemented_algorithm() {
        let directory = TestDirectory::new();
        let minhash = directory.join("m.aiips");
        let other = directory.join("o.aiips");
        let frac = directory.join("f.aiips");
        write_sketch(&minhash, b">a\nACGTACGT\n", 3);
        fs::write(&other, b"aiips-sketch v1 algorithm=bottomk\n").unwrap();
        let p = FracMinHashParams::new(NonZeroUsize::new(3).unwrap(), 1, 42).unwrap();
        fs::write(
            &frac,
            FracMinHashSketch::build(b">a\nACGTACGT\n".as_slice(), p)
                .unwrap()
                .to_bytes(),
        )
        .unwrap();
        for (inputs, names) in [
            ([minhash.clone(), frac.clone()], "minhash, "),
            ([frac.clone(), minhash.clone()], "fracminhash, "),
            ([frac.clone(), other.clone()], "bottomk)"),
        ] {
            let mut stdout = Vec::new();
            let error = execute(
                Command::Compare {
                    inputs,
                    options: AlgorithmOptions::default(),
                },
                &mut stdout,
            )
            .unwrap_err();
            assert!(error.contains("two sketches of one algorithm"), "{error}");
            assert!(error.contains(names), "{error}");
            assert!(stdout.is_empty());
        }
        for inputs in [
            [minhash.clone(), other.clone()],
            [other.clone(), minhash.clone()],
            [other.clone(), other.clone()],
        ] {
            let mut stdout = Vec::new();
            let error = execute(
                Command::Compare {
                    inputs,
                    options: AlgorithmOptions::default(),
                },
                &mut stdout,
            )
            .unwrap_err();
            assert!(
                error.contains("comparison requires two sketches of one algorithm"),
                "{error}"
            );
            assert!(error.contains("bottomk"), "{error}");
            assert!(stdout.is_empty());
        }
    }

    #[test]
    fn compare_halts_on_parameter_mismatch_before_printing() {
        let directory = TestDirectory::new();
        let a = directory.join("a.aiips");
        let b = directory.join("b.aiips");
        write_sketch(&a, b">a\nACGTACGT\n", 3);
        write_sketch(&b, b">a\nACGTACGT\n", 4);
        let mut stdout = Vec::new();
        let error = execute(
            Command::Compare {
                inputs: [a, b],
                options: AlgorithmOptions::default(),
            },
            &mut stdout,
        )
        .unwrap_err();
        assert!(
            error.contains("sketch parameters differ: k (3 vs 4)"),
            "{error}"
        );
        assert!(stdout.is_empty());
    }

    #[test]
    fn oversized_sketch_files_are_refused_without_parsing() {
        let directory = TestDirectory::new();
        let big = directory.join("big.aiips");
        let file = File::create(&big).unwrap();
        file.set_len(sketch_format::max_file_bytes(b"") + 1)
            .unwrap();
        let error = read_sketch(&big).unwrap_err();
        assert!(error.contains("too large"), "{error}");
    }

    #[test]
    fn fracminhash_files_are_bounded_by_their_declared_size() {
        let directory = TestDirectory::new();
        let path = directory.join("f.aiips");
        let header =
            b"aiips-sketch v1 algorithm=fracminhash k=3 scaled=1 seed=42 hash=murmur3-fmix64 canonical=true size=0\n";
        let mut bytes = header.to_vec();
        bytes.resize(sketch_format::HEADER_BYTES as usize + 1, b'0');
        fs::write(&path, &bytes).unwrap();
        assert!(read_sketch(&path).unwrap_err().contains("too large"));
        // A large declared size admits a file above the MinHash bound.
        let size = minhash::MAX_NUM_HASHES as u64 + 1000;
        let p = FracMinHashParams::new(NonZeroUsize::new(3).unwrap(), 1, 42).unwrap();
        let sketch = FracMinHashSketch {
            params: p,
            summary: Summary {
                records: 1,
                bases: 10,
                non_acgt_bases: 0,
                kmers: 8,
            },
            hashes: (0..size).map(|i| u64::MAX - size + i).collect(),
        };
        let bytes = sketch.to_bytes();
        assert!(bytes.len() as u64 > sketch_format::max_file_bytes(b""));
        fs::write(&path, &bytes).unwrap();
        match read_sketch(&path).unwrap() {
            Sketch::FracMinHash(parsed) => assert_eq!(parsed, sketch),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn compare_halts_on_an_empty_fracminhash_sketch() {
        let directory = TestDirectory::new();
        let empty = directory.join("e.aiips");
        let full = directory.join("f.aiips");
        let build = |scaled| {
            let p = FracMinHashParams::new(NonZeroUsize::new(3).unwrap(), scaled, 42).unwrap();
            FracMinHashSketch::build(b">a\nACGTAC\n".as_slice(), p).unwrap()
        };
        let mut nonempty = build(u64::MAX);
        assert!(nonempty.hashes.is_empty());
        fs::write(&empty, nonempty.to_bytes()).unwrap();
        nonempty.hashes = vec![0];
        fs::write(&full, nonempty.to_bytes()).unwrap();
        for inputs in [[empty.clone(), full.clone()], [full.clone(), empty.clone()]] {
            let mut stdout = Vec::new();
            let error = execute(
                Command::Compare {
                    inputs,
                    options: AlgorithmOptions::default(),
                },
                &mut stdout,
            )
            .unwrap_err();
            assert!(
                error.starts_with(&format!("{}: ", empty.display())),
                "{error}"
            );
            assert!(error.contains(fracminhash::EMPTY), "{error}");
            assert!(stdout.is_empty());
        }
    }

    #[test]
    fn retry_unique_is_bounded_and_retries_only_collisions() {
        let mut calls = 0;
        let result: io::Result<()> = retry_unique(|| {
            calls += 1;
            Err(io::ErrorKind::AlreadyExists.into())
        });
        assert_eq!(result.unwrap_err().kind(), io::ErrorKind::AlreadyExists);
        assert_eq!(calls, MAX_UNIQUE_ATTEMPTS);

        let mut calls = 0;
        let result = retry_unique(|| {
            calls += 1;
            if calls < 3 {
                Err(io::ErrorKind::AlreadyExists.into())
            } else {
                Ok(calls)
            }
        });
        assert_eq!(result.unwrap(), 3);

        let mut calls = 0;
        let result: io::Result<()> = retry_unique(|| {
            calls += 1;
            Err(io::ErrorKind::PermissionDenied.into())
        });
        assert_eq!(result.unwrap_err().kind(), io::ErrorKind::PermissionDenied);
        assert_eq!(calls, 1);
    }
}
