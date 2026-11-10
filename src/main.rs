mod cli;
mod commands;
mod report;

#[cfg(test)]
#[path = "../tests/common/mod.rs"]
mod test_support;

use std::io::{self, Write};
use std::process::ExitCode;

fn run() -> Result<(), String> {
    let command = cli::parse(std::env::args_os().skip(1))?;
    let mut output = io::BufWriter::new(io::stdout().lock());
    let result = commands::execute(command, &mut output);
    finish_stdout(output, result)
}

/// On failure, discard buffered stdout instead of letting BufWriter's Drop
/// retry the flush, which could print a success report for an output that
/// has already been rolled back.
fn finish_stdout<W: Write>(
    output: io::BufWriter<W>,
    result: Result<(), String>,
) -> Result<(), String> {
    if result.is_err() {
        drop(output.into_parts());
    }
    result
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("aiips: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::TestDirectory;
    use std::cell::RefCell;
    use std::num::NonZeroUsize;
    use std::rc::Rc;

    /// Fails its first write with a transient error, then accepts everything.
    struct TransientWriter {
        failed: bool,
        sink: Rc<RefCell<Vec<u8>>>,
    }

    impl Write for TransientWriter {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            if !self.failed {
                self.failed = true;
                return Err(io::ErrorKind::WouldBlock.into());
            }
            self.sink.borrow_mut().extend_from_slice(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    fn failed_sketch(discard: bool) -> (Vec<u8>, bool) {
        let directory = TestDirectory::new();
        let input = directory.join("in.fa");
        let destination = directory.join("out.aiips");
        std::fs::write(&input, b">a\nACGT\n").unwrap();
        let sink = Rc::new(RefCell::new(Vec::new()));
        let mut output = io::BufWriter::new(TransientWriter {
            failed: false,
            sink: Rc::clone(&sink),
        });
        let result = commands::execute(
            cli::Command::Sketch {
                input,
                output: destination.clone(),
                k: NonZeroUsize::new(3).unwrap(),
                force: false,
            },
            &mut output,
        );
        assert!(result.as_ref().unwrap_err().starts_with("stdout: "));
        if discard {
            assert!(finish_stdout(output, result).is_err());
        } else {
            drop(output);
        }
        let printed = sink.borrow().clone();
        (printed, destination.exists())
    }

    #[test]
    fn failed_report_is_not_printed_after_rollback() {
        let (printed, exists) = failed_sketch(true);
        assert!(!exists, "output must be rolled back");
        assert!(printed.is_empty(), "{}", String::from_utf8_lossy(&printed));
    }

    #[test]
    fn plain_bufwriter_drop_would_print_rolled_back_report() {
        // Control: demonstrates the hazard finish_stdout guards against.
        let (printed, exists) = failed_sketch(false);
        assert!(!exists);
        assert!(String::from_utf8_lossy(&printed).contains("output\t"));
    }
}
