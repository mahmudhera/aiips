use std::ffi::OsString;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

pub const HELP: &str = "Usage: aiips <COMMAND>\n\nCommands:\n  sketch   Summarize FASTA/FASTQ and create an empty .aiips placeholder\n  compare  Summarize two inputs and print placeholder metrics\n\nUse aiips <COMMAND> --help for options. Stdin is not supported.";
const SKETCH_HELP: &str = "Usage: aiips sketch INPUT -o OUTPUT [-k LENGTH] [-f]\n\n  -o, --output PATH        Output path (.aiips appended if missing)\n  -k, --kmer-size LENGTH   Positive k-mer length [default: 21]\n  -f, --force              Replace an existing output\n  -h, --help               Print help\n\nReads plain/gzip FASTA or FASTQ. Creates an empty placeholder only.\nUse -- before positional filenames beginning with -. Stdin is not supported.";
const COMPARE_HELP: &str = "Usage: aiips compare INPUT_A INPUT_B [-k LENGTH]\n\n  -k, --kmer-size LENGTH   Positive k-mer length [default: 21]\n  -h, --help               Print help\n\n.aiips inputs must be empty placeholder files; other inputs are plain/gzip\nFASTA or FASTQ. Mixed pairs are supported. All six metric values are -1\n(unavailable), with containment and AF reported in both directions.\nUse -- before positional filenames beginning with -. Stdin is not supported.";

pub enum Command {
    Help(&'static str),
    Sketch {
        input: PathBuf,
        output: PathBuf,
        k: NonZeroUsize,
        force: bool,
    },
    Compare {
        inputs: [PathBuf; 2],
        k: NonZeroUsize,
    },
}

pub fn is_sketch(path: &Path) -> bool {
    path.extension()
        .is_some_and(|extension| extension.eq_ignore_ascii_case("aiips"))
}

fn path(value: OsString) -> Result<PathBuf, String> {
    if value.is_empty() || value == "-" {
        return Err("a filesystem path is required; stdin (-) is not supported".into());
    }
    Ok(value.into())
}

pub fn parse(args: impl IntoIterator<Item = OsString>) -> Result<Command, String> {
    let mut args = args.into_iter();
    let first = args.next().ok_or_else(|| HELP.to_string())?;
    if first == "-h" || first == "--help" {
        return Ok(Command::Help(HELP));
    }
    let sketch = match first.to_str() {
        Some("sketch") => true,
        Some("compare") => false,
        _ => {
            return Err(format!(
                "unknown command: {}\n{HELP}",
                first.to_string_lossy()
            ))
        }
    };
    let help = if sketch { SKETCH_HELP } else { COMPARE_HELP };
    let mut inputs = Vec::new();
    let mut output = None;
    let mut k = NonZeroUsize::new(21).unwrap();
    let mut force = false;
    let mut positional_only = false;
    while let Some(arg) = args.next() {
        if !positional_only && (arg == "-h" || arg == "--help") {
            return Ok(Command::Help(help));
        } else if !positional_only && arg == "--" {
            positional_only = true;
        } else if !positional_only && (arg == "-k" || arg == "--kmer-size") {
            k = args
                .next()
                .and_then(|v| v.to_str().and_then(|s| s.parse().ok()))
                .ok_or_else(|| {
                    format!("k-mer length must be a positive integer fitting usize\n{help}")
                })?;
        } else if !positional_only && sketch && (arg == "-o" || arg == "--output") {
            if output.is_some() {
                return Err(format!("output may only be specified once\n{help}"));
            }
            let value = args
                .next()
                .ok_or_else(|| format!("missing output path\n{help}"))?;
            if value != "-" && value.to_string_lossy().starts_with('-') {
                return Err(format!("missing output path before {}; use ./ for output names beginning with -\n{help}", value.to_string_lossy()));
            }
            output = Some(path(value)?);
        } else if !positional_only && sketch && (arg == "-f" || arg == "--force") {
            force = true;
        } else if !positional_only && arg != "-" && arg.to_string_lossy().starts_with('-') {
            return Err(format!("unknown option: {}\n{help}", arg.to_string_lossy()));
        } else {
            inputs.push(path(arg)?);
        }
    }
    let expected = if sketch { 1 } else { 2 };
    if inputs.len() != expected {
        return Err(format!(
            "expected exactly {expected} input file(s), received {}\n{help}",
            inputs.len()
        ));
    }
    if sketch {
        let mut output = output.ok_or_else(|| format!("-o/--output is required\n{help}"))?;
        if !is_sketch(&output) {
            let mut name = output.into_os_string();
            name.push(".aiips");
            output = name.into();
        }
        Ok(Command::Sketch {
            input: inputs.remove(0),
            output,
            k,
            force,
        })
    } else {
        Ok(Command::Compare {
            inputs: inputs.try_into().unwrap(),
            k,
        })
    }
}
