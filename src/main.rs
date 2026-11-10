use std::env;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Write};
use std::num::NonZeroUsize;
use std::path::Path;
use std::process::ExitCode;

use aiips::fastx::summarize;

const HELP: &str = "Usage: aiips [-k <length>] <input>...\n\n\
    Summarize FASTA/FASTQ files (plain or gzip); use - for stdin.\n\
    -k, --kmer-size <length>  Positive k-mer length [default: 21]\n\
    -h, --help                Print help\n\n\
    Counts overlapping A/C/G/T k-mers, accepting either case.\n\
    Other bases break k-mer runs. FASTQ qualities are not filtered.";

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args_os().skip(1);
    let mut k = NonZeroUsize::new(21).unwrap();
    let mut inputs = Vec::<OsString>::new();
    let mut positional_only = false;
    while let Some(arg) = args.next() {
        if !positional_only && (arg == "-h" || arg == "--help") {
            println!("{HELP}");
            return Ok(());
        } else if !positional_only && arg == "--" {
            positional_only = true;
        } else if !positional_only && (arg == "-k" || arg == "--kmer-size") {
            k = args
                .next()
                .and_then(|value| value.to_str().and_then(|s| s.parse().ok()))
                .ok_or("k-mer length must be a positive integer fitting usize")?;
        } else if !positional_only && arg != "-" && arg.to_string_lossy().starts_with('-') {
            return Err(format!("unknown option: {}", arg.to_string_lossy()).into());
        } else {
            inputs.push(arg);
        }
    }
    if inputs.is_empty() {
        return Err(format!("no input files supplied\n{HELP}").into());
    }
    if inputs.iter().filter(|input| *input == "-").count() > 1 {
        return Err("stdin (-) can only be supplied once".into());
    }

    let mut output = io::BufWriter::new(io::stdout().lock());
    writeln!(
        output,
        "input\trecords\tbases\tnon_acgt_bases\tk\tvalid_kmers"
    )?;
    for input in inputs {
        let path = Path::new(&input);
        let result = if input == "-" {
            summarize(io::stdin(), k).map_err(|error| error.to_string())
        } else {
            File::open(path)
                .map_err(|error| error.to_string())
                .and_then(|file| summarize(file, k).map_err(|error| error.to_string()))
        };
        let summary = result.map_err(|error| format!("{}: {error}", path.display()))?;
        // Escape delimiters so unusual filenames cannot create extra TSV columns.
        let label = path
            .to_string_lossy()
            .replace('\\', "\\\\")
            .replace('\t', "\\t")
            .replace('\r', "\\r")
            .replace('\n', "\\n");
        writeln!(
            output,
            "{label}\t{}\t{}\t{}\t{k}\t{}",
            summary.records, summary.bases, summary.non_acgt_bases, summary.kmers
        )?;
    }
    output.flush()?;
    Ok(())
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
