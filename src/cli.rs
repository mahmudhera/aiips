use std::ffi::OsString;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use aiips::sketch::minhash::{
    MinHashParams, DEFAULT_K, DEFAULT_NUM_HASHES, DEFAULT_SEED, MAX_NUM_HASHES, MIN_NUM_HASHES,
};
use aiips::sketch::Algorithm;

pub const HELP: &str = "Usage: aiips <COMMAND>\n\nCommands:\n  sketch   Sketch a FASTA/FASTQ file into a .aiips file\n  compare  Estimate similarity between two .aiips sketches\n\nUse aiips <COMMAND> --help for options. Stdin is not supported.";
const SKETCH_HELP: &str = "Usage: aiips sketch INPUT -o OUTPUT [-a ALGORITHM] [-k LENGTH] [-s NUM] [--seed SEED] [-f]\n\n  -o, --output PATH         Output path (.aiips appended if missing)\n  -a, --algorithm NAME      Sketch algorithm; implemented: minhash [default: minhash]\n  -k, --kmer-size LENGTH    k-mer length, 1 to 32 for minhash [default: 21]\n  -s, --num-hashes NUM      minhash: number of hash functions, 2 to 1000000 [default: 1000]\n      --seed SEED           minhash: unsigned 64-bit hash seed [default: 42]\n  -f, --force               Replace an existing output\n  -h, --help                Print help\n\nReads plain/gzip FASTA or FASTQ; k-mers are canonical (strand-independent).\nReserved, not yet implemented: fracminhash, bottomk, minimizer, open-syncmer,\nclosed-syncmer, strobemer.\nUse -- before positional filenames beginning with -. Stdin is not supported.";
const COMPARE_HELP: &str = "Usage: aiips compare INPUT_A INPUT_B [-k LENGTH]\n\n  -k, --kmer-size LENGTH   k-mer length for FASTA/FASTQ inputs [default: 21, or the sketch's k]\n  -h, --help               Print help\n\nEach input is a .aiips sketch or a plain/gzip FASTA/FASTQ file; mixed pairs are\nsupported. Sequence inputs are sketched in memory with minhash, using the other\ninput's sketch parameters when it is a sketch (-k must then match), otherwise\nk from -k, 1000 hash functions and seed 42. Sketches must be minhash with\nidentical parameters. Prints Jaccard, containment and ANI estimates; AF is -1:\nit is not estimable from MinHash.\nUse -- before positional filenames beginning with -. Stdin is not supported.";

pub enum SketchParams {
    MinHash(MinHashParams),
}

pub enum Command {
    Help(&'static str),
    Sketch {
        input: PathBuf,
        output: PathBuf,
        params: SketchParams,
        force: bool,
    },
    Compare {
        inputs: [PathBuf; 2],
        /// For FASTA/FASTQ inputs; must match any sketch input's k.
        k: Option<NonZeroUsize>,
    },
}

pub fn is_sketch(path: &Path) -> bool {
    path.file_name()
        .is_some_and(|name| name.eq_ignore_ascii_case(".aiips"))
        || path
            .extension()
            .is_some_and(|extension| extension.eq_ignore_ascii_case("aiips"))
}

fn output_path(value: OsString) -> Result<PathBuf, String> {
    let output = path(value)?;
    // Components normalize trailing separators and '.', so check the original spelling too.
    let spelling = output.as_os_str().to_string_lossy();
    if spelling.chars().last().is_some_and(std::path::is_separator)
        || output.file_name().is_none()
        || matches!(
            spelling.rsplit(std::path::is_separator).next(),
            Some("." | "..")
        )
    {
        return Err(format!(
            "{}: output must name a file, not a directory",
            output.display()
        ));
    }
    // Existence and type are checked on the final, suffixed destination by
    // the sketch preflight; an unsuffixed directory name does not block it.
    if is_sketch(&output) {
        return Ok(output);
    }
    let mut name = output.file_name().unwrap().to_os_string();
    name.push(".aiips");
    Ok(output.with_file_name(name))
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
    let remaining: Vec<_> = args.collect();
    if remaining
        .iter()
        .take_while(|arg| *arg != "--")
        .any(|arg| arg == "-h" || arg == "--help")
    {
        return Ok(Command::Help(help));
    }
    let mut args = remaining.into_iter();
    let mut inputs = Vec::new();
    let mut output = None;
    let mut k: Option<NonZeroUsize> = None;
    let mut algorithm: Option<Algorithm> = None;
    let mut num_hashes: Option<usize> = None;
    let mut seed: Option<u64> = None;
    let mut force = false;
    let mut positional_only = false;
    while let Some(arg) = args.next() {
        let option = !positional_only && arg != "-" && arg.to_string_lossy().starts_with('-');
        let sketch_option = option && sketch;
        if !positional_only && arg == "--" {
            positional_only = true;
        } else if option && (arg == "-k" || arg == "--kmer-size") {
            once(&k, "k-mer length", help)?;
            k = Some(number(
                args.next(),
                "k-mer length must be a positive integer fitting usize",
                help,
            )?);
        } else if sketch_option && (arg == "-a" || arg == "--algorithm") {
            once(&algorithm, "algorithm", help)?;
            algorithm = Some(parse_algorithm(args.next(), help)?);
        } else if sketch_option && (arg == "-s" || arg == "--num-hashes") {
            once(&num_hashes, "number of hash functions", help)?;
            num_hashes = Some(number(
                args.next(),
                &format!("number of hash functions must be an integer from {MIN_NUM_HASHES} to {MAX_NUM_HASHES}"),
                help,
            )?);
        } else if sketch_option && arg == "--seed" {
            once(&seed, "seed", help)?;
            seed = Some(number(
                args.next(),
                "seed must be an unsigned 64-bit integer",
                help,
            )?);
        } else if sketch_option && (arg == "-o" || arg == "--output") {
            if output.is_some() {
                return Err(format!("output may only be specified once\n{help}"));
            }
            let value = args
                .next()
                .ok_or_else(|| format!("missing output path\n{help}"))?;
            let value = path(value)?;
            if value.as_os_str().to_string_lossy().starts_with('-') {
                return Err(format!("missing output path before {}; use ./ for output names beginning with -\n{help}", value.display()));
            }
            output = Some(value);
        } else if sketch_option && (arg == "-f" || arg == "--force") {
            force = true;
        } else if option {
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
    if !sketch {
        return Ok(Command::Compare {
            inputs: inputs.try_into().unwrap(),
            k,
        });
    }
    let output = output_path(
        output
            .ok_or_else(|| format!("-o/--output is required\n{help}"))?
            .into_os_string(),
    )?;
    // Omitting -a keeps pre-algorithm commands working; minhash is the default.
    let params = match algorithm.unwrap_or(Algorithm::MinHash) {
        Algorithm::MinHash => SketchParams::MinHash(MinHashParams::new(
            k.unwrap_or(NonZeroUsize::new(DEFAULT_K).unwrap()),
            num_hashes.unwrap_or(DEFAULT_NUM_HASHES),
            seed.unwrap_or(DEFAULT_SEED),
        )?),
        other => unreachable!("{} rejected while parsing", other.name()),
    };
    Ok(Command::Sketch {
        input: inputs.remove(0),
        output,
        params,
        force,
    })
}

fn once<T>(value: &Option<T>, name: &str, help: &str) -> Result<(), String> {
    match value {
        Some(_) => Err(format!("{name} may only be specified once\n{help}")),
        None => Ok(()),
    }
}

fn number<T: std::str::FromStr>(
    value: Option<OsString>,
    message: &str,
    help: &str,
) -> Result<T, String> {
    value
        .and_then(|v| v.to_str().and_then(|s| s.parse().ok()))
        .ok_or_else(|| format!("{message}\n{help}"))
}

fn parse_algorithm(value: Option<OsString>, help: &str) -> Result<Algorithm, String> {
    let value = value.ok_or_else(|| format!("missing algorithm name\n{help}"))?;
    let name = value.to_string_lossy();
    match Algorithm::from_name(&name) {
        Some(algorithm) if algorithm.is_implemented() => Ok(algorithm),
        Some(_) => Err(format!(
            "{name} is not implemented yet; implemented: minhash\n{help}"
        )),
        None => Err(format!(
            "unknown algorithm: {name}; known: {}\n{help}",
            Algorithm::names()
        )),
    }
}
