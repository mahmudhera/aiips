use std::ffi::OsString;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use aiips::sketch::aiips::{self as aiips_sketch, AiipsParams, WalkOptions};
use aiips::sketch::fracminhash::{self, FracMinHashParams};
use aiips::sketch::minhash::{self, MinHashParams, MAX_NUM_HASHES, MIN_NUM_HASHES};
use aiips::sketch::{Algorithm, Sketch};

pub const HELP: &str = "Usage: aiips <COMMAND>\n\nCommands:\n  sketch   Sketch a FASTA/FASTQ file into a .aiips file\n  compare  Estimate similarity between two .aiips sketches\n\nUse aiips <COMMAND> --help for options. Stdin is not supported.";
const SKETCH_HELP: &str = "Usage: aiips sketch INPUT -o OUTPUT [-a ALGORITHM] [-k LENGTH] [-s NUM] [--scaled S] [--seed SEED] [-f]\n\n  -o, --output PATH         Output path (.aiips appended if missing)\n  -a, --algorithm NAME      Sketch algorithm: minhash, fracminhash or aiips [default: minhash]\n  -k, --kmer-size LENGTH    k-mer length, 1 to 32 [default: 21]\n  -s, --num-hashes NUM      minhash: number of hash functions, 2 to 1000000 [default: 1000]\n      --scaled S            fracminhash: scale factor, keeps about 1/S of k-mers,\n                            1 to 18446744073709551615 [default: 1000]\n      --seed SEED           unsigned 64-bit hash seed [default: 42]\n  -f, --force               Replace an existing output\n  -h, --help                Print help\n\naiips options (instead of -k, -s and --scaled):\n      --anchor-k LENGTH     anchor k-mer length k1, 2 to 32 [default: 21]\n      --anchor-scaled S     keep about 1/S of k1-mers as anchors [default: 500]\n      --probe-k LENGTH      probe k-mer length k2, below k1 [default: 15]\n      --probe-scaled S      keep about 1/S of k2-mers as probes [default: 100]\n\nReads plain/gzip FASTA or FASTQ; k-mers are canonical (strand-independent).\nReserved, not yet implemented: bottomk, minimizer, open-syncmer, closed-syncmer,\nstrobemer.\nUse -- before positional filenames beginning with -. Stdin is not supported.";
const COMPARE_HELP: &str = "Usage: aiips compare INPUT_A INPUT_B [-a ALGORITHM] [-k LENGTH] [-s NUM] [--scaled S] [--seed SEED]\n\n  -a, --algorithm NAME     Algorithm for FASTA/FASTQ inputs: minhash, fracminhash or aiips\n                           [default: minhash, or the sketch's algorithm]\n  -k, --kmer-size LENGTH   k-mer length [default: 21, or the sketch's k]\n  -s, --num-hashes NUM     minhash: number of hash functions [default: 1000]\n      --scaled S           fracminhash: scale factor [default: 1000]\n      --seed SEED          unsigned 64-bit hash seed [default: 42]\n  -h, --help               Print help\n\naiips options: --anchor-k, --anchor-scaled, --probe-k, --probe-scaled (as for\nsketch), and\n      --max-gap G          consecutive low islands tolerated by a walk [default: 50]\n      --min-island-ani X   an island is high enough when its probe containment is\n                           at least X^k2, 0 to 1 [default: 0.8]\n      --symmetric-probes   count a probe only if it also avoids the other genome's\n                           anchors (corrects a downward ANI bias)\n      --bidirectional      walk backward from each seed as well as forward\n\nEach input is a .aiips sketch or a plain/gzip FASTA/FASTQ file; mixed pairs are\nsupported. Sequence inputs are sketched in memory, using the other input's\nalgorithm and parameters when it is a sketch (given options must then match),\notherwise the options above. Sketches must share one algorithm and identical\nparameters. minhash and fracminhash print Jaccard, containment and ANI\nestimates with AF -1 (not estimable from them); aiips prints ANI and AF.\nUse -- before positional filenames beginning with -. Stdin is not supported.";

/// Validated parameters for building one sketch.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SketchParams {
    MinHash(MinHashParams),
    FracMinHash(FracMinHashParams),
    Aiips(AiipsParams),
}

impl SketchParams {
    /// The parameters recorded in an implemented sketch.
    pub fn of(sketch: &Sketch) -> Option<Self> {
        match sketch {
            Sketch::MinHash(s) => Some(SketchParams::MinHash(s.params.clone())),
            Sketch::FracMinHash(s) => Some(SketchParams::FracMinHash(s.params.clone())),
            Sketch::Aiips(s) => Some(SketchParams::Aiips(s.params.clone())),
            Sketch::Unsupported(_) => None,
        }
    }

    pub fn k(&self) -> NonZeroUsize {
        match self {
            SketchParams::MinHash(p) => p.k,
            SketchParams::FracMinHash(p) => p.k,
            SketchParams::Aiips(p) => p.probe_k,
        }
    }
}

/// Algorithm options as given on the command line; `None` means not given.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AlgorithmOptions {
    pub algorithm: Option<Algorithm>,
    pub k: Option<NonZeroUsize>,
    pub num_hashes: Option<usize>,
    pub scaled: Option<u64>,
    pub seed: Option<u64>,
    pub anchor_k: Option<NonZeroUsize>,
    pub anchor_scaled: Option<u64>,
    pub probe_k: Option<NonZeroUsize>,
    pub probe_scaled: Option<u64>,
    /// Comparison-only aiips options.
    pub max_gap: Option<u32>,
    pub min_island_ani: Option<f64>,
    pub symmetric_probes: bool,
    pub bidirectional: bool,
}

impl AlgorithmOptions {
    /// Parameters from the given options and defaults. Omitting -a keeps
    /// pre-algorithm commands working; minhash is the default.
    pub fn params(&self) -> Result<SketchParams, String> {
        self.params_for(self.algorithm.unwrap_or(Algorithm::MinHash))
    }

    fn aiips_flag_given(&self) -> Option<&'static str> {
        [
            (self.anchor_k.is_some(), "--anchor-k"),
            (self.anchor_scaled.is_some(), "--anchor-scaled"),
            (self.probe_k.is_some(), "--probe-k"),
            (self.probe_scaled.is_some(), "--probe-scaled"),
            (self.max_gap.is_some(), "--max-gap"),
            (self.min_island_ani.is_some(), "--min-island-ani"),
            (self.symmetric_probes, "--symmetric-probes"),
            (self.bidirectional, "--bidirectional"),
        ]
        .into_iter()
        .find_map(|(given, flag)| given.then_some(flag))
    }

    /// Walk options for an aiips comparison, with defaults.
    pub fn walk(&self) -> WalkOptions {
        let defaults = WalkOptions::default();
        WalkOptions {
            max_gap: self.max_gap.unwrap_or(defaults.max_gap),
            min_island_ani: self.min_island_ani.unwrap_or(defaults.min_island_ani),
            symmetric_probes: self.symmetric_probes,
            bidirectional: self.bidirectional,
        }
    }

    /// Refuses options that do not apply to `algorithm`, then fills defaults.
    fn params_for(&self, algorithm: Algorithm) -> Result<SketchParams, String> {
        let k = self.k;
        let seed = self.seed;
        if algorithm != Algorithm::Aiips {
            if let Some(flag) = self.aiips_flag_given() {
                return Err(format!(
                    "{flag} applies only to aiips, not {}",
                    algorithm.name()
                ));
            }
        }
        match algorithm {
            Algorithm::Aiips => {
                for (given, flag) in [
                    (self.k.is_some(), "-k/--kmer-size"),
                    (self.num_hashes.is_some(), "-s/--num-hashes"),
                    (self.scaled.is_some(), "--scaled"),
                ] {
                    if given {
                        return Err(format!(
                            "{flag} does not apply to aiips; use --anchor-k, --anchor-scaled, --probe-k and --probe-scaled"
                        ));
                    }
                }
                let size = |k: Option<NonZeroUsize>, default: usize| {
                    k.unwrap_or(NonZeroUsize::new(default).unwrap())
                };
                Ok(SketchParams::Aiips(AiipsParams::new(
                    size(self.anchor_k, aiips_sketch::DEFAULT_ANCHOR_K),
                    self.anchor_scaled
                        .unwrap_or(aiips_sketch::DEFAULT_ANCHOR_SCALED),
                    size(self.probe_k, aiips_sketch::DEFAULT_PROBE_K),
                    self.probe_scaled
                        .unwrap_or(aiips_sketch::DEFAULT_PROBE_SCALED),
                    seed.unwrap_or(aiips_sketch::DEFAULT_SEED),
                )?))
            }
            Algorithm::MinHash => {
                if self.scaled.is_some() {
                    return Err("--scaled applies only to fracminhash, not minhash".into());
                }
                Ok(SketchParams::MinHash(MinHashParams::new(
                    k.unwrap_or(NonZeroUsize::new(minhash::DEFAULT_K).unwrap()),
                    self.num_hashes.unwrap_or(minhash::DEFAULT_NUM_HASHES),
                    seed.unwrap_or(minhash::DEFAULT_SEED),
                )?))
            }
            Algorithm::FracMinHash => {
                if self.num_hashes.is_some() {
                    return Err("-s/--num-hashes applies only to minhash, not fracminhash".into());
                }
                Ok(SketchParams::FracMinHash(FracMinHashParams::new(
                    k.unwrap_or(NonZeroUsize::new(fracminhash::DEFAULT_K).unwrap()),
                    self.scaled.unwrap_or(fracminhash::DEFAULT_SCALED),
                    seed.unwrap_or(fracminhash::DEFAULT_SEED),
                )?))
            }
            other => unreachable!("{} rejected while parsing", other.name()),
        }
    }

    /// The parameters of `template`, a sketch read from `source`, after
    /// checking that every given option applies to it and matches it.
    pub fn matching(&self, template: &SketchParams, source: &Path) -> Result<SketchParams, String> {
        let recorded = format!("recorded in {}", source.display());
        let algorithm = match template {
            SketchParams::MinHash(_) => Algorithm::MinHash,
            SketchParams::FracMinHash(_) => Algorithm::FracMinHash,
            SketchParams::Aiips(_) => Algorithm::Aiips,
        };
        if let Some(given) = self.algorithm.filter(|a| *a != algorithm) {
            return Err(format!(
                "-a {} does not match algorithm={} {}",
                given.name(),
                algorithm.name(),
                recorded
            ));
        }
        // Validates ranges and applicability of the given options.
        self.params_for(algorithm)?;
        let check = |flag: &str, given: Option<String>, name: &str, value: String| match given
            .filter(|given| *given != value)
        {
            Some(given) => Err(format!(
                "{flag} {given} does not match {name}={value} {recorded}"
            )),
            None => Ok(()),
        };
        if !matches!(template, SketchParams::Aiips(_)) {
            check(
                "-k",
                self.k.map(|k| k.to_string()),
                "k",
                template.k().to_string(),
            )?;
        }
        match template {
            SketchParams::Aiips(p) => {
                let pairs = [
                    (
                        "--anchor-k",
                        self.anchor_k.map(|k| k.to_string()),
                        "anchor_k",
                        p.anchor_k.to_string(),
                    ),
                    (
                        "--anchor-scaled",
                        self.anchor_scaled.map(|s| s.to_string()),
                        "anchor_scaled",
                        p.anchor_scaled.to_string(),
                    ),
                    (
                        "--probe-k",
                        self.probe_k.map(|k| k.to_string()),
                        "probe_k",
                        p.probe_k.to_string(),
                    ),
                    (
                        "--probe-scaled",
                        self.probe_scaled.map(|s| s.to_string()),
                        "probe_scaled",
                        p.probe_scaled.to_string(),
                    ),
                    (
                        "--seed",
                        self.seed.map(|s| s.to_string()),
                        "seed",
                        p.seed.to_string(),
                    ),
                ];
                for (flag, given, name, value) in pairs {
                    check(flag, given, name, value)?;
                }
            }
            SketchParams::MinHash(p) => {
                check(
                    "-s",
                    self.num_hashes.map(|n| n.to_string()),
                    "num_hashes",
                    p.num_hashes.to_string(),
                )?;
                check(
                    "--seed",
                    self.seed.map(|s| s.to_string()),
                    "seed",
                    p.seed.to_string(),
                )?;
            }
            SketchParams::FracMinHash(p) => {
                check(
                    "--scaled",
                    self.scaled.map(|s| s.to_string()),
                    "scaled",
                    p.scaled.to_string(),
                )?;
                check(
                    "--seed",
                    self.seed.map(|s| s.to_string()),
                    "seed",
                    p.seed.to_string(),
                )?;
            }
        }
        Ok(template.clone())
    }
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
        /// For FASTA/FASTQ inputs; given options must match any sketch input.
        options: AlgorithmOptions,
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
    let mut options = AlgorithmOptions::default();
    let mut force = false;
    let mut positional_only = false;
    while let Some(arg) = args.next() {
        let option = !positional_only && arg != "-" && arg.to_string_lossy().starts_with('-');
        let sketch_option = option && sketch;
        if !positional_only && arg == "--" {
            positional_only = true;
        } else if option && (arg == "-k" || arg == "--kmer-size") {
            once(&options.k, "k-mer length", help)?;
            options.k = Some(number(
                args.next(),
                "k-mer length must be a positive integer fitting usize",
                help,
            )?);
        } else if option && (arg == "-a" || arg == "--algorithm") {
            once(&options.algorithm, "algorithm", help)?;
            options.algorithm = Some(parse_algorithm(args.next(), help)?);
        } else if option && (arg == "-s" || arg == "--num-hashes") {
            once(&options.num_hashes, "number of hash functions", help)?;
            options.num_hashes = Some(number(
                args.next(),
                &format!("number of hash functions must be an integer from {MIN_NUM_HASHES} to {MAX_NUM_HASHES}"),
                help,
            )?);
        } else if option && arg == "--scaled" {
            once(&options.scaled, "scaled", help)?;
            let scaled: u64 = number(args.next(), SCALED_RANGE, help)?;
            if scaled < fracminhash::MIN_SCALED {
                return Err(format!("{SCALED_RANGE}\n{help}"));
            }
            options.scaled = Some(scaled);
        } else if option && (arg == "--anchor-k" || arg == "--probe-k") {
            let (slot, name) = if arg == "--anchor-k" {
                (&mut options.anchor_k, "anchor k-mer length")
            } else {
                (&mut options.probe_k, "probe k-mer length")
            };
            once(slot, name, help)?;
            *slot = Some(number(
                args.next(),
                &format!("{name} must be a positive integer fitting usize"),
                help,
            )?);
        } else if option && (arg == "--anchor-scaled" || arg == "--probe-scaled") {
            let (slot, name) = if arg == "--anchor-scaled" {
                (&mut options.anchor_scaled, "anchor scaled")
            } else {
                (&mut options.probe_scaled, "probe scaled")
            };
            once(slot, name, help)?;
            let scaled: u64 = number(args.next(), SCALED_RANGE, help)?;
            if scaled == 0 {
                return Err(format!("{SCALED_RANGE}\n{help}"));
            }
            *slot = Some(scaled);
        } else if option && !sketch && arg == "--max-gap" {
            once(&options.max_gap, "max gap", help)?;
            options.max_gap = Some(number(
                args.next(),
                "max gap must be an unsigned 32-bit integer",
                help,
            )?);
        } else if option && !sketch && arg == "--min-island-ani" {
            once(&options.min_island_ani, "min island ANI", help)?;
            let value: f64 = number(
                args.next(),
                "min island ANI must be a number from 0 to 1",
                help,
            )?;
            if !(0.0..=1.0).contains(&value) {
                return Err(format!(
                    "min island ANI must be a number from 0 to 1\n{help}"
                ));
            }
            options.min_island_ani = Some(value);
        } else if option && !sketch && arg == "--symmetric-probes" {
            if std::mem::replace(&mut options.symmetric_probes, true) {
                return Err(format!(
                    "--symmetric-probes may only be specified once\n{help}"
                ));
            }
        } else if option && !sketch && arg == "--bidirectional" {
            if std::mem::replace(&mut options.bidirectional, true) {
                return Err(format!(
                    "--bidirectional may only be specified once\n{help}"
                ));
            }
        } else if option && arg == "--seed" {
            once(&options.seed, "seed", help)?;
            options.seed = Some(number(
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
        // Applicability and ranges are checked here when no sketch input can
        // supply the algorithm; otherwise against the sketch in `compare`.
        if options.algorithm.is_some() {
            options.params()?;
        }
        return Ok(Command::Compare {
            inputs: inputs.try_into().unwrap(),
            options,
        });
    }
    let output = output_path(
        output
            .ok_or_else(|| format!("-o/--output is required\n{help}"))?
            .into_os_string(),
    )?;
    let params = options.params()?;
    Ok(Command::Sketch {
        input: inputs.remove(0),
        output,
        params,
        force,
    })
}

const SCALED_RANGE: &str = "scaled must be an integer from 1 to 18446744073709551615";

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
            "{name} is not implemented yet; implemented: {}\n{help}",
            Algorithm::implemented_names()
        )),
        None => Err(format!(
            "unknown algorithm: {name}; known: {}\n{help}",
            Algorithm::names()
        )),
    }
}
