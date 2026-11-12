use std::ffi::OsString;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use aiips::sketch::aiips::local::LEAST_MAPPING_HITS;
use aiips::sketch::aiips::{self as aiips_sketch, AiipsParams, LocalConfig};
use aiips::sketch::fracminhash::{self, FracMinHashParams};
use aiips::sketch::minhash::{self, MinHashParams, MAX_NUM_HASHES, MIN_NUM_HASHES};
use aiips::sketch::{Algorithm, Sketch};

pub const HELP: &str = "Usage: aiips <COMMAND>\n\nCommands:\n  sketch   Sketch a FASTA/FASTQ file into a .aiips file\n  compare  Estimate similarity between two .aiips sketches\n\nUse aiips <COMMAND> --help for options. Stdin is not supported.";
const SKETCH_HELP: &str = "Usage: aiips sketch INPUT -o OUTPUT [-a ALGORITHM] [-k LENGTH] [-s NUM] [--scaled S] [--seed SEED] [--probe-k K2] [--probe-scaled C2] [--anchor-k K1] [--anchor-scaled C1] [-f]\n\n  -o, --output PATH         Output path (.aiips appended if missing)\n  -a, --algorithm NAME      Sketch algorithm: minhash, fracminhash or aiips\n                            [default: minhash]\n  -k, --kmer-size LENGTH    minhash, fracminhash: k-mer length, 1 to 32 [default: 21]\n  -s, --num-hashes NUM      minhash: number of hash functions, 2 to 1000000 [default: 1000]\n      --scaled S            fracminhash: scale factor, keeps about 1/S of k-mers,\n                            1 to 18446744073709551615 [default: 1000]\n      --probe-k K2          aiips: probe k-mer length, odd, 1 to 31 [default: 15]\n      --probe-scaled C2     aiips: probe scale factor, keeps about 1/C2 of k-mers,\n                            1 to 18446744073709551615 [default: 125]\n      --anchor-k K1         aiips: anchor k-mer length, odd, 1 to 31 [default: 21]\n      --anchor-scaled C1    aiips: anchor scale factor, keeps about 1/C1 of k-mers,\n                            1 to 18446744073709551615 [default: 125]\n      --seed SEED           unsigned 64-bit hash seed [default: 42]\n  -f, --force               Replace an existing output\n  -h, --help                Print help\n\nReads plain/gzip FASTA or FASTQ; k-mers are canonical (strand-independent).\nReserved, not yet implemented: bottomk, minimizer, open-syncmer, closed-syncmer,\nstrobemer.\nUse -- before positional filenames beginning with -. Stdin is not supported.";
const COMPARE_HELP: &str = "Usage: aiips compare INPUT_A INPUT_B [-a ALGORITHM] [-k LENGTH] [-s NUM] [--scaled S] [--seed SEED] [--probe-k K2] [--probe-scaled C2] [--anchor-k K1] [--anchor-scaled C1] [aiips thresholds]\n\n  -a, --algorithm NAME     Algorithm for FASTA/FASTQ inputs: minhash, fracminhash or\n                           aiips [default: minhash, or the sketch's algorithm]\n  -k, --kmer-size LENGTH   minhash, fracminhash: k-mer length [default: 21, or the\n                           sketch's k]\n  -s, --num-hashes NUM     minhash: number of hash functions [default: 1000]\n      --scaled S           fracminhash: scale factor [default: 1000]\n      --probe-k K2         aiips: probe k-mer length [default: 15]\n      --probe-scaled C2    aiips: probe scale factor [default: 125]\n      --anchor-k K1        aiips: anchor k-mer length [default: 21]\n      --anchor-scaled C1   aiips: anchor scale factor [default: 125]\n      --seed SEED          unsigned 64-bit hash seed [default: 42]\n  -h, --help               Print help\n\naiips thresholds for positional ANI and AF (prototype defaults, not tuned):\n      --min-block-span N        cut a block once its span reaches N bases [default: 10000]\n      --max-anchor-gap N        break a run at a larger gap between shared anchors\n                                [default: 50000]\n      --max-gap-ratio R         break a run when gaps differ by more than R times\n                                [default: 2]\n      --min-mapping-hits N      mapping-probe hits to accept a block, at least 2\n                                [default: 3]\n      --min-order-fraction F    share of hits in collinear order, 0 to 1 [default: 0.8]\n      --min-scoring-probes N    pooled scoring probes needed for ANI [default: 50]\n\nEach input is a .aiips sketch or a plain/gzip FASTA/FASTQ file; mixed pairs are\nsupported. Sequence inputs are sketched in memory, using the other input's\nalgorithm and parameters when it is a sketch (given options must then match),\notherwise the options above. Sketches must share one algorithm and identical\nparameters. Prints Jaccard, containment and ANI estimates; AF is -1: it is not\nestimable from these sketches. For aiips, ANI and AF come from shared anchors and\nvalidated blocks and are -1 when unsupported; a table of counts and statuses follows.\nUse -- before positional filenames beginning with -. Stdin is not supported.";

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

    /// k-mer length; the probe k for aiips.
    pub fn k(&self) -> NonZeroUsize {
        match self {
            SketchParams::MinHash(p) => p.k,
            SketchParams::FracMinHash(p) => p.k,
            SketchParams::Aiips(p) => p.probe_k,
        }
    }

    pub fn algorithm(&self) -> Algorithm {
        match self {
            SketchParams::MinHash(_) => Algorithm::MinHash,
            SketchParams::FracMinHash(_) => Algorithm::FracMinHash,
            SketchParams::Aiips(_) => Algorithm::Aiips,
        }
    }
}

/// Algorithm options as given on the command line; `None` means not given.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AlgorithmOptions {
    pub algorithm: Option<Algorithm>,
    pub k: Option<NonZeroUsize>,
    pub num_hashes: Option<usize>,
    pub scaled: Option<u64>,
    pub probe_k: Option<NonZeroUsize>,
    pub probe_scaled: Option<u64>,
    pub anchor_k: Option<NonZeroUsize>,
    pub anchor_scaled: Option<u64>,
    pub seed: Option<u64>,
}

impl AlgorithmOptions {
    /// Parameters from the given options and defaults. Omitting -a keeps
    /// pre-algorithm commands working; minhash is the default.
    pub fn params(&self) -> Result<SketchParams, String> {
        self.params_for(self.algorithm.unwrap_or(Algorithm::MinHash))
    }

    /// Refuses options that do not apply to `algorithm`, then fills defaults.
    fn params_for(&self, algorithm: Algorithm) -> Result<SketchParams, String> {
        self.check_applicable(algorithm)?;
        let k = self.k;
        let seed = self.seed;
        match algorithm {
            Algorithm::MinHash => Ok(SketchParams::MinHash(MinHashParams::new(
                k.unwrap_or(NonZeroUsize::new(minhash::DEFAULT_K).unwrap()),
                self.num_hashes.unwrap_or(minhash::DEFAULT_NUM_HASHES),
                seed.unwrap_or(minhash::DEFAULT_SEED),
            )?)),
            Algorithm::FracMinHash => Ok(SketchParams::FracMinHash(FracMinHashParams::new(
                k.unwrap_or(NonZeroUsize::new(fracminhash::DEFAULT_K).unwrap()),
                self.scaled.unwrap_or(fracminhash::DEFAULT_SCALED),
                seed.unwrap_or(fracminhash::DEFAULT_SEED),
            )?)),
            Algorithm::Aiips => Ok(SketchParams::Aiips(AiipsParams::new(
                self.anchor_k
                    .unwrap_or(NonZeroUsize::new(aiips_sketch::DEFAULT_ANCHOR_K).unwrap()),
                self.anchor_scaled
                    .unwrap_or(aiips_sketch::DEFAULT_ANCHOR_SCALED),
                self.probe_k
                    .unwrap_or(NonZeroUsize::new(aiips_sketch::DEFAULT_PROBE_K).unwrap()),
                self.probe_scaled
                    .unwrap_or(aiips_sketch::DEFAULT_PROBE_SCALED),
                seed.unwrap_or(aiips_sketch::DEFAULT_SEED),
            )?)),
            other => unreachable!("{} rejected while parsing", other.name()),
        }
    }

    /// Refuses given options that `algorithm` does not take.
    fn check_applicable(&self, algorithm: Algorithm) -> Result<(), String> {
        let name = algorithm.name();
        let refuse = |given: bool, message: String| if given { Err(message) } else { Ok(()) };
        let probe_options = || -> Result<(), String> {
            for (given, flag) in [
                (self.probe_k.is_some(), "--probe-k"),
                (self.probe_scaled.is_some(), "--probe-scaled"),
                (self.anchor_k.is_some(), "--anchor-k"),
                (self.anchor_scaled.is_some(), "--anchor-scaled"),
            ] {
                refuse(given, format!("{flag} applies only to aiips, not {name}"))?;
            }
            Ok(())
        };
        match algorithm {
            Algorithm::MinHash => {
                refuse(
                    self.scaled.is_some(),
                    "--scaled applies only to fracminhash, not minhash".into(),
                )?;
                probe_options()
            }
            Algorithm::FracMinHash => {
                refuse(
                    self.num_hashes.is_some(),
                    "-s/--num-hashes applies only to minhash, not fracminhash".into(),
                )?;
                probe_options()
            }
            Algorithm::Aiips => {
                refuse(
                    self.k.is_some(),
                    "-k/--kmer-size does not apply to aiips; use --probe-k".into(),
                )?;
                refuse(
                    self.num_hashes.is_some(),
                    "-s/--num-hashes applies only to minhash, not aiips".into(),
                )?;
                refuse(
                    self.scaled.is_some(),
                    "--scaled applies only to fracminhash, not aiips; use --probe-scaled".into(),
                )
            }
            _ => Ok(()),
        }
    }

    /// The parameters of `template`, a sketch read from `source`, after
    /// checking that every given option applies to it and matches it.
    pub fn matching(&self, template: &SketchParams, source: &Path) -> Result<SketchParams, String> {
        let recorded = format!("recorded in {}", source.display());
        let algorithm = template.algorithm();
        if let Some(given) = self.algorithm.filter(|a| *a != algorithm) {
            return Err(format!(
                "-a {} does not match algorithm={} {}",
                given.name(),
                algorithm.name(),
                recorded
            ));
        }
        // Given options must apply to the sketch's algorithm and equal its
        // recorded values; those are valid, so no range check is needed (and
        // combining given options with defaults could break k1 > k2).
        self.check_applicable(algorithm)?;
        let check = |flag: &str, given: Option<String>, name: &str, value: String| match given
            .filter(|given| *given != value)
        {
            Some(given) => Err(format!(
                "{flag} {given} does not match {name}={value} {recorded}"
            )),
            None => Ok(()),
        };
        check(
            "-k",
            self.k.map(|k| k.to_string()),
            "k",
            template.k().to_string(),
        )?;
        match template {
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
            SketchParams::Aiips(p) => {
                check(
                    "--anchor-k",
                    self.anchor_k.map(|k| k.to_string()),
                    "anchor_k",
                    p.anchor_k.to_string(),
                )?;
                check(
                    "--anchor-scaled",
                    self.anchor_scaled.map(|s| s.to_string()),
                    "anchor_scaled",
                    p.anchor_scaled.to_string(),
                )?;
                check(
                    "--probe-k",
                    self.probe_k.map(|k| k.to_string()),
                    "probe_k",
                    p.probe_k.to_string(),
                )?;
                check(
                    "--probe-scaled",
                    self.probe_scaled.map(|s| s.to_string()),
                    "probe_scaled",
                    p.probe_scaled.to_string(),
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

/// aiips comparison thresholds as given on the command line.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LocalOptions {
    pub min_block_span: Option<u32>,
    pub max_anchor_gap: Option<u32>,
    pub max_gap_ratio: Option<f64>,
    pub min_mapping_hits: Option<u32>,
    pub min_order_fraction: Option<f64>,
    pub min_scoring_probes: Option<u64>,
}

impl LocalOptions {
    /// The first given flag, for refusing them with other algorithms.
    pub fn first_given(&self) -> Option<&'static str> {
        [
            (self.min_block_span.is_some(), "--min-block-span"),
            (self.max_anchor_gap.is_some(), "--max-anchor-gap"),
            (self.max_gap_ratio.is_some(), "--max-gap-ratio"),
            (self.min_mapping_hits.is_some(), "--min-mapping-hits"),
            (self.min_order_fraction.is_some(), "--min-order-fraction"),
            (self.min_scoring_probes.is_some(), "--min-scoring-probes"),
        ]
        .into_iter()
        .find_map(|(given, flag)| given.then_some(flag))
    }

    /// Refuses given thresholds unless `algorithm` is aiips.
    pub fn check_applicable(&self, algorithm: Algorithm) -> Result<(), String> {
        match self.first_given() {
            Some(flag) if algorithm != Algorithm::Aiips => Err(format!(
                "{flag} applies only to aiips, not {}",
                algorithm.name()
            )),
            _ => Ok(()),
        }
    }

    /// The given thresholds over the defaults (ranges checked while parsing).
    pub fn config(&self) -> LocalConfig {
        let d = LocalConfig::default();
        LocalConfig {
            min_block_span: self.min_block_span.unwrap_or(d.min_block_span),
            max_anchor_gap: self.max_anchor_gap.unwrap_or(d.max_anchor_gap),
            max_gap_ratio: self.max_gap_ratio.unwrap_or(d.max_gap_ratio),
            min_mapping_hits: self.min_mapping_hits.unwrap_or(d.min_mapping_hits),
            min_order_fraction: self.min_order_fraction.unwrap_or(d.min_order_fraction),
            min_scoring_probes: self.min_scoring_probes.unwrap_or(d.min_scoring_probes),
        }
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
        /// aiips comparison thresholds.
        local: LocalOptions,
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
    let mut local = LocalOptions::default();
    let mut force = false;
    let mut positional_only = false;
    while let Some(arg) = args.next() {
        let option = !positional_only && arg != "-" && arg.to_string_lossy().starts_with('-');
        let sketch_option = option && sketch;
        let compare_option = option && !sketch;
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
        } else if option && arg == "--probe-k" {
            once(&options.probe_k, "probe k-mer length", help)?;
            options.probe_k = Some(number(args.next(), PROBE_K_RANGE, help)?);
        } else if option && arg == "--probe-scaled" {
            once(&options.probe_scaled, "probe scaled", help)?;
            let scaled: u64 = number(args.next(), PROBE_SCALED_RANGE, help)?;
            if scaled < aiips_sketch::MIN_SCALED {
                return Err(format!("{PROBE_SCALED_RANGE}\n{help}"));
            }
            options.probe_scaled = Some(scaled);
        } else if option && arg == "--anchor-k" {
            once(&options.anchor_k, "anchor k-mer length", help)?;
            options.anchor_k = Some(number(args.next(), ANCHOR_K_RANGE, help)?);
        } else if option && arg == "--anchor-scaled" {
            once(&options.anchor_scaled, "anchor scaled", help)?;
            let scaled: u64 = number(args.next(), ANCHOR_SCALED_RANGE, help)?;
            if scaled < aiips_sketch::MIN_SCALED {
                return Err(format!("{ANCHOR_SCALED_RANGE}\n{help}"));
            }
            options.anchor_scaled = Some(scaled);
        } else if option && arg == "--seed" {
            once(&options.seed, "seed", help)?;
            options.seed = Some(number(
                args.next(),
                "seed must be an unsigned 64-bit integer",
                help,
            )?);
        } else if compare_option && arg == "--min-block-span" {
            once(&local.min_block_span, "min block span", help)?;
            local.min_block_span = Some(at_least(args.next(), 1, BLOCK_SPAN_RANGE, help)?);
        } else if compare_option && arg == "--max-anchor-gap" {
            once(&local.max_anchor_gap, "max anchor gap", help)?;
            local.max_anchor_gap = Some(at_least(args.next(), 1, ANCHOR_GAP_RANGE, help)?);
        } else if compare_option && arg == "--max-gap-ratio" {
            once(&local.max_gap_ratio, "max gap ratio", help)?;
            local.max_gap_ratio = Some(fraction(
                args.next(),
                1.0..=f64::MAX,
                GAP_RATIO_RANGE,
                help,
            )?);
        } else if compare_option && arg == "--min-mapping-hits" {
            once(&local.min_mapping_hits, "min mapping hits", help)?;
            local.min_mapping_hits = Some(at_least(
                args.next(),
                LEAST_MAPPING_HITS,
                MAPPING_HITS_RANGE,
                help,
            )?);
        } else if compare_option && arg == "--min-order-fraction" {
            once(&local.min_order_fraction, "min order fraction", help)?;
            local.min_order_fraction = Some(fraction(
                args.next(),
                0.0..=1.0,
                ORDER_FRACTION_RANGE,
                help,
            )?);
        } else if compare_option && arg == "--min-scoring-probes" {
            once(&local.min_scoring_probes, "min scoring probes", help)?;
            local.min_scoring_probes = Some(at_least(args.next(), 1, SCORING_PROBES_RANGE, help)?);
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
        if let Some(algorithm) = options.algorithm {
            options.params()?;
            local.check_applicable(algorithm)?;
        }
        return Ok(Command::Compare {
            inputs: inputs.try_into().unwrap(),
            options,
            local,
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
const BLOCK_SPAN_RANGE: &str = "min block span must be an integer from 1 to 4294967295";
const ANCHOR_GAP_RANGE: &str = "max anchor gap must be an integer from 1 to 4294967295";
const GAP_RATIO_RANGE: &str = "max gap ratio must be a finite number of at least 1";
const MAPPING_HITS_RANGE: &str = "min mapping hits must be an integer from 2 to 4294967295";
const ORDER_FRACTION_RANGE: &str = "min order fraction must be a number from 0 to 1";
const SCORING_PROBES_RANGE: &str =
    "min scoring probes must be an integer from 1 to 18446744073709551615";
const PROBE_K_RANGE: &str = "probe k-mer length must be an odd integer from 1 to 31";
const PROBE_SCALED_RANGE: &str = "probe scaled must be an integer from 1 to 18446744073709551615";
const ANCHOR_K_RANGE: &str = "anchor k-mer length must be an odd integer from 1 to 31";
const ANCHOR_SCALED_RANGE: &str = "anchor scaled must be an integer from 1 to 18446744073709551615";

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

/// An integer of at least `least`.
fn at_least<T: std::str::FromStr + PartialOrd>(
    value: Option<OsString>,
    least: T,
    message: &str,
    help: &str,
) -> Result<T, String> {
    let value: T = number(value, message, help)?;
    if value < least {
        return Err(format!("{message}\n{help}"));
    }
    Ok(value)
}

/// A finite number within `range`.
fn fraction(
    value: Option<OsString>,
    range: std::ops::RangeInclusive<f64>,
    message: &str,
    help: &str,
) -> Result<f64, String> {
    let value: f64 = number(value, message, help)?;
    if !(value.is_finite() && range.contains(&value)) {
        return Err(format!("{message}\n{help}"));
    }
    Ok(value)
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
