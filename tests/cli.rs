use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};
mod common;
use common::{reserve_directory, TestDirectory};

const STATS: &str = "input\trecords\tbases\tnon_acgt_bases\tk\tvalid_kmers";
const METRICS: &str = "input_a\tinput_b\tjaccard\tcontainment_a_in_b\tcontainment_b_in_a\tani\taf_a_to_b\taf_b_to_a\tani_containment_a_in_b\tani_containment_b_in_a";

struct Workspace(TestDirectory);

/// A valid k = 3 minhash sketch of the sample FASTA, built via the library.
fn fixture_sketch() -> Vec<u8> {
    use aiips::sketch::minhash::{MinHashParams, MinHashSketch};
    let params = MinHashParams::new(std::num::NonZeroUsize::new(3).unwrap(), 8, 42).unwrap();
    MinHashSketch::build(include_bytes!("fixtures/sample.fa").as_slice(), params)
        .unwrap()
        .to_bytes()
}

impl Workspace {
    fn new() -> Self {
        let workspace = Self(TestDirectory::new());
        for (name, bytes) in [
            ("good.fa", include_bytes!("fixtures/sample.fa").as_slice()),
            ("good.fq", include_bytes!("fixtures/sample.fq").as_slice()),
            (
                "good.fa.gz",
                include_bytes!("fixtures/sample.fa.gz").as_slice(),
            ),
            (
                "good.fq.gz",
                include_bytes!("fixtures/sample.fq.gz").as_slice(),
            ),
            ("empty.aiips", b"".as_slice()),
            ("sketch.aiips", fixture_sketch().as_slice()),
            ("other.aiips", b"".as_slice()),
            ("bad.aiips", b"not a placeholder".as_slice()),
            ("bad.fq", b"@bad\nACGT\n+\nIII\n".as_slice()),
            ("empty.fa", b"".as_slice()),
        ] {
            workspace.write(name, bytes);
        }
        workspace
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }

    fn write(&self, name: &str, bytes: &[u8]) {
        fs::write(self.path(name), bytes).unwrap();
    }

    fn run(&self, args: &[&str]) -> Output {
        Command::new(env!("CARGO_BIN_EXE_aiips"))
            .current_dir(&self.0)
            .args(args)
            .output()
            .unwrap()
    }

    fn success(&self, args: &[&str]) -> String {
        let output = self.run(args);
        assert!(
            output.status.success(),
            "{args:?}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(output.stderr.is_empty());
        String::from_utf8(output.stdout).unwrap()
    }

    fn failure(&self, args: &[&str], message: &str) {
        let output = self.run(args);
        assert_eq!(output.status.code(), Some(1), "{args:?}");
        assert!(output.stdout.is_empty(), "{args:?}: unexpected stdout");
        let error = String::from_utf8(output.stderr).unwrap();
        assert!(
            error.contains(message),
            "{args:?}: expected {message:?} in {error:?}"
        );
        assert!(error.starts_with("aiips: "));
    }

    /// Asserts `name` is a regular file holding a default-parameter minhash
    /// sketch with the expected number of hash values.
    fn assert_sketch(&self, name: &str) {
        let metadata = fs::metadata(self.path(name)).unwrap();
        assert!(metadata.is_file());
        let text = fs::read_to_string(self.path(name)).unwrap();
        assert!(
            text.starts_with("aiips-sketch v1 algorithm=minhash "),
            "{text:.80}"
        );
        let num_hashes: usize = text
            .split_whitespace()
            .find_map(|field| field.strip_prefix("num_hashes="))
            .unwrap()
            .parse()
            .unwrap();
        assert_eq!(text.lines().count(), 2 + num_hashes);
    }

    /// Builds a sketch of `input` at `name` with extra sketch arguments.
    fn sketch(&self, input: &str, name: &str, extra: &[&str]) {
        let mut args = vec!["sketch", "-a", "minhash", input, "-o", name];
        args.extend(extra);
        self.success(&args);
    }

    fn assert_no_staging_files(&self) {
        assert!(fs::read_dir(&self.0).unwrap().all(|entry| !entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .starts_with(".aiips-")));
    }
}

macro_rules! reject {
    ($name:ident, [$($arg:expr),*], $message:expr) => {
        #[test]
        fn $name() { Workspace::new().failure(&[$($arg),*], $message); }
    };
}

macro_rules! sketch_format {
    ($name:ident, $input:expr) => {
        #[test]
        fn $name() {
            let t = Workspace::new();
            let stdout = t.success(&["sketch", "-a", "minhash", $input, "-o", "out", "-k", "3"]);
            assert_eq!(
                stdout,
                format!("{STATS}\n{}\t2\t13\t2\t3\t7\noutput\tout.aiips\n", $input)
            );
            t.assert_sketch("out.aiips");
        }
    };
}

sketch_format!(sketch_fasta_known_counts, "good.fa");
sketch_format!(sketch_fastq_known_counts, "good.fq");
sketch_format!(sketch_gzip_fasta, "good.fa.gz");
sketch_format!(sketch_gzip_fastq, "good.fq.gz");

#[test]
fn sketch_default_parameters_and_exact_suffix() {
    let t = Workspace::new();
    t.write("long.fa", b">a\nGATTACAGGCTTAACGTTAGCCATGCAAGTCTAGGCATCG\n");
    let stdout = t.success(&["sketch", "-a", "minhash", "long.fa", "-o", "out.aiips"]);
    assert!(stdout.contains("long.fa\t1\t40\t0\t21\t20\n"), "{stdout}");
    t.assert_sketch("out.aiips");
    let text = fs::read_to_string(t.path("out.aiips")).unwrap();
    assert!(text.starts_with(
        "aiips-sketch v1 algorithm=minhash k=21 num_hashes=1000 seed=42 hash=murmur3-fmix64 canonical=true\nrecords=1 bases=40 non_acgt_bases=0 valid_kmers=20\n"
    ));
    assert!(!t.path("out.aiips.aiips").exists());
}

#[test]
fn sketch_refuses_input_without_valid_kmers_and_creates_nothing() {
    let t = Workspace::new();
    // good.fa has records of 13 bases at most: no 21-mers at the default k.
    for args in [
        vec!["sketch", "-a", "minhash", "good.fa", "-o", "out"],
        vec!["sketch", "-a", "minhash", "good.fa", "-o", "out", "-f"],
    ] {
        t.failure(&args, "no valid k-mers of length 21; nothing to sketch");
        assert!(!t.path("out.aiips").exists());
        t.assert_no_staging_files();
    }
    t.write("ambiguous.fa", b">n\nNNNNNNNNNN\n");
    t.failure(
        &[
            "sketch",
            "-a",
            "minhash",
            "-k",
            "3",
            "ambiguous.fa",
            "-o",
            "out",
        ],
        "no valid k-mers",
    );
}

#[test]
fn sketch_without_algorithm_defaults_to_minhash() {
    // Pre-algorithm commands keep working and produce the same sketch.
    let t = Workspace::new();
    let implicit = t.success(&["sketch", "-k", "3", "good.fa", "-o", "implicit"]);
    let explicit = t.success(&[
        "sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "explicit",
    ]);
    assert_eq!(
        implicit.replace("implicit", "X"),
        explicit.replace("explicit", "X")
    );
    assert_eq!(
        fs::read(t.path("implicit.aiips")).unwrap(),
        fs::read(t.path("explicit.aiips")).unwrap()
    );
    t.assert_sketch("implicit.aiips");
}

#[test]
fn sketch_rejects_unimplemented_and_unknown_algorithms() {
    let t = Workspace::new();
    for name in [
        "bottomk",
        "minimizer",
        "open-syncmer",
        "closed-syncmer",
        "strobemer",
    ] {
        t.failure(
            &["sketch", "-a", name, "-k", "3", "good.fa", "-o", "out"],
            &format!("{name} is not implemented yet; implemented: minhash, fracminhash"),
        );
    }
    for name in ["MinHash", "sha1", ""] {
        t.failure(
            &["sketch", "--algorithm", name, "good.fa", "-o", "out"],
            "unknown algorithm",
        );
    }
    t.failure(
        &["sketch", "good.fa", "-o", "out", "-a"],
        "missing algorithm name",
    );
    t.failure(
        &[
            "sketch", "-a", "minhash", "-a", "minhash", "good.fa", "-o", "out",
        ],
        "algorithm may only be specified once",
    );
    assert!(!t.path("out.aiips").exists());
}

#[test]
fn sketch_validates_minhash_parameters() {
    let t = Workspace::new();
    let base = ["sketch", "-a", "minhash", "good.fa", "-o", "out"];
    for (extra, message) in [
        (vec!["-k", "33"], "1 to 32"),
        (vec!["-s", "1"], "from 2 to 1000000"),
        (vec!["-s", "1000001"], "from 2 to 1000000"),
        (vec!["-s", "0"], "from 2 to 1000000"),
        (
            vec!["-s", "x"],
            "number of hash functions must be an integer",
        ),
        (
            vec!["--num-hashes"],
            "number of hash functions must be an integer",
        ),
        (
            vec!["--seed", "-1"],
            "seed must be an unsigned 64-bit integer",
        ),
        (
            vec!["--seed", "18446744073709551616"],
            "seed must be an unsigned",
        ),
        (
            vec!["-s", "8", "--num-hashes", "8"],
            "may only be specified once",
        ),
        (
            vec!["--seed", "1", "--seed", "1"],
            "seed may only be specified once",
        ),
    ] {
        let mut args = base.to_vec();
        args.extend(extra);
        t.failure(&args, message);
    }
    assert!(!t.path("out.aiips").exists());
}

#[test]
fn sketch_records_custom_parameters_and_is_deterministic() {
    let t = Workspace::new();
    let extra = ["-k", "3", "-s", "5", "--seed", "18446744073709551615"];
    t.sketch("good.fa", "a", &extra);
    t.sketch("good.fa", "b", &extra);
    let text = fs::read_to_string(t.path("a.aiips")).unwrap();
    assert!(text.starts_with(
        "aiips-sketch v1 algorithm=minhash k=3 num_hashes=5 seed=18446744073709551615 hash=murmur3-fmix64 canonical=true\nrecords=2 bases=13 non_acgt_bases=2 valid_kmers=7\n"
    ));
    assert_eq!(text.lines().count(), 7);
    assert_eq!(text, fs::read_to_string(t.path("b.aiips")).unwrap());
    // FASTA, FASTQ and gzip inputs with identical sequences give identical sketches.
    for (input, name) in [("good.fq", "c"), ("good.fa.gz", "d"), ("good.fq.gz", "e")] {
        t.sketch(input, name, &extra);
        assert_eq!(
            text,
            fs::read_to_string(t.path(&format!("{name}.aiips"))).unwrap()
        );
    }
}

#[test]
fn sketch_long_options_and_uppercase_suffix() {
    let t = Workspace::new();
    let stdout = t.success(&[
        "sketch",
        "--algorithm",
        "minhash",
        "--kmer-size",
        "3",
        "--output",
        "out.AIIPS",
        "good.fa",
    ]);
    assert!(stdout.ends_with("output\tout.AIIPS\n"));
    t.assert_sketch("out.AIIPS");
}

#[test]
fn sketch_appends_to_other_extension() {
    let t = Workspace::new();
    t.success(&[
        "sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "out.txt",
    ]);
    t.assert_sketch("out.txt.aiips");
}

#[test]
fn sketch_ambiguity_and_record_boundaries() {
    let t = Workspace::new();
    t.write("runs.fa", b">a\nAC\nGTNacg\n>b\nTT\n>c\nAA\n");
    let stdout = t.success(&["sketch", "-a", "minhash", "runs.fa", "-k", "3", "-o", "out"]);
    assert!(stdout.contains("runs.fa\t3\t12\t1\t3\t3\n"));
}

#[test]
fn sketch_existing_output_refused_and_preserved() {
    let t = Workspace::new();
    t.write("out.aiips", b"existing result");
    t.failure(
        &["sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "out"],
        "use -f or --force",
    );
    assert_eq!(fs::read(t.path("out.aiips")).unwrap(), b"existing result");
}

/// Spells a `/`-separated relative path with the platform separator, as the CLI
/// does when it appends `.aiips` to a path with a parent directory.
fn native(path: &str) -> String {
    path.split('/').collect::<PathBuf>().display().to_string()
}

#[cfg(unix)]
#[test]
fn sketch_force_replaces_and_cleans_staging() {
    for flag in ["-f", "--force"] {
        let t = Workspace::new();
        t.write("out.aiips", b"existing result");
        t.success(&[
            "sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "out", flag,
        ]);
        t.assert_sketch("out.aiips");
        t.assert_no_staging_files();
    }
}

// Without Unix device/inode identity checks, forced replacement is refused.
#[cfg(not(unix))]
#[test]
fn sketch_force_refuses_existing_output_without_unix_identity_checks() {
    for flag in ["-f", "--force"] {
        let t = Workspace::new();
        t.write("out.aiips", b"existing result");
        t.failure(
            &[
                "sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "out", flag,
            ],
            "forced replacement requires Unix file-identity safeguards",
        );
        assert_eq!(fs::read(t.path("out.aiips")).unwrap(), b"existing result");
        t.assert_no_staging_files();
    }
}

#[test]
fn sketch_force_creates_new_output() {
    let t = Workspace::new();
    t.success(&[
        "sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "out", "--force",
    ]);
    t.assert_sketch("out.aiips");
    t.assert_no_staging_files();
}

#[test]
fn sketch_force_preserves_output_on_parse_failure() {
    let t = Workspace::new();
    t.write("out.aiips", b"existing result");
    // Off Unix, the forced-replacement refusal precedes input parsing.
    let expected = if cfg!(unix) {
        "bad.fq"
    } else {
        "forced replacement requires Unix"
    };
    t.failure(
        &[
            "sketch", "-a", "minhash", "-k", "3", "bad.fq", "-o", "out", "--force",
        ],
        expected,
    );
    assert_eq!(fs::read(t.path("out.aiips")).unwrap(), b"existing result");
    t.assert_no_staging_files();
}

#[test]
fn sketch_failure_does_not_create_output() {
    let t = Workspace::new();
    t.failure(
        &["sketch", "-a", "minhash", "-k", "3", "bad.fq", "-o", "out"],
        "bad.fq",
    );
    assert!(!t.path("out.aiips").exists());
}

#[test]
fn sketch_missing_destination_directory() {
    let t = Workspace::new();
    for force in [vec![], vec!["--force"]] {
        let mut args = vec![
            "sketch",
            "-a",
            "minhash",
            "-k",
            "3",
            "good.fa",
            "-o",
            "missing/out",
        ];
        args.extend(force);
        t.failure(&args, &native("missing/out.aiips"));
    }
    assert!(!t.path("missing").exists());
}

#[test]
fn sketch_directory_destination_preserved_and_staging_cleaned() {
    let t = Workspace::new();
    fs::create_dir(t.path("out.aiips")).unwrap();
    t.write("out.aiips/result", b"keep");
    t.failure(
        &[
            "sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "out", "--force",
        ],
        "out.aiips",
    );
    assert_eq!(fs::read(t.path("out.aiips/result")).unwrap(), b"keep");
    t.assert_no_staging_files();
}

#[test]
fn sketch_rejects_identical_input_output() {
    let t = Workspace::new();
    let data = b">a\nACGT\n";
    t.write("input.aiips", data);
    t.failure(
        &[
            "sketch",
            "-a",
            "minhash",
            "-k",
            "3",
            "input.aiips",
            "-o",
            "./input.aiips",
            "--force",
        ],
        "output resolves to the input",
    );
    assert_eq!(fs::read(t.path("input.aiips")).unwrap(), data);
}

#[cfg(unix)]
#[test]
fn sketch_rejects_symlink_alias() {
    let t = Workspace::new();
    std::os::unix::fs::symlink(t.path("good.fa"), t.path("alias.aiips")).unwrap();
    t.failure(
        &[
            "sketch",
            "-a",
            "minhash",
            "-k",
            "3",
            "good.fa",
            "-o",
            "alias.aiips",
            "-f",
        ],
        "output resolves to the input",
    );
    assert_eq!(
        fs::read(t.path("good.fa")).unwrap(),
        include_bytes!("fixtures/sample.fa")
    );
}

#[cfg(unix)]
#[test]
fn sketch_dangling_symlink_requires_force() {
    let t = Workspace::new();
    std::os::unix::fs::symlink(t.path("missing"), t.path("out.aiips")).unwrap();
    t.failure(
        &["sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "out"],
        "--force",
    );
    t.success(&[
        "sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "out", "-f",
    ]);
    assert!(!fs::symlink_metadata(t.path("out.aiips"))
        .unwrap()
        .file_type()
        .is_symlink());
    assert!(!t.path("missing").exists());
}

#[test]
fn sketch_double_dash_input_beginning_with_dash() {
    let t = Workspace::new();
    t.write("-input.fa", b">a\nACGT\n");
    let stdout = t.success(&[
        "sketch",
        "-a",
        "minhash",
        "-k",
        "3",
        "-o",
        "out",
        "--",
        "-input.fa",
    ]);
    assert!(stdout.contains("\n-input.fa\t1\t4\t0\t3\t2\n"));
    assert!(stdout.ends_with("output\tout.aiips\n"));
    t.assert_sketch("out.aiips");
}

#[test]
fn sketch_output_label_with_parent_uses_native_separator() {
    let t = Workspace::new();
    fs::create_dir(t.path("sub")).unwrap();
    let stdout = t.success(&[
        "sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "sub/out",
    ]);
    // Labels escape backslashes, so a Windows separator prints as `\\`.
    let label = native("sub/out.aiips").replace('\\', "\\\\");
    assert!(stdout.ends_with(&format!("output\t{label}\n")), "{stdout}");
    t.assert_sketch("sub/out.aiips");
    // An already-suffixed output keeps the spelling exactly as given.
    let stdout = t.success(&[
        "sketch",
        "-a",
        "minhash",
        "-k",
        "3",
        "good.fa",
        "-o",
        "sub/kept.aiips",
    ]);
    assert!(stdout.ends_with("output\tsub/kept.aiips\n"), "{stdout}");
}

// Control characters are invalid in Windows file names, and `\` is a separator.
#[cfg(unix)]
#[test]
fn sketch_double_dash_and_escaped_paths() {
    let t = Workspace::new();
    t.write("-input\t.fa", b">a\nACGT\n");
    let stdout = t.success(&[
        "sketch",
        "-a",
        "minhash",
        "-k",
        "3",
        "-o",
        "out\nname",
        "--",
        "-input\t.fa",
    ]);
    assert!(stdout.contains("-input\\t.fa\t1\t4"));
    assert!(stdout.ends_with("output\tout\\nname.aiips\n"));
    t.assert_sketch("out\nname.aiips");
}

#[test]
fn sketch_help() {
    for flag in ["-h", "--help"] {
        assert!(Workspace::new()
            .success(&["sketch", "-a", "minhash", "-k", "3", flag])
            .contains("Usage: aiips sketch"));
    }
}

reject!(
    sketch_missing_input,
    ["sketch", "-a", "minhash", "-k", "3", "-o", "out"],
    "exactly 1"
);
reject!(
    sketch_extra_input,
    ["sketch", "-a", "minhash", "-k", "3", "good.fa", "good.fq", "-o", "out"],
    "exactly 1"
);
reject!(
    sketch_missing_output,
    ["sketch", "-a", "minhash", "-k", "3", "good.fa"],
    "--output is required"
);
reject!(
    sketch_missing_output_value,
    ["sketch", "-a", "minhash", "-k", "3", "good.fa", "-o"],
    "missing output path"
);
reject!(
    sketch_duplicate_output,
    ["sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "out", "-o", "other"],
    "only be specified once"
);
reject!(
    sketch_missing_file,
    [
        "sketch",
        "-a",
        "minhash",
        "-k",
        "3",
        "missing.fa",
        "-o",
        "out"
    ],
    "missing.fa"
);
reject!(
    sketch_malformed_fastq,
    ["sketch", "-a", "minhash", "-k", "3", "bad.fq", "-o", "out"],
    "bad.fq"
);
reject!(
    sketch_empty_input,
    ["sketch", "-a", "minhash", "-k", "3", "empty.fa", "-o", "out"],
    "empty.fa"
);
reject!(
    sketch_zero_k,
    ["sketch", "-a", "minhash", "good.fa", "-o", "out", "-k", "0"],
    "positive integer"
);
reject!(
    sketch_negative_k,
    ["sketch", "-a", "minhash", "good.fa", "-o", "out", "-k", "-1"],
    "positive integer"
);
reject!(
    sketch_nonnumeric_k,
    ["sketch", "-a", "minhash", "good.fa", "-o", "out", "-k", "abc"],
    "positive integer"
);
reject!(
    sketch_missing_k,
    ["sketch", "-a", "minhash", "good.fa", "-o", "out", "-k"],
    "positive integer"
);
reject!(
    sketch_overflow_k,
    [
        "sketch",
        "good.fa",
        "-o",
        "out",
        "-k",
        "999999999999999999999999999999"
    ],
    "positive integer"
);
reject!(
    sketch_unknown_option,
    [
        "sketch",
        "-a",
        "minhash",
        "-k",
        "3",
        "good.fa",
        "-o",
        "out",
        "--unknown"
    ],
    "unknown option"
);
reject!(
    sketch_stdin,
    ["sketch", "-a", "minhash", "-k", "3", "-", "-o", "out"],
    "stdin"
);
reject!(
    sketch_stdout_destination,
    ["sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "-"],
    "filesystem path"
);
reject!(
    sketch_directory_input,
    ["sketch", "-a", "minhash", "-k", "3", ".", "-o", "out"],
    "regular file"
);

/// Two sketches of identical sequences (FASTA vs FASTQ) at k = 3.
fn identical_pair(t: &Workspace) {
    t.sketch("good.fa", "a", &["-k", "3"]);
    t.sketch("good.fq", "b", &["-k", "3"]);
}

#[test]
fn compare_identical_sequences_reports_exact_ones() {
    let t = Workspace::new();
    identical_pair(&t);
    assert_eq!(
        t.success(&["compare", "a.aiips", "b.aiips"]),
        format!(
            "{STATS}\na.aiips\t2\t13\t2\t3\t7\nb.aiips\t2\t13\t2\t3\t7\n\n{METRICS}\na.aiips\tb.aiips\t1.000000\t1.000000\t1.000000\t1.000000\t-1\t-1\t1.000000\t1.000000\n"
        )
    );
}

#[test]
fn compare_disjoint_kmer_sets_reports_zeros() {
    let t = Workspace::new();
    t.write("a.fa", b">a\nAAAAAAAA\n");
    t.write("c.fa", b">c\nCCCCCCCC\n");
    t.sketch("a.fa", "a", &["-k", "3"]);
    t.sketch("c.fa", "c", &["-k", "3"]);
    let stdout = t.success(&["compare", "a.aiips", "c.aiips"]);
    assert!(stdout.ends_with(
        "\na.aiips\tc.aiips\t0.000000\t0.000000\t0.000000\t0.000000\t-1\t-1\t0.000000\t0.000000\n"
    ));
}

#[test]
fn compare_is_strand_independent() {
    let t = Workspace::new();
    t.write("f.fa", b">f\nGATTACAGGCTTAACGTTAG\n");
    t.write("r.fa", b">r\nCTAACGTTAAGCCTGTAATC\n");
    t.sketch("f.fa", "f", &["-k", "7"]);
    t.sketch("r.fa", "r", &["-k", "7"]);
    assert!(t
        .success(&["compare", "f.aiips", "r.aiips"])
        .contains("\t1.000000\t1.000000\t1.000000\t1.000000\t-1\t-1\t1.000000\t1.000000\n"));
}

/// Metrics row expected from the library on the same sketch files, checking
/// the CLI's parsing, guardrails and formatting against an independent path.
fn library_row(t: &Workspace, a: &str, b: &str) -> String {
    use aiips::sketch::{minhash, parse, Sketch};
    let load = |name: &str| match parse(&fs::read(t.path(name)).unwrap()).unwrap() {
        Sketch::MinHash(sketch) => sketch,
        other => panic!("{other:?}"),
    };
    let e = minhash::estimate(&load(a), &load(b)).unwrap();
    format!(
        "{a}\t{b}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t-1\t-1\t{:.6}\t{:.6}\n",
        e.jaccard,
        e.containment_a_in_b,
        e.containment_b_in_a,
        e.ani,
        e.ani_containment_a_in_b,
        e.ani_containment_b_in_a
    )
}

#[test]
fn compare_partial_overlap_matches_library_and_swaps_directions() {
    let t = Workspace::new();
    // B contains A's first 60 bases plus 60 unrelated bases.
    let shared = "GATTACAGGCTTAACGTTAGCCATGCAAGTCTAGGCATCGTTGACCTAGGATCCATGGCA";
    let extra = "TTGCAACGGTACCGATTGCATGCATCCGGAATTCGCGTATAGCTAGCATCGATCGGCTAA";
    t.write("a.fa", format!(">a\n{shared}\n").as_bytes());
    t.write("b.fa", format!(">b\n{shared}{extra}\n").as_bytes());
    t.sketch("a.fa", "a", &["-k", "11", "-s", "500"]);
    t.sketch("b.fa", "b", &["-k", "11", "-s", "500"]);
    let forward = t.success(&["compare", "a.aiips", "b.aiips"]);
    assert!(
        forward.ends_with(&library_row(&t, "a.aiips", "b.aiips")),
        "{forward}"
    );
    let backward = t.success(&["compare", "b.aiips", "a.aiips"]);
    assert!(backward.ends_with(&library_row(&t, "b.aiips", "a.aiips")));
    let fields = |out: &str| -> Vec<String> {
        out.lines()
            .last()
            .unwrap()
            .split('\t')
            .map(String::from)
            .collect()
    };
    let (f, b) = (fields(&forward), fields(&backward));
    assert_eq!(f[2], b[2], "Jaccard is symmetric");
    assert_eq!((&f[3], &f[4]), (&b[4], &b[3]), "containments swap");
    assert_eq!((&f[8], &f[9]), (&b[9], &b[8]), "containment ANIs swap");
    // A is (almost) contained in B; B is about half contained in A.
    let c_ab: f64 = f[3].parse().unwrap();
    let c_ba: f64 = f[4].parse().unwrap();
    assert!(c_ab > c_ba, "{c_ab} vs {c_ba}");
}

#[test]
fn compare_detects_sketches_by_content_not_extension() {
    let t = Workspace::new();
    identical_pair(&t);
    fs::copy(t.path("a.aiips"), t.path("renamed.txt")).unwrap();
    fs::copy(t.path("b.aiips"), t.path("upper.AiIpS")).unwrap();
    assert!(t
        .success(&["compare", "renamed.txt", "upper.AiIpS"])
        .contains("renamed.txt\tupper.AiIpS\t1.000000"));
}

#[test]
fn compare_aiips_names_must_be_sketches() {
    let t = Workspace::new();
    identical_pair(&t);
    // bad.aiips holds non-sketch text: an .aiips name is never read as FASTA.
    t.failure(
        &["compare", "bad.aiips", "a.aiips"],
        "bad.aiips: not an aiips sketch; run `aiips sketch` first",
    );
}

macro_rules! compare_sequences {
    ($name:ident, $first:expr, $second:expr) => {
        #[test]
        fn $name() {
            // Pre-sketch compare commands keep working: sequence inputs are
            // sketched in memory with minhash (k from -k, default s and seed).
            let t = Workspace::new();
            assert_eq!(
                t.success(&["compare", $first, $second, "-k", "3"]),
                format!(
                    "{STATS}\n{0}\t2\t13\t2\t3\t7\n{1}\t2\t13\t2\t3\t7\n\n{METRICS}\n{0}\t{1}\t1.000000\t1.000000\t1.000000\t1.000000\t-1\t-1\t1.000000\t1.000000\n",
                    $first, $second
                )
            );
        }
    };
}

compare_sequences!(compare_fasta_fastq, "good.fa", "good.fq");
compare_sequences!(compare_gzip_fastx, "good.fa.gz", "good.fq.gz");
compare_sequences!(compare_same_file_twice, "good.fa", "good.fa");

#[test]
fn compare_sequence_inputs_equal_explicit_default_sketches() {
    let t = Workspace::new();
    let shared = "GATTACAGGCTTAACGTTAGCCATGCAAGTCTAGGCATCGTTGACCTAGGATCCATGGCA";
    let extra = "TTGCAACGGTACCGATTGCATGCATCCGGAATTCGCGTATAGCTAGCATCGATCGGCTAA";
    t.write("a.fa", format!(">a\n{shared}\n").as_bytes());
    t.write("b.fa", format!(">b\n{shared}{extra}\n").as_bytes());
    t.sketch("a.fa", "a", &["-k", "11"]);
    t.sketch("b.fa", "b", &["-k", "11"]);
    let metrics = |out: String| {
        out.lines()
            .last()
            .unwrap()
            .split_once("\t")
            .unwrap()
            .1
            .split_once("\t")
            .unwrap()
            .1
            .to_string()
    };
    let from_sketches = metrics(t.success(&["compare", "a.aiips", "b.aiips"]));
    assert_eq!(
        metrics(t.success(&["compare", "a.fa", "b.fa", "-k", "11"])),
        from_sketches
    );
    assert_eq!(
        metrics(t.success(&["compare", "a.aiips", "b.fa"])),
        from_sketches
    );
    assert_eq!(
        metrics(t.success(&["compare", "a.fa", "b.aiips", "-k", "11"])),
        from_sketches
    );
}

#[test]
fn compare_sequence_input_adopts_sketch_parameters() {
    let t = Workspace::new();
    t.sketch("good.fa", "custom", &["-k", "3", "-s", "8", "--seed", "1"]);
    for args in [
        vec!["compare", "custom.aiips", "good.fq"],
        vec!["compare", "good.fq.gz", "custom.aiips"],
        vec!["compare", "custom.aiips", "good.fq", "--kmer-size", "3"],
    ] {
        assert!(t
            .success(&args)
            .contains("\t1.000000\t1.000000\t1.000000\t1.000000\t-1\t-1\t1.000000\t1.000000\n"));
    }
    for flag in ["-k", "--kmer-size"] {
        t.failure(
            &["compare", "custom.aiips", "good.fq", flag, "4"],
            "-k 4 does not match k=3 recorded in custom.aiips",
        );
    }
}

#[test]
fn compare_sequences_detected_by_content_not_extension() {
    let t = Workspace::new();
    identical_pair(&t);
    // A sketch with a non-.aiips name is still read as a sketch.
    fs::copy(t.path("a.aiips"), t.path("sketch.txt")).unwrap();
    assert!(t
        .success(&["compare", "sketch.txt", "good.fq", "-k", "3"])
        .contains("sketch.txt\tgood.fq\t1.000000"));
}

#[test]
fn compare_sequence_inputs_are_validated() {
    let t = Workspace::new();
    // good.fa has no 21-mers at the default k: Jaccard would be undefined.
    t.failure(
        &["compare", "good.fa", "good.fq"],
        "good.fa: no valid k-mers of length 21; nothing to sketch",
    );
    t.failure(&["compare", "good.fa", "good.fq", "-k", "33"], "1 to 32");
    t.failure(&["compare", "good.fa", "bad.fq", "-k", "3"], "bad.fq");
    t.failure(&["compare", "empty.fa", "good.fa", "-k", "3"], "empty.fa");
    t.failure(
        &["compare", "good.fa", "good.fq", "-k", "3", "-k", "3"],
        "k-mer length may only be specified once",
    );
    for bad in ["0", "-1", "abc"] {
        t.failure(
            &["compare", "good.fa", "good.fq", "-k", bad],
            "positive integer",
        );
    }
    t.failure(&["compare", "good.fa", "good.fq", "-k"], "positive integer");
}

#[test]
fn compare_rejects_old_empty_placeholders() {
    let t = Workspace::new();
    identical_pair(&t);
    t.failure(
        &["compare", "empty.aiips", "a.aiips"],
        "empty.aiips: empty file: placeholder from an older aiips version",
    );
    t.failure(
        &["compare", "a.aiips", "empty.aiips"],
        "empty.aiips: empty file",
    );
}

#[test]
fn compare_halts_unless_both_sketches_are_implemented() {
    let t = Workspace::new();
    identical_pair(&t);
    for name in [
        "bottomk",
        "minimizer",
        "open-syncmer",
        "closed-syncmer",
        "strobemer",
    ] {
        t.write(
            "other.aiips",
            format!("aiips-sketch v1 algorithm={name}\n").as_bytes(),
        );
        t.failure(
            &["compare", "a.aiips", "other.aiips"],
            &format!(
                "comparison requires two sketches of one algorithm, minhash or fracminhash (a.aiips: minhash, other.aiips: {name})"
            ),
        );
        t.failure(
            &["compare", "other.aiips", "other.aiips"],
            &format!("(other.aiips: {name}, other.aiips: {name})"),
        );
    }
    // A sequence input is never sketched to match an unimplemented sketch.
    t.write("other.aiips", b"aiips-sketch v1 algorithm=bottomk\n");
    t.failure(
        &["compare", "good.fa", "other.aiips", "-k", "3"],
        "comparison requires two sketches of one algorithm, minhash or fracminhash (good.fa: fasta/fastq, other.aiips: bottomk)",
    );
    t.write("other.aiips", b"aiips-sketch v1 algorithm=hyperloglog\n");
    t.failure(
        &["compare", "a.aiips", "other.aiips"],
        "unknown sketch algorithm: hyperloglog",
    );
}

#[test]
fn compare_halts_on_each_parameter_mismatch() {
    let t = Workspace::new();
    t.sketch("good.fa", "base", &["-k", "3", "-s", "8", "--seed", "1"]);
    for (name, extra, message) in [
        ("k", ["-k", "4", "-s", "8", "--seed", "1"], "k (3 vs 4)"),
        (
            "s",
            ["-k", "3", "-s", "9", "--seed", "1"],
            "num_hashes (8 vs 9)",
        ),
        (
            "seed",
            ["-k", "3", "-s", "8", "--seed", "2"],
            "seed (1 vs 2)",
        ),
    ] {
        t.sketch("good.fa", name, &extra);
        t.failure(
            &["compare", "base.aiips", &format!("{name}.aiips")],
            &format!("sketch parameters differ: {message}; re-sketch with identical parameters"),
        );
    }
    let text = fs::read_to_string(t.path("base.aiips")).unwrap();
    t.write(
        "family.aiips",
        text.replace(
            "hash=murmur3-fmix64 canonical=true",
            "hash=xxh3 canonical=false",
        )
        .as_bytes(),
    );
    t.failure(
        &["compare", "base.aiips", "family.aiips"],
        "hash (murmur3-fmix64 vs xxh3), canonical (true vs false)",
    );
}

#[test]
fn compare_rejects_corrupt_sketches() {
    let t = Workspace::new();
    t.sketch("good.fa", "a", &["-k", "3", "-s", "4"]);
    let text = fs::read_to_string(t.path("a.aiips")).unwrap();
    let truncated: String = text.lines().take(5).map(|l| format!("{l}\n")).collect();
    for (content, message) in [
        (truncated, "expected num_hashes=4 hash values, found 3"),
        (format!("{text}5\n"), "more than num_hashes=4"),
        (text.replace('\n', "\r\n"), "carriage return"),
        (text.trim_end().to_string(), "missing final newline"),
        (text.replace(" v1 ", " v9 "), "unsupported format version"),
        (
            text.replace("seed=42", "seed=42 seed=42"),
            "duplicate header field: seed",
        ),
        (text.replace("k=3", "k=40"), "1 to 32"),
        (
            text.replace("records=2", "records=two"),
            "records is not an unsigned integer",
        ),
    ] {
        t.write("bad.aiips", content.as_bytes());
        t.failure(&["compare", "a.aiips", "bad.aiips"], message);
    }
    t.write("binary.aiips", b"aiips-sketch v1 algorithm=minhash \xff\n");
    t.failure(&["compare", "a.aiips", "binary.aiips"], "not valid UTF-8");
}

#[test]
fn compare_kmer_option_must_match_sketch_inputs() {
    let t = Workspace::new();
    identical_pair(&t);
    assert!(t
        .success(&["compare", "a.aiips", "b.aiips", "-k", "3"])
        .contains("\na.aiips\tb.aiips\t1.000000\t"));
    t.failure(
        &["compare", "a.aiips", "b.aiips", "-k", "21"],
        "-k 21 does not match k=3 recorded in a.aiips",
    );
}

#[test]
fn compare_double_dash_input_beginning_with_dash() {
    let t = Workspace::new();
    identical_pair(&t);
    fs::copy(t.path("a.aiips"), t.path("-a.aiips")).unwrap();
    let stdout = t.success(&["compare", "--", "-a.aiips", "b.aiips"]);
    assert!(
        stdout.contains("\n-a.aiips\tb.aiips\t1.000000\t"),
        "{stdout}"
    );
}

// Control characters are invalid in Windows file names, and `\` is a separator.
#[cfg(unix)]
#[test]
fn compare_double_dash_and_escaped_paths() {
    let t = Workspace::new();
    identical_pair(&t);
    let name = "-input\\\t\r\n.aiips";
    fs::copy(t.path("a.aiips"), t.path(name)).unwrap();
    let stdout = t.success(&["compare", "--", name, "b.aiips"]);
    assert!(stdout.contains("-input\\\\\\t\\r\\n.aiips\t2\t13"));
    assert_eq!(stdout.lines().count(), 6);
    assert_eq!(stdout.lines().last().unwrap().split('\t').count(), 10);
}

#[test]
fn compare_help() {
    for flag in ["-h", "--help"] {
        let stdout = Workspace::new().success(&["compare", flag]);
        assert!(stdout.contains(
            "Usage: aiips compare INPUT_A INPUT_B [-a ALGORITHM] [-k LENGTH] [-s NUM] [--scaled S] [--seed SEED]"
        ));
        assert!(stdout.contains("AF is -1"));
    }
}

#[test]
fn compare_sketch_directory_rejected() {
    let t = Workspace::new();
    fs::create_dir(t.path("dir.aiips")).unwrap();
    t.failure(&["compare", "dir.aiips", "empty.aiips"], "regular file");
}

#[test]
fn compare_generated_sketch_round_trip() {
    let t = Workspace::new();
    t.success(&[
        "sketch",
        "-a",
        "minhash",
        "-k",
        "3",
        "good.fa",
        "-o",
        "generated",
    ]);
    let stdout = t.success(&["compare", "generated.aiips", "generated.aiips"]);
    assert!(stdout.contains("generated.aiips\t2\t13\t2\t3\t7\n"));
    assert!(stdout.contains("generated.aiips\tgenerated.aiips\t1.000000\t"));
}

reject!(compare_missing_inputs, ["compare"], "exactly 2");
reject!(compare_one_input, ["compare", "empty.aiips"], "exactly 2");
reject!(
    compare_extra_input,
    ["compare", "empty.aiips", "other.aiips", "bad.aiips"],
    "exactly 2"
);
reject!(
    compare_missing_first,
    ["compare", "missing.aiips", "empty.aiips"],
    "missing.aiips"
);
reject!(
    compare_missing_second,
    ["compare", "sketch.aiips", "missing.aiips"],
    "missing.aiips"
);
reject!(
    compare_unknown_option,
    ["compare", "empty.aiips", "other.aiips", "--unknown"],
    "unknown option"
);
reject!(
    compare_force_rejected,
    ["compare", "empty.aiips", "other.aiips", "--force"],
    "unknown option"
);
reject!(
    compare_output_rejected,
    ["compare", "empty.aiips", "other.aiips", "-o", "out"],
    "unknown option"
);
reject!(
    compare_unimplemented_algorithm_rejected,
    ["compare", "empty.aiips", "other.aiips", "-a", "bottomk"],
    "bottomk is not implemented yet"
);
reject!(
    compare_stdin_first,
    ["compare", "-", "empty.aiips"],
    "stdin"
);
reject!(
    compare_stdin_second,
    ["compare", "empty.aiips", "-"],
    "stdin"
);
reject!(
    compare_directory_input,
    ["compare", "sketch.aiips", "."],
    "regular file"
);

#[test]
fn top_level_help() {
    for flag in ["-h", "--help"] {
        let stdout = Workspace::new().success(&[flag]);
        assert!(stdout.contains("sketch"));
        assert!(stdout.contains("compare"));
    }
}

reject!(no_command, [], "Usage:");
reject!(unknown_command, ["unknown"], "unknown command");
reject!(old_interface_removed, ["good.fa"], "unknown command");
reject!(
    old_options_removed,
    ["-k", "3", "good.fa"],
    "unknown command"
);

reject!(
    sketch_missing_output_before_option,
    ["sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "--force"],
    "missing output path"
);

#[test]
fn sketch_dash_prefixed_output_path() {
    let t = Workspace::new();
    t.success(&[
        "sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "./-out",
    ]);
    t.assert_sketch("-out.aiips");
}

#[test]
fn workspace_retries_stale_directory_without_reusing_it() {
    let t = Workspace::new();
    let stale = t.path("stale");
    fs::create_dir(&stale).unwrap();
    fs::write(stale.join("keep"), b"old run").unwrap();
    let fresh = t.path("fresh");
    let mut attempts = [stale.clone(), fresh.clone()].into_iter();
    let reserved = reserve_directory(|| attempts.next().unwrap()).unwrap();
    assert_eq!(reserved, fresh);
    assert_eq!(fs::read(stale.join("keep")).unwrap(), b"old run");
}

#[test]
fn workspace_cleanup_failure_does_not_mask_original_panic() {
    let panic = std::panic::catch_unwind(|| {
        let t = Workspace::new();
        fs::remove_dir_all(&t.0).unwrap();
        panic!("original assertion");
    })
    .unwrap_err();
    assert_eq!(panic.downcast_ref::<&str>(), Some(&"original assertion"));
}

#[cfg(unix)]
#[test]
fn sketch_rejects_hardlink_alias() {
    let t = Workspace::new();
    fs::hard_link(t.path("good.fa"), t.path("alias.aiips")).unwrap();
    t.failure(
        &[
            "sketch",
            "-a",
            "minhash",
            "-k",
            "3",
            "good.fa",
            "-o",
            "alias.aiips",
            "-f",
        ],
        "output resolves to the input",
    );
    assert_eq!(
        fs::read(t.path("alias.aiips")).unwrap(),
        include_bytes!("fixtures/sample.fa")
    );
    t.assert_no_staging_files();
}

#[cfg(unix)]
#[test]
fn sketch_rejects_case_alias_on_case_insensitive_filesystem() {
    use std::os::unix::fs::MetadataExt;
    let t = Workspace::new();
    let data = b">a\nACGT\n";
    t.write("reads.aiips", data);
    match fs::metadata(t.path("READS.AIIPS")) {
        Ok(alias) if alias.ino() == fs::metadata(t.path("reads.aiips")).unwrap().ino() => {
            t.failure(
                &[
                    "sketch",
                    "-a",
                    "minhash",
                    "-k",
                    "3",
                    "reads.aiips",
                    "-o",
                    "READS.AIIPS",
                    "--force",
                ],
                "output resolves to the input",
            );
            assert_eq!(fs::read(t.path("reads.aiips")).unwrap(), data);
        }
        _ => eprintln!("case-insensitive alias scenario unavailable on this filesystem"),
    }
}

#[test]
fn sketch_rejects_directory_like_output_spellings() {
    let t = Workspace::new();
    fs::create_dir(t.path("outdir")).unwrap();
    for name in [
        "outdir/",
        "missing/",
        ".",
        "..",
        "outdir/.",
        "outdir/..",
        "/",
    ] {
        t.failure(
            &["sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", name],
            "directory",
        );
    }
    assert!(!t.path("outdir/.aiips").exists());
    assert!(!t.path("missing").exists());
}

#[test]
fn sketch_existing_directory_never_suggests_force() {
    let t = Workspace::new();
    fs::create_dir(t.path("dir.aiips")).unwrap();
    for name in ["dir", "dir.aiips"] {
        for force in [false, true] {
            let mut args = vec!["sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", name];
            if force {
                args.push("-f");
            }
            let result = t.run(&args);
            assert_eq!(result.status.code(), Some(1), "{args:?}");
            let stderr = String::from_utf8(result.stderr).unwrap();
            assert!(
                stderr.contains("dir.aiips: output is a directory"),
                "{stderr}"
            );
            assert!(!stderr.contains("use -f"), "{stderr}");
            assert!(result.stdout.is_empty());
        }
    }
    t.assert_no_staging_files();
}

#[test]
fn sketch_unsuffixed_directory_name_does_not_block_suffixed_output() {
    let t = Workspace::new();
    fs::create_dir(t.path("results")).unwrap();
    let stdout = t.success(&[
        "sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "results",
    ]);
    assert!(stdout.ends_with("output\tresults.aiips\n"));
    t.assert_sketch("results.aiips");
    assert!(fs::metadata(t.path("results")).unwrap().is_dir());
    assert_eq!(fs::read_dir(t.path("results")).unwrap().count(), 0);
}

#[test]
fn workspace_reservation_gives_up_after_bounded_collisions() {
    let t = Workspace::new();
    let taken = t.path("taken");
    fs::create_dir(&taken).unwrap();
    let mut calls = 0;
    let error = reserve_directory(|| {
        calls += 1;
        taken.clone()
    })
    .unwrap_err();
    assert_eq!(error.kind(), std::io::ErrorKind::AlreadyExists);
    assert_eq!(calls, common::MAX_RESERVE_ATTEMPTS);
}

#[test]
fn sketch_dotfile_suffix_and_compare_detection() {
    let t = Workspace::new();
    for name in [".aiips", ".AIIPS"] {
        let result = t.success(&["sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", name]);
        assert!(result.ends_with(&format!("output\t{name}\n")));
        t.assert_sketch(name);
        assert!(!t.path(&format!("{name}.aiips")).exists());
        assert!(t
            .success(&["compare", name, name])
            .contains(&format!("\n{name}\t{name}\t1.000000\t")));
        fs::remove_file(t.path(name)).unwrap(); // Names may alias on APFS.
    }
    t.success(&[
        "sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", ".hidden",
    ]);
    t.assert_sketch(".hidden.aiips");
}

#[cfg(unix)]
#[test]
fn sketch_force_preserves_permission_bits() {
    use std::os::unix::fs::PermissionsExt;
    let t = Workspace::new();
    t.write("out.aiips", b"old");
    fs::set_permissions(t.path("out.aiips"), fs::Permissions::from_mode(0o640)).unwrap();
    t.success(&[
        "sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "out", "--force",
    ]);
    assert_eq!(
        fs::metadata(t.path("out.aiips"))
            .unwrap()
            .permissions()
            .mode()
            & 0o7777,
        0o640
    );
    t.assert_sketch("out.aiips");
    t.assert_no_staging_files();
}

#[test]
fn command_help_precedes_invalid_options_and_values() {
    let t = Workspace::new();
    for subcommand in ["sketch", "compare"] {
        for args in [
            vec![subcommand, "--bogus", "--help"],
            vec![subcommand, "-k", "0", "-h"],
            vec![subcommand, "-o", "--help"],
            vec![subcommand, "--help", "--bogus"],
        ] {
            assert!(t
                .success(&args)
                .contains(&format!("Usage: aiips {subcommand}")));
        }
    }
}

#[test]
fn help_after_double_dash_is_a_filename() {
    let t = Workspace::new();
    t.write("--help", b">a\nACGT\n");
    t.success(&[
        "sketch", "-a", "minhash", "-k", "3", "-o", "out", "--", "--help",
    ]);
    // A sketch file literally named --help is a positional input after --.
    fs::copy(t.path("out.aiips"), t.path("--help")).unwrap();
    assert!(t
        .success(&["compare", "--", "--help", "out.aiips"])
        .contains("\n--help\tout.aiips\t1.000000\t"));
    t.failure(
        &["compare", "--bogus", "--", "--help", "out.aiips"],
        "unknown option",
    );
}

#[cfg(unix)]
#[test]
fn sketch_closed_stdout_rolls_back_new_and_forced_outputs() {
    use std::os::unix::fs::PermissionsExt;
    use std::process::Stdio;
    let t = Workspace::new();
    for force in [false, true] {
        if force {
            t.write("out.aiips", b"old result");
            fs::set_permissions(t.path("out.aiips"), fs::Permissions::from_mode(0o640)).unwrap();
        }
        // Drop the read end before starting the child: deterministic EPIPE, no shell race.
        let (reader, writer) = std::io::pipe().unwrap();
        drop(reader);
        let mut command = Command::new(env!("CARGO_BIN_EXE_aiips"));
        command
            .current_dir(&t.0)
            .args(["sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "out"])
            .stdout(Stdio::from(writer));
        if force {
            command.arg("-f");
        }
        let result = command.output().unwrap();
        assert_eq!(result.status.code(), Some(1));
        assert!(String::from_utf8_lossy(&result.stderr).contains("stdout:"));
        if force {
            assert_eq!(fs::read(t.path("out.aiips")).unwrap(), b"old result");
            assert_eq!(
                fs::metadata(t.path("out.aiips"))
                    .unwrap()
                    .permissions()
                    .mode()
                    & 0o7777,
                0o640
            );
        } else {
            assert!(!t.path("out.aiips").exists());
        }
        t.assert_no_staging_files();
    }
}

#[cfg(unix)]
#[test]
fn sketch_closed_stdout_restores_original_symlink_and_target() {
    use std::process::Stdio;
    let t = Workspace::new();
    t.write("original", b"keep target");
    std::os::unix::fs::symlink("original", t.path("out.aiips")).unwrap();
    let (reader, writer) = std::io::pipe().unwrap();
    drop(reader);
    let result = Command::new(env!("CARGO_BIN_EXE_aiips"))
        .current_dir(&t.0)
        .args([
            "sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", "out", "--force",
        ])
        .stdout(Stdio::from(writer))
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(1));
    assert_eq!(
        fs::read_link(t.path("out.aiips")).unwrap(),
        PathBuf::from("original")
    );
    assert_eq!(fs::read(t.path("original")).unwrap(), b"keep target");
    t.assert_no_staging_files();
}

#[test]
fn repeated_k_is_rejected_for_all_aliases() {
    let t = Workspace::new();
    for flags in [["-k", "-k"], ["-k", "--kmer-size"], ["--kmer-size", "-k"]] {
        let mut args = vec!["sketch", "-a", "minhash", "good.fa", "-o", "out"];
        args.extend([flags[0], "15", flags[1], "31"]);
        t.failure(&args, "k-mer length may only be specified once");
    }
    assert!(!t.path("out.aiips").exists());
}

#[test]
fn sketch_existing_output_is_rejected_before_parsing_invalid_input() {
    let t = Workspace::new();
    t.write("out.aiips", b"keep existing result");
    t.failure(
        &["sketch", "-a", "minhash", "-k", "3", "bad.fq", "-o", "out"],
        "output already exists",
    );
    assert_eq!(
        fs::read(t.path("out.aiips")).unwrap(),
        b"keep existing result"
    );
    t.assert_no_staging_files();
}

#[test]
fn sketch_alias_is_rejected_before_parsing_invalid_input() {
    let t = Workspace::new();
    t.write("invalid.aiips", b"not FASTX");
    t.failure(
        &[
            "sketch",
            "-a",
            "minhash",
            "-k",
            "3",
            "invalid.aiips",
            "-o",
            "invalid.aiips",
            "--force",
        ],
        "output resolves to the input",
    );
    assert_eq!(fs::read(t.path("invalid.aiips")).unwrap(), b"not FASTX");
}

#[cfg(unix)]
#[test]
fn sketch_symlink_is_rejected_before_parsing_invalid_input() {
    let t = Workspace::new();
    std::os::unix::fs::symlink("empty.aiips", t.path("out.aiips")).unwrap();
    t.failure(
        &["sketch", "-a", "minhash", "-k", "3", "bad.fq", "-o", "out"],
        "output already exists",
    );
    assert_eq!(
        fs::read_link(t.path("out.aiips")).unwrap(),
        PathBuf::from("empty.aiips")
    );
}

#[cfg(unix)]
#[test]
fn sketch_special_output_is_rejected_before_parsing_invalid_input() {
    let t = Workspace::new();
    let _socket = std::os::unix::net::UnixListener::bind(t.path("out.aiips")).unwrap();
    t.failure(
        &[
            "sketch", "-a", "minhash", "-k", "3", "bad.fq", "-o", "out", "--force",
        ],
        "regular file or symlink",
    );
    t.assert_no_staging_files();
}

#[test]
fn help_precedes_duplicate_k_but_not_after_double_dash() {
    let t = Workspace::new();
    assert!(t
        .success(&["compare", "-k", "15", "--kmer-size", "31", "--help"])
        .contains("Usage: aiips compare"));
    t.failure(
        &[
            "sketch", "-a", "minhash", "-k", "15", "-k", "31", "-o", "out", "--", "--help",
        ],
        "may only be specified once",
    );
    t.failure(
        &["sketch", "-a", "minhash", "-k", "--", "--help"],
        "positive integer",
    );
}

#[test]
fn sketch_missing_parent_is_rejected_before_parsing_invalid_input() {
    let t = Workspace::new();
    t.failure(
        &[
            "sketch",
            "-a",
            "minhash",
            "-k",
            "3",
            "bad.fq",
            "-o",
            "missing/out",
        ],
        &native("missing/out.aiips"),
    );
    assert!(!t.path("missing").exists());
}

#[test]
fn workspace_removes_retained_backup_and_directory_during_unwind() {
    let directory = TestDirectory::new();
    let root = directory.to_path_buf();
    let result = std::panic::catch_unwind(move || {
        fs::write(directory.join("retained-backup.tmp"), b"original result").unwrap();
        fs::create_dir(directory.join("failed-replacement.tmp")).unwrap();
        panic!("simulated failed assertion");
    });
    assert!(result.is_err());
    assert!(!root.exists());
}

#[cfg(windows)]
#[test]
fn sketch_windows_accepts_both_separator_styles_for_path_validation() {
    let t = Workspace::new();
    for spelling in [
        "missing/",
        "missing\\",
        "missing/.",
        "missing\\.",
        "missing/..",
        "missing\\..",
    ] {
        t.failure(
            &[
                "sketch", "-a", "minhash", "-k", "3", "good.fa", "-o", spelling,
            ],
            "directory",
        );
    }
}

// ---------------------------------------------------------------- FracMinHash

/// Library-built FracMinHash sketch bytes, an independent path to the CLI's.
fn frac_bytes(fasta: &[u8], k: usize, scaled: u64, seed: u64) -> Vec<u8> {
    use aiips::sketch::fracminhash::{FracMinHashParams, FracMinHashSketch};
    let p = FracMinHashParams::new(std::num::NonZeroUsize::new(k).unwrap(), scaled, seed).unwrap();
    FracMinHashSketch::build(fasta, p).unwrap().to_bytes()
}

fn frac_library_row(t: &Workspace, a: &str, b: &str) -> String {
    use aiips::sketch::{fracminhash, parse, Sketch};
    let load = |name: &str| match parse(&fs::read(t.path(name)).unwrap()).unwrap() {
        Sketch::FracMinHash(sketch) => sketch,
        other => panic!("{other:?}"),
    };
    let e = fracminhash::estimate(&load(a), &load(b)).unwrap();
    format!(
        "{a}\t{b}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t-1\t-1\t{:.6}\t{:.6}\n",
        e.jaccard,
        e.containment_a_in_b,
        e.containment_b_in_a,
        e.ani,
        e.ani_containment_a_in_b,
        e.ani_containment_b_in_a
    )
}

fn frac(t: &Workspace, input: &str, name: &str, extra: &[&str]) -> String {
    let mut args = vec!["sketch", "-a", "fracminhash", input, "-o", name];
    args.extend(extra);
    t.success(&args)
}

#[test]
fn sketch_fracminhash_records_parameters_and_matches_library() {
    let t = Workspace::new();
    let stdout = frac(
        &t,
        "good.fa",
        "a",
        &["-k", "3", "--scaled", "1", "--seed", "7"],
    );
    assert_eq!(
        stdout,
        format!("{STATS}\ngood.fa\t2\t13\t2\t3\t7\noutput\ta.aiips\n")
    );
    let text = fs::read_to_string(t.path("a.aiips")).unwrap();
    let expected = frac_bytes(include_bytes!("fixtures/sample.fa"), 3, 1, 7);
    assert_eq!(text.as_bytes(), expected);
    assert!(text.starts_with(
        "aiips-sketch v1 algorithm=fracminhash k=3 scaled=1 seed=7 hash=murmur3-fmix64 canonical=true size="
    ));
    let size: usize = text
        .lines()
        .next()
        .unwrap()
        .rsplit_once("size=")
        .unwrap()
        .1
        .parse()
        .unwrap();
    assert!(size > 0);
    assert_eq!(text.lines().count(), 2 + size);
    // FASTA, FASTQ and gzip inputs with identical sequences give identical sketches.
    for (input, name) in [("good.fq", "b"), ("good.fa.gz", "c"), ("good.fq.gz", "d")] {
        frac(
            &t,
            input,
            name,
            &["-k", "3", "--scaled", "1", "--seed", "7"],
        );
        assert_eq!(
            fs::read_to_string(t.path(&format!("{name}.aiips"))).unwrap(),
            text
        );
    }
}

#[test]
fn sketch_fracminhash_defaults_and_empty_sketch_warning() {
    let t = Workspace::new();
    let sequence = "GATTACAGGCTTAACGTTAGCCATGCAAGTCTAGGCATCG";
    t.write("long.fa", format!(">l\n{sequence}\n").as_bytes());
    let expected = frac_bytes(format!(">l\n{sequence}\n").as_bytes(), 21, 1000, 42);
    let header = "aiips-sketch v1 algorithm=fracminhash k=21 scaled=1000 seed=42 hash=murmur3-fmix64 canonical=true size=";
    assert!(String::from_utf8_lossy(&expected).starts_with(header));
    let output = t.run(&["sketch", "-a", "fracminhash", "long.fa", "-o", "out"]);
    assert!(output.status.success());
    assert_eq!(
        String::from_utf8(output.stdout).unwrap(),
        format!("{STATS}\nlong.fa\t1\t40\t0\t21\t20\noutput\tout.aiips\n")
    );
    assert_eq!(fs::read(t.path("out.aiips")).unwrap(), expected);
    // Fixed seed: none of the 20 distinct 21-mers falls under the threshold.
    assert!(String::from_utf8_lossy(&expected).starts_with(&format!("{header}0\n")));
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert_eq!(
        stderr,
        "aiips: warning: out.aiips: no k-mer hash is at or below the scaled=1000 threshold; the sketch is empty and cannot be compared\n"
    );
    for pair in [["out.aiips", "out.aiips"], ["long.fa", "out.aiips"]] {
        t.failure(
            &["compare", pair[0], pair[1]],
            ": fracminhash sketch has no hash values; containment is undefined; use a smaller --scaled",
        );
    }
    // Inputs without valid k-mers are still refused, creating nothing.
    t.failure(
        &["sketch", "-a", "fracminhash", "good.fa", "-o", "none"],
        "no valid k-mers of length 21",
    );
    assert!(!t.path("none.aiips").exists());
}

#[test]
fn sketch_validates_fracminhash_parameters() {
    let t = Workspace::new();
    for (args, message) in [
        (
            vec!["-a", "fracminhash", "--scaled", "0"],
            "scaled must be an integer from 1",
        ),
        (
            vec!["-a", "fracminhash", "--scaled", "-1"],
            "scaled must be an integer from 1",
        ),
        (
            vec!["-a", "fracminhash", "--scaled", "abc"],
            "scaled must be an integer",
        ),
        (
            vec!["-a", "fracminhash", "--scaled", "18446744073709551616"],
            "scaled must be",
        ),
        (
            vec!["-a", "fracminhash", "--scaled"],
            "scaled must be an integer",
        ),
        (
            vec!["-a", "fracminhash", "--scaled", "2", "--scaled", "2"],
            "scaled may only be specified once",
        ),
        (vec!["-a", "fracminhash", "-k", "33"], "1 to 32"),
        (
            vec!["-a", "fracminhash", "-s", "8"],
            "-s/--num-hashes applies only to minhash, not fracminhash",
        ),
        (
            vec!["-a", "minhash", "--scaled", "10"],
            "--scaled applies only to fracminhash, not minhash",
        ),
        (
            vec!["--scaled", "10"],
            "--scaled applies only to fracminhash",
        ),
    ] {
        let mut full = vec!["sketch", "good.fa", "-o", "out"];
        full.extend(args);
        t.failure(&full, message);
    }
    assert!(!t.path("out.aiips").exists());
}

#[test]
fn compare_fracminhash_matches_library_and_exact_cases() {
    let t = Workspace::new();
    let shared = "GATTACAGGCTTAACGTTAGCCATGCAAGTCTAGGCATCGTTGACCTAGGATCCATGGCA";
    let extra = "TTGCAACGGTACCGATTGCATGCATCCGGAATTCGCGTATAGCTAGCATCGATCGGCTAA";
    t.write("a.fa", format!(">a\n{shared}\n").as_bytes());
    t.write("b.fa", format!(">b\n{shared}{extra}\n").as_bytes());
    for (name, scaled) in [("1", "1"), ("2", "2")] {
        frac(
            &t,
            "a.fa",
            &format!("a{name}"),
            &["-k", "11", "--scaled", scaled],
        );
        frac(
            &t,
            "b.fa",
            &format!("b{name}"),
            &["-k", "11", "--scaled", scaled],
        );
        let (a, b) = (format!("a{name}.aiips"), format!("b{name}.aiips"));
        let forward = t.success(&["compare", &a, &b]);
        assert!(
            forward.ends_with(&frac_library_row(&t, &a, &b)),
            "{forward}"
        );
        let backward = t.success(&["compare", &b, &a]);
        assert!(backward.ends_with(&frac_library_row(&t, &b, &a)));
    }
    // At scaled=1 every k-mer is kept, and every 11-mer of A occurs in B.
    let row = t.success(&["compare", "a1.aiips", "b1.aiips"]);
    let fields: Vec<&str> = row.lines().last().unwrap().split('\t').collect();
    assert_eq!(fields[3], "1.000000", "A is contained in B");
    assert_eq!(fields[8], "1.000000");
    // Identical, reverse-complement inputs report exact ones.
    t.write("f.fa", b">f\nGATTACAGGCTTAACGTTAG\n");
    t.write("r.fa", b">r\nCTAACGTTAAGCCTGTAATC\n");
    frac(&t, "f.fa", "f", &["-k", "7", "--scaled", "1"]);
    frac(&t, "r.fa", "r", &["-k", "7", "--scaled", "1"]);
    assert!(t.success(&["compare", "f.aiips", "r.aiips"]).ends_with(
        "\nf.aiips\tr.aiips\t1.000000\t1.000000\t1.000000\t1.000000\t-1\t-1\t1.000000\t1.000000\n"
    ));
}

#[test]
fn compare_fracminhash_sequence_inputs_use_options_or_sketch_parameters() {
    let t = Workspace::new();
    frac(&t, "good.fa", "f", &["-k", "3", "--scaled", "1"]);
    let metrics = |out: &str| -> String {
        let line = out.lines().last().unwrap();
        line.splitn(3, '\t').nth(2).unwrap().to_owned()
    };
    let sketched = t.success(&["compare", "f.aiips", "f.aiips"]);
    for args in [
        vec![
            "compare",
            "good.fa",
            "good.fq.gz",
            "-a",
            "fracminhash",
            "-k",
            "3",
            "--scaled",
            "1",
        ],
        vec!["compare", "f.aiips", "good.fq"],
        vec![
            "compare", "good.fq", "f.aiips", "--scaled", "1", "--seed", "42", "-k", "3",
        ],
        vec!["compare", "f.aiips", "good.fa", "-a", "fracminhash"],
    ] {
        let out = t.success(&args);
        assert_eq!(metrics(&out), metrics(&sketched), "{args:?}");
        assert!(out.starts_with(&format!("{STATS}\n")));
    }
    // Two sequence inputs: options are validated against the chosen algorithm.
    for (args, message) in [
        (
            vec!["--scaled", "1"],
            "--scaled applies only to fracminhash",
        ),
        (
            vec!["-a", "fracminhash", "-s", "8"],
            "-s/--num-hashes applies only to minhash",
        ),
        (
            vec!["-a", "fracminhash", "--scaled", "0"],
            "scaled must be an integer from 1",
        ),
        (vec!["-a", "minhash", "-s", "1"], "from 2 to 1000000"),
    ] {
        let mut full = vec!["compare", "good.fa", "good.fq", "-k", "3"];
        full.extend(args);
        t.failure(&full, message);
    }
}

#[test]
fn compare_options_must_match_sketch_inputs() {
    let t = Workspace::new();
    frac(&t, "good.fa", "f", &["-k", "3", "--scaled", "1"]);
    t.sketch("good.fa", "m", &["-k", "3", "-s", "8"]);
    for (sketch, args, message) in [
        (
            "f.aiips",
            vec!["-a", "minhash"],
            "-a minhash does not match algorithm=fracminhash recorded in f.aiips",
        ),
        (
            "f.aiips",
            vec!["--scaled", "2"],
            "--scaled 2 does not match scaled=1 recorded in f.aiips",
        ),
        (
            "f.aiips",
            vec!["--seed", "1"],
            "--seed 1 does not match seed=42 recorded in f.aiips",
        ),
        (
            "f.aiips",
            vec!["-k", "4"],
            "-k 4 does not match k=3 recorded in f.aiips",
        ),
        (
            "f.aiips",
            vec!["-s", "8"],
            "-s/--num-hashes applies only to minhash",
        ),
        (
            "m.aiips",
            vec!["-a", "fracminhash"],
            "-a fracminhash does not match algorithm=minhash recorded in m.aiips",
        ),
        (
            "m.aiips",
            vec!["--scaled", "1"],
            "--scaled applies only to fracminhash",
        ),
        (
            "m.aiips",
            vec!["-s", "9"],
            "-s 9 does not match num_hashes=8 recorded in m.aiips",
        ),
        (
            "m.aiips",
            vec!["--seed", "2"],
            "--seed 2 does not match seed=42 recorded in m.aiips",
        ),
    ] {
        for pair in [[sketch, "good.fa"], ["good.fa", sketch]] {
            let mut full = vec!["compare", pair[0], pair[1]];
            full.extend(args.iter().copied());
            t.failure(&full, message);
        }
    }
    // Options are also checked when both inputs are sketches.
    t.failure(
        &["compare", "f.aiips", "f.aiips", "--scaled", "2"],
        "--scaled 2 does not match scaled=1",
    );
    t.success(&[
        "compare", "m.aiips", "good.fa", "-a", "minhash", "-s", "8", "--seed", "42",
    ]);
}

#[test]
fn compare_halts_on_mixed_algorithms_and_fracminhash_mismatches() {
    let t = Workspace::new();
    frac(
        &t,
        "good.fa",
        "f",
        &["-k", "3", "--scaled", "1", "--seed", "1"],
    );
    t.sketch("good.fa", "m", &["-k", "3"]);
    t.failure(
        &["compare", "m.aiips", "f.aiips"],
        "comparison requires two sketches of one algorithm, minhash or fracminhash (m.aiips: minhash, f.aiips: fracminhash)",
    );
    t.failure(
        &["compare", "f.aiips", "m.aiips"],
        "(f.aiips: fracminhash, m.aiips: minhash)",
    );
    for (name, extra, message) in [
        (
            "k",
            ["-k", "4", "--scaled", "1", "--seed", "1"],
            "k (3 vs 4)",
        ),
        (
            "s",
            ["-k", "3", "--scaled", "2", "--seed", "1"],
            "scaled (1 vs 2)",
        ),
        (
            "seed",
            ["-k", "3", "--scaled", "1", "--seed", "2"],
            "seed (1 vs 2)",
        ),
    ] {
        frac(&t, "good.fa", name, &extra);
        t.failure(
            &["compare", "f.aiips", &format!("{name}.aiips")],
            &format!("sketch parameters differ: {message}; re-sketch with identical parameters"),
        );
    }
}

#[test]
fn compare_rejects_corrupt_fracminhash_sketches() {
    let t = Workspace::new();
    frac(&t, "good.fa", "f", &["-k", "3", "--scaled", "1"]);
    let text = fs::read_to_string(t.path("f.aiips")).unwrap();
    let mut lines: Vec<&str> = text.lines().collect();
    let last = lines.pop().unwrap();
    let swapped = format!(
        "{}\n{last}\n{}\n",
        lines[..lines.len() - 1].join("\n"),
        lines[lines.len() - 1]
    );
    t.write("unsorted.aiips", swapped.as_bytes());
    t.failure(
        &["compare", "unsorted.aiips", "f.aiips"],
        "not strictly increasing",
    );
    t.write("short.aiips", format!("{}\n", lines.join("\n")).as_bytes());
    t.failure(&["compare", "short.aiips", "f.aiips"], "expected size=");
    t.write(
        "bounded.aiips",
        text.replace("scaled=1 ", "scaled=18446744073709551615 ")
            .as_bytes(),
    );
    t.failure(
        &["compare", "bounded.aiips", "bounded.aiips"],
        "exceeds max_hash",
    );
}

#[test]
fn compare_reports_saturated_fracminhash_size_instead_of_panicking() {
    // 21 · size + 4096 saturates the read bound to u64::MAX for these sizes;
    // the reader must still read the whole body and report the value count.
    let t = Workspace::new();
    frac(&t, "good.fa", "f", &["-k", "3", "--scaled", "1"]);
    let text = fs::read_to_string(t.path("f.aiips")).unwrap();
    let declared = text.lines().next().unwrap().rsplit_once("size=").unwrap().1;
    let found = text.lines().count() - 2;
    for size in [u64::MAX, u64::MAX / 21 + 1] {
        t.write(
            "huge.aiips",
            text.replacen(&format!("size={declared}\n"), &format!("size={size}\n"), 1)
                .as_bytes(),
        );
        let output = t.run(&["compare", "huge.aiips", "f.aiips"]);
        assert_eq!(output.status.code(), Some(1), "size={size}");
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert!(
            stderr.contains(&format!(
                "huge.aiips: invalid sketch: expected size={size} hash values, found {found}"
            )),
            "{stderr}"
        );
        assert!(output.stdout.is_empty());
    }
}
