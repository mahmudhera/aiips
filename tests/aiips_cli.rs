//! Command-line behavior of the aiips algorithm (`-a aiips`).

use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};
mod common;
use aiips::hash::SplitMix64;
use common::TestDirectory;

const AIIPS_METRICS: &str = "input_a\tinput_b\tani\taf_a_to_b\taf_b_to_a\tani_a\tani_b\tanchors_matched\tanchors\tprobes_matched\tprobes";

struct Workspace(TestDirectory);

/// Uniform random bases from a fixed seed.
fn genome(length: usize, seed: u64) -> Vec<u8> {
    let mut rng = SplitMix64::new(seed);
    (0..length)
        .map(|_| b"ACGT"[(rng.next().unwrap() % 4) as usize])
        .collect()
}

fn fasta(records: &[&[u8]]) -> Vec<u8> {
    let mut out = Vec::new();
    for (i, record) in records.iter().enumerate() {
        out.extend_from_slice(format!(">r{i}\n").as_bytes());
        for line in record.chunks(80) {
            out.extend_from_slice(line);
            out.push(b'\n');
        }
    }
    out
}

impl Workspace {
    fn new() -> Self {
        let workspace = Self(TestDirectory::new());
        let a = genome(60_000, 1);
        // B: A's second half, then unrelated sequence.
        let mut b = a[30_000..].to_vec();
        b.extend(genome(10_000, 2));
        workspace.write("a.fa", &fasta(&[&a]));
        workspace.write("b.fa", &fasta(&[&b]));
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
        assert!(output.stderr.is_empty(), "{args:?}");
        String::from_utf8(output.stdout).unwrap()
    }

    fn failure(&self, args: &[&str], message: &str) {
        let output = self.run(args);
        assert_eq!(output.status.code(), Some(1), "{args:?}");
        assert!(output.stdout.is_empty(), "{args:?}: unexpected stdout");
        let error = String::from_utf8(output.stderr).unwrap();
        assert!(error.contains(message), "{args:?}: {error}");
    }
}

/// The comparison row's fields after the two labels.
fn metrics(stdout: &str) -> Vec<String> {
    let lines: Vec<&str> = stdout.lines().collect();
    let header = lines.iter().position(|l| *l == AIIPS_METRICS).unwrap();
    lines[header + 1]
        .split('\t')
        .skip(2)
        .map(String::from)
        .collect()
}

#[test]
fn sketch_writes_a_deterministic_aiips_sketch() {
    let t = Workspace::new();
    let out = t.success(&["sketch", "a.fa", "-a", "aiips", "-o", "a"]);
    assert_eq!(
        out,
        "input\trecords\tbases\tnon_acgt_bases\tk\tvalid_kmers\na.fa\t1\t60000\t0\t15\t59986\noutput\ta.aiips\n"
    );
    let first = fs::read(t.path("a.aiips")).unwrap();
    let header = String::from_utf8_lossy(&first);
    assert!(header.starts_with(
        "aiips-sketch v1 algorithm=aiips anchor_k=21 anchor_scaled=500 probe_k=15 probe_scaled=100 seed=42 hash=murmur3-fmix64 canonical=true contigs=1 anchors="
    ));
    t.success(&["sketch", "a.fa", "-a", "aiips", "-o", "a", "-f"]);
    assert_eq!(fs::read(t.path("a.aiips")).unwrap(), first);
}

#[test]
fn compare_reports_identity_and_alignment_fraction() {
    let t = Workspace::new();
    let out = t.success(&["compare", "a.fa", "b.fa", "-a", "aiips"]);
    let m = metrics(&out);
    // Shared half is identical: ANI 1. A aligns about 30/60, B about 30/40.
    assert_eq!(m[0], "1.000000");
    let (af_a, af_b): (f64, f64) = (m[1].parse().unwrap(), m[2].parse().unwrap());
    assert!((af_a - 0.5).abs() < 0.01, "{out}");
    assert!((af_b - 0.75).abs() < 0.01, "{out}");
    assert_eq!(m[7], m[8], "every probe matches: {out}");
}

#[test]
fn sketches_and_sequences_give_identical_comparisons() {
    let t = Workspace::new();
    t.success(&["sketch", "a.fa", "-a", "aiips", "-o", "a"]);
    t.success(&["sketch", "b.fa", "-a", "aiips", "-o", "b"]);
    let from_files = metrics(&t.success(&["compare", "a.fa", "b.fa", "-a", "aiips"]));
    let from_sketches = metrics(&t.success(&["compare", "a.aiips", "b.aiips"]));
    let mixed = metrics(&t.success(&["compare", "a.aiips", "b.fa"]));
    assert_eq!(from_files, from_sketches);
    assert_eq!(from_files, mixed);
}

#[test]
fn walk_options_are_accepted_and_change_only_the_walk() {
    let t = Workspace::new();
    for extra in [
        vec!["--max-gap", "0"],
        vec!["--min-island-ani", "0.9"],
        vec!["--symmetric-probes"],
        vec!["--bidirectional"],
        vec!["--symmetric-probes", "--bidirectional", "--max-gap", "10"],
    ] {
        let mut args = vec!["compare", "a.fa", "b.fa", "-a", "aiips"];
        args.extend(extra);
        assert_eq!(metrics(&t.success(&args))[0], "1.000000");
    }
}

#[test]
fn custom_sketch_parameters_are_recorded_and_enforced() {
    let t = Workspace::new();
    let custom = [
        "--anchor-k",
        "19",
        "--anchor-scaled",
        "200",
        "--probe-k",
        "13",
        "--probe-scaled",
        "50",
        "--seed",
        "7",
    ];
    let mut args = vec!["sketch", "a.fa", "-a", "aiips", "-o", "a"];
    args.extend(custom);
    t.success(&args);
    let header = fs::read_to_string(t.path("a.aiips")).unwrap();
    assert!(header.contains("anchor_k=19 anchor_scaled=200 probe_k=13 probe_scaled=50 seed=7"));
    // Sequence inputs adopt the sketch's parameters; given options must match.
    t.success(&["compare", "a.aiips", "b.fa"]);
    t.failure(
        &["compare", "a.aiips", "b.fa", "--anchor-k", "21"],
        "--anchor-k 21 does not match anchor_k=19",
    );
    t.success(&["sketch", "b.fa", "-a", "aiips", "-o", "b"]);
    t.failure(
        &["compare", "a.aiips", "b.aiips"],
        "sketch parameters differ: anchor_k (19 vs 21), anchor_scaled (200 vs 500), probe_k (13 vs 15), probe_scaled (50 vs 100), seed (7 vs 42)",
    );
}

#[test]
fn options_are_validated() {
    let t = Workspace::new();
    for (args, message) in [
        (
            vec!["sketch", "a.fa", "-a", "aiips", "-o", "x", "-k", "21"],
            "-k/--kmer-size does not apply to aiips",
        ),
        (
            vec!["sketch", "a.fa", "-a", "aiips", "-o", "x", "--scaled", "9"],
            "--scaled does not apply to aiips",
        ),
        (
            vec!["sketch", "a.fa", "-o", "x", "--anchor-k", "21"],
            "--anchor-k applies only to aiips, not minhash",
        ),
        (
            vec![
                "sketch",
                "a.fa",
                "-a",
                "aiips",
                "-o",
                "x",
                "--probe-k",
                "21",
            ],
            "probe k must be from 1 to anchor k - 1",
        ),
        (
            vec![
                "sketch",
                "a.fa",
                "-a",
                "aiips",
                "-o",
                "x",
                "--anchor-k",
                "33",
            ],
            "anchor k must be from 2 to 32",
        ),
        (
            vec![
                "sketch",
                "a.fa",
                "-a",
                "aiips",
                "-o",
                "x",
                "--probe-scaled",
                "0",
            ],
            "scaled must be an integer",
        ),
        (
            vec!["sketch", "a.fa", "-a", "aiips", "-o", "x", "--max-gap", "3"],
            "unknown option: --max-gap",
        ),
        (
            vec![
                "compare",
                "a.fa",
                "b.fa",
                "-a",
                "aiips",
                "--min-island-ani",
                "1.5",
            ],
            "min island ANI must be a number from 0 to 1",
        ),
        (
            vec![
                "compare",
                "a.fa",
                "b.fa",
                "-a",
                "aiips",
                "--bidirectional",
                "--bidirectional",
            ],
            "--bidirectional may only be specified once",
        ),
        (
            vec!["compare", "a.fa", "b.fa", "--max-gap", "3"],
            "--max-gap applies only to aiips, not minhash",
        ),
        (
            vec![
                "compare",
                "a.fa",
                "b.fa",
                "-a",
                "fracminhash",
                "--symmetric-probes",
            ],
            "--symmetric-probes applies only to aiips, not fracminhash",
        ),
    ] {
        t.failure(&args, message);
    }
    assert!(!t.path("x.aiips").exists());
}

#[test]
fn unrelated_inputs_report_no_identity() {
    let t = Workspace::new();
    t.write("c.fa", &fasta(&[&genome(20_000, 3)]));
    let m = metrics(&t.success(&["compare", "a.fa", "c.fa", "-a", "aiips"]));
    assert_eq!(&m[..3], ["-1", "0.000000", "0.000000"]);
}

#[test]
fn help_lists_aiips_options() {
    let t = Workspace::new();
    let sketch = t.success(&["sketch", "--help"]);
    for option in [
        "--anchor-k",
        "--anchor-scaled",
        "--probe-k",
        "--probe-scaled",
    ] {
        assert!(sketch.contains(option), "{option}");
    }
    let compare = t.success(&["compare", "--help"]);
    for option in [
        "--max-gap",
        "--min-island-ani",
        "--symmetric-probes",
        "--bidirectional",
    ] {
        assert!(compare.contains(option), "{option}");
    }
}
