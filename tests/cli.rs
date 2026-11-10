use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

const STATS: &str = "input\trecords\tbases\tnon_acgt_bases\tk\tvalid_kmers";
const METRICS: &str =
    "input_a\tinput_b\tjaccard\tcontainment_a_in_b\tcontainment_b_in_a\tani\taf_a_to_b\taf_b_to_a";

struct Workspace(PathBuf);

impl Workspace {
    fn new() -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "aiips-cli-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        let workspace = Self(path);
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

    fn assert_empty(&self, name: &str) {
        let metadata = fs::metadata(self.path(name)).unwrap();
        assert!(metadata.is_file());
        assert_eq!(metadata.len(), 0);
    }

    fn assert_no_staging_files(&self) {
        assert!(fs::read_dir(&self.0).unwrap().all(|entry| !entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .starts_with(".aiips-")));
    }
}

impl Drop for Workspace {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
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
            let stdout = t.success(&["sketch", $input, "-o", "out", "-k", "3"]);
            assert_eq!(
                stdout,
                format!("{STATS}\n{}\t2\t13\t2\t3\t7\noutput\tout.aiips\n", $input)
            );
            t.assert_empty("out.aiips");
        }
    };
}

sketch_format!(sketch_fasta_known_counts, "good.fa");
sketch_format!(sketch_fastq_known_counts, "good.fq");
sketch_format!(sketch_gzip_fasta, "good.fa.gz");
sketch_format!(sketch_gzip_fastq, "good.fq.gz");

#[test]
fn sketch_default_k_and_exact_suffix() {
    let t = Workspace::new();
    let stdout = t.success(&["sketch", "good.fa", "-o", "out.aiips"]);
    assert!(stdout.contains("good.fa\t2\t13\t2\t21\t0\n"));
    t.assert_empty("out.aiips");
    assert!(!t.path("out.aiips.aiips").exists());
}

#[test]
fn sketch_long_options_and_uppercase_suffix() {
    let t = Workspace::new();
    let stdout = t.success(&[
        "sketch",
        "--kmer-size",
        "3",
        "--output",
        "out.AIIPS",
        "good.fa",
    ]);
    assert!(stdout.ends_with("output\tout.AIIPS\n"));
    t.assert_empty("out.AIIPS");
}

#[test]
fn sketch_appends_to_other_extension() {
    let t = Workspace::new();
    t.success(&["sketch", "good.fa", "-o", "out.txt"]);
    t.assert_empty("out.txt.aiips");
}

#[test]
fn sketch_ambiguity_and_record_boundaries() {
    let t = Workspace::new();
    t.write("runs.fa", b">a\nAC\nGTNacg\n>b\nTT\n>c\nAA\n");
    let stdout = t.success(&["sketch", "runs.fa", "-k", "3", "-o", "out"]);
    assert!(stdout.contains("runs.fa\t3\t12\t1\t3\t3\n"));
}

#[test]
fn sketch_existing_output_refused_and_preserved() {
    let t = Workspace::new();
    t.write("out.aiips", b"existing result");
    t.failure(&["sketch", "good.fa", "-o", "out"], "use -f or --force");
    assert_eq!(fs::read(t.path("out.aiips")).unwrap(), b"existing result");
}

#[test]
fn sketch_force_replaces_and_cleans_staging() {
    for flag in ["-f", "--force"] {
        let t = Workspace::new();
        t.write("out.aiips", b"existing result");
        t.success(&["sketch", "good.fa", "-o", "out", flag]);
        t.assert_empty("out.aiips");
        t.assert_no_staging_files();
    }
}

#[test]
fn sketch_force_creates_new_output() {
    let t = Workspace::new();
    t.success(&["sketch", "good.fa", "-o", "out", "--force"]);
    t.assert_empty("out.aiips");
    t.assert_no_staging_files();
}

#[test]
fn sketch_force_preserves_output_on_parse_failure() {
    let t = Workspace::new();
    t.write("out.aiips", b"existing result");
    t.failure(&["sketch", "bad.fq", "-o", "out", "--force"], "bad.fq");
    assert_eq!(fs::read(t.path("out.aiips")).unwrap(), b"existing result");
    t.assert_no_staging_files();
}

#[test]
fn sketch_failure_does_not_create_output() {
    let t = Workspace::new();
    t.failure(&["sketch", "bad.fq", "-o", "out"], "bad.fq");
    assert!(!t.path("out.aiips").exists());
}

#[test]
fn sketch_missing_destination_directory() {
    let t = Workspace::new();
    for force in [vec![], vec!["--force"]] {
        let mut args = vec!["sketch", "good.fa", "-o", "missing/out"];
        args.extend(force);
        t.failure(&args, "missing/out.aiips");
    }
    assert!(!t.path("missing").exists());
}

#[test]
fn sketch_directory_destination_preserved_and_staging_cleaned() {
    let t = Workspace::new();
    fs::create_dir(t.path("out.aiips")).unwrap();
    t.write("out.aiips/result", b"keep");
    t.failure(&["sketch", "good.fa", "-o", "out", "--force"], "out.aiips");
    assert_eq!(fs::read(t.path("out.aiips/result")).unwrap(), b"keep");
    t.assert_no_staging_files();
}

#[test]
fn sketch_rejects_identical_input_output() {
    let t = Workspace::new();
    let data = b">a\nACGT\n";
    t.write("input.aiips", data);
    t.failure(
        &["sketch", "input.aiips", "-o", "./input.aiips", "--force"],
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
        &["sketch", "good.fa", "-o", "alias.aiips", "-f"],
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
    t.failure(&["sketch", "good.fa", "-o", "out"], "--force");
    t.success(&["sketch", "good.fa", "-o", "out", "-f"]);
    assert!(!fs::symlink_metadata(t.path("out.aiips"))
        .unwrap()
        .file_type()
        .is_symlink());
    assert!(!t.path("missing").exists());
}

#[test]
fn sketch_double_dash_and_escaped_paths() {
    let t = Workspace::new();
    t.write("-input\t.fa", b">a\nACGT\n");
    let stdout = t.success(&["sketch", "-o", "out\nname", "--", "-input\t.fa"]);
    assert!(stdout.contains("-input\\t.fa\t1\t4"));
    assert!(stdout.ends_with("output\tout\\nname.aiips\n"));
    t.assert_empty("out\nname.aiips");
}

#[test]
fn sketch_help() {
    for flag in ["-h", "--help"] {
        assert!(Workspace::new()
            .success(&["sketch", flag])
            .contains("Usage: aiips sketch"));
    }
}

reject!(sketch_missing_input, ["sketch", "-o", "out"], "exactly 1");
reject!(
    sketch_extra_input,
    ["sketch", "good.fa", "good.fq", "-o", "out"],
    "exactly 1"
);
reject!(
    sketch_missing_output,
    ["sketch", "good.fa"],
    "--output is required"
);
reject!(
    sketch_missing_output_value,
    ["sketch", "good.fa", "-o"],
    "missing output path"
);
reject!(
    sketch_duplicate_output,
    ["sketch", "good.fa", "-o", "out", "-o", "other"],
    "only be specified once"
);
reject!(
    sketch_missing_file,
    ["sketch", "missing.fa", "-o", "out"],
    "missing.fa"
);
reject!(
    sketch_malformed_fastq,
    ["sketch", "bad.fq", "-o", "out"],
    "bad.fq"
);
reject!(
    sketch_empty_input,
    ["sketch", "empty.fa", "-o", "out"],
    "empty.fa"
);
reject!(
    sketch_zero_k,
    ["sketch", "good.fa", "-o", "out", "-k", "0"],
    "positive integer"
);
reject!(
    sketch_negative_k,
    ["sketch", "good.fa", "-o", "out", "-k", "-1"],
    "positive integer"
);
reject!(
    sketch_nonnumeric_k,
    ["sketch", "good.fa", "-o", "out", "-k", "abc"],
    "positive integer"
);
reject!(
    sketch_missing_k,
    ["sketch", "good.fa", "-o", "out", "-k"],
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
    ["sketch", "good.fa", "-o", "out", "--unknown"],
    "unknown option"
);
reject!(sketch_stdin, ["sketch", "-", "-o", "out"], "stdin");
reject!(
    sketch_stdout_destination,
    ["sketch", "good.fa", "-o", "-"],
    "filesystem path"
);
reject!(
    sketch_directory_input,
    ["sketch", ".", "-o", "out"],
    "regular file"
);

macro_rules! compare_pair {
    ($name:ident, $first:expr, $second:expr, $first_row:expr, $second_row:expr) => {
        #[test]
        fn $name() {
            let t = Workspace::new();
            let stdout = t.success(&["compare", $first, $second, "-k", "3"]);
            assert_eq!(
                stdout,
                format!(
                    "{STATS}\n{}\t{}\n{}\t{}\n\n{METRICS}\n{}\t{}\t-1\t-1\t-1\t-1\t-1\t-1\n",
                    $first, $first_row, $second, $second_row, $first, $second
                )
            );
        }
    };
}

compare_pair!(
    compare_fasta_fastq,
    "good.fa",
    "good.fq",
    "2\t13\t2\t3\t7",
    "2\t13\t2\t3\t7"
);
compare_pair!(
    compare_gzip_fastx,
    "good.fa.gz",
    "good.fq.gz",
    "2\t13\t2\t3\t7",
    "2\t13\t2\t3\t7"
);
compare_pair!(
    compare_two_sketches,
    "empty.aiips",
    "other.aiips",
    "NA\tNA\tNA\tNA\tNA",
    "NA\tNA\tNA\tNA\tNA"
);
compare_pair!(
    compare_fastx_sketch,
    "good.fa",
    "empty.aiips",
    "2\t13\t2\t3\t7",
    "NA\tNA\tNA\tNA\tNA"
);
compare_pair!(
    compare_sketch_fastx,
    "empty.aiips",
    "good.fq",
    "NA\tNA\tNA\tNA\tNA",
    "2\t13\t2\t3\t7"
);
compare_pair!(
    compare_same_file_twice,
    "good.fa",
    "good.fa",
    "2\t13\t2\t3\t7",
    "2\t13\t2\t3\t7"
);

#[test]
fn compare_default_k() {
    let stdout = Workspace::new().success(&["compare", "good.fa", "good.fq"]);
    assert!(stdout.contains("good.fa\t2\t13\t2\t21\t0\n"));
    assert!(stdout.contains("good.fq\t2\t13\t2\t21\t0\n"));
}

#[test]
fn compare_long_k_and_ambiguity_boundaries() {
    let t = Workspace::new();
    t.write("runs.fa", b">a\nac\nGTNacg\n>b\nTT\n>c\nAA\n");
    let stdout = t.success(&["compare", "runs.fa", "good.fq", "--kmer-size", "3"]);
    assert!(stdout.contains("runs.fa\t3\t12\t1\t3\t3\n"));
}

#[test]
fn compare_case_insensitive_sketch_extension() {
    let t = Workspace::new();
    t.write("empty.AiIpS", b"");
    let stdout = t.success(&["compare", "empty.AiIpS", "empty.aiips"]);
    assert!(stdout.contains("empty.AiIpS\tNA\tNA\tNA\tNA\tNA\n"));
}

#[test]
fn compare_detects_only_final_extension() {
    let t = Workspace::new();
    t.write("sequence.aiips.fa", b">a\nACGT\n");
    let stdout = t.success(&["compare", "sequence.aiips.fa", "good.fa", "-k", "3"]);
    assert!(stdout.contains("sequence.aiips.fa\t1\t4\t0\t3\t2\n"));
}

#[test]
fn compare_double_dash_and_escaped_paths() {
    let t = Workspace::new();
    let name = "-input\\\t\r\n.fa";
    t.write(name, b">a\nACGT\n");
    let stdout = t.success(&["compare", "--", name, "empty.aiips"]);
    assert!(stdout.contains("-input\\\\\\t\\r\\n.fa\t1\t4"));
    assert_eq!(stdout.lines().count(), 6);
    assert_eq!(stdout.lines().last().unwrap().split('\t').count(), 8);
}

#[test]
fn compare_help() {
    for flag in ["-h", "--help"] {
        assert!(Workspace::new()
            .success(&["compare", flag])
            .contains("Usage: aiips compare"));
    }
}

#[test]
fn compare_sketch_directory_rejected() {
    let t = Workspace::new();
    fs::create_dir(t.path("dir.aiips")).unwrap();
    t.failure(&["compare", "dir.aiips", "empty.aiips"], "regular file");
}

#[test]
fn compare_generated_placeholder_round_trip() {
    let t = Workspace::new();
    t.success(&["sketch", "good.fa", "-o", "generated"]);
    let stdout = t.success(&["compare", "generated.aiips", "good.fq"]);
    assert!(stdout.contains("generated.aiips\tNA\tNA\tNA\tNA\tNA\n"));
}

reject!(compare_missing_inputs, ["compare"], "exactly 2");
reject!(compare_one_input, ["compare", "good.fa"], "exactly 2");
reject!(
    compare_extra_input,
    ["compare", "good.fa", "good.fq", "empty.aiips"],
    "exactly 2"
);
reject!(
    compare_missing_first,
    ["compare", "missing.fa", "good.fa"],
    "missing.fa"
);
reject!(
    compare_missing_second,
    ["compare", "good.fa", "missing.fq"],
    "missing.fq"
);
reject!(
    compare_missing_sketch,
    ["compare", "empty.aiips", "missing.aiips"],
    "missing.aiips"
);
reject!(
    compare_malformed_first,
    ["compare", "bad.fq", "good.fa"],
    "bad.fq"
);
reject!(
    compare_malformed_second,
    ["compare", "good.fa", "bad.fq"],
    "bad.fq"
);
reject!(
    compare_empty_fastx,
    ["compare", "good.fa", "empty.fa"],
    "empty.fa"
);
reject!(
    compare_nonempty_sketch,
    ["compare", "empty.aiips", "bad.aiips"],
    "only empty"
);
reject!(
    compare_zero_k,
    ["compare", "good.fa", "good.fq", "-k", "0"],
    "positive integer"
);
reject!(
    compare_negative_k,
    ["compare", "good.fa", "good.fq", "-k", "-1"],
    "positive integer"
);
reject!(
    compare_nonnumeric_k,
    ["compare", "good.fa", "good.fq", "-k", "abc"],
    "positive integer"
);
reject!(
    compare_missing_k,
    ["compare", "good.fa", "good.fq", "-k"],
    "positive integer"
);
reject!(
    compare_overflow_k,
    [
        "compare",
        "good.fa",
        "good.fq",
        "-k",
        "999999999999999999999999999999"
    ],
    "positive integer"
);
reject!(
    compare_unknown_option,
    ["compare", "good.fa", "good.fq", "--unknown"],
    "unknown option"
);
reject!(
    compare_force_rejected,
    ["compare", "good.fa", "good.fq", "--force"],
    "unknown option"
);
reject!(
    compare_output_rejected,
    ["compare", "good.fa", "good.fq", "-o", "out"],
    "unknown option"
);
reject!(compare_stdin_first, ["compare", "-", "good.fa"], "stdin");
reject!(compare_stdin_second, ["compare", "good.fa", "-"], "stdin");
reject!(
    compare_directory_fastx,
    ["compare", "good.fa", "."],
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
    ["sketch", "good.fa", "-o", "--force"],
    "missing output path"
);

#[test]
fn sketch_dash_prefixed_output_path() {
    let t = Workspace::new();
    t.success(&["sketch", "good.fa", "-o", "./-out"]);
    t.assert_empty("-out.aiips");
}
