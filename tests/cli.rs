use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

fn command() -> Command {
    Command::new(env!("CARGO_BIN_EXE_aiips"))
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

#[test]
fn fasta_fastq_and_gzip_have_identical_known_counts() {
    let output = command()
        .args(["-k", "3"])
        .args(["sample.fa", "sample.fq", "sample.fa.gz", "sample.fq.gz"].map(fixture))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines[0],
        "input\trecords\tbases\tnon_acgt_bases\tk\tvalid_kmers"
    );
    assert_eq!(lines.len(), 5);
    for line in &lines[1..] {
        assert!(line.ends_with("\t2\t13\t2\t3\t7"), "{line}");
    }
}

#[test]
fn supports_stdin_and_default_k() {
    let mut child = command()
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(b">a\nACGT\n")
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success());
    assert!(String::from_utf8(output.stdout)
        .unwrap()
        .ends_with("-\t1\t4\t0\t21\t0\n"));
}

#[test]
fn invalid_arguments_fail_and_help_succeeds() {
    for args in [
        vec![],
        vec!["-k", "0"],
        vec!["-k", "-1"],
        vec!["-k"],
        vec!["-k", "abc"],
        vec!["--unknown"],
        vec!["-", "-"],
    ] {
        let output = command().args(args).output().unwrap();
        assert!(!output.status.success());
        assert!(!output.stderr.is_empty());
    }
    let output = command().arg("--help").output().unwrap();
    assert!(output.status.success());
    assert!(String::from_utf8(output.stdout).unwrap().contains("Usage:"));
}

#[test]
fn missing_file_names_input_in_error() {
    let output = command()
        .arg(fixture("does-not-exist.fa"))
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("does-not-exist.fa"));
}

#[test]
fn malformed_fastq_returns_error_without_summary_row() {
    let mut child = command()
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(b"@a\nACGT\n+\nIII\n")
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("aiips: -:"));
    assert_eq!(String::from_utf8(output.stdout).unwrap().lines().count(), 1);
}
