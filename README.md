[![CI](https://github.com/mahmudhera/aiips/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/mahmudhera/aiips/actions/workflows/ci.yml)

## Summary

AIIPS means Anchor Indexed Island Probe Sketching. This project explores nested
sketching. Currently, the CLI summarizes sequence files and implements placeholder
sketch creation and comparisons; no sketching or similarity estimation is performed.

## Build

```sh
conda env create -f environment.yml
conda activate aiips
cargo build --locked --release
```

If `aiips` already exists, use `conda env update -n aiips -f environment.yml`.
The environment includes Rust/Cargo, rustfmt, Clippy, and a C compiler/linker.
On macOS, Apple's command-line tools and SDK are also required
(`xcode-select --install` if absent). Cargo installs the Rust dependencies pinned
in `Cargo.lock`.

## Commands

```sh
./target/release/aiips --help
./target/release/aiips sketch --help
./target/release/aiips compare --help

./target/release/aiips sketch tests/fixtures/sample.fa -k 3 -o sample
./target/release/aiips compare tests/fixtures/sample.fa tests/fixtures/sample.fq.gz -k 3
./target/release/aiips compare sample.aiips tests/fixtures/sample.fa -k 3
```

The previous `aiips INPUT...` interface has been removed. A subcommand is required.
Stdin (`-`) is not supported. Use `--` before positional filenames beginning with
`-`; options must precede that separator. For output names beginning with `-`,
pass an explicit relative path, such as `-o ./-sample`.

### Create a placeholder sketch

```text
aiips sketch INPUT -o OUTPUT [-k LENGTH] [-f]
```

Accepts exactly one plain or gzip-compressed FASTA/FASTQ file, detected by content.
`-o`/`--output` is required. `.aiips` is appended unless that extension is already
present (case-insensitive): `sample` becomes `sample.aiips`, and `sample.txt`
becomes `sample.txt.aiips`. The actual destination is printed on success.
The dotfile `.aiips` is already a valid sketch name and is not suffixed again.
Directory destinations and directory-like spellings (`outdir/`, `.`, `..`) are
rejected, without suggesting that force can replace a directory.
Output existence, file type, input aliases, and the destination parent are checked
before scanning sequence content, then rechecked before installation. New output
uses atomic exclusive creation and does not require filesystem hard-link support.

The output is an **empty placeholder file**, containing no sequence information,
parameters, or actual sketch data. Existing output paths are refused unless
`-f`/`--force` is supplied. Force replacement is staged in the destination directory
after input validation; malformed input or replacement failure preserves the old
output. File identity checks on Unix reject input aliases, including case aliases,
hardlinks, and symlinks, even with force.
Destination directories must already exist. Forced replacement replaces the
output directory entry, including a symlink, rather than following it to truncate
its target.

Replacing a regular file preserves its permission bits. Replacement creates a new
inode: ownership, group, ACLs, extended attributes, and timestamps are not
preserved. Replacing a symlink does not copy metadata from its target. Forced
replacement of existing paths is supported on Unix; other platforms refuse it
because the dependency-free implementation cannot guarantee file identity there.

The old output is retained in a temporary backup until the stdout report has been
written and flushed successfully. A reported write or flush failure removes a new
output or restores the previous file/symlink, including its original metadata.
If restoration itself fails, the error reports the retained backup path for
recovery. Output may already be partially printed when stdout fails. This is
rollback for reported errors, not a crash-recovery protocol or a guarantee against
concurrent external changes to the same paths.
After a successful stdout flush, the new output is committed. A failure to remove
the old backup then leaves the new output in place, prints a warning on stderr
with the retained backup path, and keeps exit status 0.

Successful output is TSV statistics followed by the output path:

```text
input\trecords\tbases\tnon_acgt_bases\tk\tvalid_kmers
tests/fixtures/sample.fa\t2\t13\t2\t3\t7
output\tsample.aiips
```

Examples show tab characters as `\t` for readability.

### Compare two inputs

```text
aiips compare INPUT_A INPUT_B [-k LENGTH]
```

Each input whose final extension or whole filename is `.aiips` (case-insensitive) must be a readable,
zero-byte regular file. Nonempty sketches are rejected until a real format is
implemented. Other inputs are parsed as plain/gzip FASTA/FASTQ. Two sequence files,
two placeholders, or a mixed pair are supported.

Output contains a statistics table, a blank line, and a comparison table. For
placeholder inputs all five statistic fields, including `k`, are `NA`: an empty
file contains no metadata. All six metric values are `-1`, meaning unavailable,
**not measured similarity or distance**.

```text
input\trecords\tbases\tnon_acgt_bases\tk\tvalid_kmers
sample.aiips\tNA\tNA\tNA\tNA\tNA
tests/fixtures/sample.fa\t2\t13\t2\t3\t7

input_a\tinput_b\tjaccard\tcontainment_a_in_b\tcontainment_b_in_a\tani\taf_a_to_b\taf_b_to_a
sample.aiips\ttests/fixtures/sample.fa\t-1\t-1\t-1\t-1\t-1\t-1
```

Containment and AF (alignment fraction) have columns for both directions; Jaccard
and ANI each have one column. No statistical definitions or estimators are
implemented by these placeholders. `compare` does not accept `-o` or `--force`.

### Shared behavior

`-k`/`--kmer-size` defaults to 21 and accepts any positive length fitting `usize`.
Each command accepts this option once; repeated short/long spellings are errors.
It applies to sequence statistics; placeholders cannot recover or validate `k`.
`valid_kmers` counts every overlapping A/C/G/T occurrence in either case, including
repeats, rather than distinct k-mers. Other bases break valid runs. Windows never
cross record boundaries; wrapped FASTA lines belong to one record. FASTQ quality
lengths are checked but qualities do not filter bases.

Statistics and success output go to stdout. Errors identify the affected path
where applicable, go to stderr, and return status 1; success returns 0. All inputs
are validated before any success output is printed. Empty or malformed sequence
files, missing files, and non-regular inputs are errors. Filenames in TSV escape
backslash, tab, carriage return, and newline.

Command-specific help takes precedence over invalid arguments anywhere before
`--`. After `--`, `-h` and `--help` are literal filenames.

## Implementation

- `src/fastx.rs` streams one record at a time using
  [needletail](https://docs.rs/needletail/), also used by
  [skani](https://github.com/bluenote-1577/skani/blob/main/src/file_io.rs).
  Only gzip compression is enabled to keep dependencies small.
- `src/kmer.rs` iterates borrowed forward-strand windows and zero-based positions
  without per-k-mer allocation. Original case is preserved. Iteration takes
  O(sequence length) time and O(1) extra space.
- The binary separates argument parsing, command execution and output handling,
  and TSV formatting. Memory scales with parser buffering and the largest record.

Hashing, canonicalization, real sketch serialization, similarity/ANI estimation,
quality filtering, parallel processing, and compression other than gzip remain
unimplemented. The placeholder `.aiips` extension establishes no future binary
format.

## Tests

```sh
cargo test --locked
cargo fmt --check
cargo clippy --locked --all-targets -- -D warnings
cargo build --locked --release && scripts/smoke.sh target/release/aiips
```

`scripts/smoke.sh` runs every subcommand against the committed fixtures and checks
exit codes and exact stdout.

### Continuous integration

[`.github/workflows/ci.yml`](.github/workflows/ci.yml) runs on every push to `main`
(and on manual dispatch) with Rust 1.99.1. On Linux, macOS (arm64) and Windows, it
builds the release binary, runs `cargo test --locked`, and runs `scripts/smoke.sh`.
A separate Linux job checks formatting and Clippy. See [SPEC.md](SPEC.md).

The library tests compare k-mer iteration with an exhaustive naive baseline.
Integration tests exercise each command, all input combinations, gzip, known
counts, output creation and preservation, errors, help, and removal of the old
interface. They use Cargo's test executable and isolated temporary directories;
the downloaded genome and existing results are not modified. Fixtures are
synthetic; gzip fixtures have fixed timestamps. Symlink checks run on Unix.
