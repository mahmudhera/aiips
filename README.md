## Summary

AIIPS means Anchor Indexed Island Probe Sketching. In this project, we explore if using nested sketching can lead to something interesting.

## Build and run

Create and activate the Conda development environment:

```sh
conda env create -f environment.yml
conda activate aiips
```

If the environment `aiips` already exists, use `conda env update -n aiips -f environment.yml`.
The environment includes Rust/Cargo, rustfmt, Clippy, and a C compiler/linker.
On macOS, Apple's command-line tools (including the macOS SDK) are also required;
install them with `xcode-select --install` if absent.
Cargo installs the Rust dependencies using the checked-in `Cargo.lock`:

```sh
cargo build --locked --release
./target/release/aiips -k 21 data/phiX174_NC_001422.1.fasta
./target/release/aiips -k 3 tests/fixtures/sample.fa
cat genome.fa | ./target/release/aiips -
```

The CLI accepts one or more FASTA/FASTQ files, auto-detected from their contents,
including gzip compression. `-` reads stdin; `--` ends option parsing.
The default k-mer length is 21; any positive `usize` length is supported.
Output is TSV with one summary per file:

```text
input                            records bases non_acgt_bases k valid_kmers
tests/fixtures/sample.fa         2       13    2              3 7
```

`valid_kmers` counts all overlapping occurrences, including repeats, rather than
distinct k-mers. A/C/G/T in either case are valid. Every other base breaks a
k-mer run; windows never cross record boundaries. Wrapped FASTA lines belong to
the same record. FASTQ qualities are parsed and their lengths checked but do not
filter bases. Malformed or empty inputs produce an error and a nonzero exit
status. In a multi-file run, earlier successful files may already have output.
Filenames in TSV escape backslash, tab, carriage return, and newline.

## Implementation

- `src/fastx.rs`: record-at-a-time parsing and summary aggregation using
  [needletail](https://docs.rs/needletail/), also used by
  [skani](https://github.com/bluenote-1577/skani/blob/main/src/file_io.rs).
  Only gzip support is enabled to keep dependencies small.
- `src/kmer.rs`: allocation-free iterator yielding borrowed forward-strand
  windows and zero-based positions. Original case is preserved. Iteration takes
  O(sequence length) time and O(1) extra space; no hashing or canonicalization yet.
- `src/main.rs`: CLI and TSV reporting. Memory scales with parser buffering and
  the largest record, rather than the entire input file.

Sketching, similarity/ANI estimation, quality filtering, and parallel processing
are not implemented yet. Compression other than gzip is not supported.

## Tests

```sh
cargo test
cargo fmt --check
cargo clippy --all-targets -- -D warnings
```

Tests cover parsing, gzip, stdin, CLI errors, known summary counts, record
boundaries, ambiguity, and k-mer iteration against an exhaustive naive baseline.
Fixtures are synthetic; gzip fixtures have a fixed timestamp.
