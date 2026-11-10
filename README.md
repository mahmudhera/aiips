[![CI](https://github.com/mahmudhera/aiips/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/mahmudhera/aiips/actions/workflows/ci.yml)
[![Coverage](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/mahmudhera/aiips/badges/coverage.json)](https://github.com/mahmudhera/aiips/actions/workflows/ci.yml?query=branch%3Amain)

## Summary

AIIPS means Anchor Indexed Island Probe Sketching. This project explores nested
sketching. Currently, the CLI builds classic MinHash sketches of FASTA/FASTQ files
and compares two sketches, estimating Jaccard, containment, and ANI. Other sketch
algorithms (fracminhash, bottomk, minimizer, open-syncmer, closed-syncmer,
strobemer) are reserved names, not yet implemented.

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

./target/release/aiips sketch data/phiX174_NC_001422.1.fasta -a minhash -o phix
./target/release/aiips sketch tests/fixtures/sample.fa -a minhash -k 3 -o sample
./target/release/aiips sketch tests/fixtures/sample.fq.gz -a minhash -k 3 -o reads
./target/release/aiips compare sample.aiips reads.aiips
```

The previous `aiips INPUT...` interface has been removed. A subcommand is required.
Stdin (`-`) is not supported. Use `--` before positional filenames beginning with
`-`; options must precede that separator. For output names beginning with `-`,
pass an explicit relative path, such as `-o ./-sample`.

### Create a sketch

```text
aiips sketch INPUT -o OUTPUT [-a ALGORITHM] [-k LENGTH] [-s NUM] [--seed SEED] [-f]
```

`-a`/`--algorithm` defaults to `minhash`, the only implemented algorithm, so
commands without it keep working.
Reserved names (fracminhash, bottomk, minimizer, open-syncmer, closed-syncmer,
strobemer) are refused as "not implemented yet"; other names are unknown. MinHash
options: `-k` from 1 to 32 (default 21), `-s`/`--num-hashes` from 2 to 1,000,000
(default 1000), and `--seed`, any unsigned 64-bit integer (default 42). Each option
may be given once. An input without any valid k-mer is refused and creates no
output.

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

The input is fully parsed and sketched in one streaming pass before any output is
created. Existing output paths are refused unless
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

### Sketch format

`.aiips` files are text. Line 1 names the format, version, sketch type, and every
parameter; line 2 records the input statistics; then come the `num_hashes` minima,
one decimal value per line, in hash-function order:

```text
aiips-sketch v1 algorithm=minhash k=3 num_hashes=1000 seed=42 hash=murmur3-fmix64 canonical=true
records=2 bases=13 non_acgt_bases=2 valid_kmers=7
<minimum of hash function 0>
...
```

Reading is strict: unknown, missing, or duplicate fields, a wrong value count,
non-decimal values, carriage returns, a missing final newline, or another version
are errors. Zero-byte files written by earlier placeholder versions are rejected
with a hint to re-run `aiips sketch`.

### MinHash

Classic MinHash (Broder 1997): each of `s` hash functions keeps its minimum over
the input's distinct canonical k-mers. A k-mer's canonical form is the smaller
2-bit encoding (A=0, C=1, G=2, T=3) of the k-mer and its reverse complement, so
sketches are strand-independent. Hash function `i` is
`fmix64(kmer XOR seed_i)`, MurmurHash3's 64-bit finalizer with `seed_i` the `i`-th
SplitMix64 output from `--seed`. fmix64 is a bijection, so distinct k-mers never
collide within one hash function, but that does not make the family min-wise
independent (Broder et al. 1998 distinguish exact from approximate min-wise
independence). All functions derive deterministically from one seed rather than
being independently sampled, and no approximation-error bound is established for
this family; it is a practical heuristic. Sketching costs `s` hash evaluations per k-mer and uses
O(`s`) memory; sequences are never stored.

### Compare two inputs

```text
aiips compare INPUT_A INPUT_B [-k LENGTH]
```

Each input is a sketch or a plain/gzip FASTA/FASTQ file; mixed pairs are
supported, so commands written before sketches existed keep working. Files named
`.aiips` (case-insensitive) must be sketches; other files are sketches when they
start with the sketch header and are otherwise read as sequence.

Sequence inputs are sketched in memory with MinHash. When the other input is a
MinHash sketch, its parameters are used and a given `-k` must equal its `k`.
Otherwise `k` comes from `-k` (default 21) with the default 1000 hash functions
and seed 42. A sequence input without valid k-mers is an error.

Unless every sketch input is MinHash, `compare` halts and names each input's type.
Two MinHash sketches must agree in `k`, `num_hashes`, `seed`, `hash`, and
`canonical`; every differing field is reported.

Output contains a statistics table, a blank line, and a
comparison table:

```text
input\trecords\tbases\tnon_acgt_bases\tk\tvalid_kmers
sample.aiips\t2\t13\t2\t3\t7
reads.aiips\t2\t13\t2\t3\t7

input_a\tinput_b\tjaccard\tcontainment_a_in_b\tcontainment_b_in_a\tani\taf_a_to_b\taf_b_to_a\tani_containment_a_in_b\tani_containment_b_in_a
sample.aiips\treads.aiips\t1.000000\t1.000000\t1.000000\t1.000000\t-1\t-1\t1.000000\t1.000000
```

Estimators, with `s` hash functions and minima `a_i`, `b_i`. Their stated
properties hold in the ideal model of independent, uniformly random hash
functions; the seeded fmix64 family above only approximates that model.

- `jaccard`: `J = #{i : a_i = b_i} / s` (Broder 1997).
- Distinct k-mer counts: `n = (s - 1) / Σ -ln(1 - a_i / 2^64)` (Cohen 1997). In the
  ideal continuous model each term is exponential with rate `n`, which makes the
  estimator unbiased; here it is a practical approximation.
- `containment_a_in_b`: `|A ∩ B| / |A|` with `|A ∩ B| = J / (1 + J) · (|A| + |B|)`
  (Koslicki & Zabeti 2019), clamped to [0, 1].
- `ani`: `1 - D`, with Mash distance `D = -(1/k) ln(2J / (1 + J))` and `D = 1` at
  `J = 0` (Ondov et al. 2016), clamped at 0.
- `ani_containment_a_in_b`: `C^(1/k)`. Hera, Pierce-Ward & Koslicki (2023) derive
  `E[C] = (1 - p)^k` (Theorem 5, Eq. 5) and the point estimate
  `p̂ = 1 - C^(1/k)`, so ANI = `1 - p̂` = `C^(1/k)` (see also Blanca et al. 2022).
  Their model assumes independent substitutions and excludes spurious k-mer
  matches and repeats. The paper's guarantees and confidence intervals are for its
  debiased FracMinHash containment; they do not transfer to this containment,
  which is derived from MinHash Jaccard and estimated counts. The root is
  nonlinear, so it does not preserve unbiasedness. Reporting 0 when `C = 0` is a
  convention, not evidence of zero biological identity.
- `af_a_to_b`, `af_b_to_a`: always `-1` (unavailable). Alignment fraction needs
  k-mer positions, which a global MinHash discards: tools reporting AF map
  positional seeds (FastANI via Mashmap, skani via chaining). Under the simple
  mutation model containment is about `AF · q^k`, so one containment value cannot
  separate AF from identity `q`.

Values are printed with six decimals. In the ideal model, estimates have sampling
error of roughly `sqrt(J(1 - J)/s)` for Jaccard and `1/sqrt(s)` relative error for
the counts. These figures omit any bias or dependence from the hash family. The
tests confirm accuracy on their fixed cases (known sets, phiX174 mutated at known
rates), not hash-family independence or an error bound.
`compare` does not accept `-o`, `--force`, or sketch options.

### Shared behavior

`-k`/`--kmer-size` may be given once per command; repeated short/long spellings
are errors.
`valid_kmers` counts every overlapping A/C/G/T occurrence in either case, including
repeats, rather than distinct k-mers. Other bases break valid runs. Windows never
cross record boundaries; wrapped FASTA lines belong to one record. FASTQ quality
lengths are checked but qualities do not filter bases.

Statistics and success output go to stdout. Errors identify the affected path
where applicable, go to stderr, and return status 1; success returns 0. All inputs
are validated before any success output is printed. Empty or malformed sequence
files, missing files, non-regular inputs, and invalid sketch files are errors. Filenames in TSV escape
backslash, tab, carriage return, and newline.

Command-specific help takes precedence over invalid arguments anywhere before
`--`. After `--`, `-h` and `--help` are literal filenames.

## Implementation

- `src/fastx.rs` streams one record at a time using
  [needletail](https://docs.rs/needletail/), also used by
  [skani](https://github.com/bluenote-1577/skani/blob/main/src/file_io.rs).
  Only gzip compression is enabled to keep dependencies small.
- `src/kmer.rs` iterates borrowed forward-strand windows and zero-based positions
  without per-k-mer allocation, and rolling canonical 2-bit k-mer codes for
  k <= 32. Both take O(sequence length) time and O(1) extra space.
- `src/hash.rs` provides MurmurHash3's fmix64 and SplitMix64.
- `src/sketch/` holds the algorithm registry and `.aiips` format (`mod.rs`) and
  MinHash construction and estimators (`minhash.rs`).
- The binary separates argument parsing, command execution and output handling,
  and TSV formatting. Memory scales with parser buffering and the largest record.

Sketch algorithms other than MinHash, AF estimation, k > 32, confidence intervals,
quality filtering, parallel processing, and compression other than gzip remain
unimplemented.

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
A separate Linux job checks formatting and Clippy. Tests that rely on Unix-only
behavior (symlinks, file identity, forced replacement, control characters in file
names) are compiled only on Unix; Windows runs equivalent checks of the documented
Windows behavior instead.

A Linux `coverage` job measures line coverage of `src/` with
[`cargo-llvm-cov`](https://github.com/taiki-e/cargo-llvm-cov) over the full test
suite, including the `aiips` binary launched by the integration tests. Code under
`tests/` is excluded, but unit-test modules inside `src/` are counted. Coverage is
reported, never enforced. Each run's summary page lists coverage by file. On `main`,
`scripts/coverage_badge.sh` turns the result into a badge file, which is pushed to
the `badges` branch for the README badge. To measure locally (requires
`cargo-llvm-cov` and the `llvm-tools` component):

```sh
cargo llvm-cov --locked --ignore-filename-regex '(^|/)tests/'
```

The library tests compare forward and canonical k-mer iteration with exhaustive
naive baselines, and streamed MinHash sketches with a brute-force baseline over the
distinct k-mer set. Estimators are checked against exact cases (identical,
disjoint, and reverse-complement inputs), hand-computed formula values, synthetic
sets with known Jaccard, containment, and cardinality (within 4 standard errors,
fixed seeds), and the phiX174 genome in `data/` mutated in memory at known
substitution rates (1% and 5%), where containment ANI recovers 1 - p within 0.008.
SplitMix64 is checked against its published reference outputs, fmix64 by
inversion.
Integration tests exercise each command, all input combinations, gzip, known
counts, output creation and preservation, errors, help, and removal of the old
interface. They use Cargo's test executable and isolated temporary directories;
the downloaded genome and existing results are not modified. Fixtures are
synthetic; gzip fixtures have fixed timestamps. Symlink checks run on Unix.
