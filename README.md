[![CI](https://github.com/mahmudhera/aiips/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/mahmudhera/aiips/actions/workflows/ci.yml)
[![Coverage](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/mahmudhera/aiips/badges/coverage.json)](https://github.com/mahmudhera/aiips/actions/workflows/ci.yml?query=branch%3Amain)

## Summary

AIIPS means Anchor Indexed Island Probe Sketching. This project explores nested
sketching. Currently, the CLI builds classic MinHash, FracMinHash, or AIIPS
sketches of FASTA/FASTQ files and compares two sketches, estimating Jaccard,
containment, and ANI. AIIPS sketches store positioned anchor and probe k-mer
occurrences. Besides Jaccard and directional containment, an AIIPS comparison
estimates ANI and alignment fraction from positionally matched regions, without
base-level alignment; both are reported as unavailable when the evidence is
insufficient. Other sketch algorithms
(bottomk, minimizer, open-syncmer, closed-syncmer, strobemer) are reserved names,
not yet implemented.

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

./target/release/aiips sketch data/phiX174_NC_001422.1.fasta -a fracminhash --scaled 10 -o phix_frac
./target/release/aiips compare phix_frac.aiips data/phiX174_NC_001422.1.fasta

./target/release/aiips sketch data/phiX174_NC_001422.1.fasta -a aiips -o phix_aiips
./target/release/aiips compare phix_aiips.aiips data/phiX174_NC_001422.1.fasta
```

The previous `aiips INPUT...` interface has been removed. A subcommand is required.
Stdin (`-`) is not supported. Use `--` before positional filenames beginning with
`-`; options must precede that separator. For output names beginning with `-`,
pass an explicit relative path, such as `-o ./-sample`.

### Create a sketch

```text
aiips sketch INPUT -o OUTPUT [-a ALGORITHM] [-k LENGTH] [-s NUM] [--scaled S] [--seed SEED] [--probe-k K2] [--probe-scaled C2] [--anchor-k K1] [--anchor-scaled C1] [-f]
```

`-a`/`--algorithm` is `minhash` (the default, so commands without it keep
working), `fracminhash`, or `aiips`.
Reserved names (bottomk, minimizer, open-syncmer, closed-syncmer, strobemer) are
refused as "not implemented yet"; other names are unknown. Every algorithm takes
`--seed`, any unsigned 64-bit integer (default 42). MinHash and FracMinHash take
`-k` from 1 to 32 (default 21). MinHash also takes `-s`/`--num-hashes` from 2 to
1,000,000 (default 1000); FracMinHash takes `--scaled`, an integer scale factor
from 1 to 2^64 - 1 (default 1000). AIIPS takes `--probe-k` and `--anchor-k`,
odd k-mer lengths from 1 to 31 (defaults 15 and 21), and `--probe-scaled` and
`--anchor-scaled`, integer scale factors from 1 to 2^64 - 1 (default 125 each);
`--anchor-k` must exceed `--probe-k`. AIIPS records must be at most
4,294,967,295 bases long. Options that do not
apply to the chosen algorithm (`--scaled` with MinHash, `-s` with FracMinHash,
`-k`, `-s` or `--scaled` with AIIPS, the `--probe-*` and `--anchor-*` options
with the others) are errors. Each option may be given once. An input without any valid k-mer is
refused and creates no output. A FracMinHash or AIIPS sketch that keeps no hash
value is still written, with a warning on stderr, because it is a valid sample;
it cannot be compared.

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

Examples show tab characters as `\t` for readability. For AIIPS, `k` and
`valid_kmers` refer to the probe k-mers, and a table of build diagnostics
precedes the output path:

```text
input\trecords\tbases\tnon_acgt_bases\tk\tvalid_kmers
tests/fixtures/sample.fa\t2\t13\t2\t3\t7

input\tanchor_occurrences\tdistinct_anchors\tunique_anchors\tprobe_occurrences\tdistinct_probes\tambiguous_kmers\tfile_bytes
tests/fixtures/sample.fa\t1\t1\t1\t7\t3\t2\t324
output\tsample.aiips
```

(This example used `--anchor-k 3 --anchor-scaled 2 --probe-k 3 --probe-scaled 1`.)

- `anchor_occurrences` and `probe_occurrences` count every retained position,
  including repeats.
- `distinct_anchors` and `distinct_probes` count distinct retained hashes.
- `unique_anchors` counts anchor hashes that occur exactly once.
- `ambiguous_kmers` counts the probe-length windows inside records that contain
  a base other than A/C/G/T.
- `file_bytes` is the size of the written sketch.

### Sketch format

`.aiips` files are text. Line 1 names the format, version, sketch type, and every
parameter; line 2 records the input statistics; then come the hash values, one
decimal value per line. A MinHash sketch holds its `num_hashes` minima in
hash-function order:

```text
aiips-sketch v1 algorithm=minhash k=3 num_hashes=1000 seed=42 hash=murmur3-fmix64 canonical=true
records=2 bases=13 non_acgt_bases=2 valid_kmers=7
<minimum of hash function 0>
...
```

A FracMinHash sketch holds its `size` retained hash values in strictly increasing
order. `size` is the number of values, not a comparison parameter:

```text
aiips-sketch v1 algorithm=fracminhash k=3 scaled=1 seed=42 hash=murmur3-fmix64 canonical=true size=3
records=2 bases=13 non_acgt_bases=2 valid_kmers=7
<smallest retained hash value>
...
```

An AIIPS sketch keeps the two text lines (its `valid_kmers` counts probe
k-mers), followed by a binary body of exactly `body_bytes` bytes:

```text
aiips-sketch v1 algorithm=aiips anchor_k=21 anchor_scaled=125 probe_k=15 probe_scaled=125 seed=42 hash=murmur3-fmix64 canonical=true body_bytes=N
records=2 bases=13 non_acgt_bases=2 valid_kmers=7
<binary body>
```

The body is a sequence of tagged sections, in this order, followed by a
checksum. Fixed-width integers are little-endian. Every other integer is an
unsigned LEB128 varint (7 bits per byte, least significant group first) in its
shortest encoding. Bit columns are packed least significant bit first, with
zero padding:

```text
section  := tag (4 ASCII bytes) | payload length (u64) | payload
GPRB     := n | h_1 | h_2 - h_1 | ... | h_n - h_(n-1)    (distinct probe hashes, increasing)
GANC     := the same for the distinct anchor hashes
CTGS     := n | per contig: length, anchor count, probe count
ANCH     := per anchor: position - previous anchor's position on its contig (first: position)
            | per anchor: index into GANC | strand bits
PROB     := as ANCH, indices into GPRB | strand bits | role bits
checksum := FNV-1a 64 of all preceding body bytes (u64)
```

Occurrences are listed by contig, then position. Each occurrence stores an
index into its dictionary (GPRB or GANC) rather than a hash, and the dictionary
must contain exactly the hashes the occurrences use. Contigs are recorded by
length only, in input order.

Island links (which probes each anchor owns, and the next anchor) are not
stored: they follow from the positions and are rebuilt on loading.

The checksum detects corruption that would otherwise decode into plausible
values. For a 5 Mbp genome at the defaults, the file is 839 kB and holds about
80,000 positioned occurrences.

Reading is strict: unknown, missing, or duplicate fields, a wrong value count,
non-decimal values, FracMinHash values out of order, repeated, or above the
`scaled` threshold, carriage returns, a missing final newline, or another version
are errors. In AIIPS files these are also errors:
- missing, reordered, or extra sections, or a payload that is not consumed exactly;
- a body length other than `body_bytes`, or a checksum mismatch;
- non-minimal or overflowing varints;
- dictionary hashes that are not strictly increasing, or probe hashes above
  their threshold (anchor identities are full 64-bit hashes);
- dictionary indices out of range, or dictionary entries that no occurrence uses;
- positions that do not increase within a contig, or k-mers that would run past the contig end;
- non-zero padding bits;
- a contig table that disagrees with `records` or `bases`;
- more probe occurrences than `valid_kmers`, or more `valid_kmers` than windows;
- contigs longer than 4,294,967,295 bases;
- an even or out-of-range `anchor_k` or `probe_k`, or `anchor_k` not above `probe_k`.

A FracMinHash file may not exceed 21 bytes per declared value plus 4096 bytes,
and an AIIPS file its declared `body_bytes` plus 4096 bytes; every other sketch
file is bounded by the largest MinHash sketch. Zero-byte files written by earlier placeholder versions are rejected
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

### FracMinHash

FracMinHash (Irber et al. 2022; Hera, Pierce-Ward & Koslicki 2023) keeps every
distinct canonical k-mer whose hash is at most `max_hash = floor((2^64 - 1) / S)`
for scale factor `S` (`--scaled`), the convention used by sourmash. The hash is
`fmix64(kmer XOR seed_0)`, MinHash's hash function 0 for the same `--seed`, so the
hash-family caveats above apply. In the ideal model of uniformly random hash
values a k-mer is kept with probability `(max_hash + 1) / 2^64`. This is exactly
`1/S` when `S` is a power of two and otherwise within `2^-64` of `1/S`. `S = 1`
keeps every distinct k-mer, so estimates are then exact set statistics. In that
model the sketch size `|FRAC(A)|` is Binomial(`|A|`, about `1/S`), so
`S · |FRAC(A)|` estimates the number of distinct k-mers `|A|`, with standard
deviation `sqrt(|A|(S - 1))`.

Sketching costs one hash evaluation per k-mer. Retained hashes are buffered and
compacted (sorted, deduplicated) whenever the buffer doubles, so memory is
proportional to the number of distinct retained hashes, about `|A| / S`;
sequences are never stored.

### AIIPS

AIIPS samples canonical k-mers of two kinds, genome-wide. `x` is the canonical
2-bit code, `seed_i` is SplitMix64 output `i` from `--seed`, and
`T_i = floor((2^64 - 1) / c_i)`:

- *Probes* (k2-mers, `--probe-k`) measure k-mer overlap and local similarity.
  An occurrence is kept when
  `h2(x) = fmix64(x XOR seed_0) <= floor((2^64 - 1) / c2)`, where c2 is
  `--probe-scaled`. This is the FracMinHash rule and hash above, so the
  distinct probe hashes (the *global probe set*) are exactly the FracMinHash
  sketch with `-k` k2, `--scaled` c2, and the same seed. The tests check this
  bit for bit.
- Each probe has a *role* bit, `fmix64(x XOR seed_2) & 1`: 0 marks a mapping
  probe, 1 a scoring probe. The bit depends only on the k-mer, so a shared
  probe has the same role in both genomes. Which value means which role is
  this implementation's convention.
- *Anchors* (k1-mers, `--anchor-k`, which must exceed k2) locate corresponding
  regions of two genomes. A k1-mer is kept when
  `NOT(h2(u) XOR h2(v)) <= T1`, where `u` and `v` are its first and last
  k2-mers and c1 is `--anchor-scaled`: the XOR falls in the top `1/c1` band.
  Kept anchors are identified by `h1(x) = fmix64(x XOR seed_1)`.
  - Every position is already hashed with `h2` for the probes, so selecting
    anchors costs no further per-position hashing; `h1` is computed only for
    the anchors kept.
  - Under ideal hashing, the XOR of two independent uniform values is uniform
    and independent of each of them. Anchors are therefore kept at rate `1/c1`,
    independently of probe selection.
  - Reverse complementing swaps `u` and `v`, so the rule is strand-independent.
  - A k1-mer whose two k2-mers are equal (a short tandem repeat) is never kept.

  An independent hash of every k1-mer would sample the same way, but it was
  measured to cost about a third of the comparison time (see *Benchmark
  against skani*).

Every retained occurrence is stored with its contig, its zero-based start
position, and its strand: 1 when the reverse complement is the canonical form.
Probes that overlap anchors are kept as well.

Both k-mer lengths must be odd, so that no k-mer equals its reverse complement
and every strand is defined. A k-mer of at most 31 bases has a 2-bit code that
fits in a u64, and fmix64 is a bijection, so distinct k-mers never share a
hash: these sketches have no hash collisions.

*Islands.* On each contig, an anchor owns the probes that start at or after its
position and before the next anchor. The contig's last anchor owns the probes
up to the contig end, its island being closed by an END sentinel. Probes before
the first anchor form the contig's START bucket. The probes between any two
anchors on one contig, across all the islands in between, form one contiguous
slice. Ownership is computed from the positions by a single linear merge.

Sketching makes one pass per record:
- each base is decoded once, and only the k1 codes are rolled (the k2 window
  ending at the same base is a bit slice of them);
- every position gets one `h2` hash;
- the anchor rule reuses the two `h2` values.

Memory is proportional to the retained occurrences: a 32-bit position, a 32-bit
dictionary index and a strand bit per occurrence, plus the distinct hashes.
Sequences are never stored. At the defaults, a 5 Mbp genome takes 18 ms and
18 MB (an 873 kB file); a 100 Mbp record takes 0.3 s and 220 MB, mostly the
record itself.

### Compare two inputs

```text
aiips compare INPUT_A INPUT_B [-a ALGORITHM] [-k LENGTH] [-s NUM] [--scaled S] [--seed SEED] [--probe-k K2] [--probe-scaled C2] [--anchor-k K1] [--anchor-scaled C1] [aiips thresholds]
```

The aiips thresholds (`--min-block-span`, `--max-anchor-gap`, `--max-gap-ratio`,
`--min-mapping-hits`, `--min-order-fraction`, `--min-scoring-probes`) apply only
to AIIPS comparisons and are errors with other algorithms. They are described
under AIIPS estimators below.

Each input is a sketch or a plain/gzip FASTA/FASTQ file; mixed pairs are
supported, so commands written before sketches existed keep working. Files named
`.aiips` (case-insensitive) must be sketches; other files are sketches when they
start with the sketch header and are otherwise read as sequence.

Sequence inputs are sketched in memory. When an input is a sketch, its algorithm
and parameters are used, and every given option must equal the recorded value
(`-a` its algorithm, `-k` its `k`, `-s` its `num_hashes`, `--scaled` its
`scaled`, `--probe-k` its `probe_k`, `--probe-scaled` its `probe_scaled`,
`--anchor-k` its `anchor_k`, `--anchor-scaled` its `anchor_scaled`, `--seed`
its `seed`); options that do not apply to its algorithm are errors. These checks also apply when both inputs are sketches. When both inputs
are sequences, the algorithm comes from `-a` (default `minhash`) and the
parameters from the options, with the `sketch` defaults. A sequence input without
valid k-mers is an error.

Unless both inputs are MinHash, both FracMinHash, or both AIIPS, `compare` halts
and names each input's type. Two MinHash sketches must agree in `k`, `num_hashes`,
`seed`, `hash`, and `canonical`; two FracMinHash sketches in `k`, `scaled`, `seed`,
`hash`, and `canonical`; two AIIPS sketches in `anchor_k`, `anchor_scaled`,
`probe_k`, `probe_scaled`, `seed`, `hash`, and `canonical`, even though Jaccard
and containment use only the probes. Every differing field is reported. An empty FracMinHash
or AIIPS sketch halts the comparison, because its containment is undefined.

Output contains a statistics table, a blank line, and a
comparison table:

```text
input\trecords\tbases\tnon_acgt_bases\tk\tvalid_kmers
sample.aiips\t2\t13\t2\t3\t7
reads.aiips\t2\t13\t2\t3\t7

input_a\tinput_b\tjaccard\tcontainment_a_in_b\tcontainment_b_in_a\tani\taf_a_to_b\taf_b_to_a\tani_containment_a_in_b\tani_containment_b_in_a
sample.aiips\treads.aiips\t1.000000\t1.000000\t1.000000\t1.000000\t-1\t-1\t1.000000\t1.000000
```

MinHash estimators, with `s` hash functions and minima `a_i`, `b_i`. Their stated
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

FracMinHash estimators, with retained sets `FRAC(A)` and `FRAC(B)`:

- `jaccard`: `|FRAC(A) ∩ FRAC(B)| / |FRAC(A) ∪ FRAC(B)|`. The union of the
  retained sets is the retained set of the union, so this is a ratio of sampled
  counts.
- `containment_a_in_b`: `|FRAC(A) ∩ FRAC(B)| / |FRAC(A)|`, the raw ratio,
  reported only when both sketches are non-empty. Its bias depends on how empty
  sketches are treated. In the ideal model where each distinct k-mer is kept
  independently with probability `q ≈ 1/S`:
  - Hera, Pierce-Ward & Koslicki (2023) score an empty `FRAC(A)` as 0. The
    expectation is then `C · (1 - (1 - q)^|A|)`, and they divide by that factor to
    remove the bias. aiips does neither.
  - Given a non-empty `FRAC(A)`, the raw ratio is unbiased: given its size `m`,
    the kept part of A is a uniform `m`-subset, so the shared count has mean
    `m · C`.
  - aiips also refuses an empty `FRAC(B)`. An empty `FRAC(B)` would give ratio 0,
    so the reported value has expectation `C / P(FRAC(B) ≠ ∅ | FRAC(A) ≠ ∅)`. That
    is at least `C` and at most `C / (1 - (1 - q)^|B|)`, because the two
    non-emptiness events are positively correlated (Harris inequality). So the reported
    containment is biased upward, not downward. For example, with A = {x, y},
    B = {x, z}, `q = 1/2` and true containment 0.5, the three policies average
    0.375, 0.5 and 0.6.

  The last two statements are derived here, not taken from the paper. The tests
  check all three by exhaustive enumeration. The bias is negligible unless an
  expected sketch size (`|A|/S` or `|B|/S`) is small. Real hash functions only
  approximate the model.
- `ani` and `ani_containment_*`: the MinHash transforms above, Mash ANI from
  `jaccard` and `C^(1/k)` from each containment. The containment-ANI model and
  caveats above apply; the paper's guarantees and confidence intervals are for its
  debiased containment, not this raw ratio.
- `af_a_to_b`, `af_b_to_a`: always `-1`, as for MinHash: FracMinHash keeps no
  positions.

AIIPS Jaccard and containment are computed from the distinct global probe
hashes, never from occurrence counts. Let `a` and `b` be the sizes of the two
retained sets and `x` the size of their intersection:

- `jaccard`, `containment_a_in_b`, `containment_b_in_a`: `x / (a + b - x)`,
  `x / a`, and `x / b`. These are the FracMinHash ratios above, with the same
  bias caveats, and equal FracMinHash's values for the same k, scale, and seed.
- `ani_containment_a_in_b`, `ani_containment_b_in_a`: `C^(1/k2)`, the
  containment-ANI transform above. The same model and caveats apply: it is a
  screening value, invalid for incomplete or rearranged genomes.

AIIPS `ani`, `af_a_to_b` and `af_b_to_a` are positional. They are never filled
with a k-mer-set transform such as Mash ANI, and are `-1` when unsupported. The
steps involve no alignment and no dynamic programming; each threshold has a
flag (default in brackets):

1. *Shared anchors.* An anchor hash occurring exactly once in each genome
   is a shared anchor; repeated hashes are ignored. The relative orientation is
   strand_A XOR strand_B. For reverse matches, B positions become
   `len_B - (pos_B + k1)`, so collinear matches increase in both genomes.
2. *Runs.* Matches are grouped by contig pair and orientation, then sorted by
   A position. Consecutive matches stay in one run while all of these hold:
   - the B coordinate increases;
   - both gaps are at most `--max-anchor-gap` [50000];
   - the larger gap is at most `--max-gap-ratio` [2] times the smaller.
3. *Candidate blocks.* A run is cut wherever its A span reaches
   `--min-block-span` [10000] bases. Blocks next to each other share their
   boundary anchor. A shorter leftover joins the run's last block. A run that
   never reaches the span is still one candidate, because validation decides.
   Where blocks of different runs overlap on either genome, the blocks with
   more shared anchors are kept.
4. *Validation.* A block's probes are those lying entirely within
   `[first anchor, last anchor + k1)` and overlapping neither boundary anchor.
   Those probes match by construction, so they carry no evidence. Mapping-role
   probes whose hash occurs once in each block give hits. A block is accepted
   with at least `--min-mapping-hits` [3, at least 2] hits, when at least
   `--min-order-fraction` [0.8] of adjacent hits keep their order. This is a
   provisional heuristic, not a calibrated test.
5. *ANI.* The scoring-role probes of accepted blocks are pooled; any hash
   occurring more than once in either block is excluded from both. Then
   `q = 2X / (N_A + N_B)` is retained as the survival diagnostic. Mapping-only
   contrasts between up to four neighboring blocks on each side identify
   heterogeneous regions using a binomial working likelihood with a `2 ln N`
   penalty and at least 40 mapping-union hashes per side. These are heuristics,
   not calibrated significance thresholds: overlapping k-mers are dependent.
   Scoring counts are pooled within inferred regions, then regional roots are
   averaged by span. Uniform comparisons retain `ani = q^(1/k2)` to avoid
   small-window transformation bias. If any inferred region has no shared
   scoring probe or fewer than `--min-scoring-probes` [50] eligible probes, ANI
   falls back to the pooled root. Under the simple substitution
   model a homologous k2-mer is unchanged with probability `ANI^k2`.
   `ani` is `-1` in three cases: no block is accepted; fewer than
   `--min-scoring-probes` [50] probes are eligible; or none is shared. It is
   never reported as 0.
6. *AF.* The union of accepted block intervals on each genome, divided by
   its total length, ambiguous bases included. It is a positional coverage
   proxy, not a base-level alignment fraction:
   - ends outside the outermost shared anchors are not counted;
   - an unshared segment inside a block is counted.

   For fragmented inputs (either genome has more than one record), reported AF
   can be increased toward clipped `C/q` when all of these hold:
   - the mapping contrast finds no heterogeneity;
   - in each genome, distinct anchor hashes occurring more than once number
     under 2% of those occurring once;
   - scoring probes excluded as repeats are under 2% of eligible ones.

   Run breaks (order, gap, ratio) do not block it: islands and rearrangements
   break runs, and `C/q` estimates the homologous fraction regardless. The
   anchor repeat test is per genome because `repeated_anchors / shared_anchors`
   grows with divergence: a repeated hash is shared if any copy survives, a
   unique one only if its single copy does. At 100% identity the two ratios are
   equal. The correction never decreases observed span coverage. These guards
   are provisional and do not prove uniform homology. The `covered_bases_*`
   diagnostics always retain raw observed coverage, so reported AF may exceed
   their ratio to genome length.

These defaults favour divergent genomes; they are prototype values, not tuned.

**Why only the boundary anchors.** Excluding scoring probes that overlap *any*
shared anchor of a block biases ANI downward: the remaining probes then
include those overlapping anchors that substitutions destroyed, and these are
depleted of matches. Excluding only the two boundary anchors avoids this. Mean error (estimate minus truth) over three simulated
1 Mbp genomes with substitutions only, at the defaults:

| ANI | 99% | 95% | 90% | 85% |
| --- | --- | --- | --- | --- |
| any shared anchor excluded | -0.0015 | -0.0035 | -0.0039 | -0.0061 |
| boundary anchors excluded (used) | +0.0001 | +0.0002 | -0.0009 | -0.0043 |

At about 85% ANI and below, many 10 kb blocks lack three exact mapping hits.
Their rejection reduces AF (0.81-0.89 at 85% in these simulations) and selects
better-conserved blocks. The positional ANI has not been compared with
base-level alignments; treat it as uncalibrated.

`af_containment_a` and `af_containment_b` are `C / q`. For heterogeneous identity,
`q` differs from the regional mean ANI raised to k2. Under a
model where a fraction of the k2-mers is homologous and each survives with
probability ANI^k2, this estimates that fraction. It is a k-mer homologous
fraction, not an aligned-base fraction. It is not clipped, and can exceed 1
through repeats, incompleteness, or sampling noise; a value above 1 also
prints a warning on stderr. Reported AF uses its clipped value only under the
fragmentation guards above.

A diagnostics table follows the comparison row:

- `global_a`, `global_b`, `global_shared`, `max_containment` (`x / min(a, b)`),
  `putative_ani` (`max_containment^(1/k2)`): global set statistics.
- `shared_anchors`, `repeated_anchors`: matched and ignored anchor hashes.
- `breaks_order`, `breaks_gap`, `breaks_ratio`: runs broken for each reason.
- `candidate_blocks`, `overlapping_blocks`: blocks cut, and blocks dropped for
  overlap.
- `rejected_few_hits`, `rejected_disorder`, `accepted_blocks`, `mapping_hits`:
  validation results.
- `scoring_a`, `scoring_b`, `scoring_shared`: `N_A`, `N_B`, `X`.
- `excluded_boundary_overlap`, `excluded_repeats`: excluded scoring probes.
- `covered_bases_a`, `covered_bases_b`: AF numerators, 0 when no block is
  accepted.
- `af_containment_a`, `af_containment_b`: the cross-check above.
- `jaccard_status`, `containment_status`, `ani_status`, `af_status`:
  - `estimated`;
  - `insufficient_evidence` (blocks accepted, but scoring evidence too thin);
  - `no_supported_mapping` (no block accepted).

Comparing two 5 Mbp sketches takes about 9 ms and 7 MB.

Values are printed with six decimals. In the ideal model, MinHash estimates have
sampling error of roughly `sqrt(J(1 - J)/s)` for Jaccard and `1/sqrt(s)` relative
error for the counts, and FracMinHash containment roughly `sqrt(C(1 - C)/(|A|/S))`.
These figures are approximations that omit any bias or dependence from the hash
family. The tests confirm accuracy on their fixed cases (known sets, phiX174
mutated at known rates), not hash-family independence or an error bound.
`compare` does not accept `-o` or `--force`.

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

## Benchmark against skani

The original figures below precede the regional ANI and fragmented-AF
corrections. A rerun on the same 252 pairs reduces mosaic ANI MAE from 2.19 to
0.70 pp and draft AF MAE from 8.82 to 0.14 pp. Uniform-substitution ANI MAE
remains 0.10 pp. Both-sided-island ANI MAE rises from 0.16 to 0.19 pp.
In a fresh interleaved 84-pair cost run (five repeats, seed 101, both backbones),
the revised version takes 1.021x skani's time and 0.925x its peak memory;
the original binary takes 1.019x and 0.899x in the same run (per pair, the
revised/original time ratio has median 1.001). It remains slightly over skani's
time budget. These evaluation pairs informed development, so this rerun is not
an independent validation of the new heuristics. Mosaic ANI is still
overestimated: regions shorter than about 40 kb are not detected, and
low-identity blocks are likely rejected more often, so accepted blocks
over-represent high identity (not measured directly). Other limitations are the heuristic segmentation and the k-mer
assumptions behind AF.

These figures come from simulated genome pairs whose true ANI and AF are known
exactly. They compare AIIPS with skani 0.3.1 on a matched compute budget. They
are simulation results, not a calibration on real genomes.

**Data.** Genome A is a backbone, and B is derived from it by controlled
changes:
- *Backbones:* *E. coli* K-12 MG1655 (NC_000913.3, 4.64 Mb, downloaded from
  NCBI), with real repeats and base composition, and, as a control, i.i.d.
  random sequence of the same length.
- *subs:* substitutions only, ANI from 99.9% to 80%.
- *af_one:* 20 kb islands deleted from B (10-70% of A), at 99%, 95% and 90% ANI.
- *af_both:* islands deleted from B and novel 20 kb islands inserted into B
  (10-50% each).
- *indel:* substitutions plus small indels, at a tenth of the substitution
  rate, mean length 3.
- *mosaic:* alternating 50 kb regions of different identity (99/90% and 98/85%).
- *draft:* B fragmented into contigs with mean length 20 kb or 5 kb.
- *rearranged:* five inversions and five translocations of 50-200 kb.

Truth is tracked from the simulated events:
- ANI is the identity over aligned bases (gap columns excluded);
- AF of each genome is its aligned bases over its length.

The simulator checks its own bookkeeping against brute-force comparisons
(`simulate.py --self-test`).

**Protocol.**
- Both tools run on one thread, starting from FASTA, so sketching is included:
  `skani dist -t 1` and `aiips compare -a aiips`.
- AIIPS was tuned on separate tuning pairs (seed 1). The results below come
  from 252 held-out pairs (seeds 101-103, both backbones, every scenario).
- skani runs as shipped, with its learned ANI correction, and also with
  `--no-learned-ani`, which turned out to be more accurate on these
  simulations.
- Cost is measured with the tools run back to back on each pair (5 repeats,
  medians, 84 pairs).

**Cost.** AIIPS at its defaults takes 1.02x skani's wall time (median per
pair; skani's median is 28 ms) and 0.95x its peak memory.

To reach this, AIIPS uses a single fused pass per record and selects anchors
from the probe hashes; an independent anchor hash cost an extra ~9 ms per
pair, about a third more time. Its dictionaries are built with a radix sort.

**Accuracy.** Mean absolute error in percentage points; AIIPS at its defaults:

| Family | ANI: AIIPS | ANI: skani | ANI: skani, no learned correction | AF: AIIPS | AF: skani |
| --- | --- | --- | --- | --- | --- |
| subs (60 pairs) | **0.10** | 0.58 | 0.43 | **9.9** | 21.5 |
| af_one (72) | **0.06** | 0.30 | 0.11 | **4.5** | 5.2 |
| af_both (54) | **0.16** | 1.54 | 1.38 | **2.8** | 4.6 |
| indel (18) | **0.50** | 0.83 | 0.63 | **0.4** | 6.5 |
| rearranged (12) | 0.02 | 0.35 | 0.02 | **0.3** | 1.0 |
| draft (24) | 0.03 | 0.48 | 0.04 | 8.8 | **2.9** |
| mosaic (12) | 2.19 | **0.20** | 0.25 | **3.4** | 17.0 |

AIIPS reports every pair. skani reports nothing for 10% of the substitution
pairs, those near 80% ANI, below its screening threshold. Both backbones give
the same picture.

Within the matched budget, AIIPS estimates ANI more accurately than either
skani mode for:
- substitution-only genomes;
- genomes with deleted or novel islands;
- small indels.

Against skani without the learned correction, AIIPS is better on ANI in 49 of
72 one-sided, 39 of 54 both-sided and 39 of 60 substitution pairs.

AIIPS's AF is more accurate for substitution-only, both-sided-island, indel,
mosaic and rearranged genomes. The one-sided islands are mixed; AF error in
pp, AIIPS vs skani:

| ANI | 99% | 95% | 90% |
| --- | --- | --- | --- |
| AIIPS | 1.25 | 3.02 | **9.1** |
| skani | **0.39** | **2.37** | 12.9 |

The mean over all one-sided pairs favours AIIPS only because of the 90% pairs.

skani is better in three cases:
- *Mosaic genomes:* AIIPS pools k-mer survival across regions of different
  identity, and q^(1/k) of a mean survival overstates mean identity. Windowed
  estimates help only with windows small enough to bias homogeneous genomes.
- *AF of fragmented drafts:* AIIPS does not count each contig's ends beyond its
  outermost shared anchors.
- *AF of one-sided islands at high identity:* each island boundary loses the
  span to the nearest shared anchor.

Denser probes (`--probe-scaled 80`) cut the AF error to 4.5 (subs), 2.3
(af_both) and 0.8 (mosaic) at 1.10x skani's time, beyond the matched budget.

**Reproduce** (skani on `PATH`; data and results go to the git-ignored
`data/test_against_skani`):

```sh
curl -o data/test_against_skani/backbones/ecoli_K12_MG1655_NC_000913.3.fasta --create-dirs \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text"
python3 scripts/skani_benchmark/simulate.py --self-test
python3 scripts/skani_benchmark/simulate.py --seeds 101 102 103
python3 scripts/skani_benchmark/run.py --tool skani --label skani --seeds 101 102 103 --out data/test_against_skani/results/evaluation.tsv
python3 scripts/skani_benchmark/run.py --tool skani --label skani_no_learned "--options=--no-learned-ani" --seeds 101 102 103 --out data/test_against_skani/results/evaluation.tsv
python3 scripts/skani_benchmark/run.py --tool aiips --label aiips --seeds 101 102 103 --out data/test_against_skani/results/evaluation.tsv
python3 scripts/skani_benchmark/summarize.py data/test_against_skani/results/evaluation.tsv --by family_all --versus aiips skani_no_learned
cd scripts/skani_benchmark && python3 cost.py --pairs ../../data/test_against_skani/pairs --aiips ../../target/release/aiips \
  --seeds 101 --config skani=skani --config aiips=aiips --out ../../data/test_against_skani/results/cost.tsv
```

## Implementation

- `src/fastx.rs` streams one record at a time using
  [needletail](https://docs.rs/needletail/), also used by
  [skani](https://github.com/bluenote-1577/skani/blob/main/src/file_io.rs).
  Only gzip compression is enabled to keep dependencies small.
- `src/kmer.rs` iterates borrowed forward-strand windows and zero-based positions
  without per-k-mer allocation, and rolling canonical 2-bit k-mer codes for
  k <= 32, optionally with each window's position and strand. All take
  O(sequence length) time and O(1) extra space.
- `src/hash.rs` provides MurmurHash3's fmix64 and SplitMix64.
- `src/sketch/` holds the algorithm registry and `.aiips` format (`mod.rs`),
  MinHash construction and estimators (`minhash.rs`), FracMinHash construction
  and estimators (`fracminhash.rs`), the binary body encoding (`binary.rs`:
  varints, sections, FNV-1a 64), and AIIPS in `aiips/`: sampling, island
  navigation and estimators (`mod.rs`), dictionary-encoded occurrence columns
  and the contig table (`occurrences.rs`), the file layout (`format.rs`),
  shared anchors and candidate blocks (`blocks.rs`), and validation, ANI and AF
  (`local.rs`). FracMinHash and AIIPS comparisons merge two sorted sketches in
  linear time; shared anchors are found by merging the two sorted anchor
  dictionaries.
- The binary separates argument parsing, command execution and output handling,
  and TSV formatting. Memory scales with parser buffering and the largest record.
- `scripts/skani_benchmark/` simulates genome pairs with known ANI and AF
  (`simulate.py`), runs skani and AIIPS on them (`run.py`), and summarizes
  accuracy (`summarize.py`) and interleaved cost (`cost.py`); Python standard
  library only.

The following remain unimplemented:
- sketch algorithms other than MinHash, FracMinHash, and AIIPS;
- in AIIPS comparison: diagonal-consistency checks, sliding-window segmentation
  of blocks, matching of repeated anchors, per-stage timings, and calibration
  of ANI and AF against base-level alignments;
- debiased FracMinHash containment;
- AF for MinHash and FracMinHash;
- k > 32;
- confidence intervals;
- quality filtering;
- parallel processing;
- compression other than gzip.

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
naive baselines, and streamed MinHash and FracMinHash sketches with brute-force
baselines over the distinct k-mer set (for FracMinHash across several scale
factors, including inputs large enough to exercise buffer compaction).
Estimators are checked against exact cases (identical, disjoint, and
reverse-complement inputs, and exact set statistics at `--scaled 1`),
hand-computed formula values, synthetic sets with known Jaccard, containment, and
cardinality (within 4 standard errors, fixed seeds), and the phiX174 genome in
`data/` mutated in memory at known substitution rates (1% and 5%). There,
containment ANI recovers 1 - p within 0.008 for MinHash and for FracMinHash at
`--scaled` 1 and 10.
SplitMix64 is checked against its published reference outputs, fmix64 by
inversion, and FNV-1a 64 against its published test vectors.

AIIPS tests check the following:
- The global probe set equals a brute-force baseline (each window encoded from
  scratch) and is bit-identical to the FracMinHash sketch with the same
  parameters.
- At `--probe-scaled 1`, Jaccard and containment equal exact k-mer set
  statistics. Cases: identical, unrelated, subsequence, fragmented,
  reverse-complement, and mutated inputs, and phiX174 mutated at 1% and 5%,
  where `C^(1/k2)` recovers the identity within 0.008.
- At the defaults, a subsequence of a 200 kb synthetic genome is fully
  contained, and the reverse containment and Jaccard are within 4 standard
  errors of the exact values.
- Occurrence counts include repeats, while global statistics use distinct
  hashes, and no k-mer crosses a record boundary.
- Every anchor and probe occurrence (contig, position, hash, strand, role)
  equals a brute-force baseline, across scales, seeds, ambiguous bases, and
  short records; anchors are re-derived from the hashes of their two k2-mers.
- Occurrence counts on a 400 kb genome are within 4 standard deviations of
  L/c, and roles split evenly and are constant per k-mer.
- Islands:
  - every probe is assigned exactly once, to the last anchor at or before it
    or to START;
  - next-anchor links stop at each contig end;
  - for every pair of anchors up to 12 apart, the slice between them equals
    both the walk through the islands in between and a naive position filter;
  - a reverse-complemented genome mirrors every position and flips every
    strand.
- The sketch file for the sample fixture matches, byte for byte, one produced
  by an independent Python implementation of the hashing and file format.
- Every strict-parsing rule rejects files that break it, and files round-trip
  exactly.
- Shared anchors and blocks:
  - identical and reverse-complement genomes match every unique anchor in
    identical oriented coordinates;
  - repeated anchors are ignored;
  - a lost or an extra intermediate anchor leaves the block spanning it, with
    its probes retrievable;
  - each kind of run break is counted;
  - cuts fall at the first anchor reaching the span;
  - the better-supported of two overlapping blocks is kept.
- Validation and estimates:
  - blocks with unrelated or disordered content are rejected;
  - scoring and coverage counts equal a naive recomputation;
  - identical and reverse-complement inputs give ANI 1 with identical counts;
  - thin evidence gives -1 with a status, never 0 or NaN.
- At the defaults, positional ANI on 600 kb genomes with 0.1-10% substitutions,
  and on phiX174 at 1% and 5%, is within 4 standard errors (plus 0.001) of the
  realized identity.
- AF for subsequences, deletions, insertions, inversions, translocations,
  fragmented assemblies, duplications and unrelated genomes matches the
  homologous fraction expected from the construction.
Integration tests exercise each command, all input combinations, gzip, known
counts, output creation and preservation, errors, help, and removal of the old
interface. They use Cargo's test executable and isolated temporary directories;
the downloaded genome and existing results are not modified. Fixtures are
synthetic; gzip fixtures have fixed timestamps. Symlink checks run on Unix.
