[![CI](https://github.com/mahmudhera/aiips/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/mahmudhera/aiips/actions/workflows/ci.yml)
[![Coverage](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/mahmudhera/aiips/badges/coverage.json)](https://github.com/mahmudhera/aiips/actions/workflows/ci.yml?query=branch%3Amain)

## Summary

AIIPS means Anchor Indexed Island Probe Sketching. This project explores nested
sketching. The CLI builds classic MinHash or FracMinHash sketches of FASTA/FASTQ
files and compares two sketches, estimating Jaccard, containment, and ANI. It
also builds AIIPS sketches (`-a aiips`): two nested, positional FracMinHash
samples, sparse long anchors and dense short probes between them, compared by
walking from matched anchors through the islands between anchors, which
estimates both ANI and alignment fraction (AF). Other sketch algorithms (bottomk,
minimizer, open-syncmer, closed-syncmer, strobemer) are reserved names, not yet
implemented.

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

./target/release/aiips sketch genome_a.fa -a aiips -o genome_a
./target/release/aiips compare genome_a.aiips genome_b.fa
./target/release/aiips compare genome_a.fa genome_b.fa -a aiips --symmetric-probes --bidirectional
```

The previous `aiips INPUT...` interface has been removed. A subcommand is required.
Stdin (`-`) is not supported. Use `--` before positional filenames beginning with
`-`; options must precede that separator. For output names beginning with `-`,
pass an explicit relative path, such as `-o ./-sample`.

### Create a sketch

```text
aiips sketch INPUT -o OUTPUT [-a ALGORITHM] [-k LENGTH] [-s NUM] [--scaled S] [--seed SEED] [-f]
```

`-a`/`--algorithm` is `minhash` (the default, so commands without it keep
working), `fracminhash`, or `aiips`.
Reserved names (bottomk, minimizer, open-syncmer, closed-syncmer, strobemer) are
refused as "not implemented yet"; other names are unknown. MinHash and
FracMinHash take `-k` from 1 to 32 (default 21); all three take `--seed`, any
unsigned 64-bit integer (default 42). MinHash also takes `-s`/`--num-hashes` from
2 to 1,000,000 (default 1000); FracMinHash takes `--scaled`, an integer scale
factor from 1 to 2^64 - 1 (default 1000). AIIPS takes `--anchor-k` (2 to 32,
default 21), `--anchor-scaled` (default 500), `--probe-k` (1 to anchor k - 1,
default 15) and `--probe-scaled` (default 100) instead of `-k`, `-s` and
`--scaled`. Options that do not apply to the chosen algorithm (`--scaled` with
MinHash, `-s` with FracMinHash, `-k` with AIIPS, the AIIPS options with the
others) are errors. Each option may be given once. An
input without any valid k-mer is refused and creates no output. A FracMinHash
sketch that keeps no hash value is still written, with a warning on stderr,
because it is a valid sample; it cannot be compared.

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

An AIIPS sketch lists each record (`c LENGTH`), then each of its anchors in
position order (`a POS STRAND HASH` for a k1 anchor, `p POS STRAND HASH` for a
promoted k2 anchor, always the record's first), each followed by the probes of
the island after it (`q POS STRAND HASH`). Positions are zero-based starts in the
record; the strand is `+` when the k-mer's forward code is canonical, `-` when its
reverse complement is. `contigs`, `anchors` and `probes` are the line counts:

```text
aiips-sketch v1 algorithm=aiips anchor_k=21 anchor_scaled=500 probe_k=15 probe_scaled=100 seed=42 hash=murmur3-fmix64 canonical=true contigs=1 anchors=2084 probes=9353
records=1 bases=1000000 non_acgt_bases=0 valid_kmers=999986
c 1000000
p 30 + 56731525506569147
q 30 + 56731525506569147
q 44 + 120953976884060884
q 100 - 89964033739757942
...
```

Reading is strict: unknown, missing, or duplicate fields, a wrong value count,
non-decimal values, FracMinHash values out of order, repeated, or above the
`scaled` threshold, carriage returns, a missing final newline, or another version
are errors. AIIPS records must hold their k-mers in order, inside the record,
below their thresholds, with no probe overlapping a k1 anchor and the promoted
anchor (if any) first and repeated as its island's first probe. A FracMinHash
file may not exceed 21 bytes per declared value plus 4096 bytes, an AIIPS file 13
bytes per declared record plus 36 per declared k-mer plus 4096; every other sketch
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

An AIIPS sketch samples two k-mer sizes with two independent hash functions, both
by the FracMinHash rule. Per record (contig):

- Anchors: canonical k1-mers whose hash `fmix64(x XOR seed_0)` is at most
  `floor((2^64 - 1) / S1)` (`--anchor-k`, default 21; `--anchor-scaled`, default
  500), i.e. a FracMinHash sample at k1.
- Probes: canonical k2-mers, k2 < k1, whose hash `fmix64(x XOR seed_1)` is at most
  `floor((2^64 - 1) / S2)` (`--probe-k`, default 15; `--probe-scaled`, default
  100), excluding every k2-mer that overlaps an anchor. `seed_0` and `seed_1` are
  the first two SplitMix64 outputs of `--seed`.
- Islands: the stretch from an anchor to the next anchor or the record end. Every
  probe belongs to the island it lies in; island lengths follow from the stored
  positions.
- The record's first probe before its first anchor (or its first probe, if it has
  no anchor) becomes an anchor too, so the first sampled k-mer of every record is
  an anchor and no island lacks a preceding anchor. It also stays a probe of its
  island (a design decision: it is a sampled k2-mer like the other genome's copy,
  so containment counts it like any probe).

Every k-mer keeps its position and strand. With the defaults a 5 Mbp genome has
about 10,000 anchors and 50,000 probes. Sketching is one pass and two hashes per
base:

- Each record is decoded in 64 Ki-base chunks straight from the parser's buffer:
  line breaks are dropped on the fly (lines are cut where the first line's width
  predicts, checked against the record's line-break count, else found by
  search) and eight bases are encoded at a time with word arithmetic.
- A chunk is scanned as 8 segments side by side, so the rolling codes, both hashes
  and both threshold tests of 8 positions advance together in vector registers.
  Lanes read eight bases at a time, transposed with word arithmetic; each lane
  warms up on the k1 - 1 bases before its segment (carried over between chunks).
  The scan is compiled for AVX-512, AVX2 and baseline x86-64 and chosen at
  runtime (other architectures use the baseline code).
- About nine steps in ten sample nothing and cost only the threshold test; the
  others are buffered with their hashes and strands and emitted in bulk.
- Per-record buffers are kept by the thread and reused by the next sketch.

Every path gives identical sketches; the tests compare them with a brute-force
builder, line-wrapped input included. Besides the sketch, memory holds the parser's
buffer (the longest record) and the chunk buffers. `valid_kmers` counts
probe-length (k2) windows.

`compare` with AIIPS sketches, or with `-a aiips` for two sequence inputs:

1. Seeds: anchor pairs with equal hash (k1 anchors with k1 anchors; a promoted
   anchor with any sampled k2-mer of the other genome). A hash occurring more than
   32 times in either genome is skipped. Strands give the relative orientation.
2. Walks: from a seed, the island after the anchor in one genome is compared with
   the island after (same orientation) or before (opposite orientation) the
   matched k-mer in the other, each genome keeping its own islands. When the next
   anchors are the same k-mer both genomes cross them and the walk re-anchors;
   otherwise the genome whose island ends first (by island length) crosses its
   next anchor, so an island can pair with several of the other genome's. Island
   ends are compared after correcting for the offset between the genomes measured
   at the last shared probe, which indels change between re-anchorings. A seed
   walks forward from the anchor in each genome whose anchor it is (one walk if the
   strands agree, two if not), unless an earlier walk in that orientation crossed
   both its anchors.
3. An island is judged when the walk leaves it, from all the probes it shared with
   the islands it was paired with: it is high enough when that fraction is at
   least `X^k2` (`--min-island-ani X`, default 0.8, so effectively one shared
   probe). When both genomes leave an island together, either passing suffices
   ("containment in either genome"); matched next anchors always pass. The walk
   goes beyond low islands but stops after more than `--max-gap G` consecutive low
   judgements (default 50) or at a record end, where the other genome's island is
   clipped. It is then cut back to its last high judgement and trimmed to its last
   shared k-mer, except that a walk reaching a record end with no low island since
   its last high one keeps its islands up to that end.
4. Overlapping walks (repeats) are resolved best first by score (shared probes plus
   matched anchors). An island's probe counts come from the first walk claiming it
   and its aligned bases from the union of the walks' spans; each anchor counts
   once.
5. ANI: the maximum-likelihood per-base identity `q` from the anchors crossed
   inside walks (`n1` of `N1` matched, k1) and the probes of aligned islands (`n2`
   of `N2` shared, k2), pooled over both genomes, assuming each k-mer of length k
   matches independently with probability `q^k`:
   `log L(q) = n1 k1 ln q + (N1 - n1) ln(1 - q^k1) + n2 k2 ln q + (N2 - n2) ln(1 - q^k2)`.
   With one k-mer size this is the containment estimate `(n/N)^(1/k)` (Hera,
   Pierce-Ward & Koslicki 2023); `log L` is concave in `ln q`, and its maximum is
   found by bisection. Anchors and probes never overlap, so they share no bases.
   Seed anchors are left out because they were selected for matching. `ani_a` and
   `ani_b` are the same estimate from one genome's counts.
6. AF: aligned bases over all bases, in each genome.

Output after the statistics table (`-1` for an ANI when nothing aligned):

```text
input_a\tinput_b\tani\taf_a_to_b\taf_b_to_a\tani_a\tani_b\tanchors_matched\tanchors\tprobes_matched\tprobes
a.fa\tb.fa\t0.950404\t0.999604\t0.999604\t0.950592\t0.950217\t6730\t19766\t41590\t89049
```

(Here B is A with 5% of its 5 Mbp substituted, true ANI 0.94998, compared with
`--symmetric-probes --bidirectional`. The counts are summed over both genomes.)

Two opt-in flags change the comparison, not the sketch:

- `--symmetric-probes`: a probe counts only if it also avoids the other genome's
  anchors at the same walk offset (re-estimated at each shared probe to follow
  indels). Without it the exclusion rule differs between genomes wherever an
  anchor is mutated away in one of them: that genome keeps probes the other
  lacks, which lowers containment (by about 6% at 90% ANI with the defaults). On
  the benchmark below the plain estimate is 0.1-0.35 ANI points low on
  substitution-only pairs; with the flag it is within 0.05 on average.
- `--bidirectional`: seeds also walk backward, which recovers the part of each
  homologous stretch before its first matched anchor (contig starts, breakpoints).
  Second-direction walks start only from seeds no earlier walk passed.

Decisions beyond the specified algorithm, flagged here: islands are judged when
completed rather than per island pair (pairs of islands overlapping by a few
bases otherwise fail and end walks); low islands inside a walk count toward ANI
(dropping them biases ANI up); walks tolerate up to `G` low islands (with about 5
probes per island, an island has no shared probe with probability 0.35 at 90% ANI
even when homologous; `G` was tuned on the simulated benchmark, where larger
values only extend true alignments); walks are trimmed to their last shared
k-mer; seed anchors are excluded from the likelihood.

Limitations: indels destroy k-mers like substitutions, so ANI computed over
gap-free aligned columns is underestimated when indels are frequent (as for
skani); a straddling island can add a few hundred bases of AF at each breakpoint;
at low ANI with short contigs some contigs hold no matched seed and are missed;
records are limited to 2^32 - 1 bases and k1 to 32.

### Compare two inputs

```text
aiips compare INPUT_A INPUT_B [-a ALGORITHM] [-k LENGTH] [-s NUM] [--scaled S] [--seed SEED]
```

Each input is a sketch or a plain/gzip FASTA/FASTQ file; mixed pairs are
supported, so commands written before sketches existed keep working. Files named
`.aiips` (case-insensitive) must be sketches; other files are sketches when they
start with the sketch header and are otherwise read as sequence.

Sequence inputs are sketched in memory. When an input is a sketch, its algorithm
and parameters are used, and every given option must equal the recorded value
(`-a` its algorithm, `-k` its `k`, `-s` its `num_hashes`, `--scaled` its
`scaled`, `--seed` its `seed`); options that do not apply to its algorithm are
errors. These checks also apply when both inputs are sketches. When both inputs
are sequences, the algorithm comes from `-a` (default `minhash`) and the
parameters from the options, with the `sketch` defaults. A sequence input without
valid k-mers is an error.

Unless both inputs are of one implemented algorithm, `compare` halts and names
each input's type. Two MinHash sketches must agree in `k`, `num_hashes`, `seed`,
`hash`, and `canonical`; two FracMinHash sketches in `k`, `scaled`, `seed`, `hash`,
and `canonical`; two AIIPS sketches in `anchor_k`, `anchor_scaled`, `probe_k`,
`probe_scaled`, `seed`, `hash`, and `canonical`. Every differing field is
reported. The AIIPS comparison options (`--max-gap`, `--min-island-ani`,
`--symmetric-probes`, `--bidirectional`) are errors with the other algorithms;
AIIPS output is described above. An empty FracMinHash sketch
halts the comparison, because its containment is undefined.

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

## Implementation

- `src/fastx.rs` streams one record at a time using
  [needletail](https://docs.rs/needletail/), also used by
  [skani](https://github.com/bluenote-1577/skani/blob/main/src/file_io.rs).
  Only gzip compression is enabled to keep dependencies small.
- `src/kmer.rs` iterates borrowed forward-strand windows and zero-based positions
  without per-k-mer allocation, and rolling canonical 2-bit k-mer codes for
  k <= 32. Both take O(sequence length) time and O(1) extra space.
- `src/hash.rs` provides MurmurHash3's fmix64 and SplitMix64.
- `src/sketch/` holds the algorithm registry and `.aiips` format (`mod.rs`),
  MinHash construction and estimators (`minhash.rs`), and FracMinHash
  construction and estimators (`fracminhash.rs`). FracMinHash comparison merges
  the two sorted sketches in linear time.
- `src/sketch/aiips.rs` holds AIIPS sketching (streamed decoding and a
  vectorized scan with runtime instruction-set dispatch, no extra dependencies),
  the format, the walk, and the likelihood. Comparison joins the sorted k1 anchor
  hashes, matches the few promoted anchors by binary search instead of sorting all
  probes, and touches each island a bounded number of times; its memory is linear
  in the sketch sizes.
- The binary separates argument parsing, command execution and output handling,
  and TSV formatting. Memory scales with parser buffering and the largest record.

Sketch algorithms other than MinHash, FracMinHash and AIIPS, debiased FracMinHash
containment, AF for MinHash and FracMinHash, k > 32, confidence intervals, quality
filtering, parallel processing, and compression other than gzip remain
unimplemented.

## Benchmark against skani

`scripts/simulate.py OUTDIR` (standard library, fixed seeds) writes 32 pairs of
genomes (5 Mbp by default, `--length` to change) with known ANI and AF, and
`scripts/benchmark.py` runs aiips and skani on them: one thread, FASTA inputs (so
sketching is included), optionally pinned to one core with `--cpu`. Wall time is
measured around each process, user and system CPU time come from `getrusage`, and
peak memory from `/usr/bin/time`; the fastest of `--repeat` runs is kept.
`--other-aiips LABEL=PATH` adds another aiips binary, e.g. an earlier build:

```sh
python3 scripts/simulate.py sim --length 10000000
python3 scripts/benchmark.py sim/truth.tsv results.tsv --aiips target/release/aiips \
    --skani "$(command -v skani)" --repeat 3 --cpu 1 \
    --variant both="--symmetric-probes --bidirectional"
```

Truth: ANI is one minus substitutions over aligned (homologous, gap-free) columns;
AF is the fraction of a genome's bases with a homolog in the other. Scenarios:
uniform substitutions (ANI 99.9% to 80%); B holding 20, 50 or 90% of A in blocks
amid foreign sequence, at 1, 5 and 10% divergence; substitutions plus indels (a
tenth of the substitution rate, lengths mean 3, at most 30) with B in one record
or in contigs of mean 50, 10 or 3 kbp in random orientation; and A with repeat
families (40 families, 2-5 copies of 1-7 kbp, up to 1% diverged) against B
substituted, cut into about 40 blocks, shuffled and partly inverted.

Results for 10 Mbp genomes (`benchmarks/simulated-10mb.tsv`; truth in
`benchmarks/simulated-10mb-truth.tsv`), with skani 0.4.0 at its defaults and
`--min-af 0`, and the initial aiips implementation for reference. Errors are mean
absolute errors in percentage points, AF over both directions. Time and memory
are per pair, on one pinned core of a shared Xeon Gold 5416S; timings vary with
the machine's load, so all tools were run interleaved, fastest of 3:

| Scenarios | Tool | ANI MAE | AF MAE | Wall (s) | User (s) | System (s) | Memory (MiB) |
|---|---|---|---|---|---|---|---|
| Substitutions (8) | aiips | 0.14 | 1.3 | 0.103 | 0.074 | 0.028 | 26.2 |
| Substitutions (8) | aiips, both flags | 0.07 | 0.4 | 0.110 | 0.083 | 0.026 | 26.2 |
| Substitutions (8) | skani (1 pair not reported) | 0.44 | 11.2 | 0.124 | 0.083 | 0.041 | 39.9 |
| Partial overlap (9) | aiips | 0.19 | 2.2 | 0.082 | 0.060 | 0.021 | 24.3 |
| Partial overlap (9) | aiips, both flags | 0.06 | 0.5 | 0.086 | 0.064 | 0.022 | 24.3 |
| Partial overlap (9) | skani | 1.09 | 3.1 | 0.089 | 0.059 | 0.030 | 34.1 |
| Indels + draft contigs (12) | aiips | 0.63 | 10.7 | 0.084 | 0.065 | 0.018 | 24.2 |
| Indels + draft contigs (12) | aiips, both flags | 0.43 | 7.8 | 0.087 | 0.068 | 0.019 | 24.2 |
| Indels + draft contigs (12) | skani | 0.84 | 11.7 | 0.085 | 0.061 | 0.024 | 33.6 |
| Repeats + rearrangements (3) | aiips | 0.14 | 0.6 | 0.069 | 0.050 | 0.018 | 26.0 |
| Repeats + rearrangements (3) | aiips, both flags | 0.06 | 0.3 | 0.075 | 0.057 | 0.017 | 26.0 |
| Repeats + rearrangements (3) | skani | 0.29 | 4.9 | 0.076 | 0.048 | 0.027 | 39.9 |
| All 32 | aiips | 0.34 | 5.0 | 0.087 | 0.064 | 0.022 | 24.9 |
| All 32 | aiips, both flags | 0.20 | 3.2 | 0.091 | 0.070 | 0.021 | 24.9 |
| All 32 | skani (1 pair not reported) | 0.77 | 8.4 | 0.095 | 0.065 | 0.030 | 35.9 |
| All 32 | aiips, initial implementation | 0.34 | 5.0 | 0.118 | 0.086 | 0.031 | 38.2 |

Substitution sweep, ANI / AF of A:

| True ANI | aiips | aiips, both flags | skani |
|---|---|---|---|
| 0.9990 | 0.9990 / 1.000 | 0.9990 / 1.000 | 0.9991 / 1.000 |
| 0.9900 | 0.9896 / 1.000 | 0.9903 / 1.000 | 0.9892 / 0.999 |
| 0.9699 | 0.9686 / 1.000 | 0.9704 / 1.000 | 0.9677 / 0.993 |
| 0.9499 | 0.9481 / 1.000 | 0.9505 / 1.000 | 0.9438 / 0.985 |
| 0.9199 | 0.9177 / 1.000 | 0.9205 / 1.000 | 0.9153 / 0.948 |
| 0.8999 | 0.8975 / 1.000 | 0.9006 / 1.000 | 0.8985 / 0.875 |
| 0.8497 | 0.8472 / 0.999 | 0.8505 / 1.000 | 0.8342 / 0.415 |
| 0.8000 | 0.7991 / 0.896 | 0.8025 / 0.968 | not reported |

Both tools underestimate ANI when indels are frequent (k-mers spanning an indel
are lost), and both lose AF on 3 kbp contigs at 90% ANI (aiips 0.41-0.50, skani
0.50). Simulated genomes are uniform random sequence, without the compositional
structure of real genomes; skani's learned ANI regression was trained on real
assemblies. `benchmarks/simulated-5mb.tsv` holds the same benchmark on 5 Mbp
genomes from the initial implementation.

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
inversion.
AIIPS sketches are checked against a brute-force builder (every window hashed
from scratch, exclusion and promotion applied literally) on records with
ambiguous bases, lane and chunk boundaries, short and empty records, and line
breaks (LF, CRLF, ragged and blank lines), for every instruction set the machine
offers; the word encoder against the lookup table for every byte, the 8x8 byte
transpose against a naive one, and decoding against plain line-break removal
(including a mispredicted line width, which must be caught); the likelihood against its closed form and a
grid search; and comparisons against exact cases (identical genomes, reverse
complement cut into contigs, half shared with a breakpoint, unrelated genomes) and
1 Mbp genomes substituted at 1% and 5%.
Integration tests exercise each command, all input combinations, gzip, known
counts, output creation and preservation, errors, help, and removal of the old
interface. They use Cargo's test executable and isolated temporary directories;
the downloaded genome and existing results are not modified. Fixtures are
synthetic; gzip fixtures have fixed timestamps. Symlink checks run on Unix.
