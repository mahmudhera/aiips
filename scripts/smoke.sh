#!/usr/bin/env bash
# Smoke-tests every aiips subcommand against the committed fixtures.
# Usage: scripts/smoke.sh PATH_TO_AIIPS_BINARY
#
# Fixtures are copied into a temporary work directory and every path passed to the
# binary is relative and /-separated, so printed path labels are identical on
# every platform. The sketch output has no parent directory: appending .aiips to
# a parented path would print the platform separator (escaped `\\` on Windows).
# Nothing is written into the repository tree.
set -euo pipefail

[ $# -eq 1 ] || { echo "usage: $0 PATH_TO_AIIPS_BINARY" >&2; exit 2; }
bin=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
[ -x "$bin" ] || [ -x "$bin.exe" ] || { echo "not executable: $1" >&2; exit 2; }
[ -x "$bin" ] || bin=$bin.exe # Windows
repo=$(cd "$(dirname "$0")/.." && pwd)

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
mkdir -p "$work/tests"
cp -R "$repo/tests/fixtures" "$work/tests/"
cd "$work"
errors=$(mktemp)
trap 'rm -rf "$work" "$errors"' EXIT

fail() {
    echo "FAIL: $*" >&2
    exit 1
}

# expect_stdout NAME EXPECTED ARGS...: exit 0, exact stdout, empty stderr.
expect_stdout() {
    local name=$1 expected=$2 actual status=0
    shift 2
    actual=$("$bin" "$@" 2>"$errors" | tr -d '\r') || status=$?
    [ "$status" -eq 0 ] || fail "$name: exit $status: $(cat "$errors")"
    [ ! -s "$errors" ] || fail "$name: unexpected stderr: $(cat "$errors")"
    if [ "$actual" != "$expected" ]; then
        diff <(printf '%s\n' "$expected") <(printf '%s\n' "$actual") >&2 || true
        fail "$name: stdout mismatch"
    fi
    echo "ok: $name"
}

# expect_help NAME ARGS...: exit 0 and stdout starts with the usage line.
expect_help() {
    local name=$1 actual status=0
    shift
    actual=$("$bin" "$@" 2>"$errors" | tr -d '\r') || status=$?
    [ "$status" -eq 0 ] || fail "$name: exit $status: $(cat "$errors")"
    [[ $actual == "Usage: aiips"* ]] || fail "$name: missing usage line: $actual"
    echo "ok: $name"
}

# expect_error NAME MESSAGE ARGS...: exit 1, empty stdout, stderr contains MESSAGE.
expect_error() {
    local name=$1 message=$2 actual status=0
    shift 2
    actual=$("$bin" "$@" 2>"$errors") || status=$?
    [ "$status" -eq 1 ] || fail "$name: expected exit 1, got $status"
    [ -z "$actual" ] || fail "$name: unexpected stdout: $actual"
    grep -qF -- "$message" "$errors" || fail "$name: stderr lacks '$message': $(cat "$errors")"
    echo "ok: $name"
}

stats=$'input\trecords\tbases\tnon_acgt_bases\tk\tvalid_kmers'
metrics=$'input_a\tinput_b\tjaccard\tcontainment_a_in_b\tcontainment_b_in_a\tani\taf_a_to_b\taf_b_to_a\tani_containment_a_in_b\tani_containment_b_in_a'
# Known baseline for tests/fixtures/sample.* at k=3 (tests/cli.rs, sketch_format!).
counts=$'2\t13\t2\t3\t7'
ones=$'1.000000\t1.000000\t1.000000\t1.000000\t-1\t-1\t1.000000\t1.000000'
fa=tests/fixtures/sample.fa
fq=tests/fixtures/sample.fq.gz
header='aiips-sketch v1 algorithm=minhash k=3 num_hashes=1000 seed=42 hash=murmur3-fmix64 canonical=true'

expect_help "top-level help" --help
expect_help "sketch help" sketch --help
expect_help "compare help" compare --help

expect_stdout "sketch fasta" \
    "$stats"$'\n'"$fa"$'\t'"$counts"$'\noutput\tsample.aiips' \
    sketch "$fa" -a minhash -k 3 -o sample
[ "$(head -n 1 sample.aiips)" = "$header" ] ||
    fail "sketch: unexpected header: $(head -n 1 sample.aiips)"
[ "$(wc -l <sample.aiips | tr -d ' ')" -eq 1002 ] ||
    fail "sketch: expected 1002 lines (2 header lines + 1000 hash values)"

expect_stdout "sketch gzip fastq" \
    "$stats"$'\n'"$fq"$'\t'"$counts"$'\noutput\treads.aiips' \
    sketch "$fq" -a minhash -k 3 -o reads

# Identical sequences in different formats must give identical sketches.
cmp -s sample.aiips reads.aiips || fail "fasta and gzip fastq sketches differ"

expect_stdout "compare identical sketches" \
    "$stats"$'\nsample.aiips\t'"$counts"$'\nreads.aiips\t'"$counts"$'\n\n'"$metrics"$'\nsample.aiips\treads.aiips\t'"$ones" \
    compare sample.aiips reads.aiips

expect_stdout "sketch with different k" \
    "$stats"$'\n'"$fa"$'\t2\t13\t2\t4\t5\noutput\tk4.aiips' \
    sketch "$fa" -a minhash -k 4 -o k4

expect_error "compare refuses mismatched parameters" "sketch parameters differ: k (3 vs 4)" \
    compare sample.aiips k4.aiips
expect_stdout "compare fasta vs gzip fastq (sketched in memory)" \
    "$stats"$'\n'"$fa"$'\t'"$counts"$'\n'"$fq"$'\t'"$counts"$'\n\n'"$metrics"$'\n'"$fa"$'\t'"$fq"$'\t'"$ones" \
    compare "$fa" "$fq" -k 3

expect_stdout "compare sketch vs fasta" \
    "$stats"$'\nsample.aiips\t'"$counts"$'\n'"$fa"$'\t'"$counts"$'\n\n'"$metrics"$'\nsample.aiips\t'"$fa"$'\t'"$ones" \
    compare sample.aiips "$fa"

expect_error "compare refuses -k differing from a sketch" "-k 4 does not match k=3" \
    compare sample.aiips "$fa" -k 4
expect_stdout "sketch defaults to minhash" \
    "$stats"$'\n'"$fa"$'\t'"$counts"$'\noutput\tdefault.aiips' \
    sketch "$fa" -k 3 -o default
cmp -s sample.aiips default.aiips || fail "sketch without -a differs from -a minhash"
expect_error "sketch refuses reserved algorithm" "bottomk is not implemented yet" \
    sketch "$fa" -a bottomk -k 3 -o missing
[ ! -e missing.aiips ] || fail "refused sketch created missing.aiips"

# FracMinHash: scaled=1 keeps every distinct k-mer of the tiny fixtures.
frac_header='aiips-sketch v1 algorithm=fracminhash k=3 scaled=1 seed=42 hash=murmur3-fmix64 canonical=true size='
expect_stdout "sketch fracminhash fasta" \
    "$stats"$'\n'"$fa"$'\t'"$counts"$'\noutput\tfrac.aiips' \
    sketch "$fa" -a fracminhash -k 3 --scaled 1 -o frac
[[ "$(head -n 1 frac.aiips)" == "$frac_header"* ]] ||
    fail "fracminhash sketch: unexpected header: $(head -n 1 frac.aiips)"
size=$(head -n 1 frac.aiips)
size=${size##*size=}
[ "$(wc -l <frac.aiips | tr -d ' ')" -eq $((size + 2)) ] ||
    fail "fracminhash sketch: expected $((size + 2)) lines (2 header lines + size values)"
expect_stdout "sketch fracminhash gzip fastq" \
    "$stats"$'\n'"$fq"$'\t'"$counts"$'\noutput\tfracreads.aiips' \
    sketch "$fq" -a fracminhash -k 3 --scaled 1 -o fracreads
cmp -s frac.aiips fracreads.aiips || fail "fasta and gzip fastq fracminhash sketches differ"
expect_stdout "compare identical fracminhash sketches" \
    "$stats"$'\nfrac.aiips\t'"$counts"$'\nfracreads.aiips\t'"$counts"$'\n\n'"$metrics"$'\nfrac.aiips\tfracreads.aiips\t'"$ones" \
    compare frac.aiips fracreads.aiips
expect_stdout "compare fasta vs gzip fastq with fracminhash" \
    "$stats"$'\n'"$fa"$'\t'"$counts"$'\n'"$fq"$'\t'"$counts"$'\n\n'"$metrics"$'\n'"$fa"$'\t'"$fq"$'\t'"$ones" \
    compare "$fa" "$fq" -a fracminhash -k 3 --scaled 1
expect_error "compare refuses mixed algorithms" "(sample.aiips: minhash, frac.aiips: fracminhash)" \
    compare sample.aiips frac.aiips
expect_error "sketch refuses --scaled with minhash" "--scaled applies only to fracminhash" \
    sketch "$fa" -a minhash -k 3 --scaled 10 -o missing
[ ! -e missing.aiips ] || fail "refused sketch created missing.aiips"

# AIIPS: anchor_k=5, anchor_scaled=2 keeps one anchor; probe_k=3,
# probe_scaled=1 keeps every 3-mer: 7 occurrences, 3 distinct canonical probes,
# 2 windows touching N, a 325-byte file.
aiips_flags=(--anchor-k 5 --anchor-scaled 2 --probe-k 3 --probe-scaled 1)
aiips_build=$'input\tanchor_occurrences\tdistinct_anchors\tunique_anchors\tprobe_occurrences\tdistinct_probes\tambiguous_kmers\tfile_bytes'
aiips_diag=$'input_a\tinput_b\tglobal_a\tglobal_b\tglobal_shared\tmax_containment\tputative_ani\tshared_anchors\trepeated_anchors\tbreaks_order\tbreaks_gap\tbreaks_ratio\tcandidate_blocks\toverlapping_blocks\trejected_few_hits\trejected_disorder\taccepted_blocks\tmapping_hits\tscoring_a\tscoring_b\tscoring_shared\texcluded_boundary_overlap\texcluded_repeats\tcovered_bases_a\tcovered_bases_b\taf_containment_a\taf_containment_b\tjaccard_status\tcontainment_status\tani_status\taf_status'
aiips_header='aiips-sketch v1 algorithm=aiips anchor_k=5 anchor_scaled=2 probe_k=3 probe_scaled=1 seed=42 hash=murmur3-fmix64 canonical=true body_bytes=133'
aiips_ones=$'1.000000\t1.000000\t1.000000\t-1\t-1\t-1\t1.000000\t1.000000'
# One shared anchor makes no block: positional ANI and AF are unsupported.
aiips_counts=$'3\t3\t3\t1.000000\t1.000000\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t-1\t-1\testimated\testimated\tno_supported_mapping\tno_supported_mapping'
expect_stdout "sketch aiips fasta" \
    "$stats"$'\n'"$fa"$'\t'"$counts"$'\n\n'"$aiips_build"$'\n'"$fa"$'\t1\t1\t1\t7\t3\t2\t325\noutput\taiips.aiips' \
    sketch "$fa" -a aiips "${aiips_flags[@]}" -o aiips
[ "$(head -n 1 aiips.aiips)" = "$aiips_header" ] ||
    fail "aiips sketch: unexpected header: $(head -n 1 aiips.aiips)"
[ "$(wc -c <aiips.aiips | tr -d ' ')" -eq 325 ] ||
    fail "aiips sketch: expected 325 bytes"
expect_stdout "sketch aiips gzip fastq" \
    "$stats"$'\n'"$fq"$'\t'"$counts"$'\n\n'"$aiips_build"$'\n'"$fq"$'\t1\t1\t1\t7\t3\t2\t325\noutput\taiipsreads.aiips' \
    sketch "$fq" -a aiips "${aiips_flags[@]}" -o aiipsreads
cmp -s aiips.aiips aiipsreads.aiips || fail "fasta and gzip fastq aiips sketches differ"
expect_stdout "compare identical aiips sketches" \
    "$stats"$'\naiips.aiips\t'"$counts"$'\naiipsreads.aiips\t'"$counts"$'\n\n'"$metrics"$'\naiips.aiips\taiipsreads.aiips\t'"$aiips_ones"$'\n\n'"$aiips_diag"$'\naiips.aiips\taiipsreads.aiips\t'"$aiips_counts" \
    compare aiips.aiips aiipsreads.aiips
expect_stdout "compare aiips sketch vs fasta" \
    "$stats"$'\naiips.aiips\t'"$counts"$'\n'"$fa"$'\t'"$counts"$'\n\n'"$metrics"$'\naiips.aiips\t'"$fa"$'\t'"$aiips_ones"$'\n\n'"$aiips_diag"$'\naiips.aiips\t'"$fa"$'\t'"$aiips_counts" \
    compare aiips.aiips "$fa"
expect_error "compare refuses aiips vs fracminhash" "(aiips.aiips: aiips, frac.aiips: fracminhash)" \
    compare aiips.aiips frac.aiips
expect_error "sketch refuses -k with aiips" "-k/--kmer-size does not apply to aiips" \
    sketch "$fa" -a aiips -k 3 -o missing
expect_error "sketch refuses even --probe-k" "odd integer from 1 to 31" \
    sketch "$fa" -a aiips --probe-k 4 -o missing
expect_error "sketch refuses even --anchor-k" "anchor k-mer length must be an odd integer" \
    sketch "$fa" -a aiips --probe-k 3 --anchor-k 4 -o missing
expect_error "compare refuses an option differing from the sketch" \
    "--anchor-scaled 125 does not match anchor_scaled=2 recorded in aiips.aiips" \
    compare aiips.aiips "$fa" --anchor-scaled 125
[ ! -e missing.aiips ] || fail "refused sketch created missing.aiips"

before=$(cat sample.aiips)
expect_error "sketch refuses existing output" "use -f or --force" \
    sketch "$fa" -a minhash -k 3 -o sample
[ "$(cat sample.aiips)" = "$before" ] || fail "sketch refusal: sample.aiips was modified"

expect_error "unknown command" "unknown command" bogus

echo "all smoke tests passed"
