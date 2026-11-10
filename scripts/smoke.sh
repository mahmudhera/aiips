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

before=$(cat sample.aiips)
expect_error "sketch refuses existing output" "use -f or --force" \
    sketch "$fa" -a minhash -k 3 -o sample
[ "$(cat sample.aiips)" = "$before" ] || fail "sketch refusal: sample.aiips was modified"

expect_error "unknown command" "unknown command" bogus

echo "all smoke tests passed"
