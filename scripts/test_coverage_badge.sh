#!/usr/bin/env bash
# Tests scripts/coverage_badge.sh with synthetic cargo llvm-cov summaries.
# Usage: scripts/test_coverage_badge.sh
set -euo pipefail

badge=$(cd "$(dirname "$0")" && pwd)/coverage_badge.sh
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
failures=0

summary() { # summary JSON_VALUE_FOR_PERCENT
    printf '{"data":[{"totals":{"lines":{"count":10,"covered":5,"percent":%s}}}]}' "$1" \
        >"$work/in.json"
}

# expect PERCENT MESSAGE COLOR: badge succeeds with exactly this content.
expect() {
    summary "$1"
    rm -f "$work/out.json"
    local expected actual
    expected="{\"schemaVersion\":1,\"label\":\"coverage\",\"message\":\"$2\",\"color\":\"$3\"}"
    if ! "$badge" "$work/in.json" "$work/out.json" 2>"$work/err"; then
        echo "FAIL: $1: exited non-zero: $(cat "$work/err")" >&2
        failures=$((failures + 1))
        return
    fi
    actual=$(cat "$work/out.json")
    if [ "$actual" != "$expected" ]; then
        echo "FAIL: $1: expected $expected, got $actual" >&2
        failures=$((failures + 1))
        return
    fi
    echo "ok: $1 -> $2 $3"
}

# reject NAME FILE_CONTENT: badge exits non-zero and writes no output file.
reject() {
    printf '%s' "$2" >"$work/in.json"
    rm -f "$work/out.json"
    if "$badge" "$work/in.json" "$work/out.json" 2>/dev/null; then
        echo "FAIL: $1: accepted invalid input" >&2
        failures=$((failures + 1))
    elif [ -e "$work/out.json" ]; then
        echo "FAIL: $1: wrote output despite failing" >&2
        failures=$((failures + 1))
    else
        echo "ok: rejects $1"
    fi
}

expect 0 "0.0%" red
expect 59.9 "59.9%" red
expect 60 "60.0%" yellow
expect 74.9 "74.9%" yellow
expect 75 "75.0%" green
expect 89.9 "89.9%" green
expect 90 "90.0%" brightgreen
expect 100 "100.0%" brightgreen
expect 87.34 "87.3%" green
expect 87.36 "87.4%" green
# The color follows the rounded value shown on the badge.
expect 89.96 "90.0%" brightgreen
expect 1e-05 "0.0%" red

reject "missing field" '{"data":[{"totals":{"lines":{}}}]}'
reject "empty data" '{"data":[]}'
reject "string percent" '{"data":[{"totals":{"lines":{"percent":"87"}}}]}'
reject "null percent" '{"data":[{"totals":{"lines":{"percent":null}}}]}'
reject "negative percent" '{"data":[{"totals":{"lines":{"percent":-1}}}]}'
reject "percent above 100" '{"data":[{"totals":{"lines":{"percent":100.1}}}]}'
reject "invalid JSON" '{"data":'

if "$badge" "$work/in.json" 2>/dev/null; then
    echo "FAIL: accepted a missing argument" >&2
    failures=$((failures + 1))
else
    echo "ok: rejects missing argument"
fi

[ "$failures" -eq 0 ] || { echo "$failures failure(s)" >&2; exit 1; }
echo "all coverage badge tests passed"
