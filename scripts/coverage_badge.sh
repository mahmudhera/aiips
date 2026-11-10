#!/usr/bin/env bash
# Converts a `cargo llvm-cov --json --summary-only` report into a shields.io
# endpoint badge file showing line coverage.
# Usage: scripts/coverage_badge.sh SUMMARY_JSON OUTPUT_JSON
#
# Colors are display-only (coverage is reported, never enforced) and use the
# rounded value shown on the badge, so the color always matches the label.
set -euo pipefail

[ $# -eq 2 ] || { echo "usage: $0 SUMMARY_JSON OUTPUT_JSON" >&2; exit 2; }

# Refuse anything but a percentage, so a broken report never publishes a badge.
if ! percent=$(jq -er '.data[0].totals.lines.percent
        | select(type == "number" and . >= 0 and . <= 100)' "$1" 2>/dev/null); then
    echo "$1: no line coverage percentage at .data[0].totals.lines.percent" >&2
    exit 1
fi

message=$(awk -v p="$percent" 'BEGIN { printf "%.1f", p }')
color=$(awk -v m="$message" 'BEGIN {
    if (m >= 90) print "brightgreen"
    else if (m >= 75) print "green"
    else if (m >= 60) print "yellow"
    else print "red"
}')

jq -cn --arg message "$message%" --arg color "$color" \
    '{schemaVersion: 1, label: "coverage", message: $message, color: $color}' >"$2"
