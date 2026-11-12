#!/usr/bin/env python3
"""Interleaved cost comparison: wall time and peak memory per pair.

For each pair, every configuration is run back to back, `--repeats` times in
rotation (after one warm-up each), so all configurations see the same machine
state. Reports the median wall time and the largest peak resident memory per
configuration and pair, then the per-pair ratios to a reference configuration.
"""

import argparse
import os
import platform
import re
import statistics
import subprocess
import time
from collections import defaultdict

from run import pairs


def peak_rss_mb(command):
    flag = "-l" if platform.system() == "Darwin" else "-v"
    err = subprocess.run(["/usr/bin/time", flag] + command, capture_output=True, text=True).stderr
    if platform.system() == "Darwin":
        return int(re.search(r"(\d+)\s+maximum resident set size", err).group(1)) / 2**20
    return int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", err).group(1)) / 1024


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--pairs", default="data/test_against_skani/pairs")
    parser.add_argument("--skani", default="skani")
    parser.add_argument("--aiips", default="target/release/aiips")
    parser.add_argument("--seeds", nargs="*", type=int, default=[])
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument(
        "--config", action="append", default=[], metavar="LABEL=TOOL[:OPTIONS]",
        help="e.g. skani=skani or aiips_p80=aiips:--probe-scaled 80 (repeatable)",
    )
    parser.add_argument("--reference", default="skani")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    configs = {}
    for spec in args.config:
        label, _, rest = spec.partition("=")
        tool, _, options = rest.partition(":")
        configs[label] = (tool, options.split())

    def command(tool, options, a, b):
        if tool == "skani":
            return [args.skani, "dist", "-t", "1", "-r", a, "-q", b] + options
        return [args.aiips, "compare", a, b, "-a", "aiips"] + options

    rows = []
    jobs = list(pairs(args.pairs, args.seeds, [], None))
    for done, (a, path) in enumerate(jobs, 1):
        b = os.path.join(path, "b.fa")
        times = defaultdict(list)
        commands = {label: command(*configs[label], a, b) for label in configs}
        for cmd in commands.values():
            subprocess.run(cmd, capture_output=True, check=True)  # warm-up
        for _ in range(args.repeats):
            for label, cmd in commands.items():
                start = time.perf_counter()
                subprocess.run(cmd, capture_output=True, check=True)
                times[label].append(time.perf_counter() - start)
        scenario = os.path.relpath(path, args.pairs)
        for label, cmd in commands.items():
            rows.append((scenario, label, statistics.median(times[label]), peak_rss_mb(cmd)))
        print(f"\rcost: {done}/{len(jobs)}", end="", flush=True)
    print()

    with open(args.out, "w") as out:
        out.write("pair\tlabel\tmedian_seconds\tmax_rss_mb\n")
        for row in rows:
            out.write("\t".join(map(str, row)) + "\n")

    by_pair = defaultdict(dict)
    for pair, label, seconds, rss in rows:
        by_pair[pair][label] = (seconds, rss)
    print(f"| configuration | median time / {args.reference} | pairs at or below | "
          f"median memory / {args.reference} | pairs at or below |")
    print("|---|---|---|---|---|")
    for label in configs:
        if label == args.reference:
            continue
        t = [v[label][0] / v[args.reference][0] for v in by_pair.values()]
        m = [v[label][1] / v[args.reference][1] for v in by_pair.values()]
        print(f"| {label} | {statistics.median(t):.3f} | {sum(x <= 1 for x in t)}/{len(t)} "
              f"| {statistics.median(m):.3f} | {sum(x <= 1 for x in m)}/{len(m)} |")
    ref = [v[args.reference][0] for v in by_pair.values()]
    print(f"\n{args.reference}: median {statistics.median(ref) * 1000:.1f} ms per pair")


if __name__ == "__main__":
    main()
