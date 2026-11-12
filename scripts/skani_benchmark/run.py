#!/usr/bin/env python3
"""Run skani and aiips on simulated pairs and record estimates and costs.

Both tools start from FASTA (sketching included) on one thread:
    skani dist -t 1 -r A.fa -q B.fa
    aiips compare A.fa B.fa -a aiips [sketch options] [threshold options]
Wall time and peak resident memory come from /usr/bin/time. Results are
appended to a TSV, one row per (tool configuration, pair).
"""

import argparse
import json
import os
import platform
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

COLUMNS = [
    "label", "tool", "backbone", "seed", "scenario",
    "true_ani", "true_af_a", "true_af_b",
    "ani", "af_a", "af_b", "status", "seconds", "max_rss_mb",
]


def timed(command):
    """Runs `command` under /usr/bin/time; returns (stdout, seconds, rss_mb)."""
    if platform.system() == "Darwin":
        wrapper = ["/usr/bin/time", "-l"]
    else:
        wrapper = ["/usr/bin/time", "-v"]
    start = time.perf_counter()
    done = subprocess.run(wrapper + command, capture_output=True, text=True)
    seconds = time.perf_counter() - start
    if done.returncode != 0:
        raise RuntimeError(f"{' '.join(command)} failed: {done.stderr[-2000:]}")
    err = done.stderr
    # Wall time from a monotonic clock (finer than /usr/bin/time's 10 ms);
    # peak resident memory from /usr/bin/time.
    if platform.system() == "Darwin":
        rss = int(re.search(r"(\d+)\s+maximum resident set size", err).group(1)) / 2**20
    else:
        rss = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", err).group(1)) / 1024
    return done.stdout, seconds, rss


def run_skani(skani, a, b, options):
    out, seconds, rss = timed([skani, "dist", "-t", "1", "-r", a, "-q", b] + options)
    rows = [line.split("\t") for line in out.strip().splitlines()[1:]]
    if not rows:  # below skani's screening or AF thresholds
        return None, None, None, "not_reported", seconds, rss
    ani, af_ref, af_query = (float(x) / 100 for x in rows[0][2:5])
    return ani, af_ref, af_query, "reported", seconds, rss


def run_aiips(aiips, a, b, options):
    out, seconds, rss = timed([aiips, "compare", a, b, "-a", "aiips"] + options)
    lines = out.splitlines()
    row = lines[lines.index(next(l for l in lines if l.startswith("input_a\tinput_b\tjaccard"))) + 1]
    fields = row.split("\t")
    value = lambda x: None if x == "-1" else float(x)
    status = lines[-1].split("\t")[-2]  # ani_status
    return value(fields[5]), value(fields[6]), value(fields[7]), status, seconds, rss


def pairs(root, seeds, backbones, only):
    for group in sorted(os.listdir(root)):
        backbone, _, seed = group.rpartition("_s")
        if backbones and backbone not in backbones:
            continue
        if seeds and int(seed) not in seeds:
            continue
        for scenario in sorted(os.listdir(os.path.join(root, group))):
            path = os.path.join(root, group, scenario)
            if not os.path.isfile(os.path.join(path, "truth.json")):
                continue
            if only and only not in scenario:
                continue
            yield os.path.join(root, group, "a.fa"), path


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--pairs", default="data/test_against_skani/pairs")
    parser.add_argument("--out", required=True, help="TSV to append to")
    parser.add_argument("--label", required=True, help="name of this tool configuration")
    parser.add_argument("--tool", choices=["skani", "aiips"], required=True)
    parser.add_argument("--skani", default="skani")
    parser.add_argument("--aiips", default="target/release/aiips")
    parser.add_argument("--options", default="", help="extra tool options, space-separated")
    parser.add_argument("--seeds", nargs="*", type=int, default=[])
    parser.add_argument("--backbones", nargs="*", default=[])
    parser.add_argument("--only", help="run only scenarios whose name contains this")
    parser.add_argument("--workers", type=int, default=1,
                        help="pairs run concurrently; use 1 when measuring time")
    args = parser.parse_args()

    jobs = list(pairs(args.pairs, args.seeds, args.backbones, args.only))
    options = args.options.split()

    def one(job):
        a, path = job
        b = os.path.join(path, "b.fa")
        with open(os.path.join(path, "truth.json")) as handle:
            truth = json.load(handle)
        if args.tool == "skani":
            result = run_skani(args.skani, a, b, options)
        else:
            result = run_aiips(args.aiips, a, b, options)
        return [args.label, args.tool, truth["backbone"], truth["seed"], truth["scenario"],
                truth["ani"], truth["af_a"], truth["af_b"], *result]

    new = not os.path.exists(args.out)
    with open(args.out, "a") as out, ThreadPoolExecutor(args.workers) as pool:
        if new:
            out.write("\t".join(COLUMNS) + "\n")
        for done, row in enumerate(pool.map(one, jobs), 1):
            out.write("\t".join("NA" if x is None else str(x) for x in row) + "\n")
            out.flush()
            print(f"\r{args.label}: {done}/{len(jobs)}", end="", file=sys.stderr, flush=True)
    print(file=sys.stderr)


if __name__ == "__main__":
    main()
