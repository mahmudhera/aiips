#!/usr/bin/env python3
"""Benchmarks aiips against skani on simulated pairs with known ANI and AF.

Usage:
  scripts/benchmark.py TRUTH_TSV OUT_TSV --aiips BIN --skani BIN
      [--variant LABEL=ARGS ...] [--skani-args ARGS] [--repeat N]

TRUTH_TSV comes from scripts/simulate.py. Each tool compares the two FASTA
files directly (sketching included) on one thread; wall time and peak resident
memory come from /usr/bin/time. A missing estimate (e.g. skani screening a
pair out) is recorded as NA. Standard library only.

--variant adds an aiips run with extra `compare` arguments, e.g.
  --variant gap0="--max-gap 0"; the default variant has no extra arguments.
"""

import argparse
import csv
import os
import shlex
import statistics
import subprocess
import tempfile


def timed(command):
    """Runs command; returns (stdout, wall seconds, max RSS in KiB)."""
    with tempfile.NamedTemporaryFile("r", suffix=".time") as t:
        result = subprocess.run(
            ["/usr/bin/time", "-f", "%e %M", "-o", t.name] + command,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"{shlex.join(command)} failed: {result.stderr}")
        wall, rss = t.read().split()[-2:]
    return result.stdout, float(wall), int(rss)


def run_skani(skani, extra, a, b):
    out, wall, rss = timed([skani, "dist", "-t", "1", "--min-af", "0", *extra, "-q", b, "-r", a])
    rows = [line.split("\t") for line in out.strip().split("\n")[1:] if line]
    if not rows:
        return None, None, None, wall, rss
    # Columns: Ref_file Query_file ANI Align_fraction_ref Align_fraction_query ...
    ani, af_ref, af_query = (float(x) / 100 for x in rows[0][2:5])
    return ani, af_ref, af_query, wall, rss


def run_aiips(aiips, extra, a, b):
    out, wall, rss = timed([aiips, "compare", a, b, "-a", "aiips", *extra])
    fields = out.strip().split("\n")[-1].split("\t")
    ani = float(fields[2])
    return (None if ani < 0 else ani), float(fields[3]), float(fields[4]), wall, rss


def fmt(x, digits=6):
    return "NA" if x is None else f"{x:.{digits}f}"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("truth")
    parser.add_argument("out")
    parser.add_argument("--aiips", required=True)
    parser.add_argument("--skani", required=True)
    parser.add_argument("--skani-args", default="")
    parser.add_argument("--variant", action="append", default=[])
    parser.add_argument("--no-skani", action="store_true")
    parser.add_argument("--repeat", type=int, default=1,
                        help="runs per tool and pair; the fastest is kept")
    args = parser.parse_args()
    variants = [("aiips", [])] + [
        (f"aiips:{label}", shlex.split(extra))
        for label, extra in (v.split("=", 1) for v in args.variant)
    ]
    tools = [(label, lambda a, b, extra=extra: run_aiips(args.aiips, extra, a, b))
             for label, extra in variants]
    if not args.no_skani:
        skani_extra = shlex.split(args.skani_args)
        tools.append(("skani", lambda a, b: run_skani(args.skani, skani_extra, a, b)))
    if os.path.exists(args.out):
        raise SystemExit(f"{args.out} exists; refusing to overwrite results")
    with open(args.truth) as f:
        truth = list(csv.DictReader(f, delimiter="\t"))
    rows = []
    with open(args.out, "w") as out:
        out.write("scenario\tcategory\ttool\ttrue_ani\ttrue_af_a\ttrue_af_b\tani\taf_a\taf_b\twall_s\tmax_rss_kb\n")
        for t in truth:
            for label, run in tools:
                runs = [run(t["a"], t["b"]) for _ in range(args.repeat)]
                ani, af_a, af_b, _, _ = runs[0]
                wall = min(r[3] for r in runs)
                rss = min(r[4] for r in runs)
                row = dict(scenario=t["scenario"], category=t["category"], tool=label,
                           true_ani=float(t["ani"]), true_af_a=float(t["af_a"]),
                           true_af_b=float(t["af_b"]), ani=ani, af_a=af_a, af_b=af_b,
                           wall=wall, rss=rss)
                rows.append(row)
                out.write(f"{t['scenario']}\t{t['category']}\t{label}\t{t['ani']}\t{t['af_a']}\t{t['af_b']}\t"
                          f"{fmt(ani)}\t{fmt(af_a)}\t{fmt(af_b)}\t{wall:.2f}\t{rss}\n")
                out.flush()
                print(f"{t['scenario']:22s} {label:14s} ani {fmt(ani, 4)} (true {float(t['ani']):.4f})  "
                      f"af {fmt(af_a, 3)}/{fmt(af_b, 3)} (true {float(t['af_a']):.3f}/{float(t['af_b']):.3f})  "
                      f"{wall:.2f}s {rss / 1024:.0f}MiB", flush=True)
    summarize(rows, [label for label, _ in tools])


def summarize(rows, tools):
    """Mean/max absolute errors in percentage points, per category and tool."""
    print("\ncategory  tool            n  ani_mae  ani_max  ani_bias  af_mae   af_max   missing  wall_s  rss_MiB")
    for category in sorted({r["category"] for r in rows}) + ["ALL"]:
        for tool in tools:
            sel = [r for r in rows if r["tool"] == tool and category in ("ALL", r["category"])]
            have = [r for r in sel if r["ani"] is not None]
            ani_err = [100 * (r["ani"] - r["true_ani"]) for r in have]
            af_err = [100 * abs(r[k] - r["true_" + k]) for r in have for k in ("af_a", "af_b")]
            mean = lambda xs: statistics.mean(xs) if xs else float("nan")
            print(f"{category:9s} {tool:14s} {len(sel):3d}  {mean([abs(e) for e in ani_err]):7.3f}  "
                  f"{max((abs(e) for e in ani_err), default=float('nan')):7.3f}  {mean(ani_err):8.3f}  "
                  f"{mean(af_err):7.2f}  {max(af_err, default=float('nan')):7.2f}  {len(sel) - len(have):7d}  "
                  f"{mean([r['wall'] for r in sel]):6.2f}  {mean([r['rss'] / 1024 for r in sel]):7.1f}")


if __name__ == "__main__":
    main()
