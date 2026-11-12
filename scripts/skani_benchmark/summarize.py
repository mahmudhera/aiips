#!/usr/bin/env python3
"""Summarize run.py results: accuracy, missing estimates and cost per scenario.

Errors are in percentage points (pp) of ANI or AF. MAE and bias use only
pairs where the tool reported a value; `reported` gives that share. AF error
averages the two directions (AF of A and of B). With --versus A B, pairs run
by both configurations are compared head to head.
"""

import argparse
import csv
import statistics
from collections import defaultdict


def family(scenario):
    return scenario.split("_p")[0].split("_0")[0]


def number(x):
    return None if x in ("NA", "") else float(x)


def load(paths):
    rows = []
    for path in paths:
        with open(path) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                for key in ("true_ani", "true_af_a", "true_af_b", "ani", "af_a", "af_b",
                            "seconds", "max_rss_mb"):
                    row[key] = number(row[key])
                rows.append(row)
    return rows


def stats(rows):
    ani_err = [r["ani"] - r["true_ani"] for r in rows if r["ani"] is not None]
    af_err = [
        (abs(r["af_a"] - r["true_af_a"]) + abs(r["af_b"] - r["true_af_b"])) / 2
        for r in rows if r["af_a"] is not None and r["af_b"] is not None
    ]
    mean = lambda v: statistics.fmean(v) if v else float("nan")
    return {
        "n": len(rows),
        "reported": len(ani_err) / len(rows),
        "ani_mae": 100 * mean([abs(e) for e in ani_err]),
        "ani_bias": 100 * mean(ani_err),
        "af_mae": 100 * mean(af_err),
        "seconds": statistics.median(r["seconds"] for r in rows),
        "rss": max(r["max_rss_mb"] for r in rows),
    }


def table(rows, key):
    groups = defaultdict(list)
    for r in rows:
        groups[key(r)].append(r)
    header = ("| group | n | reported | ANI MAE (pp) | ANI bias (pp) | AF MAE (pp) "
              "| median s | max MB |")
    lines = [header, "|" + "---|" * 8]
    for name in sorted(groups):
        s = stats(groups[name])
        lines.append(
            f"| {' '.join(map(str, name))} | {s['n']} | {s['reported']:.2f} | {s['ani_mae']:.3f} "
            f"| {s['ani_bias']:+.3f} | {s['af_mae']:.2f} | {s['seconds']:.2f} | {s['rss']:.0f} |"
        )
    return "\n".join(lines)


def versus(rows, first, second, key):
    """Head-to-head over pairs both configurations ran."""
    by = defaultdict(dict)
    for r in rows:
        by[(r["backbone"], r["seed"], r["scenario"])][r["label"]] = r
    groups = defaultdict(lambda: defaultdict(int))
    for pair, labels in by.items():
        if first not in labels or second not in labels:
            continue
        a, b = labels[first], labels[second]
        g = groups[key(a)]
        g["pairs"] += 1
        g.setdefault("time_ratio", []).append(a["seconds"] / b["seconds"])
        g.setdefault("rss_ratio", []).append(a["max_rss_mb"] / b["max_rss_mb"])
        for metric, err in (
            ("ani", lambda r: abs(r["ani"] - r["true_ani"]) if r["ani"] is not None else None),
            ("af", lambda r: (abs(r["af_a"] - r["true_af_a"]) + abs(r["af_b"] - r["true_af_b"])) / 2
             if r["af_a"] is not None else None),
        ):
            ea, eb = err(a), err(b)
            if ea is None and eb is None:
                g[f"{metric}_both_missing"] += 1
            elif eb is None or (ea is not None and ea < eb):
                g[f"{metric}_{first}"] += 1
            elif ea is None or eb < ea:
                g[f"{metric}_{second}"] += 1
            else:
                g[f"{metric}_tie"] += 1
    lines = [f"| group | pairs | ANI better: {first} | {second} | tie | both missing "
             f"| AF better: {first} | {second} | tie | time {first}/{second} "
             f"| memory {first}/{second} |", "|" + "---|" * 11]
    for name in sorted(groups):
        g = groups[name]
        lines.append(
            f"| {' '.join(map(str, name))} | {g['pairs']} | {g[f'ani_{first}']} | {g[f'ani_{second}']} "
            f"| {g['ani_tie']} | {g['ani_both_missing']} | {g[f'af_{first}']} | {g[f'af_{second}']} "
            f"| {g['af_tie']} | {statistics.median(g['time_ratio']):.2f} "
            f"| {statistics.median(g['rss_ratio']):.2f} |"
        )
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("results", nargs="+")
    parser.add_argument("--by", choices=["family", "family_all", "scenario"], default="family",
                        help="family_all pools both backbones")
    parser.add_argument("--versus", nargs=2, metavar=("LABEL_A", "LABEL_B"))
    parser.add_argument("--labels", nargs="*", help="only these configurations")
    args = parser.parse_args()
    rows = load(args.results)
    if args.labels:
        rows = [r for r in rows if r["label"] in args.labels]
    if args.by == "family":
        key = lambda r: (r["label"], r["backbone"], family(r["scenario"]))
    elif args.by == "family_all":
        key = lambda r: (r["label"], family(r["scenario"]))
    else:
        key = lambda r: (r["label"], r["backbone"], r["scenario"])
    print(table(rows, key))
    if args.versus:
        drop_label = lambda r: key(r)[1:]
        print()
        print(versus(rows, *args.versus, drop_label))


if __name__ == "__main__":
    main()
