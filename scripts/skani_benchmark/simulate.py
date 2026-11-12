#!/usr/bin/env python3
"""Simulate genome pairs with known ANI and alignment fraction (AF).

Genome A is a backbone (a real genome or i.i.d. random sequence). Genome B is
derived from A by, in order:

1. deleting A-only islands (fixed-length segments; lowers AF of A),
2. mutating the remaining sequence: substitutions at rate p (or per-region
   rates for a mosaic), and optional small insertions/deletions,
3. inserting B-only islands of random sequence (lowers AF of B),
4. optional rearrangements (inversions, translocations) and fragmentation
   into contigs (neither changes the truth).

Truth is tracked exactly from the simulated events:
    aligned  = A bases that have a homologous base in B (excludes all
               deleted A bases and all inserted B bases)
    ANI      = 1 - substitutions / aligned          (identity over aligned bases)
    AF_A     = aligned / len(A),  AF_B = aligned / len(B)

Standard library only; deterministic for a given seed.
"""

import argparse
import json
import math
import os
import random

BASES = "ACGT"
OTHERS = {b: [c for c in BASES if c != b] for b in BASES}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def read_fasta(path):
    with open(path) as handle:
        return "".join(line.strip() for line in handle if not line.startswith(">")).upper()


def write_fasta(path, records, width=80):
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i : i + width] + "\n")


def random_sequence(rng, length):
    return "".join(rng.choices(BASES, k=length))


def event_positions(rng, start, end, rate):
    """Positions in [start, end) hit by independent per-base events."""
    if rate <= 0:
        return []
    out, pos = [], start
    log_q = math.log1p(-rate)
    while True:
        pos += int(math.log(1.0 - rng.random()) / log_q)
        if pos >= end:
            return out
        out.append(pos)
        pos += 1


def choose_islands(rng, length, fraction, island_len, forbidden=()):
    """Non-overlapping [start, end) islands covering about `fraction` of length."""
    count = round(fraction * length / island_len)
    taken = sorted(forbidden)
    islands = []
    attempts = 0
    while len(islands) < count:
        attempts += 1
        if attempts > 100 * count + 1000:
            raise RuntimeError("could not place islands; lower the fraction")
        start = rng.randrange(0, length - island_len)
        end = start + island_len
        if any(s < end and start < e for s, e in taken):
            continue
        taken.append((start, end))
        islands.append((start, end))
    return sorted(islands)


def mutate(rng, a, keep, rate_of, indel_rate, indel_mean):
    """Mutate the kept intervals of `a`.

    rate_of(start, end) splits [start, end) into (lo, hi, rate) regions.
    Returns (pieces, aligned, subs, deleted, inserted), where deleted and
    inserted count the bases removed and added by small indels.
    """
    pieces, aligned, subs, deleted, inserted = [], 0, 0, 0, 0
    for start, end in keep:
        events = []
        for lo, hi, rate in rate_of(start, end):
            events += [(p, "s") for p in event_positions(rng, lo, hi, rate)]
        if indel_rate > 0:
            events += [(p, "i") for p in event_positions(rng, start, end, indel_rate)]
        events.sort()
        cursor = start
        for pos, kind in events:
            if pos < cursor:
                continue  # inside a deletion just applied
            pieces.append(a[cursor:pos])
            aligned += pos - cursor
            if kind == "s":
                pieces.append(rng.choice(OTHERS[a[pos]]))
                aligned += 1
                subs += 1
                cursor = pos + 1
            else:
                length = 1 + int(math.log(1.0 - rng.random()) / math.log1p(-1.0 / indel_mean))
                if rng.random() < 0.5:  # deletion of `length` A bases
                    deleted += min(end, pos + length) - pos
                    cursor = min(end, pos + length)
                else:  # insertion of `length` random bases before pos
                    pieces.append(random_sequence(rng, length))
                    inserted += length
                    cursor = pos
        pieces.append(a[cursor:end])
        aligned += end - cursor
    return pieces, aligned, subs, deleted, inserted


def simulate(a, spec, seed):
    rng = random.Random(seed)
    island_len = spec.get("island_len", 20_000)
    # 1. A-only islands.
    removed = choose_islands(rng, len(a), spec.get("a_unique", 0.0), island_len)
    keep, cursor = [], 0
    for s, e in removed:
        keep.append((cursor, s))
        cursor = e
    keep.append((cursor, len(a)))

    # 2. Substitutions (uniform or mosaic) and small indels.
    mosaic = spec.get("mosaic")  # list of rates, cycled over region_len windows
    region_len = spec.get("region_len", 50_000)

    def rate_of(start, end):
        if not mosaic:
            return [(start, end, spec.get("sub_rate", 0.0))]
        out, pos = [], start
        while pos < end:
            region = pos // region_len
            hi = min(end, (region + 1) * region_len)
            out.append((pos, hi, mosaic[region % len(mosaic)]))
            pos = hi
        return out

    pieces, aligned, subs, deleted, inserted = mutate(
        rng, a, keep, rate_of, spec.get("indel_rate", 0.0), spec.get("indel_mean", 3.0)
    )
    b = "".join(pieces)
    island_bases = sum(e - s for s, e in removed)
    # Conservation: every A base is aligned, in an A-only island, or deleted;
    # every B base so far is aligned or inserted.
    assert aligned + island_bases + deleted == len(a)
    assert aligned + inserted == len(b)

    # 3. B-only islands of random sequence.
    b_unique = spec.get("b_unique", 0.0)
    if b_unique > 0:
        total = round(b_unique / (1 - b_unique) * len(b))
        count = max(1, round(total / island_len))
        cuts = sorted(rng.randrange(0, len(b)) for _ in range(count))
        parts, cursor = [], 0
        for cut in cuts:
            parts.append(b[cursor:cut])
            parts.append(random_sequence(rng, island_len))
            cursor = cut
        parts.append(b[cursor:])
        b = "".join(parts)

    # 4. Rearrangements: inversions, then translocations.
    for _ in range(spec.get("inversions", 0)):
        length = rng.randrange(50_000, 200_000)
        start = rng.randrange(0, len(b) - length)
        b = b[:start] + b[start : start + length][::-1].translate(COMPLEMENT) + b[start + length :]
    for _ in range(spec.get("translocations", 0)):
        length = rng.randrange(50_000, 200_000)
        start = rng.randrange(0, len(b) - length)
        segment, rest = b[start : start + length], b[:start] + b[start + length :]
        target = rng.randrange(0, len(rest))
        b = rest[:target] + segment + rest[target:]

    # Fragmentation into contigs with exponential lengths (at least 1 kb).
    mean = spec.get("contig_mean", 0)
    if mean:
        contigs, pos = [], 0
        while pos < len(b):
            length = max(1_000, int(rng.expovariate(1.0 / mean)))
            contigs.append(b[pos : pos + length])
            pos += length
    else:
        contigs = [b]

    assert sum(map(len, contigs)) == len(b)
    truth = {
        "ani": 1 - subs / aligned,
        "af_a": aligned / len(a),
        "af_b": aligned / len(b),
        "aligned": aligned,
        "substitutions": subs,
        "len_a": len(a),
        "len_b": len(b),
        "contigs_b": len(contigs),
    }
    return contigs, truth


def scenarios():
    """The benchmark grid: (name, spec) pairs."""
    grid = []
    for p in (0.001, 0.01, 0.03, 0.05, 0.08, 0.10, 0.13, 0.15, 0.18, 0.20):
        grid.append((f"subs_p{p}", {"sub_rate": p}))
    for p in (0.01, 0.05, 0.10):
        for f in (0.1, 0.3, 0.5, 0.7):
            grid.append((f"af_one_p{p}_f{f}", {"sub_rate": p, "a_unique": f}))
        for f in (0.1, 0.3, 0.5):
            grid.append((f"af_both_p{p}_f{f}", {"sub_rate": p, "a_unique": f, "b_unique": f}))
        grid.append((f"indel_p{p}", {"sub_rate": p, "indel_rate": p / 10}))
    grid.append(("mosaic_0.01_0.10", {"mosaic": [0.01, 0.10]}))
    grid.append(("mosaic_0.02_0.15", {"mosaic": [0.02, 0.15]}))
    for p in (0.01, 0.05):
        for mean in (20_000, 5_000):
            grid.append((f"draft_p{p}_n{mean}", {"sub_rate": p, "contig_mean": mean}))
        grid.append(
            (f"rearranged_p{p}", {"sub_rate": p, "inversions": 5, "translocations": 5})
        )
    return grid


def self_test():
    """Checks the tracked truth against brute force on small genomes."""
    a = random_sequence(random.Random(1), 200_000)
    # Substitutions only: B aligns position by position.
    contigs, truth = simulate(a, {"sub_rate": 0.05}, seed=2)
    b = contigs[0]
    assert len(b) == len(a)
    assert truth["substitutions"] == sum(x != y for x, y in zip(a, b))
    # A-only islands: the kept A bases, in order, align to B.
    rng = random.Random(3)
    removed = choose_islands(rng, len(a), 0.3, 20_000)
    contigs, truth = simulate(a, {"sub_rate": 0.02, "a_unique": 0.3}, seed=3)
    kept = "".join(a[s:e] for s, e in _complement(removed, len(a)))
    b = contigs[0]
    assert truth["aligned"] == len(kept) == len(b)
    assert truth["substitutions"] == sum(x != y for x, y in zip(kept, b))
    assert abs(truth["af_a"] - 0.7) < 1e-9 and truth["af_b"] == 1.0
    # B-only islands: aligned bases unchanged, B longer by the insertions.
    contigs, truth = simulate(a, {"sub_rate": 0.02, "b_unique": 0.3}, seed=4)
    assert truth["aligned"] == len(a) and abs(truth["af_b"] - 0.7) < 0.02
    # Indels: conservation is asserted inside simulate; rearrangements and
    # fragments keep the truth.
    _, with_indels = simulate(a, {"sub_rate": 0.05, "indel_rate": 0.005}, seed=5)
    assert with_indels["aligned"] < len(a)
    _, plain = simulate(a, {"sub_rate": 0.05}, seed=6)
    _, moved = simulate(a, {"sub_rate": 0.05, "inversions": 1, "translocations": 1, "contig_mean": 5_000}, seed=6)
    assert (plain["ani"], plain["af_a"], plain["af_b"]) == (moved["ani"], moved["af_a"], moved["af_b"])
    assert moved["contigs_b"] > 10
    # Mosaic: half the windows at each rate.
    _, mosaic = simulate(a, {"mosaic": [0.01, 0.2]}, seed=7)
    assert abs((1 - mosaic["ani"]) - 0.105) < 0.005
    print("self-test passed")


def _complement(intervals, length):
    out, cursor = [], 0
    for s, e in intervals:
        out.append((cursor, s))
        cursor = e
    out.append((cursor, length))
    return out


def backbone(name, seed, ecoli_path, length):
    if name == "ecoli":
        return read_fasta(ecoli_path)
    return random_sequence(random.Random(seed * 7919 + 1), length)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", default="data/test_against_skani/pairs")
    parser.add_argument(
        "--ecoli",
        default="data/test_against_skani/backbones/ecoli_K12_MG1655_NC_000913.3.fasta",
    )
    parser.add_argument("--backbones", nargs="+", default=["ecoli", "random"])
    parser.add_argument("--seeds", nargs="+", type=int, default=[])
    parser.add_argument("--self-test", action="store_true", help="check the truth tracking and exit")
    parser.add_argument("--only", help="simulate only scenarios whose name contains this")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    for name in args.backbones:
        for seed in args.seeds:
            a = None
            base = os.path.join(args.out, f"{name}_s{seed}")
            os.makedirs(base, exist_ok=True)
            for scenario, spec in scenarios():
                if args.only and args.only not in scenario:
                    continue
                target = os.path.join(base, scenario)
                if os.path.exists(os.path.join(target, "truth.json")):
                    continue
                if a is None:
                    a = backbone(name, seed, args.ecoli, 4_641_652)
                    a_path = os.path.join(base, "a.fa")
                    if not os.path.exists(a_path):
                        write_fasta(a_path, [(f"{name}_backbone", a)])
                # Each scenario has its own stream, independent of the grid order.
                contigs, truth = simulate(a, spec, seed=_seed(name, seed, scenario))
                os.makedirs(target, exist_ok=True)
                write_fasta(
                    os.path.join(target, "b.fa"),
                    [(f"b_contig{i}", c) for i, c in enumerate(contigs)],
                )
                truth.update({"backbone": name, "seed": seed, "scenario": scenario, "spec": spec})
                with open(os.path.join(target, "truth.json"), "w") as out:
                    json.dump(truth, out, indent=1)
                print(f"{name} s{seed} {scenario}: ANI {truth['ani']:.5f} "
                      f"AF {truth['af_a']:.4f}/{truth['af_b']:.4f}", flush=True)


def _seed(name, seed, scenario):
    """A stable per-scenario seed (Python's hash() is salted per process)."""
    text = f"{name}|{seed}|{scenario}"
    value = 0
    for ch in text.encode():
        value = (value * 1_000_003 + ch) % (1 << 61)
    return value


if __name__ == "__main__":
    main()
