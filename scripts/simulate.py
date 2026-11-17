#!/usr/bin/env python3
"""Simulates genome pairs with known ANI and alignment fraction (AF).

Usage: scripts/simulate.py OUTDIR [--length N] [--seed S]

Writes the shared A genomes (OUTDIR/a.fa, OUTDIR/a_repeats.fa), one
OUTDIR/<scenario>_b.fa per scenario, and OUTDIR/truth.tsv.
Standard library only; every scenario uses a fixed seed.

Truth definitions (computed from the simulation, not estimated):
- ANI = 1 - substitutions / aligned columns, where aligned columns are the
  homologous base pairs (indel gaps excluded).
- AF of A = A bases with a homologous base in B / len(A); likewise for B.

Scenarios:
- sub:     B = A with uniform substitutions (AF = 1).
- partial: B = a random subset of A's blocks, substituted, interleaved with
           foreign random sequence.
- draft:   substitutions plus small indels; B cut into contigs of a given mean
           length, each in random orientation.
- rearr:   A carries repeat families; B = A substituted, then cut into blocks
           that are shuffled (translocations) and reverse-complemented
           (inversions).
"""

import argparse
import math
import os
import random

BASES = "ACGT"
OTHERS = {b: [c for c in BASES if c != b] for b in BASES}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def random_genome(length, rng):
    return "".join(rng.choices(BASES, k=length))


def geometric(mean, rng, cap):
    """Length >= 1 with the given mean, capped."""
    p = 1.0 / mean
    return min(cap, 1 + int(math.log(1.0 - rng.random()) / math.log(1.0 - p)))


def mutate(seq, p_sub, p_indel, rng, indel_mean=3, indel_cap=30):
    """Applies substitutions and indels at independent per-base rates.

    Returns (mutated, aligned_columns, substitutions)."""
    rate = p_sub + p_indel
    n = len(seq)
    if rate <= 0:
        return seq, n, 0
    log_keep = math.log(1.0 - rate)
    out, pos, subs, deleted = [], 0, 0, 0
    while True:
        nxt = pos + int(math.log(1.0 - rng.random()) / log_keep)
        if nxt >= n:
            out.append(seq[pos:])
            break
        out.append(seq[pos:nxt])
        if rng.random() < p_sub / rate:
            out.append(rng.choice(OTHERS[seq[nxt]]))
            subs += 1
            pos = nxt + 1
        else:
            length = geometric(indel_mean, rng, indel_cap)
            if rng.random() < 0.5:
                # Insertion before base nxt, which is kept unchanged.
                out.append(random_genome(length, rng))
                out.append(seq[nxt])
                pos = nxt + 1
            else:
                length = min(length, n - nxt)
                deleted += length
                pos = nxt + length
    return "".join(out), n - deleted, subs


def cut(seq, mean, rng, minimum=1000):
    """Cuts seq into pieces with exponential lengths (at least `minimum`)."""
    pieces, pos = [], 0
    while pos < len(seq):
        length = max(minimum, int(rng.expovariate(1.0 / mean)))
        if len(seq) - (pos + length) < minimum:
            length = len(seq) - pos
        pieces.append(seq[pos:pos + length])
        pos += length
    return pieces


def add_repeats(genome, rng, families=40, copies=(2, 5), size=(1000, 7000), divergence=0.01):
    """Overwrites random positions with diverged copies of source segments."""
    g = list(genome)
    for _ in range(families):
        length = rng.randint(*size)
        start = rng.randrange(0, len(g) - length)
        source = "".join(g[start:start + length])
        for _ in range(rng.randint(*copies) - 1):
            copy, _, _ = mutate(source, divergence * rng.random(), 0.0, rng)
            if rng.random() < 0.5:
                copy = revcomp(copy)
            target = rng.randrange(0, len(g) - length)
            g[target:target + length] = copy
    return "".join(g)


def write_fasta(path, records):
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80])
                out.write("\n")


class Writer:
    def __init__(self, outdir):
        self.outdir = outdir
        self.truth = open(os.path.join(outdir, "truth.tsv"), "w")
        self.truth.write("scenario\tcategory\tparams\ta\tb\tani\taf_a\taf_b\tlen_a\tlen_b\tcontigs_b\n")

    def emit(self, name, category, params, a_records, b_records, aligned_a, aligned_b, columns, subs):
        a_path = os.path.join(self.outdir, f"{a_records[0][0]}.fa")
        b_path = os.path.join(self.outdir, f"{name}_b.fa")
        if not os.path.exists(a_path):
            write_fasta(a_path, a_records)
        write_fasta(b_path, b_records)
        len_a = sum(len(s) for _, s in a_records)
        len_b = sum(len(s) for _, s in b_records)
        ani = 1.0 - subs / columns if columns else 0.0
        self.truth.write(
            f"{name}\t{category}\t{params}\t{a_path}\t{b_path}\t{ani:.6f}\t"
            f"{aligned_a / len_a:.6f}\t{aligned_b / len_b:.6f}\t{len_a}\t{len_b}\t{len(b_records)}\n"
        )
        self.truth.flush()
        print(f"{name}: ani={ani:.4f} af_a={aligned_a / len_a:.4f} af_b={aligned_b / len_b:.4f}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("outdir")
    parser.add_argument("--length", type=int, default=5_000_000)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    w = Writer(args.outdir)
    L = args.length
    base = random_genome(L, random.Random(args.seed))
    a = [("a", base)]

    # 1. Substitution sweep: ANI from 99.9% to 80%.
    for i, p in enumerate([0.001, 0.01, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20]):
        rng = random.Random(args.seed * 100 + i)
        b, columns, subs = mutate(base, p, 0.0, rng)
        w.emit(f"sub_p{p}", "sub", f"p={p}", a, [("b", b)], L, len(b), columns, subs)

    # 2. Partial overlap: B holds a fraction f of A (in 20-100 kb blocks, kept
    # in order) interleaved with foreign sequence totalling half the rest.
    for i, f in enumerate([0.2, 0.5, 0.9]):
        for j, p in enumerate([0.01, 0.05, 0.10]):
            rng = random.Random(args.seed * 200 + 10 * i + j)
            blocks, pos = [], 0
            while pos < L:
                length = min(rng.randint(20_000, 100_000), L - pos)
                blocks.append((pos, pos + length))
                pos += length
            rng.shuffle(blocks)
            chosen, total = [], 0
            for s, e in blocks:
                if total >= f * L:
                    break
                chosen.append((s, e))
                total += e - s
            chosen.sort()
            foreign_total = int(0.5 * (L - total))
            pieces, columns, subs = [], 0, 0
            for s, e in chosen:
                seq, cols, n_subs = mutate(base[s:e], p, 0.0, rng)
                pieces.append(seq)
                columns += cols
                subs += n_subs
                pieces.append(random_genome(foreign_total // len(chosen), rng))
            b = "".join(pieces)
            w.emit(f"partial_f{f}_p{p}", "partial", f"f={f},p={p}", a, [("b", b)],
                   columns, columns, columns, subs)

    # 3. Draft assemblies: substitutions + indels (rate p/10), B fragmented.
    for i, p in enumerate([0.01, 0.05, 0.10]):
        for j, mean in enumerate([0, 50_000, 10_000, 3_000]):
            rng = random.Random(args.seed * 300 + 10 * i + j)
            b, columns, subs = mutate(base, p, p / 10, rng)
            if mean:
                contigs = [s if rng.random() < 0.5 else revcomp(s) for s in cut(b, mean, rng)]
            else:
                contigs = [b]
            records = [(f"b{k}", s) for k, s in enumerate(contigs)]
            w.emit(f"draft_p{p}_n{mean}", "draft", f"p={p},indel={p / 10},contig_mean={mean}",
                   a, records, columns, columns, columns, subs)

    # 4. Repeats and rearrangements.
    rep = add_repeats(base, random.Random(args.seed + 7))
    a_rep = [("a_repeats", rep)]
    for i, p in enumerate([0.01, 0.05, 0.10]):
        rng = random.Random(args.seed * 400 + i)
        b, columns, subs = mutate(rep, p, 0.0, rng)
        blocks = cut(b, len(b) // 40, rng, minimum=5000)
        rng.shuffle(blocks)
        blocks = [s if rng.random() < 0.5 else revcomp(s) for s in blocks]
        w.emit(f"rearr_p{p}", "rearr", f"p={p},blocks={len(blocks)}", a_rep,
               [("b", "".join(blocks))], columns, columns, columns, subs)


if __name__ == "__main__":
    main()
