//! Base identity from overlaps of decoded sampled k-mers.
//!
//! The versioned fmix64 hash is invertible, so retained hashes reveal the
//! canonical k-mer words without rereading FASTA. Exact unique shared seeds
//! establish constant-diagonal correspondence. The measurement then counts
//! literal base mismatches in overlaps of independently selected features.
//!
//! Under ideal hashing and substitutions independent of the selection hashes,
//! cross probe/anchor pairs sample base identity linearly even when identity
//! varies along a genome. Same-family equal-position, equal-word, or equal
//! anchor-endpoint-pair events are coupled and are excluded. Repeats and sparse
//! mapping still affect variance and coverage; this is not a base alignment.
//! Internal indels can invalidate a projected diagonal. Callers must apply the
//! geometry diagnostics and retain a fallback for unsupported comparisons.
//!
//! Only the two outer exact support seeds of each constant-diagonal run are
//! masked. Shared and unshared interior features are measured together, since
//! excluding all shared seeds would condition the sample on mutations. Optional
//! terminal projection has separate evidence: a caller must check terminal
//! homology before interpreting projected sequence as aligned sequence.
use crate::sketch::aiips::{blocks::Candidates, AiipsSketch};

#[derive(Clone, Copy, Debug)]
pub struct Geometry {
    pub ca: usize,
    pub cb: usize,
    pub reverse: bool,
    pub lo: i64,
    pub hi: i64,
    pub diagonal: i64,
    pub terminal: bool,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct GeometryFlags {
    pub small_steps: u64,
    pub large_steps: u64,
    pub comparisons: u64,
}

/// Read all shared unique-anchor matches, before block validation. Small
/// changes of diagonal suggest indels; islands usually make much larger
/// changes. Counts are diagnostics, not a calibrated indel-rate estimate.
pub fn geometry_flags(candidates: &Candidates) -> GeometryFlags {
    let mut out = GeometryFlags::default();
    for pair in candidates.matches.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        if (a.contig_a, a.contig_b, a.reverse) != (b.contig_a, b.contig_b, b.reverse)
            || b.x <= a.x
            || b.y <= a.y
            || b.x - a.x > 15_000
        {
            continue;
        }
        out.comparisons += 1;
        let step = ((b.y as i64 - b.x as i64) - (a.y as i64 - a.x as i64)).abs();
        if (1..=100).contains(&step) {
            out.small_steps += 1;
        } else if step > 100 {
            out.large_steps += 1;
        }
    }
    out
}

fn dispersion_counts(n: &[u32], k: &[u64]) -> f64 {
    let total = n.iter().map(|&n| u64::from(n)).sum::<u64>();
    let hits = k.iter().sum::<u64>();
    let groups = n.iter().filter(|&&v| v > 0).count();
    if total == 0 || hits < 20 || groups < 2 {
        return 0.0;
    }
    let p = hits as f64 / total as f64;
    if p <= 0.0 || p >= 1.0 {
        return 0.0;
    }
    let pearson = n
        .iter()
        .zip(k)
        .filter(|(n, _)| **n > 0)
        .map(|(&n, &k)| {
            let residual = k as f64 - p * n as f64;
            residual * residual / (p * (1.0 - p) * n as f64)
        })
        .sum::<f64>();
    pearson / (groups - 1) as f64
}

/// Approximate Pearson overdispersion of shared unique-anchor counts among
/// contigs, conditional on sampled unique-anchor counts. Uniform identity is
/// near1 under the working binomial model; mixtures increase it. Repeats,
/// incompleteness, overlapping words and short counts can also increase it.
/// Callers must retain their independent confirmation/fallback safeguards.
pub fn fragment_dispersion(a: &AiipsSketch, b: &AiipsSketch, candidates: &Candidates) -> f64 {
    let side = |s: &AiipsSketch, source_a: bool| {
        if s.contigs().len() < 20 {
            return 0.0;
        }
        // The anchor join already counted unique occurrences per contig.
        // Reusing those counts avoids another genome-wide occurrence scan.
        let n = &candidates.unique_by_contig[usize::from(!source_a)];
        if n.len() != s.contigs().len() {
            return 0.0;
        }
        let mut k = vec![0u64; s.contigs().len()];
        for m in &candidates.matches {
            k[if source_a { m.contig_a } else { m.contig_b }] += 1;
        }
        dispersion_counts(n, &k)
    };
    side(a, true).max(side(b, false))
}

#[derive(Clone, Copy, Debug, Default)]
pub struct RawScore {
    pub bases: u64,
    pub mismatches: u64,
    pub pairs: u64,
    pub pair_sq: u64,
}
impl RawScore {
    fn push(&mut self, n: u64, mismatches: u64) {
        self.bases += n;
        self.mismatches += mismatches;
        self.pairs += 1;
        self.pair_sq += n * n;
    }
    pub fn identity(self) -> Option<f64> {
        (self.bases > 0).then(|| 1.0 - self.mismatches as f64 / self.bases as f64)
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct FullPA {
    pub pairs: u64,
    pub mutated_pairs: u64,
    pub mismatches: u64,
    pub mismatches_sq: u64,
    pub histogram: [u64; 32],
}

#[derive(Clone, Debug, Default)]
pub struct Evidence {
    pub score_pp: RawScore,
    pub score_cross: RawScore,
    pub score_all: RawScore,
    pub terminal_pp: RawScore,
    pub terminal_cross: RawScore,
    pub terminal_all: RawScore,
    pub full_pa: FullPA,
    pub q_geometry: Option<f64>,
    pub probe_a: u64,
    pub probe_b: u64,
    pub probe_shared: u64,
    pub geometry_bases: u64,
    pub geometries: usize,
    pub mapped_bases_a: u64,
    pub mapped_bases_b: u64,
    pub diagonal_conflicts: usize,
}

#[derive(Clone, Copy, Debug)]
struct Seed {
    ca: usize,
    cb: usize,
    reverse: bool,
    x: i64,
    y: i64,
    k: i64,
}
impl Seed {
    fn class(self) -> (usize, usize, bool, i64) {
        (self.ca, self.cb, self.reverse, self.y - self.x)
    }
}

fn probe_locations(s: &AiipsSketch) -> Vec<u64> {
    let p = s.probes();
    let counts = p.counts();
    // Packed (contig,occurrence) avoids a 24-byte Option<(usize,usize)> per
    // hash on 64-bit machines. Builders already bound occurrence counts to u32.
    if p.len() > u32::MAX as usize || s.contigs().len() > u32::MAX as usize {
        return vec![u64::MAX; p.hashes.len()];
    }
    let mut out = vec![u64::MAX; p.hashes.len()];
    for c in 0..s.contigs().len() {
        for i in s.contigs().probes(c) {
            if counts[p.ids[i] as usize] == 1 {
                out[p.ids[i] as usize] = ((c as u64) << 32) | i as u64;
            }
        }
    }
    out
}

/// Exact constant-diagonal runs from anchors, optionally supplemented by unique
/// shared probes. Isolated diagonal classes are discarded before positional
/// scanning so accidental single15 matches do not split a supported true run.
pub fn geometry(
    a: &AiipsSketch,
    b: &AiipsSketch,
    candidates: &Candidates,
    auxiliary: bool,
    project: bool,
) -> Vec<Geometry> {
    let k1 = a.params.anchor_k.get() as i64;
    let k2 = a.params.probe_k.get() as i64;
    let mut seeds = candidates
        .matches
        .iter()
        .map(|m| Seed {
            ca: m.contig_a,
            cb: m.contig_b,
            reverse: m.reverse,
            x: m.x as i64,
            y: m.y as i64,
            k: k1,
        })
        .collect::<Vec<_>>();
    if auxiliary {
        let (pa, pb) = (a.probes(), b.probes());
        let (la, lb) = (probe_locations(a), probe_locations(b));
        let (mut i, mut j) = (0, 0);
        while i < pa.hashes.len() && j < pb.hashes.len() {
            match pa.hashes[i].cmp(&pb.hashes[j]) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    if la[i] != u64::MAX && lb[j] != u64::MAX {
                        let (ca, ia) = ((la[i] >> 32) as usize, la[i] as u32 as usize);
                        let (cb, ib) = ((lb[j] >> 32) as usize, lb[j] as u32 as usize);
                        let reverse = pa.strands.get(ia) != pb.strands.get(ib);
                        let y = if reverse {
                            b.contigs().length(cb) as i64 - pb.positions[ib] as i64 - k2
                        } else {
                            pb.positions[ib] as i64
                        };
                        seeds.push(Seed {
                            ca,
                            cb,
                            reverse,
                            x: pa.positions[ia] as i64,
                            y,
                            k: k2,
                        });
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
    }
    seeds.sort_unstable_by_key(|s| (s.class(), s.x));
    let mut write = 0;
    let mut start = 0;
    while start < seeds.len() {
        let mut end = start + 1;
        while end < seeds.len() && seeds[end].class() == seeds[start].class() {
            end += 1;
        }
        if seeds[end - 1].x - seeds[start].x >= 50 {
            seeds.copy_within(start..end, write);
            write += end - start;
        }
        start = end;
    }
    seeds.truncate(write);
    let mut supported = seeds;
    supported.sort_unstable_by_key(|s| (s.ca, s.cb, s.reverse, s.x, s.k));
    let mut runs = Vec::<(Seed, Seed)>::new();
    let mut start = 0;
    while start < supported.len() {
        let mut end = start + 1;
        while end < supported.len()
            && supported[end].class() == supported[start].class()
            && supported[end].x - supported[end - 1].x <= 15_000
        {
            end += 1;
        }
        if supported[end - 1].x - supported[start].x >= 50 {
            runs.push((supported[start], supported[end - 1]));
        }
        start = end;
    }
    let mut classes_a = vec![Vec::new(); a.contigs().len()];
    let mut classes_b = vec![Vec::new(); b.contigs().len()];
    let mut run_counts = std::collections::BTreeMap::new();
    for &(s, _) in &runs {
        classes_a[s.ca].push(s.class());
        classes_b[s.cb].push(s.class());
    }
    for &(s, _) in &runs {
        *run_counts.entry(s.class()).or_insert(0usize) += 1;
    }
    for v in classes_a.iter_mut().chain(classes_b.iter_mut()) {
        v.sort_unstable();
        v.dedup();
    }
    let mut out = Vec::new();
    for &(first, last) in &runs {
        let diagonal = first.y - first.x;
        let full_lo = 0.max(-diagonal);
        let full_hi = (a.contigs().length(first.ca) as i64)
            .min(b.contigs().length(first.cb) as i64 - diagonal);
        let single_run = run_counts[&first.class()] == 1;
        let extend = project
            && single_run
            && (classes_a[first.ca].len() == 1 || classes_b[first.cb].len() == 1);
        let mut push = |lo: i64, hi: i64, terminal: bool| {
            let lo = lo.max(full_lo);
            let hi = hi.min(full_hi);
            if hi > lo {
                out.push(Geometry {
                    ca: first.ca,
                    cb: first.cb,
                    reverse: first.reverse,
                    lo,
                    hi,
                    diagonal,
                    terminal,
                });
            }
        };
        if extend {
            push(full_lo, first.x, true);
        }
        push(first.x + first.k, last.x, false);
        if extend {
            push(last.x + last.k, full_hi, true);
        }
    }
    // Conservative overlap resolution: independent runs/classes must not
    // count the same source-A bases repeatedly.
    out.sort_unstable_by_key(|g| (g.ca, g.lo, g.hi));
    let mut kept = Vec::<Geometry>::new();
    for mut g in out {
        if let Some(last) = kept.last() {
            if last.ca == g.ca && last.hi > g.lo {
                g.lo = last.hi;
            }
        }
        if g.hi > g.lo {
            kept.push(g);
        }
    }
    kept
}

#[inline]
fn inverse(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0x9cb4_b2f8_1293_37db);
    x ^= x >> 33;
    x = x.wrapping_mul(0x4f74_430c_22a5_4005);
    x ^ (x >> 33)
}
#[inline]
fn mask(k: u32) -> u64 {
    (1u64 << (2 * k)) - 1
}
#[inline]
fn rc(x: u64, k: u32) -> u64 {
    let r = x.reverse_bits();
    let r = ((r & 0x5555_5555_5555_5555) << 1) | ((r & 0xaaaa_aaaa_aaaa_aaaa) >> 1);
    (r >> (64 - 2 * k)) ^ mask(k)
}
#[inline]
fn canonical(x: u64, k: u32) -> u64 {
    x.min(rc(x, k))
}
#[inline]
fn endpoints(x: u64, k: u32, probe_k: u32) -> (u64, u64) {
    let first = canonical(x >> (2 * (k - probe_k)), probe_k);
    let last = canonical(x & mask(probe_k), probe_k);
    (first.min(last), first.max(last))
}

struct Words {
    anchors: Vec<u64>,
    probes: Vec<u64>,
}
impl Words {
    fn new(s: &AiipsSketch) -> Self {
        let decode = |anchor: bool| {
            let (o, k, seed) = if anchor {
                (
                    s.anchors(),
                    s.params.anchor_k.get() as u32,
                    s.params.anchor_seed(),
                )
            } else {
                (
                    s.probes(),
                    s.params.probe_k.get() as u32,
                    s.params.probe_seed(),
                )
            };
            // Decode each dictionary word once, then scatter into occurrence order.
            let dictionary = o
                .hashes
                .iter()
                .map(|&h| inverse(h) ^ seed)
                .collect::<Vec<_>>();
            o.ids
                .iter()
                .enumerate()
                .map(|(i, &id)| {
                    let word = dictionary[id as usize];
                    if o.strands.get(i) {
                        rc(word, k)
                    } else {
                        word
                    }
                })
                .collect()
        };
        Self {
            anchors: decode(true),
            probes: decode(false),
        }
    }
}

struct View<'a> {
    positions: &'a [u32],
    codes: &'a [u64],
    start: usize,
    end: usize,
    reverse: bool,
    len: i64,
    k: u32,
    d: i64,
}
impl View<'_> {
    fn count(&self) -> usize {
        self.end - self.start
    }
    fn entry(&self, i: usize) -> (i64, u64) {
        let j = if self.reverse {
            self.end - 1 - i
        } else {
            self.start + i
        };
        let p = self.positions[j] as i64;
        if self.reverse {
            (
                self.len - p - self.k as i64 - self.d,
                rc(self.codes[j], self.k),
            )
        } else {
            (p - self.d, self.codes[j])
        }
    }
    fn position(&self, i: usize) -> i64 {
        let j = if self.reverse {
            self.end - 1 - i
        } else {
            self.start + i
        };
        let p = self.positions[j] as i64;
        if self.reverse {
            self.len - p - self.k as i64 - self.d
        } else {
            p - self.d
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn view<'a>(
    s: &'a AiipsSketch,
    words: &'a Words,
    c: usize,
    anchor: bool,
    reverse: bool,
    d: i64,
    lo: i64,
    hi: i64,
) -> View<'a> {
    let (o, k, codes, range) = if anchor {
        (
            s.anchors(),
            s.params.anchor_k.get() as u32,
            &words.anchors,
            s.contigs().anchors(c),
        )
    } else {
        (
            s.probes(),
            s.params.probe_k.get() as u32,
            &words.probes,
            s.contigs().probes(c),
        )
    };
    let len = s.contigs().length(c) as i64;
    let (b_lo, b_hi) = if reverse {
        (len - hi - d, len - lo - d)
    } else {
        (lo + d, hi + d)
    };
    let positions = &o.positions[range.clone()];
    let start = range.start + positions.partition_point(|&p| (p as i64) < b_lo - k as i64 + 1);
    let end = range.start + positions.partition_point(|&p| (p as i64) < b_hi);
    View {
        positions: &o.positions,
        codes,
        start,
        end,
        reverse,
        len,
        k,
        d,
    }
}

#[allow(clippy::too_many_arguments)]
fn scan(
    a: &View<'_>,
    b: &View<'_>,
    g: &Geometry,
    same: bool,
    cross: bool,
    probe: bool,
    probe_k: u32,
    e: &mut Evidence,
) {
    let mut j = 0;
    for i in 0..a.count() {
        let (x, wa) = a.entry(i);
        let lo = x.max(g.lo);
        let hi = (x + a.k as i64).min(g.hi);
        if hi <= lo {
            continue;
        }
        while j < b.count() && b.position(j) + b.k as i64 <= lo {
            j += 1;
        }
        let mut t = j;
        while t < b.count() {
            let (y, wb) = b.entry(t);
            if y >= hi {
                break;
            }
            t += 1;
            let left = lo.max(y);
            let right = hi.min(y + b.k as i64);
            if right <= left {
                continue;
            }
            if probe && x == y && lo == x && hi == x + a.k as i64 && wa == wb {
                e.probe_shared += 1;
            }
            if same && x == y {
                continue;
            }
            // Same canonical word has the same sampling event even at a
            // different position. Equal anchor endpoint-pairs share the XOR
            // selection event even if the full anchor identity differs.
            if same && canonical(wa, a.k) == canonical(wb, b.k) {
                continue;
            }
            if same && !probe && endpoints(wa, a.k, probe_k) == endpoints(wb, b.k, probe_k) {
                continue;
            }
            let n = (right - left) as u32;
            let aa = (wa >> (2 * (a.k as i64 - (left - x) - n as i64))) & mask(n);
            let bb = (wb >> (2 * (b.k as i64 - (left - y) - n as i64))) & mask(n);
            let z = aa ^ bb;
            let mismatches = ((z | (z >> 1)) & 0x5555_5555_5555_5555).count_ones() as u64;
            e.score_all.push(n as u64, mismatches);
            if g.terminal {
                e.terminal_all.push(n as u64, mismatches);
            }
            if cross {
                e.score_cross.push(n as u64, mismatches);
                if g.terminal {
                    e.terminal_cross.push(n as u64, mismatches);
                }
                if n == a.k.min(b.k) {
                    e.full_pa.pairs += 1;
                    e.full_pa.mutated_pairs += u64::from(mismatches > 0);
                    e.full_pa.mismatches += mismatches;
                    e.full_pa.mismatches_sq += mismatches * mismatches;
                    e.full_pa.histogram[mismatches as usize] += 1;
                }
            } else if probe {
                e.score_pp.push(n as u64, mismatches);
                if g.terminal {
                    e.terminal_pp.push(n as u64, mismatches);
                }
            }
        }
    }
}

pub fn estimate(a: &AiipsSketch, b: &AiipsSketch, geometry: &[Geometry]) -> Evidence {
    if a.params.hash != "murmur3-fmix64"
        || b.params.hash != "murmur3-fmix64"
        || !a.params.canonical
        || !b.params.canonical
    {
        return Evidence::default();
    }
    let (aw, bw) = (Words::new(a), Words::new(b));
    let mut e = Evidence::default();
    for g in geometry {
        e.geometry_bases += (g.hi - g.lo) as u64;
        e.geometries += 1;
        let ap = view(a, &aw, g.ca, false, false, 0, g.lo, g.hi);
        let aa = view(a, &aw, g.ca, true, false, 0, g.lo, g.hi);
        let bp = view(b, &bw, g.cb, false, g.reverse, g.diagonal, g.lo, g.hi);
        let ba = view(b, &bw, g.cb, true, g.reverse, g.diagonal, g.lo, g.hi);
        for i in 0..ap.count() {
            let x = ap.position(i);
            if x >= g.lo && x + ap.k as i64 <= g.hi {
                e.probe_a += 1;
            }
        }
        for i in 0..bp.count() {
            let x = bp.position(i);
            if x >= g.lo && x + bp.k as i64 <= g.hi {
                e.probe_b += 1;
            }
        }
        let pk = a.params.probe_k.get() as u32;
        scan(&ap, &bp, g, true, false, true, pk, &mut e);
        scan(&ap, &ba, g, false, true, false, pk, &mut e);
        scan(&aa, &bp, g, false, true, false, pk, &mut e);
        scan(&aa, &ba, g, true, false, false, pk, &mut e);
    }
    e.mapped_bases_a = e.geometry_bases;
    let mut intervals = geometry
        .iter()
        .map(|g| {
            let l = b.contigs().length(g.cb) as i64;
            let (s, t) = (g.lo + g.diagonal, g.hi + g.diagonal);
            if g.reverse {
                (g.cb, l - t, l - s)
            } else {
                (g.cb, s, t)
            }
        })
        .collect::<Vec<_>>();
    intervals.sort_unstable();
    let mut current = None;
    for (c, s, t) in intervals {
        match current {
            Some((cc, ss, tt)) if cc == c && s <= tt => current = Some((cc, ss, tt.max(t))),
            _ => {
                if let Some((_, ss, tt)) = current {
                    e.mapped_bases_b += (tt - ss) as u64;
                }
                current = Some((c, s, t));
            }
        }
    }
    if let Some((_, s, t)) = current {
        e.mapped_bases_b += (t - s) as u64;
    }
    let mut classes = geometry
        .iter()
        .map(|g| (g.ca, g.cb, g.reverse, g.diagonal))
        .collect::<Vec<_>>();
    classes.sort_unstable();
    classes.dedup();
    e.diagonal_conflicts = classes
        .windows(2)
        .filter(|v| v[0].0 == v[1].0 && v[0].1 == v[1].1 && v[0].2 == v[1].2)
        .count();
    if e.probe_a + e.probe_b > 0 {
        e.q_geometry = Some(2.0 * e.probe_shared as f64 / (e.probe_a + e.probe_b) as f64);
    }
    e
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mixing_and_reverse_complement() {
        for x in [0, 1, 17, 12345, (1u64 << 42) - 1] {
            assert_eq!(inverse(crate::hash::fmix64(x)), x);
        }
        assert_eq!(rc(1, 2), 11); // AC -> GT
        for k in [3, 15, 21, 31] {
            let x = 1234567 & mask(k);
            assert_eq!(rc(rc(x, k), k), x);
        }
    }
    #[test]
    fn count_dispersion_distinguishes_uniform_and_mixed_identity() {
        let n = vec![100; 40];
        assert_eq!(dispersion_counts(&n, &vec![50; 40]), 0.0);
        let k = (0..40)
            .map(|i| if i % 2 == 0 { 80 } else { 20 })
            .collect::<Vec<_>>();
        assert!(dispersion_counts(&n, &k) > 30.0);
        assert_eq!(dispersion_counts(&[], &[]), 0.0);
    }
    #[test]
    fn diagonal_flags_separate_small_indels_and_large_island_steps() {
        use crate::sketch::aiips::blocks::AnchorMatch;
        let matches = [0, 0, 3, 3, 0, 20_000]
            .iter()
            .enumerate()
            .map(|(i, &d)| {
                let x = (i * 1000) as u32;
                AnchorMatch {
                    anchor_a: i,
                    anchor_b: i,
                    contig_a: 0,
                    contig_b: 0,
                    reverse: false,
                    x,
                    pos_b: x + d,
                    y: x + d,
                }
            })
            .collect();
        let c = Candidates {
            matches,
            blocks: Vec::new(),
            counts: Default::default(),
            unique_by_contig: Default::default(),
        };
        let f = geometry_flags(&c);
        assert_eq!(f.small_steps, 2);
        assert_eq!(f.large_steps, 1);
        assert_eq!(f.comparisons, 5);
    }
    #[test]
    fn coupled_sampling_events_are_excluded_at_different_positions() {
        let g = Geometry {
            ca: 0,
            cb: 0,
            reverse: false,
            lo: 0,
            hi: 50,
            diagonal: 0,
            terminal: false,
        };
        let pa = [0];
        let pb = [1];
        let words = [0];
        let a = View {
            positions: &pa,
            codes: &words,
            start: 0,
            end: 1,
            reverse: false,
            len: 100,
            k: 15,
            d: 0,
        };
        let b = View {
            positions: &pb,
            codes: &words,
            start: 0,
            end: 1,
            reverse: false,
            len: 100,
            k: 15,
            d: 0,
        };
        let mut e = Evidence::default();
        scan(&a, &b, &g, true, false, true, 15, &mut e);
        assert_eq!(e.score_all.bases, 0);
        // A 31-mer endpoints exclude its middlebase. Different middlebases can
        // therefore leave exactly the same XOR selection event.
        let first = 0x1555_5555u64; //15 Cs at the right end,15 As at the left
        let wa = [first];
        let wb = [first | (3u64 << 30)];
        let a = View {
            positions: &pa,
            codes: &wa,
            start: 0,
            end: 1,
            reverse: false,
            len: 100,
            k: 31,
            d: 0,
        };
        let b = View {
            positions: &pb,
            codes: &wb,
            start: 0,
            end: 1,
            reverse: false,
            len: 100,
            k: 31,
            d: 0,
        };
        assert_ne!(canonical(wa[0], 31), canonical(wb[0], 31));
        scan(&a, &b, &g, true, false, false, 15, &mut e);
        assert_eq!(e.score_all.bases, 0);
    }
    fn random(seed: u64, n: usize) -> Vec<u8> {
        crate::hash::SplitMix64::new(seed)
            .take(n)
            .map(|x| b"ACGT"[(x >> 62) as usize])
            .collect()
    }
    fn sketch(records: &[&[u8]]) -> AiipsSketch {
        use std::num::NonZeroUsize;
        let mut fa = Vec::new();
        for seq in records {
            fa.extend_from_slice(b">c\n");
            fa.extend_from_slice(seq);
            fa.push(b'\n');
        }
        let p = crate::sketch::aiips::AiipsParams::new(
            NonZeroUsize::new(21).unwrap(),
            125,
            NonZeroUsize::new(15).unwrap(),
            125,
            42,
        )
        .unwrap();
        AiipsSketch::build(fa.as_slice(), p).unwrap()
    }
    fn mutate(seq: &[u8], seed: u64, mosaic: bool) -> (Vec<u8>, usize) {
        let mut rng = crate::hash::SplitMix64::new(seed);
        let mut n = 0;
        let b = seq
            .iter()
            .enumerate()
            .map(|(i, &c)| {
                let p = if mosaic {
                    if i / 50_000 % 2 == 0 {
                        0.02
                    } else {
                        0.15
                    }
                } else {
                    0.05
                };
                if (rng.next().unwrap() >> 11) as f64 / ((1u64 << 53) as f64) < p {
                    n += 1;
                    let j = b"ACGT".iter().position(|&x| x == c).unwrap();
                    b"ACGT"[(j + 1 + (rng.next().unwrap() % 3) as usize) % 4]
                } else {
                    c
                }
            })
            .collect();
        (b, n)
    }
    #[test]
    fn decoded_words_equal_original_sequences() {
        let seq = random(1, 50_000);
        let s = sketch(&[&seq]);
        let w = Words::new(&s);
        let code = |seq: &[u8]| {
            seq.iter().fold(0, |x, &c| {
                (x << 2) | b"ACGT".iter().position(|&b| b == c).unwrap() as u64
            })
        };
        for (o, words, k) in [(s.anchors(), &w.anchors, 21), (s.probes(), &w.probes, 15)] {
            for (i, &p) in o.positions.iter().enumerate() {
                assert_eq!(words[i], code(&seq[p as usize..p as usize + k]));
            }
        }
    }
    #[test]
    fn unknown_hash_family_has_no_direct_evidence() {
        let seq = random(7, 10_000);
        let mut a = sketch(&[&seq]);
        let b = sketch(&[&seq]);
        a.params.hash = "unknown".into();
        let g = [Geometry {
            ca: 0,
            cb: 0,
            reverse: false,
            lo: 0,
            hi: 10_000,
            diagonal: 0,
            terminal: false,
        }];
        assert_eq!(estimate(&a, &b, &g).score_all.bases, 0);
    }
    #[test]
    fn identical_and_reverse_genomes_have_exact_identity() {
        let seq = random(2, 200_000);
        let rev = seq
            .iter()
            .rev()
            .map(|&c| b"TGCA"[b"ACGT".iter().position(|&x| x == c).unwrap()])
            .collect::<Vec<_>>();
        let a = sketch(&[&seq]);
        let mut scores = Vec::new();
        for bs in [&seq, &rev] {
            let b = sketch(&[bs]);
            let c = crate::sketch::aiips::blocks::candidate_blocks(&a, &b, &Default::default());
            let g = geometry(&a, &b, &c, true, true);
            let e = estimate(&a, &b, &g);
            assert!(e.score_all.bases > 10_000);
            assert_eq!(e.score_all.mismatches, 0);
            assert_eq!(e.q_geometry, Some(1.0));
            assert_eq!(e.full_pa.mutated_pairs, 0);
            scores.push(e.score_all.bases);
        }
        assert_eq!(scores[0], scores[1]);
    }
    #[test]
    fn direct_identity_recovers_mosaic_mean() {
        let seq = random(3, 2_000_000);
        let (bs, n) = mutate(&seq, 4, true);
        let (a, b) = (sketch(&[&seq]), sketch(&[&bs]));
        let c = crate::sketch::aiips::blocks::candidate_blocks(&a, &b, &Default::default());
        let g = geometry(&a, &b, &c, true, true);
        let e = estimate(&a, &b, &g);
        let truth = 1.0 - n as f64 / seq.len() as f64;
        assert_eq!(e.full_pa.histogram.iter().sum::<u64>(), e.full_pa.pairs);
        assert_eq!(
            e.full_pa
                .histogram
                .iter()
                .enumerate()
                .map(|(k, &n)| k as u64 * n)
                .sum::<u64>(),
            e.full_pa.mismatches
        );
        assert_eq!(
            e.full_pa
                .histogram
                .iter()
                .enumerate()
                .map(|(k, &n)| (k * k) as u64 * n)
                .sum::<u64>(),
            e.full_pa.mismatches_sq
        );
        assert_eq!(
            e.full_pa.pairs - e.full_pa.histogram[0],
            e.full_pa.mutated_pairs
        );
        assert!(
            (e.score_cross.identity().unwrap() - truth).abs() < 0.004,
            "{e:?}, truth {truth}"
        );
        assert!(
            (e.score_all.identity().unwrap() - truth).abs() < 0.004,
            "{e:?}, truth {truth}"
        );
        assert!(e.q_geometry.unwrap().powf(1.0 / 15.0) > truth + 0.02);
    }
    #[test]
    fn fragments_are_projected_to_reciprocal_bounds() {
        let seq = random(5, 400_000);
        let (bs, n) = mutate(&seq, 6, false);
        let chunks = bs.chunks(5000).collect::<Vec<_>>();
        let (a, b) = (sketch(&[&seq]), sketch(&chunks));
        let c = crate::sketch::aiips::blocks::candidate_blocks(&a, &b, &Default::default());
        assert!(fragment_dispersion(&a, &b, &c) < 3.0);
        let g = geometry(&a, &b, &c, true, true);
        let e = estimate(&a, &b, &g);
        assert!(e.terminal_all.bases > 0);
        assert!(e.mapped_bases_b as f64 / bs.len() as f64 > 0.95);
        let truth = 1.0 - n as f64 / seq.len() as f64;
        assert!(
            (e.score_all.identity().unwrap() - truth).abs() < 0.006,
            "{e:?}, truth {truth}"
        );
    }
    #[test]
    fn fragmented_mosaics_show_anchor_overdispersion() {
        let seq = random(21, 2_000_000);
        let (bs, _) = mutate(&seq, 22, true);
        let chunks = bs.chunks(2000).collect::<Vec<_>>();
        let (a, b) = (sketch(&[&seq]), sketch(&chunks));
        let c = crate::sketch::aiips::blocks::candidate_blocks(&a, &b, &Default::default());
        assert!(fragment_dispersion(&a, &b, &c) > 2.0);
    }
}
