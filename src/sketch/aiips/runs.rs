//! Positional comparison by run extension (`--local-method runs`): shared
//! anchors start runs, and the anchors and probes that follow extend them,
//! across further anchors, until a gap without support ends them. No
//! alignment or dynamic programming is involved.
//!
//! 1. **Seeds.** An anchor hash occurring exactly once on contig a of A and
//!    exactly once on contig b of B seeds a run for that contig pair, even if
//!    it also occurs on other contigs. Orientation is strand_A XOR strand_B;
//!    for reverse runs B positions use the oriented coordinate
//!    y = len_B − (pos_B + k), as for blocks. Seeds are processed per (contig
//!    A, contig B, orientation) group in A order.
//! 2. **Supports.** From the seed, A's anchors and probes are visited in
//!    position order, first to the right, then to the left. A k-mer is a
//!    support when its hash has exactly one occurrence on contig b whose
//!    strand relation matches the run and whose offset y − x lies within
//!    `run_band` of the last support's offset, with y moving in the walk's
//!    direction. B k-mers absent from A neither extend nor stop a run. A
//!    shared anchor reached as a support is consumed and seeds nothing.
//! 3. **Stop.** A walk ends at the contig end, at the end of the group's
//!    previous run (left walks), or when the gap since the last support, in
//!    A or in B, exceeds G. While the run has only its seed, G = `run_max_gap`;
//!    otherwise G = min(`run_max_gap`, `run_gap_factor` · s · walked /
//!    supports), where s = 1 / (1/c1 + 1/c2) is the expected number of A bases
//!    per sampled k-mer and walked counts the k-mers visited up to the run's
//!    latest support, seed included. s · walked / supports estimates the run's
//!    spacing between supports from its hit rate over its supported part;
//!    misses after the latest support do not enlarge G, so a run cannot walk
//!    far past a break and pick up its old offset again. It is a heuristic
//!    stopping rule, not a calibrated test.
//! 4. **Accept.** A run needs `run_min_supports` supports, seed included.
//!    Runs that overlap on either genome are resolved by descending support
//!    count: a weaker run keeps only its supports outside kept runs, split into
//!    pieces where kept runs intervene, and pieces with too few supports are
//!    dropped.
//! 5. **ANI.** Every probe occurrence inside an accepted run's interval counts
//!    on each side, including the run's end supports, probes overlapping
//!    anchors and repeated hashes. A probe matches *positionally*: when its hash
//!    occurs, with the run's strand relation, in the other genome near the
//!    place the run maps it to, that is with an offset within `run_band` of the
//!    offset of either flanking support. Copies elsewhere never match, so a
//!    repeat cannot match through another copy. X_A and X_B count the matched
//!    probes of each side; pooled over runs, q = (X_A + X_B) / (N_A + N_B) and
//!    ANI = q^(1/k2).
//! 6. **AF.** Each run interval, [first support, last support end), is padded
//!    by `run_pad` bases on both sides (default: the probe scale c2, as skani
//!    pads chains by its c), clipped to the contig. AF is the union over runs
//!    divided by the genome's bases.

use std::cmp::{Ordering, Reverse};
use std::collections::BTreeMap;
use std::ops::Range;

use super::local::{covered, LocalConfig, LocalCounts, LocalEstimates, LocalStatus};
use super::occurrences::Occurrences;
use super::AiipsSketch;

/// No B dictionary entry.
const NONE: u32 = u32::MAX;

/// Counters of the run method, for diagnostics and tuning.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RunCounts {
    /// Seed pairs: anchors unique on one contig of each genome.
    pub seeds: u64,
    /// Seeds reached as supports by an earlier run of their group.
    pub seeds_absorbed: u64,
    /// Runs extended from seeds, before acceptance.
    pub runs_built: u64,
    /// Walks (two per run) ended by the gap rule.
    pub stops_gap: u64,
    /// Walks ended by a contig end or by the group's previous run.
    pub stops_end: u64,
    /// K-mers with two or more B copies inside the band: not supports.
    pub ambiguous: u64,
    /// Runs with fewer than `run_min_supports` supports.
    pub runs_short: u64,
    /// Runs that lost supports to better-supported overlapping runs.
    pub runs_trimmed: u64,
    /// Runs left with no piece of `run_min_supports` supports.
    pub runs_dropped: u64,
    /// Accepted runs (pieces), and their supports.
    pub accepted_runs: u64,
    pub supports: u64,
    /// Probe occurrences in accepted runs (N_A, N_B), and those matched
    /// positionally in the other genome (X_A, X_B).
    pub probes_a: u64,
    pub probes_b: u64,
    pub matched_a: u64,
    pub matched_b: u64,
    /// Bases covered by the union of padded run intervals on each genome.
    pub covered_bases_a: u64,
    pub covered_bases_b: u64,
}

/// One support: A start, oriented B start and k-mer length.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Support {
    pub x: u32,
    pub y: u32,
    pub k: u32,
}

/// Supports of one run, ordered by A position.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Run {
    pub contig_a: usize,
    pub contig_b: usize,
    pub reverse: bool,
    pub supports: Vec<Support>,
}

impl Run {
    /// Half-open interval on A: [first support, last support end).
    pub fn span_a(&self) -> (u64, u64) {
        span(self.supports.iter().map(|s| (s.x, s.k)))
    }

    /// Half-open interval on B in forward coordinates.
    pub fn span_b(&self, len_b: u32) -> (u64, u64) {
        let (low, high) = span(self.supports.iter().map(|s| (s.y, s.k)));
        forward(low, high, len_b, self.reverse)
    }
}

/// [min start, max end) of (start, length) pairs.
fn span(items: impl Iterator<Item = (u32, u32)>) -> (u64, u64) {
    items.fold((u64::MAX, 0), |(low, high), (start, k)| {
        (
            low.min(start.into()),
            high.max(u64::from(start) + u64::from(k)),
        )
    })
}

/// An oriented B interval in forward coordinates.
fn forward(low: u64, high: u64, len_b: u32, reverse: bool) -> (u64, u64) {
    if reverse {
        (u64::from(len_b) - high, u64::from(len_b) - low)
    } else {
        (low, high)
    }
}

#[derive(Clone, Copy, Debug)]
struct Seed {
    contig_a: usize,
    contig_b: usize,
    reverse: bool,
    anchor_a: usize,
    x: u32,
    y: u32,
}

impl Seed {
    fn group(&self) -> (usize, usize, bool) {
        (self.contig_a, self.contig_b, self.reverse)
    }
}

/// (hash, occurrence) of every anchor whose hash occurs once on its contig,
/// sorted by hash.
fn contig_unique(s: &AiipsSketch) -> Vec<(u64, usize)> {
    let anchors = s.anchors();
    let mut out = Vec::new();
    let mut ids = Vec::new();
    for contig in 0..s.contigs().len() {
        ids.clear();
        ids.extend(s.contigs().anchors(contig).map(|i| (anchors.ids[i], i)));
        ids.sort_unstable();
        let mut i = 0;
        while i < ids.len() {
            let end = i + ids[i..].iter().take_while(|e| e.0 == ids[i].0).count();
            if end == i + 1 {
                out.push((anchors.hashes[ids[i].0 as usize], ids[i].1));
            }
            i = end;
        }
    }
    out.sort_unstable();
    out
}

/// Seed pairs sorted by (contig A, contig B, orientation, x).
fn seeds(a: &AiipsSketch, b: &AiipsSketch) -> Vec<Seed> {
    let k1 = a.params.anchor_k.get() as u32;
    let (unique_a, unique_b) = (contig_unique(a), contig_unique(b));
    let (anchors_a, anchors_b) = (a.anchors(), b.anchors());
    let mut seeds = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < unique_a.len() && j < unique_b.len() {
        match unique_a[i].0.cmp(&unique_b[j].0) {
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
            Ordering::Equal => {
                let hash = unique_a[i].0;
                let end_a = i + unique_a[i..].iter().take_while(|e| e.0 == hash).count();
                let end_b = j + unique_b[j..].iter().take_while(|e| e.0 == hash).count();
                for &(_, anchor_a) in &unique_a[i..end_a] {
                    for &(_, anchor_b) in &unique_b[j..end_b] {
                        let contig_b = b.contigs().of_anchor(anchor_b);
                        let reverse =
                            anchors_a.strands.get(anchor_a) != anchors_b.strands.get(anchor_b);
                        let pos_b = anchors_b.positions[anchor_b];
                        // Stored anchors end inside their contig: no underflow.
                        let y = if reverse {
                            b.contigs().length(contig_b) - k1 - pos_b
                        } else {
                            pos_b
                        };
                        seeds.push(Seed {
                            contig_a: a.contigs().of_anchor(anchor_a),
                            contig_b,
                            reverse,
                            anchor_a,
                            x: anchors_a.positions[anchor_a],
                            y,
                        });
                    }
                }
                i = end_a;
                j = end_b;
            }
        }
    }
    seeds.sort_unstable_by_key(|s| (s.group(), s.x));
    seeds
}

enum Found {
    None,
    One(i64),
    Many,
}

/// B occurrences of the hashes A shares with B, for one kind of k-mer.
struct Lookup<'s> {
    a: &'s Occurrences,
    b: &'s Occurrences,
    k: i64,
    /// Per A dictionary entry: the B entry with the same hash, or NONE.
    to_b: Vec<u32>,
    /// B occurrence indices grouped by dictionary entry: entry e owns
    /// `order[starts[e]..starts[e + 1]]`, in occurrence (contig, position) order.
    starts: Vec<u32>,
    order: Vec<u32>,
}

impl<'s> Lookup<'s> {
    fn new(a: &'s Occurrences, b: &'s Occurrences, k: u32) -> Self {
        let mut to_b = vec![NONE; a.hashes.len()];
        let (mut i, mut j) = (0, 0);
        while i < a.hashes.len() && j < b.hashes.len() {
            match a.hashes[i].cmp(&b.hashes[j]) {
                Ordering::Less => i += 1,
                Ordering::Greater => j += 1,
                Ordering::Equal => {
                    to_b[i] = j as u32;
                    i += 1;
                    j += 1;
                }
            }
        }
        // Counting sort of B occurrences by entry. Occurrence counts fit u32
        // for genomes below about 5e11 bases at c >= 125.
        let mut starts = vec![0u32; b.hashes.len() + 1];
        for &id in &b.ids {
            starts[id as usize + 1] += 1;
        }
        for e in 0..b.hashes.len() {
            starts[e + 1] += starts[e];
        }
        let mut fill = starts.clone();
        let mut order = vec![0u32; b.ids.len()];
        for (occurrence, &id) in b.ids.iter().enumerate() {
            order[fill[id as usize] as usize] = occurrence as u32;
            fill[id as usize] += 1;
        }
        Self {
            a,
            b,
            k: i64::from(k),
            to_b,
            starts,
            order,
        }
    }

    /// The oriented B start of A occurrence `i`'s hash when exactly one copy
    /// in `range_b` (one contig of length `len_b`) has the run's strand
    /// relation and an oriented start in [lo, hi].
    fn find(
        &self,
        i: usize,
        range_b: &Range<usize>,
        len_b: i64,
        reverse: bool,
        lo: i64,
        hi: i64,
    ) -> Found {
        self.find_in(i, range_b, len_b, reverse, reverse, lo, hi)
    }

    /// As [`Self::find`], with the window [lo, hi] in B's oriented coordinates
    /// when `flip` (else forward), and copies kept when strand_A XOR strand_B
    /// equals `reverse`.
    #[allow(clippy::too_many_arguments)]
    fn find_in(
        &self,
        i: usize,
        range_b: &Range<usize>,
        len_b: i64,
        flip: bool,
        reverse: bool,
        lo: i64,
        hi: i64,
    ) -> Found {
        let entry = self.to_b[self.a.ids[i] as usize];
        if entry == NONE {
            return Found::None;
        }
        let list = &self.order
            [self.starts[entry as usize] as usize..self.starts[entry as usize + 1] as usize];
        let from = list.partition_point(|&o| (o as usize) < range_b.start);
        let to = list.partition_point(|&o| (o as usize) < range_b.end);
        let list = &list[from..to];
        let (low, high) = if flip {
            (len_b - self.k - hi, len_b - self.k - lo)
        } else {
            (lo, hi)
        };
        let position = |o: u32| i64::from(self.b.positions[o as usize]);
        let start = list.partition_point(|&o| position(o) < low);
        let strand_a = self.a.strands.get(i);
        let mut found = Found::None;
        for &o in &list[start..] {
            let p = position(o);
            if p > high {
                break;
            }
            if (self.b.strands.get(o as usize) != strand_a) != reverse {
                continue;
            }
            let y = if flip { len_b - self.k - p } else { p };
            found = match found {
                Found::None => Found::One(y),
                _ => return Found::Many,
            };
        }
        found
    }
}

/// Visits occurrence indices from `next` towards `stop` (exclusive).
struct Cursor {
    next: isize,
    stop: isize,
    step: isize,
}

impl Cursor {
    fn new(range: Range<usize>, right: bool) -> Self {
        let (start, end) = (range.start as isize, range.end as isize);
        if right {
            Self {
                next: start,
                stop: end,
                step: 1,
            }
        } else {
            Self {
                next: end - 1,
                stop: start - 1,
                step: -1,
            }
        }
    }

    fn peek(&self) -> Option<usize> {
        (self.next != self.stop).then_some(self.next as usize)
    }

    fn advance(&mut self) {
        self.next += self.step;
    }
}

/// The run's hit-rate bookkeeping, shared by both walks.
struct Progress {
    /// K-mers visited, seed included.
    walked: u64,
    supports: u64,
    /// `walked` when the latest support was found: the hit rate covers the
    /// supported part of the run only.
    walked_at_support: u64,
}

struct Walker<'s> {
    a: &'s AiipsSketch,
    b: &'s AiipsSketch,
    anchors: Lookup<'s>,
    probes: Lookup<'s>,
    k1: u32,
    k2: u32,
    band: i64,
    gap_factor: f64,
    max_gap: f64,
    /// Expected A bases per sampled k-mer, anchors and probes together.
    spacing: f64,
}

impl<'s> Walker<'s> {
    fn new(a: &'s AiipsSketch, b: &'s AiipsSketch, config: &LocalConfig) -> Self {
        let (k1, k2) = (
            a.params.anchor_k.get() as u32,
            a.params.probe_k.get() as u32,
        );
        let rate = 1.0 / a.params.anchor_scaled as f64 + 1.0 / a.params.probe_scaled as f64;
        Self {
            a,
            b,
            anchors: Lookup::new(a.anchors(), b.anchors(), k1),
            probes: Lookup::new(a.probes(), b.probes(), k2),
            k1,
            k2,
            band: i64::from(config.run_band),
            gap_factor: config.run_gap_factor,
            max_gap: f64::from(config.run_max_gap),
            spacing: 1.0 / rate,
        }
    }

    /// The largest gap a run may cross without support.
    fn gap(&self, p: &Progress) -> f64 {
        if p.supports <= 1 {
            self.max_gap
        } else {
            self.max_gap.min(
                self.gap_factor * self.spacing * p.walked_at_support as f64 / p.supports as f64,
            )
        }
    }

    /// The run grown from `seed`. Left walks stop at A starts at or before
    /// `limit`. Anchor occurrences that became supports go to `absorbed`.
    fn extend(
        &self,
        seed: &Seed,
        limit: Option<u32>,
        absorbed: &mut Vec<usize>,
        counts: &mut RunCounts,
    ) -> Run {
        let start = Support {
            x: seed.x,
            y: seed.y,
            k: self.k1,
        };
        let mut progress = Progress {
            walked: 1,
            supports: 1,
            walked_at_support: 1,
        };
        let mut right = Vec::new();
        self.walk(
            seed,
            true,
            None,
            &mut right,
            absorbed,
            &mut progress,
            counts,
        );
        let mut left = Vec::new();
        self.walk(
            seed,
            false,
            limit,
            &mut left,
            absorbed,
            &mut progress,
            counts,
        );
        left.reverse();
        left.push(start);
        left.extend(right);
        Run {
            contig_a: seed.contig_a,
            contig_b: seed.contig_b,
            reverse: seed.reverse,
            supports: left,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn walk(
        &self,
        seed: &Seed,
        right: bool,
        limit: Option<u32>,
        out: &mut Vec<Support>,
        absorbed: &mut Vec<usize>,
        progress: &mut Progress,
        counts: &mut RunCounts,
    ) {
        let (a, b) = (self.a, self.b);
        let anchors = a.contigs().anchors(seed.contig_a);
        let probes = a.contigs().probes(seed.contig_a);
        // The first probe at or after the seed starts the right walk.
        let split = a.owned_probes(seed.anchor_a).start;
        let (mut next_anchor, mut next_probe) = if right {
            (
                Cursor::new(seed.anchor_a + 1..anchors.end, true),
                Cursor::new(split..probes.end, true),
            )
        } else {
            (
                Cursor::new(anchors.start..seed.anchor_a, false),
                Cursor::new(probes.start..split, false),
            )
        };
        let (anchor_positions, probe_positions) = (&a.anchors().positions, &a.probes().positions);
        let ranges_b = (
            b.contigs().anchors(seed.contig_b),
            b.contigs().probes(seed.contig_b),
        );
        let len_b = i64::from(b.contigs().length(seed.contig_b));
        let mut frontier = Support {
            x: seed.x,
            y: seed.y,
            k: self.k1,
        };
        loop {
            let anchor = next_anchor.peek().map(|i| anchor_positions[i]);
            let probe = next_probe.peek().map(|i| probe_positions[i]);
            // Anchors first on equal positions.
            let take_anchor = match (anchor, probe) {
                (None, None) => {
                    counts.stops_end += 1;
                    return;
                }
                (Some(_), None) => true,
                (None, Some(_)) => false,
                (Some(pa), Some(pp)) => {
                    if right {
                        pa <= pp
                    } else {
                        pa >= pp
                    }
                }
            };
            let (x, i) = if take_anchor {
                (anchor.unwrap(), next_anchor.peek().unwrap())
            } else {
                (probe.unwrap(), next_probe.peek().unwrap())
            };
            if !right && limit.is_some_and(|limit| x <= limit) {
                counts.stops_end += 1;
                return;
            }
            let gap = self.gap(progress);
            if f64::from(x.abs_diff(frontier.x)) > gap {
                counts.stops_gap += 1;
                return;
            }
            if take_anchor {
                next_anchor.advance();
            } else {
                next_probe.advance();
            }
            progress.walked += 1;
            let center = i64::from(x) + i64::from(frontier.y) - i64::from(frontier.x);
            let fy = i64::from(frontier.y);
            let (lo, hi) = if right {
                ((center - self.band).max(fy), center + self.band)
            } else {
                (center - self.band, (center + self.band).min(fy))
            };
            if lo > hi {
                continue;
            }
            let (lookup, range_b, k) = if take_anchor {
                (&self.anchors, &ranges_b.0, self.k1)
            } else {
                (&self.probes, &ranges_b.1, self.k2)
            };
            match lookup.find(i, range_b, len_b, seed.reverse, lo, hi) {
                Found::One(y) => {
                    if (y - fy).unsigned_abs() as f64 > gap {
                        counts.stops_gap += 1;
                        return;
                    }
                    frontier = Support { x, y: y as u32, k };
                    out.push(frontier);
                    progress.supports += 1;
                    progress.walked_at_support = progress.walked;
                    if take_anchor {
                        absorbed.push(i);
                    }
                }
                Found::Many => counts.ambiguous += 1,
                Found::None => {}
            }
        }
    }
}

/// Disjoint covered intervals on one contig: start -> end.
#[derive(Default)]
struct Taken(BTreeMap<u64, u64>);

impl Taken {
    fn free(&self, start: u64, end: u64) -> bool {
        // Intervals are disjoint, so only the last one starting before `end`
        // can reach past `start`.
        self.0
            .range(..end)
            .next_back()
            .is_none_or(|(_, &e)| e <= start)
    }

    fn insert(&mut self, mut start: u64, mut end: u64) {
        while let Some((&s, &e)) = self.0.range(..=end).next_back() {
            if e < start {
                break;
            }
            start = start.min(s);
            end = end.max(e);
            self.0.remove(&s);
        }
        self.0.insert(start, end);
    }
}

/// Keeps runs by descending support count; a weaker run keeps the pieces of
/// its supports that lie outside kept runs on both genomes.
fn resolve_overlaps(
    b: &AiipsSketch,
    runs: Vec<Run>,
    min_supports: usize,
    counts: &mut RunCounts,
) -> Vec<Run> {
    let mut order: Vec<usize> = (0..runs.len()).collect();
    order.sort_by_key(|&i| (Reverse(runs[i].supports.len()), i));
    let mut taken_a: BTreeMap<usize, Taken> = BTreeMap::new();
    let mut taken_b: BTreeMap<usize, Taken> = BTreeMap::new();
    let mut kept = Vec::new();
    for i in order {
        let run = &runs[i];
        let len_b = b.contigs().length(run.contig_b);
        let free = |x0: u32, x1: u64, y0: u32, y1: u64| {
            let (b0, b1) = forward(y0.into(), y1, len_b, run.reverse);
            taken_a
                .get(&run.contig_a)
                .is_none_or(|t| t.free(x0.into(), x1))
                && taken_b.get(&run.contig_b).is_none_or(|t| t.free(b0, b1))
        };
        let mut pieces: Vec<Vec<Support>> = Vec::new();
        let mut piece: Vec<Support> = Vec::new();
        for &s in &run.supports {
            let end_x = u64::from(s.x) + u64::from(s.k);
            let end_y = u64::from(s.y) + u64::from(s.k);
            if !free(s.x, end_x, s.y, end_y) {
                pieces.push(std::mem::take(&mut piece));
                continue;
            }
            if let Some(last) = piece.last() {
                // A kept run between two supports splits the piece.
                let last_x = u64::from(last.x) + u64::from(last.k);
                let last_y = u64::from(last.y) + u64::from(last.k);
                if !free(
                    last.x,
                    end_x.max(last_x),
                    last.y.min(s.y),
                    end_y.max(last_y),
                ) {
                    pieces.push(std::mem::take(&mut piece));
                }
            }
            piece.push(s);
        }
        pieces.push(piece);
        let whole = pieces.len() == 1 && pieces[0].len() == run.supports.len();
        let pieces: Vec<Run> = pieces
            .into_iter()
            .filter(|p| p.len() >= min_supports)
            .map(|supports| Run {
                supports,
                ..run.clone()
            })
            .collect();
        if pieces.is_empty() {
            counts.runs_dropped += 1;
            continue;
        }
        if !whole {
            counts.runs_trimmed += 1;
        }
        for piece in pieces {
            let (a0, a1) = piece.span_a();
            let (b0, b1) = piece.span_b(len_b);
            taken_a.entry(piece.contig_a).or_default().insert(a0, a1);
            taken_b.entry(piece.contig_b).or_default().insert(b0, b1);
            kept.push(piece);
        }
    }
    kept.sort_by_key(|r| (r.contig_a, r.contig_b, r.reverse, r.supports[0].x));
    kept
}

/// Accepted runs of `a` against `b`, with counters.
#[cfg(test)]
pub(crate) fn runs(
    a: &AiipsSketch,
    b: &AiipsSketch,
    config: &LocalConfig,
) -> (Vec<Run>, RunCounts) {
    runs_with(&Walker::new(a, b, config), config)
}

fn runs_with(walker: &Walker, config: &LocalConfig) -> (Vec<Run>, RunCounts) {
    let (a, b) = (walker.a, walker.b);
    let mut counts = RunCounts::default();
    let seeds = seeds(a, b);
    counts.seeds = seeds.len() as u64;
    let min_supports = config.run_min_supports as usize;
    let mut built = Vec::new();
    let mut absorbed_anchors = Vec::new();
    let mut start = 0;
    while start < seeds.len() {
        let group = seeds[start].group();
        let end = start
            + seeds[start..]
                .iter()
                .take_while(|s| s.group() == group)
                .count();
        let mut absorbed = vec![false; end - start];
        // A start of the group's previous accepted-size run's last support.
        let mut previous: Option<u32> = None;
        for s in start..end {
            if absorbed[s - start] {
                counts.seeds_absorbed += 1;
                continue;
            }
            absorbed_anchors.clear();
            let run = walker.extend(&seeds[s], previous, &mut absorbed_anchors, &mut counts);
            counts.runs_built += 1;
            for &anchor in &absorbed_anchors {
                let x = a.anchors().positions[anchor];
                let at = start + seeds[start..end].partition_point(|t| t.x < x);
                if at < end && seeds[at].anchor_a == anchor && at > s {
                    absorbed[at - start] = true;
                }
            }
            if run.supports.len() >= min_supports {
                let last = run.supports[run.supports.len() - 1].x;
                previous = Some(previous.map_or(last, |p| p.max(last)));
                built.push(run);
            } else {
                counts.runs_short += 1;
            }
        }
        start = end;
    }
    let kept = resolve_overlaps(b, built, min_supports, &mut counts);
    counts.accepted_runs = kept.len() as u64;
    counts.supports = kept.iter().map(|r| r.supports.len() as u64).sum();
    (kept, counts)
}

/// Offsets (oriented B start minus A start) of the supports flanking a
/// position, found by `key` (A start, or oriented B start, of each support).
fn flanking(supports: &[Support], at: i64, key: impl Fn(&Support) -> i64) -> (i64, i64) {
    let next = supports.partition_point(|s| key(s) <= at);
    let offset = |s: &Support| i64::from(s.y) - i64::from(s.x);
    let before = supports[next.saturating_sub(1)];
    let after = supports[next.min(supports.len() - 1)];
    let (d0, d1) = (offset(&before), offset(&after));
    (d0.min(d1), d0.max(d1))
}

/// Probe occurrences of `contig` lying wholly inside [start, end).
fn probes_within(s: &AiipsSketch, contig: usize, start: u64, end: u64) -> Range<usize> {
    let range = s.contigs().probes(contig);
    let k = s.params.probe_k.get() as u64;
    let positions = &s.probes().positions[range.clone()];
    let from = positions.partition_point(|&p| u64::from(p) < start);
    let to = positions.partition_point(|&p| u64::from(p) + k <= end);
    range.start + from..range.start + to.max(from)
}

/// Positional ANI and AF from accepted runs.
pub fn estimate(a: &AiipsSketch, b: &AiipsSketch, config: &LocalConfig) -> LocalEstimates {
    let k2 = a.params.probe_k.get() as u32;
    let walker = Walker::new(a, b, config);
    let (kept, mut counts) = runs_with(&walker, config);
    let pad = config
        .run_pad
        .map_or(a.params.probe_scaled, u64::from)
        .min(u64::from(u32::MAX));
    let (a_to_b, b_to_a) = (&walker.probes, Lookup::new(b.probes(), a.probes(), k2));
    let (k, band) = (i64::from(k2), i64::from(config.run_band));
    let (mut covered_a, mut covered_b) = (Vec::new(), Vec::new());
    for run in &kept {
        let (a0, a1) = run.span_a();
        let (b0, b1) = run.span_b(b.contigs().length(run.contig_b));
        let len_a = u64::from(a.contigs().length(run.contig_a));
        let len_b = u64::from(b.contigs().length(run.contig_b));
        let (range_a, range_b) = (
            a.contigs().probes(run.contig_a),
            b.contigs().probes(run.contig_b),
        );
        // A probe at x maps to oriented B starts x + d, d between the
        // flanking supports' offsets.
        for i in probes_within(a, run.contig_a, a0, a1) {
            let x = i64::from(a.probes().positions[i]);
            let (low, high) = flanking(&run.supports, x, |s| i64::from(s.x));
            let found = a_to_b.find(
                i,
                &range_b,
                len_b as i64,
                run.reverse,
                x + low - band,
                x + high + band,
            );
            counts.probes_a += 1;
            counts.matched_a += u64::from(!matches!(found, Found::None));
        }
        // A B probe at oriented y maps to A starts y - d.
        for j in probes_within(b, run.contig_b, b0, b1) {
            let p = i64::from(b.probes().positions[j]);
            let y = if run.reverse { len_b as i64 - k - p } else { p };
            let (low, high) = flanking(&run.supports, y, |s| i64::from(s.y));
            let found = b_to_a.find_in(
                j,
                &range_a,
                len_a as i64,
                false,
                run.reverse,
                y - high - band,
                y - low + band,
            );
            counts.probes_b += 1;
            counts.matched_b += u64::from(!matches!(found, Found::None));
        }
        covered_a.push((run.contig_a, a0.saturating_sub(pad), (a1 + pad).min(len_a)));
        covered_b.push((run.contig_b, b0.saturating_sub(pad), (b1 + pad).min(len_b)));
    }
    counts.covered_bases_a = covered(covered_a);
    counts.covered_bases_b = covered(covered_b);
    let local = LocalCounts {
        covered_bases_a: counts.covered_bases_a,
        covered_bases_b: counts.covered_bases_b,
        ..LocalCounts::default()
    };
    if kept.is_empty() {
        return LocalEstimates {
            ani: None,
            q: None,
            af_a: None,
            af_b: None,
            ani_status: LocalStatus::NoSupportedMapping,
            af_status: LocalStatus::NoSupportedMapping,
            counts: local,
            runs: Some(counts),
        };
    }
    let pooled = counts.probes_a + counts.probes_b;
    let shared = counts.matched_a + counts.matched_b;
    let q =
        (pooled >= config.min_scoring_probes && shared > 0).then(|| shared as f64 / pooled as f64);
    LocalEstimates {
        ani: q.map(|q| q.powf(1.0 / f64::from(k2))),
        q,
        af_a: Some(counts.covered_bases_a as f64 / a.summary.bases as f64),
        af_b: Some(counts.covered_bases_b as f64 / b.summary.bases as f64),
        ani_status: if q.is_some() {
            LocalStatus::Estimated
        } else {
            LocalStatus::InsufficientEvidence
        },
        af_status: LocalStatus::Estimated,
        counts: local,
        runs: Some(counts),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sketch::aiips::blocks::{candidate_blocks, tests::dense};
    use crate::sketch::aiips::local::LocalMethod;
    use crate::sketch::aiips::tests::{build, fasta, params, random_bases, reverse_complement};
    use crate::sketch::aiips::AiipsParams;

    /// Default k-mer sizes and scales (k1 = 21, c1 = 125, k2 = 15, c2 = 125).
    fn standard() -> AiipsParams {
        params(21, 125, 15, 125, 42)
    }

    fn config() -> LocalConfig {
        LocalConfig {
            method: LocalMethod::Runs,
            ..LocalConfig::default()
        }
    }

    fn sketches(a: &[&[u8]], b: &[&[u8]], p: &AiipsParams) -> (AiipsSketch, AiipsSketch) {
        (build(&fasta(a), p.clone()), build(&fasta(b), p.clone()))
    }

    fn compare(a: &[&[u8]], b: &[&[u8]], p: &AiipsParams, c: &LocalConfig) -> LocalEstimates {
        let (a, b) = sketches(a, b, p);
        estimate(&a, &b, c)
    }

    /// `count` substitutions at distinct pseudo-random positions.
    fn substitute(sequence: &[u8], count: usize, seed: u64) -> Vec<u8> {
        let mut out = sequence.to_vec();
        let mut hit = vec![false; out.len()];
        let mut state = seed;
        let mut done = 0;
        while done < count {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let at = (state % out.len() as u64) as usize;
            if !hit[at] {
                hit[at] = true;
                out[at] = match out[at] {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                };
                done += 1;
            }
        }
        out
    }

    fn concat(parts: &[&[u8]]) -> Vec<u8> {
        parts.concat()
    }

    #[test]
    fn identical_genomes_give_one_run_with_full_identity_and_coverage() {
        let genome = random_bases(1, 100_000);
        let (a, b) = sketches(&[&genome], &[&genome], &standard());
        let e = estimate(&a, &b, &config());
        let r = e.runs.unwrap();
        assert_eq!(r.accepted_runs, 1);
        assert_eq!(r.seeds_absorbed, r.seeds - 1);
        assert_eq!((e.ani, e.q), (Some(1.0), Some(1.0)));
        assert_eq!((e.af_a, e.af_b), (Some(1.0), Some(1.0)));
        assert_eq!(
            (r.matched_a, r.probes_a),
            (r.probes_a, a.probes().len() as u64)
        );
        assert_eq!(e.ani_status, LocalStatus::Estimated);
    }

    #[test]
    fn substitutions_give_identity_within_sampling_error() {
        let genome = random_bases(2, 300_000);
        for (count, tolerance) in [(3_000, 0.002), (15_000, 0.004), (30_000, 0.008)] {
            let mutated = substitute(&genome, count, 9);
            let truth = 1.0 - count as f64 / genome.len() as f64;
            let e = compare(&[&genome], &[&mutated], &standard(), &config());
            let ani = e.ani.unwrap();
            assert!((ani - truth).abs() < tolerance, "{count}: {ani} vs {truth}");
            assert!(e.af_a.unwrap() > 0.995 && e.af_b.unwrap() > 0.995, "{e:?}");
        }
    }

    #[test]
    fn two_sided_novel_island_breaks_the_run_and_lowers_af() {
        let (left, right) = (random_bases(3, 60_000), random_bases(4, 60_000));
        let a = concat(&[&left, &random_bases(5, 5_000), &right]);
        let b = concat(&[&left, &random_bases(6, 5_000), &right]);
        let (sa, sb) = sketches(&[&a], &[&b], &standard());
        let (kept, _) = runs(&sa, &sb, &config());
        assert!(kept.len() >= 2);
        for run in &kept {
            let (start, end) = run.span_a();
            assert!(end <= 60_000 || start >= 65_000, "{:?}", (start, end));
        }
        let e = estimate(&sa, &sb, &config());
        for af in [e.af_a.unwrap(), e.af_b.unwrap()] {
            assert!((af - 120.0 / 125.0).abs() < 0.005, "{af}");
        }
        assert_eq!(e.ani, Some(1.0));
    }

    #[test]
    fn one_sided_insertion_breaks_the_run_through_the_band() {
        let (left, right) = (random_bases(7, 60_000), random_bases(8, 60_000));
        let a = concat(&[&left, &right]);
        let b = concat(&[&left, &random_bases(9, 3_000), &right]);
        let (sa, sb) = sketches(&[&a], &[&b], &standard());
        let (kept, _) = runs(&sa, &sb, &config());
        assert!(kept.len() >= 2);
        // No run crosses the insertion: offsets are 0 before it, 3000 after.
        for run in &kept {
            let (start, end) = run.span_a();
            let offsets: Vec<i64> = run
                .supports
                .iter()
                .map(|s| i64::from(s.y) - i64::from(s.x))
                .collect();
            assert!(end <= 60_000 || start >= 60_000, "{:?}", (start, end));
            assert!(offsets.iter().all(|&d| d == offsets[0]));
        }
        let e = estimate(&sa, &sb, &config());
        assert!(e.af_a.unwrap() > 0.995, "{e:?}");
        let af_b = e.af_b.unwrap();
        assert!((af_b - 120.0 / 123.0).abs() < 0.005, "{af_b}");
    }

    #[test]
    fn small_indels_stay_in_one_run_and_a_long_deletion_splits_it() {
        let genome = random_bases(10, 100_000);
        // One-base deletions every 997 bases shift the offset by about 100
        // over the genome, but by at most 1 between neighbouring supports.
        let small: Vec<u8> = genome
            .iter()
            .enumerate()
            .filter(|(i, _)| i % 997 != 500)
            .map(|(_, &base)| base)
            .collect();
        let e = compare(&[&genome], &[&small], &standard(), &config());
        assert_eq!(e.runs.unwrap().accepted_runs, 1, "{e:?}");
        assert!(e.af_a.unwrap() > 0.995);
        assert!(e.runs.unwrap().stops_gap <= 1);
        let long = concat(&[&genome[..50_000], &genome[50_500..]]);
        let (sa, sb) = sketches(&[&genome], &[&long], &standard());
        let (kept, _) = runs(&sa, &sb, &config());
        assert!(kept.len() >= 2);
        for run in &kept {
            let (start, end) = run.span_a();
            assert!(end <= 50_000 || start >= 50_500, "{:?}", (start, end));
        }
        let e = estimate(&sa, &sb, &config());
        let af_a = e.af_a.unwrap();
        assert!(af_a < 0.998 && af_a > 0.99, "{af_a}");
    }

    #[test]
    fn reverse_complemented_b_gives_the_same_estimates() {
        let genome = random_bases(11, 100_000);
        let mutated = substitute(&genome, 5_000, 3);
        let forward = compare(&[&genome], &[&mutated], &standard(), &config());
        let reverse = compare(
            &[&genome],
            &[&reverse_complement(&mutated)],
            &standard(),
            &config(),
        );
        assert_eq!(forward, reverse);
        assert!(forward.runs.unwrap().accepted_runs >= 1);
    }

    #[test]
    fn runs_stay_within_one_contig_pair() {
        let genome = random_bases(12, 100_000);
        let (a, b) = sketches(
            &[&genome],
            &[&genome[..50_000], &genome[50_000..]],
            &standard(),
        );
        let (kept, counts) = runs(&a, &b, &config());
        assert!(counts.accepted_runs >= 2);
        for run in &kept {
            let (start, end) = run.span_b(b.contigs().length(run.contig_b));
            assert!(start < end && end <= u64::from(b.contigs().length(run.contig_b)));
        }
        let contigs: std::collections::BTreeSet<usize> = kept.iter().map(|r| r.contig_b).collect();
        assert_eq!(contigs.len(), 2);
        let e = estimate(&a, &b, &config());
        assert!(e.af_a.unwrap() > 0.99 && e.af_b.unwrap() > 0.99, "{e:?}");
    }

    #[test]
    fn contig_unique_anchors_seed_on_each_contig_and_overlaps_are_trimmed() {
        let genome = random_bases(13, 80_000);
        let duplicate = genome[10_000..30_000].to_vec();
        let (a, b) = sketches(&[&genome], &[&genome, &duplicate], &standard());
        let genome_unique = candidate_blocks(&a, &b, &LocalConfig::default())
            .counts
            .shared_anchors;
        let e = estimate(&a, &b, &config());
        let r = e.runs.unwrap();
        // Anchors of the duplicated stretch seed against both B contigs.
        assert!(
            r.seeds > genome_unique + 50,
            "{} vs {genome_unique}",
            r.seeds
        );
        assert_eq!(r.accepted_runs, 1);
        assert_eq!(r.runs_dropped, 1);
        assert!(e.af_a.unwrap() > 0.995 && r.covered_bases_a <= 80_000);
        let af_b = e.af_b.unwrap();
        assert!((af_b - 0.8).abs() < 0.005, "{af_b}");
    }

    #[test]
    fn left_walk_covers_sequence_before_the_first_shared_anchor() {
        let genome = random_bases(14, 30_000);
        let sparse = params(21, 5_000, 15, 125, 42);
        let (a, b) = sketches(&[&genome], &[&genome], &sparse);
        let first_anchor = a.anchors().positions[0];
        assert!(first_anchor > 1_000, "{first_anchor}");
        let (kept, _) = runs(&a, &b, &config());
        assert_eq!(kept.len(), 1);
        assert!(kept[0].span_a().0 < 200, "{:?}", kept[0].span_a());
        assert_eq!(estimate(&a, &b, &config()).af_a, Some(1.0));
    }

    #[test]
    fn runs_need_the_configured_supports() {
        let shared = random_bases(15, 600);
        let a = concat(&[
            &random_bases(16, 20_000),
            &shared,
            &random_bases(17, 20_000),
        ]);
        let b = concat(&[
            &random_bases(18, 30_000),
            &shared,
            &random_bases(19, 10_000),
        ]);
        let (a, b) = sketches(&[&a], &[&b], &dense());
        let two = LocalConfig {
            run_min_supports: 2,
            ..config()
        };
        let e = estimate(&a, &b, &two);
        let r = e.runs.unwrap();
        assert_eq!(r.accepted_runs, 1);
        let enough = LocalConfig {
            run_min_supports: r.supports as u32,
            ..config()
        };
        assert_eq!(estimate(&a, &b, &enough).runs.unwrap().accepted_runs, 1);
        let too_many = LocalConfig {
            run_min_supports: r.supports as u32 + 1,
            ..config()
        };
        let e = estimate(&a, &b, &too_many);
        assert_eq!(e.runs.unwrap().accepted_runs, 0);
        assert!(e.runs.unwrap().runs_short >= 1);
        assert_eq!((e.ani, e.af_a, e.af_b), (None, None, None));
        assert_eq!(e.ani_status, LocalStatus::NoSupportedMapping);
        assert_eq!(e.af_status, LocalStatus::NoSupportedMapping);
    }

    #[test]
    fn stop_gap_follows_the_run_hit_rate() {
        let genome = random_bases(20, 5_000);
        let (a, b) = sketches(&[&genome], &[&genome], &standard());
        let walker = Walker::new(&a, &b, &config());
        // 1 / (1/125 + 1/125) = 62.5 bases per sampled k-mer.
        assert_eq!(walker.spacing, 62.5);
        let gap = |walked, supports, walked_at_support| {
            walker.gap(&Progress {
                walked,
                supports,
                walked_at_support,
            })
        };
        assert_eq!(gap(1, 1, 1), 20_000.0);
        assert_eq!(gap(50, 1, 1), 20_000.0);
        assert_eq!(gap(20, 10, 20), 10.0 * 62.5 * 2.0);
        assert_eq!(gap(100, 100, 100), 625.0);
        assert_eq!(gap(10_000, 10, 10_000), 20_000.0);
        // Misses after the latest support leave the gap unchanged.
        assert_eq!(gap(5_000, 100, 100), 625.0);
        assert_eq!(gap(50, 10, 20), gap(20, 10, 20));
    }

    #[test]
    fn padding_extends_af_but_never_past_contig_ends() {
        let genome = random_bases(21, 20_000);
        let (a, b) = sketches(&[&genome], &[&genome], &standard());
        let bare = estimate(
            &a,
            &b,
            &LocalConfig {
                run_pad: Some(0),
                ..config()
            },
        );
        let (kept, _) = runs(&a, &b, &config());
        let (start, end) = kept[0].span_a();
        assert_eq!(bare.runs.unwrap().covered_bases_a, end - start);
        assert!(bare.af_a.unwrap() < 1.0);
        let padded = estimate(
            &a,
            &b,
            &LocalConfig {
                run_pad: Some(1_000),
                ..config()
            },
        );
        assert_eq!(padded.runs.unwrap().covered_bases_a, 20_000);
        assert_eq!(padded.af_a, Some(1.0));
        // ANI ignores padding.
        assert_eq!(bare.q, padded.q);
    }

    #[test]
    fn probes_match_only_near_the_place_the_run_maps_them_to() {
        // B's copy of A[40000, 40300) is replaced by a copy of A[70000, 70300):
        // those B probes occur in A, but 30 kb away from where the run maps them.
        let genome = random_bases(28, 100_000);
        let mut b = genome.clone();
        b.copy_within(70_000..70_300, 40_000);
        let (sa, sb) = sketches(&[&genome], &[&b], &dense());
        // A loose gap rule, so that one run spans the replaced stretch.
        let loose = LocalConfig {
            run_gap_factor: 100.0,
            ..config()
        };
        let (kept, _) = runs(&sa, &sb, &loose);
        assert!(kept
            .iter()
            .any(|r| r.span_a().0 < 40_000 && r.span_a().1 > 40_300));
        let moved: Vec<usize> = (0..sb.probes().len())
            .filter(|&j| (40_000..40_300 - 15).contains(&sb.probes().positions[j]))
            .collect();
        assert!(moved.len() > 10);
        assert!(moved.iter().all(|&j| sa
            .global_probes()
            .binary_search(&sb.probes().hash(j))
            .is_ok()));
        let replaced = (0..sa.probes().len())
            .filter(|&i| (40_000..40_300 - 15).contains(&sa.probes().positions[i]))
            .count() as u64;
        let r = estimate(&sa, &sb, &loose).runs.unwrap();
        // None of the moved B probes match; nor do A's probes they replaced.
        assert!(r.matched_b + moved.len() as u64 <= r.probes_b, "{r:?}");
        assert!(r.matched_a + replaced <= r.probes_a, "{r:?}");
    }

    #[test]
    fn runs_do_not_cross_long_stretches_at_another_offset() {
        // 100-bp insertions into B, each undone by a 100-bp deletion 6 kb later,
        // at 95% identity: the offset leaves and returns to the original one.
        // With a 50-bp band, runs break at each indel and must not walk on and
        // pick up the original offset 6 kb later, which would score the
        // stretch in between at the wrong offset.
        let genome = random_bases(29, 120_000);
        let mutated = substitute(&genome, 6_000, 6);
        let (mut b, mut aligned, mut same) = (Vec::new(), 0usize, 0usize);
        let mut i = 0;
        let mut event = 0;
        while i < genome.len() {
            if i % 6_000 == 3_000 && i + 3_000 < genome.len() {
                if event % 2 == 0 {
                    b.extend(random_bases(30 + event as u64, 100));
                } else {
                    i += 100;
                }
                event += 1;
            }
            b.push(mutated[i]);
            aligned += 1;
            same += usize::from(mutated[i] == genome[i]);
            i += 1;
        }
        let truth = same as f64 / aligned as f64;
        let narrow = LocalConfig {
            run_gap_factor: 20.0,
            run_max_gap: 50_000,
            run_band: 50,
            run_min_supports: 2,
            ..config()
        };
        let (sa, sb) = sketches(&[&genome], &[&b], &standard());
        let (kept, _) = runs(&sa, &sb, &narrow);
        let longest = kept
            .iter()
            .flat_map(|r| r.supports.windows(2).map(|w| w[1].x - w[0].x))
            .max()
            .unwrap();
        assert!(
            longest < 4_000,
            "a run crossed {longest} bases without support"
        );
        let ani = estimate(&sa, &sb, &narrow).ani.unwrap();
        assert!((ani - truth).abs() < 0.004, "{ani} vs {truth}");
    }

    #[test]
    fn repeats_do_not_bias_identity() {
        // 12 copies of a 3 kb element in 200 kb, every copy mutated
        // independently, 6 simulations: counting a repeated hash as matched
        // through any copy read about 1.7 pp high here. One 200 kb genome at
        // 90% identity has a sampling error of about 0.3 pp, so the mean of
        // 6 is tested.
        let mut errors = Vec::new();
        for seed in 0..6u64 {
            let element = random_bases(40 + seed, 3_000);
            let mut genome = random_bases(50 + seed, 200_000);
            for copy in 0..12 {
                let at = 5_000 + copy * 16_000;
                genome[at..at + element.len()].copy_from_slice(&element);
            }
            let mutated = substitute(&genome, 20_000, 60 + seed);
            let e = compare(&[&genome], &[&mutated], &standard(), &config());
            errors.push(e.ani.unwrap() - 0.9);
        }
        let mean = errors.iter().sum::<f64>() / errors.len() as f64;
        assert!(mean.abs() < 0.0025, "{errors:?}");
    }

    #[test]
    fn unrelated_genomes_have_no_runs() {
        let e = compare(
            &[&random_bases(22, 50_000)],
            &[&random_bases(23, 50_000)],
            &standard(),
            &config(),
        );
        assert_eq!((e.ani, e.q, e.af_a, e.af_b), (None, None, None, None));
        assert_eq!(e.ani_status, LocalStatus::NoSupportedMapping);
        assert_eq!(e.runs.unwrap().accepted_runs, 0);
    }

    #[test]
    fn thin_evidence_gives_null_ani_but_keeps_af() {
        let genome = random_bases(24, 30_000);
        let e = compare(
            &[&genome],
            &[&genome],
            &standard(),
            &LocalConfig {
                min_scoring_probes: 1_000_000,
                ..config()
            },
        );
        assert_eq!((e.ani, e.q), (None, None));
        assert_eq!(e.ani_status, LocalStatus::InsufficientEvidence);
        assert_eq!(e.af_status, LocalStatus::Estimated);
        assert!(e.af_a.unwrap() > 0.995);
    }

    #[test]
    fn estimates_are_deterministic() {
        let genome = random_bases(25, 60_000);
        let mutated = substitute(&genome, 3_000, 4);
        let first = compare(&[&genome], &[&mutated], &standard(), &config());
        let second = compare(&[&genome], &[&mutated], &standard(), &config());
        assert_eq!(first, second);
    }

    #[test]
    fn taken_intervals_merge_and_answer_overlaps() {
        let mut t = Taken::default();
        t.insert(10, 20);
        t.insert(30, 40);
        assert!(t.free(0, 10) && t.free(20, 30) && t.free(40, 50));
        assert!(!t.free(15, 16) && !t.free(5, 11) && !t.free(39, 45));
        t.insert(18, 32);
        assert_eq!(t.0.len(), 1);
        assert_eq!(t.0.get(&10), Some(&40));
        assert!(!t.free(25, 26) && t.free(40, 41));
    }

    #[test]
    fn run_settings_are_checked() {
        assert_eq!(config().check(), Ok(()));
        for (bad, message) in [
            (
                LocalConfig {
                    run_gap_factor: 0.0,
                    ..config()
                },
                "run_gap_factor",
            ),
            (
                LocalConfig {
                    run_gap_factor: f64::NAN,
                    ..config()
                },
                "run_gap_factor",
            ),
            (
                LocalConfig {
                    run_max_gap: 0,
                    ..config()
                },
                "run_max_gap",
            ),
            (
                LocalConfig {
                    run_min_supports: 1,
                    ..config()
                },
                "run_min_supports",
            ),
        ] {
            assert!(bad.check().unwrap_err().contains(message), "{bad:?}");
        }
    }
}
