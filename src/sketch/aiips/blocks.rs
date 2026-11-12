//! Positional correspondence from shared anchors, without alignment or DP.
//!
//! 1. A *shared anchor* is a hash that occurs exactly once in each genome;
//!    repeated hashes are ignored. Both anchor dictionaries are sorted, so
//!    shared anchors are found by one merge, with no hash table.
//! 2. Relative orientation is strand_A XOR strand_B. For reverse matches the B
//!    position is converted to the oriented coordinate y = len_B − (pos_B + k1),
//!    so collinear matches increase in both x (A position) and y.
//! 3. Matches are grouped by (contig A, contig B, orientation) and sorted by x;
//!    blocks never cross contigs.
//! 4. Consecutive matches extend a run when y increases, both gaps are at most
//!    `max_anchor_gap`, and max(gap) / min(gap) <= `max_gap_ratio`; otherwise the
//!    run breaks, and the reason is counted.
//! 5. A run is cut into blocks at the first match where the A span (last x + k1
//!    − first x) reaches `min_block_span`. Neighbouring blocks of one run share
//!    their boundary anchor. A shorter leftover joins the run's last block; a
//!    run that never reaches the span is still one candidate, left to
//!    validation.
//! 6. Blocks of different runs that overlap on either genome are resolved by
//!    keeping blocks with more shared anchors first (ties: earlier in A).

use std::cmp::Reverse;
use std::collections::BTreeMap;
use std::ops::Range;

use super::local::LocalConfig;
use super::occurrences::Occurrences;
use super::AiipsSketch;

/// One shared anchor: the single occurrence of a hash in each genome.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AnchorMatch {
    /// Anchor occurrence indices in A and in B.
    pub anchor_a: usize,
    pub anchor_b: usize,
    pub contig_a: usize,
    pub contig_b: usize,
    /// strand_A XOR strand_B: B is reverse-complemented relative to A.
    pub reverse: bool,
    /// Position in A.
    pub x: u32,
    /// Position in B (forward strand).
    pub pos_b: u32,
    /// Oriented B coordinate: `pos_b`, or len_B − (pos_b + k1) when reversed.
    pub y: u32,
}

/// Consecutive collinear matches of one run.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Block {
    /// Indices into [`Candidates::matches`].
    pub matches: Range<usize>,
    pub run: usize,
}

/// Counters for the matching stage.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MatchCounts {
    /// Hashes occurring exactly once in each genome.
    pub shared_anchors: u64,
    /// Hashes in both genomes but repeated in at least one; ignored.
    pub repeated_anchors: u64,
    /// Runs broken because y did not increase.
    pub breaks_order: u64,
    /// Runs broken because a gap exceeded `max_anchor_gap`.
    pub breaks_gap: u64,
    /// Runs broken because the gaps differed by more than `max_gap_ratio`.
    pub breaks_ratio: u64,
    /// Blocks cut from runs, before overlap resolution.
    pub candidate_blocks: u64,
    /// Candidate blocks dropped because they overlapped better-supported ones.
    pub overlapping_blocks: u64,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Candidates {
    /// Shared anchors sorted by (contig A, contig B, orientation, x).
    pub matches: Vec<AnchorMatch>,
    /// Non-overlapping candidate blocks, in match order.
    pub blocks: Vec<Block>,
    pub counts: MatchCounts,
}

impl Candidates {
    pub fn block_matches(&self, block: &Block) -> &[AnchorMatch] {
        &self.matches[block.matches.clone()]
    }

    /// Half-open interval of `block` on its A contig: [first x, last x + k1).
    pub fn interval_a(&self, block: &Block, k1: u32) -> (usize, u64, u64) {
        interval_a(self.block_matches(block), k1)
    }

    /// Half-open interval of `block` on its B contig, in forward coordinates.
    pub fn interval_b(&self, block: &Block, k1: u32) -> (usize, u64, u64) {
        interval_b(self.block_matches(block), k1)
    }
}

fn interval_a(m: &[AnchorMatch], k1: u32) -> (usize, u64, u64) {
    let (first, last) = (m[0], m[m.len() - 1]);
    (
        first.contig_a,
        first.x.into(),
        u64::from(last.x) + u64::from(k1),
    )
}

fn interval_b(m: &[AnchorMatch], k1: u32) -> (usize, u64, u64) {
    let (first, last) = (m[0], m[m.len() - 1]);
    let (low, high) = if first.reverse {
        (last.pos_b, first.pos_b)
    } else {
        (first.pos_b, last.pos_b)
    };
    (first.contig_b, low.into(), u64::from(high) + u64::from(k1))
}

/// Per dictionary entry: its occurrence index when it occurs exactly once.
fn single_occurrences(occurrences: &Occurrences) -> Vec<Option<usize>> {
    let mut seen = vec![0u8; occurrences.hashes.len()];
    let mut single = vec![None; occurrences.hashes.len()];
    for (index, &id) in occurrences.ids.iter().enumerate() {
        let id = id as usize;
        seen[id] = seen[id].saturating_add(1);
        single[id] = (seen[id] == 1).then_some(index);
    }
    single
}

/// Shared anchors of `a` and `b`, unsorted, and the count of repeated ones.
fn shared_anchors(a: &AiipsSketch, b: &AiipsSketch) -> (Vec<AnchorMatch>, u64) {
    let k1 = a.params.anchor_k.get() as u32;
    let (anchors_a, anchors_b) = (a.anchors(), b.anchors());
    let (single_a, single_b) = (single_occurrences(anchors_a), single_occurrences(anchors_b));
    let (hashes_a, hashes_b) = (&anchors_a.hashes, &anchors_b.hashes);
    let (mut i, mut j) = (0, 0);
    let mut matches = Vec::new();
    let mut repeated = 0;
    while i < hashes_a.len() && j < hashes_b.len() {
        match hashes_a[i].cmp(&hashes_b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                match (single_a[i], single_b[j]) {
                    (Some(anchor_a), Some(anchor_b)) => {
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
                        matches.push(AnchorMatch {
                            anchor_a,
                            anchor_b,
                            contig_a: a.contigs().of_anchor(anchor_a),
                            contig_b,
                            reverse,
                            x: anchors_a.positions[anchor_a],
                            pos_b,
                            y,
                        });
                    }
                    _ => repeated += 1,
                }
                i += 1;
                j += 1;
            }
        }
    }
    (matches, repeated)
}

/// Shared anchors, collinear runs and non-overlapping candidate blocks.
/// The sketches must have identical parameters.
pub fn candidate_blocks(a: &AiipsSketch, b: &AiipsSketch, config: &LocalConfig) -> Candidates {
    let k1 = a.params.anchor_k.get() as u32;
    let (mut matches, repeated) = shared_anchors(a, b);
    matches.sort_unstable_by_key(|m| (m.contig_a, m.contig_b, m.reverse, m.x));
    let mut counts = MatchCounts {
        shared_anchors: matches.len() as u64,
        repeated_anchors: repeated,
        ..MatchCounts::default()
    };
    let mut blocks = Vec::new();
    let mut run_start = 0;
    let mut run = 0;
    for next in 1..=matches.len() {
        let extends = next < matches.len() && {
            let (left, right) = (matches[next - 1], matches[next]);
            let same_group = (left.contig_a, left.contig_b, left.reverse)
                == (right.contig_a, right.contig_b, right.reverse);
            same_group && collinear(left, right, config, &mut counts)
        };
        if !extends {
            cut_run(&matches, run_start..next, run, k1, config, &mut blocks);
            run_start = next;
            run += 1;
        }
    }
    counts.candidate_blocks = blocks.len() as u64;
    let blocks = resolve_overlaps(&matches, blocks, k1);
    counts.overlapping_blocks = counts.candidate_blocks - blocks.len() as u64;
    Candidates {
        matches,
        blocks,
        counts,
    }
}

/// Whether `right` extends the run ending at `left`; counts the break reason.
fn collinear(
    left: AnchorMatch,
    right: AnchorMatch,
    config: &LocalConfig,
    counts: &mut MatchCounts,
) -> bool {
    // x strictly increases: distinct anchors of one contig differ in position.
    let gap_a = u64::from(right.x - left.x);
    if right.y <= left.y {
        counts.breaks_order += 1;
        return false;
    }
    let gap_b = u64::from(right.y - left.y);
    let (low, high) = (gap_a.min(gap_b), gap_a.max(gap_b));
    if high > u64::from(config.max_anchor_gap) {
        counts.breaks_gap += 1;
        return false;
    }
    if high as f64 > config.max_gap_ratio * low as f64 {
        counts.breaks_ratio += 1;
        return false;
    }
    true
}

/// Cuts one run of matches into blocks (see the module documentation).
fn cut_run(
    matches: &[AnchorMatch],
    range: Range<usize>,
    run: usize,
    k1: u32,
    config: &LocalConfig,
    blocks: &mut Vec<Block>,
) {
    if range.len() < 2 {
        return; // A lone shared anchor spans no interval.
    }
    let first_block = blocks.len();
    let mut start = range.start;
    for end in range.start + 1..range.end {
        let span = u64::from(matches[end].x) + u64::from(k1) - u64::from(matches[start].x);
        if span >= u64::from(config.min_block_span) {
            blocks.push(Block {
                matches: start..end + 1,
                run,
            });
            start = end;
        }
    }
    if start + 1 < range.end {
        if blocks.len() > first_block {
            blocks.last_mut().unwrap().matches.end = range.end;
        } else {
            blocks.push(Block {
                matches: start..range.end,
                run,
            });
        }
    }
}

/// Intervals kept so far on one contig: (start, run) -> end. Same-run blocks
/// may share a boundary anchor; different-run blocks never overlap.
#[derive(Default)]
struct Kept {
    intervals: BTreeMap<(u64, usize), u64>,
    longest: u64,
}

impl Kept {
    fn overlaps(&self, start: u64, end: u64, run: usize) -> bool {
        self.intervals
            .range(..(end, 0))
            .rev()
            .take_while(|((other_start, _), _)| other_start + self.longest > start)
            .any(|(&(_, other_run), &other_end)| other_end > start && other_run != run)
    }

    fn insert(&mut self, start: u64, end: u64, run: usize) {
        self.intervals.insert((start, run), end);
        self.longest = self.longest.max(end - start);
    }
}

/// Keeps blocks by descending number of shared anchors, skipping any that
/// overlaps a kept block of another run on either genome.
fn resolve_overlaps(matches: &[AnchorMatch], blocks: Vec<Block>, k1: u32) -> Vec<Block> {
    let mut order: Vec<usize> = (0..blocks.len()).collect();
    order.sort_by_key(|&i| (Reverse(blocks[i].matches.len()), i));
    let mut kept_a: BTreeMap<usize, Kept> = BTreeMap::new();
    let mut kept_b: BTreeMap<usize, Kept> = BTreeMap::new();
    let mut keep = vec![false; blocks.len()];
    for i in order {
        let block = &blocks[i];
        let (contig_a, start_a, end_a) = interval_a(&matches[block.matches.clone()], k1);
        let (contig_b, start_b, end_b) = interval_b(&matches[block.matches.clone()], k1);
        let clash = kept_a
            .get(&contig_a)
            .is_some_and(|kept| kept.overlaps(start_a, end_a, block.run))
            || kept_b
                .get(&contig_b)
                .is_some_and(|kept| kept.overlaps(start_b, end_b, block.run));
        if !clash {
            kept_a
                .entry(contig_a)
                .or_default()
                .insert(start_a, end_a, block.run);
            kept_b
                .entry(contig_b)
                .or_default()
                .insert(start_b, end_b, block.run);
            keep[i] = true;
        }
    }
    blocks
        .into_iter()
        .zip(keep)
        .filter_map(|(block, keep)| keep.then_some(block))
        .collect()
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::hash::fmix64;
    use crate::kmer::CanonicalKmers;
    use crate::sketch::aiips::tests::{build, fasta, params, random_bases, reverse_complement};
    use crate::sketch::aiips::AiipsParams;
    use std::num::NonZeroUsize;

    /// Dense anchors (about one per 20 bases) and probes for small genomes.
    pub(crate) fn dense() -> AiipsParams {
        params(21, 20, 15, 10, 42)
    }

    pub(crate) fn config(min_block_span: u32, max_anchor_gap: u32) -> LocalConfig {
        LocalConfig {
            min_block_span,
            max_anchor_gap,
            ..LocalConfig::default()
        }
    }

    fn run(
        a: &[&[u8]],
        b: &[&[u8]],
        p: &AiipsParams,
        c: &LocalConfig,
    ) -> (AiipsSketch, AiipsSketch, Candidates) {
        let (a, b) = (build(&fasta(a), p.clone()), build(&fasta(b), p.clone()));
        let candidates = candidate_blocks(&a, &b, c);
        (a, b, candidates)
    }

    fn span(matches: &[AnchorMatch]) -> u64 {
        u64::from(matches[matches.len() - 1].x) + 21 - u64::from(matches[0].x)
    }

    /// Whether `p` samples the k1-mer `kmer` as an anchor: NOT(h2(first k2-mer)
    /// XOR h2(last k2-mer)) <= T1.
    pub(crate) fn is_anchor(p: &AiipsParams, kmer: &[u8]) -> bool {
        let (k1, k2) = (p.anchor_k.get(), p.probe_k.get());
        let code = |window: &[u8]| {
            let k = NonZeroUsize::new(window.len()).unwrap();
            CanonicalKmers::new(window, k).unwrap().next().unwrap()
        };
        let first = fmix64(code(&kmer[..k2]) ^ p.probe_seed());
        let last = fmix64(code(&kmer[k1 - k2..]) ^ p.probe_seed());
        !(first ^ last) <= p.max_anchor_hash()
    }

    /// A random k-mer that `p` samples as an anchor.
    pub(crate) fn anchor_kmer(p: &AiipsParams, mut seed: u64) -> Vec<u8> {
        loop {
            let kmer = random_bases(seed, p.anchor_k.get());
            if is_anchor(p, &kmer) {
                return kmer;
            }
            seed += 1;
        }
    }

    #[test]
    fn identical_genomes_match_every_unique_anchor_in_one_run() {
        let genome = random_bases(1, 20_000);
        let c = config(2_000, 5_000);
        let (a, _, found) = run(&[&genome], &[&genome], &dense(), &c);
        let anchors = a.anchors();
        assert_eq!(found.counts.shared_anchors as usize, anchors.unique());
        assert_eq!(
            found.counts.repeated_anchors as usize,
            anchors.hashes.len() - anchors.unique()
        );
        assert!(found.matches.len() > 800);
        for m in &found.matches {
            assert!(!m.reverse);
            assert_eq!((m.x, m.y), (m.pos_b, m.x));
        }
        let counts = found.counts;
        assert_eq!(
            (counts.breaks_order, counts.breaks_gap, counts.breaks_ratio),
            (0, 0, 0)
        );
        assert_eq!(counts.overlapping_blocks, 0);
        // One run, cut into blocks that tile the matches and share boundaries.
        let blocks = &found.blocks;
        assert!(blocks.len() >= 8, "{} blocks", blocks.len());
        assert_eq!(blocks[0].matches.start, 0);
        assert_eq!(blocks[blocks.len() - 1].matches.end, found.matches.len());
        for pair in blocks.windows(2) {
            assert_eq!(pair[0].matches.end - 1, pair[1].matches.start);
        }
        // Each cut is at the first match reaching the span; the last block
        // also absorbs the shorter leftover.
        for block in &blocks[..blocks.len() - 1] {
            let m = found.block_matches(block);
            assert!(span(m) >= 2_000 && span(&m[..m.len() - 1]) < 2_000);
        }
        assert!(span(found.block_matches(&blocks[blocks.len() - 1])) >= 2_000);
    }

    #[test]
    fn reverse_complement_matches_have_identical_oriented_coordinates() {
        let genome = random_bases(2, 20_000);
        let c = config(2_000, 5_000);
        let (_, _, forward) = run(&[&genome], &[&genome], &dense(), &c);
        let (_, _, reverse) = run(&[&genome], &[&reverse_complement(&genome)], &dense(), &c);
        assert_eq!(forward.matches.len(), reverse.matches.len());
        for (f, r) in forward.matches.iter().zip(&reverse.matches) {
            assert!(r.reverse);
            assert_eq!((r.x, r.y), (f.x, f.y));
            assert_eq!(r.pos_b, 20_000 - 21 - f.x);
        }
        assert_eq!(forward.blocks, reverse.blocks);
        // B intervals mirror the A intervals.
        for block in &reverse.blocks {
            let (_, start_a, end_a) = reverse.interval_a(block, 21);
            let (_, start_b, end_b) = reverse.interval_b(block, 21);
            assert_eq!((start_b, end_b), (20_000 - end_a, 20_000 - start_a));
        }
    }

    #[test]
    fn repeated_anchor_hashes_are_ignored() {
        let repeat = random_bases(3, 3_000);
        let unique = random_bases(4, 3_000);
        let genome = [repeat.as_slice(), &unique, &repeat].concat();
        let (_, _, found) = run(&[&genome], &[&genome], &dense(), &config(1_000, 5_000));
        assert!(found.counts.repeated_anchors > 100);
        assert!(!found.matches.is_empty());
        // Only anchors overlapping the unique middle segment are matched
        // (those straddling a junction are unique too).
        assert!(found
            .matches
            .iter()
            .all(|m| (3_000 - 20..6_000).contains(&m.x)));
    }

    #[test]
    fn a_lost_intermediate_anchor_keeps_the_block_and_its_probes() {
        let genome = random_bases(5, 20_000);
        let c = config(2_000, 5_000);
        let (_, _, before) = run(&[&genome], &[&genome], &dense(), &c);
        let block = &before.blocks[2];
        let lost = before.matches[block.matches.start + block.matches.len() / 2];
        let mut mutated = genome.clone();
        let at = lost.x as usize + 10;
        mutated[at] = if mutated[at] == b'A' { b'C' } else { b'A' };
        let (a, _, after) = run(&[&genome], &[&mutated], &dense(), &c);
        assert!(
            after.matches.iter().all(|m| m.anchor_a != lost.anchor_a),
            "anchor lost"
        );
        let spanning = after
            .blocks
            .iter()
            .find(|block| {
                let (_, start, end) = after.interval_a(block, 21);
                start < u64::from(lost.x) && u64::from(lost.x) + 21 < end
            })
            .expect("a block spans the lost anchor");
        let m = after.block_matches(spanning);
        let slice = a
            .probes_between(m[0].anchor_a, m[m.len() - 1].anchor_a)
            .unwrap();
        // The slice holds the lost anchor's whole island, and exactly the
        // probes between the block's outer anchors.
        let island = a.owned_probes(lost.anchor_a);
        assert!(slice.start <= island.start && island.end <= slice.end);
        let positions = &a.probes().positions;
        let naive = positions
            .iter()
            .filter(|&&p| p >= m[0].x && p < m[m.len() - 1].x)
            .count();
        assert_eq!(slice.len(), naive);
    }

    #[test]
    fn an_inserted_anchor_between_conserved_anchors_keeps_the_block() {
        let genome = random_bases(6, 20_000);
        let p = dense();
        let c = config(2_000, 5_000);
        let (_, _, before) = run(&[&genome], &[&genome], &p, &c);
        // Insert into a gap of at least 40 bases inside a block, so the gap
        // ratio stays below 2.
        let block = &before.blocks[2];
        let m = before.block_matches(block);
        let gap = (1..m.len() - 1)
            .find(|&i| m[i + 1].x - m[i].x >= 40)
            .expect("a wide gap");
        let at = m[gap].x as usize + 21;
        let inserted = anchor_kmer(&p, 1_000);
        let b_genome = [&genome[..at], inserted.as_slice(), &genome[at..]].concat();
        let (_, b, after) = run(&[&genome], &[&b_genome], &p, &c);
        let new_anchor = (0..b.anchors().len())
            .find(|&i| b.anchors().positions[i] == at as u32)
            .expect("the inserted k-mer is a B anchor");
        assert!(
            after.matches.iter().all(|m| m.anchor_b != new_anchor),
            "not shared"
        );
        let spanning = after
            .blocks
            .iter()
            .find(|block| {
                let (_, start, end) = after.interval_b(block, 21);
                start < at as u64 && (at as u64) + 21 < end
            })
            .expect("a block spans the inserted anchor");
        let m = after.block_matches(spanning);
        let slice = b
            .probes_between(m[0].anchor_b, m[m.len() - 1].anchor_b)
            .unwrap();
        let island = b.owned_probes(new_anchor);
        assert!(slice.start <= island.start && island.end <= slice.end);
    }

    #[test]
    fn runs_break_for_order_gap_and_ratio_reasons() {
        let p = dense();
        let c = config(1_000, 5_000);
        let (left, right) = (random_bases(7, 6_000), random_bases(8, 6_000));
        let genome = [left.as_slice(), &right].concat();
        let count = |b: Vec<u8>| run(&[&genome], &[&b], &p, &c).2.counts;
        // Transposition: y falls once, at the junction.
        let swapped = count([right.as_slice(), &left].concat());
        assert_eq!(
            (
                swapped.breaks_order,
                swapped.breaks_gap,
                swapped.breaks_ratio
            ),
            (1, 0, 0)
        );
        // An 8 kb insertion exceeds the 5 kb gap limit.
        let long = count([left.as_slice(), &random_bases(9, 8_000), &right].concat());
        assert_eq!(
            (long.breaks_order, long.breaks_gap, long.breaks_ratio),
            (0, 1, 0)
        );
        // A 300-base insertion makes the gaps differ by far more than 2x.
        let short = count([left.as_slice(), &random_bases(10, 300), &right].concat());
        assert_eq!(
            (short.breaks_order, short.breaks_gap, short.breaks_ratio),
            (0, 0, 1)
        );
    }

    #[test]
    fn a_run_shorter_than_the_span_is_one_candidate() {
        let genome = random_bases(11, 1_500);
        let (_, _, found) = run(&[&genome], &[&genome], &dense(), &config(2_000, 5_000));
        assert_eq!(found.blocks.len(), 1);
        assert_eq!(found.blocks[0].matches, 0..found.matches.len());
        // A lone shared anchor makes no block.
        let k = NonZeroUsize::new(21).unwrap();
        let lone = anchor_kmer(&params(21, 100_000, 15, 10, 42), 5);
        let p = AiipsParams::new(k, 100_000, NonZeroUsize::new(15).unwrap(), 10, 42).unwrap();
        let a = [
            random_bases(12, 500).as_slice(),
            &lone,
            &random_bases(13, 500),
        ]
        .concat();
        let b = [
            random_bases(14, 500).as_slice(),
            &lone,
            &random_bases(15, 500),
        ]
        .concat();
        let (_, _, found) = run(&[&a], &[&b], &p, &config(10, 5_000));
        assert_eq!((found.counts.shared_anchors, found.blocks.len()), (1, 0));
    }

    #[test]
    fn overlapping_runs_keep_the_better_supported_block() {
        // A interleaves x1 y1 x2 y2. B contig 0 holds x1, filler, x2 and B
        // contig 1 holds y1, filler, y2, so each run's A interval spans the
        // other's segment.
        let (x1, x2) = (random_bases(16, 4_000), random_bases(17, 4_000));
        let (y1, y2) = (random_bases(18, 2_000), random_bases(19, 2_000));
        let a = [x1.as_slice(), &y1, &x2, &y2].concat();
        let b0 = [x1.as_slice(), &random_bases(20, 2_000), &x2].concat();
        let b1 = [y1.as_slice(), &random_bases(21, 4_000), &y2].concat();
        let (_, _, found) = run(&[&a], &[&b0, &b1], &dense(), &config(50_000, 5_000));
        assert_eq!(found.counts.candidate_blocks, 2);
        assert_eq!(found.counts.overlapping_blocks, 1);
        assert_eq!(found.blocks.len(), 1);
        let kept = found.block_matches(&found.blocks[0]);
        assert_eq!(kept[0].contig_b, 0, "the x run has more shared anchors");
    }
}
