//! Local comparison inside candidate blocks: validation with mapping probes,
//! ANI from scoring probes, and an alignment-fraction proxy. No base-level
//! alignment or dynamic programming is involved.
//!
//! **Validation** (each candidate block). On each genome, take the
//! mapping-role probes between the block's first and last shared anchor,
//! except those overlapping either of these two anchors: they match by
//! construction and carry no evidence (keeping them would also make the
//! result depend on orientation). Mask hashes that occur more than once in
//! either block. The remaining hashes
//! found in both give hits (x_A, y_B), with y_B oriented like the anchors:
//! len_B − (pos_B + k2) when the block is reversed. A block is accepted when
//! it has at least `min_mapping_hits` hits, and at least `min_order_fraction`
//! of adjacent hit pairs (in A order) increase in y. This is a provisional
//! heuristic, not a calibrated significance test.
//!
//! **ANI** (accepted blocks only). A scoring-role probe in the block is
//! eligible unless it overlaps the block's first or last shared anchor: the
//! block exists only because those two matched exactly, so overlapping probes
//! would be enriched for matches. Hashes occurring more than once in either
//! block are excluded from both sides.
//!
//! Excluding probes that overlap *any* of the block's shared anchors would
//! over-correct: the probes left behind then include those
//! overlapping anchors destroyed by a substitution (sampled but not shared),
//! which are depleted of matches. In simulations of 1 Mbp genomes with
//! substitutions only (3 replicates), it biased ANI by −0.0015, −0.0035 and
//! −0.0039 at 99%, 95% and 90% ANI; excluding only the two delimiting anchors
//! gave +0.0001, +0.0002 and −0.0009. Interior probes are counted whether or
//! not they overlap anchors, so enrichment and depletion cancel. With x shared eligible hashes and n_A, n_B
//! eligible hashes per side, counts are pooled over blocks:
//! q = 2 X / (N_A + N_B). Mapping-only local contrasts split heterogeneous
//! regions; ANI averages regional roots by span. Uniform comparisons retain
//! the pooled root, avoiding noise bias from small windows. Under the substitution
//! model a homologous k2-mer is unchanged with probability ANI^k2, so this is
//! a model-based estimate, not an alignment identity.
//!
//! ANI is null (never 0) in three cases: no block is accepted; fewer than
//! `min_scoring_probes` eligible probes remain; or no scoring probe is shared.
//!
//! **AF**. The union of the accepted blocks' intervals ([first anchor,
//! last anchor + k1)) on each genome, divided by that genome's total length,
//! ambiguous bases included. It is a positional coverage proxy: terminal
//! regions outside the outermost anchors are not counted, and a block holding
//! an unshared segment counts it. For fragmented inputs with no detected
//! heterogeneity and few repeats, a guarded C/q correction can increase AF.
//! The raw covered-base counters always retain observed span coverage.

use std::ops::Range;

use super::blocks::{candidate_blocks, Candidates, MatchCounts};
use super::runs::RunCounts;
use super::{AiipsSketch, Role};

pub const DEFAULT_MIN_BLOCK_SPAN: u32 = 10_000;
pub const DEFAULT_MAX_ANCHOR_GAP: u32 = 50_000;
pub const DEFAULT_MAX_GAP_RATIO: f64 = 2.0;
pub const DEFAULT_MIN_MAPPING_HITS: u32 = 3;
pub const DEFAULT_MIN_ORDER_FRACTION: f64 = 0.8;
pub const DEFAULT_MIN_SCORING_PROBES: u64 = 50;
pub const DEFAULT_RUN_GAP_FACTOR: f64 = 10.0;
pub const DEFAULT_RUN_MAX_GAP: u32 = 20_000;
pub const DEFAULT_RUN_BAND: u32 = 100;
pub const DEFAULT_RUN_MIN_SUPPORTS: u32 = 3;
/// One support (the seed alone) cannot establish homology.
pub const LEAST_RUN_SUPPORTS: u32 = 2;

/// How positional ANI and AF are computed.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum LocalMethod {
    /// Collinear shared anchors cut into blocks, validated by mapping probes.
    #[default]
    Blocks,
    /// Runs seeded by shared anchors and extended by anchors and probes
    /// (see [`super::runs`]).
    Runs,
}

impl LocalMethod {
    pub fn name(self) -> &'static str {
        match self {
            LocalMethod::Blocks => "blocks",
            LocalMethod::Runs => "runs",
        }
    }
}
/// One or two matching probes cannot establish homology.
pub const LEAST_MAPPING_HITS: u32 = 2;

/// Thresholds for positional comparison; prototype defaults, not tuned.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LocalConfig {
    /// Cut a block once its span on A reaches this many bases.
    pub min_block_span: u32,
    /// Break a run when either gap between consecutive shared anchors exceeds this.
    pub max_anchor_gap: u32,
    /// Break a run when max(gap_A, gap_B) / min(gap_A, gap_B) exceeds this.
    pub max_gap_ratio: f64,
    /// Unique mapping-probe hits needed to accept a block (at least 2).
    pub min_mapping_hits: u32,
    /// Share of adjacent hit pairs that must keep A's order in B.
    pub min_order_fraction: f64,
    /// Pooled eligible scoring probes (N_A + N_B) needed for ANI.
    pub min_scoring_probes: u64,
    pub method: LocalMethod,
    /// Runs: a gap may reach this many times the run's expected spacing
    /// between supports.
    pub run_gap_factor: f64,
    /// Runs: the largest gap without support, in bases.
    pub run_max_gap: u32,
    /// Runs: largest change of offset (B − A position) between supports.
    pub run_band: u32,
    /// Runs: supports needed to accept a run, seed included (at least 2).
    pub run_min_supports: u32,
    /// Runs: AF padding on each side of a run; `None` is the probe scale.
    pub run_pad: Option<u32>,
}

impl Default for LocalConfig {
    fn default() -> Self {
        Self {
            min_block_span: DEFAULT_MIN_BLOCK_SPAN,
            max_anchor_gap: DEFAULT_MAX_ANCHOR_GAP,
            max_gap_ratio: DEFAULT_MAX_GAP_RATIO,
            min_mapping_hits: DEFAULT_MIN_MAPPING_HITS,
            min_order_fraction: DEFAULT_MIN_ORDER_FRACTION,
            min_scoring_probes: DEFAULT_MIN_SCORING_PROBES,
            method: LocalMethod::Blocks,
            run_gap_factor: DEFAULT_RUN_GAP_FACTOR,
            run_max_gap: DEFAULT_RUN_MAX_GAP,
            run_band: DEFAULT_RUN_BAND,
            run_min_supports: DEFAULT_RUN_MIN_SUPPORTS,
            run_pad: None,
        }
    }
}

impl LocalConfig {
    /// Refuses values the stages cannot use meaningfully.
    pub fn check(&self) -> Result<(), String> {
        if self.min_block_span == 0 || self.max_anchor_gap == 0 {
            return Err("min_block_span and max_anchor_gap must be at least 1".into());
        }
        if !(self.max_gap_ratio.is_finite() && self.max_gap_ratio >= 1.0) {
            return Err(format!(
                "max_gap_ratio must be a finite number of at least 1, got {}",
                self.max_gap_ratio
            ));
        }
        if self.min_mapping_hits < LEAST_MAPPING_HITS {
            return Err(format!(
                "min_mapping_hits must be at least {LEAST_MAPPING_HITS}, got {}",
                self.min_mapping_hits
            ));
        }
        if !(0.0..=1.0).contains(&self.min_order_fraction) {
            return Err(format!(
                "min_order_fraction must be from 0 to 1, got {}",
                self.min_order_fraction
            ));
        }
        if self.min_scoring_probes == 0 {
            return Err("min_scoring_probes must be at least 1".into());
        }
        if !(self.run_gap_factor.is_finite() && self.run_gap_factor > 0.0) {
            return Err(format!(
                "run_gap_factor must be a finite number above 0, got {}",
                self.run_gap_factor
            ));
        }
        if self.run_max_gap == 0 {
            return Err("run_max_gap must be at least 1".into());
        }
        if self.run_min_supports < LEAST_RUN_SUPPORTS {
            return Err(format!(
                "run_min_supports must be at least {LEAST_RUN_SUPPORTS}, got {}",
                self.run_min_supports
            ));
        }
        Ok(())
    }
}

/// Why a positional estimate is or is not reported.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LocalStatus {
    Estimated,
    /// Blocks were accepted, but scoring evidence is too thin for ANI.
    InsufficientEvidence,
    /// No candidate block passed validation.
    NoSupportedMapping,
}

impl LocalStatus {
    pub fn name(self) -> &'static str {
        match self {
            LocalStatus::Estimated => "estimated",
            LocalStatus::InsufficientEvidence => "insufficient_evidence",
            LocalStatus::NoSupportedMapping => "no_supported_mapping",
        }
    }
}

/// Counters for every stage, for diagnostics and tuning.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct LocalCounts {
    pub matching: MatchCounts,
    /// Candidate blocks with fewer than `min_mapping_hits` hits.
    pub rejected_few_hits: u64,
    /// Candidate blocks whose hits were too often out of order.
    pub rejected_disorder: u64,
    pub accepted_blocks: u64,
    /// Mapping-probe hits in accepted blocks.
    pub mapping_hits: u64,
    /// Pooled eligible scoring hashes: N_A, N_B and the shared X.
    pub scoring_a: u64,
    pub scoring_b: u64,
    pub scoring_shared: u64,
    /// Scoring probes excluded for overlapping a block's first or last anchor.
    pub excluded_boundary_overlap: u64,
    /// Scoring probes excluded because their hash repeats within a block.
    pub excluded_repeats: u64,
    /// Bases covered by the union of accepted blocks on each genome.
    pub covered_bases_a: u64,
    pub covered_bases_b: u64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LocalEstimates {
    /// Regional mean identity (pooled root for uniform inputs), or `None`.
    pub ani: Option<f64>,
    /// 2 X / (N_A + N_B), k2-mer survival; differs from mean ANI^k2 in mixtures.
    pub q: Option<f64>,
    /// Estimated AF, including guarded fragmentation correction when applicable.
    pub af_a: Option<f64>,
    pub af_b: Option<f64>,
    pub ani_status: LocalStatus,
    pub af_status: LocalStatus,
    pub counts: LocalCounts,
    /// Run-method counters; `None` for blocks.
    pub runs: Option<RunCounts>,
}

/// Every stage's result: candidate blocks, which were accepted, and the estimates.
pub struct Analysis {
    pub candidates: Candidates,
    /// Per candidate block: accepted by validation.
    pub accepted: Vec<bool>,
    pub estimates: LocalEstimates,
}

/// Positional ANI and AF of two sketches with identical parameters.
pub fn estimate_local(a: &AiipsSketch, b: &AiipsSketch, config: &LocalConfig) -> LocalEstimates {
    analyze(a, b, config).estimates
}

/// One genome's view of a block: its first and last shared anchor, in
/// position order on that genome.
struct Side<'s> {
    sketch: &'s AiipsSketch,
    first: usize,
    last: usize,
}

/// (hash, position) entries of one block's probes.
type Entries = Vec<(u64, u32)>;

impl Side<'_> {
    /// Probes lying entirely within the block interval [first anchor, last
    /// anchor + k1): the islands from the first to the last anchor, plus the
    /// start of the last anchor's island. Symmetric under reverse complement.
    fn probes(&self, k1: u32, k2: u32) -> Range<usize> {
        let between = self
            .sketch
            .probes_between(self.first, self.last)
            .expect("a block's anchors lie on one contig, in order");
        let contig_end = self
            .sketch
            .contigs()
            .probes(self.sketch.contigs().of_anchor(self.first))
            .end;
        let end_of_block = u64::from(self.sketch.anchors().positions[self.last]) + u64::from(k1);
        let positions = &self.sketch.probes().positions;
        let inside = positions[between.start..contig_end]
            .partition_point(|&p| u64::from(p) + u64::from(k2) <= end_of_block);
        between.start..between.start + inside
    }

    /// Splits the block's probes that overlap neither boundary anchor into
    /// mapping and scoring entries (buffers are cleared first); returns how
    /// many scoring probes overlapped a boundary anchor.
    fn split_interior(
        &self,
        k1: u32,
        k2: u32,
        mapping: &mut Entries,
        scoring: &mut Entries,
    ) -> u64 {
        mapping.clear();
        scoring.clear();
        let probes = self.sketch.probes();
        let anchors = &self.sketch.anchors().positions;
        let (first, last) = (anchors[self.first], anchors[self.last]);
        let overlaps = |start: u32, anchor: u32| start < anchor + k1 && anchor < start + k2;
        let mut overlapping = 0;
        for p in self.probes(k1, k2) {
            let start = probes.positions[p];
            let is_scoring = self.sketch.probe_role(p) == Role::Scoring;
            if overlaps(start, first) || overlaps(start, last) {
                overlapping += u64::from(is_scoring);
            } else if is_scoring {
                scoring.push((probes.hash(p), start));
            } else {
                mapping.push((probes.hash(p), start));
            }
        }
        overlapping
    }
}

/// Sorts by hash and keeps only entries whose hash occurs exactly once.
fn keep_singles(entries: &mut Entries) {
    entries.sort_unstable();
    let (mut write, mut i) = (0, 0);
    while i < entries.len() {
        let run = entries[i..]
            .iter()
            .take_while(|e| e.0 == entries[i].0)
            .count();
        if run == 1 {
            entries[write] = entries[i];
            write += 1;
        }
        i += run;
    }
    entries.truncate(write);
}

/// Mapping hits between two blocks' single-copy mapping entries (sorted by
/// hash): their number and how many adjacent pairs (in A order) increase in
/// oriented B coordinate.
fn mapping_hits(
    a: &[(u64, u32)],
    b: &[(u64, u32)],
    reverse: bool,
    k2: u32,
    length_b: u32,
    hits: &mut Vec<(u32, u32)>,
) -> (u64, u64) {
    hits.clear();
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].0.cmp(&b[j].0) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                // Stored probes end inside their contig: no underflow.
                let y = if reverse {
                    length_b - k2 - b[j].1
                } else {
                    b[j].1
                };
                hits.push((a[i].1, y));
                i += 1;
                j += 1;
            }
        }
    }
    hits.sort_unstable();
    let ordered = hits.windows(2).filter(|w| w[1].1 > w[0].1).count();
    (hits.len() as u64, ordered as u64)
}

/// (shared, n_A, n_B) over two blocks' eligible scoring entries, excluding
/// hashes repeated on either side from both.
fn score(a: &mut Entries, b: &mut Entries, repeats: &mut u64) -> (u64, u64, u64) {
    a.sort_unstable();
    b.sort_unstable();
    let (mut i, mut j) = (0, 0);
    let (mut shared, mut n_a, mut n_b) = (0, 0, 0);
    while i < a.len() || j < b.len() {
        let hash = match (a.get(i), b.get(j)) {
            (Some(x), Some(y)) => x.0.min(y.0),
            (Some(x), None) | (None, Some(x)) => x.0,
            (None, None) => unreachable!(),
        };
        let count_a = a[i..].iter().take_while(|e| e.0 == hash).count() as u64;
        let count_b = b[j..].iter().take_while(|e| e.0 == hash).count() as u64;
        i += count_a as usize;
        j += count_b as usize;
        if count_a > 1 || count_b > 1 {
            *repeats += count_a + count_b;
        } else {
            n_a += count_a;
            n_b += count_b;
            shared += count_a * count_b;
        }
    }
    (shared, n_a, n_b)
}

/// Total length of the union of half-open (contig, start, end) intervals.
pub(crate) fn covered(mut intervals: Vec<(usize, u64, u64)>) -> u64 {
    intervals.sort_unstable();
    let mut total = 0;
    let mut current: Option<(usize, u64, u64)> = None;
    for (contig, start, end) in intervals {
        match current {
            Some((c, s, e)) if c == contig && start <= e => current = Some((c, s, e.max(end))),
            _ => {
                if let Some((_, s, e)) = current {
                    total += e - s;
                }
                current = Some((contig, start, end));
            }
        }
    }
    total + current.map_or(0, |(_, s, e)| e - s)
}

/// q = 2 X / (N_A + N_B) from pooled scoring counts; `None` when fewer than
/// `min_scoring_probes` are eligible or none is shared (never q = 0).
fn survival(counts: &LocalCounts, config: &LocalConfig) -> Option<f64> {
    let pooled = counts.scoring_a + counts.scoring_b;
    (pooled >= config.min_scoring_probes && counts.scoring_shared > 0)
        .then(|| 2.0 * counts.scoring_shared as f64 / pooled as f64)
}

#[derive(Clone, Copy, Default)]
struct Region {
    mapping_shared: f64,
    mapping_union: f64,
    shared: u64,
    pooled: u64,
    span: u64,
}

impl Region {
    fn add(self, b: Self) -> Self {
        Self {
            mapping_shared: self.mapping_shared + b.mapping_shared,
            mapping_union: self.mapping_union + b.mapping_union,
            shared: self.shared + b.shared,
            pooled: self.pooled + b.pooled,
            span: self.span + b.span,
        }
    }
    fn likelihood(self) -> f64 {
        let x = self.mapping_shared;
        let n = self.mapping_union;
        if x == 0.0 || x == n {
            return 0.0;
        }
        x * (x / n).ln() + (n - x) * ((n - x) / n).ln()
    }
}

/// Binary segmentation using only mapping counts. Conditional on the union
/// count, independent sampled shared/unshared hashes give a binomial working
/// likelihood. Overlapping k-mers violate independence: the conservative
/// 2*ln(N) penalty is a heuristic, not a significance level. Scoring counts
/// never choose a split. Prefix sums bound work and recursion to O(B log B)
/// for balanced splits, with a hard depth limit for adversarial input.
/// A leaf needs at least `min_pooled` eligible scoring probes
/// (`--min-scoring-probes`), as the pooled estimate does.
fn regional_ani(regions: &[Region], k: u32, min_pooled: u64) -> (f64, bool) {
    fn split(
        prefix: &[Region],
        lo: usize,
        hi: usize,
        depth: usize,
        k: u32,
        min_pooled: u64,
    ) -> (f64, bool) {
        let difference = |a: Region, b: Region| Region {
            mapping_shared: a.mapping_shared - b.mapping_shared,
            mapping_union: a.mapping_union - b.mapping_union,
            shared: a.shared - b.shared,
            pooled: a.pooled - b.pooled,
            span: a.span - b.span,
        };
        let total = difference(prefix[hi], prefix[lo]);
        let mut best = (0.0, lo);
        if depth < 16 {
            for mid in lo + 1..hi {
                // Local contrast also detects alternating regions whose
                // global left/right averages are almost identical.
                let left = difference(prefix[mid], prefix[mid.saturating_sub(4).max(lo)]);
                let right = difference(prefix[(mid + 4).min(hi)], prefix[mid]);
                if left.mapping_union < 40.0 || right.mapping_union < 40.0 {
                    continue;
                }
                let gain = left.likelihood() + right.likelihood()
                    - left.add(right).likelihood()
                    - 2.0 * (left.mapping_union + right.mapping_union).max(2.0).ln();
                if gain > best.0 {
                    best = (gain, mid);
                }
            }
        }
        if best.1 != lo {
            let a = split(prefix, lo, best.1, depth + 1, k, min_pooled);
            let b = split(prefix, best.1, hi, depth + 1, k, min_pooled);
            (a.0 + b.0, true)
        } else {
            // A segment with no scoring evidence uses the pooled parent
            // estimate at the caller, rather than selecting positive windows.
            if total.shared == 0 || total.pooled < min_pooled {
                return (f64::NAN, false);
            }
            let q = 2.0 * total.shared as f64 / total.pooled as f64;
            (total.span as f64 * q.powf(1.0 / f64::from(k)), false)
        }
    }
    let mut prefix = vec![Region::default()];
    for &r in regions {
        prefix.push(prefix.last().unwrap().add(r));
    }
    let (sum, heterogeneous) = split(&prefix, 0, regions.len(), 0, k, min_pooled);
    (
        sum / prefix.last().unwrap().span.max(1) as f64,
        heterogeneous,
    )
}

/// Whether fewer than 2% of a genome's distinct anchor hashes repeat, relative
/// to those occurring once. This is the shared-anchor ratio
/// `repeated_anchors / shared_anchors` at 100% identity, measured on one
/// genome: shared counts depend on divergence (a repeated hash is kept if any
/// copy survives, a unique one only if its single copy does), this does not.
fn few_repeats(sketch: &AiipsSketch, unique: usize) -> bool {
    let anchors = sketch.anchors();
    (anchors.hashes.len() - unique) * 50 < unique
}

pub fn analyze(a: &AiipsSketch, b: &AiipsSketch, config: &LocalConfig) -> Analysis {
    let global = super::estimate(a, b).ok();
    analyze_impl(a, b, config, global.as_ref())
}

pub(crate) fn analyze_with_global(
    a: &AiipsSketch,
    b: &AiipsSketch,
    config: &LocalConfig,
    global: &super::AiipsEstimates,
) -> Analysis {
    analyze_impl(a, b, config, Some(global))
}

fn analyze_impl(
    a: &AiipsSketch,
    b: &AiipsSketch,
    config: &LocalConfig,
    global: Option<&super::AiipsEstimates>,
) -> Analysis {
    let k1 = a.params.anchor_k.get() as u32;
    let k2 = a.params.probe_k.get() as u32;
    let candidates = candidate_blocks(a, b, config);
    let mut counts = LocalCounts {
        matching: candidates.counts,
        ..LocalCounts::default()
    };
    let mut accepted = vec![false; candidates.blocks.len()];
    let (mut covered_a, mut covered_b) = (Vec::new(), Vec::new());
    // Buffers reused across blocks: [A, B] mapping and scoring entries, hits.
    let (mut mapping, mut scoring): ([Entries; 2], [Entries; 2]) = Default::default();
    let mut hits = Vec::new();
    let mut regions = Vec::new();
    for (index, block) in candidates.blocks.iter().enumerate() {
        let matches = candidates.block_matches(block);
        let (first, last) = (matches[0], matches[matches.len() - 1]);
        let reverse = first.reverse;
        let side_a = Side {
            sketch: a,
            first: first.anchor_a,
            last: last.anchor_a,
        };
        // Reversed blocks run backwards along B.
        let (low, high) = if reverse {
            (last.anchor_b, first.anchor_b)
        } else {
            (first.anchor_b, last.anchor_b)
        };
        let side_b = Side {
            sketch: b,
            first: low,
            last: high,
        };
        let [mapping_a, mapping_b] = &mut mapping;
        let [scoring_a, scoring_b] = &mut scoring;
        let overlapping = side_a.split_interior(k1, k2, mapping_a, scoring_a)
            + side_b.split_interior(k1, k2, mapping_b, scoring_b);
        keep_singles(mapping_a);
        keep_singles(mapping_b);
        let length_b = b.contigs().length(first.contig_b);
        let (hit_count, ordered) =
            mapping_hits(mapping_a, mapping_b, reverse, k2, length_b, &mut hits);
        if hit_count < u64::from(config.min_mapping_hits) {
            counts.rejected_few_hits += 1;
            continue;
        }
        // min_mapping_hits >= 2, so there is at least one adjacent pair.
        if (ordered as f64) < config.min_order_fraction * (hit_count - 1) as f64 {
            counts.rejected_disorder += 1;
            continue;
        }
        accepted[index] = true;
        counts.accepted_blocks += 1;
        counts.mapping_hits += hit_count;
        counts.excluded_boundary_overlap += overlapping;
        let (shared, n_a, n_b) = score(scoring_a, scoring_b, &mut counts.excluded_repeats);
        counts.scoring_shared += shared;
        counts.scoring_a += n_a;
        counts.scoring_b += n_b;
        regions.push(Region {
            mapping_shared: hit_count as f64,
            mapping_union: (mapping_a.len() + mapping_b.len()) as f64 - hit_count as f64,
            shared,
            pooled: n_a + n_b,
            span: u64::from(last.x - first.x),
        });
        covered_a.push(candidates.interval_a(block, k1));
        covered_b.push(candidates.interval_b(block, k1));
    }
    counts.covered_bases_a = covered(covered_a);
    counts.covered_bases_b = covered(covered_b);

    let mut estimates = if counts.accepted_blocks == 0 {
        LocalEstimates {
            ani: None,
            q: None,
            af_a: None,
            af_b: None,
            ani_status: LocalStatus::NoSupportedMapping,
            af_status: LocalStatus::NoSupportedMapping,
            counts,
            runs: None,
        }
    } else {
        let q = survival(&counts, config);
        let (regional, heterogeneous) = regional_ani(&regions, k2, config.min_scoring_probes);
        let mut af_a = counts.covered_bases_a as f64 / a.summary.bases as f64;
        let mut af_b = counts.covered_bases_b as f64 / b.summary.bases as f64;
        // Correct fragmentation loss only when mapping finds no identity
        // heterogeneity and neither genome is repetitive. Islands and
        // rearrangements break runs but are what C/q measures, so run breaks
        // do not block the correction. Preserve raw span numerators as
        // diagnostics. C/q remains a k-mer model, not a base alignment.
        if (a.contigs().len() > 1 || b.contigs().len() > 1)
            && !heterogeneous
            && few_repeats(
                a,
                candidates.unique_by_contig[0]
                    .iter()
                    .map(|&n| n as usize)
                    .sum(),
            )
            && few_repeats(
                b,
                candidates.unique_by_contig[1]
                    .iter()
                    .map(|&n| n as usize)
                    .sum(),
            )
            && counts.excluded_repeats * 50 < counts.scoring_a + counts.scoring_b
        {
            if let (Some(q), Some(global)) = (q, global) {
                af_a = (global.containment_a_in_b / q).clamp(af_a, 1.0);
                af_b = (global.containment_b_in_a / q).clamp(af_b, 1.0);
            }
        }
        LocalEstimates {
            ani: q.map(|q| {
                if heterogeneous && regional.is_finite() {
                    regional
                } else {
                    q.powf(1.0 / f64::from(k2))
                }
            }),
            q,
            af_a: Some(af_a),
            af_b: Some(af_b),
            ani_status: if q.is_some() {
                LocalStatus::Estimated
            } else {
                LocalStatus::InsufficientEvidence
            },
            af_status: LocalStatus::Estimated,
            counts,
            runs: None,
        }
    };
    if let Some(global) = global {
        super::heterogeneous::adjust(a, b, &candidates, global, &mut estimates);
    }
    Analysis {
        candidates,
        accepted,
        estimates,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sketch::aiips::blocks::tests::{anchor_kmer, config, dense};
    use crate::sketch::aiips::tests::{build, fasta, params, random_bases, reverse_complement};
    use std::collections::HashMap;

    /// A block's counts at identity `p`, k2 = 15, with expected values.
    fn region(p: f64) -> Region {
        let q = p.powi(15);
        Region {
            mapping_shared: 1000.0 * q,
            mapping_union: 1000.0 * (2.0 - q),
            shared: (10000.0 * q).round() as u64,
            pooled: 20000,
            span: 10000,
        }
    }

    /// 40 blocks alternating between 98% and 85% identity every 5 blocks.
    fn alternating() -> Vec<Region> {
        (0..40)
            .map(|i| region(if i / 5 % 2 == 0 { 0.98 } else { 0.85 }))
            .collect()
    }

    #[test]
    fn regional_estimator_preserves_uniform_pool_and_separates_alternating_identity() {
        let uniform = vec![region(0.95); 40];
        let (ani, changed) = regional_ani(&uniform, 15, DEFAULT_MIN_SCORING_PROBES);
        assert!(!changed);
        assert!((ani - 0.95).abs() < 0.0001);
        let (ani, changed) = regional_ani(&alternating(), 15, DEFAULT_MIN_SCORING_PROBES);
        assert!(changed);
        assert!((ani - 0.915).abs() < 0.0001, "{ani}");
    }

    #[test]
    fn regional_leaves_need_the_configured_scoring_probes() {
        // One unsplit leaf pools 40 * 20000 eligible probes.
        let uniform = vec![region(0.95); 40];
        assert!(regional_ani(&uniform, 15, 800_000).0.is_finite());
        let (ani, changed) = regional_ani(&uniform, 15, 800_001);
        assert!(ani.is_nan() && !changed);
        // Split leaves of 5 blocks pool 100000: a minimum above that leaves
        // the split flagged without a regional estimate.
        assert!(regional_ani(&alternating(), 15, 100_000).0.is_finite());
        let (ani, changed) = regional_ani(&alternating(), 15, 100_001);
        assert!(ani.is_nan() && changed);
    }

    #[test]
    fn few_repeats_is_the_shared_anchor_repeat_ratio_at_full_identity() {
        let p = params(21, 20, 15, 10, 42);
        let mut genome = random_bases(8, 200_000);
        let element = random_bases(7, 2_000);
        for at in [10_000, 90_000, 150_000] {
            genome[at..at + element.len()].copy_from_slice(&element);
        }
        let s = build(&fasta(&[&genome]), p.clone());
        let anchors = s.anchors();
        let repeated = anchors.hashes.len() - anchors.unique();
        let m = candidate_blocks(&s, &s, &LocalConfig::default()).counts;
        assert_eq!(m.shared_anchors as usize, anchors.unique());
        assert_eq!(m.repeated_anchors as usize, repeated);
        // About 100 repeated hashes against 9700 unique: under 2%.
        assert!(
            repeated > 50 && few_repeats(&s, s.anchors().unique()),
            "{repeated}"
        );
        // An 8 kb element in 3 copies: about 400 repeated, over 2%.
        let element = random_bases(9, 8_000);
        for at in [10_000, 90_000, 150_000] {
            genome[at..at + element.len()].copy_from_slice(&element);
        }
        let s = build(&fasta(&[&genome]), p);
        assert!(!few_repeats(&s, s.anchors().unique()));
    }

    fn substitute(sequence: &[u8], every: usize) -> Vec<u8> {
        let mut out = sequence.to_vec();
        for (i, base) in out.iter_mut().enumerate() {
            if i % every == every / 2 {
                *base = if *base == b'A' { b'G' } else { b'A' };
            }
        }
        out
    }

    fn analyze_pair(
        a: &[&[u8]],
        b: &[&[u8]],
        p: &AiipsParams,
        c: &LocalConfig,
    ) -> (AiipsSketch, AiipsSketch, Analysis) {
        let (a, b) = (build(&fasta(a), p.clone()), build(&fasta(b), p.clone()));
        let analysis = analyze(&a, &b, c);
        (a, b, analysis)
    }

    use crate::sketch::aiips::AiipsParams;

    #[test]
    fn identical_and_reverse_complement_genomes_give_ani_one() {
        let genome = random_bases(1, 30_000);
        let c = config(2_000, 5_000);
        let (_, _, same) = analyze_pair(&[&genome], &[&genome], &dense(), &c);
        let (_, _, rc) = analyze_pair(&[&genome], &[&reverse_complement(&genome)], &dense(), &c);
        let e = same.estimates;
        assert_eq!((e.ani, e.q), (Some(1.0), Some(1.0)));
        assert_eq!(
            (e.ani_status, e.af_status),
            (LocalStatus::Estimated, LocalStatus::Estimated)
        );
        let n = e.counts;
        assert_eq!(
            (n.scoring_a, n.scoring_b),
            (n.scoring_shared, n.scoring_shared)
        );
        assert_eq!((n.rejected_few_hits, n.rejected_disorder), (0, 0));
        assert_eq!(e.af_a, e.af_b);
        assert!(e.af_a.unwrap() > 0.99, "{:?}", e.af_a);
        // Reverse complement: the same blocks, counts and estimates.
        assert_eq!(rc.estimates, e);
    }

    #[test]
    fn unrelated_genomes_have_no_supported_mapping() {
        let (a, b) = (random_bases(2, 30_000), random_bases(3, 30_000));
        let (_, _, analysis) = analyze_pair(&[&a], &[&b], &dense(), &config(2_000, 5_000));
        let e = analysis.estimates;
        assert_eq!((e.ani, e.q, e.af_a, e.af_b), (None, None, None, None));
        assert_eq!(e.ani_status, LocalStatus::NoSupportedMapping);
        assert_eq!(e.af_status, LocalStatus::NoSupportedMapping);
        assert_eq!(e.counts.accepted_blocks, 0);
        assert_eq!((e.counts.covered_bases_a, e.counts.covered_bases_b), (0, 0));
    }

    #[test]
    fn counts_match_a_naive_recomputation() {
        // Two contigs, substitutions every 37 bases (about 97% identity), and a
        // reverse-complemented contig in B, so every path is exercised.
        let (one, two) = (random_bases(4, 25_000), random_bases(5, 15_000));
        let b_two = reverse_complement(&substitute(&two, 37));
        let p = dense();
        let (a, b, analysis) = analyze_pair(
            &[&one, &two],
            &[&substitute(&one, 37), &b_two],
            &p,
            &config(2_000, 5_000),
        );
        let candidates = &analysis.candidates;
        assert!(analysis.accepted.iter().filter(|&&x| x).count() > 10);
        assert!(candidates.matches.iter().any(|m| m.reverse));
        let (k1, k2) = (21u32, 15u32);
        let (mut shared, mut n_a, mut n_b) = (0u64, 0u64, 0u64);
        let mut cover = [
            vec![false; 25_000],
            vec![false; 15_000],
            vec![false; 25_000],
            vec![false; 15_000],
        ];
        for (block, _) in candidates
            .blocks
            .iter()
            .zip(&analysis.accepted)
            .filter(|(_, &ok)| ok)
        {
            let m = candidates.block_matches(block);
            // Naive: scan every probe of the contig by position.
            let eligible =
                |sketch: &AiipsSketch, contig: usize, ends: [u32; 2]| -> HashMap<u64, u64> {
                    let (low, high) = (ends[0].min(ends[1]), ends[0].max(ends[1]));
                    let probes = sketch.probes();
                    let mut counts = HashMap::new();
                    for i in sketch.contigs().probes(contig) {
                        let start = probes.positions[i];
                        let touches = |anchor: u32| start < anchor + k1 && anchor < start + k2;
                        if start >= low
                            && start < high
                            && sketch.probe_role(i) == Role::Scoring
                            && !touches(low)
                            && !touches(high)
                        {
                            *counts.entry(probes.hash(i)).or_insert(0) += 1;
                        }
                    }
                    counts
                };
            let (first, last) = (m[0], m[m.len() - 1]);
            let side_a = eligible(&a, first.contig_a, [first.x, last.x]);
            let side_b = eligible(&b, first.contig_b, [first.pos_b, last.pos_b]);
            for (hash, &count) in &side_a {
                let other = side_b.get(hash).copied().unwrap_or(0);
                if count == 1 && other <= 1 {
                    n_a += 1;
                    shared += other;
                }
            }
            n_b += side_b
                .iter()
                .filter(|(hash, &count)| count == 1 && side_a.get(*hash).copied().unwrap_or(0) <= 1)
                .count() as u64;
            let (contig, start, end) = candidates.interval_a(block, k1);
            cover[contig][start as usize..end as usize].fill(true);
            let (contig, start, end) = candidates.interval_b(block, k1);
            cover[2 + contig][start as usize..end as usize].fill(true);
        }
        let n = analysis.estimates.counts;
        assert_eq!(
            (n.scoring_shared, n.scoring_a, n.scoring_b),
            (shared, n_a, n_b)
        );
        let covered = |v: &[Vec<bool>]| v.iter().flatten().filter(|&&c| c).count() as u64;
        assert_eq!(n.covered_bases_a, covered(&cover[..2]));
        assert_eq!(n.covered_bases_b, covered(&cover[2..]));
        let e = analysis.estimates;
        assert_eq!(e.af_a, Some(n.covered_bases_a as f64 / 40_000.0));
        let q = 2.0 * shared as f64 / (n_a + n_b) as f64;
        assert_eq!(e.ani, Some(q.powf(1.0 / 15.0)));
    }

    /// Two planted anchors around `a_middle` / `b_middle`, with anchors so
    /// sparse (1 in 100,000) that no other is sampled, and dense probes.
    fn planted(a_middle: &[u8], b_middle: &[u8]) -> Analysis {
        let p = params(21, 100_000, 15, 4, 42);
        let (left, right) = (anchor_kmer(&p, 100), anchor_kmer(&p, 900_000));
        let a = [left.as_slice(), a_middle, &right].concat();
        let b = [left.as_slice(), b_middle, &right].concat();
        let (sa, sb, analysis) = analyze_pair(&[&a], &[&b], &p, &config(2_000, 10_000));
        assert_eq!(
            (sa.anchors().len(), sb.anchors().len()),
            (2, 2),
            "only the planted anchors"
        );
        assert_eq!(analysis.candidates.blocks.len(), 1);
        analysis
    }

    #[test]
    fn blocks_without_mapping_support_are_rejected() {
        let analysis = planted(&random_bases(30, 6_000), &random_bases(31, 6_000));
        let e = analysis.estimates;
        assert_eq!(e.counts.rejected_few_hits, 1);
        assert_eq!(e.counts.accepted_blocks, 0);
        assert_eq!((e.ani, e.af_a), (None, None));
        assert_eq!(e.ani_status, LocalStatus::NoSupportedMapping);
    }

    #[test]
    fn blocks_with_disordered_hits_are_rejected() {
        // Same content, but 30-base chunks in reverse order: plenty of hits,
        // about half of adjacent pairs out of order.
        let middle = random_bases(32, 6_000);
        let shuffled: Vec<u8> = middle.chunks(30).rev().flatten().copied().collect();
        let e = planted(&middle, &shuffled).estimates;
        assert_eq!(e.counts.rejected_disorder, 1);
        assert_eq!(e.counts.accepted_blocks, 0);
        // In order, the same content is accepted.
        let e = planted(&middle, &middle).estimates;
        assert_eq!(e.counts.accepted_blocks, 1);
        assert!(e.counts.mapping_hits > 100);
    }

    #[test]
    fn thin_scoring_evidence_gives_null_ani_but_keeps_af() {
        let genome = random_bases(6, 30_000);
        let c = LocalConfig {
            min_scoring_probes: 1_000_000,
            ..config(2_000, 5_000)
        };
        let (_, _, analysis) = analyze_pair(&[&genome], &[&genome], &dense(), &c);
        let e = analysis.estimates;
        assert_eq!((e.ani, e.q), (None, None));
        assert_eq!(e.ani_status, LocalStatus::InsufficientEvidence);
        assert_eq!(e.af_status, LocalStatus::Estimated);
        assert!(e.af_a.unwrap() > 0.99);
        // No shared scoring probe is never reported as q = 0.
        let counts = LocalCounts {
            scoring_a: 400,
            scoring_b: 400,
            ..LocalCounts::default()
        };
        assert_eq!(survival(&counts, &LocalConfig::default()), None);
        let counts = LocalCounts {
            scoring_shared: 300,
            ..counts
        };
        assert_eq!(survival(&counts, &LocalConfig::default()), Some(0.75));
    }

    #[test]
    fn tiny_inputs_give_statuses_not_nan() {
        let (_, _, analysis) = analyze_pair(
            &[b"ACGTACGTACGTACGTACGTAC"],
            &[b"ACGTACGTACGTACGTACGTAC"],
            &dense(),
            &LocalConfig::default(),
        );
        let e = analysis.estimates;
        assert_eq!((e.ani, e.af_a, e.af_b), (None, None, None));
        assert_eq!(e.af_status, LocalStatus::NoSupportedMapping);
    }

    #[test]
    fn config_checks_refuse_unusable_thresholds() {
        assert_eq!(LocalConfig::default().check(), Ok(()));
        for (bad, message) in [
            (
                LocalConfig {
                    min_block_span: 0,
                    ..LocalConfig::default()
                },
                "at least 1",
            ),
            (
                LocalConfig {
                    max_gap_ratio: 0.5,
                    ..LocalConfig::default()
                },
                "max_gap_ratio",
            ),
            (
                LocalConfig {
                    max_gap_ratio: f64::NAN,
                    ..LocalConfig::default()
                },
                "max_gap_ratio",
            ),
            (
                LocalConfig {
                    min_mapping_hits: 1,
                    ..LocalConfig::default()
                },
                "min_mapping_hits",
            ),
            (
                LocalConfig {
                    min_order_fraction: 1.5,
                    ..LocalConfig::default()
                },
                "min_order_fraction",
            ),
            (
                LocalConfig {
                    min_scoring_probes: 0,
                    ..LocalConfig::default()
                },
                "min_scoring_probes",
            ),
        ] {
            assert!(bad.check().unwrap_err().contains(message), "{bad:?}");
        }
        assert_eq!(
            LocalStatus::InsufficientEvidence.name(),
            "insufficient_evidence"
        );
    }

    #[test]
    fn interval_unions_count_each_base_once() {
        assert_eq!(covered(vec![]), 0);
        assert_eq!(
            covered(vec![
                (0, 10, 20),
                (0, 15, 30),
                (0, 30, 35),
                (1, 0, 5),
                (0, 40, 41)
            ]),
            26 + 5
        );
    }
}
