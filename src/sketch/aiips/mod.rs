//! AIIPS (Anchor Indexed Island Probe Sketching): a genome sketch made of two
//! sampled sets of canonical k-mer occurrences.
//!
//! - **Probes** (k2-mers) measure k-mer overlap and local similarity. An
//!   occurrence is kept iff h2(x) = fmix64(x XOR seed_0) <= T2, with
//!   T_i = floor((2^64 − 1) / c_i). Each probe also has a role bit,
//!   fmix64(x XOR seed_2) & 1: 0 = mapping, 1 = scoring. The bit depends only
//!   on the k-mer, so a shared probe has the same role in both genomes.
//! - **Anchors** (k1-mers, k1 > k2) locate corresponding regions of two
//!   genomes. A k1-mer is kept iff NOT(h2(u) XOR h2(v)) <= T1, where u and v
//!   are its first and last k2-mers: the XOR lies in the top 1/c1 band. h2 is
//!   computed at every position for the probes anyway, so selecting anchors
//!   costs no extra hashing. Under ideal hashing the XOR of the two
//!   independent uniform values is uniform and independent of either one, so
//!   anchors are sampled at rate 1/c1 independently of probe selection. The
//!   rule is strand-symmetric (reverse complementing swaps u and v), and a
//!   k1-mer whose two k2-mers are equal (a short tandem repeat) is never an
//!   anchor. Kept anchors are identified by h1(x) = fmix64(x XOR seed_1),
//!   which is computed for those alone and never collides (fmix64 is a
//!   bijection), so matching compares whole k1-mers.
//!
//! seed_i is SplitMix64 output i of the user seed; x is the canonical 2-bit
//! code. Probe sampling is this crate's FracMinHash rule and hash, so the
//! distinct probe hashes (the global probe set used for Jaccard and
//! containment) equal a FracMinHash sketch with k = k2, scaled = c2.
//!
//! Every retained occurrence is kept with its contig, position and strand
//! (1 when the reverse complement is canonical), including probes that
//! overlap anchors. k1 and k2 must be odd, so no k-mer equals its reverse
//! complement and every strand is defined; at most 31 bases fit one u64 code,
//! and fmix64 is a bijection, so distinct k-mers never share a hash.
//!
//! **Islands.** On each contig, anchor `a` owns the probes that start at or
//! after it and before the next anchor; the contig's last anchor owns the
//! rest (its island is closed by END). Probes before the first anchor belong
//! to the contig's START bucket. Ownership is derived from positions by one
//! linear merge, both after building and after loading a file.
//!
//! **Comparison.** Jaccard and containment use the global probe sets
//! ([`estimate`]). Positional ANI and AF come from shared anchors, candidate
//! blocks validated by mapping probes, and scoring probes ([`blocks`],
//! [`local`]); [`compare`] reports both.

pub mod blocks;
mod format;
pub mod local;
pub mod occurrences;

use std::io::Read;
use std::num::NonZeroUsize;
use std::ops::Range;

use super::fracminhash::{intersection_size, max_hash, HASH};
use super::minhash::{containment_ani, BuildError};
use super::{invalid, ParseError};
use crate::fastx::{summarize_counted, Summary};
use crate::hash::{fmix64, SplitMix64};
use crate::kmer::CODES;
pub use local::{LocalConfig, LocalEstimates, LocalStatus};
use occurrences::link_islands;
use occurrences::OccurrenceBuilder;
pub use occurrences::{Bits, Contigs, Occurrences};

pub const DEFAULT_ANCHOR_K: usize = 21;
pub const DEFAULT_ANCHOR_SCALED: u64 = 125;
pub const DEFAULT_PROBE_K: usize = 15;
pub const DEFAULT_PROBE_SCALED: u64 = 125;
pub const DEFAULT_SEED: u64 = 42;
/// Largest odd k whose 2-bit encoding fits in a u64.
pub const MAX_K: usize = 31;
pub const MIN_SCALED: u64 = 1;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AiipsParams {
    /// Anchor k-mer length k1: odd, 1 to [`MAX_K`].
    pub anchor_k: NonZeroUsize,
    /// Anchor scale factor c1: about 1/c1 of the k1-mers are kept.
    pub anchor_scaled: u64,
    /// Probe k-mer length k2: odd, 1 to [`MAX_K`].
    pub probe_k: NonZeroUsize,
    /// Probe scale factor c2: about 1/c2 of the k2-mers are kept.
    pub probe_scaled: u64,
    pub seed: u64,
    /// Hash family identifier; sketches compare only within one family.
    pub hash: String,
    pub canonical: bool,
}

impl AiipsParams {
    /// Parameters for building a sketch with this crate's hash function.
    pub fn new(
        anchor_k: NonZeroUsize,
        anchor_scaled: u64,
        probe_k: NonZeroUsize,
        probe_scaled: u64,
        seed: u64,
    ) -> Result<Self, String> {
        let params = Self {
            anchor_k,
            anchor_scaled,
            probe_k,
            probe_scaled,
            seed,
            hash: HASH.into(),
            canonical: true,
        };
        params.check_ranges().map_err(|e| e.to_string())?;
        Ok(params)
    }

    /// SplitMix64 output `index` of the user seed.
    fn derived_seed(&self, index: usize) -> u64 {
        SplitMix64::new(self.seed).nth(index).unwrap()
    }

    /// Seed of the probe hash h2: SplitMix64 output 0.
    pub fn probe_seed(&self) -> u64 {
        self.derived_seed(0)
    }

    /// Seed of the anchor identity hash h1: SplitMix64 output 1.
    pub fn anchor_seed(&self) -> u64 {
        self.derived_seed(1)
    }

    /// Seed of the probe role bit: SplitMix64 output 2.
    pub fn role_seed(&self) -> u64 {
        self.derived_seed(2)
    }

    /// Anchor selection threshold T1 = floor((2^64 − 1) / anchor_scaled),
    /// applied to NOT(h2(first k2-mer) XOR h2(last k2-mer)).
    pub fn max_anchor_hash(&self) -> u64 {
        max_hash(self.anchor_scaled)
    }

    /// Largest retained probe hash, floor((2^64 − 1) / probe_scaled).
    pub fn max_probe_hash(&self) -> u64 {
        max_hash(self.probe_scaled)
    }

    /// Names of the parameters that differ; sketches compare only if empty.
    pub fn mismatches(&self, other: &Self) -> Vec<String> {
        let mut differences = Vec::new();
        let mut check = |name: &str, a: String, b: String| {
            if a != b {
                differences.push(format!("{name} ({a} vs {b})"));
            }
        };
        check(
            "anchor_k",
            self.anchor_k.to_string(),
            other.anchor_k.to_string(),
        );
        check(
            "anchor_scaled",
            self.anchor_scaled.to_string(),
            other.anchor_scaled.to_string(),
        );
        check(
            "probe_k",
            self.probe_k.to_string(),
            other.probe_k.to_string(),
        );
        check(
            "probe_scaled",
            self.probe_scaled.to_string(),
            other.probe_scaled.to_string(),
        );
        check("seed", self.seed.to_string(), other.seed.to_string());
        check("hash", self.hash.clone(), other.hash.clone());
        check(
            "canonical",
            self.canonical.to_string(),
            other.canonical.to_string(),
        );
        differences
    }

    fn check_ranges(&self) -> Result<(), ParseError> {
        check_k("anchor", self.anchor_k.get() as u64)?;
        check_scaled("anchor_scaled", self.anchor_scaled)?;
        check_k("probe", self.probe_k.get() as u64)?;
        check_scaled("probe_scaled", self.probe_scaled)?;
        check_k_order(self.anchor_k, self.probe_k)
    }

    /// Checks that the builder can produce exactly what these parameters
    /// describe; fields are public, so construction must not trust them.
    fn check_buildable(&self) -> Result<(), BuildError> {
        self.check_ranges().map_err(|e| match e {
            ParseError::Invalid(message) => BuildError::Unsupported(message),
            other => BuildError::Unsupported(other.to_string()),
        })?;
        if self.hash != HASH {
            return Err(BuildError::Unsupported(format!(
                "aiips can only build hash={HASH}, not hash={}",
                self.hash
            )));
        }
        if !self.canonical {
            return Err(BuildError::Unsupported(
                "aiips can only build canonical=true sketches".into(),
            ));
        }
        Ok(())
    }
}

fn check_k(kind: &str, k: u64) -> Result<NonZeroUsize, ParseError> {
    usize::try_from(k)
        .ok()
        .filter(|&k| k <= MAX_K && k % 2 == 1)
        .and_then(NonZeroUsize::new)
        .ok_or_else(|| {
            invalid(format!(
                "aiips {kind} k-mer length must be an odd integer from 1 to {MAX_K}, got {k}"
            ))
        })
}

/// Anchors are selected from their first and last k2-mers, so k1 > k2.
fn check_k_order(anchor_k: NonZeroUsize, probe_k: NonZeroUsize) -> Result<(), ParseError> {
    if anchor_k > probe_k {
        Ok(())
    } else {
        Err(invalid(format!(
            "aiips anchor k-mer length ({anchor_k}) must exceed the probe k-mer length ({probe_k})"
        )))
    }
}

fn check_scaled(field: &str, scaled: u64) -> Result<(), ParseError> {
    if scaled >= MIN_SCALED {
        Ok(())
    } else {
        Err(invalid(format!(
            "{field} must be an integer from {MIN_SCALED} to {}, got {scaled}",
            u64::MAX
        )))
    }
}

/// One pass over a record that samples probes and anchors (see the module
/// documentation). Only the k1 codes are rolled; the k2 window ending at the
/// same base is the low 2·k2 bits of the forward code and the high 2·k2 bits
/// of the reverse-complement code (the newest base enters it at the top).
struct Scanner {
    k1: usize,
    k2: usize,
    /// h2 seed, h1 (anchor identity) seed, role seed.
    seeds: [u64; 3],
    /// T2 (probes), T1 (anchors).
    thresholds: [u64; 2],
}

impl Scanner {
    /// Appends the record's probes and anchors; returns (valid k2-mers,
    /// non-ACGT bases). Positions fit in u32: records are at most u32::MAX.
    #[inline]
    fn scan(
        &self,
        sequence: &[u8],
        anchors: &mut OccurrenceBuilder,
        probes: &mut OccurrenceBuilder,
        roles: &mut Bits,
    ) -> (u64, u64) {
        let (k1, k2) = (self.k1, self.k2);
        let [probe_seed, anchor_seed, role_seed] = self.seeds;
        let [max_probe, max_anchor] = self.thresholds;
        let mask = |k: usize| {
            if k == 32 {
                u64::MAX
            } else {
                (1u64 << (2 * k)) - 1
            }
        };
        let (long_mask, short_mask) = (mask(k1), mask(k2));
        let (shift, drop) = (2 * (k1 as u32 - 1), 2 * (k1 - k2) as u32);
        let gap = k1 - k2; // < 32: both are odd and at most 31
        let (mut forward, mut reverse) = (0u64, 0u64);
        let (mut run, mut windows, mut non_acgt) = (0usize, 0u64, 0u64);
        // h2 of the k2-windows ending at the last 32 bases.
        let mut recent = [0u64; 32];
        for (index, &byte) in sequence.iter().enumerate() {
            let code = CODES[byte as usize];
            if code > 3 {
                run = 0;
                non_acgt += 1;
                continue;
            }
            let code = u64::from(code);
            run += 1;
            forward = ((forward << 2) | code) & long_mask;
            reverse = (reverse >> 2) | ((3 - code) << shift);
            if run < k2 {
                continue;
            }
            windows += 1;
            let (f2, r2) = (forward & short_mask, reverse >> drop);
            let probe = f2.min(r2);
            let hash = fmix64(probe ^ probe_seed);
            recent[index % 32] = hash;
            if hash <= max_probe {
                probes.push(hash, (index + 1 - k2) as u32, r2 < f2);
                roles.push(fmix64(probe ^ role_seed) & 1 == 1);
            }
            // The anchor ending here spans the k2-windows ending `gap` bases
            // ago and here; both are valid once run >= k1.
            if run >= k1 && !(recent[(index - gap) % 32] ^ hash) <= max_anchor {
                let anchor = forward.min(reverse);
                anchors.push(
                    fmix64(anchor ^ anchor_seed),
                    (index + 1 - k1) as u32,
                    reverse < forward,
                );
            }
        }
        (windows, non_acgt)
    }
}

/// What a probe is used for in local comparison.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Role {
    /// Role bit 0: establishes positional correspondence.
    Mapping,
    /// Role bit 1: scores similarity inside established regions.
    Scoring,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AiipsSketch {
    pub params: AiipsParams,
    /// Input statistics; `valid_kmers` counts k2-mer occurrences.
    pub summary: Summary,
    contigs: Contigs,
    anchors: Occurrences,
    /// Probe occurrences; their dictionary is the global probe set.
    probes: Occurrences,
    probe_roles: Bits,
    /// Derived: index of the first probe each anchor owns.
    anchor_probe_starts: Vec<usize>,
}

impl AiipsSketch {
    /// Builds the sketch and summary in one streaming pass over FASTA/FASTQ.
    /// Unsupported parameters are refused before any input is read; records
    /// longer than u32::MAX bases are refused after reading.
    pub fn build<R: Read + Send>(input: R, params: AiipsParams) -> Result<Self, BuildError> {
        params.check_buildable()?;
        let (anchor_k, probe_k) = (params.anchor_k, params.probe_k);
        let (anchor_seed, probe_seed, role_seed) = (
            params.anchor_seed(),
            params.probe_seed(),
            params.role_seed(),
        );
        let (max_anchor, max_probe) = (params.max_anchor_hash(), params.max_probe_hash());
        let mut anchors = OccurrenceBuilder::default();
        let mut probes = OccurrenceBuilder::default();
        let mut roles = Bits::default();
        let mut contigs = Contigs::default();
        let mut too_long = None;
        let scanner = Scanner {
            k1: anchor_k.get(),
            k2: probe_k.get(),
            seeds: [probe_seed, anchor_seed, role_seed],
            thresholds: [max_probe, max_anchor],
        };
        let summary = summarize_counted(input, |sequence| {
            if too_long.is_some() {
                return (0, 0);
            }
            let Ok(length) = u32::try_from(sequence.len()) else {
                too_long = Some(sequence.len());
                return (0, 0);
            };
            let counts = scanner.scan(sequence, &mut anchors, &mut probes, &mut roles);
            contigs.push(length, anchors.len(), probes.len());
            counts
        })
        .map_err(BuildError::Input)?;
        if let Some(length) = too_long {
            return Err(BuildError::Unsupported(format!(
                "aiips supports records of at most {} bases; found one of {length}",
                u32::MAX
            )));
        }
        let too_many = || {
            BuildError::Unsupported(
                "aiips supports at most 2^32 distinct anchor or probe hashes".into(),
            )
        };
        let anchors = anchors.finish().ok_or_else(too_many)?;
        let probes = probes.finish().ok_or_else(too_many)?;
        Ok(Self::assemble(
            params, summary, contigs, anchors, probes, roles,
        ))
    }

    /// Links islands; the parts must already be consistent.
    fn assemble(
        params: AiipsParams,
        summary: Summary,
        contigs: Contigs,
        anchors: Occurrences,
        probes: Occurrences,
        probe_roles: Bits,
    ) -> Self {
        let anchor_probe_starts = link_islands(&contigs, &anchors.positions, &probes.positions);
        Self {
            params,
            summary,
            contigs,
            anchors,
            probes,
            probe_roles,
            anchor_probe_starts,
        }
    }

    /// Distinct retained probe hashes, strictly increasing: the global set.
    pub fn global_probes(&self) -> &[u64] {
        &self.probes.hashes
    }

    pub fn contigs(&self) -> &Contigs {
        &self.contigs
    }

    /// Anchor occurrences, ordered by contig, then position.
    pub fn anchors(&self) -> &Occurrences {
        &self.anchors
    }

    /// Probe occurrences, ordered by contig, then position.
    pub fn probes(&self) -> &Occurrences {
        &self.probes
    }

    pub fn probe_role(&self, probe: usize) -> Role {
        if self.probe_roles.get(probe) {
            Role::Scoring
        } else {
            Role::Mapping
        }
    }

    /// k2-length windows inside contigs that contain a non-ACGT base.
    pub fn ambiguous_kmers(&self) -> u64 {
        self.contigs.windows(self.params.probe_k.get()) - self.summary.kmers
    }

    /// Probes of `contig` before its first anchor (all of them when the
    /// contig has no anchor).
    pub fn start_bucket(&self, contig: usize) -> Range<usize> {
        let probes = self.contigs.probes(contig);
        let anchors = self.contigs.anchors(contig);
        let end = if anchors.is_empty() {
            probes.end
        } else {
            self.anchor_probe_starts[anchors.start]
        };
        probes.start..end
    }

    /// The next anchor on the same contig; `None` is END.
    pub fn next_anchor(&self, anchor: usize) -> Option<usize> {
        let contig = self.contigs.of_anchor(anchor);
        (anchor + 1 < self.contigs.anchors(contig).end).then_some(anchor + 1)
    }

    /// The probes in `anchor`'s island: from its position up to the next
    /// anchor, or to the end of the contig.
    pub fn owned_probes(&self, anchor: usize) -> Range<usize> {
        let start = self.anchor_probe_starts[anchor];
        let end = match self.next_anchor(anchor) {
            Some(next) => self.anchor_probe_starts[next],
            None => self.contigs.probes(self.contigs.of_anchor(anchor)).end,
        };
        start..end
    }

    /// Every probe in the islands of anchors `left..right`: those starting at
    /// or after `left`'s position and before `right`'s. One contiguous slice,
    /// identical to following next-anchor links from `left` to `right`.
    /// `None` unless `left <= right` lie on one contig.
    pub fn probes_between(&self, left: usize, right: usize) -> Option<Range<usize>> {
        (left <= right
            && right < self.anchors.len()
            && self.contigs.of_anchor(left) == self.contigs.of_anchor(right))
        .then(|| self.anchor_probe_starts[left]..self.anchor_probe_starts[right])
    }
}

/// Global set statistics of two probe sketches; positional ANI and AF are in
/// [`LocalEstimates`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AiipsEstimates {
    /// Distinct retained probe hashes of A, of B, and shared by both.
    pub global_a: u64,
    pub global_b: u64,
    pub global_shared: u64,
    pub jaccard: f64,
    pub containment_a_in_b: f64,
    pub containment_b_in_a: f64,
    /// shared / min(|A|, |B|), the larger of the two containments.
    pub max_containment: f64,
    pub ani_containment_a_in_b: f64,
    pub ani_containment_b_in_a: f64,
    /// max_containment^(1/k2): a screening value, not a positional ANI.
    pub putative_ani: f64,
}

/// Compares two non-empty sketches built with identical parameters, using
/// the distinct global probe hashes (never occurrence counts):
///
/// - Jaccard: x / (a + b − x), with a = |A|, b = |B|, x = |A ∩ B|.
/// - Containment of A in B: x / a; of B in A: x / b; maximum: x / min(a, b).
///   These are the raw FracMinHash ratios (no debiasing).
/// - ANI columns: C^(1/k2) for each containment and for the maximum. Under the
///   model of Hera, Pierce-Ward & Koslicki (2023) this assumes homologous,
///   fully contained sequence with independent substitutions and no repeats;
///   it is for screening only and is invalid for incomplete or rearranged
///   genomes.
pub fn estimate(a: &AiipsSketch, b: &AiipsSketch) -> Result<AiipsEstimates, String> {
    let differences = a.params.mismatches(&b.params);
    if !differences.is_empty() {
        return Err(format!(
            "sketch parameters differ: {}",
            differences.join(", ")
        ));
    }
    let (global_a, global_b) = (a.global_probes(), b.global_probes());
    if global_a.is_empty() || global_b.is_empty() {
        return Err(EMPTY.into());
    }
    let k = a.params.probe_k.get();
    let shared = intersection_size(global_a, global_b) as u64;
    let (size_a, size_b) = (global_a.len() as u64, global_b.len() as u64);
    let ratio = |n: u64, d: u64| n as f64 / d as f64;
    let containment_a_in_b = ratio(shared, size_a);
    let containment_b_in_a = ratio(shared, size_b);
    let max_containment = ratio(shared, size_a.min(size_b));
    Ok(AiipsEstimates {
        global_a: size_a,
        global_b: size_b,
        global_shared: shared,
        jaccard: ratio(shared, size_a + size_b - shared),
        containment_a_in_b,
        containment_b_in_a,
        max_containment,
        ani_containment_a_in_b: containment_ani(containment_a_in_b, k),
        ani_containment_b_in_a: containment_ani(containment_b_in_a, k),
        putative_ani: containment_ani(max_containment, k),
    })
}

/// Global and positional comparison of two sketches.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AiipsComparison {
    pub global: AiipsEstimates,
    pub local: LocalEstimates,
    /// Containment / q, i.e. C / ANI^k2, for A in B and B in A: under a model
    /// where a fraction f of k2-mers is homologous and each survives with
    /// probability ANI^k2, C = f · ANI^k2. A k-mer homologous fraction, not an
    /// aligned-base fraction; unclipped, so it can exceed 1 under repeats or
    /// incompleteness. `None` when ANI is not estimated.
    pub af_containment_a: Option<f64>,
    pub af_containment_b: Option<f64>,
}

/// [`estimate`] plus positional ANI and AF ([`local`]). Refuses the same
/// inputs as [`estimate`].
pub fn compare(
    a: &AiipsSketch,
    b: &AiipsSketch,
    config: &LocalConfig,
) -> Result<AiipsComparison, String> {
    let global = estimate(a, b)?;
    let local = local::estimate_local(a, b, config);
    let proxy = |containment: f64| local.q.map(|q| containment / q);
    Ok(AiipsComparison {
        af_containment_a: proxy(global.containment_a_in_b),
        af_containment_b: proxy(global.containment_b_in_a),
        global,
        local,
    })
}

/// Reason an empty probe sketch cannot be compared.
pub const EMPTY: &str =
    "aiips sketch has no probe hash values; containment is undefined; use a smaller --probe-scaled";

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::sketch::fracminhash::{FracMinHashParams, FracMinHashSketch};
    use std::collections::{BTreeSet, HashSet};

    pub(crate) fn params(
        anchor_k: usize,
        anchor_scaled: u64,
        probe_k: usize,
        probe_scaled: u64,
        seed: u64,
    ) -> AiipsParams {
        let k = |k| NonZeroUsize::new(k).unwrap();
        AiipsParams::new(k(anchor_k), anchor_scaled, k(probe_k), probe_scaled, seed).unwrap()
    }

    /// Probe parameters, with default anchors (or k1 = k2 + 2 for long k2).
    fn probe_params(k: usize, scaled: u64, seed: u64) -> AiipsParams {
        params(
            DEFAULT_ANCHOR_K.max(k + 2),
            DEFAULT_ANCHOR_SCALED,
            k,
            scaled,
            seed,
        )
    }

    pub(crate) fn build(fasta: &[u8], p: AiipsParams) -> AiipsSketch {
        AiipsSketch::build(fasta, p).unwrap()
    }

    pub(crate) fn fasta(records: &[&[u8]]) -> Vec<u8> {
        let mut out = Vec::new();
        for (index, record) in records.iter().enumerate() {
            out.extend_from_slice(format!(">r{index}\n").as_bytes());
            out.extend_from_slice(record);
            out.push(b'\n');
        }
        out
    }

    /// Deterministic xorshift64* stream for synthetic data in tests.
    pub(crate) fn random_bases(seed: u64, len: usize) -> Vec<u8> {
        let mut state = seed;
        (0..len)
            .map(|_| {
                state ^= state >> 12;
                state ^= state << 25;
                state ^= state >> 27;
                b"ACGT"[(state.wrapping_mul(0x2545_f491_4f6c_dd1d) >> 62) as usize]
            })
            .collect()
    }

    pub(crate) fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
        sequence
            .iter()
            .rev()
            .map(|b| match b.to_ascii_uppercase() {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                other => other,
            })
            .collect()
    }

    /// Independent baseline: (position, canonical code, strand) of every valid
    /// window of one record, encoded from scratch, and its ambiguous windows.
    fn naive_windows(record: &[u8], k: usize) -> (Vec<(u32, u64, bool)>, u64) {
        let encode = |s: &[u8]| {
            s.iter().fold(0u64, |code, b| {
                (code << 2)
                    | match b {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        _ => 3,
                    }
            })
        };
        let upper = record.to_ascii_uppercase();
        let mut hits = Vec::new();
        let mut ambiguous = 0;
        for (position, window) in upper.windows(k).enumerate() {
            if window.iter().all(|b| b"ACGT".contains(b)) {
                let (forward, reverse) = (encode(window), encode(&reverse_complement(window)));
                hits.push((position as u32, forward.min(reverse), reverse < forward));
            } else {
                ambiguous += 1;
            }
        }
        (hits, ambiguous)
    }

    /// Brute-force sampling of one record: h = fmix64(x XOR SplitMix64[seed
    /// index]) kept iff h <= floor((2^64 − 1) / c), computed in u128. Returns
    /// (position, hash, strand, role) with role from SplitMix64[2].
    fn naive_sample(
        record: &[u8],
        k: usize,
        scaled: u64,
        seed: u64,
        index: usize,
    ) -> Vec<(u32, u64, bool, bool)> {
        let seeds: Vec<u64> = SplitMix64::new(seed).take(3).collect();
        let threshold = (u64::MAX as u128 / scaled as u128) as u64;
        naive_windows(record, k)
            .0
            .into_iter()
            .filter_map(|(position, code, strand)| {
                let hash = fmix64(code ^ seeds[index]);
                (hash <= threshold).then(|| {
                    let role = fmix64(code ^ seeds[2]) & 1 == 1;
                    (position, hash, strand, role)
                })
            })
            .collect()
    }

    /// Brute-force anchors of one record: the k1-window at p is kept iff
    /// NOT(h2(k2-window at p) XOR h2(k2-window at p + k1 - k2)) <= T1 (in
    /// u128), identified by fmix64(canonical k1 code XOR SplitMix64[1]).
    fn naive_anchors(
        record: &[u8],
        (k1, c1): (usize, u64),
        k2: usize,
        seed: u64,
    ) -> Vec<(u32, u64, bool)> {
        let seeds: Vec<u64> = SplitMix64::new(seed).take(2).collect();
        let threshold = (u64::MAX as u128 / c1 as u128) as u64;
        let short: std::collections::HashMap<u32, u64> = naive_windows(record, k2)
            .0
            .into_iter()
            .map(|(position, code, _)| (position, code))
            .collect();
        naive_windows(record, k1)
            .0
            .into_iter()
            .filter(|&(position, _, _)| {
                let first = fmix64(short[&position] ^ seeds[0]);
                let last = fmix64(short[&(position + (k1 - k2) as u32)] ^ seeds[0]);
                !(first ^ last) <= threshold
            })
            .map(|(position, code, strand)| (position, fmix64(code ^ seeds[1]), strand))
            .collect()
    }

    fn anchors_of(sketch: &AiipsSketch, contig: usize) -> Vec<(u32, u64, bool)> {
        let a = sketch.anchors();
        sketch
            .contigs()
            .anchors(contig)
            .map(|i| (a.positions[i], a.hash(i), a.strands.get(i)))
            .collect()
    }

    fn probes_of(sketch: &AiipsSketch, contig: usize) -> Vec<(u32, u64, bool, bool)> {
        let p = sketch.probes();
        sketch
            .contigs()
            .probes(contig)
            .map(|i| {
                let role = sketch.probe_role(i) == Role::Scoring;
                (p.positions[i], p.hash(i), p.strands.get(i), role)
            })
            .collect()
    }

    #[test]
    fn occurrences_match_brute_force_baseline() {
        // The long record is large enough for thousands of occurrences.
        let long = random_bases(11, 60_000);
        let records: Vec<&[u8]> = vec![
            b"ACGTacgtNNacgtRYACGTTGCAAGGCTTAACGTTAGCC",
            b"GATTACA",
            b"nnnn",
            &long,
            b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        ];
        let file = fasta(&records);
        for (k1, c1, k2, c2, seed) in [
            (3, 1, 1, 1, 0),
            (5, 3, 3, 1, 42),
            (21, 125, 15, 125, 42),
            (21, 2, 15, 1, 9),
            (31, 1, 29, 3, u64::MAX),
            (11, u64::MAX, 7, 4, 5),
        ] {
            let sketch = build(&file, params(k1, c1, k2, c2, seed));
            let context = format!("k1={k1} c1={c1} k2={k2} c2={c2} seed={seed}");
            assert_eq!(sketch.contigs().len(), records.len(), "{context}");
            let (mut distinct, mut valid, mut ambiguous) = (BTreeSet::new(), 0, 0);
            for (contig, record) in records.iter().enumerate() {
                assert_eq!(sketch.contigs().length(contig) as usize, record.len());
                let expected = naive_anchors(record, (k1, c1), k2, seed);
                assert_eq!(
                    anchors_of(&sketch, contig),
                    expected,
                    "{context} anchors {contig}"
                );
                let expected = naive_sample(record, k2, c2, seed, 0);
                assert_eq!(
                    probes_of(&sketch, contig),
                    expected,
                    "{context} probes {contig}"
                );
                distinct.extend(expected.iter().map(|o| o.1));
                let (windows, bad) = naive_windows(record, k2);
                valid += windows.len() as u64;
                ambiguous += bad;
            }
            let distinct: Vec<u64> = distinct.into_iter().collect();
            assert_eq!(sketch.global_probes(), distinct, "{context}");
            assert_eq!(sketch.summary.kmers, valid, "{context}");
            assert_eq!(sketch.ambiguous_kmers(), ambiguous, "{context}");
        }
    }

    #[test]
    fn seeds_are_successive_splitmix64_outputs() {
        let p = params(21, 125, 15, 125, 1234567);
        // Published reference outputs for seed 1234567 (Vigna's splitmix64.c).
        assert_eq!(p.probe_seed(), 6457827717110365317);
        assert_eq!(p.anchor_seed(), 3203168211198807973);
        assert_eq!(p.role_seed(), 9817491932198370423);
    }

    #[test]
    fn global_probes_are_bit_identical_to_fracminhash() {
        let file = fasta(&[&random_bases(3, 20_000), b"ACGTNNacgtgattaca"]);
        for (k, scaled, seed) in [(15, 125, 42), (15, 1, 42), (21, 7, 3), (3, 1, 0)] {
            let aiips = build(&file, probe_params(k, scaled, seed));
            let frac = FracMinHashSketch::build(
                file.as_slice(),
                FracMinHashParams::new(NonZeroUsize::new(k).unwrap(), scaled, seed).unwrap(),
            )
            .unwrap();
            assert_eq!(aiips.global_probes(), frac.hashes, "k={k} scaled={scaled}");
            assert_eq!(aiips.summary, frac.summary);
        }
    }

    /// Exact Jaccard and both containments of the distinct canonical k-mer sets.
    fn exact(a: &[&[u8]], b: &[&[u8]], k: usize) -> (f64, f64, f64) {
        let set = |records: &[&[u8]]| -> HashSet<u64> {
            records
                .iter()
                .flat_map(|record| naive_windows(record, k).0)
                .map(|(_, code, _)| code)
                .collect()
        };
        let (a, b) = (set(a), set(b));
        let shared = a.intersection(&b).count() as f64;
        (
            shared / a.union(&b).count() as f64,
            shared / a.len() as f64,
            shared / b.len() as f64,
        )
    }

    #[test]
    fn scaled_one_estimates_are_exact_set_statistics() {
        let base = random_bases(1, 5_000);
        let other = random_bases(2, 5_000);
        let rc = reverse_complement(&base);
        let mut mutated = base.clone();
        for base in mutated.iter_mut().step_by(37) {
            *base = if *base == b'A' { b'C' } else { b'A' };
        }
        // (name, records of A, records of B, known [J, C(A in B), C(B in A)])
        type Case<'a> = (&'a str, Vec<&'a [u8]>, Vec<&'a [u8]>, Option<[f64; 3]>);
        let cases: Vec<Case> = vec![
            ("identical", vec![&base], vec![&base], Some([1.0; 3])),
            ("reverse complement", vec![&base], vec![&rc], Some([1.0; 3])),
            ("unrelated", vec![&base], vec![&other], None),
            ("subsequence", vec![&base[1000..3000]], vec![&base], None),
            ("mutated", vec![&base], vec![&mutated], None),
            (
                "fragmented",
                vec![&base[..1700], &base[1700..4000], &base[4000..]],
                vec![&base],
                None,
            ),
        ];
        for k in [7, 15] {
            for (name, a, b, known) in &cases {
                let p = probe_params(k, 1, 42);
                let (sa, sb) = (build(&fasta(a), p.clone()), build(&fasta(b), p));
                let e = estimate(&sa, &sb).unwrap();
                let (jaccard, a_in_b, b_in_a) = exact(a, b, k);
                let context = format!("{name} k={k}");
                assert_eq!(e.jaccard, jaccard, "{context}");
                assert_eq!(e.containment_a_in_b, a_in_b, "{context}");
                assert_eq!(e.containment_b_in_a, b_in_a, "{context}");
                if let Some(known) = known {
                    assert_eq!([e.jaccard, a_in_b, b_in_a], *known, "{context}");
                }
            }
        }
    }

    #[test]
    fn occurrences_count_repeats_and_windows_stay_inside_records() {
        // One repeated 3-mer: five occurrences, one distinct probe.
        let repeat = build(b">r\nAAAAAAA\n", probe_params(3, 1, 42));
        assert_eq!(
            (repeat.probes().len(), repeat.global_probes().len()),
            (5, 1)
        );
        assert_eq!(repeat.probes().counts(), [5]);
        // ACG and CGT are reverse complements: four occurrences, one probe;
        // the three windows touching N are ambiguous; the 2-base record has none.
        let ambiguous = build(b">a\nACGTNACGT\n>b\nAC\n", probe_params(3, 1, 42));
        assert_eq!(ambiguous.summary.kmers, 4);
        assert_eq!(ambiguous.probes().len(), 4);
        assert_eq!(ambiguous.global_probes().len(), 1);
        assert_eq!(ambiguous.ambiguous_kmers(), 3);
        assert_eq!(ambiguous.probes().positions, [0, 1, 5, 6]);
        // AACGTT holds AAC, ACG, CGT, GTT: canonical {AAC, ACG}. Split into
        // two records, the boundary-spanning ACG and CGT must disappear.
        let joined = build(b">a\nAACGTT\n", probe_params(3, 1, 42));
        let split = build(b">a\nAAC\n>b\nGTT\n", probe_params(3, 1, 42));
        assert_eq!(joined.global_probes().len(), 2);
        assert_eq!(split.global_probes().len(), 1);
        assert_eq!(split.summary.kmers, 2);
        assert_eq!(split.contigs().probes(0), 0..1);
        assert_eq!(split.contigs().probes(1), 1..2);
        assert!(joined.global_probes().contains(&split.global_probes()[0]));
    }

    /// A multi-contig genome with dense anchors (c1 = 20) and probes (c2 = 4),
    /// a contig too short for any anchor, and one too short for anything.
    pub(crate) fn island_genome() -> AiipsSketch {
        let a = random_bases(21, 4_000);
        let b = random_bases(22, 2_500);
        let short = random_bases(23, 18); // < k1, >= k2: probes, no anchors
        let tiny = b"ACGTACG"; // < k2: nothing
        build(&fasta(&[&a, &short, &b, tiny]), params(21, 20, 15, 4, 42))
    }

    #[test]
    fn every_probe_is_assigned_once_to_its_preceding_anchor() {
        let sketch = island_genome();
        let (anchors, probes, contigs) = (sketch.anchors(), sketch.probes(), sketch.contigs());
        assert!(anchors.len() > 200 && probes.len() > 1000);
        assert!(contigs.anchors(1).is_empty() && !contigs.probes(1).is_empty());
        assert!(contigs.probes(3).is_empty());
        let mut owner = vec![None::<Option<usize>>; probes.len()];
        for contig in 0..contigs.len() {
            // Walk START, then each island by next-anchor links until END.
            let mut walked = vec![sketch.start_bucket(contig)];
            let mut anchor = contigs.anchors(contig).next();
            for p in sketch.start_bucket(contig) {
                assert!(owner[p].replace(None).is_none(), "probe {p} assigned twice");
            }
            while let Some(a) = anchor {
                for p in sketch.owned_probes(a) {
                    assert!(
                        owner[p].replace(Some(a)).is_none(),
                        "probe {p} assigned twice"
                    );
                }
                walked.push(sketch.owned_probes(a));
                anchor = sketch.next_anchor(a);
            }
            // The buckets tile the contig's probes in order.
            let mut cursor = contigs.probes(contig).start;
            for range in walked {
                assert_eq!(range.start, cursor);
                cursor = range.end;
            }
            assert_eq!(cursor, contigs.probes(contig).end);
            // Naive owner: the last anchor on the contig at or before the probe.
            for p in contigs.probes(contig) {
                let expected = contigs
                    .anchors(contig)
                    .rev()
                    .find(|&a| anchors.positions[a] <= probes.positions[p]);
                assert_eq!(owner[p], Some(expected), "probe {p}");
            }
        }
        assert!(owner.iter().all(Option::is_some), "every probe assigned");
    }

    #[test]
    fn next_anchor_stops_at_each_contig_end() {
        let sketch = island_genome();
        let contigs = sketch.contigs();
        for contig in 0..contigs.len() {
            let anchors = contigs.anchors(contig);
            for a in anchors.clone() {
                let expected = (a + 1 < anchors.end).then_some(a + 1);
                assert_eq!(sketch.next_anchor(a), expected);
                assert_eq!(contigs.of_anchor(a), contig);
            }
            if anchors.is_empty() {
                assert_eq!(sketch.start_bucket(contig), contigs.probes(contig));
            }
        }
    }

    #[test]
    fn probes_between_nonadjacent_anchors_match_the_island_walk() {
        let sketch = island_genome();
        let (anchors, probes, contigs) = (sketch.anchors(), sketch.probes(), sketch.contigs());
        let mut spanning_several = 0;
        for contig in 0..contigs.len() {
            let range = contigs.anchors(contig);
            for left in range.clone() {
                // Pairs up to 12 anchors apart cover several inner islands.
                for right in left..range.end.min(left + 12) {
                    let slice = sketch.probes_between(left, right).unwrap();
                    // Concatenation of islands left, ..., right - 1.
                    let mut walked = Vec::new();
                    let mut a = left;
                    while a < right {
                        walked.extend(sketch.owned_probes(a));
                        a = sketch.next_anchor(a).unwrap();
                    }
                    assert_eq!(slice.clone().collect::<Vec<_>>(), walked);
                    // Naive: every probe on the contig starting in [left, right).
                    let (lo, hi) = (anchors.positions[left], anchors.positions[right]);
                    let naive: Vec<usize> = contigs
                        .probes(contig)
                        .filter(|&p| (lo..hi).contains(&probes.positions[p]))
                        .collect();
                    assert_eq!(slice.collect::<Vec<_>>(), naive);
                    if right >= left + 3 && !walked.is_empty() {
                        spanning_several += 1;
                    }
                }
            }
        }
        assert!(spanning_several > 100);
        // Across contigs or backwards there is no slice.
        let first_b = contigs.anchors(2).start;
        assert_eq!(sketch.probes_between(0, first_b), None);
        assert_eq!(sketch.probes_between(5, 4), None);
        assert_eq!(sketch.probes_between(0, anchors.len()), None);
    }

    #[test]
    fn reverse_complement_mirrors_positions_and_flips_strands() {
        let genome = random_bases(31, 20_000);
        let p = params(21, 10, 15, 5, 42);
        let (forward, reverse) = (
            build(&fasta(&[&genome]), p.clone()),
            build(&fasta(&[&reverse_complement(&genome)]), p),
        );
        let length = genome.len() as u32;
        let mirror = |k: usize, occurrences: Vec<(u32, u64, bool)>| -> Vec<(u32, u64, bool)> {
            let mut mirrored: Vec<_> = occurrences
                .into_iter()
                .map(|(position, hash, strand)| (length - k as u32 - position, hash, !strand))
                .collect();
            mirrored.reverse();
            mirrored
        };
        assert_eq!(mirror(21, anchors_of(&forward, 0)), anchors_of(&reverse, 0));
        let strip = |o: Vec<(u32, u64, bool, bool)>| -> Vec<(u32, u64, bool)> {
            o.into_iter().map(|(p, h, s, _)| (p, h, s)).collect()
        };
        assert_eq!(
            mirror(15, strip(probes_of(&forward, 0))),
            strip(probes_of(&reverse, 0))
        );
        assert_eq!(forward.global_probes(), reverse.global_probes());
        assert_eq!(forward.anchors().hashes, reverse.anchors().hashes);
    }

    #[test]
    fn roles_depend_only_on_the_kmer_and_split_about_evenly() {
        let genome = random_bases(41, 50_000);
        // Two copies: every probe hash occurs twice, with one role per hash.
        let sketch = build(&fasta(&[&genome, &genome]), probe_params(15, 4, 42));
        let probes = sketch.probes();
        let mut role_of = vec![None; probes.hashes.len()];
        let mut scoring = 0;
        for i in 0..probes.len() {
            let role = sketch.probe_role(i);
            let first = *role_of[probes.ids[i] as usize].get_or_insert(role);
            assert_eq!(role, first);
            scoring += usize::from(role == Role::Scoring);
        }
        // Binomial(n, 1/2): within 4 standard deviations of n/2.
        let n = probes.len() as f64;
        assert!(
            ((scoring as f64) - n / 2.0).abs() < 4.0 * (n / 4.0).sqrt(),
            "{scoring} of {n}"
        );
    }

    /// One contig holding one probe per given (sorted, distinct) hash.
    fn with_probes(hashes: Vec<u64>, params: AiipsParams) -> AiipsSketch {
        let n = hashes.len();
        let mut contigs = Contigs::default();
        contigs.push(n as u32 + 100, 0, n);
        let probes = Occurrences {
            ids: (0..n as u32).collect(),
            positions: (0..n as u32).collect(),
            strands: Bits::from_bytes(&vec![0; n.div_ceil(8)], n).unwrap(),
            hashes,
        };
        let roles = probes.strands.clone();
        let summary = Summary {
            records: 1,
            bases: n as u64 + 100,
            non_acgt_bases: 0,
            kmers: n as u64 + 86,
        };
        AiipsSketch::assemble(
            params,
            summary,
            contigs,
            Occurrences::default(),
            probes,
            roles,
        )
    }

    #[test]
    fn estimates_match_hand_computed_values() {
        let p = probe_params(15, 125, 42);
        let a = with_probes(vec![1, 2, 3, 4], p.clone());
        let b = with_probes(vec![3, 4, 5], p.clone());
        let e = estimate(&a, &b).unwrap();
        assert_eq!((e.global_a, e.global_b, e.global_shared), (4, 3, 2));
        assert_eq!(e.jaccard, 0.4);
        assert_eq!(e.containment_a_in_b, 0.5);
        assert_eq!(e.containment_b_in_a, 2.0 / 3.0);
        assert_eq!(e.max_containment, 2.0 / 3.0);
        assert_eq!(e.ani_containment_a_in_b, 0.5f64.powf(1.0 / 15.0));
        assert_eq!(e.ani_containment_b_in_a, (2.0f64 / 3.0).powf(1.0 / 15.0));
        assert_eq!(e.putative_ani, e.ani_containment_b_in_a);
        let reverse = estimate(&b, &a).unwrap();
        assert_eq!(reverse.containment_a_in_b, e.containment_b_in_a);
        assert_eq!(reverse.max_containment, e.max_containment);
        assert_eq!(reverse.jaccard, e.jaccard);
        // No shared probe: zero set similarity, and the C^(1/k) convention 0.
        let c = with_probes(vec![6, 7], p.clone());
        let none = estimate(&a, &c).unwrap();
        assert_eq!(
            (none.jaccard, none.max_containment, none.putative_ani),
            (0.0, 0.0, 0.0)
        );
    }

    #[test]
    fn estimates_refuse_mismatched_parameters_and_empty_sketches() {
        let p = params(21, 125, 15, 125, 42);
        let a = with_probes(vec![1, 2], p.clone());
        let mut q = params(19, 100, 13, 100, 1);
        q.hash = "other".into();
        q.canonical = false;
        let error = estimate(&a, &with_probes(vec![1], q)).unwrap_err();
        for field in [
            "anchor_k (21 vs 19)",
            "anchor_scaled (125 vs 100)",
            "probe_k (15 vs 13)",
            "probe_scaled (125 vs 100)",
            "seed (42 vs 1)",
            "hash (murmur3-fmix64 vs other)",
            "canonical (true vs false)",
        ] {
            assert!(error.contains(field), "{error}");
        }
        // Anchor parameters alone also block comparison.
        let anchors_differ = with_probes(vec![1, 2], params(23, 125, 15, 125, 42));
        assert!(estimate(&a, &anchors_differ)
            .unwrap_err()
            .contains("anchor_k (21 vs 23)"));
        let empty = with_probes(vec![], p);
        assert_eq!(estimate(&empty, &a).unwrap_err(), EMPTY);
        assert_eq!(estimate(&a, &empty).unwrap_err(), EMPTY);
    }

    #[test]
    fn parameters_are_validated_before_building() {
        let k = |k| NonZeroUsize::new(k).unwrap();
        for (k1, c1, k2, c2, message) in [
            (
                2,
                1,
                15,
                1,
                "anchor k-mer length must be an odd integer from 1 to 31, got 2",
            ),
            (
                33,
                1,
                15,
                1,
                "anchor k-mer length must be an odd integer from 1 to 31, got 33",
            ),
            (21, 0, 15, 1, "anchor_scaled must be an integer from 1"),
            (
                21,
                1,
                16,
                1,
                "probe k-mer length must be an odd integer from 1 to 31, got 16",
            ),
            (21, 1, 15, 0, "probe_scaled must be an integer from 1"),
        ] {
            let error = AiipsParams::new(k(k1), c1, k(k2), c2, 0).unwrap_err();
            assert!(error.contains(message), "{error}");
        }
        let error = AiipsParams::new(k(15), 1, k(15), 1, 0).unwrap_err();
        assert!(error.contains("anchor k-mer length (15) must exceed the probe k-mer length (15)"));
        assert!(AiipsParams::new(k(17), 1, k(21), 1, 0).is_err());
        assert!(AiipsParams::new(k(31), 1, k(29), u64::MAX, 0).is_ok());
        assert!(AiipsParams::new(k(3), u64::MAX, k(1), 1, 0).is_ok());
        let message = |p: AiipsParams| match AiipsSketch::build(b">a\nACGT\n".as_slice(), p) {
            Err(BuildError::Unsupported(message)) => message,
            other => panic!("{other:?}"),
        };
        let mut p = params(21, 125, 15, 125, 42);
        p.hash = "other".into();
        assert!(message(p).contains("can only build hash=murmur3-fmix64"));
        let mut p = params(21, 125, 15, 125, 42);
        p.canonical = false;
        assert!(message(p).contains("canonical=true"));
        let mut p = params(21, 125, 15, 125, 42);
        p.anchor_k = k(20);
        assert!(message(p).contains("anchor k-mer length must be an odd"));
        // Unsupported parameters are refused before the input is read.
        let mut p = params(21, 125, 15, 125, 42);
        p.probe_scaled = 0;
        assert!(matches!(
            AiipsSketch::build(b"not fasta".as_slice(), p),
            Err(BuildError::Unsupported(_))
        ));
    }
}
