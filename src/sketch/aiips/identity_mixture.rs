//! Conditional mismatch severity and local Poisson mixture deconvolution.
//! These estimates require reliable positional correspondence; callers gate
//! their use and retain the pooled estimator when evidence is insufficient.
//!
//! The mixture is a working model, not a calibrated confidence procedure:
//! homogeneous equal-length bins, rare ideal-hash samples, and reliable fixed
//! geometry give C~Poi(lambda*q), U~Poi(2lambda*(1-q)), q=p^k. Overlapping words,
//! repeats, Ns, and data-dependent geometry violate independence/density.

use crate::sketch::aiips::{blocks::Candidates, AiipsSketch, Role};
use std::collections::BTreeMap;

/// Mutation fraction from independently sampled full-word mismatches:
/// E[K/k]=(1-q)*E[K/k | K>0]. No iid substitution law is needed for this
/// identity, but q and the changed-word severity must concern the same word
/// population, and the projected base correspondence must be correct.
pub(crate) fn conditional_mismatch(
    q: f64,
    mutated_pairs: u64,
    mismatches: u64,
    k: u32,
) -> Option<f64> {
    if !q.is_finite() || !(0.0..=1.0).contains(&q) || k == 0 {
        return None;
    }
    if mutated_pairs == 0 {
        return (q == 1.0).then_some(0.0);
    }
    if mismatches < mutated_pairs {
        return None;
    }
    let severity = mismatches as f64 / (f64::from(k) * mutated_pairs as f64);
    (severity <= 1.0).then_some((1.0 - q) * severity)
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct MixtureFit {
    pub ani: f64,
    log_likelihood: f64,
    pub bic: f64,
}

impl Default for MixtureFit {
    fn default() -> Self {
        Self {
            ani: f64::NAN,
            log_likelihood: f64::NEG_INFINITY,
            bic: f64::INFINITY,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub(crate) struct MixtureEvidence {
    pub fits: [MixtureFit; 3],
    pub selected_components: usize,
}

#[derive(Clone, Copy, Debug)]
struct Bin {
    span: u64,
    na: u64,
    nb: u64,
    shared: u64,
    lambda: f64,
}

#[derive(Clone, Copy)]
struct Observation {
    bin: Bin,
    frequency: f64,
}

/// Exact sufficient-statistic compression: only X, U, lambda, and span enter
/// any fit or posterior-span calculation. Keep lambda's exact bit pattern.
fn histogram(bins: &[Bin]) -> Vec<Observation> {
    let mut counts: BTreeMap<(u64, u64, u64, u64), (Bin, u64)> = BTreeMap::new();
    for &b in bins {
        let key = (
            b.span,
            b.shared,
            b.na + b.nb - 2 * b.shared,
            b.lambda.to_bits(),
        );
        let entry = counts.entry(key).or_insert((b, 0));
        entry.1 += 1;
    }
    counts
        .into_values()
        .map(|(bin, n)| Observation {
            bin,
            frequency: n as f64,
        })
        .collect()
}

impl Bin {
    fn x(self) -> f64 {
        self.shared as f64
    }
    fn u(self) -> f64 {
        (self.na + self.nb - 2 * self.shared) as f64
    }
}

fn bounded_q(q: f64) -> f64 {
    q.clamp(1e-12, 1.0 - 1e-12)
}

fn maximize_q(x: f64, u: f64, lambda: f64) -> f64 {
    if lambda <= 0.0 {
        return bounded_q(x / (x + u).max(1.0));
    }
    let b = lambda - x - u;
    let discriminant = (b * b + 4.0 * lambda * x).sqrt();
    let q = if b >= 0.0 {
        (b + discriminant) / (2.0 * lambda)
    } else {
        2.0 * x / (discriminant - b)
    };
    bounded_q(q)
}

/// Fit a fixed-component mixture with deterministic starts. Constants common
/// to all components and component-count models are omitted from likelihoods.
#[cfg(test)]
fn fit(bins: &[Bin], k: u32, components: usize) -> MixtureFit {
    fit_observations(&histogram(bins), bins.len(), k, components)
}

fn fit_observations(
    observations: &[Observation],
    original_bins: usize,
    k: u32,
    components: usize,
) -> MixtureFit {
    if observations.is_empty() {
        return MixtureFit::default();
    }
    let mut totals = (0.0, 0.0, 0.0);
    for o in observations {
        let b = o.bin;
        totals.0 += b.x() * o.frequency;
        totals.1 += b.u() * o.frequency;
        totals.2 += b.lambda * o.frequency;
    }
    let q0 = maximize_q(totals.0, totals.1, totals.2);
    let p0 = q0.powf(1.0 / f64::from(k));
    let mut best = MixtureFit::default();
    let initializations = if components == 1 { 1 } else { 3 };
    for start in 0..initializations {
        let offsets = match (components, start) {
            (1, _) => [0.0, 0.0, 0.0],
            (2, 0) => [-0.10, 0.04, 0.0],
            (2, 1) => [-0.04, 0.06, 0.0],
            (2, _) => [-0.18, 0.015, 0.0],
            (3, 0) => [-0.12, -0.02, 0.04],
            (3, 1) => [-0.18, -0.05, 0.025],
            _ => [-0.08, 0.01, 0.06],
        };
        let mut qs = [q0; 3];
        let mut weights = [0.0; 3];
        for j in 0..components {
            qs[j] = bounded_q((p0 + offsets[j]).clamp(0.5, 0.9999).powf(f64::from(k)));
            weights[j] = 1.0 / components as f64;
        }
        if components == 1 {
            qs[0] = q0;
        }
        let mut previous = f64::NEG_INFINITY;
        for _ in 0..80 {
            let mut sx = [0.0; 3];
            let mut su = [0.0; 3];
            let mut sl = [0.0; 3];
            let mut sw = [0.0; 3];
            let mut ll = 0.0;
            let mut logq = [0.0; 3];
            let mut logprivate = [0.0; 3];
            let mut logweight = [0.0; 3];
            for j in 0..components {
                logq[j] = qs[j].ln();
                logprivate[j] = (-qs[j]).ln_1p();
                logweight[j] = weights[j].max(1e-15).ln();
            }
            for o in observations {
                let b = o.bin;
                let mut logs = [f64::NEG_INFINITY; 3];
                for j in 0..components {
                    logs[j] =
                        logweight[j] + b.x() * logq[j] + b.u() * logprivate[j] + b.lambda * qs[j];
                }
                let maximum = logs[..components]
                    .iter()
                    .copied()
                    .fold(f64::NEG_INFINITY, f64::max);
                let mut values = [0.0; 3];
                for j in 0..components {
                    values[j] = (logs[j] - maximum).exp();
                }
                let denominator = values[..components].iter().sum::<f64>();
                ll += (maximum + denominator.ln()) * o.frequency;
                for j in 0..components {
                    let r = values[j] / denominator * o.frequency;
                    sw[j] += r;
                    sx[j] += r * b.x();
                    su[j] += r * b.u();
                    sl[j] += r * b.lambda;
                }
            }
            if (ll - previous).abs() < 1e-7 {
                break;
            }
            previous = ll;
            for j in 0..components {
                weights[j] = (sw[j] / original_bins as f64).max(1e-12);
                qs[j] = maximize_q(sx[j], su[j], sl[j]);
            }
            let ws = weights[..components].iter().sum::<f64>();
            for w in &mut weights[..components] {
                *w /= ws;
            }
        }
        // Posterior-weight by covered span, not by number of partial bins.
        let mut posterior_span = [0.0; 3];
        let mut ll = 0.0;
        let mut total_span = 0.0;
        let mut logq = [0.0; 3];
        let mut logprivate = [0.0; 3];
        let mut logweight = [0.0; 3];
        for j in 0..components {
            logq[j] = qs[j].ln();
            logprivate[j] = (-qs[j]).ln_1p();
            logweight[j] = weights[j].max(1e-15).ln();
        }
        for o in observations {
            let b = o.bin;
            let mut logs = [f64::NEG_INFINITY; 3];
            for j in 0..components {
                logs[j] = logweight[j] + b.x() * logq[j] + b.u() * logprivate[j] + b.lambda * qs[j];
            }
            let maximum = logs[..components]
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            let mut values = [0.0; 3];
            for j in 0..components {
                values[j] = (logs[j] - maximum).exp();
            }
            let denominator = values[..components].iter().sum::<f64>();
            ll += (maximum + denominator.ln()) * o.frequency;
            for j in 0..components {
                posterior_span[j] += values[j] / denominator * b.span as f64 * o.frequency;
            }
            total_span += b.span as f64 * o.frequency;
        }
        let mut identities = [f64::NAN; 3];
        for j in 0..components {
            identities[j] = qs[j].powf(1.0 / f64::from(k));
        }
        let ani = (0..components)
            .map(|j| posterior_span[j] / total_span * identities[j])
            .sum();
        let candidate = MixtureFit {
            ani,
            log_likelihood: ll,
            bic: -2.0 * ll + (2 * components - 1) as f64 * (original_bins as f64).ln(),
        };
        if candidate.log_likelihood > best.log_likelihood {
            best = candidate;
        }
    }
    best
}

fn starts(sketch: &AiipsSketch, contig: usize, lo: i64, hi: i64) -> std::ops::Range<usize> {
    let range = sketch.contigs().probes(contig);
    let positions = &sketch.probes().positions[range.clone()];
    let begin = positions.partition_point(|&p| i64::from(p) < lo);
    let end = positions.partition_point(|&p| i64::from(p) < hi);
    range.start + begin..range.start + end
}

fn evidence_from(observations: &[Bin], k: u32) -> MixtureEvidence {
    let mut evidence = MixtureEvidence::default();
    let compressed = histogram(observations);
    for components in 1..=3 {
        evidence.fits[components - 1] =
            fit_observations(&compressed, observations.len(), k, components);
    }
    evidence.selected_components = (0..3)
        .min_by(|&i, &j| evidence.fits[i].bic.total_cmp(&evidence.fits[j].bic))
        .unwrap()
        + 1;
    evidence
}

/// Variant retaining matches across interpolated bin edges: shared hash pairs
/// are placed in A's bin, including their B count; unmatched B words use the
/// linear projection. This avoids losing exact matches to small indel offsets,
/// but assumes the broad block matches orthologous sequence without repeats.
pub(crate) fn estimate_blocks_linked(
    a: &AiipsSketch,
    b: &AiipsSketch,
    candidates: &Candidates,
    bin_size: u32,
) -> MixtureEvidence {
    evidence_from(
        &block_bins(a, b, candidates, bin_size),
        a.params.probe_k.get() as u32,
    )
}

fn count_words(a: &mut [(u64, usize)], b: &mut [(u64, usize)], bins: &mut [Bin]) {
    a.sort_unstable();
    b.sort_unstable();
    let (mut i, mut j) = (0, 0);
    while i < a.len() || j < b.len() {
        let hash = match (a.get(i), b.get(j)) {
            (Some(x), Some(y)) => x.0.min(y.0),
            (Some(x), None) | (None, Some(x)) => x.0,
            _ => unreachable!(),
        };
        let ca = a[i..].iter().take_while(|x| x.0 == hash).count();
        let cb = b[j..].iter().take_while(|x| x.0 == hash).count();
        if ca <= 1 && cb <= 1 {
            if ca == 1 {
                bins[a[i].1].na += 1;
            }
            if ca == 1 && cb == 1 {
                let n = a[i].1;
                bins[n].nb += 1;
                bins[n].shared += 1;
            } else if cb == 1 {
                bins[b[j].1].nb += 1;
            }
        }
        i += ca;
        j += cb;
    }
}

fn block_bins(
    a: &AiipsSketch,
    b: &AiipsSketch,
    candidates: &Candidates,
    bin_size: u32,
) -> Vec<Bin> {
    let k = a.params.probe_k.get() as i64;
    let k1 = a.params.anchor_k.get() as i64;
    let size = i64::from(bin_size.max(k as u32));
    let rho = 1.0 / (2.0 * a.params.probe_scaled as f64);
    let mut out = Vec::new();
    for block in &candidates.blocks {
        let matches = candidates.block_matches(block);
        let (first, last) = (matches[0], matches[matches.len() - 1]);
        let (alo, ahi) = (i64::from(first.x) + k1, i64::from(last.x));
        let (blo, bhi) = (i64::from(first.y) + k1, i64::from(last.y));
        if ahi - alo < k || bhi - blo < k {
            continue;
        }
        let span_a = ahi - alo;
        let span_b = bhi - blo;
        let count = ((span_a + size - 1) / size) as usize;
        let mut bins = Vec::with_capacity(count);
        let mut boundaries = Vec::with_capacity(count);
        for n in 0..count {
            let lo = alo + n as i64 * size;
            let hi = (lo + size).min(ahi);
            let bl = blo + (lo - alo) * span_b / span_a;
            let bh = blo + (hi - alo) * span_b / span_a;
            let positions = (hi.min(ahi - k + 1) - lo).max(0) + (bh.min(bhi - k + 1) - bl).max(0);
            bins.push(Bin {
                span: (hi - lo) as u64,
                na: 0,
                nb: 0,
                shared: 0,
                lambda: positions as f64 * 0.5 * rho,
            });
            boundaries.push((lo, hi, bl, bh));
        }
        let mut wa = vec![Vec::<(u64, usize)>::new(); count];
        let mut wb = wa.clone();
        for i in starts(a, first.contig_a, alo, ahi - k + 1) {
            if a.probe_role(i) != Role::Scoring {
                continue;
            }
            let x = i64::from(a.probes().positions[i]);
            let n = ((x - alo) / size) as usize;
            wa[n].push((a.probes().hash(i), n));
        }
        let blen = i64::from(b.contigs().length(first.contig_b));
        let (low, high) = if first.reverse {
            (blen - bhi, blen - k - blo + 1)
        } else {
            (blo, bhi - k + 1)
        };
        for i in starts(b, first.contig_b, low, high) {
            if b.probe_role(i) != Role::Scoring {
                continue;
            }
            let raw = i64::from(b.probes().positions[i]);
            let y = if first.reverse { blen - k - raw } else { raw };
            // The integer projection assigns starts exactly at a B edge to
            // its following bin. Adjust roundoff in the rational inversion.
            let mut n =
                (((y - blo) * span_a) / (span_b * size)).clamp(0, count as i64 - 1) as usize;
            while n + 1 < count && y >= boundaries[n].3 {
                n += 1;
            }
            while n > 0 && y < boundaries[n].2 {
                n -= 1;
            }
            if y >= boundaries[n].2 {
                wb[n].push((b.probes().hash(i), n));
            }
        }
        let mut aa = wa.into_iter().flatten().collect::<Vec<_>>();
        let mut bb = wb.into_iter().flatten().collect::<Vec<_>>();
        count_words(&mut aa, &mut bb, &mut bins);
        out.extend(bins.into_iter().filter(|b| b.lambda > 0.0));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    fn expected(p: f64) -> Bin {
        let lambda = 10_000.0;
        let q = p.powi(15);
        let shared = (lambda * q).round() as u64;
        let private = (lambda * (1.0 - q)).round() as u64;
        Bin {
            span: 10_000,
            na: shared + private,
            nb: shared + private,
            shared,
            lambda,
        }
    }
    #[test]
    fn conditional_identity_matches_known_counts() {
        assert_eq!(conditional_mismatch(1.0, 0, 0, 15), Some(0.0));
        assert_eq!(conditional_mismatch(0.9, 0, 0, 15), None);
        assert!(
            (conditional_mismatch(0.6, 100, 200, 15).unwrap() - 0.4 * 2.0 / 15.0).abs() < 1e-12
        );
        assert_eq!(conditional_mismatch(0.9, 1, 16, 15), None);
    }
    #[test]
    fn poisson_mstep_recovers_expected_homogeneous_survival() {
        for p in [0.8, 0.85, 0.95, 0.99] {
            let b = expected(p);
            let q = maximize_q(b.x(), b.u(), b.lambda);
            assert!((q.powf(1.0 / 15.0) - p).abs() < 0.0001);
        }
    }
    #[test]
    fn mixture_recovers_exact_two_identity_mean() {
        let bs = (0..100)
            .map(|i| expected(if i % 2 == 0 { 0.85 } else { 0.98 }))
            .collect::<Vec<_>>();
        let one = fit(&bs, 15, 1);
        let two = fit(&bs, 15, 2);
        assert!((two.ani - 0.915).abs() < 0.0001, "{two:?}");
        assert!(two.bic < one.bic);
    }
    #[test]
    fn uniform_expected_counts_do_not_need_extra_components() {
        let bs = vec![expected(0.95); 100];
        let one = fit(&bs, 15, 1);
        let two = fit(&bs, 15, 2);
        let three = fit(&bs, 15, 3);
        assert!(one.bic < two.bic && one.bic < three.bic);
    }
    #[test]
    fn sufficient_statistic_histogram_preserves_fits_and_original_bic_sample_size() {
        let bs = (0..200)
            .map(|i| expected(if i % 3 == 0 { 0.85 } else { 0.98 }))
            .collect::<Vec<_>>();
        let compressed = histogram(&bs);
        assert_eq!(compressed.len(), 2);
        let full = bs
            .iter()
            .map(|&bin| Observation {
                bin,
                frequency: 1.0,
            })
            .collect::<Vec<_>>();
        for components in 1..=3 {
            let a = fit_observations(&compressed, bs.len(), 15, components);
            let b = fit_observations(&full, bs.len(), 15, components);
            assert!((a.ani - b.ani).abs() < 1e-8, "{a:?} versus {b:?}");
            assert!((a.log_likelihood - b.log_likelihood).abs() < 1e-6);
            assert!((a.bic - b.bic).abs() < 1e-6);
        }
    }
    #[test]
    fn linked_word_counts_reassign_shared_target_words_and_mask_repeats() {
        let empty = Bin {
            span: 1000,
            na: 0,
            nb: 0,
            shared: 0,
            lambda: 4.0,
        };
        let mut bins = vec![empty; 2];
        let mut a = vec![(1, 0), (2, 1), (4, 0), (4, 1)];
        let mut b = vec![(1, 1), (2, 1), (3, 0), (4, 0)];
        count_words(&mut a, &mut b, &mut bins);
        assert_eq!((bins[0].na, bins[0].nb, bins[0].shared), (1, 2, 1));
        assert_eq!((bins[1].na, bins[1].nb, bins[1].shared), (1, 1, 1));
        assert_eq!(bins.iter().map(|b| b.shared).sum::<u64>(), 2);
    }
    #[test]
    fn linked_bins_preserve_forward_reverse_and_offset_coordinates() {
        use crate::hash::SplitMix64;
        use crate::sketch::aiips::AiipsParams;
        use std::num::NonZeroUsize;
        let mut rng = SplitMix64::new(20261114);
        let dna = (0..20_037)
            .map(|_| b"ACGT"[(rng.next().unwrap() >> 62) as usize])
            .collect::<Vec<_>>();
        let rc = dna
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect::<Vec<_>>();
        let params = AiipsParams::new(
            NonZeroUsize::new(21).unwrap(),
            125,
            NonZeroUsize::new(15).unwrap(),
            1,
            42,
        )
        .unwrap();
        let build = |d: &[u8]| {
            let mut fa = b">test\n".to_vec();
            fa.extend_from_slice(d);
            fa.push(b'\n');
            AiipsSketch::build(&fa[..], params.clone()).unwrap()
        };
        let a = build(&dna);
        for reverse in [false, true] {
            let mut bdna = b"ACGTACGTACGTACGTACGT".to_vec();
            bdna.extend_from_slice(if reverse { &rc } else { &dna });
            bdna.extend_from_slice(b"AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAA");
            let suffix = bdna.len() - 20 - dna.len();
            let b = build(&bdna);
            let _ = suffix;
            let candidates = crate::sketch::aiips::blocks::candidate_blocks(
                &a,
                &b,
                &crate::sketch::aiips::LocalConfig::default(),
            );
            for size in [1000, 2000, 4000] {
                let bins = block_bins(&a, &b, &candidates, size);
                assert!(!bins.is_empty());
                for b in bins {
                    assert_eq!(b.na, b.nb);
                    assert_eq!(b.na, b.shared);
                }
            }
        }
    }
}
