//! Binomial mixtures of full decoded-word Hamming counts for separating
//! nonhomologous terminal sequence from supported homologous regions.
//! Correct projected correspondence and representative word sampling are
//! essential. Per-word binomial kernels assume iid substitutions at a latent
//! local rate. They are not appropriate for unmodeled indels or repeats.
//! BIC is a heuristic because mixture tests are nonregular and words overlap.

#[derive(Clone, Copy, Debug)]
pub(crate) struct WordFit {
    #[cfg(test)]
    pub ani_homologous: f64,
    pub ani_hybrid: f64,
    pub homologous_fraction: f64,
    log_likelihood: f64,
    pub bic: f64,
}
impl Default for WordFit {
    fn default() -> Self {
        Self {
            #[cfg(test)]
            ani_homologous: f64::NAN,
            ani_hybrid: f64::NAN,
            homologous_fraction: 0.0,
            log_likelihood: f64::NEG_INFINITY,
            bic: f64::INFINITY,
        }
    }
}
#[derive(Clone, Debug, Default)]
pub(crate) struct WordMixtureEvidence {
    /// Fits use respectively one, two, three, and four latent identities.
    pub fits: [WordFit; 4],
    pub selected_components: usize,
    pub observations: u64,
}

fn fit(
    hist: &[u64],
    k: u32,
    components: usize,
    q_dense: Option<f64>,
    minimum_identity: f64,
) -> WordFit {
    let total = hist[..=k as usize].iter().sum::<u64>() as f64;
    if total == 0.0 {
        return WordFit::default();
    }
    let raw_theta = hist[..=k as usize]
        .iter()
        .enumerate()
        .map(|(j, &n)| j as f64 * n as f64)
        .sum::<f64>()
        / (total * f64::from(k));
    let mut best = WordFit::default();
    for start in 0..if components == 1 { 1 } else { 3 } {
        let mut theta = match (components, start) {
            (1, _) => [raw_theta; 4],
            (2, 0) => [
                (raw_theta * 0.3).max(0.001),
                (raw_theta * 2.0).min(0.9),
                0.0,
                0.0,
            ],
            (2, 1) => [0.05, 0.75, 0.0, 0.0],
            (2, _) => [0.02, 0.15, 0.0, 0.0],
            (3, 0) => [0.02, 0.15, 0.75, 0.0],
            (3, 1) => [
                (raw_theta * 0.3).max(0.001),
                (raw_theta * 1.5).min(0.5),
                0.75,
                0.0,
            ],
            (3, _) => [0.01, 0.05, 0.20, 0.0],
            (4, 0) => [0.01, 0.05, 0.15, 0.75],
            (4, 1) => [
                (raw_theta * 0.25).max(0.001),
                (raw_theta * 0.75).min(0.3),
                (raw_theta * 1.5).min(0.5),
                0.75,
            ],
            _ => [0.02, 0.10, 0.25, 0.7],
        };
        let mut weights = [0.0; 4];
        for j in 0..components {
            theta[j] = theta[j].clamp(1e-8, 1.0 - 1e-8);
            weights[j] = 1.0 / components as f64;
        }
        let mut previous = f64::NEG_INFINITY;
        for _ in 0..120 {
            let mut sw = [0.0; 4];
            let mut sx = [0.0; 4];
            let mut ll = 0.0;
            for (x, &n) in hist[..=k as usize]
                .iter()
                .enumerate()
                .filter(|(_, n)| **n > 0)
            {
                let mut logs = [f64::NEG_INFINITY; 4];
                for j in 0..components {
                    logs[j] = weights[j].max(1e-15).ln()
                        + x as f64 * theta[j].ln()
                        + (f64::from(k) - x as f64) * (-theta[j]).ln_1p();
                }
                let maximum = logs[..components]
                    .iter()
                    .copied()
                    .fold(f64::NEG_INFINITY, f64::max);
                let sum = logs[..components]
                    .iter()
                    .map(|v| (v - maximum).exp())
                    .sum::<f64>();
                ll += (maximum + sum.ln()) * n as f64;
                for j in 0..components {
                    let r = (logs[j] - maximum).exp() / sum * n as f64;
                    sw[j] += r;
                    sx[j] += r * x as f64;
                }
            }
            if (ll - previous).abs() < 1e-7 {
                break;
            }
            previous = ll;
            for j in 0..components {
                weights[j] = (sw[j] / total).max(1e-12);
                theta[j] = (sx[j] / (f64::from(k) * sw[j].max(1e-15))).clamp(1e-8, 1.0 - 1e-8);
            }
            let ws = weights[..components].iter().sum::<f64>();
            for w in &mut weights[..components] {
                *w /= ws;
            }
        }
        let mut ll = 0.0;
        let mut changed_homologous = 0.0;
        let mut mismatches_homologous = 0.0;
        for (x, &n) in hist[..=k as usize]
            .iter()
            .enumerate()
            .filter(|(_, n)| **n > 0)
        {
            let mut logs = [f64::NEG_INFINITY; 4];
            for j in 0..components {
                logs[j] = weights[j].max(1e-15).ln()
                    + x as f64 * theta[j].ln()
                    + (f64::from(k) - x as f64) * (-theta[j]).ln_1p();
            }
            let maximum = logs[..components]
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            let sum = logs[..components]
                .iter()
                .map(|v| (v - maximum).exp())
                .sum::<f64>();
            ll += (maximum + sum.ln()) * n as f64;
            if x > 0 {
                for j in 0..components {
                    if 1.0 - theta[j] >= minimum_identity {
                        let r = (logs[j] - maximum).exp() / sum * n as f64;
                        changed_homologous += r;
                        mismatches_homologous += r * x as f64;
                    }
                }
            }
        }
        let mut identities = [f64::NAN; 4];
        let mut fraction = 0.0;
        for j in 0..components {
            identities[j] = 1.0 - theta[j];
            if identities[j] >= minimum_identity {
                fraction += weights[j];
            }
        }
        #[cfg(test)]
        let ani_homologous = if fraction > 0.0 {
            weights
                .iter()
                .zip(identities)
                .filter(|(_, p)| *p >= minimum_identity)
                .map(|(w, p)| w * p)
                .sum::<f64>()
                / fraction
        } else {
            f64::NAN
        };
        let ani_hybrid =
            if let Some(q) = q_dense.filter(|q| q.is_finite() && (0.0..=1.0).contains(q)) {
                if fraction > 0.0 && changed_homologous > 0.0 {
                    let q_homologous = (q / fraction).clamp(0.0, 1.0);
                    let severity = mismatches_homologous / (f64::from(k) * changed_homologous);
                    1.0 - (1.0 - q_homologous) * severity
                } else if q == 1.0 && fraction > 0.0 {
                    1.0
                } else {
                    f64::NAN
                }
            } else {
                f64::NAN
            };
        let candidate = WordFit {
            #[cfg(test)]
            ani_homologous,
            ani_hybrid,
            homologous_fraction: fraction,
            log_likelihood: ll,
            bic: -2.0 * ll + (2 * components - 1) as f64 * total.ln(),
        };
        if candidate.log_likelihood > best.log_likelihood {
            best = candidate;
        }
    }
    best
}

pub(crate) fn estimate(
    hist: &[u64],
    k: u32,
    q_dense: Option<f64>,
    minimum_identity: f64,
) -> WordMixtureEvidence {
    let mut result = WordMixtureEvidence::default();
    if k == 0 || k > 31 || hist.len() <= k as usize {
        return result;
    }
    result.observations = hist[..=k as usize].iter().sum();
    if result.observations == 0 {
        return result;
    }
    for components in 1..=4 {
        result.fits[components - 1] = fit(hist, k, components, q_dense, minimum_identity);
    }
    result.selected_components = (0..4)
        .min_by(|&i, &j| result.fits[i].bic.total_cmp(&result.fits[j].bic))
        .unwrap()
        + 1;
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    fn expected_hist(parts: &[(f64, f64)], n: f64) -> [u64; 32] {
        let mut hist = [0; 32];
        for (x, value) in hist.iter_mut().enumerate().take(16) {
            let choose = (0..x).fold(1.0, |c, i| c * (15 - i) as f64 / (i + 1) as f64);
            *value = parts
                .iter()
                .map(|&(weight, p)| {
                    n * weight * choose * (1.0 - p).powi(x as i32) * p.powi(15 - x as i32)
                })
                .sum::<f64>()
                .round() as u64;
        }
        hist
    }
    #[test]
    fn exact_uniform_histogram_needs_one_component() {
        let h = expected_hist(&[(1.0, 0.95)], 10_000_000.0);
        let e = estimate(&h, 15, Some(0.95_f64.powi(15)), 0.65);
        assert_eq!(e.selected_components, 1);
        assert!((e.fits[0].ani_hybrid - 0.95).abs() < 1e-6);
    }
    #[test]
    fn novel_component_is_removed_from_identity_and_hybrid_prevalence() {
        let parts = [(0.4, 0.98), (0.4, 0.85), (0.2, 0.25)];
        let h = expected_hist(&parts, 10_000_000.0);
        let q = parts.iter().map(|&(w, p)| w * p.powi(15)).sum();
        let e = estimate(&h, 15, Some(q), 0.65);
        let f = e.fits[2];
        assert!((f.homologous_fraction - 0.8).abs() < 0.0001, "{f:?}");
        assert!((f.ani_homologous - 0.915).abs() < 0.0001, "{f:?}");
        assert!((f.ani_hybrid - 0.915).abs() < 0.0001, "{f:?}");
    }
    #[test]
    fn no_words_are_not_identity_zero() {
        assert_eq!(estimate(&[0; 32], 15, None, 0.65).selected_components, 0);
        assert!(estimate(&[0; 32], 15, None, 0.65).fits[0]
            .ani_hybrid
            .is_nan());
    }
}
