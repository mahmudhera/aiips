//! Conservative dispatch for heterogeneous identity. Statistical working models
//! and thresholds were selected on development simulations, not calibrated as
//! confidence tests. The original estimator remains the fallback.
use super::{
    blocks::Candidates, decoded, identity_mixture, word_mixture, AiipsEstimates, AiipsSketch,
    LocalEstimates,
};

pub(super) fn adjust(
    a: &AiipsSketch,
    b: &AiipsSketch,
    candidates: &Candidates,
    global: &AiipsEstimates,
    base: &mut LocalEstimates,
) {
    let Some(q) = base.q else {
        return;
    };
    let Some(previous) = base.ani else {
        return;
    };
    if !a.params.canonical || a.params.hash != "murmur3-fmix64" {
        return;
    }
    let shared = base.counts.matching.shared_anchors;
    // Small samples cannot support this extra estimator. The approximation
    // below also avoids re-counting unique anchors on ordinary comparisons.
    if shared < 1000 || a.summary.bases.min(b.summary.bases) < 500_000 {
        return;
    }
    let exponent =
        (a.params.anchor_k.get() - a.params.probe_k.get()) as f64 / a.params.probe_k.get() as f64;
    let scale = global.max_containment * q.powf(exponent);
    if !scale.is_finite() || scale <= 0. {
        return;
    }
    let fragmented = a.contigs().len() > 1 || b.contigs().len() > 1;
    let r_min = if fragmented { 1.015 } else { 1.008 };
    let gap_min = if fragmented { 0.0025 } else { 0.0015 };
    let raw_c21 = (shared + base.counts.matching.repeated_anchors) as f64
        / a.anchors()
            .hashes
            .len()
            .min(b.anchors().hashes.len())
            .max(1) as f64;
    let dispersion = if fragmented {
        decoded::fragment_dispersion(a, b, candidates)
    } else {
        0.
    };
    if raw_c21 / scale <= 1.008 && dispersion <= 2. {
        return;
    }
    let unique = candidates
        .unique_by_contig
        .iter()
        .map(|v| v.iter().map(|&n| n as usize).sum::<usize>())
        .min()
        .unwrap_or(0);
    let r = shared as f64 / unique.max(1) as f64 / scale;
    let flags = decoded::geometry_flags(candidates);
    let indels =
        flags.small_steps >= 20 && flags.small_steps as f64 > 0.005 * flags.comparisons as f64;
    if indels {
        // Broad blocks keep the displaced matches that exact-diagonal
        // intervals would preferentially lose in divergent indel-rich regions.
        if r <= 1.008 {
            return;
        }
        let result = identity_mixture::estimate_blocks_linked(a, b, candidates, 1000);
        let fit = &result.fits[result.selected_components - 1];
        if result.selected_components > 1 && fit.ani.is_finite() && previous - fit.ani > 0.0005 {
            base.ani = Some(fit.ani);
        }
        return;
    }
    if r <= r_min && dispersion <= 2. {
        return;
    }
    let geometry = decoded::geometry(a, b, candidates, true, true);
    let projected = decoded::estimate(a, b, &geometry);
    let terminal_ok = projected.terminal_all.bases < 500
        || projected.terminal_all.identity().unwrap_or(0.) >= 0.70;
    let interior;
    let evidence = if terminal_ok {
        &projected
    } else {
        let geometry = decoded::geometry(a, b, candidates, true, false);
        interior = decoded::estimate(a, b, &geometry);
        &interior
    };
    if evidence.score_all.bases < 25_000 {
        return;
    }
    let Some(qg) = evidence.q_geometry.filter(|v| *v > 0.) else {
        return;
    };
    let rg = r * (q / qg).powf(exponent);
    if rg <= r_min {
        return;
    }
    let k = a.params.probe_k.get() as u32;
    let Some(raw) = evidence.score_all.identity() else {
        return;
    };
    let f = &evidence.full_pa;
    let Some(theta) = identity_mixture::conditional_mismatch(qg, f.mutated_pairs, f.mismatches, k)
    else {
        return;
    };
    let mut ani = 0.75 * raw + 0.25 * (1. - theta);
    let mut survival = qg;
    let mut homologous = 1.;
    if !terminal_ok {
        // A low-identity terminal population is distinct from weak homology.
        // Fit only here; ordinary comparisons do not pay for this model.
        let result =
            word_mixture::estimate(&projected.full_pa.histogram, k, projected.q_geometry, 0.65);
        let fit = &result.fits[result.selected_components.saturating_sub(1).min(3)];
        if fit.homologous_fraction > 0.5
            && fit.homologous_fraction < 0.97
            && result.observations >= 1000
        {
            ani = fit.ani_hybrid;
            survival = projected.q_geometry.unwrap_or(f64::NAN);
            homologous = fit.homologous_fraction;
        }
    }
    let pool = q.powf(1. / f64::from(k));
    if !ani.is_finite() || pool - ani <= gap_min || previous - ani <= 0.0005 {
        return;
    }
    base.ani = Some(ani);
    if fragmented && survival > 0. {
        base.af_a = base
            .af_a
            .map(|v| (homologous * global.containment_a_in_b / survival).clamp(v, 1.));
        base.af_b = base
            .af_b
            .map(|v| (homologous * global.containment_b_in_a / survival).clamp(v, 1.));
    }
}
