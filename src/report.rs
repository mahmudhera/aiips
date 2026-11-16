use std::io::{self, Write};
use std::num::NonZeroUsize;
use std::path::Path;

use aiips::fastx::Summary;
use aiips::sketch::aiips::{AiipsComparison, AiipsSketch};
use aiips::sketch::Estimates;

pub fn label(path: &Path) -> String {
    path.to_string_lossy()
        .replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\r', "\\r")
        .replace('\n', "\\n")
}

pub fn statistics(
    output: &mut impl Write,
    inputs: &[(&Path, &Summary, NonZeroUsize)],
) -> io::Result<()> {
    writeln!(
        output,
        "input\trecords\tbases\tnon_acgt_bases\tk\tvalid_kmers"
    )?;
    for (path, summary, k) in inputs {
        writeln!(
            output,
            "{}\t{}\t{}\t{}\t{k}\t{}",
            label(path),
            summary.records,
            summary.bases,
            summary.non_acgt_bases,
            summary.kmers
        )?;
    }
    Ok(())
}

/// AF is not estimable from a global sketch, so both AF columns are -1.
/// Containment-ANI columns are appended so existing column positions hold.
pub fn comparison(
    output: &mut impl Write,
    inputs: &[impl AsRef<Path>; 2],
    e: &Estimates,
) -> io::Result<()> {
    comparison_row(
        output,
        inputs,
        [e.jaccard, e.containment_a_in_b, e.containment_b_in_a],
        [Some(e.ani), None, None],
        [e.ani_containment_a_in_b, e.ani_containment_b_in_a],
    )
}

/// An estimate with six decimals, or -1 when it is not available.
fn value(estimate: Option<f64>) -> String {
    estimate.map_or_else(|| "-1".into(), |v| format!("{v:.6}"))
}

/// The comparison table shared by every algorithm. `None` in `ani`, `af_a_to_b`
/// or `af_b_to_a` (not estimated) prints as -1.
fn comparison_row(
    output: &mut impl Write,
    inputs: &[impl AsRef<Path>; 2],
    [jaccard, containment_a_in_b, containment_b_in_a]: [f64; 3],
    [ani, af_a, af_b]: [Option<f64>; 3],
    [ani_containment_a_in_b, ani_containment_b_in_a]: [f64; 2],
) -> io::Result<()> {
    writeln!(output, "\ninput_a\tinput_b\tjaccard\tcontainment_a_in_b\tcontainment_b_in_a\tani\taf_a_to_b\taf_b_to_a\tani_containment_a_in_b\tani_containment_b_in_a")?;
    writeln!(
        output,
        "{}\t{}\t{jaccard:.6}\t{containment_a_in_b:.6}\t{containment_b_in_a:.6}\t{}\t{}\t{}\t{ani_containment_a_in_b:.6}\t{ani_containment_b_in_a:.6}",
        label(inputs[0].as_ref()),
        label(inputs[1].as_ref()),
        value(ani),
        value(af_a),
        value(af_b),
    )
}

/// AIIPS build diagnostics, printed after the statistics table.
pub fn aiips_build(
    output: &mut impl Write,
    input: &Path,
    sketch: &AiipsSketch,
    file_bytes: u64,
) -> io::Result<()> {
    writeln!(
        output,
        "\ninput\tanchor_occurrences\tdistinct_anchors\tunique_anchors\tprobe_occurrences\tdistinct_probes\tambiguous_kmers\tfile_bytes"
    )?;
    let anchors = sketch.anchors();
    writeln!(
        output,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{file_bytes}",
        label(input),
        anchors.len(),
        anchors.hashes.len(),
        anchors.unique(),
        sketch.probes().len(),
        sketch.global_probes().len(),
        sketch.ambiguous_kmers()
    )
}

/// Header of the AIIPS diagnostics table.
pub const AIIPS_DIAGNOSTICS: &str = "input_a\tinput_b\tglobal_a\tglobal_b\tglobal_shared\tmax_containment\tputative_ani\tshared_anchors\trepeated_anchors\tbreaks_order\tbreaks_gap\tbreaks_ratio\tcandidate_blocks\toverlapping_blocks\trejected_few_hits\trejected_disorder\taccepted_blocks\tmapping_hits\tscoring_a\tscoring_b\tscoring_shared\texcluded_boundary_overlap\texcluded_repeats\tcovered_bases_a\tcovered_bases_b\taf_containment_a\taf_containment_b\tjaccard_status\tcontainment_status\tani_status\taf_status";

/// Header of the AIIPS diagnostics table for `--local-method runs`.
pub const AIIPS_RUN_DIAGNOSTICS: &str = "input_a\tinput_b\tglobal_a\tglobal_b\tglobal_shared\tmax_containment\tputative_ani\tlocal_method\tseeds\tseeds_absorbed\truns_built\tstops_gap\tstops_end\tambiguous_supports\truns_short\truns_trimmed\truns_dropped\taccepted_runs\tsupports\tprobes_a\tprobes_b\tmatched_a\tmatched_b\tcovered_bases_a\tcovered_bases_b\taf_containment_a\taf_containment_b\tjaccard_status\tcontainment_status\tani_status\taf_status";

/// AIIPS comparison: the shared table with positional ANI and AF (-1 when
/// unsupported), then counts for every stage and per-estimate statuses.
pub fn aiips_comparison(
    output: &mut impl Write,
    inputs: &[impl AsRef<Path>; 2],
    c: &AiipsComparison,
) -> io::Result<()> {
    let (g, l) = (&c.global, &c.local);
    comparison_row(
        output,
        inputs,
        [g.jaccard, g.containment_a_in_b, g.containment_b_in_a],
        [l.ani, l.af_a, l.af_b],
        [g.ani_containment_a_in_b, g.ani_containment_b_in_a],
    )?;
    if let Some(r) = &l.runs {
        writeln!(output, "\n{AIIPS_RUN_DIAGNOSTICS}")?;
        return writeln!(
            output,
            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\truns\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\testimated\testimated\t{}\t{}",
            label(inputs[0].as_ref()),
            label(inputs[1].as_ref()),
            g.global_a,
            g.global_b,
            g.global_shared,
            g.max_containment,
            g.putative_ani,
            r.seeds,
            r.seeds_absorbed,
            r.runs_built,
            r.stops_gap,
            r.stops_end,
            r.ambiguous,
            r.runs_short,
            r.runs_trimmed,
            r.runs_dropped,
            r.accepted_runs,
            r.supports,
            r.probes_a,
            r.probes_b,
            r.matched_a,
            r.matched_b,
            r.covered_bases_a,
            r.covered_bases_b,
            value(c.af_containment_a),
            value(c.af_containment_b),
            l.ani_status.name(),
            l.af_status.name(),
        );
    }
    let n = &l.counts;
    let m = &n.matching;
    writeln!(output, "\n{AIIPS_DIAGNOSTICS}")?;
    writeln!(
        output,
        "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\testimated\testimated\t{}\t{}",
        label(inputs[0].as_ref()),
        label(inputs[1].as_ref()),
        g.global_a,
        g.global_b,
        g.global_shared,
        g.max_containment,
        g.putative_ani,
        m.shared_anchors,
        m.repeated_anchors,
        m.breaks_order,
        m.breaks_gap,
        m.breaks_ratio,
        m.candidate_blocks,
        m.overlapping_blocks,
        n.rejected_few_hits,
        n.rejected_disorder,
        n.accepted_blocks,
        n.mapping_hits,
        n.scoring_a,
        n.scoring_b,
        n.scoring_shared,
        n.excluded_boundary_overlap,
        n.excluded_repeats,
        n.covered_bases_a,
        n.covered_bases_b,
        value(c.af_containment_a),
        value(c.af_containment_b),
        l.ani_status.name(),
        l.af_status.name(),
    )
}
