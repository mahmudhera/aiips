use std::io::{self, Write};
use std::num::NonZeroUsize;
use std::path::Path;

use aiips::fastx::Summary;
use aiips::sketch::minhash::Estimates;

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

/// AF is not estimable from a global MinHash, so both AF columns are -1.
/// Containment-ANI columns are appended so existing column positions hold.
pub fn comparison(
    output: &mut impl Write,
    inputs: &[impl AsRef<Path>; 2],
    e: &Estimates,
) -> io::Result<()> {
    writeln!(output, "\ninput_a\tinput_b\tjaccard\tcontainment_a_in_b\tcontainment_b_in_a\tani\taf_a_to_b\taf_b_to_a\tani_containment_a_in_b\tani_containment_b_in_a")?;
    writeln!(
        output,
        "{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t-1\t-1\t{:.6}\t{:.6}",
        label(inputs[0].as_ref()),
        label(inputs[1].as_ref()),
        e.jaccard,
        e.containment_a_in_b,
        e.containment_b_in_a,
        e.ani,
        e.ani_containment_a_in_b,
        e.ani_containment_b_in_a
    )
}
