use std::io::{self, Write};
use std::num::NonZeroUsize;
use std::path::Path;

use aiips::fastx::Summary;

pub fn label(path: &Path) -> String {
    path.to_string_lossy()
        .replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\r', "\\r")
        .replace('\n', "\\n")
}

pub fn statistics(
    output: &mut impl Write,
    inputs: &[(&Path, Option<&Summary>)],
    k: NonZeroUsize,
) -> io::Result<()> {
    writeln!(
        output,
        "input\trecords\tbases\tnon_acgt_bases\tk\tvalid_kmers"
    )?;
    for (path, summary) in inputs {
        let name = label(path);
        if let Some(summary) = summary {
            writeln!(
                output,
                "{name}\t{}\t{}\t{}\t{k}\t{}",
                summary.records, summary.bases, summary.non_acgt_bases, summary.kmers
            )?;
        } else {
            writeln!(output, "{name}\tNA\tNA\tNA\tNA\tNA")?;
        }
    }
    Ok(())
}

pub fn comparison(output: &mut impl Write, inputs: &[impl AsRef<Path>; 2]) -> io::Result<()> {
    writeln!(output, "\ninput_a\tinput_b\tjaccard\tcontainment_a_in_b\tcontainment_b_in_a\tani\taf_a_to_b\taf_b_to_a")?;
    writeln!(
        output,
        "{}\t{}\t-1\t-1\t-1\t-1\t-1\t-1",
        label(inputs[0].as_ref()),
        label(inputs[1].as_ref())
    )
}
