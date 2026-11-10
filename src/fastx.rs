use std::io::Read;
use std::num::NonZeroUsize;

use needletail::{errors::ParseError, parse_fastx_reader};

use crate::kmer::Kmers;

#[derive(Debug, Default, PartialEq, Eq)]
pub struct Summary {
    pub records: u64,
    pub bases: u64,
    /// Bases other than A/C/G/T, accepting either case.
    pub non_acgt_bases: u64,
    /// Overlapping valid occurrences, not distinct k-mers.
    pub kmers: u64,
}

/// Auto-detect FASTA/FASTQ and gzip. Process one record at a time.
/// FASTQ qualities are validated by the parser but are not used for filtering.
/// Empty input and malformed records return errors; no partial summary is returned.
pub fn summarize<R: Read + Send>(input: R, k: NonZeroUsize) -> Result<Summary, ParseError> {
    let mut reader = parse_fastx_reader(input)?;
    let mut summary = Summary::default();
    while let Some(record) = reader.next() {
        let record = record?;
        let sequence = record.seq();
        summary.records += 1;
        summary.bases += sequence.len() as u64;
        summary.non_acgt_bases += sequence
            .iter()
            .filter(|base| !b"ACGTacgt".contains(base))
            .count() as u64;
        for _kmer in Kmers::new(&sequence, k) {
            summary.kmers += 1;
        }
    }
    Ok(summary)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stats(input: &[u8], k: usize) -> Result<Summary, ParseError> {
        summarize(input, NonZeroUsize::new(k).unwrap())
    }

    #[test]
    fn multiline_fasta_crlf_and_record_boundaries() {
        assert_eq!(
            stats(
                b">one description\r\nACG\r\nTac\r\n>two\r\nNNACG\r\n>short\r\nAC\r\n",
                3
            )
            .unwrap(),
            Summary {
                records: 3,
                bases: 13,
                non_acgt_bases: 2,
                kmers: 5
            }
        );
    }

    #[test]
    fn fastq_counts_match_fasta_without_quality_filtering() {
        let fasta = stats(b">a\nACGTNacg\n>b\nTT\n", 3).unwrap();
        let fastq = stats(b"@a\nACGTNacg\n+\n!!!!!!!!\n@b\nTT\n+\nII\n", 3).unwrap();
        assert_eq!(fasta, fastq);
        assert_eq!(
            fastq,
            Summary {
                records: 2,
                bases: 10,
                non_acgt_bases: 1,
                kmers: 3
            }
        );
    }

    #[test]
    fn accepts_no_final_newline() {
        assert_eq!(stats(b">a\nACGT", 3).unwrap().kmers, 2);
        assert_eq!(stats(b"@a\nACGT\n+\nIIII", 3).unwrap().kmers, 2);
    }

    #[test]
    fn rejects_empty_unknown_and_malformed_inputs() {
        for input in [
            b"".as_slice(),
            b"not fasta",
            b"@a\nACGT\n+\nIII\n",
            b"@a\nACGT\n",
            b"@a\nACGT\n-\nIIII\n",
        ] {
            assert!(stats(input, 3).is_err(), "accepted {input:?}");
        }
    }
}
