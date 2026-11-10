use std::num::NonZeroUsize;

/// A forward-strand k-mer borrowed from its source record.
/// Bases retain their original case; positions are zero-based within the record.
#[derive(Debug, PartialEq, Eq)]
pub struct Kmer<'a> {
    pub position: usize,
    pub bases: &'a [u8],
}

/// Visits every A/C/G/T-only window in O(sequence length) time and O(1) space.
/// Both cases are accepted. Every other byte breaks the valid run.
/// Repeated k-mers are emitted repeatedly; no canonicalization is performed.
pub struct Kmers<'a> {
    sequence: &'a [u8],
    k: usize,
    next_base: usize,
    valid_run: usize,
}

impl<'a> Kmers<'a> {
    pub fn new(sequence: &'a [u8], k: NonZeroUsize) -> Self {
        Self {
            sequence,
            k: k.get(),
            next_base: 0,
            valid_run: 0,
        }
    }
}

impl<'a> Iterator for Kmers<'a> {
    type Item = Kmer<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(&base) = self.sequence.get(self.next_base) {
            self.next_base += 1;
            if matches!(base, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
                self.valid_run += 1;
                if self.valid_run >= self.k {
                    let position = self.next_base - self.k;
                    return Some(Kmer {
                        position,
                        bases: &self.sequence[position..self.next_base],
                    });
                }
            } else {
                self.valid_run = 0;
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn windows(sequence: &[u8], k: usize) -> Vec<(usize, &[u8])> {
        Kmers::new(sequence, NonZeroUsize::new(k).unwrap())
            .map(|kmer| (kmer.position, kmer.bases))
            .collect()
    }

    #[test]
    fn visits_all_overlapping_windows_including_repeats() {
        assert_eq!(
            windows(b"AAAAC", 3),
            vec![(0, &b"AAA"[..]), (1, &b"AAA"[..]), (2, &b"AAC"[..])]
        );
    }

    #[test]
    fn ambiguity_breaks_runs_and_case_is_preserved() {
        assert_eq!(
            windows(b"acGTNACRttt", 3),
            vec![(0, &b"acG"[..]), (1, &b"cGT"[..]), (8, &b"ttt"[..])]
        );
    }

    #[test]
    fn handles_empty_short_and_single_base_windows() {
        assert!(windows(b"", 3).is_empty());
        assert!(windows(b"AC", 3).is_empty());
        assert_eq!(windows(b"ANt", 1), vec![(0, &b"A"[..]), (2, &b"t"[..])]);
        assert_eq!(windows(b"ACG", 3), vec![(0, &b"ACG"[..])]);
    }

    #[test]
    fn agrees_with_naive_window_baseline() {
        // Exhaust all short strings over valid, lowercase, and ambiguous bases.
        let alphabet = b"AcGN";
        for code in 0..4_usize.pow(6) {
            let mut value = code;
            let mut sequence = [0; 6];
            for base in &mut sequence {
                *base = alphabet[value % 4];
                value /= 4;
            }
            for k in 1..=7 {
                let expected: Vec<_> = sequence
                    .windows(k)
                    .enumerate()
                    .filter(|(_, window)| window.iter().all(|b| b"ACGTacgt".contains(b)))
                    .collect();
                assert_eq!(windows(&sequence, k), expected);
            }
        }
    }
}
