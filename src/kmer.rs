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

/// Longest k for which a 2-bit-encoded k-mer fits in a u64.
pub const MAX_PACKED_K: usize = 32;

#[inline]
fn base_code(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Visits the canonical 2-bit code (A=0, C=1, G=2, T=3, first base most
/// significant) of every A/C/G/T-only window: the smaller of the forward and
/// reverse-complement codes. Same windows as [`Kmers`]; O(1) work per base.
pub struct CanonicalKmers<'a> {
    sequence: &'a [u8],
    k: usize,
    mask: u64,
    shift: u32,
    next_base: usize,
    valid_run: usize,
    forward: u64,
    reverse: u64,
}

impl<'a> CanonicalKmers<'a> {
    /// Returns `None` when `k` exceeds [`MAX_PACKED_K`].
    pub fn new(sequence: &'a [u8], k: NonZeroUsize) -> Option<Self> {
        let k = k.get();
        (k <= MAX_PACKED_K).then(|| Self {
            sequence,
            k,
            mask: if k == MAX_PACKED_K {
                u64::MAX
            } else {
                (1 << (2 * k)) - 1
            },
            shift: 2 * (k as u32 - 1),
            next_base: 0,
            valid_run: 0,
            forward: 0,
            reverse: 0,
        })
    }
}

impl CanonicalKmers<'_> {
    /// Moves to the next valid window and returns its forward and
    /// reverse-complement codes; the window starts at `next_base - k`.
    #[inline]
    fn advance(&mut self) -> Option<(u64, u64)> {
        while let Some(&base) = self.sequence.get(self.next_base) {
            self.next_base += 1;
            let Some(code) = base_code(base) else {
                self.valid_run = 0;
                continue;
            };
            self.forward = ((self.forward << 2) | code) & self.mask;
            self.reverse = (self.reverse >> 2) | ((3 - code) << self.shift);
            self.valid_run += 1;
            if self.valid_run >= self.k {
                return Some((self.forward, self.reverse));
            }
        }
        None
    }
}

impl Iterator for CanonicalKmers<'_> {
    type Item = u64;

    #[inline]
    fn next(&mut self) -> Option<u64> {
        self.advance()
            .map(|(forward, reverse)| forward.min(reverse))
    }
}

/// One canonical k-mer occurrence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StrandedKmer {
    /// Zero-based start within the record.
    pub position: usize,
    /// Canonical 2-bit code, as yielded by [`CanonicalKmers`].
    pub code: u64,
    /// True when the reverse complement is the canonical form (strand 1).
    /// A reverse-complement palindrome (even k only) counts as forward.
    pub reverse: bool,
}

/// [`CanonicalKmers`] with each window's position and strand.
pub struct StrandedKmers<'a>(CanonicalKmers<'a>);

impl<'a> StrandedKmers<'a> {
    /// Returns `None` when `k` exceeds [`MAX_PACKED_K`].
    pub fn new(sequence: &'a [u8], k: NonZeroUsize) -> Option<Self> {
        CanonicalKmers::new(sequence, k).map(Self)
    }
}

impl Iterator for StrandedKmers<'_> {
    type Item = StrandedKmer;

    #[inline]
    fn next(&mut self) -> Option<StrandedKmer> {
        let (forward, reverse) = self.0.advance()?;
        Some(StrandedKmer {
            position: self.0.next_base - self.0.k,
            code: forward.min(reverse),
            reverse: reverse < forward,
        })
    }
}

/// 2-bit code of each byte (A/C/G/T in either case); 4 for every other byte.
pub const CODES: [u8; 256] = {
    let mut codes = [4u8; 256];
    codes[b'A' as usize] = 0;
    codes[b'a' as usize] = 0;
    codes[b'C' as usize] = 1;
    codes[b'c' as usize] = 1;
    codes[b'G' as usize] = 2;
    codes[b'g' as usize] = 2;
    codes[b'T' as usize] = 3;
    codes[b't' as usize] = 3;
    codes
};

#[cfg(test)]
mod tests {
    use super::*;

    /// Independent baseline: encode each window and its reverse complement
    /// from scratch.
    fn naive_canonical(sequence: &[u8], k: usize) -> Vec<u64> {
        let encode = |window: &[u8]| -> u64 {
            window
                .iter()
                .fold(0, |code, &b| (code << 2) | base_code(b).unwrap())
        };
        sequence
            .windows(k)
            .filter(|window| window.iter().all(|b| base_code(*b).is_some()))
            .map(|window| {
                let reverse_complement: Vec<u8> = window
                    .iter()
                    .rev()
                    .map(|b| match b.to_ascii_uppercase() {
                        b'A' => b'T',
                        b'C' => b'G',
                        b'G' => b'C',
                        _ => b'A',
                    })
                    .collect();
                encode(window).min(encode(&reverse_complement))
            })
            .collect()
    }

    /// Independent baseline: (position, canonical code, reverse) of every
    /// valid window, each encoded from scratch.
    fn naive_stranded(sequence: &[u8], k: usize) -> Vec<(usize, u64, bool)> {
        let encode = |window: &[u8]| -> u64 {
            window
                .iter()
                .fold(0, |code, &b| (code << 2) | base_code(b).unwrap())
        };
        sequence
            .windows(k)
            .enumerate()
            .filter(|(_, window)| window.iter().all(|b| base_code(*b).is_some()))
            .map(|(position, window)| {
                let reverse_complement: Vec<u8> = window
                    .iter()
                    .rev()
                    .map(|b| match b.to_ascii_uppercase() {
                        b'A' => b'T',
                        b'C' => b'G',
                        b'G' => b'C',
                        _ => b'A',
                    })
                    .collect();
                let (forward, reverse) = (encode(window), encode(&reverse_complement));
                (position, forward.min(reverse), reverse < forward)
            })
            .collect()
    }

    fn stranded(sequence: &[u8], k: usize) -> Vec<(usize, u64, bool)> {
        StrandedKmers::new(sequence, NonZeroUsize::new(k).unwrap())
            .unwrap()
            .map(|hit| (hit.position, hit.code, hit.reverse))
            .collect()
    }

    #[test]
    fn stranded_kmers_agree_with_naive_baseline() {
        let alphabet = b"AcGtN";
        for code in 0..5_usize.pow(7) {
            let mut value = code;
            let mut sequence = [0; 7];
            for base in &mut sequence {
                *base = alphabet[value % 5];
                value /= 5;
            }
            for k in 1..=8 {
                assert_eq!(
                    stranded(&sequence, k),
                    naive_stranded(&sequence, k),
                    "{sequence:?} k={k}"
                );
            }
        }
        let long = b"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAGGNACGTTGCAACGTTGCAACGTTGCAACGTTGCAT";
        for k in [31, 32] {
            assert_eq!(stranded(long, k), naive_stranded(long, k));
        }
        assert!(StrandedKmers::new(long, NonZeroUsize::new(33).unwrap()).is_none());
    }

    #[test]
    fn strand_flips_and_positions_mirror_under_reverse_complement() {
        // Odd k: no palindromes, so every occurrence flips strand.
        let forward = b"GATTACAGGCTTAACGTTAGCCATG";
        let reverse_complement: Vec<u8> = forward
            .iter()
            .rev()
            .map(|b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect();
        let k = 5;
        let a = stranded(forward, k);
        let mut b = stranded(&reverse_complement, k);
        b.reverse();
        for ((pa, ca, ra), (pb, cb, rb)) in a.into_iter().zip(b) {
            assert_eq!(ca, cb);
            assert_eq!(ra, !rb);
            assert_eq!(pa, forward.len() - k - pb);
        }
    }

    #[test]
    fn code_table_matches_base_code() {
        for byte in 0..=255u8 {
            assert_eq!(
                base_code(byte),
                (CODES[byte as usize] < 4).then(|| u64::from(CODES[byte as usize]))
            );
        }
    }

    fn canonical(sequence: &[u8], k: usize) -> Vec<u64> {
        CanonicalKmers::new(sequence, NonZeroUsize::new(k).unwrap())
            .unwrap()
            .collect()
    }

    #[test]
    fn canonical_codes_agree_with_naive_baseline() {
        let alphabet = b"AcGtN";
        for code in 0..5_usize.pow(7) {
            let mut value = code;
            let mut sequence = [0; 7];
            for base in &mut sequence {
                *base = alphabet[value % 5];
                value /= 5;
            }
            for k in 1..=8 {
                assert_eq!(
                    canonical(&sequence, k),
                    naive_canonical(&sequence, k),
                    "{sequence:?} k={k}"
                );
            }
        }
    }

    #[test]
    fn canonical_codes_handle_full_width_k() {
        let sequence = b"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAGGNACGTTGCAACGTTGCAACGTTGCAACGTTGCAT";
        assert_eq!(canonical(sequence, 32), naive_canonical(sequence, 32));
        assert_eq!(canonical(sequence, 31), naive_canonical(sequence, 31));
        assert!(CanonicalKmers::new(sequence, NonZeroUsize::new(33).unwrap()).is_none());
    }

    #[test]
    fn reverse_complement_yields_the_same_codes() {
        let forward = b"GATTACAGGCTTAACG";
        let reverse_complement = b"CGTTAAGCCTGTAATC";
        let mut a = canonical(forward, 5);
        let mut b = canonical(reverse_complement, 5);
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b);
    }

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
