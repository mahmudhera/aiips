//! FracMinHash (Irber et al. 2022; Hera, Pierce-Ward & Koslicki 2023): keep
//! every distinct canonical k-mer whose hash falls at or below a threshold, so
//! about a 1/S fraction of the k-mers is retained for scale factor S.
//!
//! The hash is h(x) = fmix64(x XOR seed_0), with seed_0 the first SplitMix64
//! output from the user seed (MinHash's hash function 0). A k-mer is kept iff
//! h(x) <= floor((2^64 − 1) / S), the convention used by sourmash. The keep
//! probability under ideal uniform hashing is (max_hash + 1) / 2^64: exactly
//! 1/S when S is a power of two, and within 2^-64 of 1/S otherwise.

use std::io::Read;
use std::num::NonZeroUsize;

use super::minhash::{containment_ani, mash_ani, BuildError};
use super::{invalid, parse_u64, Estimates, Fields, ParseError};
use crate::fastx::{summarize_with, Summary};
use crate::hash::{fmix64, SplitMix64};
use crate::kmer::{CanonicalKmers, MAX_PACKED_K};

pub const HASH: &str = "murmur3-fmix64";
pub const DEFAULT_K: usize = 21;
pub const DEFAULT_SCALED: u64 = 1000;
pub const DEFAULT_SEED: u64 = 42;
pub const MIN_SCALED: u64 = 1;

/// Buffered hashes are compacted (sorted and deduplicated) once the buffer
/// reaches twice its size after the previous compaction, but not below this.
const MIN_COMPACTION: usize = 4096;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FracMinHashParams {
    pub k: NonZeroUsize,
    /// Scale factor S: about 1/S of the distinct k-mers are kept.
    pub scaled: u64,
    pub seed: u64,
    /// Hash family identifier; sketches compare only within one family.
    pub hash: String,
    pub canonical: bool,
}

impl FracMinHashParams {
    /// Parameters for building a sketch with this crate's hash function.
    pub fn new(k: NonZeroUsize, scaled: u64, seed: u64) -> Result<Self, String> {
        check_k(k.get()).map_err(|e| e.to_string())?;
        check_scaled(scaled).map_err(|e| e.to_string())?;
        Ok(Self {
            k,
            scaled,
            seed,
            hash: HASH.into(),
            canonical: true,
        })
    }

    /// Largest retained hash value, floor((2^64 − 1) / scaled).
    pub fn max_hash(&self) -> u64 {
        max_hash(self.scaled)
    }

    /// Names of the parameters that differ; sketches compare only if empty.
    pub fn mismatches(&self, other: &Self) -> Vec<String> {
        let mut differences = Vec::new();
        let mut check = |name: &str, a: String, b: String| {
            if a != b {
                differences.push(format!("{name} ({a} vs {b})"));
            }
        };
        check("k", self.k.to_string(), other.k.to_string());
        check("scaled", self.scaled.to_string(), other.scaled.to_string());
        check("seed", self.seed.to_string(), other.seed.to_string());
        check("hash", self.hash.clone(), other.hash.clone());
        check(
            "canonical",
            self.canonical.to_string(),
            other.canonical.to_string(),
        );
        differences
    }

    /// Checks that the builder can produce exactly what these parameters
    /// describe; fields are public, so construction must not trust them.
    fn check_buildable(&self) -> Result<(), BuildError> {
        let unsupported = |e: ParseError| match e {
            ParseError::Invalid(message) => BuildError::Unsupported(message),
            other => BuildError::Unsupported(other.to_string()),
        };
        check_k(self.k.get()).map_err(unsupported)?;
        check_scaled(self.scaled).map_err(unsupported)?;
        if self.hash != HASH {
            return Err(BuildError::Unsupported(format!(
                "fracminhash can only build hash={HASH}, not hash={}",
                self.hash
            )));
        }
        if !self.canonical {
            return Err(BuildError::Unsupported(
                "fracminhash can only build canonical=true sketches".into(),
            ));
        }
        Ok(())
    }
}

fn max_hash(scaled: u64) -> u64 {
    u64::MAX / scaled
}

fn check_k(k: usize) -> Result<(), ParseError> {
    if (1..=MAX_PACKED_K).contains(&k) {
        Ok(())
    } else {
        Err(invalid(format!(
            "fracminhash supports k-mer lengths 1 to {MAX_PACKED_K}, got {k}"
        )))
    }
}

fn check_scaled(scaled: u64) -> Result<(), ParseError> {
    if scaled >= MIN_SCALED {
        Ok(())
    } else {
        Err(invalid(format!(
            "scaled must be an integer from {MIN_SCALED} to {}, got {scaled}",
            u64::MAX
        )))
    }
}

/// Retained hashes for one input, compacted periodically so memory stays
/// proportional to the number of distinct retained hashes.
pub struct FracMinHasher {
    seed: u64,
    max_hash: u64,
    hashes: Vec<u64>,
    compact_at: usize,
}

impl FracMinHasher {
    pub fn new(scaled: u64, seed: u64) -> Self {
        Self {
            seed: SplitMix64::new(seed).next().unwrap(),
            max_hash: max_hash(scaled),
            hashes: Vec::new(),
            compact_at: MIN_COMPACTION,
        }
    }

    #[inline]
    pub fn insert(&mut self, kmer: u64) {
        let hash = fmix64(kmer ^ self.seed);
        if hash <= self.max_hash {
            self.hashes.push(hash);
            if self.hashes.len() >= self.compact_at {
                self.compact();
                self.compact_at = (2 * self.hashes.len()).max(MIN_COMPACTION);
            }
        }
    }

    fn compact(&mut self) {
        self.hashes.sort_unstable();
        self.hashes.dedup();
    }

    /// Sorted, distinct retained hashes.
    pub fn into_hashes(mut self) -> Vec<u64> {
        self.compact();
        self.hashes.shrink_to_fit();
        self.hashes
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct FracMinHashSketch {
    pub params: FracMinHashParams,
    pub summary: Summary,
    /// Strictly increasing, each at most `params.max_hash()`.
    pub hashes: Vec<u64>,
}

impl FracMinHashSketch {
    /// Builds the sketch and summary in one streaming pass over FASTA/FASTQ.
    /// Unsupported parameters are refused before any input is read.
    pub fn build<R: Read + Send>(input: R, params: FracMinHashParams) -> Result<Self, BuildError> {
        params.check_buildable()?;
        let mut hasher = FracMinHasher::new(params.scaled, params.seed);
        let k = params.k;
        let summary = summarize_with(input, k, |sequence| {
            // check_buildable guarantees k <= MAX_PACKED_K.
            for kmer in CanonicalKmers::new(sequence, k).expect("validated k") {
                hasher.insert(kmer);
            }
        })
        .map_err(BuildError::Input)?;
        Ok(Self {
            params,
            summary,
            hashes: hasher.into_hashes(),
        })
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        use std::fmt::Write;
        let p = &self.params;
        let mut out = String::with_capacity(200 + 21 * self.hashes.len());
        let _ = writeln!(
            out,
            "{} {} algorithm=fracminhash k={} scaled={} seed={} hash={} canonical={} size={}",
            super::MAGIC,
            super::VERSION,
            p.k,
            p.scaled,
            p.seed,
            p.hash,
            p.canonical,
            self.hashes.len()
        );
        super::write_summary(&mut out, &self.summary);
        for hash in &self.hashes {
            let _ = writeln!(out, "{hash}");
        }
        out.into_bytes()
    }

    pub(crate) fn parse_body<'a>(
        mut header: Fields<'a>,
        mut lines: impl Iterator<Item = &'a str>,
    ) -> Result<Self, ParseError> {
        let raw_k = header.take_u64("k")?;
        let k = usize::try_from(raw_k)
            .ok()
            .and_then(NonZeroUsize::new)
            .filter(|k| k.get() <= MAX_PACKED_K)
            .ok_or_else(|| {
                invalid(format!(
                    "fracminhash supports k-mer lengths 1 to {MAX_PACKED_K}, got {raw_k}"
                ))
            })?;
        let scaled = header.take_u64("scaled")?;
        check_scaled(scaled)?;
        let seed = header.take_u64("seed")?;
        let hash = header.take("hash")?.to_owned();
        let canonical = match header.take("canonical")? {
            "true" => true,
            "false" => false,
            other => {
                return Err(invalid(format!(
                    "canonical must be true or false: {other:?}"
                )))
            }
        };
        let size = header.take_u64("size")?;
        header.finish()?;
        let summary = super::parse_summary(lines.next())?;
        let max_hash = max_hash(scaled);
        // The header's size is untrusted: never preallocate more than 1 Mi values.
        let mut hashes = Vec::with_capacity(size.min(1 << 20) as usize);
        for (index, line) in lines.enumerate() {
            if index as u64 == size {
                return Err(invalid(format!("more than size={size} hash values")));
            }
            let value = parse_u64(line).ok_or_else(|| {
                invalid(format!(
                    "hash value {} is not an unsigned integer: {line:?}",
                    index + 1
                ))
            })?;
            if value > max_hash {
                return Err(invalid(format!(
                    "hash value {} exceeds max_hash={max_hash} for scaled={scaled}",
                    index + 1
                )));
            }
            if hashes.last().is_some_and(|&last| last >= value) {
                return Err(invalid(format!(
                    "hash value {} is not strictly increasing",
                    index + 1
                )));
            }
            hashes.push(value);
        }
        if hashes.len() as u64 != size {
            return Err(invalid(format!(
                "expected size={size} hash values, found {}",
                hashes.len()
            )));
        }
        Ok(Self {
            params: FracMinHashParams {
                k,
                scaled,
                seed,
                hash,
                canonical,
            },
            summary,
            hashes,
        })
    }
}

/// Size of the intersection of two strictly increasing sequences.
pub fn intersection_size(a: &[u64], b: &[u64]) -> usize {
    let (mut i, mut j, mut shared) = (0, 0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                shared += 1;
                i += 1;
                j += 1;
            }
        }
    }
    shared
}

/// Compares two non-empty sketches built with identical parameters.
///
/// - Jaccard: |FRAC(A) ∩ FRAC(B)| / |FRAC(A) ∪ FRAC(B)|.
/// - Containment of A in B: |FRAC(A) ∩ FRAC(B)| / |FRAC(A)|, the raw ratio,
///   reported only when both sketches are non-empty. Hera et al. (2023) score
///   an empty FRAC(A) as zero and divide by 1 − (1 − 1/S)^|A|; neither is done
///   here. In the ideal independent-retention model the ratio is unbiased given
///   a non-empty FRAC(A), and also requiring a non-empty FRAC(B) makes its
///   expectation C / P(FRAC(B) ≠ ∅ | FRAC(A) ≠ ∅) ≥ C (derived here, checked by
///   exhaustive enumeration in the tests).
/// - ANI columns reuse the MinHash transforms (Mash ANI from J, C^(1/k)).
pub fn estimate(a: &FracMinHashSketch, b: &FracMinHashSketch) -> Result<Estimates, String> {
    let differences = a.params.mismatches(&b.params);
    if !differences.is_empty() {
        return Err(format!(
            "sketch parameters differ: {}",
            differences.join(", ")
        ));
    }
    if a.hashes.is_empty() || b.hashes.is_empty() {
        return Err(EMPTY.into());
    }
    let k = a.params.k.get();
    let shared = intersection_size(&a.hashes, &b.hashes) as f64;
    let size_a = a.hashes.len() as f64;
    let size_b = b.hashes.len() as f64;
    let jaccard = shared / (size_a + size_b - shared);
    let containment_a_in_b = shared / size_a;
    let containment_b_in_a = shared / size_b;
    Ok(Estimates {
        jaccard,
        containment_a_in_b,
        containment_b_in_a,
        ani: mash_ani(jaccard, k),
        ani_containment_a_in_b: containment_ani(containment_a_in_b, k),
        ani_containment_b_in_a: containment_ani(containment_b_in_a, k),
    })
}

/// Reason an empty sketch cannot be compared.
pub const EMPTY: &str =
    "fracminhash sketch has no hash values; containment is undefined; use a smaller --scaled";

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    fn params(k: usize, scaled: u64, seed: u64) -> FracMinHashParams {
        FracMinHashParams::new(NonZeroUsize::new(k).unwrap(), scaled, seed).unwrap()
    }

    fn sketch(fasta: &[u8], p: FracMinHashParams) -> FracMinHashSketch {
        FracMinHashSketch::build(fasta, p).unwrap()
    }

    fn complement(b: u8) -> u8 {
        match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        }
    }

    /// Independent baseline: distinct canonical k-mers encoded from scratch.
    fn naive_kmers(sequences: &[&[u8]], k: usize) -> HashSet<u64> {
        let encode = |s: &[u8]| {
            s.iter().fold(0u64, |c, b| {
                (c << 2)
                    | match b {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        _ => 3,
                    }
            })
        };
        let mut distinct = HashSet::new();
        for sequence in sequences {
            let upper = sequence.to_ascii_uppercase();
            for window in upper.windows(k) {
                if window.iter().all(|b| b"ACGT".contains(b)) {
                    let rc: Vec<u8> = window.iter().rev().map(|&b| complement(b)).collect();
                    distinct.insert(encode(window).min(encode(&rc)));
                }
            }
        }
        distinct
    }

    /// Brute-force FracMinHash: hash every distinct k-mer, filter, sort.
    fn naive_hashes(kmers: &HashSet<u64>, scaled: u64, seed: u64) -> Vec<u64> {
        let seed_0 = SplitMix64::new(seed).next().unwrap();
        let threshold = (u64::MAX as u128 / scaled as u128) as u64;
        let mut hashes: Vec<u64> = kmers
            .iter()
            .map(|&x| fmix64(x ^ seed_0))
            .filter(|&h| h <= threshold)
            .collect();
        hashes.sort();
        hashes
    }

    /// Deterministic xorshift64* stream for synthetic data in tests.
    fn random_bases(seed: u64, len: usize) -> Vec<u8> {
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

    fn fasta(parts: &[&[u8]]) -> Vec<u8> {
        let mut out = b">r\n".to_vec();
        for p in parts {
            out.extend_from_slice(p);
            out.push(b'N');
        }
        out.push(b'\n');
        out
    }

    #[test]
    fn streamed_sketch_equals_naive_baseline() {
        let records: [&[u8]; 3] = [
            b"ACGTTGCAnnACGgattacaACGT",
            b"TTTTCCCCAAAAGGGG",
            b"acgtNACGT",
        ];
        let file = b">a\nACGTTGCAnnACG\ngattacaACGT\n>b\nTTTTCCCCAAAAGGGG\n>c\nacgtNACGT\n";
        for (k, scaled, seed) in [
            (1, 1, 0),
            (3, 1, 42),
            (3, 2, 42),
            (5, 3, 7),
            (11, 1, u64::MAX),
            (11, u64::MAX, 5),
        ] {
            let built = sketch(file, params(k, scaled, seed));
            let kmers = naive_kmers(&records, k);
            assert_eq!(
                built.hashes,
                naive_hashes(&kmers, scaled, seed),
                "k={k} scaled={scaled}"
            );
            if scaled == 1 {
                assert_eq!(built.hashes.len(), kmers.len(), "scaled=1 keeps all");
            }
        }
    }

    #[test]
    fn compaction_preserves_the_exact_set() {
        // 50k random 21-mers at scaled=1 cross the compaction threshold many times.
        let bases = random_bases(9, 50_000);
        let built = sketch(&fasta(&[&bases]), params(21, 1, 3));
        let kmers = naive_kmers(&[&bases], 21);
        assert!(kmers.len() > 10 * MIN_COMPACTION);
        assert_eq!(built.hashes, naive_hashes(&kmers, 1, 3));
        for scaled in [2, 7, 100] {
            let built = sketch(&fasta(&[&bases]), params(21, scaled, 3));
            assert_eq!(built.hashes, naive_hashes(&kmers, scaled, 3));
        }
    }

    #[test]
    fn threshold_matches_hand_computed_values() {
        assert_eq!(max_hash(1), u64::MAX);
        assert_eq!(max_hash(2), (1 << 63) - 1);
        assert_eq!(max_hash(3), 6_148_914_691_236_517_205);
        assert_eq!(max_hash(u64::MAX), 1);
        assert_eq!(params(3, 1000, 0).max_hash(), 18_446_744_073_709_551);
    }

    #[test]
    fn identical_disjoint_and_reverse_complement_cases_are_exact() {
        let p = params(5, 1, 42);
        let a = sketch(b">a\nGATTACAGGCTTAACGTTAG\n", p.clone());
        let e = estimate(&a, &a).unwrap();
        assert_eq!(
            e,
            Estimates {
                jaccard: 1.0,
                containment_a_in_b: 1.0,
                containment_b_in_a: 1.0,
                ani: 1.0,
                ani_containment_a_in_b: 1.0,
                ani_containment_b_in_a: 1.0,
            }
        );
        let rc = sketch(b">rc\nCTAACGTTAAGCCTGTAATC\n", p.clone());
        assert_eq!(rc.hashes, a.hashes, "canonical k-mers ignore strand");

        let poly_a = sketch(b">a\nAAAAAAAAAA\n", p.clone());
        let poly_c = sketch(b">c\nCCCCCCCCCC\n", p);
        let disjoint = estimate(&poly_a, &poly_c).unwrap();
        assert_eq!(disjoint.jaccard, 0.0);
        assert_eq!(disjoint.containment_a_in_b, 0.0);
        assert_eq!(disjoint.containment_b_in_a, 0.0);
        assert_eq!(disjoint.ani, 0.0);
        assert_eq!(disjoint.ani_containment_a_in_b, 0.0);
    }

    #[test]
    fn scaled_one_gives_exact_set_statistics() {
        let x = random_bases(1, 3000);
        let y = random_bases(2, 2000);
        let z = random_bases(3, 5000);
        let k = 15;
        let p = params(k, 1, 77);
        let a = sketch(&fasta(&[&x, &y]), p.clone());
        let b = sketch(&fasta(&[&y, &z]), p);
        let sa = naive_kmers(&[&x, &y], k);
        let sb = naive_kmers(&[&y, &z], k);
        let inter = sa.intersection(&sb).count() as f64;
        let union = sa.union(&sb).count() as f64;
        let e = estimate(&a, &b).unwrap();
        assert_eq!(e.jaccard, inter / union);
        assert_eq!(e.containment_a_in_b, inter / sa.len() as f64);
        assert_eq!(e.containment_b_in_a, inter / sb.len() as f64);
    }

    #[test]
    fn formulas_match_hand_computed_values() {
        let with = |hashes: Vec<u64>| FracMinHashSketch {
            params: params(21, 1, 0),
            summary: Summary {
                records: 1,
                bases: 100,
                non_acgt_bases: 0,
                kmers: 80,
            },
            hashes,
        };
        // |A| = 4, |B| = 6, |A ∩ B| = 2: J = 2/8, C(A,B) = 1/2, C(B,A) = 1/3.
        let a = with(vec![1, 2, 3, 4]);
        let b = with(vec![3, 4, 5, 6, 7, 8]);
        assert_eq!(intersection_size(&a.hashes, &b.hashes), 2);
        let e = estimate(&a, &b).unwrap();
        assert_eq!(e.jaccard, 0.25);
        assert_eq!(e.containment_a_in_b, 0.5);
        assert_eq!(e.containment_b_in_a, 1.0 / 3.0);
        assert_eq!(e.ani, mash_ani(0.25, 21));
        assert_eq!(e.ani_containment_a_in_b, 0.5f64.powf(1.0 / 21.0));
        assert_eq!(e.ani_containment_b_in_a, (1.0f64 / 3.0).powf(1.0 / 21.0));
        let swapped = estimate(&b, &a).unwrap();
        assert_eq!(swapped.containment_a_in_b, e.containment_b_in_a);
        assert_eq!(swapped.jaccard, e.jaccard);
    }

    #[test]
    fn empty_sketches_round_trip_but_cannot_be_compared() {
        // Four distinct k-mers at the largest scale keep (almost surely) nothing.
        let p = params(3, u64::MAX, 42);
        let empty = sketch(b">a\nACGTAC\n", p.clone());
        assert!(empty.hashes.is_empty());
        let bytes = empty.to_bytes();
        assert!(String::from_utf8(bytes.clone())
            .unwrap()
            .ends_with(" size=0\nrecords=1 bases=6 non_acgt_bases=0 valid_kmers=4\n"));
        match super::super::parse(&bytes).unwrap() {
            super::super::Sketch::FracMinHash(parsed) => assert_eq!(parsed, empty),
            other => panic!("{other:?}"),
        }
        let full = FracMinHashSketch {
            hashes: vec![1],
            ..empty.clone()
        };
        for (a, b) in [(&empty, &full), (&full, &empty), (&empty, &empty)] {
            assert_eq!(estimate(a, b).unwrap_err(), EMPTY);
        }
    }

    #[test]
    fn mismatched_parameters_are_refused_by_name() {
        let fa = b">a\nACGTACGTAC\n";
        let base = sketch(fa, params(3, 1, 1));
        for (other, field) in [
            (sketch(fa, params(4, 1, 1)), "k (3 vs 4)"),
            (sketch(fa, params(3, 2, 1)), "scaled (1 vs 2)"),
            (sketch(fa, params(3, 1, 2)), "seed (1 vs 2)"),
        ] {
            assert!(estimate(&base, &other).unwrap_err().contains(field));
        }
        let mut other = base.clone();
        other.params.hash = "other".into();
        other.params.canonical = false;
        let error = estimate(&base, &other).unwrap_err();
        assert!(
            error.contains("hash") && error.contains("canonical"),
            "{error}"
        );
    }

    #[test]
    fn build_refuses_parameters_it_cannot_honor() {
        let fa = b">a\nAAAAAA\n>b\nTTTTTT\n";
        let valid = params(3, 1, 1);
        let mut cases = Vec::new();
        let mut p = valid.clone();
        p.canonical = false;
        cases.push((p, "canonical=true"));
        let mut p = valid.clone();
        p.hash = "another-hash".into();
        cases.push((p, "hash=another-hash"));
        let mut p = valid.clone();
        p.k = NonZeroUsize::new(33).unwrap();
        cases.push((p, "1 to 32"));
        let mut p = valid.clone();
        p.scaled = 0;
        cases.push((p, "scaled must be"));
        for (p, expected) in cases {
            match FracMinHashSketch::build(fa.as_slice(), p) {
                Err(BuildError::Unsupported(message)) => {
                    assert!(message.contains(expected), "{message}")
                }
                other => panic!("expected refusal for {expected}: {other:?}"),
            }
        }
        assert!(matches!(
            FracMinHashSketch::build(b"".as_slice(), valid),
            Err(BuildError::Input(_))
        ));
    }

    #[test]
    fn parameters_are_validated() {
        let k = |k| NonZeroUsize::new(k).unwrap();
        assert!(FracMinHashParams::new(k(32), 1, 0).is_ok());
        assert!(FracMinHashParams::new(k(21), u64::MAX, 0).is_ok());
        assert!(FracMinHashParams::new(k(33), 1, 0)
            .unwrap_err()
            .contains("1 to 32"));
        assert!(FracMinHashParams::new(k(21), 0, 0)
            .unwrap_err()
            .contains("scaled must be an integer from 1"));
    }

    #[test]
    fn file_round_trip_and_strict_body() {
        let original = sketch(b">a\nGATTACAGGCTTAACGTTAG\n", params(5, 2, 9));
        assert!(original.hashes.len() >= 3, "{}", original.hashes.len());
        let bytes = original.to_bytes();
        let text = String::from_utf8(bytes.clone()).unwrap();
        let n = original.hashes.len();
        assert!(text.starts_with(&format!(
            "aiips-sketch v1 algorithm=fracminhash k=5 scaled=2 seed=9 hash=murmur3-fmix64 canonical=true size={n}\nrecords=1 bases=20 non_acgt_bases=0 valid_kmers=16\n"
        )));
        match super::super::parse(&bytes).unwrap() {
            super::super::Sketch::FracMinHash(parsed) => assert_eq!(parsed, original),
            other => panic!("{other:?}"),
        }
        let header = text.lines().next().unwrap();
        let summary = text.lines().nth(1).unwrap();
        let values: Vec<&str> = text.lines().skip(2).collect();
        let with = |header: &str, summary: &str, values: &[&str]| {
            let mut file = format!("{header}\n{summary}\n");
            for v in values {
                file.push_str(v);
                file.push('\n');
            }
            super::super::parse(file.as_bytes())
                .map(|_| ())
                .unwrap_err()
                .to_string()
        };
        let size = format!("size={n}");
        let fewer = format!("size={}", n - 1);
        let too_big = (max_hash(2) + 1).to_string();
        let first = values[0];
        let cases = [
            (with(header, summary, &values[..n - 1]), "found"),
            (
                with(header, summary, &[&values[..], &["1"]].concat()),
                "more than",
            ),
            (
                with(&header.replace(&size, &fewer), summary, &values),
                "more than",
            ),
            (
                with(header, summary, &[&values[..n - 1], &["x"]].concat()),
                &format!("hash value {n} is not an unsigned"),
            ),
            (
                with(header, summary, &[&values[..n - 1], &[&too_big]].concat()),
                "exceeds max_hash",
            ),
            (
                with(header, summary, &[&values[..n - 1], &[first]].concat()),
                "not strictly increasing",
            ),
            (
                with(
                    header,
                    summary,
                    &[&values[..n - 1], &[values[n - 2]]].concat(),
                ),
                "not strictly increasing",
            ),
            (
                with(&header.replace("k=5", "k=33"), summary, &values),
                "1 to 32",
            ),
            (
                with(&header.replace("k=5", "k=0"), summary, &values),
                "1 to 32",
            ),
            (
                with(&header.replace("scaled=2", "scaled=0"), summary, &values),
                "scaled must be",
            ),
            (
                with(&header.replace(" scaled=2", ""), summary, &values),
                "missing header field: scaled",
            ),
            (
                with(&header.replace(&format!(" {size}"), ""), summary, &values),
                "missing header field: size",
            ),
            (
                with(&header.replace("seed=9", "seed=+9"), summary, &values),
                "seed",
            ),
            (
                with(
                    &header.replace("canonical=true", "canonical=yes"),
                    summary,
                    &values,
                ),
                "true or false",
            ),
            (
                with(&format!("{header} scaled=2"), summary, &values),
                "duplicate header field: scaled",
            ),
            (
                with(&format!("{header} extra=1"), summary, &values),
                "unknown header field: extra",
            ),
            (
                with(
                    header,
                    &summary.replace("valid_kmers=16", "valid_kmers=0"),
                    &values,
                ),
                "no valid k-mers",
            ),
        ];
        for (error, expected) in cases {
            assert!(
                error.contains(expected),
                "expected {expected:?} in {error:?}"
            );
        }
        // An untrusted huge size does not cause a huge allocation.
        let huge = header.replace(&size, &format!("size={}", u64::MAX));
        assert!(with(&huge, summary, &values).contains("found"));
    }

    /// Exhaustive enumeration in the ideal model, where each distinct k-mer is
    /// kept independently with probability q. Checks the containment
    /// expectations stated in the README under three policies for empty
    /// sketches: zero for an empty FRAC(A) (Hera et al. 2023), FRAC(A)
    /// non-empty only, and both non-empty (what `estimate` accepts).
    #[test]
    fn containment_expectation_depends_on_the_empty_sketch_policy() {
        let summary = Summary {
            records: 1,
            bases: 1,
            non_acgt_bases: 0,
            kmers: 1,
        };
        let frac = |hashes: Vec<u64>| FracMinHashSketch {
            params: params(21, 2, 0),
            summary,
            hashes,
        };
        for (a, b, q) in [
            // The worked example: A = {x, y}, B = {x, z}, q = 1/2.
            (vec![0u64, 1], vec![0u64, 2], 0.5f64),
            (vec![0, 1, 2], vec![1, 2, 3, 4], 0.5),
            (vec![0, 1, 2], vec![1, 2, 3, 4], 1.0 / 3.0),
            (vec![0, 1], vec![2, 3], 0.25),
        ] {
            let universe: Vec<u64> = {
                let mut u: Vec<u64> = a.iter().chain(&b).copied().collect();
                u.sort();
                u.dedup();
                u
            };
            let truth = intersection_size(&a, &b) as f64 / a.len() as f64;
            let (mut zero_policy, mut a_mass, mut a_sum) = (0.0, 0.0, 0.0);
            let (mut both_mass, mut both_sum) = (0.0, 0.0);
            for mask in 0u32..1 << universe.len() {
                let kept = |set: &[u64]| -> Vec<u64> {
                    set.iter()
                        .copied()
                        .filter(|x| mask >> universe.iter().position(|u| u == x).unwrap() & 1 == 1)
                        .collect()
                };
                let n = mask.count_ones() as i32;
                let weight = q.powi(n) * (1.0 - q).powi(universe.len() as i32 - n);
                let (fa, fb) = (kept(&a), kept(&b));
                if fa.is_empty() {
                    continue; // contributes zero under the paper's policy
                }
                let ratio = intersection_size(&fa, &fb) as f64 / fa.len() as f64;
                zero_policy += weight * ratio;
                a_mass += weight;
                a_sum += weight * ratio;
                match estimate(&frac(fa), &frac(fb.clone())) {
                    Ok(e) => {
                        assert_eq!(e.containment_a_in_b, ratio);
                        both_mass += weight;
                        both_sum += weight * ratio;
                    }
                    Err(error) => {
                        assert!(fb.is_empty());
                        assert_eq!(error, EMPTY);
                    }
                }
            }
            let close = |x: f64, y: f64| (x - y).abs() < 1e-12;
            let empty_a = (1.0 - q).powi(a.len() as i32);
            assert!(close(zero_policy, truth * (1.0 - empty_a)), "paper policy");
            assert!(
                close(a_sum / a_mass, truth),
                "unbiased given FRAC(A) non-empty"
            );
            let reported = both_sum / both_mass;
            assert!(close(reported, truth * a_mass / both_mass), "C / P(B | A)");
            let empty_b = (1.0 - q).powi(b.len() as i32);
            assert!(reported >= truth - 1e-12);
            assert!(reported <= truth / (1.0 - empty_b) + 1e-12);
            if b == [0, 2] {
                // The worked example's three means.
                assert!(close(zero_policy, 0.375) && close(a_sum / a_mass, 0.5));
                assert!(close(reported, 0.6));
            }
        }
    }

    #[test]
    fn estimates_are_within_four_sigma_of_known_truth() {
        // A = X + Y, B = Y + Z; truth is computed exactly from the k-mer sets.
        // Under ideal uniform hashing |FRAC(A)| ~ Binomial(|A|, 1/S); the
        // containment's spread is approximated by sqrt(C(1 − C) / (|A| / S)).
        let k = 21;
        let scaled = 100;
        let x = random_bases(11, 60_000);
        let y = random_bases(12, 40_000);
        let z = random_bases(13, 100_000);
        let p = params(k, scaled, 2024);
        let a = sketch(&fasta(&[&x, &y]), p.clone());
        let b = sketch(&fasta(&[&y, &z]), p);
        let sa = naive_kmers(&[&x, &y], k);
        let sb = naive_kmers(&[&y, &z], k);
        let inter = sa.intersection(&sb).count() as f64;
        let union = sa.union(&sb).count() as f64;
        let s = scaled as f64;

        for (sketch, truth) in [(&a, sa.len() as f64), (&b, sb.len() as f64)] {
            let estimate = sketch.hashes.len() as f64 * s;
            let sigma = (truth * (s - 1.0)).sqrt();
            assert!(
                (estimate - truth).abs() < 4.0 * sigma,
                "cardinality {estimate} vs {truth}"
            );
        }
        let e = estimate(&a, &b).unwrap();
        let check = |name: &str, value: f64, truth: f64, expected_size: f64| {
            let sigma = (truth * (1.0 - truth) / expected_size).sqrt();
            assert!(
                (value - truth).abs() < 4.0 * sigma,
                "{name} {value} vs {truth} (sigma {sigma})"
            );
        };
        check("J", e.jaccard, inter / union, union / s);
        check(
            "C(A,B)",
            e.containment_a_in_b,
            inter / sa.len() as f64,
            sa.len() as f64 / s,
        );
        check(
            "C(B,A)",
            e.containment_b_in_a,
            inter / sb.len() as f64,
            sb.len() as f64 / s,
        );
    }
}
