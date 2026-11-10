//! Classic MinHash (Broder 1997): `s` hash functions, keeping each one's
//! minimum over the distinct canonical k-mers of the input.
//!
//! Hash function i is h_i(x) = fmix64(x XOR seed_i), with seed_i the i-th
//! SplitMix64 output from the user seed. This approximates, but does not
//! guarantee, min-wise independence.

use std::io::Read;
use std::num::NonZeroUsize;

use super::{invalid, parse_u64, Fields, ParseError};
use crate::fastx::{summarize_with, Summary};
use crate::hash::{fmix64, SplitMix64};
use crate::kmer::{CanonicalKmers, MAX_PACKED_K};

pub const HASH: &str = "murmur3-fmix64";
pub const MIN_NUM_HASHES: usize = 2;
pub const DEFAULT_K: usize = 21;
pub const DEFAULT_NUM_HASHES: usize = 1000;
pub const DEFAULT_SEED: u64 = 42;
pub const MAX_NUM_HASHES: usize = 1_000_000;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MinHashParams {
    pub k: NonZeroUsize,
    pub num_hashes: usize,
    pub seed: u64,
    /// Hash family identifier; sketches compare only within one family.
    pub hash: String,
    pub canonical: bool,
}

impl MinHashParams {
    /// Parameters for building a sketch with this crate's hash family.
    pub fn new(k: NonZeroUsize, num_hashes: usize, seed: u64) -> Result<Self, String> {
        check_k(k.get()).map_err(|e| e.to_string())?;
        check_num_hashes(num_hashes).map_err(|e| e.to_string())?;
        Ok(Self {
            k,
            num_hashes,
            seed,
            hash: HASH.into(),
            canonical: true,
        })
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
        check(
            "num_hashes",
            self.num_hashes.to_string(),
            other.num_hashes.to_string(),
        );
        check("seed", self.seed.to_string(), other.seed.to_string());
        check("hash", self.hash.clone(), other.hash.clone());
        check(
            "canonical",
            self.canonical.to_string(),
            other.canonical.to_string(),
        );
        differences
    }
}

fn check_k(k: usize) -> Result<(), ParseError> {
    if (1..=MAX_PACKED_K).contains(&k) {
        Ok(())
    } else {
        Err(invalid(format!(
            "minhash supports k-mer lengths 1 to {MAX_PACKED_K}, got {k}"
        )))
    }
}

fn check_num_hashes(num_hashes: usize) -> Result<(), ParseError> {
    if (MIN_NUM_HASHES..=MAX_NUM_HASHES).contains(&num_hashes) {
        Ok(())
    } else {
        Err(invalid(format!(
            "number of hash functions must be from {MIN_NUM_HASHES} to {MAX_NUM_HASHES}, got {num_hashes}"
        )))
    }
}

/// Running minima for one input.
pub struct MinHasher {
    seeds: Vec<u64>,
    mins: Vec<u64>,
}

impl MinHasher {
    pub fn new(num_hashes: usize, seed: u64) -> Self {
        Self {
            seeds: SplitMix64::new(seed).take(num_hashes).collect(),
            mins: vec![u64::MAX; num_hashes],
        }
    }

    #[inline]
    pub fn insert(&mut self, kmer: u64) {
        for (min, seed) in self.mins.iter_mut().zip(&self.seeds) {
            *min = (*min).min(fmix64(kmer ^ seed));
        }
    }

    pub fn into_mins(self) -> Vec<u64> {
        self.mins
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MinHashSketch {
    pub params: MinHashParams,
    pub summary: Summary,
    pub mins: Vec<u64>,
}

/// Why a sketch could not be built.
#[derive(Debug)]
pub enum BuildError {
    /// Parameters this builder cannot honor; building would mislabel the data.
    Unsupported(String),
    Input(needletail::errors::ParseError),
}

impl std::fmt::Display for BuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BuildError::Unsupported(message) => write!(f, "{message}"),
            BuildError::Input(e) => write!(f, "{e}"),
        }
    }
}

impl MinHashParams {
    /// Checks that the builder can produce exactly what these parameters
    /// describe. Fields are public (parsed sketches carry arbitrary values),
    /// so construction must not trust them.
    fn check_buildable(&self) -> Result<(), BuildError> {
        let unsupported = |e: ParseError| match e {
            ParseError::Invalid(message) => BuildError::Unsupported(message),
            other => BuildError::Unsupported(other.to_string()),
        };
        check_k(self.k.get()).map_err(unsupported)?;
        check_num_hashes(self.num_hashes).map_err(unsupported)?;
        if self.hash != HASH {
            return Err(BuildError::Unsupported(format!(
                "minhash can only build hash={HASH}, not hash={}",
                self.hash
            )));
        }
        if !self.canonical {
            return Err(BuildError::Unsupported(
                "minhash can only build canonical=true sketches".into(),
            ));
        }
        Ok(())
    }
}

impl MinHashSketch {
    /// Builds the sketch and summary in one streaming pass over FASTA/FASTQ.
    /// Unsupported parameters are refused before any input is read.
    pub fn build<R: Read + Send>(input: R, params: MinHashParams) -> Result<Self, BuildError> {
        params.check_buildable()?;
        let mut hasher = MinHasher::new(params.num_hashes, params.seed);
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
            mins: hasher.into_mins(),
        })
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        use std::fmt::Write;
        let p = &self.params;
        let mut out = String::with_capacity(160 + 21 * self.mins.len());
        let _ = writeln!(
            out,
            "{} {} algorithm=minhash k={} num_hashes={} seed={} hash={} canonical={}",
            super::MAGIC,
            super::VERSION,
            p.k,
            p.num_hashes,
            p.seed,
            p.hash,
            p.canonical
        );
        super::write_summary(&mut out, &self.summary);
        for min in &self.mins {
            let _ = writeln!(out, "{min}");
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
                    "minhash supports k-mer lengths 1 to {MAX_PACKED_K}, got {raw_k}"
                ))
            })?;
        let num_hashes = header.take_u64("num_hashes")?;
        let num_hashes = usize::try_from(num_hashes).unwrap_or(usize::MAX);
        check_num_hashes(num_hashes)?;
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
        header.finish()?;
        let summary = super::parse_summary(lines.next())?;
        let mut mins = Vec::with_capacity(num_hashes);
        for (index, line) in lines.enumerate() {
            if index == num_hashes {
                return Err(invalid(format!(
                    "more than num_hashes={num_hashes} hash values"
                )));
            }
            mins.push(parse_u64(line).ok_or_else(|| {
                invalid(format!(
                    "hash value {} is not an unsigned integer: {line:?}",
                    index + 1
                ))
            })?);
        }
        if mins.len() != num_hashes {
            return Err(invalid(format!(
                "expected num_hashes={num_hashes} hash values, found {}",
                mins.len()
            )));
        }
        Ok(Self {
            params: MinHashParams {
                k,
                num_hashes,
                seed,
                hash,
                canonical,
            },
            summary,
            mins,
        })
    }
}

/// Point estimates for one sketch pair. AF is not estimable from a global
/// MinHash (positions are discarded), so it is not part of this result.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Estimates {
    pub jaccard: f64,
    pub containment_a_in_b: f64,
    pub containment_b_in_a: f64,
    pub ani: f64,
    pub ani_containment_a_in_b: f64,
    pub ani_containment_b_in_a: f64,
}

/// Compares two sketches built with identical parameters.
pub fn estimate(a: &MinHashSketch, b: &MinHashSketch) -> Result<Estimates, String> {
    let differences = a.params.mismatches(&b.params);
    if !differences.is_empty() {
        return Err(format!(
            "sketch parameters differ: {}",
            differences.join(", ")
        ));
    }
    let k = a.params.k.get();
    let jaccard = jaccard(&a.mins, &b.mins);
    let (containment_a_in_b, containment_b_in_a) =
        containments(jaccard, cardinality(&a.mins), cardinality(&b.mins));
    Ok(Estimates {
        jaccard,
        containment_a_in_b,
        containment_b_in_a,
        ani: mash_ani(jaccard, k),
        ani_containment_a_in_b: containment_ani(containment_a_in_b, k),
        ani_containment_b_in_a: containment_ani(containment_b_in_a, k),
    })
}

/// Fraction of hash functions whose minima agree (Broder 1997).
pub fn jaccard(a: &[u64], b: &[u64]) -> f64 {
    assert_eq!(a.len(), b.len(), "sketch sizes must match");
    let matches = a.iter().zip(b).filter(|(x, y)| x == y).count();
    matches as f64 / a.len() as f64
}

/// Distinct-element estimate n̂ = (s − 1) / Σ −ln(1 − u_i), with u_i the
/// minima scaled to [0, 1). Each −ln(1 − u_i) is Exp(n), so the estimator is
/// unbiased (Cohen 1997). Requires s >= 2.
pub fn cardinality(mins: &[u64]) -> f64 {
    assert!(mins.len() >= MIN_NUM_HASHES, "need at least two minima");
    const SCALE: f64 = 1.0 / 18_446_744_073_709_551_616.0; // 2^-64
    let sum: f64 = mins.iter().map(|&m| -(-(m as f64) * SCALE).ln_1p()).sum();
    (mins.len() - 1) as f64 / sum
}

/// Containment in both directions from |A ∩ B| = J / (1 + J) · (|A| + |B|)
/// (Koslicki & Zabeti 2019), clamped to [0, 1] against estimation noise.
pub fn containments(jaccard: f64, size_a: f64, size_b: f64) -> (f64, f64) {
    let intersection = jaccard / (1.0 + jaccard) * (size_a + size_b);
    let ratio = |size: f64| {
        if size > 0.0 && size.is_finite() {
            (intersection / size).clamp(0.0, 1.0)
        } else {
            0.0
        }
    };
    (ratio(size_a), ratio(size_b))
}

/// 1 − D with Mash distance D = −(1/k) ln(2J / (1 + J)) and D = 1 at J = 0
/// (Ondov et al. 2016), clamped at 0 because 1 − D is negative for tiny J.
pub fn mash_ani(jaccard: f64, k: usize) -> f64 {
    if jaccard <= 0.0 {
        return 0.0;
    }
    let distance = -(2.0 * jaccard / (1.0 + jaccard)).ln() / k as f64;
    (1.0 - distance).max(0.0)
}

/// Containment ANI point estimate C^(1/k) under the simple mutation model
/// (Blanca et al. 2022; Hera et al. 2023); 0 when C = 0.
pub fn containment_ani(containment: f64, k: usize) -> f64 {
    if containment <= 0.0 {
        0.0
    } else {
        containment.powf(1.0 / k as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    fn params(k: usize, num_hashes: usize, seed: u64) -> MinHashParams {
        MinHashParams::new(NonZeroUsize::new(k).unwrap(), num_hashes, seed).unwrap()
    }

    fn sketch(fasta: &[u8], p: MinHashParams) -> MinHashSketch {
        MinHashSketch::build(fasta, p).unwrap()
    }

    /// Independent baseline: distinct canonical k-mers as strings, encoded from
    /// scratch, then each hash function's minimum by brute force.
    fn naive_mins(sequences: &[&[u8]], k: usize, num_hashes: usize, seed: u64) -> Vec<u64> {
        let complement = |b: u8| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        };
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
        SplitMix64::new(seed)
            .take(num_hashes)
            .map(|s| distinct.iter().map(|&x| fmix64(x ^ s)).min().unwrap())
            .collect()
    }

    #[test]
    fn streamed_sketch_equals_naive_baseline() {
        let records: [&[u8]; 3] = [
            b"ACGTTGCAnnACGgattacaACGT",
            b"TTTTCCCCAAAAGGGG",
            b"acgtNACGT",
        ];
        let fasta = b">a\nACGTTGCAnnACG\ngattacaACGT\n>b\nTTTTCCCCAAAAGGGG\n>c\nacgtNACGT\n";
        for (k, seed) in [(1, 0), (3, 42), (5, 7), (11, u64::MAX)] {
            let built = sketch(fasta, params(k, 64, seed));
            assert_eq!(built.mins, naive_mins(&records, k, 64, seed), "k={k}");
        }
    }

    #[test]
    fn identical_disjoint_and_reverse_complement_cases_are_exact() {
        let p = params(5, 200, 42);
        let a = sketch(b">a\nGATTACAGGCTTAACGTTAG\n", p.clone());
        let estimates = estimate(&a, &a).unwrap();
        assert_eq!(estimates.jaccard, 1.0);
        assert_eq!(estimates.containment_a_in_b, 1.0);
        assert_eq!(estimates.containment_b_in_a, 1.0);
        assert_eq!(estimates.ani, 1.0);
        assert_eq!(estimates.ani_containment_a_in_b, 1.0);

        let rc = sketch(b">rc\nCTAACGTTAAGCCTGTAATC\n", p.clone());
        assert_eq!(rc.mins, a.mins, "canonical k-mers ignore strand");

        let poly_a = sketch(b">a\nAAAAAAAAAA\n", p.clone());
        let poly_c = sketch(b">c\nCCCCCCCCCC\n", p);
        let disjoint = estimate(&poly_a, &poly_c).unwrap();
        assert_eq!(disjoint.jaccard, 0.0);
        assert_eq!(disjoint.containment_a_in_b, 0.0);
        assert_eq!(disjoint.ani, 0.0);
        assert_eq!(disjoint.ani_containment_b_in_a, 0.0);
    }

    #[test]
    fn formulas_match_hand_computed_values() {
        let d = -(2.0f64 / 3.0).ln() / 21.0;
        assert!((mash_ani(0.5, 21) - (1.0 - d)).abs() < 1e-15);
        assert_eq!(mash_ani(1.0, 21), 1.0);
        assert_eq!(mash_ani(0.0, 21), 0.0);
        assert_eq!(mash_ani(1e-300, 21), 0.0, "clamped at zero");
        assert!((containment_ani(0.5, 21) - 0.5f64.powf(1.0 / 21.0)).abs() < 1e-15);
        assert_eq!(containment_ani(0.0, 21), 0.0);
        // |A| = 100, |B| = 300, |A ∩ B| = 50: J = 50 / 350.
        let (ab, ba) = containments(50.0 / 350.0, 100.0, 300.0);
        assert!((ab - 0.5).abs() < 1e-12 && (ba - 50.0 / 300.0).abs() < 1e-12);
        // Inconsistent noisy inputs: |A ∩ B| = 0.5 · 150 = 75 > |B| = 50.
        assert_eq!(
            containments(1.0, 100.0, 50.0),
            (0.75, 1.0),
            "clamped at one"
        );
        assert_eq!(containments(0.5, 0.0, 10.0).0, 0.0);
    }

    #[test]
    fn mismatched_parameters_are_refused_by_name() {
        let base = sketch(b">a\nACGTACGTAC\n", params(3, 16, 1));
        for (other, field) in [
            (sketch(b">a\nACGTACGTAC\n", params(4, 16, 1)), "k (3 vs 4)"),
            (sketch(b">a\nACGTACGTAC\n", params(3, 17, 1)), "num_hashes"),
            (sketch(b">a\nACGTACGTAC\n", params(3, 16, 2)), "seed"),
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
        let fasta = b">a\nAAAAAA\n>b\nTTTTTT\n";
        let valid = params(3, 4, 1);
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
        p.num_hashes = 1;
        cases.push((p, "from 2"));
        for (p, expected) in cases {
            match MinHashSketch::build(fasta.as_slice(), p) {
                Err(BuildError::Unsupported(message)) => {
                    assert!(message.contains(expected), "{message}")
                }
                other => panic!("expected refusal for {expected}: {other:?}"),
            }
        }
        // Input errors stay distinguishable from parameter errors.
        assert!(matches!(
            MinHashSketch::build(b"".as_slice(), valid),
            Err(BuildError::Input(_))
        ));
    }

    #[test]
    fn parameters_are_validated() {
        let k = |k| NonZeroUsize::new(k).unwrap();
        assert!(MinHashParams::new(k(32), 2, 0).is_ok());
        assert!(MinHashParams::new(k(33), 2, 0)
            .unwrap_err()
            .contains("1 to 32"));
        assert!(MinHashParams::new(k(21), 1, 0).is_err());
        assert!(MinHashParams::new(k(21), MAX_NUM_HASHES + 1, 0).is_err());
        assert!(MinHashParams::new(k(21), MAX_NUM_HASHES, 0).is_ok());
    }

    #[test]
    fn file_round_trip_and_strict_body() {
        let original = sketch(b">a\nGATTACAGGCTTAACGTTAG\n", params(5, 8, 9));
        let bytes = original.to_bytes();
        let text = String::from_utf8(bytes.clone()).unwrap();
        assert!(text.starts_with(
            "aiips-sketch v1 algorithm=minhash k=5 num_hashes=8 seed=9 hash=murmur3-fmix64 canonical=true\nrecords=1 bases=20 non_acgt_bases=0 valid_kmers=16\n"
        ));
        match super::super::parse(&bytes).unwrap() {
            super::super::Sketch::MinHash(parsed) => assert_eq!(parsed, original),
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
        let cases = [
            (with(header, summary, &values[..7]), "found 7"),
            (
                with(header, summary, &[&values[..], &["1"]].concat()),
                "more than",
            ),
            (
                with(header, summary, &[&values[..7], &["x"]].concat()),
                "hash value 8",
            ),
            (
                with(header, summary, &[&values[..7], &["-1"]].concat()),
                "hash value 8",
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
                with(
                    &header.replace("num_hashes=8", "num_hashes=1"),
                    summary,
                    &values[..1],
                ),
                "from 2",
            ),
            (
                with(&header.replace("seed=9", "seed=+9"), summary, &values),
                "seed",
            ),
            (
                with(&header.replace(" canonical=true", ""), summary, &values),
                "canonical",
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
                with(&format!("{header} extra=1"), summary, &values),
                "unknown header field: extra",
            ),
            (
                with(header, "records=1 bases=20 non_acgt_bases=0", &values),
                "valid_kmers",
            ),
            (
                with(
                    header,
                    &summary.replace("valid_kmers=16", "valid_kmers=0"),
                    &values,
                ),
                "no valid k-mers",
            ),
            (
                with(header, &format!("{summary} x=1"), &values),
                "unknown summary field",
            ),
        ];
        for (error, expected) in cases {
            assert!(
                error.contains(expected),
                "expected {expected:?} in {error:?}"
            );
        }
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

    #[test]
    fn estimates_are_within_four_sigma_of_known_truth() {
        // A = X + Y, B = Y + Z with long random segments, so distinct k-mer sets
        // are close to the segment k-mer counts. Truth is computed exactly from
        // the k-mer sets, not from segment lengths.
        let k = 21;
        let s = 2000;
        let x = random_bases(1, 30_000);
        let y = random_bases(2, 20_000);
        let z = random_bases(3, 50_000);
        let fasta = |parts: &[&[u8]]| {
            let mut out = b">r\n".to_vec();
            for p in parts {
                out.extend_from_slice(p);
                out.extend_from_slice(b"N");
            }
            out.push(b'\n');
            out
        };
        let a_fa = fasta(&[&x, &y]);
        let b_fa = fasta(&[&y, &z]);
        let p = params(k, s, 2024);
        let a = sketch(&a_fa, p.clone());
        let b = sketch(&b_fa, p);

        let set = |parts: &[&[u8]]| -> HashSet<Vec<u8>> {
            let mut set = HashSet::new();
            for part in parts {
                for w in part.windows(k) {
                    let rc: Vec<u8> = w
                        .iter()
                        .rev()
                        .map(|b| match b {
                            b'A' => b'T',
                            b'C' => b'G',
                            b'G' => b'C',
                            _ => b'A',
                        })
                        .collect();
                    set.insert(w.to_vec().min(rc));
                }
            }
            set
        };
        let sa = set(&[&x, &y]);
        let sb = set(&[&y, &z]);
        let inter = sa.intersection(&sb).count() as f64;
        let union = sa.union(&sb).count() as f64;
        let true_j = inter / union;

        let e = estimate(&a, &b).unwrap();
        let sigma_j = (true_j * (1.0 - true_j) / s as f64).sqrt();
        assert!(
            (e.jaccard - true_j).abs() < 4.0 * sigma_j,
            "J {} vs {true_j}",
            e.jaccard
        );

        for (sketch, truth) in [(&a, sa.len() as f64), (&b, sb.len() as f64)] {
            let relative = (cardinality(&sketch.mins) - truth).abs() / truth;
            assert!(
                relative < 4.0 / ((s - 2) as f64).sqrt(),
                "cardinality {relative}"
            );
        }
        let true_c_ab = inter / sa.len() as f64;
        let true_c_ba = inter / sb.len() as f64;
        // Containment combines J and two cardinality errors; allow 4σ of each.
        let tolerance = 4.0 * (sigma_j / true_j + 2.0 / ((s - 2) as f64).sqrt());
        assert!((e.containment_a_in_b / true_c_ab - 1.0).abs() < tolerance);
        assert!((e.containment_b_in_a / true_c_ba - 1.0).abs() < tolerance);
    }
}
