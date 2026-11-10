//! Library-level FracMinHash accuracy on a real genome with known mutation rates.
//!
//! The phiX174 genome (read-only, from data/) is mutated in memory by i.i.d.
//! substitutions at rate p with a fixed seed. Under the simple mutation model,
//! containment is about (1 − p)^k, so containment ANI should recover 1 − p.

use std::num::NonZeroUsize;

use aiips::sketch::fracminhash::{estimate, FracMinHashParams, FracMinHashSketch};

const GENOME: &[u8] = include_bytes!("../data/phiX174_NC_001422.1.fasta");

fn genome_sequence() -> Vec<u8> {
    GENOME
        .split(|&b| b == b'\n')
        .filter(|line| !line.starts_with(b">"))
        .flatten()
        .copied()
        .filter(|b| !b.is_ascii_whitespace())
        .collect()
}

/// Substitutes each base with probability `rate` by one of the other three,
/// using a fixed-seed xorshift64* stream. Returns the mutated sequence and the
/// realized number of substitutions.
fn mutate(sequence: &[u8], rate: f64, seed: u64) -> (Vec<u8>, usize) {
    let mut state = seed;
    let mut next = || {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        state.wrapping_mul(0x2545_f491_4f6c_dd1d)
    };
    let mut substitutions = 0;
    let mutated = sequence
        .iter()
        .map(|&base| {
            let draw = (next() >> 11) as f64 / (1u64 << 53) as f64;
            if draw >= rate {
                return base;
            }
            substitutions += 1;
            let others: Vec<u8> = b"ACGT".iter().copied().filter(|&b| b != base).collect();
            others[(next() % 3) as usize]
        })
        .collect();
    (mutated, substitutions)
}

fn sketch(sequence: &[u8], params: &FracMinHashParams) -> FracMinHashSketch {
    let mut fasta = b">s\n".to_vec();
    fasta.extend_from_slice(sequence);
    fasta.push(b'\n');
    FracMinHashSketch::build(fasta.as_slice(), params.clone()).unwrap()
}

fn params(scaled: u64) -> FracMinHashParams {
    FracMinHashParams::new(NonZeroUsize::new(21).unwrap(), scaled, 42).unwrap()
}

#[test]
fn containment_ani_recovers_known_substitution_rates() {
    let genome = genome_sequence();
    assert_eq!(genome.len(), 5386, "phiX174 NC_001422.1 length");
    for scaled in [1, 10] {
        let p = params(scaled);
        let original = sketch(&genome, &p);
        // phiX174 has 5366 distinct canonical 21-mers (scaled=1 keeps all).
        if scaled == 1 {
            assert_eq!(original.hashes.len(), 5366);
        }
        // Under ideal hashing the size is Binomial(5366, 1/scaled).
        let expected = 5366.0 / scaled as f64;
        let size = original.hashes.len() as f64;
        assert!(
            (size - expected).abs() < 4.0 * expected.sqrt() + 1.0,
            "scaled={scaled}: {size} hashes"
        );
        for (rate, seed) in [(0.01, 11), (0.05, 12)] {
            let (mutated, substitutions) = mutate(&genome, rate, seed);
            let realized = substitutions as f64 / genome.len() as f64;
            let e = estimate(&original, &sketch(&mutated, &p)).unwrap();
            for ani in [e.ani_containment_a_in_b, e.ani_containment_b_in_a] {
                assert!(
                    (ani - (1.0 - realized)).abs() < 0.008,
                    "scaled={scaled} p={rate}: containment ANI {ani} vs {}",
                    1.0 - realized
                );
            }
            assert!(
                (e.ani - (1.0 - realized)).abs() < 0.01,
                "scaled={scaled} p={rate}: Mash ANI {}",
                e.ani
            );
            assert!(e.jaccard > 0.0 && e.jaccard < 1.0);
        }
    }
}

#[test]
fn unrelated_random_sequence_shares_nothing() {
    let genome = genome_sequence();
    let p = params(10);
    let mut reversed = genome.clone();
    reversed.reverse(); // Reversal (not reverse complement) destroys shared 21-mers.
    let e = estimate(&sketch(&genome, &p), &sketch(&reversed, &p)).unwrap();
    assert!(e.jaccard < 0.01, "J {}", e.jaccard);
    assert!(e.containment_a_in_b < 0.02, "C {}", e.containment_a_in_b);
    assert!(e.ani < 0.85, "ANI {}", e.ani);
}
