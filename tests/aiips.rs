//! Library-level accuracy of the AIIPS sketch: the global probe set and the
//! positional anchor/probe occurrences.
//!
//! Checks against known baselines: exact distinct canonical k-mer sets
//! (computed independently from strings), the FracMinHash sketch with the
//! same k, scale and seed (which the global probe set must equal bit for bit),
//! expected sampling rates, and the phiX174 genome mutated in memory at known
//! substitution rates. Synthetic sequences come from a fixed-seed generator.

use std::collections::HashSet;
use std::num::NonZeroUsize;

use aiips::sketch::aiips::{
    compare, estimate, AiipsComparison, AiipsParams, AiipsSketch, LocalConfig, LocalStatus,
};
use aiips::sketch::fracminhash::{self, FracMinHashParams, FracMinHashSketch};

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

/// xorshift64* stream with a fixed seed.
fn generator(seed: u64) -> impl FnMut() -> u64 {
    let mut state = seed;
    move || {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        state.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }
}

fn random_bases(seed: u64, len: usize) -> Vec<u8> {
    let mut next = generator(seed);
    (0..len).map(|_| b"ACGT"[(next() >> 62) as usize]).collect()
}

/// Substitutes each base with probability `rate` by one of the other three.
/// Returns the mutated sequence and the realized number of substitutions.
fn mutate(sequence: &[u8], rate: f64, seed: u64) -> (Vec<u8>, usize) {
    let mut next = generator(seed);
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

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

fn fasta(records: &[&[u8]]) -> Vec<u8> {
    let mut out = Vec::new();
    for (index, record) in records.iter().enumerate() {
        out.extend_from_slice(format!(">contig{index}\n").as_bytes());
        out.extend_from_slice(record);
        out.push(b'\n');
    }
    out
}

/// Probe parameters with default anchors (k1 = 21, c1 = 125).
fn params(k: usize, scaled: u64) -> AiipsParams {
    let anchor_k = NonZeroUsize::new(21).unwrap();
    AiipsParams::new(anchor_k, 125, NonZeroUsize::new(k).unwrap(), scaled, 42).unwrap()
}

fn sketch(records: &[&[u8]], p: &AiipsParams) -> AiipsSketch {
    AiipsSketch::build(fasta(records).as_slice(), p.clone()).unwrap()
}

/// Distinct canonical k-mers as strings: an independent exact baseline.
fn kmer_set(records: &[&[u8]], k: usize) -> HashSet<Vec<u8>> {
    let mut set = HashSet::new();
    for record in records {
        for window in record.windows(k) {
            let rc = reverse_complement(window);
            set.insert(window.to_vec().min(rc));
        }
    }
    set
}

/// Exact (Jaccard, containment of A in B, containment of B in A).
fn exact(a: &[&[u8]], b: &[&[u8]], k: usize) -> (f64, f64, f64) {
    let (a, b) = (kmer_set(a, k), kmer_set(b, k));
    let shared = a.intersection(&b).count() as f64;
    (
        shared / a.union(&b).count() as f64,
        shared / a.len() as f64,
        shared / b.len() as f64,
    )
}

/// Asserts `estimate` is within 4 standard errors of the exact ratio `truth`,
/// treating the denominator's sketch as Binomial(n, 1/scaled) draws.
fn within_four_se(name: &str, estimate: f64, truth: f64, denominator: usize, scaled: u64) {
    let sampled = denominator as f64 / scaled as f64;
    let se = (truth * (1.0 - truth) / sampled).sqrt();
    assert!(
        (estimate - truth).abs() <= 4.0 * se + 1e-12,
        "{name}: estimate {estimate} vs exact {truth} (se {se})"
    );
}

#[test]
fn global_estimates_equal_fracminhash_on_mutated_phix() {
    let genome = genome_sequence();
    assert_eq!(genome.len(), 5386, "phiX174 NC_001422.1 length");
    for scaled in [1, 10, 125] {
        let p = params(15, scaled);
        let frac = FracMinHashParams::new(NonZeroUsize::new(15).unwrap(), scaled, 42).unwrap();
        let frac_sketch = |sequence: &[u8]| {
            FracMinHashSketch::build(fasta(&[sequence]).as_slice(), frac.clone()).unwrap()
        };
        let original = sketch(&[&genome], &p);
        assert_eq!(original.global_probes(), frac_sketch(&genome).hashes);
        for (rate, seed) in [(0.01, 11), (0.05, 12)] {
            let (mutated, _) = mutate(&genome, rate, seed);
            let a = estimate(&original, &sketch(&[&mutated], &p)).unwrap();
            let f = fracminhash::estimate(&frac_sketch(&genome), &frac_sketch(&mutated)).unwrap();
            let context = format!("scaled={scaled} p={rate}");
            assert_eq!(a.jaccard, f.jaccard, "{context}");
            assert_eq!(a.containment_a_in_b, f.containment_a_in_b, "{context}");
            assert_eq!(a.containment_b_in_a, f.containment_b_in_a, "{context}");
            assert_eq!(
                a.ani_containment_a_in_b, f.ani_containment_a_in_b,
                "{context}"
            );
            assert_eq!(
                a.ani_containment_b_in_a, f.ani_containment_b_in_a,
                "{context}"
            );
        }
    }
}

#[test]
fn scaled_one_matches_exact_sets_and_recovers_mutation_rate_on_phix() {
    let genome = genome_sequence();
    let p = params(15, 1);
    let original = sketch(&[&genome], &p);
    for (rate, seed) in [(0.01, 11), (0.05, 12)] {
        let (mutated, substitutions) = mutate(&genome, rate, seed);
        let e = estimate(&original, &sketch(&[&mutated], &p)).unwrap();
        let (jaccard, a_in_b, b_in_a) = exact(&[&genome], &[&mutated], 15);
        assert_eq!(
            (e.jaccard, e.containment_a_in_b, e.containment_b_in_a),
            (jaccard, a_in_b, b_in_a),
            "p={rate}"
        );
        // Screening ANI C^(1/k2) under the simple substitution model.
        let identity = 1.0 - substitutions as f64 / genome.len() as f64;
        for ani in [
            e.ani_containment_a_in_b,
            e.ani_containment_b_in_a,
            e.putative_ani,
        ] {
            assert!(
                (ani - identity).abs() < 0.008,
                "p={rate}: C^(1/15) {ani} vs identity {identity}"
            );
        }
    }
}

#[test]
fn subsequence_is_contained_and_reverse_containment_tracks_length() {
    let scaled = 125;
    let p = params(15, scaled);
    let genome = random_bases(7, 200_000);
    let part = &genome[50_000..100_000];
    let (full, sub) = (sketch(&[&genome], &p), sketch(&[part], &p));
    let e = estimate(&sub, &full).unwrap();
    // Every 15-mer of the subsequence occurs in the genome.
    assert_eq!(e.containment_a_in_b, 1.0);
    assert_eq!(e.max_containment, 1.0);
    let (jaccard, _, full_in_sub) = exact(&[part], &[&genome], 15);
    let union = kmer_set(&[&genome], 15).len();
    within_four_se(
        "C(full in sub)",
        e.containment_b_in_a,
        full_in_sub,
        union,
        scaled,
    );
    within_four_se("J", e.jaccard, jaccard, union, scaled);
    assert!(e.containment_b_in_a < 0.4, "{}", e.containment_b_in_a);
}

#[test]
fn fragmented_contigs_stay_contained_without_cross_contig_kmers() {
    let scaled = 125;
    let p = params(15, scaled);
    let genome = random_bases(8, 200_000);
    let contigs: Vec<&[u8]> = genome.chunks(10_000).collect();
    let (full, fragmented) = (sketch(&[&genome], &p), sketch(&contigs, &p));
    let e = estimate(&fragmented, &full).unwrap();
    assert_eq!(e.containment_a_in_b, 1.0);
    // 19 breaks remove 19 · 14 boundary-spanning 15-mers from the fragments.
    let (_, _, full_in_fragments) = exact(&contigs, &[&genome], 15);
    assert!(full_in_fragments < 1.0);
    within_four_se(
        "C(full in fragments)",
        e.containment_b_in_a,
        full_in_fragments,
        kmer_set(&[&genome], 15).len(),
        scaled,
    );
}

#[test]
fn unrelated_sequences_share_almost_nothing() {
    let p = params(15, 125);
    let a = sketch(&[&random_bases(9, 200_000)], &p);
    let b = sketch(&[&random_bases(10, 200_000)], &p);
    assert!(a.global_probes().len() > 1000 && b.global_probes().len() > 1000);
    let e = estimate(&a, &b).unwrap();
    assert!(e.jaccard < 0.01, "J {}", e.jaccard);
    assert!(e.max_containment < 0.01, "C {}", e.max_containment);
}

#[test]
fn reverse_complement_gives_identical_global_sets() {
    let p = params(15, 125);
    let genome = random_bases(12, 200_000);
    let rc = reverse_complement(&genome);
    let (a, b) = (sketch(&[&genome], &p), sketch(&[&rc], &p));
    assert_eq!(a.global_probes(), b.global_probes());
    let e = estimate(&a, &b).unwrap();
    assert_eq!(
        (e.jaccard, e.containment_a_in_b, e.containment_b_in_a),
        (1.0, 1.0, 1.0)
    );
}

#[test]
fn global_sets_use_distinct_hashes_while_occurrences_count_repeats() {
    let p = params(15, 125);
    let genome = random_bases(13, 100_000);
    let once = sketch(&[&genome], &p);
    let twice = sketch(&[&genome, &genome], &p);
    // Same distinct set, twice the occurrences.
    assert_eq!(once.global_probes(), twice.global_probes());
    assert_eq!(twice.probes().len(), 2 * once.probes().len());
    assert_eq!(twice.summary.kmers, 2 * once.summary.kmers);
    let e = estimate(&once, &twice).unwrap();
    assert_eq!(
        (e.jaccard, e.containment_a_in_b, e.containment_b_in_a),
        (1.0, 1.0, 1.0)
    );
    // Under ideal hashing occurrences are Binomial(L − k + 1, 1/c).
    let expected = once.summary.kmers as f64 / 125.0;
    let observed = once.probes().len() as f64;
    assert!(
        (observed - expected).abs() < 4.0 * expected.sqrt(),
        "{observed} occurrences vs {expected} expected"
    );
}

#[test]
fn occurrence_counts_follow_the_sampling_rates() {
    // Long low-repeat genome: about L/c1 anchors and L/c2 probes, nearly all
    // distinct. Under ideal hashing each count is Binomial(windows, 1/c).
    let genome = random_bases(14, 400_000);
    let k = |k| NonZeroUsize::new(k).unwrap();
    let p = AiipsParams::new(k(21), 125, k(15), 100, 42).unwrap();
    let sketch = sketch(&[&genome], &p);
    for (name, count, windows, scaled) in [
        ("anchors", sketch.anchors().len(), genome.len() - 20, 125.0),
        ("probes", sketch.probes().len(), genome.len() - 14, 100.0),
    ] {
        let expected = windows as f64 / scaled;
        let sd = (expected * (1.0 - 1.0 / scaled)).sqrt();
        assert!(
            (count as f64 - expected).abs() < 4.0 * sd,
            "{name}: {count} vs {expected}"
        );
    }
    // Chance repeats are rare among 21-mers and 15-mers of 400 kb: about
    // n^2 / 4^k colliding window pairs, so a few sampled 15-mer repeats.
    assert!(sketch.anchors().unique() as f64 > 0.99 * sketch.anchors().len() as f64);
    let (distinct, occurrences) = (sketch.global_probes().len(), sketch.probes().len());
    assert!(distinct <= occurrences && distinct as f64 > 0.99 * occurrences as f64);
}

#[test]
fn phix_islands_partition_every_probe() {
    let genome = genome_sequence();
    let k = |k| NonZeroUsize::new(k).unwrap();
    let p = AiipsParams::new(k(21), 50, k(15), 10, 42).unwrap();
    let sketch = sketch(&[&genome], &p);
    let anchors = sketch.contigs().anchors(0);
    assert!(anchors.len() > 50, "{} anchors", anchors.len());
    // START plus every island, in order, covers each probe exactly once.
    let mut covered = sketch.start_bucket(0).len();
    for anchor in anchors.clone() {
        covered += sketch.owned_probes(anchor).len();
    }
    assert_eq!(covered, sketch.probes().len());
    // The slice from the first to the last anchor equals the probes between
    // their positions, across every inner island.
    let (first, last) = (anchors.start, anchors.end - 1);
    let slice = sketch.probes_between(first, last).unwrap();
    let positions = &sketch.anchors().positions;
    let inside = sketch
        .probes()
        .positions
        .iter()
        .filter(|&&q| q >= positions[first] && q < positions[last])
        .count();
    assert_eq!(slice.len(), inside);
}

// ---------------------------------------------------------------------------
// Positional ANI and AF (shared anchors, validated blocks, scoring probes).

/// Default AIIPS parameters: k1 = 21, c1 = 125, k2 = 15, c2 = 125.
fn defaults() -> AiipsParams {
    let k = |k| NonZeroUsize::new(k).unwrap();
    AiipsParams::new(k(21), 125, k(15), 125, 42).unwrap()
}

fn positional(a: &[&[u8]], b: &[&[u8]], p: &AiipsParams) -> AiipsComparison {
    compare(&sketch(a, p), &sketch(b, p), &LocalConfig::default()).unwrap()
}

/// Asserts the positional ANI is within 4 model standard errors of `identity`,
/// plus 0.001 for the residual bias seen at 90% ANI in simulations. The
/// standard error treats q = 2X/(N_A + N_B) as a binomial proportion over
/// (N_A + N_B)/2 probes and propagates it through q^(1/k2).
fn assert_ani_near(c: &AiipsComparison, identity: f64, context: &str) {
    let l = &c.local;
    assert_eq!(l.ani_status, LocalStatus::Estimated, "{context}");
    let (ani, q) = (l.ani.unwrap(), l.q.unwrap());
    let n = (l.counts.scoring_a + l.counts.scoring_b) as f64 / 2.0;
    let se = (q * (1.0 - q) / n).sqrt() * ani / (15.0 * q);
    assert!(
        (ani - identity).abs() < 4.0 * se + 0.001,
        "{context}: ANI {ani} vs identity {identity} (se {se})"
    );
}

fn af(c: &AiipsComparison) -> (f64, f64) {
    (c.local.af_a.unwrap(), c.local.af_b.unwrap())
}

#[test]
fn positional_ani_tracks_known_substitution_rates() {
    let genome = random_bases(40, 600_000);
    let p = defaults();
    let a = sketch(&[&genome], &p);
    for (rate, seed) in [(0.001, 41), (0.01, 42), (0.05, 43), (0.10, 44)] {
        let (mutated, substitutions) = mutate(&genome, rate, seed);
        let identity = 1.0 - substitutions as f64 / genome.len() as f64;
        let c = compare(&a, &sketch(&[&mutated], &p), &LocalConfig::default()).unwrap();
        assert_ani_near(&c, identity, &format!("p={rate}"));
        // The whole genome is homologous; AF falls short only at the ends
        // (outside the outermost shared anchors) and, at 10% divergence where
        // shared anchors are about 1.1 kb apart, where a run breaks.
        let least = if rate <= 0.05 { 0.99 } else { 0.95 };
        let (af_a, af_b) = af(&c);
        assert!(af_a > least && af_b > least, "p={rate}: AF {af_a} {af_b}");
        assert_eq!(c.local.af_status, LocalStatus::Estimated);
    }
}

#[test]
fn divergent_genomes_are_estimated_or_left_unresolved() {
    // About 80% identity: many blocks lack exact mapping hits. ANI is either
    // reported close to the truth or null, never a confident wrong value.
    let genome = random_bases(45, 600_000);
    let (mutated, substitutions) = mutate(&genome, 0.2, 46);
    let identity = 1.0 - substitutions as f64 / genome.len() as f64;
    let c = positional(&[&genome], &[&mutated], &defaults());
    assert!(c.local.counts.rejected_few_hits > 0);
    match c.local.ani {
        Some(ani) => assert!((ani - identity).abs() < 0.02, "ANI {ani} vs {identity}"),
        None => assert_ne!(c.local.ani_status, LocalStatus::Estimated),
    }
    let (af_a, _) = af(&c);
    assert!(af_a > 0.0 && af_a < 1.0);
}

#[test]
fn phix_positional_ani_recovers_known_rates() {
    // phiX174 is short, so sample densely (c1 = c2 = 10) for enough probes.
    let genome = genome_sequence();
    let k = |k| NonZeroUsize::new(k).unwrap();
    let p = AiipsParams::new(k(21), 10, k(15), 10, 42).unwrap();
    for (rate, seed) in [(0.01, 11), (0.05, 12)] {
        let (mutated, substitutions) = mutate(&genome, rate, seed);
        let identity = 1.0 - substitutions as f64 / genome.len() as f64;
        let c = positional(&[&genome], &[&mutated], &p);
        assert_ani_near(&c, identity, &format!("phiX p={rate}"));
        // One block; the ends beyond the outermost shared anchors (about 100
        // bases at 5% divergence) are a few percent of this short genome.
        let (af_a, af_b) = af(&c);
        assert!(
            af_a > 0.95 && af_b > 0.95,
            "phiX p={rate}: AF {af_a} {af_b}"
        );
    }
}

#[test]
fn subsequence_af_reflects_the_shared_fraction() {
    let genome = random_bases(47, 600_000);
    let part = &genome[150_000..450_000];
    let c = positional(&[&genome], &[part], &defaults());
    let (af_full, af_sub) = af(&c);
    assert!(af_sub > 0.99, "AF of the subsequence {af_sub}");
    assert!(
        (af_full - 0.5).abs() < 0.01,
        "AF of the full genome {af_full}"
    );
    assert!(af_sub > af_full);
    assert_eq!(c.local.ani, Some(1.0));
}

#[test]
fn deletions_and_insertions_are_left_out_of_af() {
    let genome = random_bases(48, 600_000);
    // 30% of A deleted from B.
    let deleted = [&genome[..200_000], &genome[380_000..]].concat();
    let c = positional(&[&genome], &[&deleted], &defaults());
    let (af_a, af_b) = af(&c);
    assert!((af_a - 0.7).abs() < 0.01 && af_b > 0.99, "AF {af_a} {af_b}");
    assert_eq!(c.local.ani, Some(1.0));
    // 30 kb of unrelated sequence inserted into B: within the gap limit, but
    // the gap ratio breaks the run around it.
    let inserted = [
        &genome[..300_000],
        &random_bases(49, 30_000),
        &genome[300_000..],
    ]
    .concat();
    let c = positional(&[&genome], &[&inserted], &defaults());
    let (af_a, af_b) = af(&c);
    assert!(
        af_a > 0.99 && (af_b - 600.0 / 630.0).abs() < 0.01,
        "AF {af_a} {af_b}"
    );
    assert!(c.local.counts.matching.breaks_ratio >= 1);
    assert_eq!(c.local.ani, Some(1.0));
}

#[test]
fn inversions_translocations_and_fragments_stay_covered() {
    let genome = random_bases(50, 600_000);
    let inverted = [
        &genome[..200_000],
        &reverse_complement(&genome[200_000..300_000]),
        &genome[300_000..],
    ]
    .concat();
    let translocated = [&genome[300_000..], &genome[..300_000]].concat();
    let fragments: Vec<&[u8]> = genome.chunks(20_000).collect();
    for (name, b) in [
        ("inversion", vec![inverted.as_slice()]),
        ("translocation", vec![translocated.as_slice()]),
        ("fragments", fragments),
    ] {
        let c = positional(&[&genome], &b, &defaults());
        let (af_a, af_b) = af(&c);
        assert!(af_a > 0.97 && af_b > 0.97, "{name}: AF {af_a} {af_b}");
        assert_eq!(c.local.ani, Some(1.0), "{name}");
    }
}

#[test]
fn duplicated_regions_are_not_mapped() {
    // Anchors in the duplicated 100 kb occur twice in B, so they are ignored:
    // the region stays uncovered rather than being assigned to one copy.
    let genome = random_bases(51, 600_000);
    let duplicated = [genome.as_slice(), &genome[100_000..200_000]].concat();
    let c = positional(&[&genome], &[&duplicated], &defaults());
    let (af_a, af_b) = af(&c);
    assert!((af_a - 5.0 / 6.0).abs() < 0.01, "AF of A {af_a}");
    assert!((af_b - 5.0 / 7.0).abs() < 0.01, "AF of B {af_b}");
    assert!(c.local.counts.matching.repeated_anchors > 500);
    assert_eq!(c.local.ani, Some(1.0));
}

#[test]
fn unrelated_genomes_report_no_supported_mapping() {
    let c = positional(
        &[&random_bases(52, 300_000)],
        &[&random_bases(53, 300_000)],
        &defaults(),
    );
    assert_eq!(
        (c.local.ani, c.local.af_a, c.local.af_b),
        (None, None, None)
    );
    assert_eq!(c.local.ani_status, LocalStatus::NoSupportedMapping);
    assert_eq!(c.local.af_status, LocalStatus::NoSupportedMapping);
    assert_eq!((c.af_containment_a, c.af_containment_b), (None, None));
}

#[test]
fn af_containment_cross_check_is_near_one_for_full_homologs() {
    // C / q estimates the homologous k-mer fraction: about 1 when the whole
    // genome is homologous, about 0.5 for the full genome against a half.
    let genome = random_bases(54, 600_000);
    let (mutated, _) = mutate(&genome, 0.02, 55);
    let c = positional(&[&genome], &[&mutated], &defaults());
    for proxy in [c.af_containment_a, c.af_containment_b] {
        assert!((proxy.unwrap() - 1.0).abs() < 0.03, "{proxy:?}");
    }
    let c = positional(&[&genome], &[&genome[..300_000]], &defaults());
    assert!((c.af_containment_a.unwrap() - 0.5).abs() < 0.02);
    assert!((c.af_containment_b.unwrap() - 1.0).abs() < 0.02);
}
