use super::*;
use crate::hash::SplitMix64;

fn params(k1: usize, s1: u64, k2: usize, s2: u64) -> AiipsParams {
    AiipsParams::new(
        NonZeroUsize::new(k1).unwrap(),
        s1,
        NonZeroUsize::new(k2).unwrap(),
        s2,
        DEFAULT_SEED,
    )
    .unwrap()
}

fn default_params() -> AiipsParams {
    params(
        DEFAULT_ANCHOR_K,
        DEFAULT_ANCHOR_SCALED,
        DEFAULT_PROBE_K,
        DEFAULT_PROBE_SCALED,
    )
}

fn random_sequence(length: usize, seed: u64) -> Vec<u8> {
    let mut rng = SplitMix64::new(seed);
    let mut sequence = Vec::with_capacity(length);
    while sequence.len() < length {
        let mut word = rng.next().unwrap();
        for _ in 0..32 {
            if sequence.len() == length {
                break;
            }
            sequence.push(b"ACGT"[(word & 3) as usize]);
            word >>= 2;
        }
    }
    sequence
}

/// Substitutes each base independently with probability `rate` (by a
/// different base); returns the mutated copy and the substitution count.
fn substitute(sequence: &[u8], rate: f64, seed: u64) -> (Vec<u8>, usize) {
    let mut rng = SplitMix64::new(seed);
    let mut changed = 0;
    let mutated = sequence
        .iter()
        .map(|&base| {
            let draw = rng.next().unwrap();
            if (draw >> 11) as f64 / (1u64 << 53) as f64 >= rate {
                return base;
            }
            changed += 1;
            let others: Vec<u8> = b"ACGT".iter().copied().filter(|&b| b != base).collect();
            others[(rng.next().unwrap() % 3) as usize]
        })
        .collect();
    (mutated, changed)
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

fn sketch(records: &[&[u8]], p: &AiipsParams) -> AiipsSketch {
    AiipsSketch::from_records(records.iter().copied(), p.clone()).unwrap()
}

fn compare(a: &[&[u8]], b: &[&[u8]]) -> AiipsEstimates {
    let p = default_params();
    estimate(&sketch(a, &p), &sketch(b, &p), WalkOptions::default())
}

/// Independent baseline: hashes every window from scratch, then applies the
/// exclusion and promotion rules literally.
fn brute_force(records: &[&[u8]], p: &AiipsParams) -> AiipsSketch {
    let (seed1, seed2) = p.hash_seeds();
    let encode = |window: &[u8]| -> (u64, bool) {
        let code = |b: u8| {
            b"ACGT"
                .iter()
                .position(|&x| x == b.to_ascii_uppercase())
                .unwrap() as u64
        };
        let forward = window.iter().fold(0, |c, &b| (c << 2) | code(b));
        let reverse = window
            .iter()
            .rev()
            .fold(0, |c, &b| (c << 2) | (3 - code(b)));
        (forward.min(reverse), reverse < forward)
    };
    let sample = |sequence: &[u8], k: usize, seed: u64, max: u64| -> Vec<Sampled> {
        sequence
            .windows(k)
            .enumerate()
            .filter(|(_, w)| w.iter().all(|b| b"ACGTacgt".contains(b)))
            .filter_map(|(pos, w)| {
                let (canonical, reverse) = encode(w);
                let hash = fmix64(canonical ^ seed);
                (hash <= max).then_some(Sampled {
                    pos: pos as u32,
                    reverse,
                    hash,
                })
            })
            .collect()
    };
    let (k1, k2) = (p.anchor_k.get() as u32, p.probe_k.get() as u32);
    let mut out = AiipsSketch {
        params: p.clone(),
        summary: Summary::default(),
        contigs: vec![],
        anchors: vec![],
        anchor_contig: vec![],
        island_start: vec![],
        probes: vec![],
    };
    for sequence in records {
        let mut anchors = sample(sequence, k1 as usize, seed1, p.anchor_max_hash());
        let probes: Vec<Sampled> = sample(sequence, k2 as usize, seed2, p.probe_max_hash())
            .into_iter()
            .filter(|q| {
                !anchors
                    .iter()
                    .any(|a| a.pos < q.pos + k2 && q.pos < a.pos + k1)
            })
            .collect();
        let first_anchor = anchors.first().map_or(u32::MAX, |a| a.pos);
        let promoted = probes.first().is_some_and(|q| q.pos < first_anchor);
        // The promoted probe is the first anchor and stays a probe.
        if promoted {
            anchors.insert(0, probes[0]);
        }
        let first = out.anchors.len() as u32;
        for (i, anchor) in anchors.iter().enumerate() {
            out.island_start.push(out.probes.len() as u32);
            out.anchors.push(*anchor);
            out.anchor_contig.push(out.contigs.len() as u32);
            let next = anchors.get(i + 1).map_or(u32::MAX, |a| a.pos);
            out.probes.extend(
                probes
                    .iter()
                    .filter(|q| q.pos >= anchor.pos && q.pos < next),
            );
        }
        out.contigs.push(Contig {
            length: sequence.len() as u32,
            first,
            end: out.anchors.len() as u32,
            promoted,
        });
    }
    out.island_start.push(out.probes.len() as u32);
    out
}

fn assert_same_structure(actual: &AiipsSketch, expected: &AiipsSketch) {
    assert_eq!(actual.contigs, expected.contigs);
    assert_eq!(actual.anchors, expected.anchors);
    assert_eq!(actual.anchor_contig, expected.anchor_contig);
    assert_eq!(actual.island_start, expected.island_start);
    assert_eq!(actual.probes, expected.probes);
}

#[test]
fn builder_agrees_with_brute_force_baseline() {
    // Long enough to be scanned in interleaved lanes; the ambiguous bases
    // and the short records exercise lane warm-up and the single-lane path.
    let mut long = random_sequence(20_003, 1);
    // Ambiguous bases break k-mer runs.
    for i in (777..20_000).step_by(3001) {
        long[i] = b'N';
    }
    long[5000..5040].copy_from_slice(&[b'n'; 40]);
    // An N just before a lane boundary (20_003 / 4 = 5000.75).
    long[10_000] = b'N';
    let lower: Vec<u8> = random_sequence(3000, 2).to_ascii_lowercase();
    let short = random_sequence(12, 3);
    let records: Vec<&[u8]> = vec![&long, b"", b"NNNN", &lower, &short, b"ACGTAC"];
    for p in [
        params(7, 8, 4, 3),
        params(9, 30, 5, 4),
        params(21, 50, 15, 10),
        params(32, 20, 31, 5),
        params(2, 1, 1, 1),
    ] {
        let built = sketch(&records, &p);
        assert_same_structure(&built, &brute_force(&records, &p));
        assert_eq!(built.summary.records, records.len() as u64);
        // Every probe lies outside every k1 anchor of its record.
        let (k1, k2) = (p.anchor_k.get() as u32, p.probe_k.get() as u32);
        for t in 0..built.anchors.len() as u32 {
            let (lo, hi) = built.island_probes(t);
            for q in &built.probes[lo as usize..hi as usize] {
                let contig = built.contigs[built.anchor_contig[t as usize] as usize];
                for u in contig.first..contig.end {
                    if !built.is_promoted(u) {
                        let a = built.anchors[u as usize];
                        assert!(!(a.pos < q.pos + k2 && q.pos < a.pos + k1));
                    }
                }
            }
        }
    }
}

#[test]
fn first_sampled_kmer_of_every_record_is_an_anchor() {
    let p = params(9, 40, 5, 6);
    let records: Vec<Vec<u8>> = (0..50)
        .map(|i| random_sequence(300 + 7 * i, 100 + i as u64))
        .collect();
    let refs: Vec<&[u8]> = records.iter().map(Vec::as_slice).collect();
    let s = sketch(&refs, &p);
    for contig in &s.contigs {
        if contig.first == contig.end {
            continue;
        }
        let first_anchor = s.anchors[contig.first as usize].pos;
        let (lo, hi) = (
            s.island_start[contig.first as usize],
            s.island_start[contig.end as usize],
        );
        assert!(s.probes[lo as usize..hi as usize]
            .iter()
            .all(|q| q.pos >= first_anchor));
    }
    assert!(s.contigs.iter().any(|c| c.promoted));
}

#[test]
fn sketch_file_round_trips_and_is_bounded() {
    let records: Vec<&[u8]> = vec![b"", &[b'A'; 10]];
    let long = random_sequence(50_000, 9);
    let mut records = records;
    records.push(&long);
    let s = sketch(&records, &params(21, 50, 15, 10));
    let bytes = s.to_bytes();
    let parsed = match super::super::parse(&bytes).unwrap() {
        super::super::Sketch::Aiips(parsed) => parsed,
        other => panic!("{other:?}"),
    };
    assert_eq!(parsed, s);
    assert!(bytes.len() as u64 <= super::super::max_file_bytes(&bytes));
}

#[test]
fn parser_rejects_inconsistent_sketches() {
    let s = sketch(&[&random_sequence(5000, 4)], &params(9, 20, 5, 5));
    let text = String::from_utf8(s.to_bytes()).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    let first_probe = lines.iter().position(|l| l.starts_with("q ")).unwrap();
    let first_anchor = lines.iter().position(|l| l.starts_with("a ")).unwrap();
    let rebuild = |lines: &[String]| format!("{}\n", lines.join("\n"));
    let owned: Vec<String> = lines.iter().map(|l| l.to_string()).collect();
    let mut cases = Vec::new();
    // Probe moved before every anchor.
    let mut v = owned.clone();
    let probe = v.remove(first_probe);
    v.insert(3, probe);
    cases.push((v, "probe before"));
    // Wrong declared count.
    let mut v = owned.clone();
    v[0] = v[0].replace(" probes=", " probes=1");
    cases.push((v, "expected"));
    // Hash above the anchor threshold.
    let mut v = owned.clone();
    let fields: Vec<&str> = v[first_anchor].split(' ').collect();
    v[first_anchor] = format!("a {} {} {}", fields[1], fields[2], u64::MAX);
    cases.push((v, "threshold"));
    // Bad strand.
    let mut v = owned.clone();
    v[first_anchor] = v[first_anchor].replace(" + ", " * ").replace(" - ", " * ");
    cases.push((v, "strand"));
    // Unknown tag.
    let mut v = owned.clone();
    v[first_anchor] = v[first_anchor].replacen('a', "x", 1);
    cases.push((v, "unknown line tag"));
    // Probe overlapping the anchor before it.
    let mut v = owned.clone();
    let anchor_pos: u64 = v[first_anchor].split(' ').nth(1).unwrap().parse().unwrap();
    let mut probe: Vec<String> = v[first_probe].split(' ').map(String::from).collect();
    probe[1] = (anchor_pos + 1).to_string();
    v[first_probe] = probe.join(" ");
    cases.push((v, "overlapping"));
    for (lines, message) in cases {
        let error = super::super::parse(rebuild(&lines).as_bytes())
            .unwrap_err()
            .to_string();
        assert!(error.contains(message), "{message}: {error}");
    }
}

#[test]
fn identity_mle_matches_closed_form_and_grid_search() {
    assert_eq!(identity_mle(&[]), None);
    assert_eq!(identity_mle(&[(21, 0, 0), (15, 0, 0)]), None);
    assert_eq!(identity_mle(&[(21, 5, 5), (15, 7, 7)]), Some(1.0));
    assert_eq!(identity_mle(&[(21, 0, 5), (15, 0, 7)]), Some(0.0));
    // One k: the containment estimate (n/N)^(1/k).
    for (k, n, total) in [(15, 1, 2), (21, 37, 100), (15, 999, 1000), (4, 1, 1000)] {
        let q = identity_mle(&[(k, n, total), (21, 0, 0)]).unwrap();
        let expected = (n as f64 / total as f64).powf(1.0 / k as f64);
        assert!(
            (q - expected).abs() < 1e-12,
            "{k} {n} {total}: {q} vs {expected}"
        );
    }
    // Two k: maximum of the log-likelihood on a fine grid.
    let log_likelihood = |q: f64, obs: &[(usize, u64, u64)]| {
        obs.iter()
            .map(|&(k, n, total)| {
                let qk = q.powi(k as i32);
                n as f64 * qk.ln() + (total - n) as f64 * (1.0 - qk).ln()
            })
            .sum::<f64>()
    };
    for obs in [
        [(21, 30, 100), (15, 200, 500)],
        [(21, 1, 50), (15, 90, 400)],
        [(21, 99, 100), (15, 480, 500)],
    ] {
        let q = identity_mle(&obs).unwrap();
        let best = (1..1_000_000)
            .map(|i| i as f64 / 1_000_000.0)
            .max_by(|x, y| log_likelihood(*x, &obs).total_cmp(&log_likelihood(*y, &obs)))
            .unwrap();
        assert!((q - best).abs() < 2e-6, "{obs:?}: {q} vs {best}");
    }
}

#[test]
fn identical_genomes_have_full_identity_and_alignment() {
    let genome = random_sequence(400_000, 11);
    let e = compare(&[&genome], &[&genome]);
    assert_eq!(e.ani, Some(1.0));
    assert!(e.af_a > 0.995 && e.af_b > 0.995, "{e:?}");
    assert_eq!(e.a.probes, e.a.probes_matched);
    assert!(e.a.anchors > 100);
}

#[test]
fn strand_and_fragmentation_do_not_change_identity() {
    let genome = random_sequence(400_000, 12);
    // B: the reverse complement, cut into contigs, every other one flipped back.
    let rc = reverse_complement(&genome);
    let mut pieces: Vec<Vec<u8>> = Vec::new();
    let mut start = 0;
    for (i, length) in [50_000, 80_000, 20_000, 120_000, 60_000, 70_000]
        .into_iter()
        .enumerate()
    {
        let piece = &rc[start..start + length];
        pieces.push(if i % 2 == 0 {
            piece.to_vec()
        } else {
            reverse_complement(piece)
        });
        start += length;
    }
    let refs: Vec<&[u8]> = pieces.iter().map(Vec::as_slice).collect();
    let e = compare(&[&genome], &refs);
    assert_eq!(e.ani, Some(1.0), "{e:?}");
    // Each contig loses at most part of one island at each end.
    assert!(e.af_a > 0.99 && e.af_b > 0.99, "{e:?}");
}

#[test]
fn substitutions_are_recovered_by_the_identity_estimate() {
    let genome = random_sequence(1_000_000, 13);
    for (rate, tolerance) in [(0.01, 0.002), (0.05, 0.006)] {
        let (mutated, changed) = substitute(&genome, rate, 14);
        let truth = 1.0 - changed as f64 / genome.len() as f64;
        let e = compare(&[&genome], &[&mutated]);
        let ani = e.ani.unwrap();
        assert!(
            (ani - truth).abs() < tolerance,
            "rate {rate}: {ani} vs {truth} {e:?}"
        );
        assert!(e.af_a > 0.9 && e.af_b > 0.9, "rate {rate}: {e:?}");
    }
}

#[test]
fn unrelated_genomes_do_not_align() {
    let e = compare(
        &[&random_sequence(300_000, 15)],
        &[&random_sequence(300_000, 16)],
    );
    assert!(e.af_a < 0.005 && e.af_b < 0.005, "{e:?}");
}

#[test]
fn partial_overlap_gives_the_shared_fraction() {
    let genome = random_sequence(600_000, 17);
    let mut other = genome[..300_000].to_vec();
    other.extend(random_sequence(100_000, 18));
    let e = compare(&[&genome], &[&other]);
    // Only the islands around the breakpoint hold unshared probes.
    assert!(e.ani.unwrap() > 0.9995, "{e:?}");
    assert!((e.af_a - 0.5).abs() < 0.005, "{e:?}");
    assert!((e.af_b - 0.75).abs() < 0.005, "{e:?}");
}

#[test]
fn every_instruction_set_gives_the_brute_force_sketch_and_summary() {
    let mut long = random_sequence(100_003, 21);
    for i in [0, 12_499, 12_500, 12_501, 50_000, 99_990] {
        long[i] = b'N';
    }
    long[30_000..30_100].make_ascii_lowercase();
    let records: Vec<Vec<u8>> = vec![
        long,
        random_sequence(7, 22),
        random_sequence(25, 23),
        b"NNNNACGTNNNN".to_vec(),
        random_sequence(4_100, 24),
        Vec::new(),
    ];
    let refs: Vec<&[u8]> = records.iter().map(Vec::as_slice).collect();
    for p in [
        params(21, 50, 15, 10),
        params(32, 20, 31, 5),
        params(2, 1, 1, 1),
    ] {
        let expected = brute_force(&refs, &p);
        let k2 = p.probe_k;
        let summary = Summary {
            records: refs.len() as u64,
            bases: refs.iter().map(|r| r.len() as u64).sum(),
            non_acgt_bases: refs
                .iter()
                .map(|r| r.iter().filter(|b| !b"ACGTacgt".contains(b)).count() as u64)
                .sum(),
            kmers: refs
                .iter()
                .map(|r| crate::kmer::Kmers::new(r, k2).count() as u64)
                .sum(),
        };
        for isa in Isa::available() {
            let mut builder = Builder::new(&p);
            builder.isa = isa;
            for record in &refs {
                builder.add_record(record);
            }
            let built = builder.finish(p.clone()).unwrap();
            assert_same_structure(&built, &expected);
            assert_eq!(built.summary, summary, "{isa:?}");
        }
    }
}

#[test]
fn transpose8_transposes_byte_matrices() {
    let mut rng = SplitMix64::new(51);
    for _ in 0..1000 {
        let rows: [u64; 8] = std::array::from_fn(|_| rng.next().unwrap());
        let mut fast = rows;
        transpose8(&mut fast);
        for (j, column) in fast.iter().enumerate() {
            for (l, row) in rows.iter().enumerate() {
                assert_eq!(column.to_le_bytes()[l], row.to_le_bytes()[j]);
            }
        }
    }
}

#[test]
fn word_encoding_matches_the_lookup_table_for_every_byte() {
    for byte in 0..=255u8 {
        for position in 0..8 {
            let mut bytes = *b"ACGTacgt";
            bytes[position] = byte;
            let encoded = encode_word(u64::from_le_bytes(bytes)).to_le_bytes();
            let expected: Vec<u8> = bytes.iter().map(|&b| BASE_CODE[b as usize]).collect();
            assert_eq!(encoded.to_vec(), expected, "byte {byte} at {position}");
        }
    }
}

/// Wraps `seq` into lines of the given widths (cycled), joined by `newline`.
fn wrap(seq: &[u8], widths: &[usize], newline: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    let (mut at, mut i) = (0, 0);
    while at < seq.len() {
        let w = widths[i % widths.len()].min(seq.len() - at);
        out.extend_from_slice(&seq[at..at + w]);
        at += w;
        if at < seq.len() {
            out.extend_from_slice(newline);
        }
        i += 1;
    }
    out
}

/// Decodes a record as `sample_record` does (mode choice and redo), in
/// chunks of `chunk` bytes.
fn decode(raw: &[u8], chunk: usize) -> (Vec<u8>, bool) {
    let newlines = raw.iter().filter(|&&b| b == b'\n').count();
    let returns = raw.iter().filter(|&&b| b == b'\r').count();
    let mut mode = match (newlines, returns) {
        (0, 0) => Breaks::None,
        (_, 0) => Breaks::Predicted,
        _ => Breaks::Any,
    };
    let mut redone = false;
    loop {
        let mut decoder = Decoder::new(raw, mode);
        let (mut out, mut buf) = (Vec::new(), vec![0; chunk]);
        loop {
            let len = decoder.fill(&mut buf);
            if len == 0 {
                break;
            }
            out.extend_from_slice(&buf[..len]);
        }
        if mode == Breaks::Predicted && decoder.newlines != newlines {
            mode = Breaks::Any;
            redone = true;
            continue;
        }
        return (out, redone);
    }
}

#[test]
fn decoding_drops_line_breaks_like_needletail() {
    let reference = |raw: &[u8]| -> Vec<u8> {
        raw.iter()
            .filter(|&&b| b != b'\n' && b != b'\r')
            .map(|&b| BASE_CODE[b as usize])
            .collect()
    };
    let mut seq = random_sequence(5_000, 52);
    seq[100] = b'N';
    seq[2_000..2_010].make_ascii_lowercase();
    let mut cases: Vec<Vec<u8>> = vec![
        Vec::new(),
        b"\n".to_vec(),
        b"\n\nACGT\n".to_vec(),
        seq.clone(),
    ];
    for widths in [
        vec![1],
        vec![7],
        vec![8],
        vec![60],
        vec![80],
        vec![61, 80],
        vec![80, 80, 3],
        vec![5000],
    ] {
        cases.push(wrap(&seq, &widths, b"\n"));
        cases.push(wrap(&seq, &widths, b"\r\n"));
    }
    // A line longer than the first, a lone carriage return, empty lines.
    let mut odd = wrap(&seq, &[80], b"\n");
    odd.insert(500, b'\n');
    odd.insert(1_000, b'\r');
    odd.splice(1_500..1_500, b"\n\n".iter().copied());
    cases.push(odd);
    // Width predicted from a short first line while a later line holds a
    // break before the predicted position: decoding must be redone.
    // (After "ACGTACG", the break at index 7 of "AC\nACGT\nAC" matches the
    // predicted width but the line ends at index 2.)
    let tricky = b"ACGTACG\nAC\nACGT\nAC".to_vec();
    assert!(decode(&tricky, 64).1, "the misprediction must be caught");
    cases.push(tricky);
    for raw in &cases {
        for chunk in [1, 7, 64, 1 << 16] {
            let (codes, _) = decode(raw, chunk);
            assert_eq!(
                codes,
                reference(raw),
                "chunk {chunk}: {:?}",
                String::from_utf8_lossy(&raw[..raw.len().min(40)])
            );
        }
    }
}

#[test]
fn line_wrapped_records_give_the_brute_force_sketch() {
    // Longer than three chunks, with ambiguous bases around chunk ends.
    let mut long = random_sequence(3 * CHUNK + 1_234, 61);
    for i in [
        CHUNK - 21,
        CHUNK - 1,
        CHUNK,
        2 * CHUNK + 5,
        3 * CHUNK + 1_000,
    ] {
        long[i] = b'N';
    }
    let plain: Vec<Vec<u8>> = vec![long, random_sequence(70, 62), random_sequence(9_000, 63)];
    let refs: Vec<&[u8]> = plain.iter().map(Vec::as_slice).collect();
    let p = params(21, 50, 15, 10);
    let expected = brute_force(&refs, &p);
    for (widths, newline) in [
        (vec![60], &b"\n"[..]),
        (vec![61, 80], b"\n"),
        (vec![80], b"\r\n"),
        (vec![70, 70, 9], b"\n"),
        (vec![CHUNK + 3], b"\n"),
    ] {
        let wrapped: Vec<Vec<u8>> = plain.iter().map(|r| wrap(r, &widths, newline)).collect();
        for isa in Isa::available() {
            let mut builder = Builder::new(&p);
            builder.isa = isa;
            for record in &wrapped {
                builder.add_record(record);
            }
            let built = builder.finish(p.clone()).unwrap();
            assert_same_structure(&built, &expected);
            assert_eq!(
                built.summary.bases,
                plain.iter().map(|r| r.len() as u64).sum::<u64>()
            );
        }
    }
}
