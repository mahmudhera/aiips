//! The aiips file: the shared header and summary text lines, then a binary
//! body of tagged sections (see [`crate::sketch::binary`]) in this order:
//!
//! ```text
//! GPRB  distinct probe hashes (the global set and probe dictionary):
//!       n, then h_1, h_2 - h_1, ..., h_n - h_(n-1)
//! GANC  distinct anchor identity hashes (anchor dictionary), same encoding;
//!       these are full 64-bit hashes, not bounded by a threshold
//! CTGS  n, then per contig: length, anchor count, probe count
//! ANCH  per anchor: position minus the previous anchor's on its contig
//!       (the first from 0); per anchor: dictionary index; strand bits
//! PROB  as ANCH, then role bits
//! ```
//!
//! Integers are minimal LEB128 varints; bit columns are packed least
//! significant bit first with zero padding. Island links are not stored:
//! they are rebuilt from positions on load.

use super::occurrences::{Bits, Contigs, Occurrences};
use super::{check_k, check_k_order, check_scaled, AiipsParams, AiipsSketch};
use crate::fastx::Summary;
use crate::sketch::binary::{put_varint, BodyReader, BodyWriter, Section};
use crate::sketch::fracminhash::max_hash;
use crate::sketch::{invalid, Fields, ParseError};

const GLOBAL_PROBES: &[u8; 4] = b"GPRB";
const ANCHOR_HASHES: &[u8; 4] = b"GANC";
const CONTIGS: &[u8; 4] = b"CTGS";
const ANCHORS: &[u8; 4] = b"ANCH";
const PROBES: &[u8; 4] = b"PROB";

impl AiipsSketch {
    pub fn to_bytes(&self) -> Vec<u8> {
        use std::fmt::Write;
        let mut contigs = Vec::with_capacity(1 + 9 * self.contigs.len());
        put_varint(&mut contigs, self.contigs.len() as u64);
        for contig in 0..self.contigs.len() {
            put_varint(&mut contigs, u64::from(self.contigs.length(contig)));
            put_varint(&mut contigs, self.contigs.anchors(contig).len() as u64);
            put_varint(&mut contigs, self.contigs.probes(contig).len() as u64);
        }
        let mut body = BodyWriter::default();
        body.section(GLOBAL_PROBES, &dictionary_bytes(&self.probes.hashes));
        body.section(ANCHOR_HASHES, &dictionary_bytes(&self.anchors.hashes));
        body.section(CONTIGS, &contigs);
        body.section(
            ANCHORS,
            &occurrence_bytes(&self.contigs, Kind::Anchor, &self.anchors, None),
        );
        body.section(
            PROBES,
            &occurrence_bytes(
                &self.contigs,
                Kind::Probe,
                &self.probes,
                Some(&self.probe_roles),
            ),
        );
        let body = body.finish();

        let p = &self.params;
        let mut text = String::with_capacity(300);
        let _ = writeln!(
            text,
            "{} {} algorithm=aiips anchor_k={} anchor_scaled={} probe_k={} probe_scaled={} seed={} hash={} canonical={} body_bytes={}",
            crate::sketch::MAGIC,
            crate::sketch::VERSION,
            p.anchor_k,
            p.anchor_scaled,
            p.probe_k,
            p.probe_scaled,
            p.seed,
            p.hash,
            p.canonical,
            body.len()
        );
        crate::sketch::write_summary(&mut text, &self.summary);
        let mut out = text.into_bytes();
        out.extend_from_slice(&body);
        out
    }

    /// Parses everything after the header line: the summary line, then the
    /// binary body, which must be exactly `body_bytes` long.
    pub(crate) fn parse_body(mut header: Fields<'_>, rest: &[u8]) -> Result<Self, ParseError> {
        let anchor_k = check_k("anchor", header.take_u64("anchor_k")?)?;
        let anchor_scaled = header.take_u64("anchor_scaled")?;
        check_scaled("anchor_scaled", anchor_scaled)?;
        let probe_k = check_k("probe", header.take_u64("probe_k")?)?;
        let probe_scaled = header.take_u64("probe_scaled")?;
        check_scaled("probe_scaled", probe_scaled)?;
        check_k_order(anchor_k, probe_k)?;
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
        let body_bytes = header.take_u64("body_bytes")?;
        header.finish()?;
        let params = AiipsParams {
            anchor_k,
            anchor_scaled,
            probe_k,
            probe_scaled,
            seed,
            hash,
            canonical,
        };

        let newline = rest
            .iter()
            .position(|&b| b == b'\n')
            .ok_or_else(|| invalid("missing summary line"))?;
        let line = std::str::from_utf8(&rest[..newline])
            .map_err(|_| invalid("summary line is not valid UTF-8"))?;
        if line.contains('\r') {
            return Err(invalid("carriage return found; lines must end with \\n"));
        }
        let summary = crate::sketch::parse_summary(Some(line))?;
        let body = &rest[newline + 1..];
        if body.len() as u64 != body_bytes {
            return Err(invalid(format!(
                "expected body_bytes={body_bytes}, found {}",
                body.len()
            )));
        }

        let mut reader = BodyReader::new(body)?;
        let probe_hashes = read_dictionary(
            reader.section(GLOBAL_PROBES)?,
            Kind::Probe,
            Some(probe_scaled),
        )?;
        let anchor_hashes = read_dictionary(reader.section(ANCHOR_HASHES)?, Kind::Anchor, None)?;
        let contigs = read_contigs(reader.section(CONTIGS)?, &summary)?;
        let (anchors, _) = read_occurrences(
            reader.section(ANCHORS)?,
            Kind::Anchor,
            &contigs,
            anchor_k.get(),
            anchor_hashes,
        )?;
        let (probes, roles) = read_occurrences(
            reader.section(PROBES)?,
            Kind::Probe,
            &contigs,
            probe_k.get(),
            probe_hashes,
        )?;
        reader.finish()?;

        let windows = contigs.windows(probe_k.get());
        if summary.kmers > windows {
            return Err(invalid(format!(
                "valid_kmers={} exceeds the {windows} probe-length windows of the contigs",
                summary.kmers
            )));
        }
        if probes.len() as u64 > summary.kmers {
            return Err(invalid(format!(
                "{} probe occurrences exceed valid_kmers={}",
                probes.len(),
                summary.kmers
            )));
        }
        Ok(Self::assemble(
            params,
            summary,
            contigs,
            anchors,
            probes,
            roles.expect("probe sections carry roles"),
        ))
    }
}

#[derive(Clone, Copy)]
enum Kind {
    Anchor,
    Probe,
}

impl Kind {
    fn name(self) -> &'static str {
        match self {
            Kind::Anchor => "anchor",
            Kind::Probe => "probe",
        }
    }

    fn range(self, contigs: &Contigs, contig: usize) -> std::ops::Range<usize> {
        match self {
            Kind::Anchor => contigs.anchors(contig),
            Kind::Probe => contigs.probes(contig),
        }
    }

    fn total(self, contigs: &Contigs) -> usize {
        let ends = match self {
            Kind::Anchor => &contigs.anchor_ends,
            Kind::Probe => &contigs.probe_ends,
        };
        ends.last().copied().unwrap_or(0)
    }
}

fn dictionary_bytes(hashes: &[u64]) -> Vec<u8> {
    let mut out = Vec::with_capacity(10 + 8 * hashes.len());
    put_varint(&mut out, hashes.len() as u64);
    let mut previous = 0;
    for &hash in hashes {
        put_varint(&mut out, hash - previous);
        previous = hash;
    }
    out
}

fn occurrence_bytes(
    contigs: &Contigs,
    kind: Kind,
    occurrences: &Occurrences,
    roles: Option<&Bits>,
) -> Vec<u8> {
    let mut out = Vec::with_capacity(5 * occurrences.len());
    for contig in 0..contigs.len() {
        let mut previous = 0;
        for &position in &occurrences.positions[kind.range(contigs, contig)] {
            put_varint(&mut out, u64::from(position - previous));
            previous = position;
        }
    }
    for &id in &occurrences.ids {
        put_varint(&mut out, u64::from(id));
    }
    out.extend_from_slice(occurrences.strands.as_bytes());
    if let Some(roles) = roles {
        out.extend_from_slice(roles.as_bytes());
    }
    out
}

/// A strictly increasing delta-encoded hash list; with `scaled`, each value
/// must also be at most floor((2^64 − 1) / scaled).
fn read_dictionary(
    mut section: Section<'_>,
    kind: Kind,
    scaled: Option<u64>,
) -> Result<Vec<u64>, ParseError> {
    let name = kind.name();
    let size = section.varint()?;
    // The count is untrusted: every value takes at least one byte.
    if size > section.remaining() as u64 {
        return Err(invalid(format!(
            "{name} dictionary declares {size} hashes but holds at most {}",
            section.remaining()
        )));
    }
    let max = scaled.map_or(u64::MAX, max_hash);
    let mut hashes: Vec<u64> = Vec::with_capacity(size as usize);
    for index in 1..=size {
        let delta = section.varint()?;
        let value = match hashes.last() {
            None => Some(delta),
            Some(_) if delta == 0 => {
                return Err(invalid(format!(
                    "{name} hash {index} is not strictly increasing"
                )))
            }
            Some(&previous) => previous.checked_add(delta),
        }
        .filter(|&value| value <= max)
        .ok_or_else(|| {
            invalid(format!(
                "{name} hash {index} exceeds max_hash={max} for {name}_scaled={}",
                scaled.unwrap_or(1)
            ))
        })?;
        hashes.push(value);
    }
    section.finish()?;
    // Occurrences index the dictionary with u32.
    if hashes.len() as u64 > u64::from(u32::MAX) + 1 {
        return Err(invalid(format!("more than 2^32 distinct {name} hashes")));
    }
    Ok(hashes)
}

fn read_contigs(mut section: Section<'_>, summary: &Summary) -> Result<Contigs, ParseError> {
    let count = section.varint()?;
    if count != summary.records {
        return Err(invalid(format!(
            "section CTGS lists {count} contigs, but records={}",
            summary.records
        )));
    }
    // Every contig takes at least three bytes.
    if count > (section.remaining() / 3) as u64 {
        return Err(invalid(format!(
            "section CTGS declares {count} contigs but holds at most {}",
            section.remaining() / 3
        )));
    }
    let overflow = || invalid("section CTGS: occurrence counts overflow");
    let mut contigs = Contigs::default();
    let (mut anchors, mut probes, mut bases) = (0usize, 0usize, 0u64);
    for contig in 1..=count {
        let length = section.varint()?;
        let length = u32::try_from(length).map_err(|_| {
            invalid(format!(
                "contig {contig} has {length} bases; at most {} are supported",
                u32::MAX
            ))
        })?;
        for total in [&mut anchors, &mut probes] {
            let n = section.varint()?;
            *total = usize::try_from(n)
                .ok()
                .and_then(|n| total.checked_add(n))
                .ok_or_else(overflow)?;
        }
        bases += u64::from(length);
        contigs.push(length, anchors, probes);
    }
    section.finish()?;
    if bases != summary.bases {
        return Err(invalid(format!(
            "contig lengths sum to {bases}, but bases={}",
            summary.bases
        )));
    }
    Ok(contigs)
}

/// One occurrence column set; probes also return their role bits.
fn read_occurrences(
    mut section: Section<'_>,
    kind: Kind,
    contigs: &Contigs,
    k: usize,
    hashes: Vec<u64>,
) -> Result<(Occurrences, Option<Bits>), ParseError> {
    let name = kind.name();
    let total = kind.total(contigs);
    // Every position takes at least one byte; bound before allocating.
    if total > section.remaining() {
        return Err(invalid(format!(
            "{total} {name} occurrences declared, but their section holds {} bytes",
            section.remaining()
        )));
    }
    let mut positions = Vec::with_capacity(total);
    for contig in 0..contigs.len() {
        let length = u64::from(contigs.length(contig));
        let mut previous: Option<u64> = None;
        for _ in kind.range(contigs, contig) {
            let index = positions.len() + 1;
            let delta = section.varint()?;
            let position = match previous {
                None => Some(delta),
                Some(_) if delta == 0 => {
                    return Err(invalid(format!(
                        "{name} {index} is not after the previous {name} on contig {}",
                        contig + 1
                    )))
                }
                Some(previous) => previous.checked_add(delta),
            }
            .filter(|&position| position.saturating_add(k as u64) <= length)
            .ok_or_else(|| {
                invalid(format!(
                    "{name} {index}: a {k}-mer there does not fit in contig {} of {length} bases",
                    contig + 1
                ))
            })?;
            positions.push(position as u32);
            previous = Some(position);
        }
    }
    let mut used = vec![false; hashes.len()];
    let mut ids = Vec::with_capacity(total);
    for index in 1..=total {
        let id = section.varint()?;
        let slot = usize::try_from(id)
            .ok()
            .filter(|&slot| slot < hashes.len())
            .ok_or_else(|| {
                invalid(format!(
                    "{name} {index} refers to dictionary entry {id} of {}",
                    hashes.len()
                ))
            })?;
        used[slot] = true;
        ids.push(slot as u32);
    }
    if let Some(unused) = used.iter().position(|&used| !used) {
        return Err(invalid(format!(
            "{name} hash {} is in the dictionary but has no occurrence",
            unused + 1
        )));
    }
    let mut bits = |column: &str| -> Result<Bits, ParseError> {
        Bits::from_bytes(section.bytes(total.div_ceil(8))?, total).ok_or_else(|| {
            invalid(format!(
                "{name} {column} bits: padding after the last bit must be zero"
            ))
        })
    };
    let strands = bits("strand")?;
    let roles = match kind {
        Kind::Probe => Some(bits("role")?),
        Kind::Anchor => None,
    };
    section.finish()?;
    Ok((
        Occurrences {
            hashes,
            ids,
            positions,
            strands,
        },
        roles,
    ))
}

#[cfg(test)]
mod tests {
    use super::super::tests::{build, fasta, island_genome, params, random_bases};
    use super::super::Role;
    use super::*;
    use crate::sketch::{parse, Sketch};

    fn hex(bytes: &[u8]) -> String {
        bytes.iter().map(|b| format!("{b:02x}")).collect()
    }

    /// tests/fixtures/sample.fa at anchor_k=5, anchor_scaled=2, probe_k=3,
    /// probe_scaled=1, seed=42, as computed by an independent Python
    /// implementation of SplitMix64, fmix64, stranded 2-bit canonical k-mers,
    /// LEB128, bit packing and FNV-1a 64, written from their definitions and
    /// this layout. Guards against hash or format drift.
    const GOLDEN: &str = concat!(
        "61696970732d736b6574636820763120616c676f726974686d3d616969707320616e63686f725f6b",
        "3d3520616e63686f725f7363616c65643d322070726f62655f6b3d332070726f62655f7363616c65",
        "643d3120736565643d343220686173683d6d75726d7572332d666d697836342063616e6f6e696361",
        "6c3d7472756520626f64795f62797465733d3133330a7265636f7264733d322062617365733d3133",
        "206e6f6e5f616367745f62617365733d322076616c69645f6b6d6572733d370a475052421c000000",
        "0000000003dbc4aeabe0a1b2e126ecff9bf79892df8a65e6b6f7cfbef4e3ac5947414e430b000000",
        "0000000001fcd5d8ec9588ac82ba0143544753070000000000000002080004050103414e43480300",
        "00000000000000000150524f42100000000000000000010101000101020201010001023a53b9cd92",
        "f8ed800690"
    );

    #[test]
    fn sample_fixture_matches_independent_golden_file() {
        let fixture = include_bytes!("../../../tests/fixtures/sample.fa");
        let sketch = build(fixture, params(5, 2, 3, 1, 42));
        assert_eq!(
            sketch.global_probes(),
            [
                2793015780980793947,
                10076880717932044871,
                16515215794403048877
            ]
        );
        assert_eq!(sketch.anchors().hashes, [13404032185750924028]);
        assert_eq!(sketch.contigs().anchors(1), 0..1);
        let bytes = sketch.to_bytes();
        assert_eq!(hex(&bytes), GOLDEN);
        // Independent runs give identical files.
        assert_eq!(build(fixture, params(5, 2, 3, 1, 42)).to_bytes(), bytes);
    }

    fn parsed(bytes: &[u8]) -> AiipsSketch {
        match parse(bytes).unwrap() {
            Sketch::Aiips(sketch) => *sketch,
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn serialization_round_trips_exactly() {
        let long = random_bases(5, 30_000);
        let file = fasta(&[&long, b"ACGTNNNNacgt", b"GA", &long[..500]]);
        let empty = build(b">a\nACGTACGT\n", params(5, u64::MAX, 3, u64::MAX, 1));
        assert!(empty.global_probes().is_empty() && empty.anchors().is_empty());
        for sketch in [
            build(&file, params(21, 125, 15, 125, 42)),
            build(&file, params(17, 3, 15, 1, 7)),
            build(&file, params(31, 1, 29, 3, u64::MAX)),
            island_genome(),
            empty,
        ] {
            let bytes = sketch.to_bytes();
            let back = parsed(&bytes);
            assert_eq!(back, sketch);
            assert_eq!(back.to_bytes(), bytes);
        }
    }

    fn varints(values: &[u64]) -> Vec<u8> {
        let mut out = Vec::new();
        for &value in values {
            put_varint(&mut out, value);
        }
        out
    }

    fn cat(parts: &[&[u8]]) -> Vec<u8> {
        parts.concat()
    }

    const FIELDS: &str = "anchor_k=5 anchor_scaled=1 probe_k=3 probe_scaled=1 seed=42 hash=murmur3-fmix64 canonical=true";
    /// Two 10-base contigs: 8 + 8 windows of length 3.
    const SUMMARY: &str = "records=2 bases=20 non_acgt_bases=0 valid_kmers=16";

    /// A file with a correct checksum and body_bytes, so each malformed
    /// variant reaches the rule it targets.
    fn file(fields: &str, summary: &str, sections: &[(&[u8; 4], Vec<u8>)]) -> Vec<u8> {
        let mut writer = BodyWriter::default();
        for (tag, payload) in sections {
            writer.section(tag, payload);
        }
        let body = writer.finish();
        let mut out = format!(
            "aiips-sketch v1 algorithm=aiips {fields} body_bytes={}\n{summary}\n",
            body.len()
        )
        .into_bytes();
        out.extend_from_slice(&body);
        out
    }

    /// Contig 0: anchor at 4, probes at 0 (START) and 6 (owned by the
    /// anchor); contig 1: no anchor, one probe at 7.
    fn sections() -> Vec<(&'static [u8; 4], Vec<u8>)> {
        vec![
            (GLOBAL_PROBES, varints(&[2, 5, 7])),
            (ANCHOR_HASHES, varints(&[1, 9])),
            (CONTIGS, varints(&[2, 10, 1, 2, 10, 0, 1])),
            (ANCHORS, cat(&[&varints(&[4, 0]), &[0b1]])),
            (
                PROBES,
                cat(&[&varints(&[0, 6, 7, 0, 1, 1]), &[0b010], &[0b101]]),
            ),
        ]
    }

    fn with(tag: &[u8; 4], payload: Vec<u8>) -> Vec<(&'static [u8; 4], Vec<u8>)> {
        let mut sections = sections();
        sections.iter_mut().find(|(t, _)| *t == tag).unwrap().1 = payload;
        sections
    }

    fn rejects(bytes: &[u8], message: &str) {
        let error = parse(bytes).unwrap_err().to_string();
        assert!(error.contains(message), "expected {message:?} in {error:?}");
    }

    #[test]
    fn parsing_accepts_a_hand_built_file_and_links_islands() {
        let sketch = parsed(&file(FIELDS, SUMMARY, &sections()));
        assert_eq!(sketch.global_probes(), [5, 12]);
        assert_eq!(sketch.anchors().hashes, [9]);
        assert_eq!(sketch.anchors().positions, [4]);
        assert!(sketch.anchors().strands.get(0));
        assert_eq!(sketch.probes().positions, [0, 6, 7]);
        assert_eq!(
            (0..3).map(|p| sketch.probes().hash(p)).collect::<Vec<_>>(),
            [5, 12, 12]
        );
        assert_eq!(
            (0..3).map(|p| sketch.probe_role(p)).collect::<Vec<_>>(),
            [Role::Scoring, Role::Mapping, Role::Scoring]
        );
        assert_eq!(sketch.start_bucket(0), 0..1);
        assert_eq!(sketch.owned_probes(0), 1..2);
        assert_eq!(sketch.next_anchor(0), None);
        assert_eq!(sketch.start_bucket(1), 2..3);
        assert_eq!(sketch.ambiguous_kmers(), 0);
        assert_eq!(sketch.params, params(5, 1, 3, 1, 42));
    }

    #[test]
    fn parsing_rejects_inconsistent_sections() {
        let max = u64::MAX;
        for (sections, message) in [
            // Dictionaries.
            (
                with(ANCHOR_HASHES, varints(&[2, 9, 0])),
                "anchor hash 2 is not strictly increasing",
            ),
            (
                with(GLOBAL_PROBES, varints(&[2, max, 1])),
                "probe hash 2 exceeds max_hash",
            ),
            (
                with(GLOBAL_PROBES, varints(&[3, 5])),
                "probe dictionary declares 3 hashes but holds at most 1",
            ),
            // Contig table.
            (
                with(CONTIGS, varints(&[1, 20, 1, 3])),
                "section CTGS lists 1 contigs, but records=2",
            ),
            (
                with(CONTIGS, varints(&[2, 10, 1, 2, 11, 0, 1])),
                "contig lengths sum to 21, but bases=20",
            ),
            (
                with(CONTIGS, varints(&[2, 1 << 32, 1, 2, 10, 0, 1])),
                "contig 1 has 4294967296 bases; at most 4294967295 are supported",
            ),
            (
                with(CONTIGS, varints(&[2, 10, max, 2, 10, max, 1])),
                "section CTGS: occurrence counts overflow",
            ),
            (
                with(CONTIGS, varints(&[2, 10, 9, 2, 10, 0, 1])),
                "9 anchor occurrences declared, but their section holds 3 bytes",
            ),
            (
                with(CONTIGS, varints(&[2, 10, 1, 2, 10, 0, 1, 0])),
                "section CTGS: 1 unread bytes",
            ),
            // Occurrence positions.
            (
                with(ANCHORS, cat(&[&varints(&[8, 0]), &[0b1]])),
                "anchor 1: a 5-mer there does not fit in contig 1 of 10 bases",
            ),
            (
                with(PROBES, cat(&[&varints(&[0, 0, 7, 0, 1, 1]), &[0], &[0]])),
                "probe 2 is not after the previous probe on contig 1",
            ),
            (
                with(PROBES, cat(&[&varints(&[0, 6, 8, 0, 1, 1]), &[0], &[0]])),
                "probe 3: a 3-mer there does not fit in contig 2 of 10 bases",
            ),
            // Dictionary indices.
            (
                with(PROBES, cat(&[&varints(&[0, 6, 7, 0, 1, 2]), &[0], &[0]])),
                "probe 3 refers to dictionary entry 2 of 2",
            ),
            (
                with(PROBES, cat(&[&varints(&[0, 6, 7, 0, 0, 0]), &[0], &[0]])),
                "probe hash 2 is in the dictionary but has no occurrence",
            ),
            // Bit columns.
            (
                with(ANCHORS, cat(&[&varints(&[4, 0]), &[0b11]])),
                "anchor strand bits: padding after the last bit must be zero",
            ),
            (
                with(
                    PROBES,
                    cat(&[&varints(&[0, 6, 7, 0, 1, 1]), &[0], &[0b1000]]),
                ),
                "probe role bits: padding after the last bit must be zero",
            ),
            (
                with(PROBES, cat(&[&varints(&[0, 6, 7, 0, 1, 1]), &[0]])),
                "needs 1 more bytes, found 0",
            ),
            (
                with(ANCHORS, cat(&[&varints(&[4, 0]), &[0b1, 0]])),
                "section ANCH: 1 unread bytes",
            ),
        ] {
            rejects(&file(FIELDS, SUMMARY, &sections), message);
        }
        // Section order and presence.
        let mut missing = sections();
        missing.pop();
        rejects(&file(FIELDS, SUMMARY, &missing), "missing section PROB");
        let mut swapped = sections();
        swapped.swap(0, 1);
        rejects(
            &file(FIELDS, SUMMARY, &swapped),
            "expected section GPRB, found GANC",
        );
        let mut extra = sections();
        extra.push((b"MORE", vec![]));
        rejects(
            &file(FIELDS, SUMMARY, &extra),
            "unexpected bytes after the last section",
        );
        // Anchor identities are full 64-bit hashes, bounded by no threshold.
        let full = parsed(&file(
            FIELDS,
            SUMMARY,
            &with(ANCHOR_HASHES, varints(&[1, u64::MAX])),
        ));
        assert_eq!(full.anchors().hashes, [u64::MAX]);
        // Probe hashes stay bounded by their threshold.
        rejects(
            &file(
                &FIELDS.replace("probe_scaled=1", "probe_scaled=2"),
                SUMMARY,
                &with(GLOBAL_PROBES, varints(&[2, 5, 1 << 63])),
            ),
            "probe hash 2 exceeds max_hash=9223372036854775807 for probe_scaled=2",
        );
        // Anchors must be longer than probes.
        rejects(
            &file(
                &FIELDS.replace("anchor_k=5", "anchor_k=3"),
                SUMMARY,
                &sections(),
            ),
            "aiips anchor k-mer length (3) must exceed the probe k-mer length (3)",
        );
    }

    #[test]
    fn parsing_checks_counts_against_the_summary_line() {
        rejects(
            &file(
                FIELDS,
                &SUMMARY.replace("valid_kmers=16", "valid_kmers=17"),
                &sections(),
            ),
            "valid_kmers=17 exceeds the 16 probe-length windows of the contigs",
        );
        rejects(
            &file(
                FIELDS,
                &SUMMARY.replace("valid_kmers=16", "valid_kmers=2"),
                &sections(),
            ),
            "3 probe occurrences exceed valid_kmers=2",
        );
        rejects(
            &file(
                FIELDS,
                &SUMMARY.replace("records=2", "records=1000"),
                &sections(),
            ),
            "section CTGS lists 2 contigs, but records=1000",
        );
    }

    #[test]
    fn parsing_rejects_corrupted_or_resized_bodies() {
        let good = file(FIELDS, SUMMARY, &sections());
        // The body starts after the second newline (header, summary).
        let summary_end = good
            .iter()
            .enumerate()
            .filter(|(_, &b)| b == b'\n')
            .nth(1)
            .unwrap()
            .0
            + 1;
        let body_bytes = good.len() - summary_end;
        for index in summary_end..good.len() {
            let mut bad = good.clone();
            bad[index] ^= 0x01;
            rejects(&bad, "checksum mismatch");
        }
        let mut longer = good.clone();
        longer.push(0);
        rejects(
            &longer,
            &format!("expected body_bytes={body_bytes}, found {}", body_bytes + 1),
        );
        rejects(
            &good[..good.len() - 1],
            &format!("expected body_bytes={body_bytes}, found {}", body_bytes - 1),
        );
    }

    #[test]
    fn parsing_rejects_malformed_header_and_summary_lines() {
        let sections = sections();
        for (fields, message) in [
            (
                FIELDS.replace("anchor_k=5", "anchor_k=4"),
                "aiips anchor k-mer length must be an odd integer from 1 to 31, got 4",
            ),
            (FIELDS.replace("anchor_k=5", "anchor_k=33"), "got 33"),
            (
                FIELDS.replace("anchor_scaled=1", "anchor_scaled=0"),
                "anchor_scaled must be an integer from 1",
            ),
            (
                FIELDS.replace("probe_k=3", "probe_k=4"),
                "aiips probe k-mer length must be an odd integer from 1 to 31, got 4",
            ),
            (FIELDS.replace("probe_k=3", "probe_k=0"), "got 0"),
            (
                FIELDS.replace("probe_k=3", "probe_k=18446744073709551615"),
                "got 18446744073709551615",
            ),
            (
                FIELDS.replace("probe_scaled=1", "probe_scaled=0"),
                "probe_scaled must be an integer from 1",
            ),
            (
                FIELDS.replace("anchor_k=5 ", ""),
                "missing header field: anchor_k",
            ),
            (FIELDS.replace("seed=42 ", ""), "missing header field: seed"),
            (format!("{FIELDS} seed=1"), "duplicate header field: seed"),
            (format!("{FIELDS} extra=1"), "unknown header field: extra"),
            (
                FIELDS.replace("canonical=true", "canonical=yes"),
                "canonical must be true or false",
            ),
            (
                FIELDS.replace("seed=42", "seed=-1"),
                "seed is not an unsigned integer",
            ),
        ] {
            rejects(&file(&fields, SUMMARY, &sections), message);
        }
        let good = file(FIELDS, SUMMARY, &sections);
        let text = String::from_utf8_lossy(&good);
        let declared = text
            .split("body_bytes=")
            .nth(1)
            .unwrap()
            .split('\n')
            .next()
            .unwrap();
        let mut no_body_bytes = good.clone();
        let field = format!(" body_bytes={declared}");
        let at = good
            .windows(field.len())
            .position(|w| w == field.as_bytes())
            .unwrap();
        no_body_bytes.drain(at..at + field.len());
        rejects(&no_body_bytes, "missing header field: body_bytes");
        rejects(
            &file(
                FIELDS,
                &SUMMARY.replace("valid_kmers=16", "valid_kmers=0"),
                &sections,
            ),
            "no valid k-mers",
        );
        rejects(
            &file(FIELDS, &format!("{SUMMARY}\r"), &sections),
            "carriage return",
        );
        rejects(
            &file(FIELDS, "records=2", &sections),
            "missing summary field",
        );

        let newline = good.iter().position(|&b| b == b'\n').unwrap();
        let mut crlf = good.clone();
        crlf.insert(newline, b'\r');
        rejects(&crlf, "carriage return");
        let mut not_utf8 = good.clone();
        not_utf8.insert(newline, 0xff);
        rejects(&not_utf8, "not valid UTF-8");
        rejects(&good[..newline], "header line must end with \\n");
        rejects(&good[..newline + 1], "missing summary line");
        rejects(
            &[b"aiips-sketch v2", &good[15..]].concat(),
            "unsupported format version",
        );
    }
}
