//! AIIPS (Anchor-Indexed Island-Probe Sketching): a positional two-level
//! FracMinHash sketch that estimates ANI and alignment fraction (AF).
//!
//! Sketch (per record, i.e. per contig):
//! - Anchors: canonical k1-mers with h1(x) = fmix64(x XOR seed_0) at or below
//!   floor((2^64 − 1) / S1), i.e. FracMinHash at k1 with scale S1.
//! - Probes: canonical k2-mers (k2 < k1) with h2(x) = fmix64(x XOR seed_1) at or
//!   below floor((2^64 − 1) / S2), excluding every k2-mer that overlaps an
//!   anchor. seed_0 and seed_1 are the first two SplitMix64 outputs of `seed`,
//!   so the two samplings use independent hash functions.
//! - Islands: the stretch from one anchor to the next anchor (or the record
//!   end). Every probe lies in exactly one island.
//! - The record's first probe that precedes its first anchor (or its first
//!   probe, if it has no anchor) is promoted to an anchor, so the first
//!   sampled k-mer of every record is an anchor and every probe follows one.
//!
//! Every sampled k-mer keeps its record position and strand; island lengths
//! follow from the anchor positions.
//!
//! Comparison: see [`estimate`].

use std::collections::HashSet;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::Read;
use std::num::NonZeroUsize;

use super::minhash::BuildError;
use super::{invalid, parse_u64, Fields, ParseError};
use crate::fastx::{for_each_raw_sequence, Summary};
use crate::hash::{fmix64, SplitMix64};
use crate::kmer::MAX_PACKED_K;

pub const HASH: &str = "murmur3-fmix64";
pub const DEFAULT_ANCHOR_K: usize = 21;
pub const DEFAULT_ANCHOR_SCALED: u64 = 500;
pub const DEFAULT_PROBE_K: usize = 15;
pub const DEFAULT_PROBE_SCALED: u64 = 100;
pub const DEFAULT_SEED: u64 = 42;
pub const DEFAULT_MAX_GAP: u32 = 50;
pub const DEFAULT_MIN_ISLAND_ANI: f64 = 0.80;
/// A seed hash occurring more often than this in either sketch is skipped:
/// highly repetitive k-mers would otherwise start many futile walks.
pub const MAX_SEED_OCCURRENCES: usize = 32;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AiipsParams {
    /// k1, the anchor length.
    pub anchor_k: NonZeroUsize,
    /// S1: about 1/S1 of the k1-mers are anchors.
    pub anchor_scaled: u64,
    /// k2, the probe length; shorter than k1.
    pub probe_k: NonZeroUsize,
    /// S2: about 1/S2 of the k2-mers are sampled (before anchor exclusion).
    pub probe_scaled: u64,
    pub seed: u64,
    pub hash: String,
    pub canonical: bool,
}

fn check(
    anchor_k: usize,
    anchor_scaled: u64,
    probe_k: usize,
    probe_scaled: u64,
) -> Result<(), ParseError> {
    if !(2..=MAX_PACKED_K).contains(&anchor_k) {
        return Err(invalid(format!(
            "aiips anchor k must be from 2 to {MAX_PACKED_K}, got {anchor_k}"
        )));
    }
    if !(1..anchor_k).contains(&probe_k) {
        return Err(invalid(format!(
            "aiips probe k must be from 1 to anchor k - 1 ({}), got {probe_k}",
            anchor_k - 1
        )));
    }
    for (name, scaled) in [("anchor", anchor_scaled), ("probe", probe_scaled)] {
        if scaled == 0 {
            return Err(invalid(format!(
                "aiips {name} scaled must be an integer from 1 to {}",
                u64::MAX
            )));
        }
    }
    Ok(())
}

impl AiipsParams {
    pub fn new(
        anchor_k: NonZeroUsize,
        anchor_scaled: u64,
        probe_k: NonZeroUsize,
        probe_scaled: u64,
        seed: u64,
    ) -> Result<Self, String> {
        check(anchor_k.get(), anchor_scaled, probe_k.get(), probe_scaled)
            .map_err(|e| e.to_string())?;
        Ok(Self {
            anchor_k,
            anchor_scaled,
            probe_k,
            probe_scaled,
            seed,
            hash: HASH.into(),
            canonical: true,
        })
    }

    pub fn anchor_max_hash(&self) -> u64 {
        u64::MAX / self.anchor_scaled
    }

    pub fn probe_max_hash(&self) -> u64 {
        u64::MAX / self.probe_scaled
    }

    /// (anchor seed, probe seed): the first two SplitMix64 outputs.
    fn hash_seeds(&self) -> (u64, u64) {
        let mut seeds = SplitMix64::new(self.seed);
        (seeds.next().unwrap(), seeds.next().unwrap())
    }

    /// Names of the parameters that differ; sketches compare only if empty.
    pub fn mismatches(&self, other: &Self) -> Vec<String> {
        let mut differences = Vec::new();
        let mut field = |name: &str, a: String, b: String| {
            if a != b {
                differences.push(format!("{name} ({a} vs {b})"));
            }
        };
        field(
            "anchor_k",
            self.anchor_k.to_string(),
            other.anchor_k.to_string(),
        );
        field(
            "anchor_scaled",
            self.anchor_scaled.to_string(),
            other.anchor_scaled.to_string(),
        );
        field(
            "probe_k",
            self.probe_k.to_string(),
            other.probe_k.to_string(),
        );
        field(
            "probe_scaled",
            self.probe_scaled.to_string(),
            other.probe_scaled.to_string(),
        );
        field("seed", self.seed.to_string(), other.seed.to_string());
        field("hash", self.hash.clone(), other.hash.clone());
        field(
            "canonical",
            self.canonical.to_string(),
            other.canonical.to_string(),
        );
        differences
    }

    fn check_buildable(&self) -> Result<(), BuildError> {
        check(
            self.anchor_k.get(),
            self.anchor_scaled,
            self.probe_k.get(),
            self.probe_scaled,
        )
        .map_err(|e| BuildError::Unsupported(e.to_string()))?;
        if self.hash != HASH || !self.canonical {
            return Err(BuildError::Unsupported(format!(
                "aiips can only build hash={HASH} canonical=true sketches"
            )));
        }
        Ok(())
    }
}

/// One sampled k-mer: zero-based start in its record, strand, and hash.
/// `reverse` is true when the canonical form is the reverse complement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Sampled {
    pub pos: u32,
    pub reverse: bool,
    pub hash: u64,
}

/// A record: its length and its anchors `first..end`. When `promoted`, anchor
/// `first` is the promoted probe (a k2-mer, hashed with the probe hash).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Contig {
    pub length: u32,
    pub first: u32,
    pub end: u32,
    pub promoted: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AiipsSketch {
    pub params: AiipsParams,
    pub summary: Summary,
    pub contigs: Vec<Contig>,
    /// Anchors of all records, in record order then position order.
    pub anchors: Vec<Sampled>,
    /// Record index of each anchor.
    pub anchor_contig: Vec<u32>,
    /// `island_start[i]..island_start[i + 1]` are the probes of the island
    /// after anchor i; length `anchors.len() + 1`.
    pub island_start: Vec<u32>,
    pub probes: Vec<Sampled>,
}

impl AiipsSketch {
    /// Builds the sketch in one streaming pass over FASTA/FASTQ. The summary
    /// (with `valid_kmers` counting probe-length windows) is gathered by the
    /// same scan that samples k-mers.
    pub fn build<R: Read + Send>(input: R, params: AiipsParams) -> Result<Self, BuildError> {
        params.check_buildable()?;
        let mut builder = Builder::new(&params);
        for_each_raw_sequence(input, |sequence| builder.add_record(sequence))
            .map_err(BuildError::Input)?;
        builder.finish(params).map_err(BuildError::Unsupported)
    }

    /// Builds from in-memory records, as `build` would from a file.
    pub fn from_records<'a>(
        records: impl IntoIterator<Item = &'a [u8]>,
        params: AiipsParams,
    ) -> Result<Self, String> {
        params.check_buildable().map_err(|e| e.to_string())?;
        let mut builder = Builder::new(&params);
        for sequence in records {
            builder.add_record(sequence);
        }
        builder.finish(params)
    }

    pub fn is_promoted(&self, anchor: u32) -> bool {
        let contig = &self.contigs[self.anchor_contig[anchor as usize] as usize];
        contig.promoted && contig.first == anchor
    }

    /// Whether probe `probe` is a promoted anchor (its island's first probe,
    /// at the anchor's position).
    fn promoted_probe(&self, probe: u32) -> bool {
        let island = self.island_of_probe(probe);
        self.is_promoted(island)
            && self.island_start[island as usize] == probe
            && self.anchors[island as usize].pos == self.probes[probe as usize].pos
    }

    fn anchor_len(&self, anchor: u32) -> u32 {
        if self.is_promoted(anchor) {
            self.params.probe_k.get() as u32
        } else {
            self.params.anchor_k.get() as u32
        }
    }

    fn island_probes(&self, anchor: u32) -> (u32, u32) {
        (
            self.island_start[anchor as usize],
            self.island_start[anchor as usize + 1],
        )
    }

    /// End (exclusive) of the island after `anchor`: the next anchor's start,
    /// or the record end.
    fn island_end(&self, anchor: u32) -> u32 {
        let contig = &self.contigs[self.anchor_contig[anchor as usize] as usize];
        if anchor + 1 < contig.end {
            self.anchors[anchor as usize + 1].pos
        } else {
            contig.length
        }
    }

    /// Island (anchor index) holding probe `probe`.
    fn island_of_probe(&self, probe: u32) -> u32 {
        (self.island_start.partition_point(|&start| start <= probe) - 1) as u32
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        use std::fmt::Write;
        let p = &self.params;
        let mut out = String::with_capacity(
            300 + 24 * self.contigs.len() + 36 * (self.anchors.len() + self.probes.len()),
        );
        let _ = writeln!(
            out,
            "{} {} algorithm=aiips anchor_k={} anchor_scaled={} probe_k={} probe_scaled={} seed={} hash={} canonical={} contigs={} anchors={} probes={}",
            super::MAGIC,
            super::VERSION,
            p.anchor_k,
            p.anchor_scaled,
            p.probe_k,
            p.probe_scaled,
            p.seed,
            p.hash,
            p.canonical,
            self.contigs.len(),
            self.anchors.len(),
            self.probes.len()
        );
        super::write_summary(&mut out, &self.summary);
        let strand = |s: &Sampled| if s.reverse { '-' } else { '+' };
        for contig in &self.contigs {
            let _ = writeln!(out, "c {}", contig.length);
            for anchor in contig.first..contig.end {
                let a = &self.anchors[anchor as usize];
                let tag = if contig.promoted && anchor == contig.first {
                    'p'
                } else {
                    'a'
                };
                let _ = writeln!(out, "{tag} {} {} {}", a.pos, strand(a), a.hash);
                let (lo, hi) = self.island_probes(anchor);
                for probe in &self.probes[lo as usize..hi as usize] {
                    let _ = writeln!(out, "q {} {} {}", probe.pos, strand(probe), probe.hash);
                }
            }
        }
        out.into_bytes()
    }

    pub(crate) fn parse_body<'a>(
        mut header: Fields<'a>,
        mut lines: impl Iterator<Item = &'a str>,
    ) -> Result<Self, ParseError> {
        let small = |name: &str, value: u64| {
            usize::try_from(value)
                .ok()
                .and_then(NonZeroUsize::new)
                .filter(|k| k.get() <= MAX_PACKED_K)
                .ok_or_else(|| invalid(format!("{name} must be from 1 to {MAX_PACKED_K}")))
        };
        let anchor_k = small("anchor_k", header.take_u64("anchor_k")?)?;
        let anchor_scaled = header.take_u64("anchor_scaled")?;
        let probe_k = small("probe_k", header.take_u64("probe_k")?)?;
        let probe_scaled = header.take_u64("probe_scaled")?;
        check(anchor_k.get(), anchor_scaled, probe_k.get(), probe_scaled)?;
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
        let n_contigs = header.take_u64("contigs")?;
        let n_anchors = header.take_u64("anchors")?;
        let n_probes = header.take_u64("probes")?;
        header.finish()?;
        if n_anchors > u32::MAX as u64 || n_probes > u32::MAX as u64 {
            return Err(invalid("too many anchors or probes"));
        }
        let summary = super::parse_summary(lines.next())?;
        let params = AiipsParams {
            anchor_k,
            anchor_scaled,
            probe_k,
            probe_scaled,
            seed,
            hash,
            canonical,
        };
        let (k1, k2) = (anchor_k.get() as u64, probe_k.get() as u64);
        let (max1, max2) = (params.anchor_max_hash(), params.probe_max_hash());
        // Header counts are untrusted: cap preallocation.
        let cap = |n: u64| n.min(1 << 20) as usize;
        let mut sketch = Self {
            params,
            summary,
            contigs: Vec::with_capacity(cap(n_contigs)),
            anchors: Vec::with_capacity(cap(n_anchors)),
            anchor_contig: Vec::with_capacity(cap(n_anchors)),
            island_start: Vec::with_capacity(cap(n_anchors) + 1),
            probes: Vec::with_capacity(cap(n_probes)),
        };
        // Position after which the next anchor or probe must start.
        let mut anchor_floor = 0u64;
        let mut probe_floor = 0u64;
        for (number, line) in lines.enumerate() {
            let line_no = number + 3;
            let bad = |what: &str| invalid(format!("line {line_no}: {what}: {line:?}"));
            let (tag, rest) = line.split_once(' ').ok_or_else(|| bad("malformed line"))?;
            if tag == "c" {
                let length = parse_u64(rest)
                    .filter(|&l| l <= u32::MAX as u64)
                    .ok_or_else(|| bad("record length is not a 32-bit unsigned integer"))?;
                if sketch.contigs.len() as u64 == n_contigs {
                    return Err(invalid(format!("more than contigs={n_contigs} records")));
                }
                let first = sketch.anchors.len() as u32;
                sketch.contigs.push(Contig {
                    length: length as u32,
                    first,
                    end: first,
                    promoted: false,
                });
                anchor_floor = 0;
                probe_floor = 0;
                continue;
            }
            let mut fields = rest.split(' ');
            let (Some(pos), Some(strand), Some(value), None) =
                (fields.next(), fields.next(), fields.next(), fields.next())
            else {
                return Err(bad("expected position, strand and hash"));
            };
            let pos = parse_u64(pos).ok_or_else(|| bad("position is not an unsigned integer"))?;
            let reverse = match strand {
                "+" => false,
                "-" => true,
                _ => return Err(bad("strand must be + or -")),
            };
            let value = parse_u64(value).ok_or_else(|| bad("hash is not an unsigned integer"))?;
            let contig = sketch
                .contigs
                .last_mut()
                .ok_or_else(|| bad("k-mer before the first record line"))?;
            let length = contig.length as u64;
            let sampled = Sampled {
                pos: pos as u32,
                reverse,
                hash: value,
            };
            match tag {
                "a" | "p" => {
                    let (k, max) = if tag == "a" { (k1, max1) } else { (k2, max2) };
                    if tag == "p" && contig.end != contig.first {
                        return Err(bad("a promoted anchor must be its record's first anchor"));
                    }
                    if value > max {
                        return Err(bad("hash exceeds the sampling threshold"));
                    }
                    if pos < anchor_floor || pos + k > length {
                        return Err(bad(
                            "anchor out of order, overlapping a probe, or past the record end",
                        ));
                    }
                    if sketch.anchors.len() as u64 == n_anchors {
                        return Err(invalid(format!("more than anchors={n_anchors} anchors")));
                    }
                    contig.promoted |= tag == "p";
                    contig.end += 1;
                    let contig_index = (sketch.contigs.len() - 1) as u32;
                    sketch.island_start.push(sketch.probes.len() as u32);
                    sketch.anchors.push(sampled);
                    sketch.anchor_contig.push(contig_index);
                    // Probes must not overlap a k1 anchor. The promoted k2
                    // anchor is also its island's first probe (same line
                    // position), and must not overlap the k1 anchor after it.
                    if tag == "a" {
                        anchor_floor = pos + 1;
                        probe_floor = pos + k1;
                    } else {
                        anchor_floor = pos + k2;
                        probe_floor = pos;
                    }
                }
                "q" => {
                    if contig.end == contig.first {
                        return Err(bad("probe before its record's first anchor"));
                    }
                    if value > max2 {
                        return Err(bad("hash exceeds the sampling threshold"));
                    }
                    if pos < probe_floor || pos + k2 > length {
                        return Err(bad(
                            "probe out of order, overlapping an anchor, or past the record end",
                        ));
                    }
                    if sketch.probes.len() as u64 == n_probes {
                        return Err(invalid(format!("more than probes={n_probes} probes")));
                    }
                    sketch.probes.push(sampled);
                    // The next k1 anchor must not overlap this probe.
                    anchor_floor = anchor_floor.max(pos + k2);
                    probe_floor = pos + 1;
                }
                _ => return Err(bad("unknown line tag")),
            }
        }
        sketch.island_start.push(sketch.probes.len() as u32);
        if sketch.contigs.len() as u64 != n_contigs
            || sketch.anchors.len() as u64 != n_anchors
            || sketch.probes.len() as u64 != n_probes
        {
            return Err(invalid(format!(
                "expected contigs={n_contigs} anchors={n_anchors} probes={n_probes}, found {} {} {}",
                sketch.contigs.len(),
                sketch.anchors.len(),
                sketch.probes.len()
            )));
        }
        Ok(sketch)
    }
}

/// Streams records into the sketch, one pass and two rolling hashes per base.
struct Builder {
    sampler: Sampler,
    contigs: Vec<Contig>,
    anchors: Vec<Sampled>,
    anchor_contig: Vec<u32>,
    island_start: Vec<u32>,
    probes: Vec<Sampled>,
    record_anchors: Vec<Sampled>,
    record_probes: Vec<Sampled>,
    lanes: [LaneOutput; LANES],
    /// The chunk being sampled, after the codes warming it up.
    codes: Vec<u8>,
    block: Box<HitBlock>,
    isa: Isa,
    summary: Summary,
    error: Option<String>,
}

/// Per-record buffers, kept by the thread between sketches: they grow to the
/// longest record, and reusing them spares the second sketch of a comparison
/// faulting in fresh memory.
#[derive(Default)]
struct Scratch {
    record_anchors: Vec<Sampled>,
    record_probes: Vec<Sampled>,
    lanes: [LaneOutput; LANES],
    codes: Vec<u8>,
    block: Box<HitBlock>,
}

thread_local! {
    static SCRATCH: std::cell::RefCell<Scratch> = std::cell::RefCell::default();
}

/// Records are scanned as this many segments side by side, one per vector
/// lane: the rolling codes, hashes and threshold tests of all segments
/// advance together in vector registers.
const LANES: usize = 8;
/// Code for positions outside a lane: breaks k-mers like N but is not a base.
const PAD: u8 = 5;
const ONES: u64 = 0x0101_0101_0101_0101;

/// What one segment sampled, in position order.
#[derive(Default)]
struct LaneOutput {
    anchors: Vec<Sampled>,
    probes: Vec<Sampled>,
}

/// Summary counts gathered by the scan: bases other than A/C/G/T, and
/// probe-length windows of valid bases.
#[derive(Clone, Copy, Default)]
struct Counts {
    invalid: u64,
    windows: u64,
}

/// Hashing and sampling constants.
struct Sampler {
    k1: usize,
    k2: usize,
    m1: u64,
    m2: u64,
    s1: u32,
    s2: u32,
    seed1: u64,
    seed2: u64,
    max1: u64,
    max2: u64,
}

/// 2-bit codes of A, C, G, T in either case; 4 for every other byte. A table
/// lookup avoids a data-dependent branch per base.
static BASE_CODE: [u8; 256] = {
    let mut table = [4u8; 256];
    let mut i = 0;
    while i < 4 {
        table[b"ACGT"[i] as usize] = i as u8;
        table[b"acgt"[i] as usize] = i as u8;
        i += 1;
    }
    table
};

fn mask(k: usize) -> u64 {
    if k == MAX_PACKED_K {
        u64::MAX
    } else {
        (1 << (2 * k)) - 1
    }
}

/// Instruction sets the scanner is compiled for; chosen once at runtime.
/// Every variant runs the same code and gives identical sketches.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[cfg_attr(not(target_arch = "x86_64"), allow(dead_code))]
enum Isa {
    Avx512,
    Avx2,
    Baseline,
}

impl Isa {
    fn detect() -> Self {
        #[cfg(target_arch = "x86_64")]
        {
            if std::arch::is_x86_feature_detected!("avx512f")
                && std::arch::is_x86_feature_detected!("avx512dq")
                && std::arch::is_x86_feature_detected!("avx512vl")
                && std::arch::is_x86_feature_detected!("avx512bw")
            {
                return Isa::Avx512;
            }
            if std::arch::is_x86_feature_detected!("avx2") {
                return Isa::Avx2;
            }
        }
        Isa::Baseline
    }

    /// Variants usable on this machine (for tests).
    #[cfg(test)]
    fn available() -> Vec<Self> {
        let best = Self::detect();
        [Isa::Avx512, Isa::Avx2, Isa::Baseline]
            .into_iter()
            .filter(|&isa| match best {
                Isa::Avx512 => true,
                Isa::Avx2 => isa != Isa::Avx512,
                Isa::Baseline => isa == Isa::Baseline,
            })
            .collect()
    }
}

/// 0x80 in each byte of `word` equal to `byte`, 0 elsewhere (exact: no
/// borrow crosses bytes).
#[inline(always)]
fn bytes_equal(word: u64, byte: u8) -> u64 {
    let low = 0x7f * ONES;
    let z = word ^ (byte as u64 * ONES);
    !(((z & low) + low) | z | low)
}

/// [`BASE_CODE`] of eight bytes at once. A, C, G and T (either case) have
/// bits 1-2 equal to 0, 1, 3 and 2; one Gray-code step maps them to 0-3.
#[inline(always)]
fn encode_word(word: u64) -> u64 {
    let lower = word | (0x20 * ONES);
    let valid = (bytes_equal(lower, b'a')
        | bytes_equal(lower, b'c')
        | bytes_equal(lower, b'g')
        | bytes_equal(lower, b't'))
        >> 7;
    let e = (word >> 1) & (3 * ONES);
    let code = e ^ ((e >> 1) & ONES);
    (code & (valid * 0xff)) | ((valid ^ ONES) << 2)
}

/// Encodes `src` into `dst` (equal lengths), eight bytes per step.
#[inline(always)]
fn encode_into(src: &[u8], dst: &mut [u8]) {
    let (words, tail) = src.as_chunks::<8>();
    let (slots, tail_slots) = dst.as_chunks_mut::<8>();
    for (word, slot) in words.iter().zip(slots) {
        *slot = encode_word(u64::from_le_bytes(*word)).to_le_bytes();
    }
    for (&byte, slot) in tail.iter().zip(tail_slots) {
        *slot = BASE_CODE[byte as usize];
    }
}

/// How a record's line breaks are removed while decoding it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Breaks {
    /// The record has none.
    None,
    /// Only `\n`: lines are cut where the first line's width predicts,
    /// unchecked; the caller compares [`Decoder::newlines`] with the record's
    /// count and decodes again with `Any` if they differ.
    Predicted,
    /// `\n` and `\r` anywhere, found by search.
    Any,
}

/// Streams a record's raw sequence as base codes ([`BASE_CODE`]), dropping
/// the line breaks a FASTA record may hold exactly as needletail's `seq`
/// does.
struct Decoder<'a> {
    rest: &'a [u8],
    mode: Breaks,
    /// Bytes left in the current line.
    line: usize,
    width: Option<usize>,
    /// Line breaks consumed as separators.
    newlines: usize,
}

impl<'a> Decoder<'a> {
    fn new(raw: &'a [u8], mode: Breaks) -> Self {
        Self {
            rest: raw,
            mode,
            line: 0,
            width: None,
            newlines: 0,
        }
    }

    /// Writes the next codes into `dst`; returns how many (0 at the end).
    #[inline(always)]
    fn fill(&mut self, dst: &mut [u8]) -> usize {
        let mut written = 0;
        while written < dst.len() && !self.rest.is_empty() {
            if self.line == 0 {
                match self.mode {
                    Breaks::None => self.line = self.rest.len(),
                    Breaks::Predicted => {
                        if self.rest[0] == b'\n' {
                            self.rest = &self.rest[1..];
                            self.newlines += 1;
                            continue;
                        }
                        self.line = match self.width {
                            Some(w) if self.rest.len() > w && self.rest[w] == b'\n' => w,
                            _ => {
                                let found = self
                                    .rest
                                    .iter()
                                    .position(|&b| b == b'\n')
                                    .unwrap_or(self.rest.len());
                                self.width.get_or_insert(found);
                                found
                            }
                        };
                    }
                    Breaks::Any => {
                        if matches!(self.rest[0], b'\n' | b'\r') {
                            self.rest = &self.rest[1..];
                            continue;
                        }
                        self.line = self
                            .rest
                            .iter()
                            .position(|&b| b == b'\n' || b == b'\r')
                            .unwrap_or(self.rest.len());
                    }
                }
            }
            let take = self.line.min(dst.len() - written);
            encode_into(&self.rest[..take], &mut dst[written..written + take]);
            written += take;
            self.rest = &self.rest[take..];
            self.line -= take;
            if self.line == 0 && self.mode == Breaks::Predicted && !self.rest.is_empty() {
                // The separator found or predicted: a '\n'.
                self.rest = &self.rest[1..];
                self.newlines += 1;
            }
        }
        written
    }
}

/// Transposes an 8x8 byte matrix held as rows: afterwards byte l of `x[j]`
/// is what byte j of `x[l]` was.
#[inline(always)]
fn transpose8(x: &mut [u64; 8]) {
    for l in 0..4 {
        let t = ((x[l] >> 32) ^ x[l + 4]) & 0x0000_0000_ffff_ffff;
        x[l] ^= t << 32;
        x[l + 4] ^= t;
    }
    for l in [0, 1, 4, 5] {
        let t = ((x[l] >> 16) ^ x[l + 2]) & 0x0000_ffff_0000_ffff;
        x[l] ^= t << 16;
        x[l + 2] ^= t;
    }
    for l in [0, 2, 4, 6] {
        let t = ((x[l] >> 8) ^ x[l + 1]) & 0x00ff_00ff_00ff_00ff;
        x[l] ^= t << 8;
        x[l + 1] ^= t;
    }
}

/// Steps buffered before their samples are appended to the lane outputs.
const HIT_BLOCK: usize = 256;

/// Scan steps at which some lane sampled, with the lane state needed to
/// emit the samples: bit l of `anchors` (`probes`) is set when lane l's
/// k1-mer (k2-mer) ending at the step is sampled, of `reverse1`
/// (`reverse2`) when its canonical form is the reverse complement.
struct HitBlock {
    len: usize,
    offset: [u32; HIT_BLOCK],
    anchors: [u8; HIT_BLOCK],
    probes: [u8; HIT_BLOCK],
    reverse1: [u8; HIT_BLOCK],
    reverse2: [u8; HIT_BLOCK],
    h1: [[u64; LANES]; HIT_BLOCK],
    h2: [[u64; LANES]; HIT_BLOCK],
}

impl Default for HitBlock {
    fn default() -> Self {
        Self {
            len: 0,
            offset: [0; HIT_BLOCK],
            anchors: [0; HIT_BLOCK],
            probes: [0; HIT_BLOCK],
            reverse1: [0; HIT_BLOCK],
            reverse2: [0; HIT_BLOCK],
            h1: [[0; LANES]; HIT_BLOCK],
            h2: [[0; LANES]; HIT_BLOCK],
        }
    }
}

impl Sampler {
    /// Samples one chunk of a record. `buf` holds the k1 − 1 codes before
    /// the chunk (`PAD` at the record start) and then the chunk, whose first
    /// base is at record position `base`. The chunk is cut into `LANES`
    /// segments of `seg` bases; each lane first warms up on the k1 − 1 codes
    /// before its segment, so it samples exactly the k-mers ending in its
    /// segment that a sequential scan would. Lanes read eight codes at a
    /// time, transposed so each step holds one code per lane. A step whose
    /// hashes all exceed both thresholds (about nine in ten) costs only the
    /// test; the others are buffered in `block` and emitted in bulk, so the
    /// loop's state stays in registers. Lane outputs are in position order;
    /// concatenated they are the chunk's.
    #[inline(always)]
    fn scan_chunk(
        &self,
        buf: &[u8],
        base: usize,
        block: &mut HitBlock,
        out: &mut [LaneOutput; LANES],
    ) {
        let warm = self.k1 - 1;
        let n = buf.len() - warm;
        let seg = n.div_ceil(LANES);
        let steps = seg + warm;
        // Lane l reads buf[l·seg..ends[l]]: its warm-up, then its segment.
        let ends: [usize; LANES] = std::array::from_fn(|l| warm + ((l + 1) * seg).min(n));
        let codes = buf;
        let (k1, k2) = (self.k1 as u64, self.k2 as u64);
        let mut f1 = [0u64; LANES];
        let mut r1 = [0u64; LANES];
        let mut f2 = [0u64; LANES];
        let mut r2 = [0u64; LANES];
        let mut run = [0u64; LANES];
        block.len = 0;
        let mut group = 0;
        while group < steps {
            let mut rows = [0u64; LANES];
            for l in 0..LANES {
                let start = l * seg + group;
                rows[l] = if start + 8 <= ends[l] {
                    u64::from_le_bytes(codes[start..start + 8].try_into().unwrap())
                } else {
                    let mut bytes = [PAD; 8];
                    for (j, byte) in bytes.iter_mut().enumerate() {
                        if start + j < ends[l] {
                            *byte = codes[start + j];
                        }
                    }
                    u64::from_le_bytes(bytes)
                };
            }
            transpose8(&mut rows);
            // Steps past the last are all PAD: they sample nothing.
            for (j, row) in rows.iter().enumerate() {
                let step = group + j;
                let lane_codes = row.to_le_bytes();
                let mut h1 = [0u64; LANES];
                let mut h2 = [0u64; LANES];
                let mut any = false;
                for l in 0..LANES {
                    let code = lane_codes[l] as u64;
                    let c = code & 3;
                    f1[l] = ((f1[l] << 2) | c) & self.m1;
                    r1[l] = (r1[l] >> 2) | ((3 - c) << self.s1);
                    f2[l] = ((f2[l] << 2) | c) & self.m2;
                    r2[l] = (r2[l] >> 2) | ((3 - c) << self.s2);
                    run[l] = if code > 3 { 0 } else { run[l] + 1 };
                    let c1 = if f1[l] < r1[l] { f1[l] } else { r1[l] };
                    let c2 = if f2[l] < r2[l] { f2[l] } else { r2[l] };
                    h1[l] = fmix64(c1 ^ self.seed1);
                    h2[l] = fmix64(c2 ^ self.seed2);
                    any |= (h1[l] <= self.max1) | (h2[l] <= self.max2);
                }
                if !any || step < warm {
                    continue;
                }
                let (mut anchors, mut probes, mut reverse1, mut reverse2) = (0u8, 0u8, 0u8, 0u8);
                for l in 0..LANES {
                    anchors |= (((run[l] >= k1) & (h1[l] <= self.max1)) as u8) << l;
                    probes |= (((run[l] >= k2) & (h2[l] <= self.max2)) as u8) << l;
                    reverse1 |= ((r1[l] < f1[l]) as u8) << l;
                    reverse2 |= ((r2[l] < f2[l]) as u8) << l;
                }
                if anchors | probes == 0 {
                    continue;
                }
                let i = block.len;
                block.offset[i] = (step - warm) as u32;
                block.anchors[i] = anchors;
                block.probes[i] = probes;
                block.reverse1[i] = reverse1;
                block.reverse2[i] = reverse2;
                block.h1[i] = h1;
                block.h2[i] = h2;
                block.len += 1;
                if block.len == HIT_BLOCK {
                    self.emit(block, base, seg, out);
                }
            }
            group += 8;
        }
        self.emit(block, base, seg, out);
    }

    /// Appends the buffered samples to the lane outputs and empties `block`.
    #[cold]
    #[inline(never)]
    fn emit(&self, block: &mut HitBlock, base: usize, seg: usize, out: &mut [LaneOutput; LANES]) {
        let (k1, k2) = (self.k1 as u32, self.k2 as u32);
        for i in 0..block.len {
            for (l, lane) in out.iter_mut().enumerate() {
                let end = (base + l * seg) as u32 + block.offset[i] + 1;
                if block.probes[i] >> l & 1 == 1 {
                    lane.probes.push(Sampled {
                        pos: end - k2,
                        reverse: block.reverse2[i] >> l & 1 == 1,
                        hash: block.h2[i][l],
                    });
                }
                if block.anchors[i] >> l & 1 == 1 {
                    lane.anchors.push(Sampled {
                        pos: end - k1,
                        reverse: block.reverse1[i] >> l & 1 == 1,
                        hash: block.h1[i][l],
                    });
                }
            }
        }
        block.len = 0;
    }
}

/// Adds a chunk's summary counts: bases other than A/C/G/T, and windows of
/// `k` valid bases. `run` carries the valid bases ending the previous chunk.
/// Without ambiguous bases (the common case) the windows follow from the
/// length.
#[inline(always)]
fn chunk_counts(codes: &[u8], k: usize, run: &mut usize, counts: &mut Counts) {
    let invalid = codes.iter().filter(|&&c| c == 4).count();
    counts.invalid += invalid as u64;
    if invalid == 0 {
        let skipped = (k - 1).saturating_sub(*run).min(codes.len());
        counts.windows += (codes.len() - skipped) as u64;
        *run += codes.len();
        return;
    }
    for &code in codes {
        *run = if code > 3 { 0 } else { *run + 1 };
        counts.windows += (*run >= k) as u64;
    }
}

/// Bases decoded per chunk.
const CHUNK: usize = 1 << 16;

/// Per-record working state of the sampler.
struct RecordBuffers<'b> {
    buf: &'b mut Vec<u8>,
    block: &'b mut HitBlock,
    lanes: &'b mut [LaneOutput; LANES],
    anchors: &'b mut Vec<Sampled>,
    probes: &'b mut Vec<Sampled>,
}

/// Decodes and samples one record chunk by chunk, appending its anchors and
/// candidate probes in position order. Returns the record length and its
/// counts, or only the length when it is too long for 32-bit positions.
#[inline(always)]
fn sample_record(s: &Sampler, raw: &[u8], b: &mut RecordBuffers) -> (usize, Option<Counts>) {
    let (newlines, returns) = raw.iter().fold((0usize, 0usize), |(n, r), &byte| {
        (n + (byte == b'\n') as usize, r + (byte == b'\r') as usize)
    });
    let length = raw.len() - newlines - returns;
    if length > u32::MAX as usize {
        return (length, None);
    }
    let mut mode = match (newlines, returns) {
        (0, 0) => Breaks::None,
        (_, 0) => Breaks::Predicted,
        _ => Breaks::Any,
    };
    let warm = s.k1 - 1;
    loop {
        b.anchors.clear();
        b.probes.clear();
        b.buf.clear();
        b.buf.resize(warm + CHUNK, PAD);
        let mut decoder = Decoder::new(raw, mode);
        let (mut base, mut run, mut counts) = (0, 0, Counts::default());
        loop {
            let len = decoder.fill(&mut b.buf[warm..]);
            if len == 0 {
                break;
            }
            for lane in b.lanes.iter_mut() {
                lane.anchors.clear();
                lane.probes.clear();
            }
            s.scan_chunk(&b.buf[..warm + len], base, b.block, b.lanes);
            for lane in b.lanes.iter() {
                b.anchors.extend_from_slice(&lane.anchors);
                b.probes.extend_from_slice(&lane.probes);
            }
            chunk_counts(&b.buf[warm..warm + len], s.k2, &mut run, &mut counts);
            base += len;
            // The chunk's last k1 − 1 codes warm the next chunk's lanes.
            b.buf.copy_within(len..len + warm, 0);
        }
        if mode == Breaks::Predicted && decoder.newlines != newlines {
            // A predicted line held a line break: decode by search.
            mode = Breaks::Any;
            continue;
        }
        debug_assert_eq!(base, length);
        return (length, Some(counts));
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512dq,avx512vl,avx512bw")]
fn sample_avx512(s: &Sampler, raw: &[u8], b: &mut RecordBuffers) -> (usize, Option<Counts>) {
    sample_record(s, raw, b)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
fn sample_avx2(s: &Sampler, raw: &[u8], b: &mut RecordBuffers) -> (usize, Option<Counts>) {
    sample_record(s, raw, b)
}

fn sample_baseline(s: &Sampler, raw: &[u8], b: &mut RecordBuffers) -> (usize, Option<Counts>) {
    sample_record(s, raw, b)
}

impl Builder {
    fn new(params: &AiipsParams) -> Self {
        let scratch = SCRATCH.with(|s| s.take());
        let (seed1, seed2) = params.hash_seeds();
        let (k1, k2) = (params.anchor_k.get(), params.probe_k.get());
        Self {
            sampler: Sampler {
                k1,
                k2,
                m1: mask(k1),
                m2: mask(k2),
                s1: 2 * (k1 as u32 - 1),
                s2: 2 * (k2 as u32 - 1),
                seed1,
                seed2,
                max1: params.anchor_max_hash(),
                max2: params.probe_max_hash(),
            },
            contigs: Vec::new(),
            anchors: Vec::new(),
            anchor_contig: Vec::new(),
            island_start: Vec::new(),
            probes: Vec::new(),
            record_anchors: scratch.record_anchors,
            record_probes: scratch.record_probes,
            lanes: scratch.lanes,
            codes: scratch.codes,
            block: scratch.block,
            isa: Isa::detect(),
            summary: Summary::default(),
            error: None,
        }
    }

    /// Adds one record given as raw sequence (FASTA line breaks allowed).
    fn add_record(&mut self, raw: &[u8]) {
        self.summary.records += 1;
        if self.error.is_some() {
            return;
        }
        let mut buffers = RecordBuffers {
            buf: &mut self.codes,
            block: &mut self.block,
            lanes: &mut self.lanes,
            anchors: &mut self.record_anchors,
            probes: &mut self.record_probes,
        };
        let sampler = &self.sampler;
        let (length, counts) = match self.isa {
            // SAFETY: these variants are selected only when `Isa::detect`
            // found the instruction sets they are compiled for.
            #[cfg(target_arch = "x86_64")]
            Isa::Avx512 => unsafe { sample_avx512(sampler, raw, &mut buffers) },
            #[cfg(target_arch = "x86_64")]
            Isa::Avx2 => unsafe { sample_avx2(sampler, raw, &mut buffers) },
            _ => sample_baseline(sampler, raw, &mut buffers),
        };
        self.summary.bases += length as u64;
        let Some(counts) = counts else {
            self.error = Some(format!(
                "records longer than {} bases are not supported",
                u32::MAX
            ));
            return;
        };
        self.summary.non_acgt_bases += counts.invalid;
        self.summary.kmers += counts.windows;
        self.store(length as u32);
    }

    /// Drops probes overlapping an anchor, promotes the leading probe, and
    /// appends the record's anchors and islands.
    fn store(&mut self, length: u32) {
        let (k1, k2) = (self.sampler.k1 as u32, self.sampler.k2 as u32);
        let anchors = &self.record_anchors;
        let mut next = 0;
        self.record_probes.retain(|probe| {
            // First anchor ending after the probe starts; later anchors start
            // no earlier, so it alone decides overlap.
            while next < anchors.len() && anchors[next].pos + k1 <= probe.pos {
                next += 1;
            }
            !(next < anchors.len() && anchors[next].pos < probe.pos + k2)
        });
        let first_anchor = self.record_anchors.first().map_or(u32::MAX, |a| a.pos);
        let promoted = self
            .record_probes
            .first()
            .is_some_and(|probe| probe.pos < first_anchor);
        // The promoted probe becomes the record's first anchor and also stays
        // the first probe of the island it starts, so containment counts it
        // like the other genome's copy of the same k2-mer.
        if promoted {
            self.record_anchors.insert(0, self.record_probes[0]);
        }
        let mut probe = 0;
        let contig = self.contigs.len() as u32;
        let first = self.anchors.len() as u32;
        for index in 0..self.record_anchors.len() {
            self.island_start.push(self.probes.len() as u32);
            self.anchors.push(self.record_anchors[index]);
            self.anchor_contig.push(contig);
            let next_pos = self
                .record_anchors
                .get(index + 1)
                .map_or(u32::MAX, |a| a.pos);
            while probe < self.record_probes.len() && self.record_probes[probe].pos < next_pos {
                self.probes.push(self.record_probes[probe]);
                probe += 1;
            }
        }
        debug_assert_eq!(probe, self.record_probes.len());
        if self.probes.len() > u32::MAX as usize || self.anchors.len() > u32::MAX as usize {
            self.error = Some("too many sampled k-mers for one sketch".into());
        }
        self.contigs.push(Contig {
            length,
            first,
            end: self.anchors.len() as u32,
            promoted,
        });
    }

    fn finish(mut self, params: AiipsParams) -> Result<AiipsSketch, String> {
        let summary = self.summary;
        SCRATCH.with(|s| {
            s.replace(Scratch {
                record_anchors: std::mem::take(&mut self.record_anchors),
                record_probes: std::mem::take(&mut self.record_probes),
                lanes: std::mem::take(&mut self.lanes),
                codes: std::mem::take(&mut self.codes),
                block: std::mem::take(&mut self.block),
            })
        });
        if let Some(error) = self.error {
            return Err(error);
        }
        self.island_start.push(self.probes.len() as u32);
        self.anchors.shrink_to_fit();
        self.probes.shrink_to_fit();
        Ok(AiipsSketch {
            params,
            summary,
            contigs: self.contigs,
            anchors: self.anchors,
            anchor_contig: self.anchor_contig,
            island_start: self.island_start,
            probes: self.probes,
        })
    }
}

/// Maximum-likelihood per-base identity q from k-mer match outcomes, each
/// observation `(k, matched, total)` assuming every k-mer of length k matches
/// independently with probability q^k:
/// log L(q) = Σ matched·k·ln q + (total − matched)·ln(1 − q^k).
/// With a single k this is the containment estimate (matched/total)^(1/k).
/// log L is concave in ln q, so the root of its derivative is found by
/// bisection. `None` without observations.
pub fn identity_mle(observations: &[(usize, u64, u64)]) -> Option<f64> {
    let total: u64 = observations.iter().map(|o| o.2).sum();
    let matched: u64 = observations.iter().map(|o| o.1.min(o.2)).sum();
    if total == 0 {
        return None;
    }
    if matched == total {
        return Some(1.0);
    }
    if matched == 0 {
        return Some(0.0);
    }
    // q · d(log L)/dq, strictly decreasing on (0, 1).
    let slope = |q: f64| {
        let ln_q = q.ln();
        observations
            .iter()
            .map(|&(k, n, total)| {
                let (k, n, total) = (k as f64, n.min(total) as f64, total as f64);
                let qk = (k * ln_q).exp();
                let one_minus = -(k * ln_q).exp_m1();
                k * n - (total - n) * k * qk / one_minus
            })
            .sum::<f64>()
    };
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            break;
        }
        if slope(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WalkOptions {
    /// Consecutive low islands tolerated before a walk stops.
    pub max_gap: u32,
    /// An island is high enough when its probe containment in either genome
    /// is at least `min_island_ani^k2`.
    pub min_island_ani: f64,
    /// Count a probe only if it also avoids the other genome's anchors (at
    /// the same walk offset), so both genomes apply one exclusion rule.
    pub symmetric_probes: bool,
    /// Also walk backward from every seed, not only forward.
    pub bidirectional: bool,
}

impl Default for WalkOptions {
    fn default() -> Self {
        Self {
            max_gap: DEFAULT_MAX_GAP,
            min_island_ani: DEFAULT_MIN_ISLAND_ANI,
            symmetric_probes: false,
            bidirectional: false,
        }
    }
}

/// Aligned-region counts for one genome.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SideCounts {
    /// k1 anchors crossed inside walks, and how many matched the other genome.
    pub anchors: u64,
    pub anchors_matched: u64,
    /// Probes in aligned islands, and how many the other genome's paired
    /// islands share.
    pub probes: u64,
    pub probes_matched: u64,
    /// Bases inside aligned islands.
    pub aligned_bases: u64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AiipsEstimates {
    /// Joint maximum-likelihood ANI from both genomes' anchor and probe
    /// counts; `None` when nothing aligned.
    pub ani: Option<f64>,
    /// The same estimator from one genome's counts only.
    pub ani_a: Option<f64>,
    pub ani_b: Option<f64>,
    pub af_a: f64,
    pub af_b: f64,
    pub a: SideCounts,
    pub b: SideCounts,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Origin {
    Anchor(u32),
    Probe(u32),
}

#[derive(Clone, Copy, Debug)]
struct Seed {
    a: Origin,
    b: Origin,
    a_forward: bool,
    b_forward: bool,
    /// Also walk the opposite way (both directions reversed) if this seed
    /// starts a walk: seeds passed by an earlier walk have the stretch before
    /// them covered by that walk already.
    reverse_too: bool,
}

/// A walk's position in one genome: the current (possibly partial) island,
/// its probes, its bases, and the next anchor in the walk direction.
#[derive(Clone, Copy, Debug)]
struct Cursor {
    forward: bool,
    contig: u32,
    lo: u32,
    hi: u32,
    /// Island id: an anchor index for the island after it, or
    /// `anchors.len() + contig` for a record's leading island.
    island: u32,
    probes: (u32, u32),
    span: (u32, u32),
    next: Option<u32>,
    origin_pos: u32,
    origin_len: u32,
}

impl Cursor {
    fn at_anchor(sketch: &AiipsSketch, anchor: u32, forward: bool) -> Self {
        let contig_index = sketch.anchor_contig[anchor as usize];
        let contig = &sketch.contigs[contig_index as usize];
        let mut cursor = Self {
            forward,
            contig: contig_index,
            lo: contig.first,
            hi: contig.end,
            island: 0,
            probes: (0, 0),
            span: (0, 0),
            next: None,
            origin_pos: sketch.anchors[anchor as usize].pos,
            origin_len: sketch.anchor_len(anchor),
        };
        if forward {
            cursor.enter(sketch, Some(anchor), contig_index);
            // A promoted anchor is also its island's first probe: the seed
            // itself is not an observation.
            if sketch.is_promoted(anchor) {
                cursor.probes.0 += 1;
            }
        } else {
            cursor.enter(
                sketch,
                anchor.checked_sub(1).filter(|&u| u >= contig.first),
                contig_index,
            );
        }
        cursor
    }

    fn at_probe(sketch: &AiipsSketch, probe: u32, forward: bool) -> Self {
        let island = sketch.island_of_probe(probe);
        let contig_index = sketch.anchor_contig[island as usize];
        let contig = &sketch.contigs[contig_index as usize];
        let k2 = sketch.params.probe_k.get() as u32;
        let pos = sketch.probes[probe as usize].pos;
        let mut cursor = Self {
            forward,
            contig: contig_index,
            lo: contig.first,
            hi: contig.end,
            island: 0,
            probes: (0, 0),
            span: (0, 0),
            next: None,
            origin_pos: pos,
            origin_len: k2,
        };
        cursor.enter(sketch, Some(island), contig_index);
        if forward {
            cursor.probes.0 = probe + 1;
            cursor.span.0 = pos;
        } else {
            cursor.probes.1 = probe;
            cursor.span.1 = pos + k2;
        }
        cursor
    }

    /// Makes `island` (None: the leading island) current.
    fn enter(&mut self, sketch: &AiipsSketch, island: Option<u32>, contig_index: u32) {
        match island {
            Some(u) => {
                self.island = u;
                self.probes = sketch.island_probes(u);
                self.span = (sketch.anchors[u as usize].pos, sketch.island_end(u));
                self.next = if self.forward {
                    (u + 1 < self.hi).then_some(u + 1)
                } else {
                    Some(u)
                };
            }
            None => {
                let contig = &sketch.contigs[contig_index as usize];
                self.island = sketch.anchors.len() as u32 + contig_index;
                self.probes = (0, 0);
                let end = if contig.first < contig.end {
                    sketch.anchors[contig.first as usize].pos
                } else {
                    contig.length
                };
                self.span = (0, end);
                self.next = (self.forward && contig.first < contig.end).then_some(contig.first);
            }
        }
    }

    /// Distance from the walk origin to anchor `t`'s first base in walk
    /// direction.
    fn rel(&self, sketch: &AiipsSketch, t: u32) -> i64 {
        let pos = sketch.anchors[t as usize].pos as i64;
        if self.forward {
            pos - self.origin_pos as i64
        } else {
            (self.origin_pos + self.origin_len) as i64 - (pos + sketch.anchor_len(t) as i64)
        }
    }

    /// Crosses the next anchor into the following island.
    fn advance(&mut self, sketch: &AiipsSketch) {
        let t = self.next.expect("advance past the record end");
        let contig_index = sketch.anchor_contig[t as usize];
        if self.forward {
            self.enter(sketch, Some(t), contig_index);
        } else {
            self.enter(sketch, (t > self.lo).then(|| t - 1), contig_index);
        }
    }

    /// Offset of the record end in walk direction.
    fn end_rel(&self, sketch: &AiipsSketch) -> i64 {
        if self.forward {
            sketch.contigs[self.contig as usize].length as i64 - self.origin_pos as i64
        } else {
            (self.origin_pos + self.origin_len) as i64
        }
    }

    /// Restricts the current island to walk offsets below `limit`.
    fn clip(&mut self, sketch: &AiipsSketch, limit: i64, k2: u32) {
        let probes = &sketch.probes;
        let clamp = |x: i64| x.clamp(0, u32::MAX as i64) as u32;
        if self.forward {
            let end = clamp(self.origin_pos as i64 + limit);
            while self.probes.1 > self.probes.0 && probes[self.probes.1 as usize - 1].pos + k2 > end
            {
                self.probes.1 -= 1;
            }
            self.span.1 = self.span.1.min(end).max(self.span.0);
        } else {
            let start = clamp((self.origin_pos + self.origin_len) as i64 - limit);
            while self.probes.0 < self.probes.1 && probes[self.probes.0 as usize].pos < start {
                self.probes.0 += 1;
            }
            self.span.0 = self.span.0.max(start).min(self.span.1);
        }
    }

    fn resync(&mut self, sketch: &AiipsSketch, t: u32) {
        self.origin_pos = sketch.anchors[t as usize].pos;
        self.origin_len = sketch.anchor_len(t);
    }
}

#[derive(Clone, Copy, Debug)]
struct Step {
    a_island: u32,
    b_island: u32,
    a_probes: u32,
    b_probes: u32,
    a_span: (u32, u32),
    b_span: (u32, u32),
    shared: u32,
    /// k1 anchors crossed after this step: (anchor, matched).
    a_cross: Option<(u32, bool)>,
    b_cross: Option<(u32, bool)>,
    /// Probe index ranges of the step's (possibly partial) islands and, for
    /// symmetric probes, which of them (bit i: range start + i) are excluded.
    a_range: (u32, u32),
    b_range: (u32, u32),
    a_excluded: u64,
    b_excluded: u64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Status {
    Pass,
    Fail,
    /// No probe on either side: no evidence either way.
    Neutral,
}

/// One walk's contribution, aggregated per island and anchor.
#[derive(Default)]
struct Walks {
    /// (island, shared probes, probes, span) per genome.
    a_islands: Vec<(u32, u32, u32, (u32, u32))>,
    b_islands: Vec<(u32, u32, u32, (u32, u32))>,
    /// (anchor, matched) per genome.
    a_anchors: Vec<(u32, bool)>,
    b_anchors: Vec<(u32, bool)>,
    /// Per walk: score and the end offsets of the four vectors above.
    walks: Vec<(u64, [usize; 4])>,
}

/// The confirmed end of a walk in one genome: the furthest base, in walk
/// direction, of a shared k-mer.
#[derive(Clone, Copy, Debug)]
struct Frontier {
    forward: bool,
    /// Forward: exclusive end; backward: inclusive start.
    limit: u32,
}

impl Frontier {
    fn at_origin(cursor: &Cursor) -> Self {
        Self {
            forward: cursor.forward,
            limit: if cursor.forward {
                cursor.origin_pos + cursor.origin_len
            } else {
                cursor.origin_pos
            },
        }
    }

    /// A frontier that trims nothing.
    fn unbounded(cursor: &Cursor) -> Self {
        Self {
            forward: cursor.forward,
            limit: if cursor.forward { u32::MAX } else { 0 },
        }
    }

    /// Extends the frontier over the k-mer at `pos` of length `len`.
    fn extend(&mut self, pos: u32, len: u32) {
        self.limit = if self.forward {
            self.limit.max(pos + len)
        } else {
            self.limit.min(pos)
        };
    }

    /// Trims a step's island to the frontier. Returns the number of its
    /// non-excluded probes lying wholly inside the frontier when the island
    /// reaches past it, `None` when it lies wholly inside.
    fn trim(
        &self,
        sketch: &AiipsSketch,
        span: &mut (u32, u32),
        range: (u32, u32),
        excluded: u64,
    ) -> Option<u32> {
        let k2 = sketch.params.probe_k.get() as u32;
        if self.forward {
            if span.1 <= self.limit {
                return None;
            }
            span.1 = self.limit.max(span.0);
        } else {
            if span.0 >= self.limit {
                return None;
            }
            span.0 = self.limit.min(span.1);
        }
        let inside = |pos: u32| {
            if self.forward {
                pos + k2 <= self.limit
            } else {
                pos >= self.limit
            }
        };
        Some(
            (range.0..range.1)
                .filter(|&p| {
                    let bit = p - range.0;
                    inside(sketch.probes[p as usize].pos)
                        && (bit >= 64 || (excluded >> bit) & 1 == 0)
                })
                .count() as u32,
        )
    }
}

/// Excluded probes of the cursor's current island as a bit mask over its
/// probe range (probes past the 64th are never marked).
fn exclusion_mask(cursor: &Cursor, excluded: &[u32]) -> u64 {
    excluded
        .iter()
        .filter(|&&p| p >= cursor.probes.0 && p < cursor.probes.1 && p - cursor.probes.0 < 64)
        .fold(0, |mask, &p| mask | 1 << (p - cursor.probes.0))
}

struct Comparer<'s> {
    a: &'s AiipsSketch,
    b: &'s AiipsSketch,
    max_gap: u32,
    threshold: f64,
    symmetric: bool,
    bidirectional: bool,
}

/// Walk offset of probe `probe` in the cursor's walk direction.
fn probe_rel(sketch: &AiipsSketch, cursor: &Cursor, probe: u32) -> i64 {
    let pos = sketch.probes[probe as usize].pos as i64;
    if cursor.forward {
        pos - cursor.origin_pos as i64
    } else {
        (cursor.origin_pos + cursor.origin_len) as i64 - (pos + sketch.params.probe_k.get() as i64)
    }
}

/// k1 anchors bounding the cursor's current island.
fn bounding_anchors(sketch: &AiipsSketch, cursor: &Cursor) -> [Option<u32>; 2] {
    let k1 = |t: u32| (t < cursor.hi && !sketch.is_promoted(t)).then_some(t);
    if cursor.island < sketch.anchors.len() as u32 {
        [k1(cursor.island), k1(cursor.island + 1)]
    } else {
        [
            (cursor.lo < cursor.hi).then_some(cursor.lo).and_then(k1),
            None,
        ]
    }
}

/// Records in `excluded` the unshared probes of `this` genome's current
/// island that overlap a k1 anchor bounding the `other` genome's current
/// island. Offsets of the other genome map to this one by subtracting
/// `shift`.
fn exclude_overlapping(
    this: &AiipsSketch,
    this_cursor: &Cursor,
    other: &AiipsSketch,
    other_cursor: &Cursor,
    other_probes: &[Sampled],
    shift: i64,
    excluded: &mut Vec<u32>,
) {
    let (k1, k2) = (
        this.params.anchor_k.get() as i64,
        this.params.probe_k.get() as i64,
    );
    // The other genome's bounding anchors as offset ranges in this genome's
    // walk frame; overlap is rare, so it is tested first.
    let spans = bounding_anchors(other, other_cursor).map(|t| {
        t.map(|t| {
            let start = other_cursor.rel(other, t) - shift;
            (start - k2, start + k1)
        })
    });
    if spans == [None, None] {
        return;
    }
    for probe in this_cursor.probes.0..this_cursor.probes.1 {
        let start = probe_rel(this, this_cursor, probe);
        // Overlap: start < anchor + k1 and anchor < start + k2.
        if !spans
            .iter()
            .flatten()
            .any(|&(lo, hi)| lo < start && start < hi)
        {
            continue;
        }
        let hash = this.probes[probe as usize].hash;
        if other_probes.iter().any(|q| q.hash == hash) || excluded.contains(&probe) {
            continue;
        }
        excluded.push(probe);
    }
}

/// Probes of the cursor's current island that are not excluded.
fn eligible(cursor: &Cursor, excluded: &[u32]) -> u32 {
    let (lo, hi) = cursor.probes;
    hi - lo - excluded.iter().filter(|&&p| p >= lo && p < hi).count() as u32
}

impl<'s> Comparer<'s> {
    /// Status of a completed island: high enough when the fraction of its
    /// probes shared with the other genome's paired islands reaches the
    /// threshold.
    fn island_status(&self, shared: u32, probes: u32) -> Status {
        if probes == 0 {
            Status::Neutral
        } else if shared > 0 && shared as f64 >= self.threshold * probes as f64 {
            Status::Pass
        } else {
            Status::Fail
        }
    }

    fn cursor(sketch: &AiipsSketch, origin: Origin, forward: bool) -> Cursor {
        match origin {
            Origin::Anchor(t) => Cursor::at_anchor(sketch, t, forward),
            Origin::Probe(p) => Cursor::at_probe(sketch, p, forward),
        }
    }

    /// Walks from `seed`, leaving the retained steps in `steps`. Marks
    /// matched anchor pairs crossed on the way as consumed seeds.
    ///
    /// Each iteration pairs the current island of each genome, then crosses
    /// the next anchor of the genome whose island ends first (both when they
    /// end together or at the same k-mer). An island is judged when it is
    /// completed, from all the shared probes it accumulated; when both
    /// genomes complete an island at once, either passing suffices.
    fn walk(&self, seed: &Seed, steps: &mut Vec<Step>, consumed: &mut Consumed) {
        let (a, b) = (self.a, self.b);
        let k2 = a.params.probe_k.get() as u32;
        let mut ca = Self::cursor(a, seed.a, seed.a_forward);
        let mut cb = Self::cursor(b, seed.b, seed.b_forward);
        steps.clear();
        let mut fails = 0;
        let mut keep: Option<usize> = None;
        // Shared probes accumulated by each genome's current island.
        let (mut sa, mut sb) = (0u32, 0u32);
        // Symmetric probes: probes of each current island overlapping the
        // other genome's anchors, and the offset of B relative to A at the
        // last shared probe (indels shift it between anchors).
        let (mut excluded_a, mut excluded_b) = (Vec::new(), Vec::new());
        let mut shift = 0i64;
        // How far each genome is confirmed homologous (in walk direction) by
        // the seed, shared probes and matched anchors; the walk is trimmed
        // there when it stops.
        let mut hit_a = Frontier::at_origin(&ca);
        let mut hit_b = Frontier::at_origin(&cb);
        let mut keep_hits = (hit_a, hit_b);
        loop {
            // Past a record end nothing is homologous: clip the other
            // genome's island there.
            let (mut ra, mut rb) = (ca, cb);
            if ca.next.is_none() {
                rb.clip(b, ca.end_rel(a) + shift, k2);
            }
            if cb.next.is_none() {
                ra.clip(a, cb.end_rel(b) - shift, k2);
            }
            let pa = &a.probes[ra.probes.0 as usize..ra.probes.1 as usize];
            let pb = &b.probes[rb.probes.0 as usize..rb.probes.1 as usize];
            let mut shared = 0;
            for (i, x) in pa.iter().enumerate() {
                if let Some(j) = pb.iter().position(|y| y.hash == x.hash) {
                    shared += 1;
                    hit_a.extend(x.pos, k2);
                    hit_b.extend(pb[j].pos, k2);
                    // Indels shift B's walk offsets against A's between
                    // re-anchorings; the last shared probe measures it.
                    shift = probe_rel(b, &rb, rb.probes.0 + j as u32)
                        - probe_rel(a, &ra, ra.probes.0 + i as u32);
                }
            }
            sa += shared;
            sb += shared;
            let (mut na_probes, mut nb_probes) = (pa.len() as u32, pb.len() as u32);
            if self.symmetric {
                exclude_overlapping(a, &ra, b, &rb, pb, shift, &mut excluded_a);
                exclude_overlapping(b, &rb, a, &ra, pa, -shift, &mut excluded_b);
                na_probes = eligible(&ra, &excluded_a);
                nb_probes = eligible(&rb, &excluded_b);
            }
            let index = steps.len();
            steps.push(Step {
                a_island: ra.island,
                b_island: rb.island,
                a_probes: na_probes,
                b_probes: nb_probes,
                a_span: ra.span,
                b_span: rb.span,
                shared,
                a_cross: None,
                b_cross: None,
                a_range: ra.probes,
                b_range: rb.probes,
                a_excluded: exclusion_mask(&ra, &excluded_a),
                b_excluded: exclusion_mask(&rb, &excluded_b),
            });
            let next = ca.next.zip(cb.next);
            let (done_a, done_b, matched) = match next {
                Some((na, nb)) => {
                    let promoted = (a.is_promoted(na), b.is_promoted(nb));
                    if a.anchors[na as usize].hash == b.anchors[nb as usize].hash
                        && promoted.0 == promoted.1
                    {
                        (true, true, true)
                    } else {
                        // Island ends compared in A's offsets.
                        let (ea, eb) = (ca.rel(a, na), cb.rel(b, nb) - shift);
                        (ea <= eb, eb <= ea, false)
                    }
                }
                None => (true, true, false),
            };
            let judged = [
                done_a.then(|| self.island_status(sa, na_probes)),
                done_b.then(|| self.island_status(sb, nb_probes)),
            ];
            // The same k-mer bounding both islands confirms homology.
            let event = if matched || judged.contains(&Some(Status::Pass)) {
                Some(Status::Pass)
            } else if judged.contains(&Some(Status::Fail)) {
                Some(Status::Fail)
            } else {
                None
            };
            if matched {
                let (na, nb) = next.unwrap();
                hit_a.extend(a.anchors[na as usize].pos, a.anchor_len(na));
                hit_b.extend(b.anchors[nb as usize].pos, b.anchor_len(nb));
            }
            match event {
                Some(Status::Pass) => {
                    fails = 0;
                    keep = Some(index);
                    keep_hits = (hit_a, hit_b);
                }
                Some(Status::Fail) => fails += 1,
                _ => {}
            }
            let Some((na, nb)) = next else {
                // A record end, not divergence, stopped the walk: with no low
                // island since the last confirmation, homology is taken to
                // run to that end (the islands are already clipped there).
                if fails == 0 && keep.is_some() {
                    keep = Some(index);
                    keep_hits = (Frontier::unbounded(&ca), Frontier::unbounded(&cb));
                }
                break;
            };
            if fails > self.max_gap {
                break;
            }
            if matched {
                // From here the walk equals the one seeded at this pair.
                if !a.is_promoted(na) {
                    steps[index].a_cross = Some((na, true));
                    steps[index].b_cross = Some((nb, true));
                }
                ca.advance(a);
                cb.advance(b);
                ca.resync(a, na);
                cb.resync(b, nb);
                shift = 0;
                consumed.insert((na, nb, ca.forward, cb.forward));
            } else {
                if done_a {
                    if !a.is_promoted(na) {
                        steps[index].a_cross = Some((na, false));
                    }
                    ca.advance(a);
                }
                if done_b {
                    if !b.is_promoted(nb) {
                        steps[index].b_cross = Some((nb, false));
                    }
                    cb.advance(b);
                }
            }
            if done_a {
                sa = 0;
                excluded_a.clear();
            }
            if done_b {
                sb = 0;
                excluded_b.clear();
            }
        }
        match keep {
            Some(last) => {
                steps.truncate(last + 1);
                // Anchors crossed after the last kept step border it, like
                // seeds; only interior crossings are observations.
                steps[last].a_cross = None;
                steps[last].b_cross = None;
                // Homology is confirmed only up to the last shared k-mer:
                // trim the islands reaching past it.
                let (limit_a, limit_b) = keep_hits;
                for step in steps.iter_mut().rev() {
                    let trimmed_a =
                        limit_a.trim(a, &mut step.a_span, step.a_range, step.a_excluded);
                    let trimmed_b =
                        limit_b.trim(b, &mut step.b_span, step.b_range, step.b_excluded);
                    if trimmed_a.is_none() && trimmed_b.is_none() {
                        break;
                    }
                    if let Some(probes) = trimmed_a {
                        step.a_probes = probes;
                    }
                    if let Some(probes) = trimmed_b {
                        step.b_probes = probes;
                    }
                }
            }
            None => steps.clear(),
        }
    }

    fn seeds(&self) -> Vec<Seed> {
        let (a, b) = (self.a, self.b);
        let mut seeds = Vec::new();
        let anchor_keys = |s: &AiipsSketch, promoted: bool| {
            let key = |t: u32| {
                let anchor = &s.anchors[t as usize];
                (anchor.hash, Origin::Anchor(t), anchor.reverse)
            };
            let mut keys: Vec<(u64, Origin, bool)> = if promoted {
                s.contigs
                    .iter()
                    .filter(|c| c.promoted)
                    .map(|c| key(c.first))
                    .collect()
            } else {
                s.contigs
                    .iter()
                    .flat_map(|c| (c.first + c.promoted as u32)..c.end)
                    .map(key)
                    .collect()
            };
            keys.sort_unstable_by_key(|k| k.0);
            keys
        };
        // The sampled k2-mers of `s` (promoted anchors are also probes)
        // whose hash occurs in `wanted`, a sorted list of promoted anchors of
        // the other genome: there are few of those, so this avoids sorting
        // every probe.
        let k2_keys = |s: &AiipsSketch, wanted: &[(u64, Origin, bool)]| {
            let mut keys: Vec<(u64, Origin, bool)> = if wanted.is_empty() {
                Vec::new()
            } else {
                s.probes
                    .iter()
                    .enumerate()
                    .filter(|(_, p)| wanted.binary_search_by_key(&p.hash, |k| k.0).is_ok())
                    .map(|(i, p)| (p.hash, Origin::Probe(i as u32), p.reverse))
                    .collect()
            };
            keys.sort_unstable_by_key(|k| k.0);
            keys
        };
        // k1 anchor pairs: walk forward from the anchor in each genome
        // (one walk when the strands agree); with `bidirectional`, every pair
        // walks both ways.
        let both = self.bidirectional;
        join(&anchor_keys(a, false), &anchor_keys(b, false), |x, y| {
            let same = x.2 == y.2;
            seeds.push(Seed {
                a: x.1,
                b: y.1,
                a_forward: true,
                b_forward: same,
                reverse_too: !same || both,
            });
        });
        // A's promoted anchors against B's k2-mers, walking forward in A.
        let (a_promoted, b_promoted) = (anchor_keys(a, true), anchor_keys(b, true));
        join(&a_promoted, &k2_keys(b, &a_promoted), |x, y| {
            let same = x.2 == y.2;
            seeds.push(Seed {
                a: x.1,
                b: y.1,
                a_forward: true,
                b_forward: same,
                reverse_too: false,
            });
        });
        // B's promoted anchors against A's k2-mers, walking forward in B.
        // (When A's k2-mer is itself a promoted anchor and the strands
        // agree, the walk was seeded above.)
        join(&k2_keys(a, &b_promoted), &b_promoted, |x, y| {
            let same = x.2 == y.2;
            let promoted = matches!(x.1, Origin::Probe(p) if a.promoted_probe(p));
            if !(same && promoted) {
                seeds.push(Seed {
                    a: x.1,
                    b: y.1,
                    a_forward: same,
                    b_forward: true,
                    reverse_too: false,
                });
            }
        });
        seeds.sort_by_key(|seed| seed_order(a, seed));
        seeds
    }

    fn run(&self) -> Walks {
        let mut walks = Walks::default();
        let mut consumed = Consumed::default();
        let mut steps = Vec::new();
        let mut seeds = self.seeds();
        let mut reversed = Vec::new();
        // k1 anchors crossed inside kept walks, per walk orientation
        // (a_forward, b_forward): a seed whose two anchors an earlier walk of
        // its orientation crossed lies on that walk, even when indels kept the
        // pair from being crossed in the same step.
        let mut crossed_a = vec![vec![false; self.a.anchors.len()]; 4];
        let mut crossed_b = vec![vec![false; self.b.anchors.len()]; 4];
        let class = |seed: &Seed| 2 * seed.a_forward as usize + seed.b_forward as usize;
        let mut next = 0;
        while next < seeds.len() || !reversed.is_empty() {
            if next == seeds.len() {
                // Opposite-direction walks of the seeds that started walks,
                // upstream (in their direction) first.
                reversed.sort_by_key(|seed: &Seed| seed_order(self.a, seed));
                seeds.append(&mut reversed);
                continue;
            }
            let seed = seeds[next];
            next += 1;
            if let (Origin::Anchor(x), Origin::Anchor(y)) = (seed.a, seed.b) {
                let c = class(&seed);
                if consumed.contains(&(x, y, seed.a_forward, seed.b_forward))
                    || (crossed_a[c][x as usize] && crossed_b[c][y as usize])
                {
                    continue;
                }
            }
            if seed.reverse_too {
                reversed.push(Seed {
                    a_forward: !seed.a_forward,
                    b_forward: !seed.b_forward,
                    reverse_too: false,
                    ..seed
                });
            }
            self.walk(&seed, &mut steps, &mut consumed);
            let c = class(&seed);
            for step in steps.iter() {
                if let Some((anchor, _)) = step.a_cross {
                    crossed_a[c][anchor as usize] = true;
                }
                if let Some((anchor, _)) = step.b_cross {
                    crossed_b[c][anchor as usize] = true;
                }
            }
            if steps.is_empty() {
                continue;
            }
            let mut score = 0u64;
            let (a_first, b_first) = (walks.a_islands.len(), walks.b_islands.len());
            for step in steps.iter() {
                score += step.shared as u64;
                aggregate(
                    &mut walks.a_islands,
                    a_first,
                    step.a_island,
                    step.shared,
                    step.a_probes,
                    step.a_span,
                );
                aggregate(
                    &mut walks.b_islands,
                    b_first,
                    step.b_island,
                    step.shared,
                    step.b_probes,
                    step.b_span,
                );
                if let Some(cross) = step.a_cross {
                    score += cross.1 as u64;
                    walks.a_anchors.push(cross);
                }
                if let Some(cross) = step.b_cross {
                    walks.b_anchors.push(cross);
                }
            }
            walks.walks.push((
                score,
                [
                    walks.a_islands.len(),
                    walks.b_islands.len(),
                    walks.a_anchors.len(),
                    walks.b_anchors.len(),
                ],
            ));
        }
        walks
    }
}

/// Multiply-rotate hasher for the small integer keys of [`Consumed`]
/// (indices into one comparison's sketches; no untrusted keys).
#[derive(Default)]
struct IndexHasher(u64);

impl Hasher for IndexHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.write_u64(byte as u64);
        }
    }

    fn write_u8(&mut self, i: u8) {
        self.write_u64(i as u64);
    }

    fn write_u32(&mut self, i: u32) {
        self.write_u64(i as u64);
    }

    fn write_u64(&mut self, i: u64) {
        self.0 = (self.0.rotate_left(5) ^ i).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
    }
}

/// Anchor pairs (with walk directions) already walked through.
type Consumed = HashSet<(u32, u32, bool, bool), BuildHasherDefault<IndexHasher>>;

/// Upstream seeds first (in A, along the walk direction), so walks consume
/// the seeds they pass.
fn seed_order(a: &AiipsSketch, seed: &Seed) -> (u8, u64) {
    let key = match seed.a {
        Origin::Anchor(t) => 2 * t as u64,
        Origin::Probe(p) => 2 * a.island_of_probe(p) as u64 + 1,
    };
    if seed.a_forward {
        (0, key)
    } else {
        (1, u64::MAX - key)
    }
}

/// Adds a step's island to the walk record, merging consecutive steps on the
/// same island: its probes count once, shared probes are summed, and the
/// latest extent wins (only the walk's final island can be clipped).
fn aggregate(
    islands: &mut Vec<(u32, u32, u32, (u32, u32))>,
    walk_start: usize,
    island: u32,
    shared: u32,
    probes: u32,
    span: (u32, u32),
) {
    // Only records of the current walk (from `walk_start`) are merged.
    if let Some(last) = islands[walk_start..].last_mut() {
        if last.0 == island {
            last.1 += shared;
            last.2 = probes;
            last.3 = span;
            return;
        }
    }
    islands.push((island, shared, probes, span));
}

/// Calls `pair` for every cross pair of entries with equal hash, skipping
/// hashes that occur more than [`MAX_SEED_OCCURRENCES`] times on either side.
fn join<T: Copy>(
    xs: &[(u64, T, bool)],
    ys: &[(u64, T, bool)],
    mut pair: impl FnMut(&(u64, T, bool), &(u64, T, bool)),
) {
    let (mut i, mut j) = (0, 0);
    while i < xs.len() && j < ys.len() {
        match xs[i].0.cmp(&ys[j].0) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                let hash = xs[i].0;
                let i_end = i + xs[i..].partition_point(|x| x.0 == hash);
                let j_end = j + ys[j..].partition_point(|y| y.0 == hash);
                if i_end - i <= MAX_SEED_OCCURRENCES && j_end - j <= MAX_SEED_OCCURRENCES {
                    for x in &xs[i..i_end] {
                        for y in &ys[j..j_end] {
                            pair(x, y);
                        }
                    }
                }
                i = i_end;
                j = j_end;
            }
        }
    }
}

/// Compares two sketches built with identical parameters.
///
/// 1. Seeds: every anchor pair with equal hash (k1 anchors with k1 anchors;
///    a promoted anchor with any sampled k2-mer of the other genome). The two
///    strands give the relative orientation. A seed walks forward from the
///    anchor in each genome whose anchor it is (one walk when the strands
///    agree, two when they differ); with `bidirectional`, every k1 pair walks
///    both ways.
/// 2. Walk: compare the island after the anchor in one genome with the
///    island after (same orientation) or before (opposite orientation) the
///    matched k-mer in the other, using each genome's own islands. The walk
///    crosses the next anchor on both sides when the next anchors are the
///    same k-mer (re-anchoring), otherwise on the side whose island ends first
///    (from island lengths), so one island can pair with several of the
///    other's. An island is judged when completed, from all the probes it
///    shared: high enough when that fraction reaches `min_island_ani^k2`; when
///    both genomes complete an island at once, either suffices, and matched
///    next anchors always do. The walk stops after more than `max_gap`
///    consecutive low judgements or at a record end (islands are clipped
///    there). It is then cut back to its last high judgement and trimmed to
///    the last shared k-mer, except that a walk reaching a record end with no
///    low island since its last high one keeps its islands up to that end.
/// 3. Overlapping walks (repeats) are resolved best first by score (shared
///    probes plus matched anchors): an island's probe counts come from the
///    first walk claiming it, its aligned bases from the union of the walks'
///    spans; each crossed anchor counts once.
/// 4. ANI: [`identity_mle`] over the k1 anchors crossed inside walks and the
///    probes of aligned islands (k2), pooled over both genomes. Seed anchors
///    are excluded because they were selected for matching. With
///    `symmetric_probes`, probes overlapping the other genome's anchors are
///    left out of the counts, so the exclusion rule is the same on both sides.
/// 5. AF: aligned bases over all bases, per genome.
pub fn estimate(a: &AiipsSketch, b: &AiipsSketch, options: WalkOptions) -> AiipsEstimates {
    let comparer = Comparer {
        a,
        b,
        max_gap: options.max_gap,
        threshold: options.min_island_ani.powi(a.params.probe_k.get() as i32),
        symmetric: options.symmetric_probes,
        bidirectional: options.bidirectional,
    };
    let walks = comparer.run();
    let mut order: Vec<usize> = (0..walks.walks.len()).collect();
    order.sort_by_key(|&w| std::cmp::Reverse(walks.walks[w].0));
    let mut a_side = SideCounts::default();
    let mut b_side = SideCounts::default();
    let mut a_islands = vec![None; a.anchors.len() + a.contigs.len()];
    let mut b_islands = vec![None; b.anchors.len() + b.contigs.len()];
    let mut a_anchors = vec![false; a.anchors.len()];
    let mut b_anchors = vec![false; b.anchors.len()];
    for w in order {
        let start = if w == 0 { [0; 4] } else { walks.walks[w - 1].1 };
        let end = walks.walks[w].1;
        claim_islands(
            &walks.a_islands[start[0]..end[0]],
            &mut a_islands,
            &mut a_side,
        );
        claim_islands(
            &walks.b_islands[start[1]..end[1]],
            &mut b_islands,
            &mut b_side,
        );
        claim_anchors(
            &walks.a_anchors[start[2]..end[2]],
            &mut a_anchors,
            &mut a_side,
        );
        claim_anchors(
            &walks.b_anchors[start[3]..end[3]],
            &mut b_anchors,
            &mut b_side,
        );
    }
    let (k1, k2) = (a.params.anchor_k.get(), a.params.probe_k.get());
    let observations = |sides: &[&SideCounts]| {
        let sum = |f: fn(&SideCounts) -> u64| sides.iter().map(|s| f(s)).sum::<u64>();
        identity_mle(&[
            (k1, sum(|s| s.anchors_matched), sum(|s| s.anchors)),
            (k2, sum(|s| s.probes_matched), sum(|s| s.probes)),
        ])
    };
    let fraction = |aligned: u64, total: u64| {
        if total == 0 {
            0.0
        } else {
            (aligned as f64 / total as f64).min(1.0)
        }
    };
    AiipsEstimates {
        ani: observations(&[&a_side, &b_side]),
        ani_a: observations(&[&a_side]),
        ani_b: observations(&[&b_side]),
        af_a: fraction(a_side.aligned_bases, a.summary.bases),
        af_b: fraction(b_side.aligned_bases, b.summary.bases),
        a: a_side,
        b: b_side,
    }
}

/// Claims a walk's islands. The first walk to reach an island contributes
/// its probe counts; later walks only add the bases they align beyond the
/// span already claimed (an island can straddle a record end or breakpoint
/// of the other genome and be aligned in parts). The claimed extent of an
/// island is kept as one interval, the hull of its claims.
fn claim_islands(
    islands: &[(u32, u32, u32, (u32, u32))],
    claimed: &mut [Option<(u32, u32)>],
    side: &mut SideCounts,
) {
    let length = |span: (u32, u32)| span.1.saturating_sub(span.0) as u64;
    for &(island, shared, probes, span) in islands {
        if span.1 <= span.0 {
            continue;
        }
        match &mut claimed[island as usize] {
            slot @ None => {
                *slot = Some(span);
                side.probes += probes as u64;
                side.probes_matched += shared.min(probes) as u64;
                side.aligned_bases += length(span);
            }
            Some(hull) => {
                let overlap = length((span.0.max(hull.0), span.1.min(hull.1)));
                side.aligned_bases += length(span) - overlap;
                *hull = (hull.0.min(span.0), hull.1.max(span.1));
            }
        }
    }
}

fn claim_anchors(anchors: &[(u32, bool)], claimed: &mut [bool], side: &mut SideCounts) {
    for &(anchor, matched) in anchors {
        if !std::mem::replace(&mut claimed[anchor as usize], true) {
            side.anchors += 1;
            side.anchors_matched += matched as u64;
        }
    }
}

#[cfg(test)]
mod tests;
