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
use std::io::Read;
use std::num::NonZeroUsize;

use super::minhash::BuildError;
use super::{invalid, parse_u64, Fields, ParseError};
use crate::fastx::{for_each_sequence, Summary};
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
        for_each_sequence(input, |sequence| builder.add_record(sequence))
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
    block: Vec<u8>,
    isa: Isa,
    summary: Summary,
    error: Option<String>,
}

/// Records are scanned as this many segments side by side, one per vector
/// lane: the rolling codes, hashes and threshold tests of all segments
/// advance together, so the compiler can keep them in vector registers.
const LANES: usize = 8;
/// Bases per lane transposed at a time (`BLOCK * LANES` bytes of scratch).
const BLOCK: usize = 2048;
/// Code for positions outside a lane: breaks k-mers like N but is not a base.
const PAD: u8 = 5;

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

impl Sampler {
    /// Samples one record: it is cut into `LANES` segments of `seg` bases,
    /// lane l covering `[l·seg, (l+1)·seg)`. Each lane first warms up on the
    /// k1 − 1 bases before its segment (outside the record: `PAD`), so it
    /// emits exactly the k-mers ending in its segment that a sequential scan
    /// would. Lane outputs are in position order; concatenated they are the
    /// record's.
    #[inline(always)]
    fn scan_lanes(
        &self,
        sequence: &[u8],
        block: &mut [u8],
        out: &mut [LaneOutput; LANES],
    ) -> Counts {
        let n = sequence.len();
        let seg = n.div_ceil(LANES);
        let warm = self.k1 - 1;
        let steps = seg + warm;
        let (k1, k2) = (self.k1 as u64, self.k2 as u64);
        let mut f1 = [0u64; LANES];
        let mut r1 = [0u64; LANES];
        let mut f2 = [0u64; LANES];
        let mut r2 = [0u64; LANES];
        let mut run = [0u64; LANES];
        let mut invalid = [0u64; LANES];
        let mut windows = [0u64; LANES];
        let mut first = 0;
        while first < steps {
            let len = BLOCK.min(steps - first);
            // Transpose: block[j·LANES + l] is lane l's code at step first + j.
            for l in 0..LANES {
                let end = ((l + 1) * seg).min(n);
                for j in 0..len {
                    let pos = (l * seg + first + j).wrapping_sub(warm);
                    block[j * LANES + l] = if pos < end {
                        BASE_CODE[sequence[pos] as usize]
                    } else {
                        PAD
                    };
                }
            }
            for j in 0..len {
                let step = first + j;
                let codes: &[u8; LANES] = block[j * LANES..(j + 1) * LANES].try_into().unwrap();
                let mut h1 = [0u64; LANES];
                let mut h2 = [0u64; LANES];
                let mut hit = false;
                for l in 0..LANES {
                    let code = codes[l] as u64;
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
                    hit |= (h1[l] <= self.max1) | (h2[l] <= self.max2);
                }
                if step < warm {
                    continue;
                }
                for l in 0..LANES {
                    invalid[l] += (codes[l] == 4) as u64;
                    windows[l] += (run[l] >= k2) as u64;
                }
                if hit {
                    for l in 0..LANES {
                        let pos = l * seg + step - warm;
                        if run[l] >= k2 && h2[l] <= self.max2 {
                            out[l].probes.push(Sampled {
                                pos: (pos + 1 - self.k2) as u32,
                                reverse: r2[l] < f2[l],
                                hash: h2[l],
                            });
                        }
                        if run[l] >= k1 && h1[l] <= self.max1 {
                            out[l].anchors.push(Sampled {
                                pos: (pos + 1 - self.k1) as u32,
                                reverse: r1[l] < f1[l],
                                hash: h1[l],
                            });
                        }
                    }
                }
            }
            first += len;
        }
        Counts {
            invalid: invalid.iter().sum(),
            windows: windows.iter().sum(),
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512dq,avx512vl,avx512bw")]
fn scan_avx512(
    s: &Sampler,
    sequence: &[u8],
    block: &mut [u8],
    out: &mut [LaneOutput; LANES],
) -> Counts {
    s.scan_lanes(sequence, block, out)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
fn scan_avx2(
    s: &Sampler,
    sequence: &[u8],
    block: &mut [u8],
    out: &mut [LaneOutput; LANES],
) -> Counts {
    s.scan_lanes(sequence, block, out)
}

fn scan_baseline(
    s: &Sampler,
    sequence: &[u8],
    block: &mut [u8],
    out: &mut [LaneOutput; LANES],
) -> Counts {
    s.scan_lanes(sequence, block, out)
}

impl Builder {
    fn new(params: &AiipsParams) -> Self {
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
            record_anchors: Vec::new(),
            record_probes: Vec::new(),
            lanes: Default::default(),
            block: vec![PAD; BLOCK * LANES],
            isa: Isa::detect(),
            summary: Summary::default(),
            error: None,
        }
    }

    fn add_record(&mut self, sequence: &[u8]) {
        self.summary.records += 1;
        self.summary.bases += sequence.len() as u64;
        if self.error.is_some() {
            return;
        }
        let Ok(length) = u32::try_from(sequence.len()) else {
            self.error = Some(format!(
                "records longer than {} bases are not supported",
                u32::MAX
            ));
            return;
        };
        self.scan(sequence);
        self.store(length);
    }

    /// Collects the record's anchors and candidate probes in position order.
    fn scan(&mut self, sequence: &[u8]) {
        for out in &mut self.lanes {
            out.anchors.clear();
            out.probes.clear();
        }
        let (sampler, block, lanes) = (&self.sampler, &mut self.block[..], &mut self.lanes);
        let counts = match self.isa {
            // SAFETY: these variants are selected only when `Isa::detect`
            // found the instruction sets they are compiled for.
            #[cfg(target_arch = "x86_64")]
            Isa::Avx512 => unsafe { scan_avx512(sampler, sequence, block, lanes) },
            #[cfg(target_arch = "x86_64")]
            Isa::Avx2 => unsafe { scan_avx2(sampler, sequence, block, lanes) },
            _ => scan_baseline(sampler, sequence, block, lanes),
        };
        self.summary.non_acgt_bases += counts.invalid;
        self.summary.kmers += counts.windows;
        self.record_anchors.clear();
        self.record_probes.clear();
        for out in &self.lanes {
            self.record_anchors.extend_from_slice(&out.anchors);
            self.record_probes.extend_from_slice(&out.probes);
        }
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

fn shared_probes(a: &[Sampled], b: &[Sampled]) -> u32 {
    let mut shared = 0;
    for x in a {
        if b.iter().any(|y| y.hash == x.hash) {
            shared += 1;
        }
    }
    shared
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
    let anchors = bounding_anchors(other, other_cursor);
    if anchors == [None, None] {
        return;
    }
    for probe in this_cursor.probes.0..this_cursor.probes.1 {
        let hash = this.probes[probe as usize].hash;
        if other_probes.iter().any(|q| q.hash == hash) || excluded.contains(&probe) {
            continue;
        }
        let start = probe_rel(this, this_cursor, probe);
        let overlaps = anchors.iter().flatten().any(|&t| {
            let anchor = other_cursor.rel(other, t) - shift;
            start < anchor + k1 && anchor < start + k2
        });
        if overlaps {
            excluded.push(probe);
        }
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
    fn walk(
        &self,
        seed: &Seed,
        steps: &mut Vec<Step>,
        consumed: &mut HashSet<(u32, u32, bool, bool)>,
    ) {
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
                rb.clip(b, ca.end_rel(a), k2);
            }
            if cb.next.is_none() {
                ra.clip(a, cb.end_rel(b), k2);
            }
            let pa = &a.probes[ra.probes.0 as usize..ra.probes.1 as usize];
            let pb = &b.probes[rb.probes.0 as usize..rb.probes.1 as usize];
            let shared = shared_probes(pa, pb);
            sa += shared;
            sb += shared;
            if shared > 0 {
                for x in pa {
                    if let Some(y) = pb.iter().find(|y| y.hash == x.hash) {
                        hit_a.extend(x.pos, k2);
                        hit_b.extend(y.pos, k2);
                    }
                }
            }
            let (mut na_probes, mut nb_probes) = (pa.len() as u32, pb.len() as u32);
            if self.symmetric {
                if shared > 0 {
                    // Re-estimate the offset from the last shared probe.
                    for (i, x) in pa.iter().enumerate() {
                        if let Some(j) = pb.iter().position(|y| y.hash == x.hash) {
                            shift = probe_rel(b, &rb, rb.probes.0 + j as u32)
                                - probe_rel(a, &ra, ra.probes.0 + i as u32);
                        }
                    }
                }
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
                        let (ea, eb) = (ca.rel(a, na), cb.rel(b, nb));
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
            let mut keys: Vec<(u64, Origin, bool)> = (0..s.anchors.len() as u32)
                .filter(|&t| s.is_promoted(t) == promoted)
                .map(|t| {
                    (
                        s.anchors[t as usize].hash,
                        Origin::Anchor(t),
                        s.anchors[t as usize].reverse,
                    )
                })
                .collect();
            keys.sort_unstable_by_key(|k| k.0);
            keys
        };
        // Every sampled k2-mer (promoted anchors are also probes).
        let k2_keys = |s: &AiipsSketch| {
            let mut keys: Vec<(u64, Origin, bool)> = s
                .probes
                .iter()
                .enumerate()
                .map(|(i, p)| (p.hash, Origin::Probe(i as u32), p.reverse))
                .collect();
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
        join(&anchor_keys(a, true), &k2_keys(b), |x, y| {
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
        join(&k2_keys(a), &anchor_keys(b, true), |x, y| {
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
        let mut consumed = HashSet::new();
        let mut steps = Vec::new();
        let mut seeds = self.seeds();
        let mut reversed = Vec::new();
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
                if consumed.contains(&(x, y, seed.a_forward, seed.b_forward)) {
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
