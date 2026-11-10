//! Sketch algorithms and the self-describing `.aiips` text format.
//!
//! Line 1 is `aiips-sketch v1 algorithm=NAME key=value...`, so every file
//! records its sketch type and parameters. The remaining lines belong to the
//! algorithm. Parsing is strict: malformed files are rejected, never guessed.

pub mod fracminhash;
pub mod minhash;

use std::num::NonZeroUsize;

use crate::fastx::Summary;

pub const MAGIC: &str = "aiips-sketch";
pub const VERSION: &str = "v1";

/// Every algorithm name the CLI and the format recognize. `MinHash` and
/// `FracMinHash` are implemented; the others are reserved so files and errors
/// stay unambiguous.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Algorithm {
    MinHash,
    FracMinHash,
    BottomK,
    Minimizer,
    OpenSyncmer,
    ClosedSyncmer,
    Strobemer,
}

impl Algorithm {
    pub const ALL: [Algorithm; 7] = [
        Algorithm::MinHash,
        Algorithm::FracMinHash,
        Algorithm::BottomK,
        Algorithm::Minimizer,
        Algorithm::OpenSyncmer,
        Algorithm::ClosedSyncmer,
        Algorithm::Strobemer,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Algorithm::MinHash => "minhash",
            Algorithm::FracMinHash => "fracminhash",
            Algorithm::BottomK => "bottomk",
            Algorithm::Minimizer => "minimizer",
            Algorithm::OpenSyncmer => "open-syncmer",
            Algorithm::ClosedSyncmer => "closed-syncmer",
            Algorithm::Strobemer => "strobemer",
        }
    }

    pub fn from_name(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|a| a.name() == name)
    }

    pub fn is_implemented(self) -> bool {
        matches!(self, Algorithm::MinHash | Algorithm::FracMinHash)
    }

    /// Comma-separated names of the implemented algorithms.
    pub fn implemented_names() -> String {
        Self::ALL
            .into_iter()
            .filter(|a| a.is_implemented())
            .map(Self::name)
            .collect::<Vec<_>>()
            .join(", ")
    }

    pub fn names() -> String {
        Self::ALL.map(Self::name).join(", ")
    }
}

/// A parsed sketch file. Recognized but unimplemented algorithms are kept by
/// name only, so callers can report the type without trusting the body.
#[derive(Debug)]
pub enum Sketch {
    MinHash(minhash::MinHashSketch),
    FracMinHash(fracminhash::FracMinHashSketch),
    Unsupported(Algorithm),
}

impl Sketch {
    pub fn algorithm(&self) -> Algorithm {
        match self {
            Sketch::MinHash(_) => Algorithm::MinHash,
            Sketch::FracMinHash(_) => Algorithm::FracMinHash,
            Sketch::Unsupported(algorithm) => *algorithm,
        }
    }

    /// Input statistics; `None` for unsupported sketches, whose body is not read.
    pub fn summary(&self) -> Option<&Summary> {
        match self {
            Sketch::MinHash(s) => Some(&s.summary),
            Sketch::FracMinHash(s) => Some(&s.summary),
            Sketch::Unsupported(_) => None,
        }
    }

    /// k-mer length; `None` for unsupported sketches.
    pub fn k(&self) -> Option<NonZeroUsize> {
        match self {
            Sketch::MinHash(s) => Some(s.params.k),
            Sketch::FracMinHash(s) => Some(s.params.k),
            Sketch::Unsupported(_) => None,
        }
    }
}

/// Point estimates for one sketch pair. AF is not estimable from a global
/// sketch (positions are discarded), so it is not part of this result.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Estimates {
    pub jaccard: f64,
    pub containment_a_in_b: f64,
    pub containment_b_in_a: f64,
    pub ani: f64,
    pub ani_containment_a_in_b: f64,
    pub ani_containment_b_in_a: f64,
}

/// Bytes reserved for the two header lines when bounding sketch file sizes.
pub const HEADER_BYTES: u64 = 4096;

/// Largest file accepted for a sketch whose file starts with `prefix`, so
/// files are bounded before being read. FracMinHash files are bounded by the
/// `size` their header declares (each value at most 20 digits plus newline);
/// every other file by the largest MinHash sketch. Parsing still validates.
pub fn max_file_bytes(prefix: &[u8]) -> u64 {
    let minhash = 21 * minhash::MAX_NUM_HASHES as u64 + HEADER_BYTES;
    let header = prefix.split(|&b| b == b'\n').next().unwrap_or_default();
    let Ok(header) = std::str::from_utf8(header) else {
        return minhash;
    };
    let mut fields = header.split(' ');
    if fields.next() != Some(MAGIC) || fields.next() != Some(VERSION) {
        return minhash;
    }
    let fields: Vec<&str> = fields.collect();
    if !fields.contains(&"algorithm=fracminhash") {
        return minhash;
    }
    match fields
        .iter()
        .filter_map(|f| f.strip_prefix("size="))
        .collect::<Vec<_>>()[..]
    {
        [size] => parse_u64(size).map_or(minhash, |size| {
            size.saturating_mul(21).saturating_add(HEADER_BYTES)
        }),
        _ => minhash,
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum ParseError {
    /// Zero-byte file, as written by aiips versions before real sketches.
    Empty,
    /// Does not start with the sketch magic, e.g. a FASTA/FASTQ file.
    NotSketch,
    Invalid(String),
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::Empty => write!(
                f,
                "empty file: placeholder from an older aiips version, not a sketch; re-run `aiips sketch`"
            ),
            ParseError::NotSketch => write!(f, "not an aiips sketch; run `aiips sketch` first"),
            ParseError::Invalid(message) => write!(f, "invalid sketch: {message}"),
        }
    }
}

pub(crate) fn invalid(message: impl Into<String>) -> ParseError {
    ParseError::Invalid(message.into())
}

/// Parses a complete sketch file.
pub fn parse(bytes: &[u8]) -> Result<Sketch, ParseError> {
    if bytes.is_empty() {
        return Err(ParseError::Empty);
    }
    if !bytes.starts_with(MAGIC.as_bytes()) {
        return Err(ParseError::NotSketch);
    }
    let text = std::str::from_utf8(bytes).map_err(|_| invalid("not valid UTF-8"))?;
    if text.contains('\r') {
        return Err(invalid("carriage return found; lines must end with \\n"));
    }
    let body = text
        .strip_suffix('\n')
        .ok_or_else(|| invalid("missing final newline"))?;
    let mut lines = body.split('\n');
    let mut header = lines.next().unwrap_or_default().split(' ');
    if header.next() != Some(MAGIC) {
        return Err(invalid(format!("first field must be exactly {MAGIC}")));
    }
    if header.next() != Some(VERSION) {
        return Err(invalid(format!(
            "unsupported format version; expected {VERSION}"
        )));
    }
    let mut params = Fields::parse(header, "header")?;
    let name = params.take("algorithm")?;
    let algorithm = Algorithm::from_name(name)
        .ok_or_else(|| invalid(format!("unknown sketch algorithm: {name}")))?;
    match algorithm {
        Algorithm::MinHash => {
            minhash::MinHashSketch::parse_body(params, lines).map(Sketch::MinHash)
        }
        Algorithm::FracMinHash => {
            fracminhash::FracMinHashSketch::parse_body(params, lines).map(Sketch::FracMinHash)
        }
        other => Ok(Sketch::Unsupported(other)),
    }
}

/// `key=value` fields of one line: no duplicates, no empty keys or values,
/// and every field must be consumed by [`Fields::finish`].
pub(crate) struct Fields<'a> {
    line: &'static str,
    pairs: Vec<(&'a str, &'a str)>,
}

impl<'a> Fields<'a> {
    pub(crate) fn parse(
        tokens: impl Iterator<Item = &'a str>,
        line: &'static str,
    ) -> Result<Self, ParseError> {
        let mut pairs: Vec<(&str, &str)> = Vec::new();
        for token in tokens {
            let (key, value) = token
                .split_once('=')
                .filter(|(k, v)| !k.is_empty() && !v.is_empty())
                .ok_or_else(|| invalid(format!("malformed {line} field: {token:?}")))?;
            if pairs.iter().any(|(k, _)| *k == key) {
                return Err(invalid(format!("duplicate {line} field: {key}")));
            }
            pairs.push((key, value));
        }
        Ok(Self { line, pairs })
    }

    pub(crate) fn take(&mut self, key: &str) -> Result<&'a str, ParseError> {
        let index = self
            .pairs
            .iter()
            .position(|(k, _)| *k == key)
            .ok_or_else(|| invalid(format!("missing {} field: {key}", self.line)))?;
        Ok(self.pairs.remove(index).1)
    }

    pub(crate) fn take_u64(&mut self, key: &str) -> Result<u64, ParseError> {
        let value = self.take(key)?;
        parse_u64(value)
            .ok_or_else(|| invalid(format!("{key} is not an unsigned integer: {value:?}")))
    }

    pub(crate) fn finish(self) -> Result<(), ParseError> {
        match self.pairs.first() {
            Some((key, _)) => Err(invalid(format!("unknown {} field: {key}", self.line))),
            None => Ok(()),
        }
    }
}

/// Strict decimal u64: ASCII digits only (no sign, spaces, or empty string).
pub(crate) fn parse_u64(value: &str) -> Option<u64> {
    if value.is_empty() || !value.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    value.parse().ok()
}

pub(crate) fn write_summary(out: &mut String, summary: &Summary) {
    use std::fmt::Write;
    let _ = writeln!(
        out,
        "records={} bases={} non_acgt_bases={} valid_kmers={}",
        summary.records, summary.bases, summary.non_acgt_bases, summary.kmers
    );
}

pub(crate) fn parse_summary(line: Option<&str>) -> Result<Summary, ParseError> {
    let line = line.ok_or_else(|| invalid("missing summary line"))?;
    let mut fields = Fields::parse(line.split(' '), "summary")?;
    let summary = Summary {
        records: fields.take_u64("records")?,
        bases: fields.take_u64("bases")?,
        non_acgt_bases: fields.take_u64("non_acgt_bases")?,
        kmers: fields.take_u64("valid_kmers")?,
    };
    fields.finish()?;
    if summary.kmers == 0 {
        return Err(invalid("sketch of an input with no valid k-mers"));
    }
    Ok(summary)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn algorithm_names_round_trip_and_minhash_and_fracminhash_are_implemented() {
        for algorithm in Algorithm::ALL {
            assert_eq!(Algorithm::from_name(algorithm.name()), Some(algorithm));
            assert_eq!(
                algorithm.is_implemented(),
                matches!(algorithm, Algorithm::MinHash | Algorithm::FracMinHash)
            );
        }
        assert_eq!(Algorithm::implemented_names(), "minhash, fracminhash");
        assert_eq!(Algorithm::from_name("MinHash"), None);
        assert_eq!(Algorithm::from_name(""), None);
    }

    #[test]
    fn strict_u64_rejects_signs_spaces_and_overflow() {
        assert_eq!(parse_u64("0"), Some(0));
        assert_eq!(parse_u64("18446744073709551615"), Some(u64::MAX));
        for bad in [
            "",
            "+1",
            "-1",
            " 1",
            "1 ",
            "1.0",
            "18446744073709551616",
            "0x1",
        ] {
            assert_eq!(parse_u64(bad), None, "{bad:?}");
        }
    }

    #[test]
    fn classifies_empty_foreign_and_unsupported_files() {
        assert_eq!(parse(b"").unwrap_err(), ParseError::Empty);
        assert_eq!(parse(b">a\nACGT\n").unwrap_err(), ParseError::NotSketch);
        for algorithm in Algorithm::ALL.into_iter().filter(|a| !a.is_implemented()) {
            let file = format!("{MAGIC} {VERSION} algorithm={}\n", algorithm.name());
            assert_eq!(parse(file.as_bytes()).unwrap().algorithm(), algorithm);
        }
    }

    #[test]
    fn file_size_bound_follows_the_declared_fracminhash_size() {
        let minhash = 21 * minhash::MAX_NUM_HASHES as u64 + HEADER_BYTES;
        let frac = |rest: &str| {
            max_file_bytes(
                format!("{MAGIC} {VERSION} algorithm=fracminhash {rest}\nrest").as_bytes(),
            )
        };
        assert_eq!(frac("k=3 size=10"), 210 + HEADER_BYTES);
        assert_eq!(frac("size=0 k=3"), HEADER_BYTES);
        assert_eq!(frac(&format!("size={}", u64::MAX)), u64::MAX);
        for rest in ["k=3", "size=x", "size=1 size=2", "size=+1"] {
            assert_eq!(frac(rest), minhash, "{rest}");
        }
        for prefix in [
            "aiips-sketch v1 algorithm=minhash k=3 size=10\n".as_bytes(),
            b"aiips-sketch v2 algorithm=fracminhash size=10\n",
            b">a\nACGT\n",
            b"",
            b"aiips-sketch v1 algorithm=fracminhash size=1\xff\n",
        ] {
            assert_eq!(max_file_bytes(prefix), minhash, "{prefix:?}");
        }
    }

    #[test]
    fn rejects_malformed_headers() {
        for (file, message) in [
            ("aiips-sketch v2 algorithm=minhash\n", "version"),
            ("aiips-sketch v1\n", "missing header field: algorithm"),
            (
                "aiips-sketch v1 algorithm=sha\n",
                "unknown sketch algorithm",
            ),
            (
                "aiips-sketch v1 algorithm=bottomk algorithm=bottomk\n",
                "duplicate",
            ),
            ("aiips-sketch v1 algorithm=\n", "malformed"),
            ("aiips-sketch v1  algorithm=bottomk\n", "malformed"),
            ("aiips-sketch v1 algorithm=bottomk", "final newline"),
            ("aiips-sketch v1 algorithm=bottomk\r\n", "carriage return"),
        ] {
            let error = parse(file.as_bytes()).unwrap_err().to_string();
            assert!(error.contains(message), "{file:?}: {error}");
        }
        assert!(parse(b"aiips-sketch v1 algorithm=bottomk\xff\n")
            .unwrap_err()
            .to_string()
            .contains("UTF-8"));
    }
}
