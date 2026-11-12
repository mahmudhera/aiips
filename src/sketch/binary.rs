//! Binary sketch bodies: tagged sections of canonical LEB128 varints, closed
//! by an FNV-1a 64 checksum. All fixed-width integers are little-endian.
//!
//! ```text
//! body    := section* checksum
//! section := tag (4 ASCII bytes) | payload length (u64) | payload
//! checksum = FNV-1a 64 of every body byte before it (u64)
//! ```
//!
//! Reading is strict: sections must appear in the expected order, payloads
//! must be consumed exactly, varints must use their shortest encoding, and
//! nothing may follow the last section.

use super::{invalid, ParseError};

/// Bytes in the trailing checksum.
const CHECKSUM_BYTES: usize = 8;
/// Bytes in a section header: tag plus payload length.
const SECTION_HEADER_BYTES: usize = 4 + 8;

/// Appends `value` as unsigned LEB128: seven bits per byte, least significant
/// group first, with the high bit set on every byte except the last.
pub(crate) fn put_varint(out: &mut Vec<u8>, mut value: u64) {
    while value >= 0x80 {
        out.push(value as u8 | 0x80);
        value >>= 7;
    }
    out.push(value as u8);
}

/// Why a varint could not be decoded.
#[derive(Debug, PartialEq, Eq)]
pub(crate) enum VarintError {
    Truncated,
    /// A longer encoding than necessary, e.g. `80 00` for 0.
    Overlong,
    Overflow,
}

/// Decodes one canonical LEB128 value from the front of `bytes`, advancing it.
pub(crate) fn get_varint(bytes: &mut &[u8]) -> Result<u64, VarintError> {
    let mut value = 0u64;
    for (index, &byte) in bytes.iter().enumerate() {
        // The tenth byte holds bit 63 only and must end the value.
        if index == 9 && byte > 1 {
            return Err(VarintError::Overflow);
        }
        value |= u64::from(byte & 0x7f) << (7 * index);
        if byte & 0x80 == 0 {
            if byte == 0 && index > 0 {
                return Err(VarintError::Overlong);
            }
            *bytes = &bytes[index + 1..];
            return Ok(value);
        }
    }
    Err(VarintError::Truncated)
}

/// FNV-1a 64 (Fowler, Noll & Vo), used only to detect corrupted files.
pub(crate) fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, &byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

/// Assembles a body section by section.
#[derive(Default)]
pub(crate) struct BodyWriter {
    bytes: Vec<u8>,
}

impl BodyWriter {
    pub(crate) fn section(&mut self, tag: &[u8; 4], payload: &[u8]) {
        self.bytes.reserve(SECTION_HEADER_BYTES + payload.len());
        self.bytes.extend_from_slice(tag);
        self.bytes
            .extend_from_slice(&(payload.len() as u64).to_le_bytes());
        self.bytes.extend_from_slice(payload);
    }

    /// The finished body, checksum appended.
    pub(crate) fn finish(mut self) -> Vec<u8> {
        let checksum = fnv1a64(&self.bytes);
        self.bytes.extend_from_slice(&checksum.to_le_bytes());
        self.bytes
    }
}

/// Hands out the sections of a checksum-verified body in order.
pub(crate) struct BodyReader<'a> {
    rest: &'a [u8],
}

impl<'a> BodyReader<'a> {
    /// Verifies the trailing checksum before any section is read.
    pub(crate) fn new(body: &'a [u8]) -> Result<Self, ParseError> {
        let split = body
            .len()
            .checked_sub(CHECKSUM_BYTES)
            .ok_or_else(|| invalid("binary body is shorter than its checksum"))?;
        let (content, checksum) = body.split_at(split);
        if fnv1a64(content) != u64::from_le_bytes(checksum.try_into().unwrap()) {
            return Err(invalid("checksum mismatch; the file is corrupted"));
        }
        Ok(Self { rest: content })
    }

    /// The next section, which must carry `tag`.
    pub(crate) fn section(&mut self, tag: &'static [u8; 4]) -> Result<Section<'a>, ParseError> {
        let name = tag.escape_ascii();
        if self.rest.len() < SECTION_HEADER_BYTES {
            return Err(invalid(format!("missing section {name}")));
        }
        let (found, rest) = self.rest.split_at(4);
        if found != tag {
            return Err(invalid(format!(
                "expected section {name}, found {}",
                found.escape_ascii()
            )));
        }
        let (length, rest) = rest.split_at(8);
        let length = u64::from_le_bytes(length.try_into().unwrap());
        if length > rest.len() as u64 {
            return Err(invalid(format!(
                "section {name} declares {length} bytes but only {} remain",
                rest.len()
            )));
        }
        let (payload, rest) = rest.split_at(length as usize);
        self.rest = rest;
        Ok(Section { tag, payload })
    }

    /// Succeeds only when every section has been read.
    pub(crate) fn finish(self) -> Result<(), ParseError> {
        if self.rest.is_empty() {
            Ok(())
        } else {
            Err(invalid(format!(
                "{} unexpected bytes after the last section",
                self.rest.len()
            )))
        }
    }
}

/// The payload of one section, read as a sequence of varints.
pub(crate) struct Section<'a> {
    tag: &'static [u8; 4],
    payload: &'a [u8],
}

impl<'a> Section<'a> {
    pub(crate) fn varint(&mut self) -> Result<u64, ParseError> {
        get_varint(&mut self.payload).map_err(|e| {
            let reason = match e {
                VarintError::Truncated => "truncated",
                VarintError::Overlong => "not minimally encoded",
                VarintError::Overflow => "exceeds 64 bits",
            };
            invalid(format!(
                "section {}: varint {reason}",
                self.tag.escape_ascii()
            ))
        })
    }

    /// Unread payload bytes; each remaining varint takes at least one.
    pub(crate) fn remaining(&self) -> usize {
        self.payload.len()
    }

    /// The next `n` raw bytes.
    pub(crate) fn bytes(&mut self, n: usize) -> Result<&'a [u8], ParseError> {
        if n > self.payload.len() {
            return Err(invalid(format!(
                "section {}: needs {n} more bytes, found {}",
                self.tag.escape_ascii(),
                self.payload.len()
            )));
        }
        let (bytes, rest) = self.payload.split_at(n);
        self.payload = rest;
        Ok(bytes)
    }

    /// Succeeds only when the payload has been consumed exactly.
    pub(crate) fn finish(self) -> Result<(), ParseError> {
        if self.payload.is_empty() {
            Ok(())
        } else {
            Err(invalid(format!(
                "section {}: {} unread bytes",
                self.tag.escape_ascii(),
                self.payload.len()
            )))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn encode(value: u64) -> Vec<u8> {
        let mut out = Vec::new();
        put_varint(&mut out, value);
        out
    }

    #[test]
    fn varints_round_trip_with_known_encodings() {
        for (value, bytes) in [
            (0, vec![0x00]),
            (1, vec![0x01]),
            (127, vec![0x7f]),
            (128, vec![0x80, 0x01]),
            (300, vec![0xac, 0x02]),
            (16_383, vec![0xff, 0x7f]),
            (16_384, vec![0x80, 0x80, 0x01]),
            (
                1 << 63,
                vec![0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x01],
            ),
            (
                u64::MAX,
                vec![0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x01],
            ),
        ] {
            assert_eq!(encode(value), bytes, "{value}");
            let mut slice = &bytes[..];
            assert_eq!(get_varint(&mut slice), Ok(value));
            assert!(slice.is_empty());
        }
    }

    #[test]
    fn varint_decoding_advances_past_one_value() {
        let mut bytes: &[u8] = &[0xac, 0x02, 0x05];
        assert_eq!(get_varint(&mut bytes), Ok(300));
        assert_eq!(bytes, &[0x05]);
    }

    #[test]
    fn varints_reject_truncation_overlong_encodings_and_overflow() {
        for (bytes, error) in [
            (vec![], VarintError::Truncated),
            (vec![0x80], VarintError::Truncated),
            (vec![0xff, 0xff], VarintError::Truncated),
            (vec![0x80, 0x00], VarintError::Overlong),
            (vec![0xff, 0x80, 0x00], VarintError::Overlong),
            (
                vec![0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x02],
                VarintError::Overflow,
            ),
            (
                vec![
                    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x81, 0x00,
                ],
                VarintError::Overflow,
            ),
        ] {
            let mut slice = &bytes[..];
            assert_eq!(get_varint(&mut slice), Err(error), "{bytes:02x?}");
        }
    }

    #[test]
    fn fnv1a64_matches_published_test_vectors() {
        // Reference values from the FNV authors' test suite.
        assert_eq!(fnv1a64(b""), 0xcbf2_9ce4_8422_2325);
        assert_eq!(fnv1a64(b"a"), 0xaf63_dc4c_8601_ec8c);
        assert_eq!(fnv1a64(b"foobar"), 0x8594_4171_f739_67e8);
    }

    fn body(sections: &[(&[u8; 4], &[u8])]) -> Vec<u8> {
        let mut writer = BodyWriter::default();
        for (tag, payload) in sections {
            writer.section(tag, payload);
        }
        writer.finish()
    }

    fn error(result: Result<impl Sized, ParseError>) -> String {
        match result {
            Ok(_) => panic!("accepted"),
            Err(e) => e.to_string(),
        }
    }

    #[test]
    fn sections_round_trip_in_order() {
        let bytes = body(&[(b"ONE_", &[0x05]), (b"TWO_", &[]), (b"3RD_", &[0xac, 0x02])]);
        assert_eq!(bytes.len(), 3 * SECTION_HEADER_BYTES + 3 + CHECKSUM_BYTES);
        let mut reader = BodyReader::new(&bytes).unwrap();
        let mut one = reader.section(b"ONE_").unwrap();
        assert_eq!(one.remaining(), 1);
        assert_eq!(one.varint().unwrap(), 5);
        one.finish().unwrap();
        reader.section(b"TWO_").unwrap().finish().unwrap();
        let mut third = reader.section(b"3RD_").unwrap();
        assert_eq!(third.varint().unwrap(), 300);
        third.finish().unwrap();
        reader.finish().unwrap();
    }

    #[test]
    fn body_reader_rejects_corruption_and_structure_errors() {
        let good = body(&[(b"ONE_", &[0x05, 0x06])]);
        // Any flipped bit fails the checksum.
        for index in 0..good.len() {
            let mut bad = good.clone();
            bad[index] ^= 0x10;
            assert!(error(BodyReader::new(&bad)).contains("checksum"), "{index}");
        }
        assert!(error(BodyReader::new(&good[..7])).contains("shorter than its checksum"));

        let mut reader = BodyReader::new(&good).unwrap();
        assert!(error(reader.section(b"TWO_")).contains("expected section TWO_, found ONE_"));

        let mut reader = BodyReader::new(&good).unwrap();
        let mut one = reader.section(b"ONE_").unwrap();
        one.varint().unwrap();
        assert!(error(one.finish()).contains("section ONE_: 1 unread bytes"));
        assert!(error(reader.section(b"TWO_")).contains("missing section TWO_"));

        let reader = BodyReader::new(&good).unwrap();
        assert!(error(reader.finish()).contains("unexpected bytes after the last section"));

        // A declared length running past the body.
        let mut content = b"ONE_".to_vec();
        content.extend_from_slice(&100u64.to_le_bytes());
        content.push(0x01);
        content.extend_from_slice(&fnv1a64(&content).to_le_bytes());
        let mut reader = BodyReader::new(&content).unwrap();
        assert!(error(reader.section(b"ONE_")).contains("declares 100 bytes but only 1 remain"));

        let raw = body(&[(b"RAW_", &[1, 2, 3])]);
        let mut reader = BodyReader::new(&raw).unwrap();
        let mut section = reader.section(b"RAW_").unwrap();
        assert_eq!(section.bytes(2).unwrap(), &[1, 2]);
        assert!(error(section.bytes(2)).contains("section RAW_: needs 2 more bytes, found 1"));

        let bad = body(&[(b"ONE_", &[0x80, 0x00])]);
        let mut reader = BodyReader::new(&bad).unwrap();
        assert!(error(reader.section(b"ONE_").unwrap().varint()).contains("not minimally encoded"));
    }
}
