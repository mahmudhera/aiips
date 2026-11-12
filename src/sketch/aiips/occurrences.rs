//! Columnar storage for sampled k-mer occurrences and the contig table.
//!
//! Occurrences of one kind (anchors or probes) are ordered by contig, then by
//! position. Each holds an index into a dictionary of distinct hashes, a
//! zero-based start position, and a strand bit; probes also carry a role bit.

use std::ops::Range;

/// Packed bits, least significant bit first within each byte. This is both
/// the in-memory form and the file encoding.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Bits {
    bytes: Vec<u8>,
    len: usize,
}

impl Bits {
    pub fn with_capacity(bits: usize) -> Self {
        Self {
            bytes: Vec::with_capacity(bits.div_ceil(8)),
            len: 0,
        }
    }

    #[inline]
    pub fn push(&mut self, bit: bool) {
        if self.len.is_multiple_of(8) {
            self.bytes.push(0);
        }
        if bit {
            *self.bytes.last_mut().unwrap() |= 1 << (self.len % 8);
        }
        self.len += 1;
    }

    #[inline]
    pub fn get(&self, index: usize) -> bool {
        assert!(index < self.len, "bit {index} of {}", self.len);
        self.bytes[index / 8] >> (index % 8) & 1 == 1
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// `len` bits from exactly `ceil(len / 8)` bytes; `None` when the length
    /// differs or a padding bit after the last one is set.
    pub(crate) fn from_bytes(bytes: &[u8], len: usize) -> Option<Self> {
        let padded = bytes.len() == len.div_ceil(8)
            && (len.is_multiple_of(8) || bytes.last().is_some_and(|last| last >> (len % 8) == 0));
        padded.then(|| Self {
            bytes: bytes.to_vec(),
            len,
        })
    }
}

/// Sampled occurrences of one kind, dictionary-encoded.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Occurrences {
    /// Distinct retained hashes, strictly increasing; every one is used by at
    /// least one occurrence.
    pub hashes: Vec<u64>,
    /// Per occurrence: index into `hashes`.
    pub ids: Vec<u32>,
    /// Per occurrence: zero-based start within its contig.
    pub positions: Vec<u32>,
    /// Per occurrence: true when the reverse complement is canonical.
    pub strands: Bits,
}

impl Occurrences {
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    #[inline]
    pub fn hash(&self, index: usize) -> u64 {
        self.hashes[self.ids[index] as usize]
    }

    /// Occurrences of each dictionary entry, indexed like `hashes`.
    pub fn counts(&self) -> Vec<u32> {
        let mut counts = vec![0u32; self.hashes.len()];
        for &id in &self.ids {
            counts[id as usize] += 1;
        }
        counts
    }

    /// Number of dictionary entries that occur exactly once.
    pub fn unique(&self) -> usize {
        self.counts().into_iter().filter(|&n| n == 1).count()
    }
}

/// Collects occurrences with their raw hashes, then dictionary-encodes them.
#[derive(Default)]
pub(crate) struct OccurrenceBuilder {
    hashes: Vec<u64>,
    positions: Vec<u32>,
    strands: Bits,
}

impl OccurrenceBuilder {
    #[inline]
    pub(crate) fn push(&mut self, hash: u64, position: u32, reverse: bool) {
        self.hashes.push(hash);
        self.positions.push(position);
        self.strands.push(reverse);
    }

    pub(crate) fn len(&self) -> usize {
        self.hashes.len()
    }

    /// Sorts and deduplicates the hashes into the dictionary and replaces each
    /// occurrence's hash by its index. `None` when there are more than 2^32
    /// distinct hashes.
    pub(crate) fn finish(self) -> Option<Occurrences> {
        // Sort (hash, occurrence) pairs once, then assign ids in one sweep.
        let count = u32::try_from(self.hashes.len()).ok()?;
        let mut order: Vec<(u64, u32)> = self.hashes.iter().copied().zip(0..count).collect();
        radix_sort(&mut order);
        let mut dictionary: Vec<u64> = Vec::new();
        let mut ids = vec![0u32; order.len()];
        for (hash, index) in order {
            if dictionary.last() != Some(&hash) {
                dictionary.push(hash);
            }
            ids[index as usize] = (dictionary.len() - 1) as u32;
        }
        dictionary.shrink_to_fit();
        Some(Occurrences {
            hashes: dictionary,
            ids,
            positions: self.positions,
            strands: self.strands,
        })
    }
}

/// Stable LSD radix sort by the u64 key, one byte per pass; passes where all
/// keys share the byte are skipped (e.g. the high bytes of hashes below a
/// FracMinHash threshold).
fn radix_sort(items: &mut Vec<(u64, u32)>) {
    let mut buffer = vec![(0u64, 0u32); items.len()];
    for shift in (0..64).step_by(8) {
        let mut counts = [0usize; 256];
        for &(key, _) in items.iter() {
            counts[(key >> shift) as usize & 0xff] += 1;
        }
        if counts.contains(&items.len()) {
            continue;
        }
        let mut offset = 0;
        for count in counts.iter_mut() {
            let n = *count;
            *count = offset;
            offset += n;
        }
        for &item in items.iter() {
            let digit = (item.0 >> shift) as usize & 0xff;
            buffer[counts[digit]] = item;
            counts[digit] += 1;
        }
        std::mem::swap(items, &mut buffer);
    }
}

/// Contigs (input records) in input order, with each contig's half-open
/// ranges into the anchor and probe occurrence arrays.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Contigs {
    pub(crate) lengths: Vec<u32>,
    /// Cumulative ends: contig c owns occurrences `ends[c - 1]..ends[c]`.
    pub(crate) anchor_ends: Vec<usize>,
    pub(crate) probe_ends: Vec<usize>,
}

impl Contigs {
    pub(crate) fn push(&mut self, length: u32, anchor_end: usize, probe_end: usize) {
        self.lengths.push(length);
        self.anchor_ends.push(anchor_end);
        self.probe_ends.push(probe_end);
    }

    pub fn len(&self) -> usize {
        self.lengths.len()
    }

    pub fn is_empty(&self) -> bool {
        self.lengths.is_empty()
    }

    pub fn length(&self, contig: usize) -> u32 {
        self.lengths[contig]
    }

    pub fn anchors(&self, contig: usize) -> Range<usize> {
        span(&self.anchor_ends, contig)
    }

    pub fn probes(&self, contig: usize) -> Range<usize> {
        span(&self.probe_ends, contig)
    }

    /// The contig holding anchor `anchor`.
    pub fn of_anchor(&self, anchor: usize) -> usize {
        self.anchor_ends.partition_point(|&end| end <= anchor)
    }

    /// Valid-or-not windows of length `k` inside the contigs.
    pub(crate) fn windows(&self, k: usize) -> u64 {
        self.lengths
            .iter()
            .map(|&length| (length as u64 + 1).saturating_sub(k as u64))
            .sum()
    }
}

fn span(ends: &[usize], index: usize) -> Range<usize> {
    let start = if index == 0 { 0 } else { ends[index - 1] };
    start..ends[index]
}

/// For every anchor, the index of the first probe at or after its position
/// within the same contig: a linear merge of the two position-ordered arrays.
/// Anchor `a` then owns probes `starts[a]..starts[a + 1]` (or up to the
/// contig's last probe for the contig's last anchor, whose island END closes).
pub(crate) fn link_islands(
    contigs: &Contigs,
    anchor_positions: &[u32],
    probe_positions: &[u32],
) -> Vec<usize> {
    let mut starts = Vec::with_capacity(anchor_positions.len());
    for contig in 0..contigs.len() {
        let probes = contigs.probes(contig);
        let mut probe = probes.start;
        for &anchor in &anchor_positions[contigs.anchors(contig)] {
            while probe < probes.end && probe_positions[probe] < anchor {
                probe += 1;
            }
            starts.push(probe);
        }
    }
    starts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bits_pack_least_significant_first_and_reject_bad_padding() {
        let mut bits = Bits::with_capacity(10);
        for bit in [
            true, false, true, true, false, false, false, false, false, true,
        ] {
            bits.push(bit);
        }
        assert_eq!(bits.as_bytes(), &[0b0000_1101, 0b0000_0010]);
        assert_eq!(bits.len(), 10);
        assert!(bits.get(0) && !bits.get(1) && bits.get(9));
        assert_eq!(Bits::from_bytes(bits.as_bytes(), 10), Some(bits));
        assert_eq!(Bits::from_bytes(&[], 0), Some(Bits::default()));
        assert_eq!(Bits::from_bytes(&[0xff], 8).unwrap().len(), 8);
        for (bytes, len) in [
            (&[0b0000_0100][..], 2),
            (&[0b1000_0000][..], 7),
            (&[0, 0][..], 8),
            (&[][..], 1),
        ] {
            assert_eq!(Bits::from_bytes(bytes, len), None, "{bytes:?} {len}");
        }
    }

    #[test]
    fn dictionary_encoding_preserves_every_occurrence() {
        let mut builder = OccurrenceBuilder::default();
        for (hash, position) in [(9, 0), (4, 3), (9, 5), (7, 8)] {
            builder.push(hash, position, position % 2 == 1);
        }
        let occurrences = builder.finish().unwrap();
        assert_eq!(occurrences.hashes, [4, 7, 9]);
        assert_eq!(occurrences.ids, [2, 0, 2, 1]);
        assert_eq!(
            (0..4).map(|i| occurrences.hash(i)).collect::<Vec<_>>(),
            [9, 4, 9, 7]
        );
        assert_eq!(occurrences.counts(), [1, 1, 2]);
        assert_eq!(occurrences.unique(), 2);
        assert!(occurrences.strands.get(1) && !occurrences.strands.get(0));
    }

    #[test]
    fn radix_sort_orders_by_key_and_keeps_ties_stable() {
        let mut state = 7u64;
        let mut items: Vec<(u64, u32)> = (0..5_000)
            .map(|i| {
                state = crate::hash::fmix64(state.wrapping_add(i));
                // Some shared high bytes, some duplicate keys.
                (if i % 7 == 0 { 42 } else { state >> (i % 3 * 8) }, i as u32)
            })
            .collect();
        let mut expected = items.clone();
        expected.sort_by_key(|&(key, _)| key); // stable
        radix_sort(&mut items);
        assert_eq!(items, expected);
        let mut empty: Vec<(u64, u32)> = Vec::new();
        radix_sort(&mut empty);
        assert!(empty.is_empty());
    }

    #[test]
    fn contig_spans_and_anchor_lookup() {
        let mut contigs = Contigs::default();
        contigs.push(100, 2, 3);
        contigs.push(5, 2, 3); // no occurrences
        contigs.push(50, 5, 4);
        assert_eq!(contigs.anchors(0), 0..2);
        assert_eq!(contigs.anchors(1), 2..2);
        assert_eq!(contigs.anchors(2), 2..5);
        assert_eq!(contigs.probes(2), 3..4);
        assert_eq!(
            (0..5).map(|a| contigs.of_anchor(a)).collect::<Vec<_>>(),
            [0, 0, 2, 2, 2]
        );
        // 91 + 0 (shorter than k) + 41 windows.
        assert_eq!(contigs.windows(10), 132);
    }
}
