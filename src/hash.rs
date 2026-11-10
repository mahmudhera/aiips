//! 64-bit mixing functions used to build hash families.

/// MurmurHash3's 64-bit finalizer (Austin Appleby, public domain). It is a
/// bijection on u64, so distinct inputs never collide.
#[inline]
pub fn fmix64(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51_afd7_ed55_8ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    x ^ (x >> 33)
}

/// SplitMix64 (Steele, Lea & Flood 2014; Vigna's reference implementation),
/// used only to derive per-hash-function seeds from one user seed.
pub struct SplitMix64(u64);

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        Self(seed)
    }
}

impl Iterator for SplitMix64 {
    type Item = u64;

    fn next(&mut self) -> Option<u64> {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        Some(z ^ (z >> 31))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splitmix64_matches_reference_sequence() {
        // Published reference outputs for seed 1234567 (Vigna's splitmix64.c).
        let expected = [
            6457827717110365317,
            3203168211198807973,
            9817491932198370423,
            4593380528125082431,
            16408922859458223821,
        ];
        assert_eq!(
            SplitMix64::new(1234567).take(5).collect::<Vec<_>>(),
            expected
        );
    }

    #[test]
    fn fmix64_fixes_zero_and_inverts_with_the_known_inverse() {
        assert_eq!(fmix64(0), 0);
        // Independent check of the published constants: undo each step with
        // the multiplicative inverses (mod 2^64) and xor-shift inverses.
        fn unxorshift33(x: u64) -> u64 {
            x ^ (x >> 33)
        }
        fn inverse(x: u64) -> u64 {
            let mut x = unxorshift33(x);
            x = x.wrapping_mul(0x9cb4_b2f8_1293_37db);
            x = unxorshift33(x);
            x = x.wrapping_mul(0x4f74_430c_22a5_4005);
            unxorshift33(x)
        }
        for x in SplitMix64::new(7).take(1000).chain([1, u64::MAX, 1 << 63]) {
            assert_eq!(inverse(fmix64(x)), x);
        }
    }
}
