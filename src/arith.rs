use std::cmp::Ordering;
use rand::Rng;

use rustc_serialize::{Encodable, Encoder, Decodable, Decoder};
use byteorder::{ByteOrder, BigEndian};

/// 256-bit, stack allocated biginteger for use in prime field
/// arithmetic.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct U256(pub [u64; 4]);

impl Encodable for U256 {
    fn encode<S: Encoder>(&self, s: &mut S) -> Result<(), S::Error> {
        let mut buf = [0; 32];

        BigEndian::write_u64(&mut buf[0..], self.0[3]);
        BigEndian::write_u64(&mut buf[8..], self.0[2]);
        BigEndian::write_u64(&mut buf[16..], self.0[1]);
        BigEndian::write_u64(&mut buf[24..], self.0[0]);

        for i in 0..32 {
            try!(s.emit_u8(buf[i]));
        }

        Ok(())
    }
}

impl Decodable for U256 {
    fn decode<S: Decoder>(s: &mut S) -> Result<U256, S::Error> {
        let mut buf = [0; 32];

        for i in 0..32 {
            buf[i] = try!(s.read_u8());
        }

        let mut n = [0; 4];
        n[3] = BigEndian::read_u64(&buf[0..]);
        n[2] = BigEndian::read_u64(&buf[8..]);
        n[1] = BigEndian::read_u64(&buf[16..]);
        n[0] = BigEndian::read_u64(&buf[24..]);

        Ok(U256(n))
    }
}

impl Ord for U256 {
    #[inline]
    fn cmp(&self, other: &U256) -> Ordering {
        for (a, b) in self.0.iter().zip(other.0.iter()).rev() {
            if *a < *b {
                return Ordering::Less;
            } else if *a > *b {
                return Ordering::Greater;
            }
        }

        return Ordering::Equal;
    }
}

impl PartialOrd for U256 {
    #[inline]
    fn partial_cmp(&self, other: &U256) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl U256 {
    #[inline]
    pub fn zero() -> U256 {
        U256([0, 0, 0, 0])
    }

    #[inline]
    pub fn one() -> U256 {
        U256([1, 0, 0, 0])
    }

    /// Produce a random number (mod `modulo`)
    pub fn random<R: Rng>(rng: &mut R, modulo: &U256) -> U256
    {
        let mut res;

        loop {
            res = U256(rng.gen());

            for (i, x) in modulo.bits().enumerate() {
                if !x {
                    assert!(res.set_bit(255 - i, false));
                } else {
                    break;
                }
            }

            if &res < modulo { break; }
        }

        res
    }

    pub fn is_zero(&self) -> bool {
        self.0[0] == 0 &&
        self.0[1] == 0 &&
        self.0[2] == 0 &&
        self.0[3] == 0
    }

    pub fn set_bit(&mut self, n: usize, to: bool) -> bool
    {
        if n >= 256 {
            false
        } else {
            let part = n / 64;
            let bit = n - (64 * part);

            if to {
                self.0[part] |= 1 << bit;
            } else {
                self.0[part] &= !(1 << bit);
            }

            true
        }
    }

    pub fn get_bit(&self, n: usize) -> Option<bool>
    {
        if n >= 256 {
            None
        } else {
            let part = n / 64;
            let bit = n - (64 * part);

            Some(self.0[part] & (1 << bit) > 0)
        }
    }

    /// Add `other` to `self` (mod `modulo`)
    pub fn add(&mut self, other: &U256, modulo: &U256) {
        add_nocarry(&mut self.0, &other.0);

        if *self >= *modulo {
            sub_noborrow(&mut self.0, &modulo.0);
        }
    }

    /// Subtract `other` from `self` (mod `modulo`)
    pub fn sub(&mut self, other: &U256, modulo: &U256) {
        if *self < *other {
            add_nocarry(&mut self.0, &modulo.0);
        }

        sub_noborrow(&mut self.0, &other.0);
    }

    /// Multiply `self` by `other` (mod `modulo`) via the Montgomery
    /// multiplication method.
    pub fn mul(&mut self, other: &U256, modulo: &U256, inv: u64) {
        mul_reduce(&mut self.0, &other.0, &modulo.0, inv);

        if *self >= *modulo {
            sub_noborrow(&mut self.0, &modulo.0);
        }
    }

    /// Turn `self` into its additive inverse (mod `modulo`)
    pub fn neg(&mut self, modulo: &U256) {
        if *self > Self::zero() {
            let mut tmp = modulo.0;
            sub_noborrow(&mut tmp, &self.0);
            
            self.0 = tmp;
        }
    }

    #[inline]
    pub fn is_even(&self) -> bool {
        self.0[0] & 1 == 0
    }

    /// Turn `self` into its multiplicative inverse (mod `modulo`)
    pub fn invert(&mut self, modulo: &U256) {
        // Guajardo Kumar Paar Pelzl
        // Efficient Software-Implementation of Finite Fields with Applications to Cryptography
        // Algorithm 16 (BEA for Inversion in Fp)

        let mut u = *self;
        let mut v = *modulo;
        let mut b = U256::one();
        let mut c = U256::zero();

        while u != U256::one() && v != U256::one() {
            while u.is_even() {
                div2(&mut u.0);

                if b.is_even() {
                    div2(&mut b.0);
                } else {
                    add_nocarry(&mut b.0, &modulo.0);
                    div2(&mut b.0);
                }
            }
            while v.is_even() {
                div2(&mut v.0);

                if c.is_even() {
                    div2(&mut c.0);
                } else {
                    add_nocarry(&mut c.0, &modulo.0);
                    div2(&mut c.0);
                }
            }

            if u >= v {
                sub_noborrow(&mut u.0, &v.0);
                b.sub(&c, modulo);
            } else {
                sub_noborrow(&mut v.0, &u.0);
                c.sub(&b, modulo);
            }
        }

        if u == U256::one() {
            self.0 = b.0;
        } else {
            self.0 = c.0;
        }
    }

    /// Return an Iterator<Item=bool> over all bits from
    /// MSB to LSB.
    pub fn bits(&self) -> BitIterator {
        BitIterator {
            int: &self,
            n: 256
        }
    }
}

pub struct BitIterator<'a> {
    int: &'a U256,
    n: usize
}

impl<'a> Iterator for BitIterator<'a> {
    type Item = bool;

    fn next(&mut self) -> Option<bool> {
        if self.n == 0 {
            None
        }
        else {
            self.n -= 1;

            self.int.get_bit(self.n)
        }
    }
}

/// Divide by two
#[inline]
fn div2(a: &mut [u64; 4]) {
    let mut t = a[3] << 63;
    a[3] = a[3] >> 1;
    let b = a[2] << 63;
    a[2] >>= 1;
    a[2] |= t;
    t = a[1] << 63;
    a[1] >>= 1;
    a[1] |= b;
    a[0] >>= 1;
    a[0] |= t;
}

#[inline(always)]
fn split_u64(i: u64) -> (u64, u64) {
    (i >> 32, i & 0xFFFFFFFF)
}

#[inline(always)]
fn combine_u64(hi: u64, lo: u64) -> u64 {
    (hi << 32) | lo
}

#[inline]
fn add_nocarry(a: &mut [u64; 4], b: &[u64; 4]) {
    #[inline]
    fn adc(a: u64, b: u64, carry: &mut u64) -> u64 {
        let (a1, a0) = split_u64(a);
        let (b1, b0) = split_u64(b);
        let (c, r0) = split_u64(a0 + b0 + *carry);
        let (c, r1) = split_u64(a1 + b1 + c);
        *carry = c;

        combine_u64(r1, r0)
    }

    let mut carry = 0;

    for (a, b) in a.into_iter().zip(b.iter()) {
        *a = adc(*a, *b, &mut carry);
    }

    debug_assert!(0 == carry);
}

#[inline]
fn sub_noborrow(a: &mut [u64; 4], b: &[u64; 4]) {
    #[inline]
    fn sbb(a: u64, b: u64, borrow: &mut u64) -> u64 {
        let (a1, a0) = split_u64(a);
        let (b1, b0) = split_u64(b);
        let (b, r0) = split_u64((1 << 32) + a0 - b0 - *borrow);
        let (b, r1) = split_u64((1 << 32) + a1 - b1 - ((b == 0) as u64));

        *borrow = (b == 0) as u64;

        combine_u64(r1, r0)
    }

    let mut borrow = 0;

    for (a, b) in a.into_iter().zip(b.iter()) {
        *a = sbb(*a, *b, &mut borrow);
    }

    debug_assert!(0 == borrow);
}

#[inline]
fn mul_reduce(
    this: &mut [u64; 4],
    by: &[u64; 4],
    modulus: &[u64; 4],
    inv: u64
)
{
    fn mac_digit(acc: &mut [u64], b: &[u64], c: u64)
    {
        #[inline]
        fn mac_with_carry(a: u64, b: u64, c: u64, carry: &mut u64) -> u64 {
            let (b_hi, b_lo) = split_u64(b);
            let (c_hi, c_lo) = split_u64(c);

            let (a_hi, a_lo) = split_u64(a);
            let (carry_hi, carry_lo) = split_u64(*carry);
            let (x_hi, x_lo) = split_u64(b_lo * c_lo + a_lo + carry_lo);
            let (y_hi, y_lo) = split_u64(b_lo * c_hi);
            let (z_hi, z_lo) = split_u64(b_hi * c_lo);
            let (r_hi, r_lo) = split_u64(x_hi + y_lo + z_lo + a_hi + carry_hi);

            *carry = (b_hi * c_hi) + r_hi + y_hi + z_hi;

            combine_u64(r_lo, x_lo)
        }

        if c == 0 {
            return;
        }

        let mut b_iter = b.iter();
        let mut carry = 0;

        for ai in acc.iter_mut() {
            if let Some(bi) = b_iter.next() {
                *ai = mac_with_carry(*ai, *bi, c, &mut carry);
            } else if carry != 0 {
                *ai = mac_with_carry(*ai, 0, c, &mut carry);
            } else {
                break;
            }
        }

        debug_assert!(carry == 0);
    }

    // The Montgomery reduction here is based on Algorithm 14.32 in
    // Handbook of Applied Cryptography
    // <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.

    let mut res = [0; 2*4];
    for (i, xi) in this.iter().enumerate() {
        mac_digit(&mut res[i..], by, *xi);
    }

    for i in 0..4 {
        let k = inv.wrapping_mul(res[i]);
        mac_digit(&mut res[i..], modulus, k);
    }

    this.copy_from_slice(&res[4..]);
}

#[test]
fn setting_bits() {
    let rng = &mut ::rand::thread_rng();
    let modulo = U256([0xffffffffffffffff; 4]);

    let a = U256::random(rng, &modulo);
    let mut e = U256::zero();
    for (i, b) in a.bits().enumerate() {
        assert!(e.set_bit(255 - i, b));
    }

    assert_eq!(a, e);
}
