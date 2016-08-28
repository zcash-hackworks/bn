use std::cmp::Ordering;
use rand::Rng;

mod primitives;
use self::primitives::*;

/// 256-bit, stack allocated biginteger for use in prime field
/// arithmetic.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct U256([Digit; LIMBS]);

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
    fn partial_cmp(&self, other: &U256) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl From<u32> for U256 {
    #[inline]
    fn from(f: u32) -> U256 {
        U256([f, 0, 0, 0, 0, 0, 0, 0])
    }
}

impl From<[u32; 8]> for U256 {
    #[inline]
    fn from(a: [u32; 8]) -> U256 {
        U256(a)
    }
}

impl U256 {
    #[inline]
    pub fn zero() -> U256 {
        0.into()
    }

    #[inline]
    pub fn one() -> U256 {
        1.into()
    }

    /// Produce a random number (mod `modulo`)
    pub fn rand<R: Rng>(
        rng: &mut R,
        modulo: &U256
    ) -> U256
    {
        let mut res;

        loop {
            res = U256(rng.gen());

            for (i, x) in modulo.bits().enumerate() {
                if !x {
                    res.set_bit(255 - i, false);
                } else {
                    break;
                }
            }

            if &res < modulo { break; }
        }

        res
    }

    pub fn is_zero(&self) -> bool {
        for a in self.0.iter() {
            if *a != 0 {
                return false;
            }
        }

        return true;
    }

    pub fn set_bit(&mut self, n: usize, to: bool)
    {
        assert!(n < 256);

        let part = n / 32;
        let bit = n - (32 * part);

        if to {
            self.0[part] |= 1 << bit;
        } else {
            self.0[part] &= !(1 << bit);
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
    pub fn mul(&mut self, other: &U256, modulo: &U256, inv: u32) {
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

            let part = self.n / 32;
            let bit = self.n - (32 * part);

            Some(self.int.0[part] & (1 << bit) > 0)
        }
    }
}

#[test]
fn setting_bits() {
    let rng = &mut ::rand::thread_rng();
    let modulo = U256([0xffffffff; LIMBS]);

    let a = U256::rand(rng, &modulo);
    let mut e = U256::zero();
    for (i, b) in a.bits().enumerate() {
        e.set_bit(255 - i, b);
    }

    assert_eq!(a, e);
}
