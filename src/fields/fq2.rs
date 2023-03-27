use fields::{const_fq, FieldElement, Fq};
use rand::Rng;
use std::ops::{Add, Mul, Neg, Sub};

use arith::{U256, U512};

use rustc_serialize::{Decodable, Decoder, Encodable, Encoder};

#[inline]
fn fq_non_residue() -> Fq {
    // (q - 1) is a quadratic nonresidue in Fq
    // 21888242871839275222246405745257275088696311157297823662689037894645226208582
    const_fq([
        0x68c3488912edefaa,
        0x8d087f6872aabf4f,
        0x51e1a24709081231,
        0x2259d6b14729c0fa,
    ])
}

#[inline]
pub fn fq2_nonresidue() -> Fq2 {
    Fq2::new(
        const_fq([
            0xf60647ce410d7ff7,
            0x2f3d6f4dd31bd011,
            0x2943337e3940c6d1,
            0x1d9598e8a7e39857,
        ]),
        const_fq([
            0xd35d438dc58f0d9d,
            0x0a78eb28f5c70b3d,
            0x666ea36f7879462c,
            0x0e0a77c19a07df2f,
        ]),
    )
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct Fq2 {
    c0: Fq,
    c1: Fq,
}

impl Encodable for Fq2 {
    fn encode<S: Encoder>(&self, s: &mut S) -> Result<(), S::Error> {
        let c0: U256 = self.c0.into();
        let c1: U256 = self.c1.into();

        U512::from(&c1, &c0, &Fq::modulus()).encode(s)
    }
}

impl Decodable for Fq2 {
    fn decode<S: Decoder>(s: &mut S) -> Result<Fq2, S::Error> {
        let combined = try!(U512::decode(s));

        match combined.divrem(&Fq::modulus()) {
            (Some(c1), c0) => Ok(Fq2::new(Fq::new(c0).unwrap(), Fq::new(c1).unwrap())),
            _ => Err(s.error("integer not less than modulus squared")),
        }
    }
}

impl Fq2 {
    pub fn new(c0: Fq, c1: Fq) -> Self {
        Fq2 { c0: c0, c1: c1 }
    }

    pub fn scale(&self, by: Fq) -> Self {
        Fq2 {
            c0: self.c0 * by,
            c1: self.c1 * by,
        }
    }

    pub fn mul_by_nonresidue(&self) -> Self {
        *self * fq2_nonresidue()
    }

    pub fn frobenius_map(&self, power: usize) -> Self {
        if power % 2 == 0 {
            *self
        } else {
            Fq2 {
                c0: self.c0,
                c1: self.c1 * fq_non_residue(),
            }
        }
    }
}

impl FieldElement for Fq2 {
    fn zero() -> Self {
        Fq2 {
            c0: Fq::zero(),
            c1: Fq::zero(),
        }
    }

    fn one() -> Self {
        Fq2 {
            c0: Fq::one(),
            c1: Fq::zero(),
        }
    }

    fn random<R: Rng>(rng: &mut R) -> Self {
        Fq2 {
            c0: Fq::random(rng),
            c1: Fq::random(rng),
        }
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn squared(&self) -> Self {
        // Devegili OhEig Scott Dahab
        //     Multiplication and Squaring on Pairing-Friendly Fields.pdf
        //     Section 3 (Complex squaring)

        let ab = self.c0 * self.c1;

        Fq2 {
            c0: (self.c1 * fq_non_residue() + self.c0) * (self.c0 + self.c1)
                - ab
                - ab * fq_non_residue(),
            c1: ab + ab,
        }
    }

    fn inverse(self) -> Option<Self> {
        // "High-Speed Software Implementation of the Optimal Ate Pairing
        // over Barreto–Naehrig Curves"; Algorithm 8

        match (self.c0.squared() - (self.c1.squared() * fq_non_residue())).inverse() {
            Some(t) => Some(Fq2 {
                c0: self.c0 * t,
                c1: -(self.c1 * t),
            }),
            None => None,
        }
    }
}

impl Mul for Fq2 {
    type Output = Fq2;

    fn mul(self, other: Fq2) -> Fq2 {
        // Devegili OhEig Scott Dahab
        //     Multiplication and Squaring on Pairing-Friendly Fields.pdf
        //     Section 3 (Karatsuba)

        let aa = self.c0 * other.c0;
        let bb = self.c1 * other.c1;

        Fq2 {
            c0: bb * fq_non_residue() + aa,
            c1: (self.c0 + self.c1) * (other.c0 + other.c1) - aa - bb,
        }
    }
}

impl Sub for Fq2 {
    type Output = Fq2;

    fn sub(self, other: Fq2) -> Fq2 {
        Fq2 {
            c0: self.c0 - other.c0,
            c1: self.c1 - other.c1,
        }
    }
}

impl Add for Fq2 {
    type Output = Fq2;

    fn add(self, other: Fq2) -> Fq2 {
        Fq2 {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1,
        }
    }
}

impl Neg for Fq2 {
    type Output = Fq2;

    fn neg(self) -> Fq2 {
        Fq2 {
            c0: -self.c0,
            c1: -self.c1,
        }
    }
}
