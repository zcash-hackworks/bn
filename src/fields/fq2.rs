use fields::{FieldElement, const_fp, Fq};
use std::ops::{Add, Sub, Mul, Neg};
use rand::Rng;

#[inline]
fn non_residue() -> Fq {
    // (q - 1) is a quadratic nonresidue in Fq
    // 21888242871839275222246405745257275088696311157297823662689037894645226208582
    const_fp([317583274, 1757628553, 1923792719, 2366144360, 151523889, 1373741639, 1193918714, 576313009])
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Fq2 {
    c0: Fq,
    c1: Fq
}

impl Fq2 {
    pub fn new(c0: Fq, c1: Fq) -> Self {
        Fq2 {
            c0: c0,
            c1: c1
        }
    }
}

impl FieldElement for Fq2 {
    fn zero() -> Self {
        Fq2 {
            c0: Fq::zero(),
            c1: Fq::zero()
        }
    }

    fn one() -> Self {
        Fq2 {
            c0: Fq::one(),
            c1: Fq::zero()
        }
    }
    
    fn random<R: Rng>(rng: &mut R) -> Self {
        Fq2 {
            c0: Fq::random(rng),
            c1: Fq::random(rng)
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
            c0: (self.c1 * non_residue() + self.c0) * (self.c0 + self.c1) - ab - ab * non_residue(),
            c1: ab + ab
        }
    }

    fn inverse(self) -> Self {
        // "High-Speed Software Implementation of the Optimal Ate Pairing
        // over Barreto–Naehrig Curves"; Algorithm 8

        let t = (self.c0.squared() - (self.c1.squared() * non_residue())).inverse();

        Fq2 {
            c0: self.c0 * t,
            c1: -(self.c1 * t)
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
            c0: bb * non_residue() + aa,
            c1: (self.c0 + self.c1) * (other.c0 + other.c1) - aa - bb
        }
    }
}

impl Sub for Fq2 {
    type Output = Fq2;

    fn sub(self, other: Fq2) -> Fq2 {
        Fq2 {
            c0: self.c0 - other.c0,
            c1: self.c1 - other.c1
        }
    }
}

impl Add for Fq2 {
    type Output = Fq2;

    fn add(self, other: Fq2) -> Fq2 {
        Fq2 {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1
        }
    }
}

impl Neg for Fq2 {
    type Output = Fq2;

    fn neg(self) -> Fq2 {
        Fq2 {
            c0: -self.c0,
            c1: -self.c1
        }
    }
}
