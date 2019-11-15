use fields::{const_fq, FieldElement, Fq, Fq2};
use rand::Rng;
use std::ops::{Add, Mul, Neg, Sub};

fn frobenius_coeffs_c1(n: usize) -> Fq2 {
    match n % 6 {
        0 => Fq2::one(),
        1 => Fq2::new(
            const_fq([
                13075984984163199792,
                3782902503040509012,
                8791150885551868305,
                1825854335138010348,
            ]),
            const_fq([
                7963664994991228759,
                12257807996192067905,
                13179524609921305146,
                2767831111890561987,
            ]),
        ),
        2 => Fq2::new(
            const_fq([
                3697675806616062876,
                9065277094688085689,
                6918009208039626314,
                2775033306905974752,
            ]),
            Fq::zero(),
        ),
        3 => Fq2::new(
            const_fq([
                14532872967180610477,
                12903226530429559474,
                1868623743233345524,
                2316889217940299650,
            ]),
            const_fq([
                12447993766991532972,
                4121872836076202828,
                7630813605053367399,
                740282956577754197,
            ]),
        ),
        _ => unimplemented!(),
    }
}
fn frobenius_coeffs_c2(n: usize) -> Fq2 {
    match n % 6 {
        0 => Fq2::one(),
        1 => Fq2::new(
            const_fq([
                8314163329781907090,
                11942187022798819835,
                11282677263046157209,
                1576150870752482284,
            ]),
            const_fq([
                6763840483288992073,
                7118829427391486816,
                4016233444936635065,
                2630958277570195709,
            ]),
        ),
        2 => Fq2::new(
            const_fq([
                8183898218631979349,
                12014359695528440611,
                12263358156045030468,
                3187210487005268291,
            ]),
            Fq::zero(),
        ),
        3 => Fq2::new(
            const_fq([
                4938922280314430175,
                13823286637238282975,
                15589480384090068090,
                481952561930628184,
            ]),
            const_fq([
                3105754162722846417,
                11647802298615474591,
                13057042392041828081,
                1660844386505564338,
            ]),
        ),
        _ => unimplemented!(),
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct Fq6 {
    pub c0: Fq2,
    pub c1: Fq2,
    pub c2: Fq2,
}

impl Fq6 {
    pub fn new(c0: Fq2, c1: Fq2, c2: Fq2) -> Self {
        Fq6 {
            c0: c0,
            c1: c1,
            c2: c2,
        }
    }

    pub fn mul_by_nonresidue(&self) -> Self {
        Fq6 {
            c0: self.c2.mul_by_nonresidue(),
            c1: self.c0,
            c2: self.c1,
        }
    }

    pub fn scale(&self, by: Fq2) -> Self {
        Fq6 {
            c0: self.c0 * by,
            c1: self.c1 * by,
            c2: self.c2 * by,
        }
    }

    pub fn frobenius_map(&self, power: usize) -> Self {
        Fq6 {
            c0: self.c0.frobenius_map(power),
            c1: self.c1.frobenius_map(power) * frobenius_coeffs_c1(power),
            c2: self.c2.frobenius_map(power) * frobenius_coeffs_c2(power),
        }
    }
}

impl FieldElement for Fq6 {
    fn zero() -> Self {
        Fq6 {
            c0: Fq2::zero(),
            c1: Fq2::zero(),
            c2: Fq2::zero(),
        }
    }

    fn one() -> Self {
        Fq6 {
            c0: Fq2::one(),
            c1: Fq2::zero(),
            c2: Fq2::zero(),
        }
    }

    fn random<R: Rng>(rng: &mut R) -> Self {
        Fq6 {
            c0: Fq2::random(rng),
            c1: Fq2::random(rng),
            c2: Fq2::random(rng),
        }
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero() && self.c2.is_zero()
    }

    fn squared(&self) -> Self {
        let s0 = self.c0.squared();
        let ab = self.c0 * self.c1;
        let s1 = ab + ab;
        let s2 = (self.c0 - self.c1 + self.c2).squared();
        let bc = self.c1 * self.c2;
        let s3 = bc + bc;
        let s4 = self.c2.squared();

        Fq6 {
            c0: s0 + s3.mul_by_nonresidue(),
            c1: s1 + s4.mul_by_nonresidue(),
            c2: s1 + s2 + s3 - s0 - s4,
        }
    }

    fn inverse(self) -> Option<Self> {
        let c0 = self.c0.squared() - self.c1 * self.c2.mul_by_nonresidue();
        let c1 = self.c2.squared().mul_by_nonresidue() - self.c0 * self.c1;
        let c2 = self.c1.squared() - self.c0 * self.c2;
        match ((self.c2 * c1 + self.c1 * c2).mul_by_nonresidue() + self.c0 * c0).inverse() {
            Some(t) => Some(Fq6 {
                c0: t * c0,
                c1: t * c1,
                c2: t * c2,
            }),
            None => None,
        }
    }
}

impl Mul for Fq6 {
    type Output = Fq6;

    fn mul(self, other: Fq6) -> Fq6 {
        let a_a = self.c0 * other.c0;
        let b_b = self.c1 * other.c1;
        let c_c = self.c2 * other.c2;

        Fq6 {
            c0: ((self.c1 + self.c2) * (other.c1 + other.c2) - b_b - c_c).mul_by_nonresidue() + a_a,
            c1: (self.c0 + self.c1) * (other.c0 + other.c1) - a_a - b_b + c_c.mul_by_nonresidue(),
            c2: (self.c0 + self.c2) * (other.c0 + other.c2) - a_a + b_b - c_c,
        }
    }
}

impl Sub for Fq6 {
    type Output = Fq6;

    fn sub(self, other: Fq6) -> Fq6 {
        Fq6 {
            c0: self.c0 - other.c0,
            c1: self.c1 - other.c1,
            c2: self.c2 - other.c2,
        }
    }
}

impl Add for Fq6 {
    type Output = Fq6;

    fn add(self, other: Fq6) -> Fq6 {
        Fq6 {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1,
            c2: self.c2 + other.c2,
        }
    }
}

impl Neg for Fq6 {
    type Output = Fq6;

    fn neg(self) -> Fq6 {
        Fq6 {
            c0: -self.c0,
            c1: -self.c1,
            c2: -self.c2,
        }
    }
}
