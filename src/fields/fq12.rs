use fields::{FieldElement, Fq2, Fq, Fq6, const_fq};
use std::ops::{Add, Sub, Mul, Neg};
use rand::Rng;

use arith::U256;

fn frobenius_coeffs_c1(power: usize) -> Fq2 {
    match power % 12 {
        0 => Fq2::one(),
        1 => Fq2::new(
            const_fq([12653890742059813127, 14585784200204367754, 1278438861261381767, 212598772761311868]),
            const_fq([11683091849979440498, 14992204589386555739, 15866167890766973222, 1200023580730561873])
        ),
        2 => Fq2::new(
            const_fq([14595462726357228530, 17349508522658994025, 1017833795229664280, 299787779797702374]),
            Fq::zero()
        ),
        3 => Fq2::new(
            const_fq([3914496794763385213, 790120733010914719, 7322192392869644725, 581366264293887267]),
            const_fq([12817045492518885689, 4440270538777280383, 11178533038884588256, 2767537931541304486])
        ),
        _ => unimplemented!()
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct Fq12 {
    c0: Fq6,
    c1: Fq6
}

impl Fq12 {
    pub fn new(c0: Fq6, c1: Fq6) -> Self {
        Fq12 {
            c0: c0,
            c1: c1
        }
    }

    fn final_exponentiation_first_chunk(&self) -> Option<Fq12> {
        match self.inverse() {
            Some(b) => {
                let a = self.unitary_inverse();
                let c = a * b;
                let d = c.frobenius_map(2);

                Some(d * c)
            },
            None => None
        }
    }

    fn final_exponentiation_last_chunk(&self) -> Fq12 {
        let a = self.exp_by_neg_z();
        let b = a.cyclotomic_squared();
        let c = b.cyclotomic_squared();
        let d = c * b;

        let e = d.exp_by_neg_z();
        let f = e.cyclotomic_squared();
        let g = f.exp_by_neg_z();
        let h = d.unitary_inverse();
        let i = g.unitary_inverse();

        let j = i * e;
        let k = j * h;
        let l = k * b;
        let m = k * e;
        let n = *self * m;

        let o = l.frobenius_map(1);
        let p = o * n;

        let q = k.frobenius_map(2);
        let r = q * p;

        let s = self.unitary_inverse();
        let t = s * l;
        let u = t.frobenius_map(3);
        let v = u * r;

        v
    }

    pub fn final_exponentiation(&self) -> Option<Fq12> {
        self.final_exponentiation_first_chunk().map(|a| a.final_exponentiation_last_chunk())
    }

    pub fn frobenius_map(&self, power: usize) -> Self {
        Fq12 {
            c0: self.c0.frobenius_map(power),
            c1: self.c1.frobenius_map(power).scale(frobenius_coeffs_c1(power))
        }
    }

    pub fn exp_by_neg_z(&self) -> Fq12 {
        self.cyclotomic_pow(
            U256([4965661367192848881, 0, 0, 0])
        ).unitary_inverse()
    }

    pub fn unitary_inverse(&self) -> Fq12 {
        Fq12::new(self.c0, -self.c1)
    }

    pub fn mul_by_024(&self,
                      ell_0: Fq2,
                      ell_vw: Fq2,
                      ell_vv: Fq2) -> Fq12 {
        let z0 = self.c0.c0;
        let z1 = self.c0.c1;
        let z2 = self.c0.c2;
        let z3 = self.c1.c0;
        let z4 = self.c1.c1;
        let z5 = self.c1.c2;

        let x0 = ell_0;
        let x2 = ell_vv;
        let x4 = ell_vw;

        let d0 = z0 * x0;
        let d2 = z2 * x2;
        let d4 = z4 * x4;
        let t2 = z0 + z4;
        let t1 = z0 + z2;
        let s0 = z1 + z3 + z5;

        let s1 = z1 * x2;
        let t3 = s1 + d4;
        let t4 = t3.mul_by_nonresidue() + d0;
        let z0 = t4;

        let t3 = z5 * x4;
        let s1 = s1 + t3;
        let t3 = t3 + d2;
        let t4 = t3.mul_by_nonresidue();
        let t3 = z1 * x0;
        let s1 = s1 + t3;
        let t4 = t4 + t3;
        let z1 = t4;

        let t0 = x0 + x2;
        let t3 = t1 * t0 - d0 - d2;
        let t4 = z3 * x4;
        let s1 = s1 + t4;
        let t3 = t3 + t4;

        let t0 = z2 + z4;
        let z2 = t3;

        let t1 = x2 + x4;
        let t3 = t0 * t1 - d2 - d4;
        let t4 = t3.mul_by_nonresidue();
        let t3 = z3 * x0;
        let s1 = s1 + t3;
        let t4 = t4 + t3;
        let z3 = t4;

        let t3 = z5 * x2;
        let s1 = s1 + t3;
        let t4 = t3.mul_by_nonresidue();
        let t0 = x0 + x4;
        let t3 = t2 * t0 - d0 - d4;
        let t4 = t4 + t3;
        let z4 = t4;

        let t0 = x0 + x2 + x4;
        let t3 = s0 * t0 - s1;
        let z5 = t3;

        Fq12 {
            c0: Fq6::new(z0, z1, z2),
            c1: Fq6::new(z3, z4, z5)
        }
    }

    pub fn cyclotomic_squared(&self) -> Self {
        let z0 = self.c0.c0;
        let z4 = self.c0.c1;
        let z3 = self.c0.c2;
        let z2 = self.c1.c0;
        let z1 = self.c1.c1;
        let z5 = self.c1.c2;

        let tmp = z0 * z1;
        let t0 = (z0 + z1) * (z1.mul_by_nonresidue() + z0) - tmp - tmp.mul_by_nonresidue();
        let t1 = tmp + tmp;

        let tmp = z2 * z3;
        let t2 = (z2 + z3) * (z3.mul_by_nonresidue() + z2) - tmp - tmp.mul_by_nonresidue();
        let t3 = tmp + tmp;

        let tmp = z4 * z5;
        let t4 = (z4 + z5) * (z5.mul_by_nonresidue() + z4) - tmp - tmp.mul_by_nonresidue();
        let t5 = tmp + tmp;

        let z0 = t0 - z0;
        let z0 = z0 + z0;
        let z0 = z0 + t0;

        let z1 = t1 + z1;
        let z1 = z1 + z1;
        let z1 = z1 + t1;

        let tmp = t5.mul_by_nonresidue();
        let z2 = tmp + z2;
        let z2 = z2 + z2;
        let z2 = z2 + tmp;

        let z3 = t4 - z3;
        let z3 = z3 + z3;
        let z3 = z3 + t4;

        let z4 = t2 - z4;
        let z4 = z4 + z4;
        let z4 = z4 + t2;

        let z5 = t3 + z5;
        let z5 = z5 + z5;
        let z5 = z5 + t3;

        Fq12 {
            c0: Fq6::new(z0, z4, z3),
            c1: Fq6::new(z2, z1, z5)
        }
    }

    pub fn cyclotomic_pow<I: Into<U256>>(&self, by: I) -> Self {
        let mut res = Self::one();

        let mut found_one = false;

        for i in by.into().bits() {
            if found_one {
                res = res.cyclotomic_squared();
            }

            if i {
                found_one = true;
                res = *self * res;
            }
        }

        res
    }
}

impl FieldElement for Fq12 {
    fn zero() -> Self {
        Fq12 {
            c0: Fq6::zero(),
            c1: Fq6::zero()
        }
    }

    fn one() -> Self {
        Fq12 {
            c0: Fq6::one(),
            c1: Fq6::zero()
        }
    }
    
    fn random<R: Rng>(rng: &mut R) -> Self {
        Fq12 {
            c0: Fq6::random(rng),
            c1: Fq6::random(rng)
        }
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn squared(&self) -> Self {
        let ab = self.c0 * self.c1;

        Fq12 {
            c0: (self.c1.mul_by_nonresidue() + self.c0) * (self.c0 + self.c1) - ab - ab.mul_by_nonresidue(),
            c1: ab + ab
        }
    }

    fn inverse(self) -> Option<Self> {
        match (self.c0.squared() - (self.c1.squared().mul_by_nonresidue())).inverse() {
            Some(t) => Some(Fq12 {
                c0: self.c0 * t,
                c1: -(self.c1 * t)
            }),
            None => None
        }
    }
}

impl Mul for Fq12 {
    type Output = Fq12;

    fn mul(self, other: Fq12) -> Fq12 {
        let aa = self.c0 * other.c0;
        let bb = self.c1 * other.c1;

        Fq12 {
            c0: bb.mul_by_nonresidue() + aa,
            c1: (self.c0 + self.c1) * (other.c0 + other.c1) - aa - bb
        }
    }
}

impl Sub for Fq12 {
    type Output = Fq12;

    fn sub(self, other: Fq12) -> Fq12 {
        Fq12 {
            c0: self.c0 - other.c0,
            c1: self.c1 - other.c1
        }
    }
}

impl Add for Fq12 {
    type Output = Fq12;

    fn add(self, other: Fq12) -> Fq12 {
        Fq12 {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1
        }
    }
}

impl Neg for Fq12 {
    type Output = Fq12;

    fn neg(self) -> Fq12 {
        Fq12 {
            c0: -self.c0,
            c1: -self.c1
        }
    }
}
