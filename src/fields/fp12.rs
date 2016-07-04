use ::Fr;
use ::Fq2;
use ::Fq6;
use ::Fq12;
use fields::fp::PrimeFieldParams;
use fields::fp6::Fp6Params;
use params::{FrParams,Fq6Params};
use rand::Rng;
use fields::Field;
use std::ops::{Mul,Add,Sub,Neg};
use std::cmp::{PartialEq, Eq};
use std::marker::PhantomData;
use std::fmt;

pub trait Fp12Params {
    fn non_residue() -> Fq2;
    fn name() -> &'static str;
    fn frobenius_coeffs_c1(n: usize) -> Fq2;
}

pub struct Fp12<P: Fp12Params> {
    a: Fq6,
    b: Fq6,
    _marker: PhantomData<P>
}

impl<P: Fp12Params> Fp12<P> {
    pub fn new(a: Fq6, b: Fq6) -> Self {
        Fp12 {
            a: a,
            b: b,
            _marker: PhantomData
        }
    }

    pub fn frobenius_map(&self, power: usize) -> Self {
        Fp12 {
            a: self.a.frobenius_map(power),
            b: &self.b.frobenius_map(power) * &P::frobenius_coeffs_c1(power % 12),
            _marker: PhantomData
        }
    }
}

impl Fq12 {
    pub fn unitary_inverse(&self) -> Fq12 {
         Fp12 {
            a: self.a.clone(),
            b: -&self.b,
            _marker: PhantomData
        }
    }

    pub fn cyclotomic_exp(&self, exp: &Fr) -> Self {
        let mut res = Self::one();

        let mut found_one = false;

        for i in (0..FrParams::bits()).rev() {
            if found_one {
                res = res.cyclotomic_squared();
            }

            if exp.test_bit(i) {
                found_one = true;
                res = self * &res;
            }
        }

        res
    }

    pub fn cyclotomic_squared(&self) -> Self {
        let z0 = &self.a.a;
        let z4 = &self.a.b;
        let z3 = &self.a.c;
        let z2 = &self.b.a;
        let z1 = &self.b.b;
        let z5 = &self.b.c;

        let tmp = z0 * z1;
        let t0 = (z0 + z1) * (z1 * Fq6Params::non_residue() + z0) - &tmp - &tmp * Fq6Params::non_residue();
        let t1 = &tmp + &tmp;

        let tmp = z2 * z3;
        let t2 = (z2 + z3) * (z3 * Fq6Params::non_residue() + z2) - &tmp - &tmp * Fq6Params::non_residue();
        let t3 = &tmp + &tmp;

        let tmp = z4 * z5;
        let t4 = (z4 + z5) * (z5 * Fq6Params::non_residue() + z4) - &tmp - &tmp * Fq6Params::non_residue();
        let t5 = &tmp + &tmp;

        let z0 = &t0 - z0;
        let z0 = &z0 + &z0;
        let z0 = &z0 + &t0;

        let z1 = &t1 + z1;
        let z1 = &z1 + &z1;
        let z1 = &z1 + &t1;

        let tmp = &t5 * Fq6Params::non_residue();
        let z2 = &tmp + z2;
        let z2 = &z2 + &z2;
        let z2 = &z2 + &tmp;

        let z3 = &t4 - z3;
        let z3 = &z3 + &z3;
        let z3 = &z3 + &t4;

        let z4 = &t2 - z4;
        let z4 = &z4 + &z4;
        let z4 = &z4 + &t2;

        let z5 = &t3 + z5;
        let z5 = &z5 + &z5;
        let z5 = &z5 + &t3;

        Fp12 {
            a: Fq6::new(z0, z4, z3),
            b: Fq6::new(z2, z1, z5),
            _marker: PhantomData
        }
    }

    pub fn mul_by_024(&self,
                      ell_0: &Fq2,
                      ell_vw: Fq2,
                      ell_vv: Fq2) -> Fq12 {
        let z0 = &self.a.a;
        let z1 = &self.a.b;
        let z2 = &self.a.c;
        let z3 = &self.b.a;
        let z4 = &self.b.b;
        let z5 = &self.b.c;

        let x0 = ell_0;
        let x2 = &ell_vv;
        let x4 = &ell_vw;

        let d0 = z0 * x0;
        let d2 = z2 * x2;
        let d4 = z4 * x4;
        let t2 = z0 + z4;
        let t1 = z0 + z2;
        let s0 = z1 + z3 + z5;

        let s1 = z1 * x2;
        let t3 = &s1 + &d4;
        let t4 = &t3 * Fq6Params::non_residue() + &d0;
        let z0 = t4;

        let t3 = z5 * x4;
        let s1 = &s1 + &t3;
        let t3 = &t3 + &d2;
        let t4 = &t3 * Fq6Params::non_residue();
        let t3 = z1 * x0;
        let s1 = &s1 + &t3;
        let t4 = &t4 + &t3;
        let z1 = t4;

        let t0 = x0 + x2;
        let t3 = &t1 * &t0 - &d0 - &d2;
        let t4 = z3 * x4;
        let s1 = &s1 + &t4;
        let t3 = &t3 + &t4;

        let t0 = z2 + z4;
        let z2 = t3;

        let t1 = x2 + x4;
        let t3 = &t0 * &t1 - &d2 - &d4;
        let t4 = &t3 * Fq6Params::non_residue();
        let t3 = z3 * x0;
        let s1 = &s1 + &t3;
        let t4 = &t4 + &t3;
        let z3 = t4;

        let t3 = z5 * x2;
        let s1 = &s1 + &t3;
        let t4 = &t3 * Fq6Params::non_residue();
        let t0 = x0 + x4;
        let t3 = &t2 * &t0 - &d0 - &d4;
        let t4 = &t4 + &t3;
        let z4 = t4;

        let t0 = x0 + x2 + x4;
        let t3 = &s0 * &t0 - &s1;
        let z5 = t3;

        Fq12 {
            a: Fq6::new(z0, z1, z2),
            b: Fq6::new(z3, z4, z5),
            _marker: PhantomData
        }
    }
}

impl<P: Fp12Params> fmt::Debug for Fp12<P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}({:?}, {:?})", P::name(), self.a, self.b)
    }
}

impl<P: Fp12Params> Clone for Fp12<P> {
    fn clone(&self) -> Self {
        Fp12 {
            a: self.a.clone(),
            b: self.b.clone(),
            _marker: PhantomData
        }
    }
}

impl<P: Fp12Params> Field for Fp12<P> {
    fn zero() -> Self {
        Fp12 {
            a: Fq6::zero(),
            b: Fq6::zero(),
            _marker: PhantomData
        }
    }
    fn one() -> Self {
        Fp12 {
            a: Fq6::one(),
            b: Fq6::zero(),
            _marker: PhantomData
        }
    }
    fn random<R: Rng>(rng: &mut R) -> Self {
        Fp12 {
            a: Fq6::random(rng),
            b: Fq6::random(rng),
            _marker: PhantomData
        }
    }
    
    fn inverse(&self) -> Self {
        let t = (self.a.squared() - (self.b.squared().mul_by_nonresidue(&P::non_residue()))).inverse();

        Fp12 {
            a: &self.a * &t,
            b: -(&self.b * &t),
            _marker: PhantomData
        }
    }
    
    fn squared(&self) -> Self {
        let a = &self.a; let b = &self.b;
        let ab = &(a * b);

        Fp12 {
            a: (b.mul_by_nonresidue(&P::non_residue()) + a) * (a + b) - ab - ab.mul_by_nonresidue(&P::non_residue()),
            b: ab + ab,
            _marker: PhantomData
        }
    }
    
    fn eq(&self, other: &Self) -> bool {
        self.a == other.a && self.b == other.b
    }
    fn neg(&self) -> Self {
        Fp12 {
            a: -(&self.a),
            b: -(&self.b),
            _marker: PhantomData
        }
    }
    fn mul(&self, other: &Self) -> Self {
        let a_a = &(&self.a * &other.a);
        let b_b = &(&self.b * &other.b);

        Fp12 {
            a: b_b.mul_by_nonresidue(&P::non_residue()) + a_a,
            b: (&self.a + &self.b) * (&other.a + &other.b) - a_a - b_b,
            _marker: PhantomData
        }
    }
    fn sub(&self, other: &Self) -> Self {
        Fp12 {
            a: &self.a - &other.a,
            b: &self.b - &other.b,
            _marker: PhantomData
        }
    }
    fn add(&self, other: &Self) -> Self {
        Fp12 {
            a: &self.a + &other.a,
            b: &self.b + &other.b,
            _marker: PhantomData
        }
    }
}

forward_ops_to_field_ops!(impl(P: Fp12Params) Fp12<P>);
