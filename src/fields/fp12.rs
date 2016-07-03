use ::Fq2;
use ::Fq6;
use rand::Rng;
use fields::Field;
use std::ops::{Mul,Add,Sub,Neg};
use std::cmp::{PartialEq, Eq};
use std::marker::PhantomData;
use std::fmt;

pub trait Fp12Params {
    fn non_residue() -> Fq2;
    fn name() -> &'static str;
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
