use ::Fq;
use rand::Rng;
use fields::Field;
use std::ops::{Mul,Add,Sub,Neg};
use std::cmp::{PartialEq, Eq};
use std::marker::PhantomData;
use std::fmt;

pub trait Fp2Params {
    fn non_residue() -> Fq;
    fn name() -> &'static str;
    fn frobenius_coeffs_c1(usize) -> Fq;
}

pub struct Fp2<P: Fp2Params> {
    a: Fq,
    b: Fq,
    _marker: PhantomData<P>
}

impl<P: Fp2Params> Fp2<P> {
    pub fn new(a: Fq, b: Fq) -> Self {
        Fp2 {
            a: a,
            b: b,
            _marker: PhantomData
        }
    }

    pub fn frobenius_map(&self, power: usize) -> Self {
        Fp2 {
            a: self.a.clone(),
            b: &self.b * P::frobenius_coeffs_c1(power % 2),
            _marker: PhantomData
        }
    }
}

impl<P: Fp2Params> fmt::Debug for Fp2<P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}({:?}, {:?})", P::name(), self.a, self.b)
    }
}

impl<P: Fp2Params> Clone for Fp2<P> {
    fn clone(&self) -> Self {
        Fp2 {
            a: self.a.clone(),
            b: self.b.clone(),
            _marker: PhantomData
        }
    }
}

impl<P: Fp2Params> Field for Fp2<P> {
    fn zero() -> Self {
        Fp2 {
            a: Fq::zero(),
            b: Fq::zero(),
            _marker: PhantomData
        }
    }
    fn one() -> Self {
        Fp2 {
            a: Fq::one(),
            b: Fq::zero(),
            _marker: PhantomData
        }
    }
    fn random<R: Rng>(rng: &mut R) -> Self {
        Fp2 {
            a: Fq::random(rng),
            b: Fq::random(rng),
            _marker: PhantomData
        }
    }
    
    fn inverse(&self) -> Self {
        let t = (self.a.squared() - (self.b.squared() * P::non_residue())).inverse();

        Fp2 {
            a: &self.a * &t,
            b: -(&self.b * &t),
            _marker: PhantomData
        }
    }
    fn squared(&self) -> Self {
        let a = &self.a; let b = &self.b;
        let ab = &(a * b);

        Fp2 {
            a: (b * P::non_residue() + a) * (a + b) - ab - ab * P::non_residue(),
            b: ab + ab,
            _marker: PhantomData
        }
    }

    fn eq(&self, other: &Self) -> bool {
        self.a == other.a && self.b == other.b
    }
    fn neg(&self) -> Self {
        Fp2 {
            a: -(&self.a),
            b: -(&self.b),
            _marker: PhantomData
        }
    }
    fn mul(&self, other: &Self) -> Self {
        let a_a = &(&self.a * &other.a);
        let b_b = &(&self.b * &other.b);

        Fp2 {
            a: b_b * P::non_residue() + a_a,
            b: (&self.a + &self.b) * (&other.a + &other.b) - a_a - b_b,
            _marker: PhantomData
        }
    }
    fn sub(&self, other: &Self) -> Self {
        Fp2 {
            a: &self.a - &other.a,
            b: &self.b - &other.b,
            _marker: PhantomData
        }
    }
    fn add(&self, other: &Self) -> Self {
        Fp2 {
            a: &self.a + &other.a,
            b: &self.b + &other.b,
            _marker: PhantomData
        }
    }
}

impl<'a, 'b, P: Fp2Params> Mul<&'a Fq> for &'b Fp2<P> {
    type Output = Fp2<P>;

    fn mul(self, other: &Fq) -> Fp2<P> {
        Fp2 {
            a: &self.a * other,
            b: &self.b * other,
            _marker: PhantomData
        }
    }
}

forward_ops_to_field_ops!(impl(P: Fp2Params) Fp2<P>);
