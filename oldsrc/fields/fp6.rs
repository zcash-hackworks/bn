use ::Fq2;
use rand::Rng;
use fields::Field;
use std::ops::{Mul,Add,Sub,Neg};
use std::cmp::{PartialEq, Eq};
use std::marker::PhantomData;
use std::fmt;

pub trait Fp6Params {
    fn non_residue() -> Fq2;
    fn name() -> &'static str;
    fn frobenius_coeffs_c1(usize) -> Fq2;
    fn frobenius_coeffs_c2(usize) -> Fq2;
}

pub struct Fp6<P: Fp6Params> {
    pub a: Fq2,
    pub b: Fq2,
    pub c: Fq2,
    _marker: PhantomData<P>
}

impl<P: Fp6Params> Fp6<P> {
    pub fn new(a: Fq2, b: Fq2, c: Fq2) -> Self {
        Fp6 {
            a: a,
            b: b,
            c: c,
            _marker: PhantomData
        }
    }

    pub fn frobenius_map(&self, power: usize) -> Self {
        Fp6 {
            a: self.a.frobenius_map(power),
            b: self.b.frobenius_map(power) * P::frobenius_coeffs_c1(power % 6),
            c: self.c.frobenius_map(power) * P::frobenius_coeffs_c2(power % 6),
            _marker: PhantomData
        }
    }
}

impl<P: Fp6Params> fmt::Debug for Fp6<P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}({:?}, {:?}, {:?})", P::name(), self.a, self.b, self.c)
    }
}

impl<P: Fp6Params> Clone for Fp6<P> {
    fn clone(&self) -> Self {
        Fp6 {
            a: self.a.clone(),
            b: self.b.clone(),
            c: self.c.clone(),
            _marker: PhantomData
        }
    }
}

impl<P: Fp6Params> Field for Fp6<P> {
    fn zero() -> Self {
        Fp6 {
            a: Fq2::zero(),
            b: Fq2::zero(),
            c: Fq2::zero(),
            _marker: PhantomData
        }
    }

    fn one() -> Self {
        Fp6 {
            a: Fq2::one(),
            b: Fq2::zero(),
            c: Fq2::zero(),
            _marker: PhantomData
        }
    }

    fn random<R: Rng>(rng: &mut R) -> Self {
        Fp6 {
            a: Fq2::random(rng),
            b: Fq2::random(rng),
            c: Fq2::random(rng),
            _marker: PhantomData
        }
    }

    fn inverse(&self) -> Self {
        let c0 = self.a.squared() - &self.b * &self.c * P::non_residue();
        let c1 = self.c.squared() * P::non_residue() - &self.a * &self.b;
        let c2 = self.b.squared() - &self.a * &self.c;
        let t = ((&self.c * &c1 + &self.b * &c2) * P::non_residue() + &self.a * &c0).inverse();

        Fp6 {
            a: &t * &c0,
            b: &t * &c1,
            c: &t * &c2,
            _marker: PhantomData
        }
    }

    fn eq(&self, other: &Self) -> bool {
        self.a == other.a && self.b == other.b && self.c == other.c
    }

    fn neg(&self) -> Self {
        Fp6 {
        	a: -(&self.a),
        	b: -(&self.b),
        	c: -(&self.c),
        	_marker: PhantomData
        }
    }

    fn mul(&self, other: &Self) -> Self {
        let a_a = &self.a * &other.a;
        let b_b = &self.b * &other.b;
        let c_c = &self.c * &other.c;

        Fp6 {
        	a: ((&self.b + &self.c) * (&other.b + &other.c) - &b_b - &c_c) * P::non_residue() + &a_a,
        	b: (&self.a + &self.b) * (&other.a + &other.b) - &a_a - &b_b + &c_c * P::non_residue(),
        	c: (&self.a + &self.c) * (&other.a + &other.c) - &a_a + &b_b - &c_c,
        	_marker: PhantomData
        }
    }

    fn squared(&self) -> Self {
        let a = &self.a;
        let b = &self.b;
        let c = &self.c;

        let s0 = a.squared();
        let ab = a * b;
        let s1 = &ab + &ab;
        let s2 = (a - b + c).squared();
        let bc = b * c;
        let s3 = &bc + &bc;
        let s4 = c.squared();

        Fp6 {
            a: &s0 + &s3 * P::non_residue(),
            b: &s1 + &s4 * P::non_residue(),
            c: &s1 + &s2 + &s3 - &s0 - &s4,
            _marker: PhantomData
        }
    }

    fn sub(&self, other: &Self) -> Self {
        Fp6 {
            a: &self.a - &other.a,
            b: &self.b - &other.b,
            c: &self.c - &other.c,
            _marker: PhantomData
        }
    }

    fn add(&self, other: &Self) -> Self {
        Fp6 {
            a: &self.a + &other.a,
            b: &self.b + &other.b,
            c: &self.c + &other.c,
            _marker: PhantomData
        }
    }
}

impl<P: Fp6Params> Fp6<P> {
	pub fn mul_by_nonresidue(&self, other: &Fq2) -> Fp6<P> {
		Fp6 {
			a: &self.c * other,
			b: self.a.clone(),
			c: self.b.clone(),
			_marker: PhantomData
		}
	}
}

impl<'a, 'b, P: Fp6Params> Mul<&'a Fq2> for &'b Fp6<P> {
    type Output = Fp6<P>;

    fn mul(self, other: &Fq2) -> Fp6<P> {
        Fp6 {
            a: &self.a * other,
            b: &self.b * other,
            c: &self.c * other,
            _marker: PhantomData
        }
    }
}

forward_ops_to_field_ops!(impl(P: Fp6Params) Fp6<P>);
