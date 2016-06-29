use rand::Rng;
use num::{BigUint, Num};
use std::ops::{Mul,Add,Sub,Neg};
use std::convert::From;
use std::marker::PhantomData;

pub trait PrimeFieldParams {
    fn modulus() -> BigUint;
    fn bits() -> usize;
}

pub struct Fp<P: PrimeFieldParams> {
    value: BigUint,
    _marker: PhantomData<P>
}

impl<P: PrimeFieldParams> Fp<P> {
    pub fn zero() -> Self { unimplemented!() }
    pub fn one() -> Self { unimplemented!() }
    pub fn random<R: Rng>(rng: &mut R) -> Self { unimplemented!() }

    pub fn is_zero(&self) -> bool { unimplemented!() }
    pub fn inverse(&self) -> Self { unimplemented!() }
    pub fn squared(&self) -> Self { unimplemented!() }
    pub fn pow<P2: PrimeFieldParams>(&self, exp: &Fp<P2>) -> Self { unimplemented!() }
    pub fn test_bit(&self, bit: usize) -> bool { unimplemented!() }
}

impl<'a, P: PrimeFieldParams> From<&'a str> for Fp<P> {
    fn from(s: &'a str) -> Self { unimplemented!() }
}

impl<P: PrimeFieldParams> Clone for Fp<P> {
    fn clone(&self) -> Self { unimplemented!() }
}

impl<'a, 'b, P: PrimeFieldParams> Add<&'b Fp<P>> for &'a Fp<P> {
    type Output = Fp<P>;

    fn add(self, other: &Fp<P>) -> Fp<P> { unimplemented!() }
}

impl<'a, 'b, P: PrimeFieldParams> Sub<&'b Fp<P>> for &'a Fp<P> {
    type Output = Fp<P>;

    fn sub(self, other: &Fp<P>) -> Fp<P> { unimplemented!() }
}

impl<'a, 'b, P: PrimeFieldParams> Mul<&'b Fp<P>> for &'a Fp<P> {
    type Output = Fp<P>;

    fn mul(self, other: &Fp<P>) -> Fp<P> { unimplemented!() }
}

impl<'a, P: PrimeFieldParams> Neg for &'a Fp<P> {
    type Output = Fp<P>;

    fn neg(self) -> Fp<P> { unimplemented!() }
}

impl<P: PrimeFieldParams> Neg for Fp<P> {
    type Output = Fp<P>;

    fn neg(self) -> Fp<P> {
        -(&self)
    }
}

forward_all_binop_to_ref_ref!(impl(P: PrimeFieldParams) Mul for Fp<P>, mul);
