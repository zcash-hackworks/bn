#[macro_use]
mod macros;

pub mod fp;
pub mod fp2;
pub mod fp6;
pub mod fp12;

#[cfg(test)]
pub mod tests;

use rand::Rng;
use self::fp::{Fp, PrimeFieldParams};
use std::fmt::Debug;

pub trait Field: Sized + Clone + Debug {
    fn zero() -> Self;
    fn one() -> Self;
    fn random<R: Rng>(rng: &mut R) -> Self;
    fn is_zero(&self) -> bool {
        self.eq(&Self::zero())
    }
    fn inverse(&self) -> Self;
    fn squared(&self) -> Self {
        self.mul(self)
    }
    fn pow<P: PrimeFieldParams>(&self, exp: &Fp<P>) -> Self {
        let mut res = Self::one();

        let mut found_one = false;

        for i in (0..P::bits()).rev() {
            if found_one {
                res = res.squared();
            }

            if exp.test_bit(i) {
                found_one = true;
                res = self.mul(&res);
            }
        }

        res
    }
    fn eq(&self, other: &Self) -> bool;
    fn ne(&self, other: &Self) -> bool {
        !self.eq(other)
    }

    fn neg(&self) -> Self;
    fn mul(&self, other: &Self) -> Self;
    fn sub(&self, other: &Self) -> Self;
    fn add(&self, other: &Self) -> Self;
}
