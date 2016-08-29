mod fp;
mod fq2;

use arith::U256;
use rand::Rng;
use std::ops::{Add, Sub, Mul, Neg};
use std::fmt::Debug;

pub use self::fp::{Fq,Fr,const_fp};
pub use self::fq2::Fq2;

pub trait FieldElement: Sized
                        + Copy
                        + Clone
                        + Add<Output=Self>
                        + Sub<Output=Self>
                        + Mul<Output=Self>
                        + Neg<Output=Self>
                        + PartialEq
                        + Eq
                        + Debug
{
    fn zero() -> Self;
    fn one() -> Self;
    fn random<R: Rng>(&mut R) -> Self;
    fn is_zero(&self) -> bool;
    fn squared(&self) -> Self {
        (*self) * (*self)
    }
    fn inverse(self) -> Option<Self>;
    fn pow<I: Into<U256>>(&self, by: I) -> Self {
        let mut res = Self::one();

        for i in by.into().bits() {
            res = res.squared();
            if i {
                res = *self * res;
            }
        }

        res
    }
}

#[cfg(test)]
mod tests;

#[test]
fn test_fr() {
    tests::field_trials::<Fr>();
}

#[test]
fn test_fq() {
    tests::field_trials::<Fq>();
}

#[test]
fn test_fq2() {
    tests::field_trials::<Fq2>();
}

#[test]
fn test_str() {
    assert_eq!(-Fr::one(), Fr::from_str("21888242871839275222246405745257275088548364400416034343698204186575808495616").unwrap());
    assert_eq!(-Fq::one(), Fq::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208582").unwrap());
}
