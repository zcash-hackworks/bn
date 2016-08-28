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
    fn inverse(self) -> Self;
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

    assert_eq!(
        const_fp::<self::fp::FrParams, _>([0xf0000000, 0x43e1f593, 0x79b97091, 0x2833e848, 0x8181585d, 0xb85045b6, 0xe131a029, 0x30644e72]).inverse(),
        const_fp::<self::fp::FrParams, _>([1105105498, 673779534, 2522683054, 3560287638, 767940567, 738640505, 1642290052, 776830401])
    );
}

#[test]
fn test_fq() {
    tests::field_trials::<Fq>();

    assert_eq!(
        const_fp::<self::fp::FqParams, _>([0xd87cfd46, 0x3c208c16, 0x6871ca8d, 0x97816a91, 0x8181585d, 0xb85045b6, 0xe131a029, 0x30644e72]).inverse(),
        const_fp::<self::fp::FqParams, _>([2230452926, 1223921595, 2485962897, 3784987007, 2000672870, 1889871543, 377056010, 697020161])
    );
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

#[test]
fn test_fr_squaring() {
    // Test vector.
    assert_eq!(
        Fr::from_str("21888242871839275222246405745257275088548364400416034343698204186575808495610").unwrap().pow([0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff]),
        const_fp([2572677806, 4136501712, 2816341579, 1079830324, 2386666211, 2796237710, 2115470140, 673230670])
    );
}

#[test]
fn test_fq_squaring() {
    // Test vector
    assert_eq!(
        Fq::from_str("21888242871839275222246405745257275088696311157297823662689037894645226208540").unwrap().pow([0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff]),
        const_fp([1128196772, 3636718463, 2555618142, 559034227, 1481851987, 3987127805, 2583602581, 88292310])
    );
}
