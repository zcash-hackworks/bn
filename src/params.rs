use num::{Num,BigUint};
use fields::Field;
use fields::fp::PrimeFieldParams;
use super::{Fr,Fq,G1};
use groups::*;

pub struct FrParams;

impl PrimeFieldParams for FrParams {
    fn modulus() -> BigUint {
        BigUint::from_str_radix("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap()
    }
    fn bits() -> usize { 254 }
    fn name() -> &'static str { "Fr" }
}

pub struct FqParams;

impl PrimeFieldParams for FqParams {
    fn modulus() -> BigUint {
        BigUint::from_str_radix("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap()
    }
    fn bits() -> usize { 254 }
    fn name() -> &'static str { "Fq" }
}

#[test]
fn test_fr() {
    use fields;

    fields::tests::field_trials::<super::Fr>();
}

#[test]
fn test_fq() {
    use fields;

    fields::tests::field_trials::<super::Fq>();
}

pub struct G1Params;

impl GroupParams for G1Params {
    type Base = super::Fq;

    fn name() -> &'static str {
        "G1"
    }
    fn zero() -> Jacobian<Self> {
        Jacobian::new(Fq::zero(), Fq::one(), Fq::zero()).unwrap()
    }
    fn one() -> Jacobian<Self> {
        Jacobian::new(Fq::from("1"), Fq::from("2"), Fq::one()).unwrap()
    }
    fn coeff_b() -> Self::Base {
        Fq::from("3")
    }
}

#[test]
fn test_g1() {
    use groups;

    groups::tests::group_trials::<G1Params>();
}

#[test]
fn g1_test_vector() {
    let a = G1::one() * &Fr::from("19797905000333868150253315089095386158892526856493194078073564469188852136946");
    let b = G1::one() * &Fr::from("2730506433347642574983433139433778984782882168213690554721050571242082865799");
    let e = &a + &b;

    let expect = G1::new(
        Fq::from("18450621724990678172567114131642278789161361170999664461184794604011563728206"),
        Fq::from("21329688341674583036384007811166435666174342925504675855816423131698588368496"),
        Fq::one()
    ).unwrap();

    assert!(expect.eq(&e));
}
