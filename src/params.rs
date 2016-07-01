use num::{Num,BigUint};
use fields::fp::PrimeFieldParams;

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

    fields::tests::field_trials::<fields::fp::Fp<FrParams>>();
}

#[test]
fn test_fq() {
    use fields;

    fields::tests::field_trials::<fields::fp::Fp<FqParams>>();
}
