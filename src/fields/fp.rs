use rand::Rng;
use std::ops::{Add, Sub, Mul, Neg};
use std::fmt;
use std::marker::PhantomData;
use super::FieldElement;

use arith::U256;

pub trait FpParams {
    fn name() -> &'static str;
    fn modulus() -> U256;
    fn inv() -> u64;
    fn rsquared() -> U256;
    fn rcubed() -> U256;
    fn one() -> U256;
}

pub struct Fp<P: FpParams>(U256, PhantomData<P>);
impl<P: FpParams> Copy for Fp<P> { }
impl<P: FpParams> Clone for Fp<P> {
    fn clone(&self) -> Self { *self }
}
impl<P: FpParams> PartialEq for Fp<P> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}
impl<P: FpParams> Eq for Fp<P> { }
impl<P: FpParams> fmt::Debug for Fp<P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}({:?})", P::name(), self.0)
    }
}
impl<P: FpParams> From<Fp<P>> for U256 {
    fn from(mut a: Fp<P>) -> Self {
        a.0.mul(&U256::one(), &P::modulus(), P::inv());

        a.0
    }
}

#[inline]
pub fn const_fp<P: FpParams, I: Into<U256>>(i: I) -> Fp<P> {
    Fp(i.into(), PhantomData)
}

impl<P: FpParams> Fp<P> {
    pub fn from_str(s: &str) -> Option<Self> {
        let ints: Vec<_> = {
            let mut acc = Self::zero();
            (0..11).map(|_| {let tmp = acc; acc = acc + Self::one(); tmp}).collect()
        };

        let mut res = Self::zero();
        for c in s.chars() {
            match c.to_digit(10) {
                Some(d) => {
                    res = res * ints[10];
                    res = res + ints[d as usize];
                },
                None => {
                    return None;
                }
            }
        }

        Some(res)
    }
}

impl<P: FpParams> Fp<P> {
    /// Converts a U256 to an Fp so long as it's below the modulus.
    pub fn new(mut a: U256) -> Option<Self> {
        if a < P::modulus() {
            a.mul(&P::rsquared(), &P::modulus(), P::inv());

            Some(Fp(a, PhantomData))
        } else {
            None
        }
    }
}

impl<P: FpParams> FieldElement for Fp<P> {
    fn zero() -> Self {
        const_fp(U256::zero())
    }

    fn one() -> Self {
        const_fp(P::one())
    }
    
    fn random<R: Rng>(rng: &mut R) -> Self {
        const_fp(U256::random(rng, &P::modulus()))
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    fn inverse(mut self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            self.0.invert(&P::modulus());
            self.0.mul(&P::rcubed(), &P::modulus(), P::inv());

            Some(self)
        }
    }
}

impl<P: FpParams> Add for Fp<P> {
    type Output = Fp<P>;

    fn add(mut self, other: Fp<P>) -> Fp<P> {
        self.0.add(&other.0, &P::modulus());

        self
    }
}

impl<P: FpParams> Sub for Fp<P> {
    type Output = Fp<P>;

    fn sub(mut self, other: Fp<P>) -> Fp<P> {
        self.0.sub(&other.0, &P::modulus());

        self
    }
}

impl<P: FpParams> Mul for Fp<P> {
    type Output = Fp<P>;

    fn mul(mut self, other: Fp<P>) -> Fp<P> {
        self.0.mul(&other.0, &P::modulus(), P::inv());

        self
    }
}

impl<P: FpParams> Neg for Fp<P> {
    type Output = Fp<P>;

    fn neg(mut self) -> Fp<P> {
        self.0.neg(&P::modulus());

        self
    }
}

pub struct FrParams;
pub type Fr = Fp<FrParams>;
impl FpParams for FrParams {
    fn name() -> &'static str { "Fr" }

    #[inline]
    fn modulus() -> U256 {
        // 21888242871839275222246405745257275088548364400416034343698204186575808495617
        [0x43e1f593f0000001, 0x2833e84879b97091, 0xb85045b68181585d, 0x30644e72e131a029].into()
    }

    #[inline]
    fn inv() -> u64 {
        0xc2e1f593efffffff
    }

    #[inline]
    fn rsquared() -> U256 {
        // 944936681149208446651664254269745548490766851729442924617792859073125903783
        [0x1bb8e645ae216da7, 0x53fe3ab1e35c59e3, 0x8c49833d53bb8085, 0x0216d0b17f4e44a5].into()
    }

    #[inline]
    fn rcubed() -> U256 {
        // 5866548545943845227489894872040244720403868105578784105281690076696998248512
        [0x5e94d8e1b4bf0040, 0x2a489cbe1cfbb6b8, 0x893cc664a19fcfed, 0x0cf8594b7fcc657c].into()
    }

    #[inline]
    fn one() -> U256 {
        [0xac96341c4ffffffb, 0x36fc76959f60cd29, 0x666ea36f7879462e, 0xe0a77c19a07df2f].into()
    }
}

pub struct FqParams;
pub type Fq = Fp<FqParams>;
impl FpParams for FqParams {
    fn name() -> &'static str { "Fq" }

    #[inline]
    fn modulus() -> U256 {
        // 21888242871839275222246405745257275088696311157297823662689037894645226208583
        [0x3c208c16d87cfd47, 0x97816a916871ca8d, 0xb85045b68181585d, 0x30644e72e131a029].into()
    }

    #[inline]
    fn inv() -> u64 {
        0x87d20782e4866389
    }

    #[inline]
    fn rsquared() -> U256 {
        // 3096616502983703923843567936837374451735540968419076528771170197431451843209
        [0xf32cfc5b538afa89, 0xb5e71911d44501fb, 0x47ab1eff0a417ff6, 0x06d89f71cab8351f].into()
    }

    #[inline]
    fn rcubed() -> U256 {
        // 14921786541159648185948152738563080959093619838510245177710943249661917737183
        [0xb1cd6dafda1530df, 0x62f210e6a7283db6, 0xef7f0b0c0ada0afb, 0x20fd6e902d592544].into()
    }

    #[inline]
    fn one() -> U256 {
        [0xd35d438dc58f0d9d, 0xa78eb28f5c70b3d, 0x666ea36f7879462c, 0xe0a77c19a07df2f].into()
    }
}

#[test]
fn test_rsquared() {
    let rng = &mut ::rand::thread_rng();

    for _ in 0..1000 {
        let a = Fr::random(rng);
        let b: U256 = a.into();
        let c = Fr::new(b).unwrap();

        assert_eq!(a, c);
    }

    for _ in 0..1000 {
        let a = Fq::random(rng);
        let b: U256 = a.into();
        let c = Fq::new(b).unwrap();

        assert_eq!(a, c);
    }
}
