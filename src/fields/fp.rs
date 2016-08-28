use rand::Rng;
use std::ops::{Add, Sub, Mul, Neg};
use std::fmt;
use std::marker::PhantomData;
use super::FieldElement;

use arith::U256;

pub trait FpParams {
    fn name() -> &'static str;
    fn modulus() -> U256;
    fn inv() -> u32;
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
    /// Assumes input is mod p, not exposed publicly
    fn new_checked(mut a: U256) -> Self {
        a.mul(&P::rsquared(), &P::modulus(), P::inv());

        Fp(a, PhantomData)
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
        Fp::new_checked(U256::rand(rng, &P::modulus()))
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    fn inverse(mut self) -> Self {
        assert!(!self.is_zero());

        self.0.invert(&P::modulus());
        self.0.mul(&P::rcubed(), &P::modulus(), P::inv());

        self
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
        [0xf0000001, 0x43e1f593, 0x79b97091, 0x2833e848, 0x8181585d, 0xb85045b6, 0xe131a029, 0x30644e72].into()
    }

    #[inline]
    fn inv() -> u32 {
        0xefffffff
    }

    #[inline]
    fn rsquared() -> U256 {
        // 944936681149208446651664254269745548490766851729442924617792859073125903783
        [0xae216da7, 0x1bb8e645, 0xe35c59e3, 0x53fe3ab1, 0x53bb8085, 0x8c49833d, 0x7f4e44a5, 0x0216d0b1].into()
    }

    #[inline]
    fn rcubed() -> U256 {
        // 5866548545943845227489894872040244720403868105578784105281690076696998248512
        [0xb4bf0040, 0x5e94d8e1, 0x1cfbb6b8, 0x2a489cbe, 0xa19fcfed, 0x893cc664, 0x7fcc657c, 0x0cf8594b].into()
    }

    #[inline]
    fn one() -> U256 {
        [0x4ffffffb, 0xac96341c, 0x9f60cd29, 0x36fc7695, 0x7879462e, 0x666ea36f, 0x9a07df2f, 0x0e0a77c1].into()
    }
}

pub struct FqParams;
pub type Fq = Fp<FqParams>;
impl FpParams for FqParams {
    fn name() -> &'static str { "Fq" }

    #[inline]
    fn modulus() -> U256 {
        // 21888242871839275222246405745257275088696311157297823662689037894645226208583
        [0xd87cfd47, 0x3c208c16, 0x6871ca8d, 0x97816a91, 0x8181585d, 0xb85045b6, 0xe131a029, 0x30644e72].into()
    }

    #[inline]
    fn inv() -> u32 {
        0xe4866389
    }

    #[inline]
    fn rsquared() -> U256 {
        // 3096616502983703923843567936837374451735540968419076528771170197431451843209
        [0x538afa89, 0xf32cfc5b, 0xd44501fb, 0xb5e71911, 0x0a417ff6, 0x47ab1eff, 0xcab8351f, 0x06d89f71].into()
    }

    #[inline]
    fn rcubed() -> U256 {
        // 14921786541159648185948152738563080959093619838510245177710943249661917737183
        [0xda1530df, 0xb1cd6daf, 0xa7283db6, 0x62f210e6, 0x0ada0afb, 0xef7f0b0c, 0x2d592544, 0x20fd6e90].into()
    }

    #[inline]
    fn one() -> U256 {
        [0xc58f0d9d, 0xd35d438d, 0xf5c70b3d, 0x0a78eb28, 0x7879462c, 0x666ea36f, 0x9a07df2f, 0x0e0a77c1].into()
    }
}

#[test]
fn test_rsquared() {
    let rng = &mut ::rand::thread_rng();

    for _ in 0..1000 {
        let a = Fr::random(rng);
        let b: U256 = a.into();
        let c = Fr::new_checked(b);

        assert_eq!(a, c);
    }

    for _ in 0..1000 {
        let a = Fq::random(rng);
        let b: U256 = a.into();
        let c = Fq::new_checked(b);

        assert_eq!(a, c);
    }
}
