use rand::Rng;
use std::ops::{Add, Sub, Mul, Neg};
use super::FieldElement;

use rustc_serialize::{Encodable, Encoder, Decodable, Decoder};

use arith::{U512, U256};

macro_rules! field_impl {
    ($name:ident, $modulus:expr, $rsquared:expr, $rcubed:expr, $one:expr, $inv:expr) => {
        #[derive(Copy, Clone, PartialEq, Eq, Debug)]
        #[repr(C)]
        pub struct $name(U256);

        impl From<$name> for U256 {
            #[inline]
            fn from(mut a: $name) -> Self {
                a.0.mul(&U256::one(), &U256($modulus), $inv);
                
                a.0
            }
        }

        impl Encodable for $name {
            fn encode<S: Encoder>(&self, s: &mut S) -> Result<(), S::Error> {
                let normalized = U256::from(*self);

                normalized.encode(s)
            }
        }

        impl Decodable for $name {
            fn decode<S: Decoder>(s: &mut S) -> Result<$name, S::Error> {
                $name::new(try!(U256::decode(s))).ok_or_else(|| s.error("integer is not less than modulus"))
            }
        }

        impl $name {
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

            /// Converts a U256 to an Fp so long as it's below the modulus.
            pub fn new(mut a: U256) -> Option<Self> {
                if a < U256($modulus) {
                    a.mul(&U256($rsquared), &U256($modulus), $inv);

                    Some($name(a))
                } else {
                    None
                }
            }

            pub fn interpret(buf: &[u8; 64]) -> Self {
                $name::new(U512::interpret(buf).divrem(&U256($modulus)).1).unwrap()
            }

            /// Returns the modulus
            #[inline]
            pub fn modulus() -> U256 {
                U256($modulus)
            }
        }

        impl FieldElement for $name {
            #[inline]
            fn zero() -> Self {
                $name(U256([0, 0, 0, 0]))
            }

            #[inline]
            fn one() -> Self {
                $name(U256($one))
            }
            
            fn random<R: Rng>(rng: &mut R) -> Self {
                $name(U256::random(rng, &U256($modulus)))
            }

            #[inline]
            fn is_zero(&self) -> bool {
                self.0.is_zero()
            }

            fn inverse(mut self) -> Option<Self> {
                if self.is_zero() {
                    None
                } else {
                    self.0.invert(&U256($modulus));
                    self.0.mul(&U256($rcubed), &U256($modulus), $inv);

                    Some(self)
                }
            }
        }

        impl Add for $name {
            type Output = $name;

            #[inline]
            fn add(mut self, other: $name) -> $name {
                self.0.add(&other.0, &U256($modulus));

                self
            }
        }

        impl Sub for $name {
            type Output = $name;

            #[inline]
            fn sub(mut self, other: $name) -> $name {
                self.0.sub(&other.0, &U256($modulus));

                self
            }
        }

        impl Mul for $name {
            type Output = $name;

            #[inline]
            fn mul(mut self, other: $name) -> $name {
                self.0.mul(&other.0, &U256($modulus), $inv);

                self
            }
        }

        impl Neg for $name {
            type Output = $name;

            #[inline]
            fn neg(mut self) -> $name {
                self.0.neg(&U256($modulus));

                self
            }
        }
    }
}

field_impl!(
    Fr,
    [0x43e1f593f0000001, 0x2833e84879b97091, 0xb85045b68181585d, 0x30644e72e131a029],
    [0x1bb8e645ae216da7, 0x53fe3ab1e35c59e3, 0x8c49833d53bb8085, 0x0216d0b17f4e44a5],
    [0x5e94d8e1b4bf0040, 0x2a489cbe1cfbb6b8, 0x893cc664a19fcfed, 0x0cf8594b7fcc657c],
    [0xac96341c4ffffffb, 0x36fc76959f60cd29, 0x666ea36f7879462e, 0xe0a77c19a07df2f],
    0xc2e1f593efffffff
);

field_impl!(
    Fq,
    [0x3c208c16d87cfd47, 0x97816a916871ca8d, 0xb85045b68181585d, 0x30644e72e131a029],
    [0xf32cfc5b538afa89, 0xb5e71911d44501fb, 0x47ab1eff0a417ff6, 0x06d89f71cab8351f],
    [0xb1cd6dafda1530df, 0x62f210e6a7283db6, 0xef7f0b0c0ada0afb, 0x20fd6e902d592544],
    [0xd35d438dc58f0d9d, 0xa78eb28f5c70b3d, 0x666ea36f7879462c, 0xe0a77c19a07df2f],
    0x87d20782e4866389
);

#[inline]
pub fn const_fq(i: [u64; 4]) -> Fq {
    Fq(U256(i))
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
