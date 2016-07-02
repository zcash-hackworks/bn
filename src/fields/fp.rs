use rand::Rng;
use num::{BigUint, Num};
use std::ops::{Mul,Add,Sub,Neg};
use std::cmp::{PartialEq, Eq};
use std::convert::From;
use std::fmt;
use std::marker::PhantomData;
use super::Field;

pub trait PrimeFieldParams {
    fn modulus() -> BigUint;
    fn bits() -> usize;
    fn name() -> &'static str;
}

pub struct Fp<P: PrimeFieldParams> {
    value: BigUint,
    _marker: PhantomData<P>
}

impl<P: PrimeFieldParams> fmt::Debug for Fp<P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}({})", P::name(), self.value)
    }
}

impl<P: PrimeFieldParams> Field for Fp<P> {
    fn zero() -> Self {
        use num::Zero;

        Fp {
            value: BigUint::zero(),
            _marker: PhantomData
        }
    }
    fn one() -> Self {
        use num::One;

        Fp {
            value: BigUint::one(),
            _marker: PhantomData
        }
    }
    fn random<R: Rng>(rng: &mut R) -> Self {
        use num::num_bigint::RandBigInt;
        use num::Zero;

        Fp {
            value: rng.gen_biguint_range(&BigUint::zero(), &P::modulus()),
            _marker: PhantomData
        }
    }

    fn is_zero(&self) -> bool {
        use num::Zero;

        self.value == BigUint::zero()
    }

    fn inverse(&self) -> Self {
        if self.is_zero() {
            // TODO: this should likely bleed through the abstraction layers
            panic!("cannot get the multiplicative inverse of zero")
        } else {
            let mut res = Self::one();

            let mut found_one = false;

            let exp = Self::zero() - Self::one() - Self::one();

            for i in (0..P::bits()).rev() {
                if found_one {
                    res = res.squared();
                }

                if exp.test_bit(i) {
                    found_one = true;
                    res = self * &res;
                }
            }

            res
        }
    }
    
    fn pow<P2: PrimeFieldParams>(&self, exp: &Fp<P2>) -> Self {
        let mut res = Self::one();

        let mut found_one = false;

        for i in (0..P2::bits()).rev() {
            if found_one {
                res = res.squared();
            }

            if exp.test_bit(i) {
                found_one = true;
                res = res * self;
            }
        }

        res
    }

    fn neg(&self) -> Self {
        use num::Zero;

        Fp {
            value: if self.value.is_zero() {
                self.value.clone()
            } else {
                P::modulus() - &self.value
            },
            _marker: PhantomData
        }
    }

    fn mul(&self, other: &Self) -> Self {
        Fp {
            value: (&self.value * &other.value) % &P::modulus(),
            _marker: PhantomData
        }
    }

    fn sub(&self, other: &Self) -> Self {
        if other.value > self.value {
            Fp {
                value: (&self.value + P::modulus()) - &other.value,
                _marker: PhantomData
            }
        } else {
            Fp {
                value: &self.value - &other.value,
                _marker: PhantomData
            }
        }
    }

    fn add(&self, other: &Self) -> Self {
        let tmp = &self.value + &other.value;
        if tmp >= P::modulus() {
            Fp {
                value: tmp - P::modulus(),
                _marker: PhantomData
            }
        } else {
            Fp {
                value: tmp,
                _marker: PhantomData
            }
        }
    }

    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl<P: PrimeFieldParams> Fp<P> {
    pub fn test_bit(&self, bit: usize) -> bool {
        // TODO: This is a naive approach.
        use num::{One, Zero};

        let mut b = BigUint::one();
        let two = &b + &b;
        for _ in 0..bit {
            b = &b + &b;
        }

        (&self.value / b) % two != BigUint::zero()
    }
}

impl<'a, P: PrimeFieldParams> From<&'a str> for Fp<P> {
    fn from(s: &'a str) -> Self {
        Fp {
            value: BigUint::from_str_radix(s, 10).unwrap() % P::modulus(),
            _marker: PhantomData
        }
    }
}

impl<P: PrimeFieldParams> Clone for Fp<P> {
    fn clone(&self) -> Self {
        Fp {
            value: self.value.clone(),
            _marker: PhantomData
        }
    }
}

forward_ops_to_field_ops!(impl(P: PrimeFieldParams) Fp<P>);
