use rand::Rng;
use num::{BigUint, Num};
use std::ops::{Mul,Add,Sub,Neg};
use std::cmp::{PartialEq, Eq};
use std::convert::From;
use std::fmt;
use std::marker::PhantomData;

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

impl<P: PrimeFieldParams> Fp<P> {
    pub fn zero() -> Self {
        use num::Zero;

        Fp {
            value: BigUint::zero(),
            _marker: PhantomData
        }
    }
    pub fn one() -> Self {
        use num::One;

        Fp {
            value: BigUint::one(),
            _marker: PhantomData
        }
    }
    pub fn random<R: Rng>(rng: &mut R) -> Self {
        use num::num_bigint::RandBigInt;
        use num::Zero;

        Fp {
            value: rng.gen_biguint_range(&BigUint::zero(), &P::modulus()),
            _marker: PhantomData
        }
    }

    pub fn is_zero(&self) -> bool {
        use num::Zero;

        self.value == BigUint::zero()
    }

    pub fn inverse(&self) -> Self {
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
    pub fn squared(&self) -> Self {
        self * self
    }
    pub fn pow<P2: PrimeFieldParams>(&self, exp: &Fp<P2>) -> Self {
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
    fn clone(&self) -> Self { unimplemented!() }
}

impl<'a, 'b, P: PrimeFieldParams> Add<&'b Fp<P>> for &'a Fp<P> {
    type Output = Fp<P>;

    fn add(self, other: &Fp<P>) -> Fp<P> {
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
}

impl<'a, 'b, P: PrimeFieldParams> Sub<&'b Fp<P>> for &'a Fp<P> {
    type Output = Fp<P>;

    fn sub(self, other: &Fp<P>) -> Fp<P> {
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
}

impl<'a, 'b, P: PrimeFieldParams> Mul<&'b Fp<P>> for &'a Fp<P> {
    type Output = Fp<P>;

    fn mul(self, other: &Fp<P>) -> Fp<P> {
        Fp {
            value: (&self.value * &other.value) % &P::modulus(),
            _marker: PhantomData
        }
    }
}

impl<'a, P: PrimeFieldParams> Neg for &'a Fp<P> {
    type Output = Fp<P>;

    fn neg(self) -> Fp<P> {
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
}

impl<P: PrimeFieldParams> Neg for Fp<P> {
    type Output = Fp<P>;

    fn neg(self) -> Fp<P> {
        -(&self)
    }
}

impl<P: PrimeFieldParams> PartialEq for Fp<P> {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl<P: PrimeFieldParams> Eq for Fp<P> {}

forward_all_binop_to_ref_ref!(impl(P: PrimeFieldParams) Mul for Fp<P>, mul);

#[cfg(test)]
mod large_field_tests {
    use super::*;
    use rand::{Rng,SeedableRng,StdRng};
    use num::{BigUint, Num};

    struct Small;

    impl PrimeFieldParams for Small {
        fn modulus() -> BigUint {
            BigUint::from_str_radix("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap()
        }

        fn bits() -> usize { 254 }
        fn name() -> &'static str { "Small" }
    }

    type Ft = Fp<Small>;

    #[test]
    fn rand_element_squaring() {
        let seed: [usize; 4] = [0, 0, 0, 0];
        let rng = &mut StdRng::from_seed(&seed);

        for _ in 0..100 {
            let a = Ft::random(rng);

            let mul = &a * &a;
            let sq = a.squared();

            assert!(sq == mul);
        }

        let mut cur = Ft::zero();
        for _ in 0..100 {
            let mul = &cur * &cur;
            let sq = cur.squared();

            assert!(sq == mul);

            cur = &cur + &Ft::one();
        }
    }

    #[test]
    fn rand_element_multiplication() {
        // If field is not associative under multiplication, 1/8 of all triplets a, b, c
        // will fail the test (a*b)*c = a*(b*c).

        let seed: [usize; 4] = [0, 0, 0, 0];
        let rng = &mut StdRng::from_seed(&seed);

        for _ in 0..250 {
            let a = &Ft::random(rng);
            let b = &Ft::random(rng);
            let c = &Ft::random(rng);

            assert!((a * b) * c == (b * c) * a);
        }
    }

    #[test]
    fn rand_element_inverse() {
        let seed: [usize; 4] = [0, 0, 0, 0];
        let rng = &mut StdRng::from_seed(&seed);

        for _ in 0..100 {
            let mut n = Ft::random(rng);
            n = n.inverse() * n;
            assert_eq!(n, Ft::one());
        }
        for _ in 0..100 {
            let a = Ft::random(rng);
            let b = Ft::random(rng);
            assert_eq!(&a * &b * (a.inverse()), b);
        }
    }

    #[test]
    fn bit_testing() {
        let a = Ft::from("13");
        assert!(a.test_bit(0) == true);
        assert!(a.test_bit(1) == false);
        assert!(a.test_bit(2) == true);
        assert!(a.test_bit(3) == true);

        let expected: Vec<bool> = [1,1,0,1,1,0,0,0,0,1,0,0,1,1,1,0,0,0,0,0,1,1,0,0,1,0,0,1,1]
                                  .iter().map(|a| *a == 1).rev().collect();

        let a = Ft::from("453624211");

        for (i, b) in expected.into_iter().enumerate() {
            assert!(a.test_bit(i) == b);
        }

        let expected: Vec<bool> = [1,1,1,1,0,1,0,1,1,0,1,0,0,0,1,1,1,0,1,1,1,1,0,0,0,0,1,1,0,1,1,0,1,0,0,1,1,0,1,1,0,0,1,1,0,0,1,1,1,1,0,1,1,1,1,0,1,0,1,1,1,1,1,0,1,0,1,0,0,0,1,1,1,0,1,1,1,1,0,0,0,1,1,1,0,0,1,0,1,0,0,0,0,1,1,0,0,1,1,1,1,0,1,1,1,1,0,0,0,1,1,0,0,1,1,1,0,1,1,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,1,0,1,0,1,0,0,0,1,1,0,1,1,0,0,1,0,1,1,1,1,1,0,0,0,1,0,1,0,0,0,1,0,1,0,0,1,0,1,1,0,1,0,1,1,0,1,0,0,0,1,0,0,1,0,1,1,1,0,1,1,0,0,0,0,0,1,1,0,0,0,1,1,1,1,1,1,0,0,0,0,1,0,1,1,1,0,1,1,0,1,1,0,0,0,0,1,1,1,1,1,0,0,1,1,1,1,1,1,0,1,0,0,1,0,1,1,0,0]
                                  .iter().map(|a| *a == 1).rev().collect();
        let a = Ft::from("13888242871869275222244405745257275088696211157297823662689037894645226208556");

        for (i, b) in expected.into_iter().enumerate() {
            assert!(a.test_bit(i) == b);
        }
    }
}

#[cfg(test)]
mod small_field_tests {
    use super::*;
    use num::{BigUint, Num};

    struct Small;

    impl PrimeFieldParams for Small {
        fn modulus() -> BigUint {
            BigUint::from_str_radix("13", 10).unwrap()
        }

        fn bits() -> usize { 6 }
        fn name() -> &'static str { "Small" }
    }

    type Ft = Fp<Small>;

    #[test]
    fn field_ops() {
        fn test_field_operation<C: Fn(&Ft, &Ft) -> Ft>(a: u64, b: u64, f: C, expected: u64) {
            let af = Ft::from(format!("{}", a).as_ref());
            let bf = Ft::from(format!("{}", b).as_ref());
            let expectedf = Ft::from(format!("{}", expected).as_ref());

            let res = f(&af, &bf);

            if res != expectedf {
                panic!("res={:?} != expectedf={:?} (a={}, b={}, expected={})", res, expectedf, a, b, expected);
            }
        }

        const MODULO: u64 = 13;

        for a in 0..13u64 {
            for b in 0..13u64 {
                test_field_operation(a, b, |a,b| {a * b}, (a*b)%MODULO);
                test_field_operation(a, b, |a,b| {a + b}, (a+b)%MODULO);
                test_field_operation(a, b, |a,b| {a - b}, {
                    let mut tmp = (a as i64) - (b as i64);
                    if tmp < 0 {
                        tmp += MODULO as i64;
                    }

                    tmp as u64
                });
                test_field_operation(a, b, |a,b| {a.pow(b)}, (a.pow(b as u32))%MODULO);
            }
            test_field_operation(a, 0, |a,_| {-a}, if a == 0 { 0 } else { MODULO - a });
            if a > 0 {
                test_field_operation(a, 0, |a,_| {&a.inverse() * a}, 1);
            }
        }
    }
}
