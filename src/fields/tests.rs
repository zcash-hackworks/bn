use rand::{Rng,SeedableRng,StdRng};
use fields::Field;

mod large_field {
    use fields::fp::*;
    use num::{BigUint, Num};

    struct Large;

    impl PrimeFieldParams for Large {
        fn modulus() -> BigUint {
            BigUint::from_str_radix("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap()
        }

        fn bits() -> usize { 254 }
        fn name() -> &'static str { "Large" }
    }

    type Ft = Fp<Large>;

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

mod small_field {
    use fields::fp::*;
    use fields::Field;
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


fn can_invert<F: Field>() {
    let mut a = F::one();
    for _ in 0..1000 {
        assert!(a.ne(&F::zero()));

        let inv = a.inverse();

        assert!(a.mul(&inv).eq(&F::one()));

        a = a.add(&F::one());
    }
}


fn rand_element_squaring<F: Field, R: Rng>(rng: &mut R) {
    for _ in 0..100 {
        let a = F::random(rng);

        let mul = a.mul(&a);
        let sq = a.squared();

        assert!(sq.eq(&mul));
    }

    let mut cur = F::zero();
    for _ in 0..100 {
        let mul = cur.mul(&cur);
        let sq = cur.squared();

        assert!(sq.eq(&mul));

        cur = cur.add(&F::one());
    }
}

fn rand_element_addition_and_negation<F: Field, R: Rng>(rng: &mut R) {
    for _ in 0..10 {
        let mut a = F::random(rng);
        let r = F::random(rng);
        let mut b = a.add(&r);

        for _ in 0..10 {
            let r = F::random(rng);
            a = a.add(&r);
            b = b.add(&r);

            let r = F::random(rng);
            a = a.sub(&r);
            b = b.sub(&r);

            let r = F::random(rng);
            a = a.add(&r);
            b = b.add(&r);
        }

        b = b.sub(&r);
        assert!(a.eq(&b));
    }
}

fn rand_element_inverse<F: Field, R: Rng>(rng: &mut R) {
    for _ in 0..100 {
        let mut n = F::random(rng);
        n = n.inverse().mul(&n);
        assert!(n.eq(&F::one()));
    }
    for _ in 0..100 {
        let a = F::random(rng);
        let b = F::random(rng);
        assert!(a.mul(&b).mul(&a.inverse()).eq(&b));
    }
}

fn rand_element_multiplication<F: Field, R: Rng>(rng: &mut R) {
    // If field is not associative under multiplication, 1/8 of all triplets a, b, c
    // will fail the test (a*b)*c = a*(b*c).

    for _ in 0..250 {
        let a = F::random(rng);
        let b = F::random(rng);
        let c = F::random(rng);

        assert!(a.mul(&b).mul(&c).eq(&b.mul(&c).mul(&a)));
    }
}

pub fn field_trials<F: Field>() {
    can_invert::<F>();

    let seed: [usize; 4] = [103245, 191922, 1293, 192103];
    let mut rng = StdRng::from_seed(&seed);

    rand_element_squaring::<F, StdRng>(&mut rng);
    rand_element_addition_and_negation::<F, StdRng>(&mut rng);
    rand_element_multiplication::<F, StdRng>(&mut rng);
    rand_element_inverse::<F, StdRng>(&mut rng);
}


