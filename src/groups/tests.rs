use super::*;
use ::Fr;
use fields::Field;

use rand::Rng;

pub fn group_trials<P: GroupParams>() {
    fn test_doubling<P: GroupParams>() {
        type G<P> = Jacobian<P>;

        let one = G::<P>::one();
        let two = one.add(&one);
        let three = two.add(&one);
        let four = three.add(&one);

        assert_eq!(one.double(), two);
        assert_eq!(two.double(), four);
    }

    fn random_test_addition<P: GroupParams, R: Rng>(rng: &mut R) {
        type G<P> = Jacobian<P>;

        for _ in 0..50 {
            let r1 = &G::<P>::random(rng);
            let r2 = &G::<P>::random(rng);
            let r3 = &G::<P>::random(rng);

            let s1 = (r1 + r2) + r3;
            let s2 = (r2 + r3) + r1;

            assert_eq!(s1, s2);
        }
    }
    
    fn random_test_doubling<P: GroupParams, R: Rng>(rng: &mut R) {
        type G<P> = Jacobian<P>;

        for _ in 0..50 {
            let r1 = &G::<P>::random(rng);
            let r2 = &G::<P>::random(rng);

            let a = (r1 + r2) + r1;
            let b = r1.double() + r2;

            assert!(a.eq(&b));
        }
    }

    fn random_test_dh<P: GroupParams, R: Rng>(rng: &mut R) {
        type G<P> = Jacobian<P>;

        for _ in 0..50 {
            let alice_sk = Fr::random(rng);
            let bob_sk = Fr::random(rng);

            let alice_pk = G::<P>::one() * &alice_sk;
            let bob_pk = G::<P>::one() * &bob_sk;

            let alice_shared = &bob_pk * &alice_sk;
            let bob_shared = &alice_pk * &bob_sk;

            assert_eq!(alice_shared, bob_shared);
        }
    }

    test_doubling::<P>();

    use rand::{SeedableRng,StdRng};
    let seed: [usize; 4] = [103245, 191922, 1293, 192103];
    let mut rng = StdRng::from_seed(&seed);

    random_test_addition::<P, _>(&mut rng);
    random_test_doubling::<P, _>(&mut rng);
    random_test_dh::<P, _>(&mut rng);
}

