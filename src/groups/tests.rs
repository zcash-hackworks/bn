use super::GroupElement;
use fields::{FieldElement, Fr};
use rand::Rng;

fn random_test_addition<G: GroupElement, R: Rng>(rng: &mut R) {
    for _ in 0..50 {
        let r1 = G::random(rng);
        let r2 = G::random(rng);
        let r3 = G::random(rng);

        assert_eq!((r1 + r2) + r3, r1 + (r2 + r3));
        assert!(((r1 + r2 + r3) - r2 - r3 - r1).is_zero());
    }
}

fn random_test_doubling<G: GroupElement, R: Rng>(rng: &mut R) {
    for _ in 0..50 {
        let r1 = G::random(rng);
        let r2 = G::random(rng);
        let ti = Fr::from_str("2").unwrap().inverse().unwrap();

        assert_eq!((r1 + r2) + r1, r1.double() + r2);
        assert_eq!(r1, r1.double() * ti);
    }
}

fn random_test_dh<G: GroupElement, R: Rng>(rng: &mut R) {
    for _ in 0..50 {
        let alice_sk = Fr::random(rng);
        let bob_sk = Fr::random(rng);

        let alice_pk = G::one() * alice_sk;
        let bob_pk = G::one() * bob_sk;

        let alice_shared = bob_pk * alice_sk;
        let bob_shared = alice_pk * bob_sk;

        assert_eq!(alice_shared, bob_shared);
    }
}

fn random_test_equality<G: GroupElement, R: Rng>(rng: &mut R) {
    for _ in 0..50 {
        let begin = G::random(rng);

        let mut acc = begin;

        // Do a bunch of random things.

        let a = Fr::random(rng);
        let b = G::random(rng);
        let c = Fr::random(rng);
        let d = G::random(rng);

        for _ in 0..10 {
            acc = acc * a;
            acc = -acc;
            acc = acc + b;
            acc = acc * c;
            acc = -acc;
            acc = acc - d;
            acc = acc.double();
        }

        // Then reverse the operations

        let ai = a.inverse().unwrap();
        let ci = c.inverse().unwrap();
        let ti = Fr::from_str("2").unwrap().inverse().unwrap();

        for _ in 0..10 {
            acc = acc * ti;
            acc = acc + d;
            acc = -acc;
            acc = acc * ci;
            acc = acc - b;
            acc = -acc;
            acc = acc * ai;
        }

        assert_eq!(acc, begin);
    }
}

pub fn group_trials<G: GroupElement>() {
    assert!(G::zero().is_zero());
    assert!((G::one() - G::one()).is_zero());
    assert_eq!(G::one() + G::one(), G::one() * Fr::from_str("2").unwrap());
    assert!(G::zero().double().is_zero());

    assert!((G::one() * (-Fr::one()) + G::one()).is_zero());

    use rand::{SeedableRng,StdRng};
    let seed: [usize; 4] = [103245, 191922, 1293, 192103];
    let mut rng = StdRng::from_seed(&seed);

    random_test_addition::<G, _>(&mut rng);
    random_test_doubling::<G, _>(&mut rng);
    random_test_dh::<G, _>(&mut rng);
    random_test_equality::<G, _>(&mut rng);
}

