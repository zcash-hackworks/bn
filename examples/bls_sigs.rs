    extern crate rand;
    extern crate bn;

    use bn::{Field, Scalar, G1, G2, pairing};
    mod sighash;

    fn main() {
        let rng = &mut rand::thread_rng();

        // Generate Keys
        let alice_sk = Scalar::random(rng);

        // Generate Public Keys
        let alice_pk = G2::one() * &alice_sk;
        // Generate Signature
        let msgm1 = G1::random(&mut sighash::SignatureHash::from("Hello!"));
        let sigm1_a = &msgm1 * &alice_sk;

        // Verify single signature
        assert_eq!(pairing(&sigm1_a, &G2::one()), pairing(&msgm1, &alice_pk));

    }
