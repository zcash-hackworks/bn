
extern crate rand;
extern crate bn;

use bn::{Field, Scalar, G1, G2, pairing};
mod sighash;

fn main() {
    let rng = &mut rand::thread_rng();

    // Generate Keys
    let alice_sk = Scalar::random(rng);
    let bob_sk = Scalar::random(rng);

    // Generate Public Keys
    let alice_pk = G2::one() * &alice_sk;
    let bob_pk = G2::one() * &bob_sk;
    // Generate Signatures
    let msgm1 = G1::random(&mut sighash::SignatureHash::from("Hello!"));
    let sigm1_a = &msgm1 * &alice_sk;

    let msgm2 = G1::random(&mut sighash::SignatureHash::from("Hello!2"));
    let sigm2_b = &msgm2 * &bob_sk;

    // Verify single signatures
    assert_eq!(pairing(&sigm1_a, &G2::one()), pairing(&msgm1, &alice_pk));
    assert_eq!(pairing(&sigm2_b, &G2::one()), pairing(&msgm2, &bob_pk));

    // Generate Aggregate Signature
    let sig_m1m2 = &sigm1_a + &sigm2_b;
    // Show the messages are distinct
    assert!(msgm2 != msgm1);
    // Verify he Aggregate Signature
    assert_eq!(pairing(&sig_m1m2, &G2::one()),
               pairing(&msgm1, &alice_pk) * pairing(&msgm2, &bob_pk));

}
