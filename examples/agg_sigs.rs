
extern crate rand;
extern crate bn;

use bn::{Field, Scalar, G1, G2, Gt, Fq12, pairing};
mod sighash;
use std::collections::HashMap;


fn verify_aggregate_sig(agg_sig: &G1, msg_key_pairs: &[(&str, G2)]) -> bool {
    let mut unique = HashMap::new();
    let mut agg_verifier = Gt::new(Fq12::one());
    for &(msg, ref pub_key) in msg_key_pairs {
        match unique.get(msg) {
            Some(&true) => return false, //fail on duplicate messages
            _ => {} // do nothing
        }
        unique.insert(msg, true);
        let msg_e = G1::random(&mut sighash::SignatureHash::from(msg));
        agg_verifier = agg_verifier * pairing(&msg_e, &pub_key);
    }
    if pairing(agg_sig, &G2::one()) != agg_verifier {
        return false;
    }
    return true;
}

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

    // Verify he Aggregate Signature
    assert!(verify_aggregate_sig(&sig_m1m2, &[("Hello!", alice_pk), ("Hello!2", bob_pk)]));

}
