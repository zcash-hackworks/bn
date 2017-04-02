extern crate rand;
extern crate bn;
extern crate ring;
#[macro_use]
extern crate arrayref;

use bn::{Group,Fr, G1, G2, Gt, pairing};
use ring::digest;
use std::collections::HashMap;


fn verify_aggregate_sig(agg_sig: G1, msg_key_pairs: &[(&str, G2)]) -> bool {
    let mut unique = HashMap::new();
    let mut agg_verifier = Gt::one();

    for &(msg, pub_key) in msg_key_pairs {
        match unique.get(msg) {
            Some(&true) => return false, //fail on duplicate messages
            _ => {} // do nothing
        }
        unique.insert(msg, true);
        let mut ctx = digest::Context::new(&digest::SHA512);
        ctx.update(msg.as_bytes());
        let hash = ctx.finish();
        let msg_e = Fr::interpret(array_ref![hash.as_ref(),0,64]);

        agg_verifier = agg_verifier * pairing(G1::one() * msg_e, pub_key);
    }
    if pairing(agg_sig, G2::one()) != agg_verifier {
        return false;
    }
    return true;
}

fn main() {
    let rng = &mut rand::thread_rng();

    const MSG1: &'static str = "Hello!";
    const MSG2: &'static str = "Hello!2";

    // Generate Keys
    let alice_sk = Fr::random(rng);
    let bob_sk = Fr::random(rng);

    // Generate Public Keys
    let alice_pk = G2::one() * alice_sk;
    let bob_pk = G2::one() * bob_sk;
    // Generate Signatures
    let mut ctx = digest::Context::new(&digest::SHA512);
    ctx.update(MSG1.as_bytes());
    let hash1 = ctx.finish();

    let msgm1 = Fr::interpret(array_ref![hash1.as_ref(),0,64]);
    let sigm1_a = msgm1 * alice_sk;

    let mut ctx2 = digest::Context::new(&digest::SHA512);
    ctx2.update(MSG2.as_bytes());
    let hash2 = ctx2.finish();

    let msgm2 = Fr::interpret(array_ref![hash2.as_ref(),0,64]);
    let sigm2_b = msgm2 * bob_sk;

    // Verify single signatures
    assert!(pairing(G1::one() * sigm1_a, G2::one()) == pairing(G1::one() * msgm1, alice_pk));
    assert!(pairing(G1::one() * sigm2_b, G2::one()) == pairing(G1::one() * msgm2, bob_pk));

    // Generate Aggregate Signature
    let sig_m1m2 = sigm1_a + sigm2_b;

    // Verify the Aggregate Signature
    assert!(verify_aggregate_sig(G1::one() * sig_m1m2, &[(MSG1, alice_pk), (MSG2, bob_pk)]));

    //Test duplicate messages
    //Generate bob's sig of MSG
    let sigm1_b = msgm1 * bob_sk;

    //Generate duplicate aggregate signature
    let sig_m1_dup = sigm1_a + sigm1_b;
    assert!(!verify_aggregate_sig(G1::one() * sig_m1_dup, &[(MSG1, alice_pk), (MSG1, bob_pk)]));

}