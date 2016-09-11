// This is an example of three-party Diffie-Hellman key exchange
// Requires two rounds

extern crate bn;
extern crate rand;

use bn::*;

fn main() {
    let rng = &mut rand::thread_rng();

    // Construct private keys
    let alice_sk = Fr::random(rng);
    let bob_sk = Fr::random(rng);
    let carol_sk = Fr::random(rng);

    // Construct public keys
    let alice_pk = G1::one() * alice_sk;
    let bob_pk = G1::one() * bob_sk;
    let carol_pk = G1::one() * carol_sk;

    // Round one:
    let alice_dh_1 = bob_pk * carol_sk;
    let bob_dh_1 = carol_pk * alice_sk;
    let carol_dh_1 = alice_pk * bob_sk;

    // Round two:
    let alice_dh_2 = alice_dh_1 * alice_sk;
    let bob_dh_2 = bob_dh_1 * bob_sk;
    let carol_dh_2 = carol_dh_1 * carol_sk;

    // All parties should arrive to the same shared secret
    assert!(alice_dh_2 == bob_dh_2 && bob_dh_2 == carol_dh_2);
}
