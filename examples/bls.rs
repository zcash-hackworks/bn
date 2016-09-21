    extern crate rand;
    extern crate bn;
    extern crate ring;
    #[macro_use]
    extern crate arrayref;


    use bn::{Group,Fr, G1, G2, pairing};
    use ring::digest;

    fn main() {
        let rng = &mut rand::thread_rng();

        // Generate Keys
        let alice_sk = Fr::random(rng);

        // Generate Public Keys
        let alice_pk = G2::one() * alice_sk;
        // Generate Signature
        let mut ctx = digest::Context::new(&digest::SHA512);
        ctx.update(b"Hello!");
        let hash = ctx.finish();
        let msgm1 =Fr::interpret(array_ref![hash.as_ref(),0,64]);
        let sigm1_a = msgm1 * alice_sk;

        // Verify single signature
        assert!(pairing(G1::one() * sigm1_a, G2::one())== pairing(G1::one() * msgm1, alice_pk));

    }