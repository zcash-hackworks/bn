extern crate bn;
extern crate rand;
extern crate sodiumoxide;

mod sighash;

use bn::*;
use sodiumoxide::crypto::stream::chacha20;
use sodiumoxide::crypto::hash::sha256;
use std::str;

fn main() {
  ibe();
}

fn ibe() {
    let rng       = &mut rand::thread_rng();
    let master_sk = Scalar::random(rng);
    // do we need another generator than G1::one() here?
    // we use G1 since Ppub is used in the first arg of pairing
    let master_pk = G1::one() * &master_sk;

    let id = b"test";

    let derived = G2::random(&mut sighash::SignatureHash::from(&id[..]));
    println!("derived: {:?}", derived);

    let id_sk = &derived * &master_sk;

    //encrypting with BasicIdent
    let r = Scalar::random(rng);
    let g_id = pairing(&master_pk, &derived) ^ &r;
    println!("g_id: {:?}", g_id);
    let badly_serialized = format!("{:?}", g_id);
    let hash = sha256::hash(badly_serialized.as_bytes());

    println!("hash: {:?}", hash);
    let sym_key = chacha20::Key::from_slice(&hash[..32]).unwrap();
    let nonce   = chacha20::gen_nonce();

    let plaintext  = "We propose a fully functional identity-based encryption scheme (IBE). The scheme has chosen ciphertext security in the random oracle model assuming a variant of the computational Diffie-Hellman problem";
    let ciphertext = chacha20::stream_xor(plaintext.as_bytes(), &nonce, &sym_key);

    // do we need another generator than G1::one() here?
    let result = (G1::one() * &r, ciphertext);
    println!("ciphertext: {:?}", result);

    //decrypting
    let decrypting_seed = pairing(&result.0, &id_sk);
    let badly_serialized_again = format!("{:?}", decrypting_seed);
    println!("seed: {:?}", decrypting_seed);
    assert_eq!(g_id, decrypting_seed);
    let hash2 = sha256::hash(badly_serialized_again.as_bytes());
    println!("hash2: {:?}", hash2);
    let sym_key_2 = chacha20::Key::from_slice(&hash2[..32]).unwrap();
    let decrypted = chacha20::stream_xor(&result.1, &nonce, &sym_key_2);
    println!("decrypted: \"{}\"", str::from_utf8(&decrypted).unwrap());
    assert_eq!(plaintext.as_bytes(), &decrypted[..]);

}

