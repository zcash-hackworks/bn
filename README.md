# bn [![Crates.io](https://img.shields.io/crates/v/bn.svg)](https://crates.io/crates/bn) [![Build status](https://api.travis-ci.org/zcash/bn.svg)](https://travis-ci.org/zcash/bn)

This is a [pairing cryptography](https://en.wikipedia.org/wiki/Pairing-based_cryptography) library written in pure Rust. It makes use of the Barreto-Naehrig (BN) curve construction from [[BCTV2015]](https://eprint.iacr.org/2013/879.pdf) to provide two cyclic groups **G<sub>1</sub>** and **G<sub>2</sub>**, with an efficient bilinear pairing:

*e: G<sub>1</sub> × G<sub>2</sub> → G<sub>T</sub>*

## Security warnings

This library, like other pairing cryptography libraries implementing this construction, is not resistant to side-channel attacks.

## Usage

Add the `bn` crate to your dependencies in `Cargo.toml`...

```toml
[dependencies]
bn = "0.4.3"
```

...and add an `extern crate` declaration to your crate root:

```rust
extern crate bn;
```

## API

* `Fr` is an element of F<sub>r</sub>
* `G1` is a point on the BN curve E/Fq : y^2 = x^3 + b
* `G2` is a point on the twisted BN curve E'/Fq2 : y^2 = x^3 + b/xi
* `Gt` is a group element (written multiplicatively) obtained with the `pairing` function over `G1` and `G2`.

### Examples

#### Joux's key agreement protocol

In a typical Diffie-Hellman key exchange, relying on ECDLP, a three-party key exchange requires two rounds. A single round protocol is possible through the use of a bilinear pairing: given Alice's public key *a*P<sub>1</sub> and Bob's public key *b*P<sub>2</sub>, Carol can compute the shared secret with her private key *c* by *e*(*a*P<sub>1</sub>, *b*P<sub>2</sub>)<sup>c</sup>.

(See `examples/joux.rs` for the full example.)

```rust
// Generate private keys
let alice_sk = Fr::random(rng);
let bob_sk = Fr::random(rng);
let carol_sk = Fr::random(rng);

// Generate public keys in G1 and G2
let (alice_pk1, alice_pk2) = (G1::one() * alice_sk, G2::one() * alice_sk);
let (bob_pk1, bob_pk2) = (G1::one() * bob_sk, G2::one() * bob_sk);
let (carol_pk1, carol_pk2) = (G1::one() * carol_sk, G2::one() * carol_sk);

// Each party computes the shared secret
let alice_ss = pairing(bob_pk1, carol_pk2).pow(alice_sk);
let bob_ss = pairing(carol_pk1, alice_pk2).pow(bob_sk);
let carol_ss = pairing(alice_pk1, bob_pk2).pow(carol_sk);

assert!(alice_ss == bob_ss && bob_ss == carol_ss);
```

## License

Licensed under either of

 * MIT license, ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)
 * Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)

at your option.

Copyright 2016 [Zcash Electric Coin Company](https://z.cash/). The Zcash Company promises to maintain the "bn" crate on crates.io under this MIT/Apache-2.0 dual license.

### Authors

* [Sean Bowe](https://github.com/ebfull)

### Contribution

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in the work by you, as defined in the Apache-2.0
license, shall be dual licensed as above, without any additional terms or
conditions.
