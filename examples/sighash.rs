extern crate rand;
extern crate sodiumoxide;

pub struct SignatureHash(rand::ChaChaRng);

impl rand::Rng for SignatureHash {
    fn next_u32(&mut self) -> u32 {
        self.0.next_u32()
    }
}

impl<'a> From<&'a str> for SignatureHash {
    fn from(v: &str) -> SignatureHash {
        SignatureHash::from(v.as_bytes())
    }
}

impl<'a> From<&'a [u8]> for SignatureHash {
    fn from(v: &[u8]) -> SignatureHash {
        use rand::SeedableRng;
        use std::slice;
        use std::mem;

        let hash = sodiumoxide::crypto::hash::sha256::hash(v);
        assert_eq!(hash.0.len(), 32);

        SignatureHash(rand::ChaChaRng::from_seed(unsafe {
            slice::from_raw_parts(mem::transmute::<&u8, &u32>(&hash.0[0]), 8)
        }))
    }
}
