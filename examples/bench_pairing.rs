extern crate rand;
extern crate bn;
use bn::{Group, Fr, G1, G2, Gt, pairing};
use std::time::Instant;

fn main() {
    let test_size = 1000;
    let mut random_g1_a: Vec<Fr> = vec![];
    let mut random_g2_b: Vec<Fr> = vec![];
    let mut gt_compute: Vec<Gt> = vec![];
    // let mut pairing_res: Vec<Gt> = vec![];
    let gt_generator = pairing(G1::one(), G2::one());
    let rng = &mut rand::thread_rng();

    // Generate random field elements
    for _ in 0..test_size {
        let a = Fr::random(rng);
        let b = Fr::random(rng);
        random_g1_a.push(a);
        random_g2_b.push(b);
        gt_compute.push(gt_generator.pow(a * b));
    }

    // Generate public keys in G1 and G2
    let before = Instant::now();
    for i in 0..test_size {
        // pairing_res.push(pairing(G1::one() * random_g1_a[i], G2::one() * random_g2_b[i]));
        let res = pairing(G1::one() * random_g1_a[i], G2::one() * random_g2_b[i]);
        // assert!(gt_compute[i] == pairing_res[i]);
        assert!(gt_compute[i] == res);
    }
    println!("{} times pairing need {}ms", test_size, 
        before.elapsed().as_secs() * 1000 as u64 + before.elapsed().subsec_millis() as u64);

    
}