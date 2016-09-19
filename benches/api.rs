#![feature(test)]
extern crate test;
extern crate rand;
extern crate bn;
extern crate bincode;

use bn::*;
use bincode::SizeLimit::Infinite;
use bincode::rustc_serialize::{encode, decode};

#[bench]
fn g1_deserialization(b: &mut test::Bencher) {
    const SAMPLES: usize = 1000;

    let rng = &mut rand::thread_rng();

    let serialized: Vec<_> = (0..SAMPLES).map(|_| encode(&G1::random(rng), Infinite).unwrap()).collect();

    let mut ctr = 0;

    b.iter(|| {
        ctr += 1;

        decode::<G1>(&serialized[ctr % SAMPLES]).unwrap()
    });
}

#[bench]
fn g2_deserialization(b: &mut test::Bencher) {
    const SAMPLES: usize = 1000;

    let rng = &mut rand::thread_rng();

    let serialized: Vec<_> = (0..SAMPLES).map(|_| encode(&G2::random(rng), Infinite).unwrap()).collect();

    let mut ctr = 0;

    b.iter(|| {
        ctr += 1;

        decode::<G2>(&serialized[ctr % SAMPLES]).unwrap()
    });
}

#[bench]
fn fr_addition(b: &mut test::Bencher) {
    const SAMPLES: usize = 1000;

    let rng = &mut rand::thread_rng();

    let v1: Vec<_> = (0..SAMPLES).map(|_| Fr::random(rng)).collect();
    let v2: Vec<_> = (0..SAMPLES).map(|_| Fr::random(rng)).collect();

    let mut ctr = 0;

    b.iter(|| {
        ctr += 1;

        v1[ctr % SAMPLES] + v2[ctr % SAMPLES]
    });
}

#[bench]
fn fr_subtraction(b: &mut test::Bencher) {
    const SAMPLES: usize = 1000;

    let rng = &mut rand::thread_rng();

    let v1: Vec<_> = (0..SAMPLES).map(|_| Fr::random(rng)).collect();
    let v2: Vec<_> = (0..SAMPLES).map(|_| Fr::random(rng)).collect();

    let mut ctr = 0;

    b.iter(|| {
        ctr += 1;

        v1[ctr % SAMPLES] - v2[ctr % SAMPLES]
    });
}

#[bench]
fn fr_multiplication(b: &mut test::Bencher) {
    const SAMPLES: usize = 1000;

    let rng = &mut rand::thread_rng();

    let v1: Vec<_> = (0..SAMPLES).map(|_| Fr::random(rng)).collect();
    let v2: Vec<_> = (0..SAMPLES).map(|_| Fr::random(rng)).collect();

    let mut ctr = 0;

    b.iter(|| {
        ctr += 1;

        v1[ctr % SAMPLES] * v2[ctr % SAMPLES]
    });
}

#[bench]
fn fr_inverses(b: &mut test::Bencher) {
	const SAMPLES: usize = 1000;

    let rng = &mut rand::thread_rng();

    let v1: Vec<_> = (0..SAMPLES).map(|_| Fr::random(rng)).collect();

    let mut ctr = 0;

    b.iter(|| {
        ctr += 1;

        v1[ctr % SAMPLES].inverse()
    });
}

#[bench]
fn g1_addition(b: &mut test::Bencher) {
    const SAMPLES: usize = 100;

    let rng = &mut rand::thread_rng();

    let v1: Vec<_> = (0..SAMPLES).map(|_| G1::random(rng)).collect();
    let v2: Vec<_> = (0..SAMPLES).map(|_| G1::random(rng)).collect();

    let mut ctr = 0;

    b.iter(|| {
        ctr += 1;

        v1[ctr % SAMPLES] + v2[ctr % SAMPLES]
    });
}

#[bench]
fn g1_subtraction(b: &mut test::Bencher) {
    const SAMPLES: usize = 100;

    let rng = &mut rand::thread_rng();

    let v1: Vec<_> = (0..SAMPLES).map(|_| G1::random(rng)).collect();
    let v2: Vec<_> = (0..SAMPLES).map(|_| G1::random(rng)).collect();

    let mut ctr = 0;

    b.iter(|| {
        ctr += 1;

        v1[ctr % SAMPLES] - v2[ctr % SAMPLES]
    });
}

#[bench]
fn g1_scalar_multiplication(b: &mut test::Bencher) {
    const SAMPLES: usize = 100;

    let rng = &mut rand::thread_rng();

    let v1: Vec<_> = (0..SAMPLES).map(|_| G1::random(rng)).collect();
    let v2: Vec<_> = (0..SAMPLES).map(|_| Fr::random(rng)).collect();

    let mut ctr = 0;

    b.iter(|| {
        ctr += 1;

        v1[ctr % SAMPLES] * v2[ctr % SAMPLES]
    });
}

#[bench]
fn g2_addition(b: &mut test::Bencher) {
    const SAMPLES: usize = 100;

    let rng = &mut rand::thread_rng();

    let v1: Vec<_> = (0..SAMPLES).map(|_| G2::random(rng)).collect();
    let v2: Vec<_> = (0..SAMPLES).map(|_| G2::random(rng)).collect();

    let mut ctr = 0;

    b.iter(|| {
        ctr += 1;

        v1[ctr % SAMPLES] + v2[ctr % SAMPLES]
    });
}

#[bench]
fn g2_subtraction(b: &mut test::Bencher) {
    const SAMPLES: usize = 100;

    let rng = &mut rand::thread_rng();

    let v1: Vec<_> = (0..SAMPLES).map(|_| G2::random(rng)).collect();
    let v2: Vec<_> = (0..SAMPLES).map(|_| G2::random(rng)).collect();

    let mut ctr = 0;

    b.iter(|| {
        ctr += 1;

        v1[ctr % SAMPLES] - v2[ctr % SAMPLES]
    });
}

#[bench]
fn g2_scalar_multiplication(b: &mut test::Bencher) {
    const SAMPLES: usize = 100;

    let rng = &mut rand::thread_rng();

    let v1: Vec<_> = (0..SAMPLES).map(|_| G2::random(rng)).collect();
    let v2: Vec<_> = (0..SAMPLES).map(|_| Fr::random(rng)).collect();

    let mut ctr = 0;

    b.iter(|| {
        ctr += 1;

        v1[ctr % SAMPLES] * v2[ctr % SAMPLES]
    });
}

#[bench]
fn fq12_scalar_multiplication(b: &mut test::Bencher) {
    const SAMPLES: usize = 100;

    let rng = &mut rand::thread_rng();

    let v1: Vec<_> = (0..SAMPLES).map(|_| G1::random(rng)).collect();
    let v2: Vec<_> = (0..SAMPLES).map(|_| G2::random(rng)).collect();
    let v3: Vec<_> = (0..SAMPLES).map(|i| pairing(v1[i], v2[i])).collect();

    let mut ctr = 0;

    b.iter(|| {
        ctr += 1;

        v3[(ctr + SAMPLES/50) % SAMPLES] * v3[ctr % SAMPLES]
    });
}

#[bench]
fn fq12_exponentiation(b: &mut test::Bencher) {
    const SAMPLES: usize = 100;

    let rng = &mut rand::thread_rng();

    let v1: Vec<_> = (0..SAMPLES).map(|_| G1::random(rng)).collect();
    let v2: Vec<_> = (0..SAMPLES).map(|_| G2::random(rng)).collect();
    let v3: Vec<_> = (0..SAMPLES).map(|i| pairing(v1[i], v2[i])).collect();
    let v4: Vec<_> = (0..SAMPLES).map(|_| Fr::random(rng)).collect();

    let mut ctr = 0;

    b.iter(|| {
        ctr += 1;

        v3[ctr % SAMPLES].pow(v4[ctr % SAMPLES])
    });
}

#[bench]
fn perform_pairing(b: &mut test::Bencher) {
    const SAMPLES: usize = 100;

    let rng = &mut rand::thread_rng();

    let v1: Vec<_> = (0..SAMPLES).map(|_| G1::random(rng)).collect();
    let v2: Vec<_> = (0..SAMPLES).map(|_| G2::random(rng)).collect();

    let mut ctr = 0;

    b.iter(|| {
        ctr += 1;

        pairing(v1[ctr % SAMPLES], v2[ctr % SAMPLES])
    });
}
