extern crate bn;

use bn::*;

fn main() {
    let mut a = G1::one();
    let mut b = G2::one();
    let c = Fr::from_str("1901").unwrap().inverse().unwrap();
    let d = Fr::from_str("2344").unwrap().inverse().unwrap();

    let mut acc1 = Gt::one();
    for i in 0..10000 {
        acc1 = acc1 * pairing(a, b);
        a = a * c;
        b = b * d;
    }

    a = G1::one();
    b = G2::one();

    let mut acc2 = Gt::one();
    for i in 0..10000 {
        acc2 = acc2 * pairing(a, b);
        a = a * d;
        b = b * c;
    }

    assert!(acc1 == acc2);
}
