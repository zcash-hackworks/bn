extern crate rand;
extern crate rustc_serialize;
extern crate byteorder;

mod arith;
mod fields;
mod groups;

pub use arith::U256;
pub use fields::{Fq, Fr, Fq2, Fq6, Fq12, FieldElement};
pub use groups::{G1, G2, GroupElement};

pub fn pairing(p: &G1, q: &G2) -> Fq12 {
	match (p.to_affine(), q.to_affine()) {
		(None, _) | (_, None) => Fq12::one(),
		(Some(p), Some(q)) => {
			q.precompute().miller_loop(&p).final_exponentiation().expect("miller loop cannot produce zero")
		}
	}
}

#[test]
fn test_reduced_pairing() {
    let g1 = G1::one() * Fr::from_str("18097487326282793650237947474982649264364522469319914492172746413872781676").unwrap();
    let g2 = G2::one() * Fr::from_str("20390255904278144451778773028944684152769293537511418234311120800877067946").unwrap();

    let gt = pairing(&g1, &g2);

    let expected = Fq12::new(
        Fq6::new(
            Fq2::new(Fq::from_str("7520311483001723614143802378045727372643587653754534704390832890681688842501").unwrap(), Fq::from_str("20265650864814324826731498061022229653175757397078253377158157137251452249882").unwrap()),
            Fq2::new(Fq::from_str("11942254371042183455193243679791334797733902728447312943687767053513298221130").unwrap(), Fq::from_str("759657045325139626991751731924144629256296901790485373000297868065176843620").unwrap()),
            Fq2::new(Fq::from_str("16045761475400271697821392803010234478356356448940805056528536884493606035236").unwrap(), Fq::from_str("4715626119252431692316067698189337228571577552724976915822652894333558784086").unwrap())
        ),
        Fq6::new(
            Fq2::new(Fq::from_str("14901948363362882981706797068611719724999331551064314004234728272909570402962").unwrap(), Fq::from_str("11093203747077241090565767003969726435272313921345853819385060670210834379103").unwrap()),
            Fq2::new(Fq::from_str("17897835398184801202802503586172351707502775171934235751219763553166796820753").unwrap(), Fq::from_str("1344517825169318161285758374052722008806261739116142912817807653057880346554").unwrap()),
            Fq2::new(Fq::from_str("11123896897251094532909582772961906225000817992624500900708432321664085800838").unwrap(), Fq::from_str("17453370448280081813275586256976217762629631160552329276585874071364454854650").unwrap())
        )
    );

    assert_eq!(expected, gt);
}

#[test]
fn test_binlinearity() {
    use rand::{SeedableRng,StdRng};
    let seed: [usize; 4] = [103245, 191922, 1293, 192103];
    let mut rng = StdRng::from_seed(&seed);

    for _ in 0..50 {
        let p = G1::random(&mut rng);
        let q = G2::random(&mut rng);
        let s = Fr::random(&mut rng);
        let sp = p * s;
        let sq = q * s;

        let a = pairing(&p, &q).pow(s);
        let b = pairing(&sp, &q);
        let c = pairing(&p, &sq);

        assert_eq!(a, b);
        assert_eq!(b, c);

        let t = -Fr::one();

        assert!(a != Fq12::one());
        assert_eq!((a.pow(t)) * a, Fq12::one());
    }
}
