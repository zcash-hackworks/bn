use num::{Num,BigUint};
use fields::Field;
use fields::fp::PrimeFieldParams;
use fields::fp2::Fp2Params;
use fields::fp6::Fp6Params;
use fields::fp12::Fp12Params;
use super::{Fr,Fq,Fq2,Fq6,Fq12,G1,G2};
use groups::*;

pub struct FrParams;

impl PrimeFieldParams for FrParams {
    fn modulus() -> BigUint {
        BigUint::from_str_radix("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap()
    }
    fn bits() -> usize { 254 }
    fn name() -> &'static str { "Fr" }
}

pub struct FqParams;

impl PrimeFieldParams for FqParams {
    fn modulus() -> BigUint {
        BigUint::from_str_radix("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap()
    }
    fn bits() -> usize { 254 }
    fn name() -> &'static str { "Fq" }
}

#[test]
fn test_fr() {
    use fields;

    fields::tests::field_trials::<Fr>();
}

#[test]
fn test_fq() {
    use fields;

    fields::tests::field_trials::<Fq>();
}

pub struct G1Params;

impl GroupParams for G1Params {
    type Base = Fq;

    fn name() -> &'static str {
        "G1"
    }
    fn zero() -> Jacobian<Self> {
        Jacobian::new(Fq::zero(), Fq::one(), Fq::zero()).unwrap()
    }
    fn one() -> Jacobian<Self> {
        Jacobian::new(Fq::from("1"), Fq::from("2"), Fq::one()).unwrap()
    }
    fn coeff_b() -> Self::Base {
        Fq::from("3")
    }
}

#[test]
fn test_g1() {
    use groups;

    groups::tests::group_trials::<G1Params>();

    assert_eq!(G1::zero(), G1::one() + (G1::zero() - G1::one()));
    assert_eq!(G1::zero(), G1::one() * Fr::from("21888242871839275222246405745257275088548364400416034343698204186575808495616") + G1::one());
}

#[test]
fn g1_test_vector() {
    let a = G1::one() * Fr::from("19797905000333868150253315089095386158892526856493194078073564469188852136946");
    let b = G1::one() * Fr::from("2730506433347642574983433139433778984782882168213690554721050571242082865799");
    let e = &a + &b;

    let expect = G1::new(
        Fq::from("18450621724990678172567114131642278789161361170999664461184794604011563728206"),
        Fq::from("21329688341674583036384007811166435666174342925504675855816423131698588368496"),
        Fq::one()
    ).unwrap();

    assert_eq!(expect, e);
}

pub struct Fq2Params;

impl Fp2Params for Fq2Params {
    fn non_residue() -> Fq {
        Fq::from("21888242871839275222246405745257275088696311157297823662689037894645226208582")
    }
    fn name() -> &'static str {
        "Fq2"
    }
}

#[test]
fn test_fq2() {
    use fields;

    fields::tests::field_trials::<Fq2>();
}

pub struct G2Params;

impl G2Params {
    pub fn twist() -> Fq2 {
        Fq2::new(Fq::from("9"), Fq::from("1"))
    }
}

impl GroupParams for G2Params {
    type Base = Fq2;

    fn name() -> &'static str {
        "G2"
    }
    fn zero() -> Jacobian<Self> {
        Jacobian::new(Fq2::zero(), Fq2::one(), Fq2::zero()).unwrap()
    }
    fn one() -> Jacobian<Self> {
        Jacobian::new(
            Fq2::new(
                Fq::from("10857046999023057135944570762232829481370756359578518086990519993285655852781"),
                Fq::from("11559732032986387107991004021392285783925812861821192530917403151452391805634")
            ),
            Fq2::new(
                Fq::from("8495653923123431417604973247489272438418190587263600148770280649306958101930"),
                Fq::from("4082367875863433681332203403145435568316851327593401208105741076214120093531")
            ),
            Fq2::one()
        ).unwrap()
    }
    fn coeff_b() -> Self::Base {
        &G2Params::twist().inverse() * &Fq::from("3")
    }
}

#[test]
fn test_g2() {
    use groups;

    groups::tests::group_trials::<G2Params>();

    assert_eq!(G2::zero(), G2::one() + (G2::zero() - G2::one()));
    assert_eq!(G2::zero(), G2::one() * Fr::from("21888242871839275222246405745257275088548364400416034343698204186575808495616") + G2::one());
}

#[test]
fn g2_test_vector() {
    let a = G2::one() * Fr::from("19797905000333868150253315089095386158892526856493194078073564469188852136946");
    let b = G2::one() * Fr::from("2730506433347642574983433139433778984782882168213690554721050571242082865799");   
    let e = &a + &b;

    let expect = G2::new(
        Fq2::new(
            Fq::from("10805137482603266627116066166226222153808813611856467496561473491230213987197"),
            Fq::from("11018998371825437935082073888099464993330606622517843684670450190973893289235")
        ),
        Fq2::new(
            Fq::from("371699491666579792038680273553261511891341995868329474144713691525212078012"),
            Fq::from("2123259504314265904107110265140842273706723557882599408954283209162529085097")
        ),
        Fq2::one()
    ).unwrap();

    assert_eq!(expect, e);
}

pub struct Fq6Params;

impl Fp6Params for Fq6Params {
    fn non_residue() -> Fq2 {
        Fq2::new(Fq::from("9"), Fq::from("1"))
    }
    fn name() -> &'static str {
        "Fq6"
    }
}

#[test]
fn test_fq6() {
    use fields;

    fields::tests::field_trials::<Fq6>();
}

pub struct Fq12Params;

impl Fp12Params for Fq12Params {
    fn non_residue() -> Fq2 {
        Fq2::new(Fq::from("9"), Fq::from("1"))
    }
    fn name() -> &'static str {
        "Fq12"
    }
}

#[test]
fn test_fq12() {
    use fields;

    fields::tests::field_trials::<Fq12>();
}

#[test]
fn fq12_test_vector() {
    let start = Fq12::new(
        Fq6::new(
            Fq2::new(
                Fq::from("19797905000333868150253315089095386158892526856493194078073564469188852136946"),
                Fq::from("10509658143212501778222314067134547632307419253211327938344904628569123178733")
            ),
            Fq2::new(
                Fq::from("208316612133170645758860571704540129781090973693601051684061348604461399206"),
                Fq::from("12617661120538088237397060591907161689901553895660355849494983891299803248390")
            ),
            Fq2::new(
                Fq::from("2897490589776053688661991433341220818937967872052418196321943489809183508515"),
                Fq::from("2730506433347642574983433139433778984782882168213690554721050571242082865799")
            )
        ),
        Fq6::new(
            Fq2::new(
                Fq::from("17870056122431653936196746815433147921488990391314067765563891966783088591110"),
                Fq::from("14314041658607615069703576372547568077123863812415914883625850585470406221594")
            ),
            Fq2::new(
                Fq::from("10123533891707846623287020000407963680629966110211808794181173248765209982878"),
                Fq::from("5062091880848845693514855272640141851746424235009114332841857306926659567101")
            ),
            Fq2::new(
                Fq::from("9839781502639936537333620974973645053542086898304697594692219798017709586567"),
                Fq::from("1583892292110602864638265389721494775152090720173641072176370350017825640703")
            )
        )
    );

    // Do a bunch of arbitrary stuff to the element

    let mut next = start.clone();
    for _ in 0..100 {
        next = &next * &start;
    }

    let mut cpy = next.clone();

    for _ in 0..10 {
        next = next.squared();
    }

    for _ in 0..10 {
        next = &next + &start;
        next = &next - &cpy;
        next = -&next;
    }

    next = next.squared();

    let finally = Fq12::new(
        Fq6::new(
            Fq2::new(
                Fq::from("18388750939593263065521177085001223024106699964957029146547831509155008229833"),
                Fq::from("18370529854582635460997127698388761779167953912610241447912705473964014492243")
            ),
            Fq2::new(
                Fq::from("3691824277096717481466579496401243638295254271265821828017111951446539785268"),
                Fq::from("20513494218085713799072115076991457239411567892860153903443302793553884247235")
            ),
            Fq2::new(
                Fq::from("12214155472433286415803224222551966441740960297013786627326456052558698216399"),
                Fq::from("10987494248070743195602580056085773610850106455323751205990078881956262496575")
            )
        ),
        Fq6::new(
            Fq2::new(
                Fq::from("5134522153456102954632718911439874984161223687865160221119284322136466794876"),
                Fq::from("20119236909927036376726859192821071338930785378711977469360149362002019539920")
            ),
            Fq2::new(
                Fq::from("8839766648621210419302228913265679710586991805716981851373026244791934012854"),
                Fq::from("9103032146464138788288547957401673544458789595252696070370942789051858719203")
            ),
            Fq2::new(
                Fq::from("10378379548636866240502412547812481928323945124508039853766409196375806029865"),
                Fq::from("9021627154807648093720460686924074684389554332435186899318369174351765754041")
            )
        )
    );

    assert_eq!(finally, next);
}
