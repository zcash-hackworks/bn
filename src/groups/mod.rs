use std::ops::{Add,Sub,Neg,Mul};
use fields::{FieldElement, Fq, Fq2, Fq12, Fr, const_fq, fq2_nonresidue};
use arith::U256;
use std::fmt;
use rand::Rng;

use rustc_serialize::{Encodable, Encoder, Decodable, Decoder};

pub trait GroupElement: Sized +
                    Copy +
                    Clone +
                    PartialEq +
                    Eq +
                    fmt::Debug +
                    Add<Output=Self> +
                    Sub<Output=Self> +
                    Neg<Output=Self> +
                    Mul<Fr, Output=Self>
{
    fn zero() -> Self;
    fn one() -> Self;
    fn random<R: Rng>(rng: &mut R) -> Self;
    fn is_zero(&self) -> bool;
    fn double(&self) -> Self;
}

pub trait GroupParams: Sized {
    type Base: FieldElement + Decodable + Encodable;

    fn name() -> &'static str;
    fn one() -> G<Self>;
    fn coeff_b() -> Self::Base;
    fn check_order() -> bool { false }
}

#[repr(C)]
pub struct G<P: GroupParams> {
    x: P::Base,
    y: P::Base,
    z: P::Base
}

pub struct AffineG<P: GroupParams> {
    x: P::Base,
    y: P::Base
}

impl<P: GroupParams> PartialEq for AffineG<P> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<P: GroupParams> Eq for AffineG<P> { }

impl<P: GroupParams> fmt::Debug for G<P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}({:?}, {:?}, {:?})", P::name(), self.x, self.y, self.z)
    }
}

impl<P: GroupParams> Clone for G<P> {
    fn clone(&self) -> Self {
        G {
            x: self.x,
            y: self.y,
            z: self.z
        }
    }
}
impl<P: GroupParams> Copy for G<P> {}

impl<P: GroupParams> Clone for AffineG<P> {
    fn clone(&self) -> Self {
        AffineG {
            x: self.x,
            y: self.y
        }
    }
}
impl<P: GroupParams> Copy for AffineG<P> {}

impl<P: GroupParams> PartialEq for G<P> {
    fn eq(&self, other: &Self) -> bool {
        if self.is_zero() {
            return other.is_zero()
        }

        if other.is_zero() {
            return false;
        }

        let z1_squared = self.z.squared();
        let z2_squared = other.z.squared();

        if self.x * z2_squared != other.x * z1_squared {
            return false;
        }

        let z1_cubed = self.z * z1_squared;
        let z2_cubed = other.z * z2_squared;

        if self.y * z2_cubed != other.y * z1_cubed {
            return false;
        }

        return true;
    }
}
impl<P: GroupParams> Eq for G<P> { }

impl<P: GroupParams> G<P> {
    pub fn to_affine(&self) -> Option<AffineG<P>> {
        if self.z.is_zero() {
            None
        } else if self.z == P::Base::one() {
            Some(AffineG {
                x: self.x,
                y: self.y
            })
        } else {
            let zinv = self.z.inverse().unwrap();
            let zinv_squared = zinv.squared();

            Some(AffineG {
                x: self.x * zinv_squared,
                y: self.y * (zinv_squared * zinv)
            })
        }
    }
}

impl<P: GroupParams> AffineG<P> {
    pub fn to_jacobian(&self) -> G<P> {
        G {
            x: self.x,
            y: self.y,
            z: P::Base::one()
        }
    }
}

impl<P: GroupParams> Encodable for G<P> {
    fn encode<S: Encoder>(&self, s: &mut S) -> Result<(), S::Error> {
        if self.is_zero() {
            let l: u8 = 0;
            l.encode(s)
        } else {
            let l: u8 = 4;
            try!(l.encode(s));
            self.to_affine().unwrap().encode(s)
        }
    }
}

impl<P: GroupParams> Encodable for AffineG<P> {
    fn encode<S: Encoder>(&self, s: &mut S) -> Result<(), S::Error> {
        try!(self.x.encode(s));
        try!(self.y.encode(s));

        Ok(())
    }
}

impl<P: GroupParams> Decodable for G<P> {
    fn decode<S: Decoder>(s: &mut S) -> Result<G<P>, S::Error> {
        let l = try!(u8::decode(s));
        if l == 0 {
            Ok(G::zero())
        } else if l == 4 {
            Ok(try!(AffineG::decode(s)).to_jacobian())
        } else {
            Err(s.error("invalid leading byte for uncompressed group element"))
        }
    }
}

impl<P: GroupParams> Decodable for AffineG<P> {
    fn decode<S: Decoder>(s: &mut S) -> Result<AffineG<P>, S::Error> {
        let x = try!(P::Base::decode(s));
        let y = try!(P::Base::decode(s));

        // y^2 = x^3 + b
        if y.squared() == (x.squared() * x) + P::coeff_b() {
            if P::check_order() {
                let p: G<P> = G {
                    x: x,
                    y: y,
                    z: P::Base::one()
                };

                if (p * (-Fr::one())) + p != G::zero() {
                    return Err(s.error("point is not in the subgroup"))
                }
            }

            Ok(AffineG {
                x: x,
                y: y
            })
        } else {
            Err(s.error("point is not on the curve"))
        }
    }
}

impl<P: GroupParams> GroupElement for G<P> {
    fn zero() -> Self {
        G {
            x: P::Base::zero(),
            y: P::Base::one(),
            z: P::Base::zero()
        }
    }

    fn one() -> Self {
        P::one()
    }

    fn random<R: Rng>(rng: &mut R) -> Self {
        P::one() * Fr::random(rng)
    }

    fn is_zero(&self) -> bool {
        self.z.is_zero()
    }

    fn double(&self) -> Self {
        let a = self.x.squared();
        let b = self.y.squared();
        let c = b.squared();
        let mut d = (self.x + b).squared() - a - c;
        d = d + d;
        let e = a + a + a;
        let f = e.squared();
        let x3 = f - (d + d);
        let mut eight_c = c + c;
        eight_c = eight_c + eight_c;
        eight_c = eight_c + eight_c;
        let y1z1 = self.y * self.z;

        G {
            x: x3,
            y: e * (d - x3) - eight_c,
            z: y1z1 + y1z1
        }
    }
}

impl<P: GroupParams> Mul<Fr> for G<P> {
    type Output = G<P>;

    fn mul(self, other: Fr) -> G<P> {
        let mut res = G::zero();
        let mut found_one = false;

        for i in U256::from(other).bits() {
            if found_one {
                res = res.double();
            }

            if i {
                found_one = true;
                res = res + self;
            }
        }

        res
    }
}

impl<P: GroupParams> Add<G<P>> for G<P> {
    type Output = G<P>;

    fn add(self, other: G<P>) -> G<P> {
        if self.is_zero() {
            return other;
        }

        if other.is_zero() {
            return self;
        }

        let z1_squared = self.z.squared();
        let z2_squared = other.z.squared();
        let u1 = self.x * z2_squared;
        let u2 = other.x * z1_squared;
        let z1_cubed = self.z * z1_squared;
        let z2_cubed = other.z * z2_squared;
        let s1 = self.y * z2_cubed;
        let s2 = other.y * z1_cubed;

        if u1 == u2 && s1 == s2 {
            self.double()
        } else {
            let h = u2 - u1;
            let s2_minus_s1 = s2 - s1;
            let i = (h + h).squared();
            let j = h * i;
            let r = s2_minus_s1 + s2_minus_s1;
            let v = u1 * i;
            let s1_j = s1 * j;
            let x3 = r.squared() - j - (v + v);

            G {
                x: x3,
                y: r * (v - x3) - (s1_j + s1_j),
                z: ((self.z + other.z).squared() - z1_squared - z2_squared) * h
            }
        }
    }
}

impl<P: GroupParams> Neg for G<P> {
    type Output = G<P>;

    fn neg(self) -> G<P> {
        G {
            x: self.x,
            y: -self.y,
            z: self.z
        }
    }
}

impl<P: GroupParams> Neg for AffineG<P> {
    type Output = AffineG<P>;

    fn neg(self) -> AffineG<P> {
        AffineG {
            x: self.x,
            y: -self.y
        }
    }
}

impl<P: GroupParams> Sub<G<P>> for G<P> {
    type Output = G<P>;

    fn sub(self, other: G<P>) -> G<P> {
        self + (-other)
    }
}

pub struct G1Params;

impl GroupParams for G1Params {
    type Base = Fq;

    fn name() -> &'static str { "G1" }

    fn one() -> G<Self> {
        G {
            x: Fq::one(),
            y: const_fq([0xa6ba871b8b1e1b3a, 0x14f1d651eb8e167b, 0xccdd46def0f28c58, 0x1c14ef83340fbe5e]),
            z: Fq::one()
        }
    }

    fn coeff_b() -> Fq {
        const_fq([0x7a17caa950ad28d7, 0x1f6ac17ae15521b9, 0x334bea4e696bd284, 0x2a1f6744ce179d8e])
    }
}

pub type G1 = G<G1Params>;

pub struct G2Params;

impl GroupParams for G2Params {
    type Base = Fq2;

    fn name() -> &'static str { "G2" }

    fn one() -> G<Self> {
        G {
            x: Fq2::new(
                const_fq([0x8e83b5d102bc2026, 0xdceb1935497b0172, 0xfbb8264797811adf, 0x19573841af96503b]),
                const_fq([0xafb4737da84c6140, 0x6043dd5a5802d8c4, 0x09e950fc52a02f86, 0x14fef0833aea7b6b])
            ),
            y: Fq2::new(
                const_fq([0x619dfa9d886be9f6, 0xfe7fd297f59e9b78, 0xff9e1a62231b7dfe, 0x28fd7eebae9e4206]),
                const_fq([0x64095b56c71856ee, 0xdc57f922327d3cbb, 0x55f935be33351076, 0x0da4a0e693fd6482])
            ),
            z: Fq2::one()
        }
    }

    fn coeff_b() -> Fq2 {
        Fq2::new(
            const_fq([0x3bf938e377b802a8, 0x020b1b273633535d, 0x26b7edf049755260, 0x2514c6324384a86d]),
            const_fq([0x38e7ecccd1dcff67, 0x65f0b37d93ce0d3e, 0xd749d0dd22ac00aa, 0x0141b9ce4a688d4d])
        )
    }

    fn check_order() -> bool { true }
}

pub type G2 = G<G2Params>;

#[cfg(test)]
mod tests;

#[test]
fn test_g1() {
    tests::group_trials::<G1>();
}

#[test]
fn test_g2() {
    tests::group_trials::<G2>();
}

#[test]
fn test_affine_jacobian_conversion() {
    let rng = &mut ::rand::thread_rng();

    assert!(G1::zero().to_affine().is_none());
    assert!(G2::zero().to_affine().is_none());

    for _ in 0..1000 {
        let a = G1::one() * Fr::random(rng);
        let b = a.to_affine().unwrap();
        let c = b.to_jacobian();

        assert_eq!(a, c);
    }

    for _ in 0..1000 {
        let a = G2::one() * Fr::random(rng);
        let b = a.to_affine().unwrap();
        let c = b.to_jacobian();

        assert_eq!(a, c);
    }
}

#[inline]
fn twist() -> Fq2 {
    fq2_nonresidue()
}

#[inline]
fn two_inv() -> Fq {
    const_fq([9781510331150239090, 15059239858463337189, 10331104244869713732, 2249375503248834476])
}

#[inline]
fn ate_loop_count() -> U256 {
    U256([0x9d797039be763ba8, 0x0000000000000001, 0x0000000000000000, 0x0000000000000000])
}

#[inline]
fn twist_mul_by_q_x() -> Fq2 {
    Fq2::new(
        const_fq([13075984984163199792, 3782902503040509012, 8791150885551868305, 1825854335138010348]),
        const_fq([7963664994991228759, 12257807996192067905, 13179524609921305146, 2767831111890561987])
    )
}

#[inline]
fn twist_mul_by_q_y() -> Fq2 {
    Fq2::new(
        const_fq([16482010305593259561, 13488546290961988299, 3578621962720924518, 2681173117283399901]),
        const_fq([11661927080404088775, 553939530661941723, 7860678177968807019, 3208568454732775116])
    )
}

#[derive(PartialEq, Eq)]
pub struct EllCoeffs {
    pub ell_0: Fq2,
    pub ell_vw: Fq2,
    pub ell_vv: Fq2
}

#[derive(PartialEq, Eq)]
pub struct G2Precomp {
    pub q: AffineG<G2Params>,
    pub coeffs: Vec<EllCoeffs>
}

impl G2Precomp {
    pub fn miller_loop(&self, g1: &AffineG<G1Params>) -> Fq12 {
        let mut f = Fq12::one();

        let mut idx = 0;

        let mut found_one = false;

        for i in ate_loop_count().bits() {
            if !found_one {
                // skips the first bit
                found_one = i;
                continue;
            }

            let c = &self.coeffs[idx];
            idx += 1;
            f = f.squared().mul_by_024(c.ell_0, c.ell_vw.scale(g1.y), c.ell_vv.scale(g1.x));

            if i {
                let c = &self.coeffs[idx];
                idx += 1;
                f = f.mul_by_024(c.ell_0, c.ell_vw.scale(g1.y), c.ell_vv.scale(g1.x));
            }
        }

        let c = &self.coeffs[idx];
        idx += 1;
        f = f.mul_by_024(c.ell_0, c.ell_vw.scale(g1.y), c.ell_vv.scale(g1.x));

        let c = &self.coeffs[idx];
        f = f.mul_by_024(c.ell_0, c.ell_vw.scale(g1.y), c.ell_vv.scale(g1.x));

        f
    }
}

#[test]
fn test_miller_loop() {
    use fields::Fq6;

    let g1 = G1::one() * Fr::from_str("18097487326282793650237947474982649264364522469319914492172746413872781676").unwrap();
    let g2 = G2::one() * Fr::from_str("20390255904278144451778773028944684152769293537511418234311120800877067946").unwrap();

    let g1_pre = g1.to_affine().unwrap();
    let g2_pre = g2.to_affine().unwrap().precompute();

    let gt = g2_pre.miller_loop(&g1_pre);

    assert_eq!(gt,
        Fq12::new(Fq6::new(
                Fq2::new(Fq::from_str("14551901853310307118181117653102171756020286507151693083446930124375536995872").unwrap(), Fq::from_str("9312135802322424742640599513015426415694425842442244572104764725304978020017").unwrap()), 
                Fq2::new(Fq::from_str("2008578374540014049115224515107136454624926345291695498760935593377832328658").unwrap(), Fq::from_str("19401931167387470703307774451905975977586101231060812348184567722817888018105").unwrap()), 
                Fq2::new(Fq::from_str("15835061253582829097893482726334173316772697321004871665993836763948321578465").unwrap(), Fq::from_str("2434436628082562384254182545550914004674636606111293955202388712261962820365").unwrap())
            ),
            Fq6::new(
                Fq2::new(Fq::from_str("2874440054453559166574356420729655370224872280550180463983603224123901706537").unwrap(), Fq::from_str("21199736323249863378180814900160978651989782296293186487853700340281870105680").unwrap()), 
                Fq2::new(Fq::from_str("19165582755854282767090326095669835261356341739532443976394958023142879015770").unwrap(), Fq::from_str("1381947898997178910398427566832118260186305708991760706544743699683050330259").unwrap()), 
                Fq2::new(Fq::from_str("282285618133171001983721596014922591835675934808772882476123488581876545578").unwrap(), Fq::from_str("9533292755262567365755835323107174518472361243562718718917822947506880920117").unwrap())
            )
        )
    );
}

impl AffineG<G2Params> {
    fn mul_by_q(&self) -> Self {
        AffineG {
            x: twist_mul_by_q_x() * self.x.frobenius_map(1),
            y: twist_mul_by_q_y() * self.y.frobenius_map(1)
        }
    }

    pub fn precompute(&self) -> G2Precomp {
        let mut r = self.to_jacobian();

        let mut coeffs = Vec::with_capacity(102);

        let mut found_one = false;

        for i in ate_loop_count().bits() {
            if !found_one {
                // skips the first bit
                found_one = i;
                continue;
            }

            coeffs.push(r.doubling_step_for_flipped_miller_loop());

            if i {
                coeffs.push(r.mixed_addition_step_for_flipped_miller_loop(self));
            }
        }

        let q1 = self.mul_by_q();
        let q2 = -(q1.mul_by_q());

        coeffs.push(r.mixed_addition_step_for_flipped_miller_loop(&q1));
        coeffs.push(r.mixed_addition_step_for_flipped_miller_loop(&q2));

        G2Precomp {
            q: *self,
            coeffs: coeffs
        }
    }
}

impl G2 {
    fn mixed_addition_step_for_flipped_miller_loop(&mut self, base: &AffineG<G2Params>) -> EllCoeffs {
        let d = self.x - self.z * base.x;
        let e = self.y - self.z * base.y;
        let f = d.squared();
        let g = e.squared();
        let h = d * f;
        let i = self.x * f;
        let j = self.z * g + h - (i + i);

        self.x = d * j;
        self.y = e * (i - j) - h * self.y;
        self.z = self.z * h;

        EllCoeffs {
            ell_0: twist() * (e * base.x - d * base.y),
            ell_vv: e.neg(),
            ell_vw: d
        }
    }

    fn doubling_step_for_flipped_miller_loop(&mut self) -> EllCoeffs {
        let a = (self.x * self.y).scale(two_inv());
        let b = self.y.squared();
        let c = self.z.squared();
        let d = c + c + c;
        let e = G2Params::coeff_b() * d;
        let f = e + e + e;
        let g = (b + f).scale(two_inv());
        let h = (self.y + self.z).squared() - (b + c);
        let i = e - b;
        let j = self.x.squared();
        let e_sq = e.squared();

        self.x = a * (b - f);
        self.y = g.squared() - (e_sq + e_sq + e_sq);
        self.z = b * h;

        EllCoeffs {
            ell_0: twist() * i,
            ell_vw: h.neg(),
            ell_vv: j + j + j
        }
    }
}

#[test]
fn test_prepared_g2() {
    let g2 = G2::one() * Fr::from_str("20390255904278144451778773028944684152769293537511418234311120800877067946").unwrap();

    let g2_p = g2.to_affine().unwrap().precompute();

    let expected_g2_p = G2Precomp {
        q: AffineG {
            x: Fq2::new(
                Fq::from_str("13936578204263895229092967414825041569724079982537616431212348844388899776640").unwrap(),
                Fq::from_str("6372636725635053773371996212293600406870925440022386078671828127711809436031").unwrap()
            ),
            y: Fq2::new(
                Fq::from_str("19293035970010898452454381709939058714495082898872914526540420247178075881697").unwrap(),
                Fq::from_str("13822349107533275437553410197128434338071760794202849712402800746887549494974").unwrap(),
            )
        },
        coeffs: vec![
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("2627043964130257285960335798481049684473505235282434252362906013168360965985").unwrap(), Fq::from_str("14787221188041042838526170263203159226721816513593198738527261794343476005673").unwrap()), ell_vw: Fq2::new(Fq::from_str("5190413803656753539584048070636432748402456516849818272297235294934300653772").unwrap(), Fq::from_str("16131787528611999569385991096257681501249100726189947900572474295515353427218").unwrap()), ell_vv: Fq2::new(Fq::from_str("11284217811624345285836951206193951701052480344824103057273798441033858118424").unwrap(), Fq::from_str("6392116536365389562188363539617434898082427085812013656312597845585909267447").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("827617134098165717451808940080463277390770457691666780560712143809003953598").unwrap(), Fq::from_str("6776229088211374446530917353066321938640163548858035832637439231634790575465").unwrap()), ell_vw: Fq2::new(Fq::from_str("987776078024262725561041258416387561158070255475504730561661362421251696401").unwrap(), Fq::from_str("15312963471998242334683179861466148222641884112991952428739813077336923189144").unwrap()), ell_vv: Fq2::new(Fq::from_str("2813988028633040066320201189843971639620433430176492766961373503539074898364").unwrap(), Fq::from_str("17055167212030988864288747645634552775658830115224987200613220012554982651578").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("17093477194591041266397380404112224367919571835678040889250664706433487493043").unwrap(), Fq::from_str("8643160753646550179100165428401120938543691942215161396502118430923535820500").unwrap()), ell_vw: Fq2::new(Fq::from_str("10933962898922943120024964690690003920835888964592980369975615298196656183772").unwrap(), Fq::from_str("10054435552662455933652209211749007190981947542684285471884164039871580864235").unwrap()), ell_vv: Fq2::new(Fq::from_str("1475406195060586931258578851599932495282912871465827155842007043242102046659").unwrap(), Fq::from_str("5049830924161187876168845800328902567850325970867534538120475909389784841990").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("10434276065336457961840942984114655645641355797054067910262735814196303889171").unwrap(), Fq::from_str("13878009844584745291193244072670264656592436467816383615606044714978763054890").unwrap()), ell_vw: Fq2::new(Fq::from_str("20993505060568459085534388708128029516059470166174503494622391683433030357130").unwrap(), Fq::from_str("10192769806017258272908841309051731389378033999625803563624881649401237928876").unwrap()), ell_vv: Fq2::new(Fq::from_str("17377091829483421118284926147077357677507183290871250621761583281560333718550").unwrap(), Fq::from_str("18297780639039540893260206902113850213293978821982237672423968999898103045256").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("5064242330961655837810472366173802140673816516016739528934082190598759216129").unwrap(), Fq::from_str("20172028727708469864987113767379594545870343866584998839699240724281670882460").unwrap()), ell_vw: Fq2::new(Fq::from_str("2996734698689069564052128450078797792056313730363946953355442668966438200951").unwrap(), Fq::from_str("9910941594900797404370917094355311738210439203708671856236954583027355115302").unwrap()), ell_vv: Fq2::new(Fq::from_str("6899792008443933570863502403567742315007264221354931185809349393350618040985").unwrap(), Fq::from_str("18692424201235977456836009406318755884773471215777102117879530727096807383696").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("9406196439119227017197200697795682187217614676635965071601916278141238949664").unwrap(), Fq::from_str("1156397052245190909533108891675012415349316096091404012699655652677336827483").unwrap()), ell_vw: Fq2::new(Fq::from_str("9249275523859165912800910211783530541095037634552893633043938333737241478198").unwrap(), Fq::from_str("13674086724537885208439774394455021050268654286498757861381516198880504360984").unwrap()), ell_vv: Fq2::new(Fq::from_str("5760486636714317173624828119078349011897522431614835257025421978626292083787").unwrap(), Fq::from_str("19195594701877562770296730558597966521343920474874926443609709880727520577918").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("12166638768862632131967787101186130559865354321458713212316891357299196109517").unwrap(), Fq::from_str("1517854601414954363921968422436855204960697222979994177630673862982254831193").unwrap()), ell_vw: Fq2::new(Fq::from_str("13366104107979832799705928625277388583673352228059815944738106504080223007902").unwrap(), Fq::from_str("21648874073566751720910466086303831459661520354253552089758049591613143727506").unwrap()), ell_vv: Fq2::new(Fq::from_str("3488776451414334656597670851583995239177248056342375163461521452312663358741").unwrap(), Fq::from_str("20333414955874213710027365371143864644434939642550352414403887754960702108637").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("20468104499950674906991599591024893438007901046236752159815194628487603192867").unwrap(), Fq::from_str("19069685413250623899891139234721987960230148622397590130148098858949057649761").unwrap()), ell_vw: Fq2::new(Fq::from_str("2458162095258922865086779772179745380154704066745856173311580123496018896409").unwrap(), Fq::from_str("14594824564551859025266692529845827840804005950672166057231487026855067744741").unwrap()), ell_vv: Fq2::new(Fq::from_str("9689655654286496198845571066967567466951511863880400774464590553022163229348").unwrap(), Fq::from_str("10870127206120874488753221889350215535903993427218589071224774335669836789433").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("15334990572515495832170458963151379997993357347608954164027751164770129172212").unwrap(), Fq::from_str("13576640273125978757652414940234881692674718901710647917847939367156260249620").unwrap()), ell_vw: Fq2::new(Fq::from_str("17441399177222030197667035179458353890777336233312342937564354849059149975908").unwrap(), Fq::from_str("17007488955907409465543686050546106251553784366237996525801373957978511005278").unwrap()), ell_vv: Fq2::new(Fq::from_str("11907792503811346855438657769533223907267837220830333730410940125854558637383").unwrap(), Fq::from_str("19490060942462197243379937108210110931787757078893580540054690416131017931738").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("17315779415820661602591970083599829709759601467527578787441772726289774228451").unwrap(), Fq::from_str("12080802446892182171725856453369856169039939341293608060423916935382017714476").unwrap()), ell_vw: Fq2::new(Fq::from_str("12912774602739952889552025090054804346758888115608888766725054474682313950190").unwrap(), Fq::from_str("18204955971588024842905918028829132053489100674596042021111813303048407250440").unwrap()), ell_vv: Fq2::new(Fq::from_str("15721987682182104464109212111163313043868345920774975073102123159159010623181").unwrap(), Fq::from_str("12041652399509652544977783565303273845883245418094457529808044120827920548972").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("19999070717670362479390656255660288732744476951677182489244240784338456791622").unwrap(), Fq::from_str("3906002755238658676658974095984902840159437088441472633765497396881255783772").unwrap()), ell_vw: Fq2::new(Fq::from_str("4998301920940599524538129746790252820575786197710804708434384224972978627969").unwrap(), Fq::from_str("13740819604543965458992763357705315276969968429480411530007322182403523242226").unwrap()), ell_vv: Fq2::new(Fq::from_str("12357083524385039711813704948190809675818592911528835736620903716727226601953").unwrap(), Fq::from_str("13949012942996636079549854439740332990496218971788334041233648192605019948979").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("19805154539184543843735550309592886829271314363908537964547876413114433734624").unwrap(), Fq::from_str("16481416862400234230001935040590314552205484054708993091993487440981092053599").unwrap()), ell_vw: Fq2::new(Fq::from_str("2850933690997557614680452825756662353227291804228932023648526608593113870536").unwrap(), Fq::from_str("4786423507519240726599235257249783193682886531662683554713943115802063984530").unwrap()), ell_vv: Fq2::new(Fq::from_str("1023533030117941985730558522296044345509721849960856436490999820443312747347").unwrap(), Fq::from_str("16431548648882316869128893343209930754508425074299553782883149795599907406302").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("17203075478724975791441470425150871365957382183648220405350599958349478966969").unwrap(), Fq::from_str("19566853320706390536799238411444563555185156853468371794865496424286288424312").unwrap()), ell_vw: Fq2::new(Fq::from_str("15310908404682220539815850064952575682898427316407655213706791187971849703724").unwrap(), Fq::from_str("16518878079785062271298378680376701241016697782308466806532677867880353948389").unwrap()), ell_vv: Fq2::new(Fq::from_str("14461573034674366293217641735721285869692106578982161360427147953320122006893").unwrap(), Fq::from_str("6682535969691272372531048932989583075363217464304983455197465867679730319705").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("20895553497690200905762157235991552162755238828234592755315911339400881373107").unwrap(), Fq::from_str("11947919061415103173593352950044808416347417810515010536678086292686107222044").unwrap()), ell_vw: Fq2::new(Fq::from_str("6277981470671734560229149673222400771590408919719217743512900463425397921908").unwrap(), Fq::from_str("857212248193599195410711267373529614085739023676782457331058533755337294689").unwrap()), ell_vv: Fq2::new(Fq::from_str("13885916052404740835214284721162700704805447448922857617672445286139626987191").unwrap(), Fq::from_str("19585389042100086858633036225209396447924441966778662250950999864988284098182").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("12748401205319976481200663793287165427761341490879316890008250070038707706903").unwrap(), Fq::from_str("9347025175062646190962067344106887221752685406563421042329129538495481178145").unwrap()), ell_vw: Fq2::new(Fq::from_str("12819473645128401077507234830331609583133696659610406468927174080626676619197").unwrap(), Fq::from_str("9039389340404634934394726812078086316924048399159604642340339967187032834011").unwrap()), ell_vv: Fq2::new(Fq::from_str("3576972655652865483567741810151575278761470911419381093571973994483543095285").unwrap(), Fq::from_str("15017468236028111406921228193688973260631093975412858847563400987354453569177").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("15407570760857043664675198351546637734249070940748070513828892036078049928562").unwrap(), Fq::from_str("16864837887933576789187120028615006826049916557295792205508706545038103706967").unwrap()), ell_vw: Fq2::new(Fq::from_str("10721360465919602144653520112889689786555051462198759164582568961994647563952").unwrap(), Fq::from_str("15294101887395462362997873520785213441488786364828649600395509678143054744021").unwrap()), ell_vv: Fq2::new(Fq::from_str("338877711237691732068166553384070140091720428235937017002925264628733553275").unwrap(), Fq::from_str("8565770517833678396194963042166078986348581926120249815567974171840711179000").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("17302401852599488824262121640328139426135883695023154134865582322858078206108").unwrap(), Fq::from_str("17762524466294762819341377108361075499043287383142581974755577861093917129883").unwrap()), ell_vw: Fq2::new(Fq::from_str("13015933066581331586057031655765394156438563571733887223640331819027233831808").unwrap(), Fq::from_str("758435331870656231667722761984400968773108978847112431665229252181053186864").unwrap()), ell_vv: Fq2::new(Fq::from_str("1969703955244816291771524536928710541897907643982046151796341732075930805350").unwrap(), Fq::from_str("7031491051442716548063495812915693048061021502544985435130728834982309403017").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("3256431897334799395000264533777543923057110731812322608781781329228902421476").unwrap(), Fq::from_str("12181683857561584905762347955319259913388176042023651776593530168946356566719").unwrap()), ell_vw: Fq2::new(Fq::from_str("796109596241493969908761241268545233387446802589790513588935996990189190052").unwrap(), Fq::from_str("1881011435205668659920588004752042238985141500261993259857887774967416952583").unwrap()), ell_vv: Fq2::new(Fq::from_str("8939978626202472531630769965940628317376245215184468794030552730424109408120").unwrap(), Fq::from_str("6876296539554825957959953448623382464100292373640362636397165101276891176372").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("9977776755404508431321449793657366203354342168643747704006143818738023991417").unwrap(), Fq::from_str("17287266532802977377832339048078093713694407321084566908200498605546793435820").unwrap()), ell_vw: Fq2::new(Fq::from_str("13828304001931491210449888619926027655755205497065268858009259221646932328015").unwrap(), Fq::from_str("15160507616608365341639459841002920786281000303114113782821091199951204875277").unwrap()), ell_vv: Fq2::new(Fq::from_str("11161349646944933007552237093278324942110971730799870343491465916610232191523").unwrap(), Fq::from_str("1406353163795290818414111649956135144227974028417320470897284179379367096127").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("7051926909684750402992796647032741235565177570654429988875094870973533133596").unwrap(), Fq::from_str("20900584116079647608136401104271956963096697793208792824067775068966626057321").unwrap()), ell_vw: Fq2::new(Fq::from_str("16845876804097000415961425779514720605919619149076212452284059758659650969766").unwrap(), Fq::from_str("20806461795944565303534619813753388058736049320841472653130287521152836019826").unwrap()), ell_vv: Fq2::new(Fq::from_str("14700937697004489523945953560825193972493636932341330040072362595384517591325").unwrap(), Fq::from_str("17137955424118909932502168320776219042080643724609838664247554484115994974522").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("11997142079850902731018491830701934107296522353910410046841604519639199995406").unwrap(), Fq::from_str("7846127392296015970710836702134527490160588060582987110380879354952036676012").unwrap()), ell_vw: Fq2::new(Fq::from_str("17663414534578185845724108244425572237371794484114570058559740097895281553299").unwrap(), Fq::from_str("5029140669755495687983399191536061323171057684732073778423303267496666431971").unwrap()), ell_vv: Fq2::new(Fq::from_str("18428442704115060890708987869503078260482829714024061977805185634781850796860").unwrap(), Fq::from_str("10997089106320530533400287751108615584976629751702164406696900094998369084403").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("13049677819639682137107272322562642308503411465986079606334189062279441432203").unwrap(), Fq::from_str("8874624972677935217588093490885486721673594121072678233714072789659128555367").unwrap()), ell_vw: Fq2::new(Fq::from_str("17057994029556488440657870420250293132874472170733557198308184169150330844407").unwrap(), Fq::from_str("21633949388482190843211557537195287403106536528012825991981214057356654097223").unwrap()), ell_vv: Fq2::new(Fq::from_str("6897044784874980042628318562068623805084522888972873438474139199226581764756").unwrap(), Fq::from_str("2692625602037608101594849144393949134766066123857076012132350588113519149572").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("9195356099870032305655798667367006686897880523242566329014149856866611060003").unwrap(), Fq::from_str("8052880750949660720559492348565167987253311417692073778040482170195175946428").unwrap()), ell_vw: Fq2::new(Fq::from_str("14126308928695836758114496610102001680995565158909643521688860213965125059872").unwrap(), Fq::from_str("671075569038912473223544948981506373517307914345286644404839238190008163502").unwrap()), ell_vv: Fq2::new(Fq::from_str("16828103633974871724641544770525140436938333764966095252502003672344759889143").unwrap(), Fq::from_str("11294529501746593948507624028584472676389760897716457892342249916318954501723").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("21124500310722447577044753724306337899524073055454053171030128408285122310058").unwrap(), Fq::from_str("6266274689105996173186441194758121417834460925587732781640173398761521329055").unwrap()), ell_vw: Fq2::new(Fq::from_str("21576331134605070996263178912402854747994495442070419140713471256922277734125").unwrap(), Fq::from_str("10074747339717657657921727946300324710951646905879034242434327125816600264379").unwrap()), ell_vv: Fq2::new(Fq::from_str("21765519811276677250031796638148457542841530645055777846258335421409598000534").unwrap(), Fq::from_str("21113224637771003677145667454854210585751245642107143223999187554888833454243").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("3972300323627467772646468544189830113523911494046612234673787196935654523692").unwrap(), Fq::from_str("4678843490808368914310764558657240771944434980659329233409943752329617161210").unwrap()), ell_vw: Fq2::new(Fq::from_str("17273599953138772488805612764323870175561796764345770950040891840763740358904").unwrap(), Fq::from_str("20415406312940209525918367703821373463952436739154646433333601347381723237136").unwrap()), ell_vv: Fq2::new(Fq::from_str("9160987523541118769404859219229455872777126347530590315596894278463599930633").unwrap(), Fq::from_str("9054614112960826755739076858442558227570393423663838494067070216861352843202").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("9994954547022379603539679290217274617475455164007187108524364097751858323148").unwrap(), Fq::from_str("20207771379440433457553782484784849880846158790602128386495229840158065235194").unwrap()), ell_vw: Fq2::new(Fq::from_str("4256518944443059194982602310083745017699415852751800911480635583088649250060").unwrap(), Fq::from_str("15878038811235111024146051545935959211923223574364333740058604258589674044428").unwrap()), ell_vv: Fq2::new(Fq::from_str("15953919860849853924955813198899405295748673040941543212644580663271517228860").unwrap(), Fq::from_str("19594437991878043589671466710739298733201896072505590020656735731744884281871").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("1333543700232406379832215757617068657339748039338978081735192495239856176346").unwrap(), Fq::from_str("13583045269132775344737980756385257860562307101313068058661941588822340692712").unwrap()), ell_vw: Fq2::new(Fq::from_str("11316045827391404343818336116749364562466461792138476673985290066961360847098").unwrap(), Fq::from_str("6683447190424452097068562814518005185732181765401062334661636064120276671717").unwrap()), ell_vv: Fq2::new(Fq::from_str("21483741541431414142042872578822487244163622589552506298253453995156009765340").unwrap(), Fq::from_str("16502358172728062550923345759488126380665850742066670591562640787191536290296").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("4727679401549029925466046980042210063407620223804057104870816459356196798234").unwrap(), Fq::from_str("4245208253099058296636688301996648265313738760450861302514451065904449375848").unwrap()), ell_vw: Fq2::new(Fq::from_str("19589426516750929424306774984752060788246605243046082989256140076730177865962").unwrap(), Fq::from_str("19202599366770569221216501935740579180019961851730854119509268985533547631061").unwrap()), ell_vv: Fq2::new(Fq::from_str("18329650333197074546830734977650617140322812110454751288122555715911145072440").unwrap(), Fq::from_str("18685828888894718887886440900247948352873037333999376257561751598045361934940").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("15348453159719894937627693446890053385899104841918429189733580400902437025934").unwrap(), Fq::from_str("7501190780205720975723852235835733848687299193764526538831345783610137314719").unwrap()), ell_vw: Fq2::new(Fq::from_str("12341649147852680664345690534373468938440543264264387581538764313149640739302").unwrap(), Fq::from_str("140543540021318882234026314588789476198722270107997413913717074968501481488").unwrap()), ell_vv: Fq2::new(Fq::from_str("5960374673747504135132782820167542238733653290730217772559545641742563125282").unwrap(), Fq::from_str("14507139191276441762644240084427079430815677404266073216544260783868466836469").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("1463974716527100203803811366059855506152810619047460592789012516365557382809").unwrap(), Fq::from_str("12894560712946121716686238503640546420854178898471895472462918225441749701869").unwrap()), ell_vw: Fq2::new(Fq::from_str("10508120321944855490271402587060090544503473558782513226881398929302598320386").unwrap(), Fq::from_str("20659878091446285231719184399846880718538938504006216923931157918871132219202").unwrap()), ell_vv: Fq2::new(Fq::from_str("9314687360620406591930278900314522432964075856747694070384017245464130048974").unwrap(), Fq::from_str("17694675150016443475335083504221747995503842086249559620223686305992007745569").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("14196539467892272085884829545587228956041818801561947743113049357924737022135").unwrap(), Fq::from_str("11319357657686312684765696405805456576475636756706867678717962076341276658626").unwrap()), ell_vw: Fq2::new(Fq::from_str("2758949814978246411081695103954624266457693267107844100147125044169545866221").unwrap(), Fq::from_str("9492100316900274072027501035053797502489796059345864111677781413985743864383").unwrap()), ell_vv: Fq2::new(Fq::from_str("4704107896591496617393915690357581449785399000813194972371050438639491862314").unwrap(), Fq::from_str("11531365345870313813093010071065540963321572589131340821665512078366854738990").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("18802465664727774773666122120775907428094487475952383150269507875840089541246").unwrap(), Fq::from_str("1731637922302665944958348222164323910730404298281499149800917211427254836819").unwrap()), ell_vw: Fq2::new(Fq::from_str("9297609456558817283944780415839637083190742770917810313181849518385379162642").unwrap(), Fq::from_str("16377246183633200781091121197835476043913744319877154056563766308564474143489").unwrap()), ell_vv: Fq2::new(Fq::from_str("9542461764772800142773586971354345045722377861760780709590923740133551871717").unwrap(), Fq::from_str("21640615613954317564353959148225302612886621120325752435425437529806262143694").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("19906049876176293236988182778475214525096978055900789654464603831282713330781").unwrap(), Fq::from_str("20362580925091945825523087376754169675844247448742703441518744743887060798158").unwrap()), ell_vw: Fq2::new(Fq::from_str("12267599301269062964822189844562856172284191535014283588478875529832316401461").unwrap(), Fq::from_str("16669825137620567269982790548542800570592942517936126877686519256685036586668").unwrap()), ell_vv: Fq2::new(Fq::from_str("8883573260261802165780406246211290131940260376316892224137508921505842308805").unwrap(), Fq::from_str("19231108274886583687672953546761020040720197580241947146286750565108059707992").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("634182853515641884871095612114627494589971356340786857463490962898713112220").unwrap(), Fq::from_str("11430884471074068345749191262196070058384500365559620896062382995797164195815").unwrap()), ell_vw: Fq2::new(Fq::from_str("10550902705393315932525976596381459439846414130661652807183357995272197708172").unwrap(), Fq::from_str("3961587743443277010160905820636488892023891053433292002931434257476626412826").unwrap()), ell_vv: Fq2::new(Fq::from_str("12451440882697904111698352784988882971851697245540249016187908252126835129000").unwrap(), Fq::from_str("11782592623486603129029758483031965320434451423371375539247662743178943477576").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("1404156345033137762434737426689373534725204733157905834242027185649572912742").unwrap(), Fq::from_str("19570890322851668915358612375948367153948206621537633922343603465980161569078").unwrap()), ell_vw: Fq2::new(Fq::from_str("8686503710040621382075118769097735385052600877432244428197116763766844918162").unwrap(), Fq::from_str("17223447941612107204326938316335236704702567123017447936901253775530936851086").unwrap()), ell_vv: Fq2::new(Fq::from_str("18767158861904777443587157807982315566545548385817262926280782765269843879565").unwrap(), Fq::from_str("21334543677892438295079840326280352051766440588072368601899791201007143326750").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("20453014818098756686381445228623803935714316636974320420742351477862354368007").unwrap(), Fq::from_str("11004429848286161246832723709453157691786114066215975460914111961981461575367").unwrap()), ell_vw: Fq2::new(Fq::from_str("14444837779919119651023367195508370508399423510380962187645910198528103345317").unwrap(), Fq::from_str("4021193867366929691652413553892868612329317391162172266658368471326882495507").unwrap()), ell_vv: Fq2::new(Fq::from_str("21098047915839588616030399663527845846589272989005835932893550785408588982098").unwrap(), Fq::from_str("16796487305049888740077005894512139332211008298712590271618863352993144873711").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("8799932768264091269394526725962837485495458036153083882004277055769655900574").unwrap(), Fq::from_str("18596145973176446599121174893654309094121045316907408604913712044917019537350").unwrap()), ell_vw: Fq2::new(Fq::from_str("5727901458562943136998396728651681550262281371313580881804977990767435431080").unwrap(), Fq::from_str("16665276274186006945492338822028666786935953718794303569570538737692828117635").unwrap()), ell_vv: Fq2::new(Fq::from_str("10910280368195212990489805322100030141483989387720003214110076576431959738476").unwrap(), Fq::from_str("9606982070265587786427675480691869167035030735018706315315080889814741903522").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("14453637744147393731588213406970287317078500183808667567110703622821377404467").unwrap(), Fq::from_str("2399421543794348043991999124897133983828094348779941685913680103324762089779").unwrap()), ell_vw: Fq2::new(Fq::from_str("9949473626857777965864312230577133154733097780222892008806269288857107115742").unwrap(), Fq::from_str("3174442934451895428948522429276679446276043250911620596846605340399481680631").unwrap()), ell_vv: Fq2::new(Fq::from_str("10773554875563236900680111254025179423521083612506388720894431016952344645908").unwrap(), Fq::from_str("1791500946161395076321848651097476041018873800532919907790655022966446958329").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("13228018442560564090007964181563674306630092611719086476342149061864561185280").unwrap(), Fq::from_str("12610285327624889222715450538150491139050620893563766277510163987972586901233").unwrap()), ell_vw: Fq2::new(Fq::from_str("17479662858773093405375538474745376030990506794065246603293683266124925174043").unwrap(), Fq::from_str("2084404960859519365733972471606438409123459865402210971655148862692485946826").unwrap()), ell_vv: Fq2::new(Fq::from_str("9675635536108496330861750800277390856149950003475979013330722219857375990089").unwrap(), Fq::from_str("17889926389980731198337394126857435112726851886206769936779651682292365620935").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("19205463777127243508919644186723680978776945615422050305472500178064412123473").unwrap(), Fq::from_str("10530028762478985738084339660443655112791893212694172883230610760959423970096").unwrap()), ell_vw: Fq2::new(Fq::from_str("16515559472874148379011472911124765449097459867226728232468301906190158644000").unwrap(), Fq::from_str("17379548905015313648731271429222704537289418711815257252434537403381161039222").unwrap()), ell_vv: Fq2::new(Fq::from_str("13018298656134311744614296809909600459918460239498208291306800551968894517759").unwrap(), Fq::from_str("16285303269049929546880902543736204779249639088136468614535789316688521663976").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("9506592875949693201584225006409064359870098448027454395748349563972724759341").unwrap(), Fq::from_str("1756212124305090096775291444252901778765176938185769777272403090820692916832").unwrap()), ell_vw: Fq2::new(Fq::from_str("8608069305722185875126845710591577689851258364137915830961192516727536204958").unwrap(), Fq::from_str("7879742654857051185278955345176850926952262018841407228150061698899939295057").unwrap()), ell_vv: Fq2::new(Fq::from_str("2311752966405194522531977544398652994349075269275918521122768201037004674541").unwrap(), Fq::from_str("8056192790171072186146581305753657840809572079825584316343698126369167589177").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("16894171323623324294071017663897158494527746239981369238071300334653959886935").unwrap(), Fq::from_str("21839557351524624846710837291698003965911117746801809883425361529677608757736").unwrap()), ell_vw: Fq2::new(Fq::from_str("7921995242367245298209685641602469741817826911044923939882608549077724733892").unwrap(), Fq::from_str("16195366803093832645900011617380218833645406154510352140817411408898365560241").unwrap()), ell_vv: Fq2::new(Fq::from_str("12718951609280728995721445896408716663776093135974310357265585518411722037823").unwrap(), Fq::from_str("14775730569800772809191320148533632897226925930462946258120395650683306784222").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("18324740767049924394294219212760164185884801967765967078427354810675386933254").unwrap(), Fq::from_str("18510525798971302215942990229348890627301667188256739950256532108744423958467").unwrap()), ell_vw: Fq2::new(Fq::from_str("392130232730851920199052091589353916248087256056205523658797958171832234380").unwrap(), Fq::from_str("21095230708713942590602820178155041621291598326436725983218983799781802054278").unwrap()), ell_vv: Fq2::new(Fq::from_str("3496296009938486548567151348321691013571513425434397040346440937975425900960").unwrap(), Fq::from_str("5982055369526470853931810147684497423036063399292908540001447206525852512531").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("19247021831371751344058078378442480867649167715937305884736113857856820190448").unwrap(), Fq::from_str("3844115161163747856814522390569661892955450798273343821401275758096863392646").unwrap()), ell_vw: Fq2::new(Fq::from_str("6827864046826742026186179073427918229064686953470940699125130948982461542396").unwrap(), Fq::from_str("17615043663229174389408922909891358248158892597678816385681324344031892930590").unwrap()), ell_vv: Fq2::new(Fq::from_str("17143020152984715373298117189870583206488981021188881359521310100992509008220").unwrap(), Fq::from_str("2174519777022093979454116656397980458948553322768446892932246705483799164016").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("7787149350876246869417695770873628344850130125888658086308716316049428077539").unwrap(), Fq::from_str("3645376584593179034263104864050174094110419584940335745722492256564582048708").unwrap()), ell_vw: Fq2::new(Fq::from_str("11834206166934740935662681278907472701573638397628522018255718023777042308993").unwrap(), Fq::from_str("11564024956967266505944728035150285591786947709753433733529805379472238406883").unwrap()), ell_vv: Fq2::new(Fq::from_str("15009525058621581543328262019401183048318643177731888421709672368574056014500").unwrap(), Fq::from_str("13573139489477783981960389988307434036015536834494019264597795654229345703051").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("16059660214342726778635842789601639915164060169646706775268999439826160969051").unwrap(), Fq::from_str("3592864998692836817549974752747753891215991801035643051485066473571282986623").unwrap()), ell_vw: Fq2::new(Fq::from_str("16014086539360826224382965134738311010906671679277894939278016943520015391650").unwrap(), Fq::from_str("8080101227643664909861133876167034094240488311720250596061359317310944652491").unwrap()), ell_vv: Fq2::new(Fq::from_str("3252331663236937425057834000082561720577180169369353108561319513563871033743").unwrap(), Fq::from_str("11096433584385810520286571872297862049422429286608206525566845219548956753724").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("19106172305007554707393707438702629275907204152148571359255256006833275609105").unwrap(), Fq::from_str("19930525972818199875851038877132322547334829658949554806801108404941414101103").unwrap()), ell_vw: Fq2::new(Fq::from_str("6679449617056144648556013688092503718528081567224358397706944806960707450866").unwrap(), Fq::from_str("11834339054115289581728384446510154030614495257639549098914576992681792644317").unwrap()), ell_vv: Fq2::new(Fq::from_str("1103971970080878735020200145931068396548255816344141900100578642050167835756").unwrap(), Fq::from_str("17442014287618549059909689032626628911431121434337187918457502534017974292771").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("13076805492928497748141329405065621544173515733010216253264862326898999186011").unwrap(), Fq::from_str("19346998537378535293491878821381642632408077348855336790184043880601230192418").unwrap()), ell_vw: Fq2::new(Fq::from_str("2545953185107430710067715026283426197979102699516436299879612082830465064387").unwrap(), Fq::from_str("658480144944520773976742499043453183928200527625254756411898212059863027220").unwrap()), ell_vv: Fq2::new(Fq::from_str("10599092041977168661053673645861514974472253862403262653806335823439984547520").unwrap(), Fq::from_str("17162893574285047796735108280435654463748265136801600105571927738662909784323").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("7085319511141515608176313666085289847826695130590757892136226409211984372940").unwrap(), Fq::from_str("4640142023499829801003091193578633557246618475866500623404375352964089302137").unwrap()), ell_vw: Fq2::new(Fq::from_str("1378711095367589966539356394028311641404402942969200900226036867512926984525").unwrap(), Fq::from_str("20227188036459114276006891788262806649925222994898646338447622388573680807347").unwrap()), ell_vv: Fq2::new(Fq::from_str("10909351851886459521617902618958119972478357132140936115187854789408351897716").unwrap(), Fq::from_str("4270359828288113001419972416654771076387260764523167288282363565700959902034").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("10664974041382313497332907041293813168318012018754207087217051246266695094381").unwrap(), Fq::from_str("7224602128989495563505314962160716046811485526332201590284671097356858321790").unwrap()), ell_vw: Fq2::new(Fq::from_str("3934655983500232748821065857705540333423053242147172231346384664554964864345").unwrap(), Fq::from_str("11454752219285395886386371667205283572988335397450041936477253847740701967591").unwrap()), ell_vv: Fq2::new(Fq::from_str("8613140983903038954690786727135433848146709090100659195246608245554768285578").unwrap(), Fq::from_str("11683437859485859733100235618517710813205908191940301821880216475492534630964").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("9260311336006363265451970156091146121161647548618044834929651955256342089015").unwrap(), Fq::from_str("13937163080647677520680042908210664749449075029135461146729576039103333691760").unwrap()), ell_vw: Fq2::new(Fq::from_str("17902331142943185999552054453100857674678810044830315483268693319046882556841").unwrap(), Fq::from_str("1544882902322119030227919062381439448920928851067359863906156836436873153173").unwrap()), ell_vv: Fq2::new(Fq::from_str("16995331248700240765237826468302734126993772776070338276810983166039248498016").unwrap(), Fq::from_str("14969634758504522458179971123440531767135906355134596040598769037804654587702").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("9691457346326988270440356562380914538881384898220972529069554878570918299675").unwrap(), Fq::from_str("4212749851190062070411861546581519573149862920541429109749681953529770836138").unwrap()), ell_vw: Fq2::new(Fq::from_str("10890917529063563080488636605724946896609903254351684402546841801646043334166").unwrap(), Fq::from_str("17064796604657821037034221187126605208939910043822292844975397354283757360178").unwrap()), ell_vv: Fq2::new(Fq::from_str("10053764157025690527213061522668163518874554311426109722786780729745637570153").unwrap(), Fq::from_str("12784239716790558610229592219566905741360808124013978959238298462172061891010").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("18178610037180803777372696468852749823139836169202607246443680867419301660001").unwrap(), Fq::from_str("2203553150397509140176760042861998991151797164185565920046676990283238338489").unwrap()), ell_vw: Fq2::new(Fq::from_str("20994248737185978263276750124680792452522185817020429107744966692073225126536").unwrap(), Fq::from_str("6349838002511478763817816496894821255234760390494115897123421446389963162622").unwrap()), ell_vv: Fq2::new(Fq::from_str("1263928993764270656816573911502202456032289420152723320064785370753332332564").unwrap(), Fq::from_str("3046827108505433561924898780283267608115953987844810798949064271834105292424").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("8956928906602747002792065663685826131559810003765184492377190892700326559543").unwrap(), Fq::from_str("9256558483224710959389211180674778886236784676421453067453985103303008910100").unwrap()), ell_vw: Fq2::new(Fq::from_str("19258483215045500292709619645410473292886118206783068420253408398264757067719").unwrap(), Fq::from_str("4451865829956913287893327659699467549695832858717473490066802599791690802339").unwrap()), ell_vv: Fq2::new(Fq::from_str("2764660330371796211645912465865514165647712230266940655139955982219479135487").unwrap(), Fq::from_str("13240152986607724183470740198192310673584785226612774164092563912390438227662").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("681767682327814473892054649470454096768283412403007527891528077929769397010").unwrap(), Fq::from_str("7668193584311505399758378551993462119934308304509766677185196010163712481402").unwrap()), ell_vw: Fq2::new(Fq::from_str("5655243760535315520783349710797795360471266639432050916814563868068635911845").unwrap(), Fq::from_str("6563255051651012450772847952915771842363624186837705815170827122750163703821").unwrap()), ell_vv: Fq2::new(Fq::from_str("649940975249127008407571724305690292046357719924151085455652491670338045916").unwrap(), Fq::from_str("21579311737676937819033214124612551530127167774002946532496743344036322675381").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("21175778737847401044104466077056088822729568143966863788586142950892805137666").unwrap(), Fq::from_str("5595589116740977138584720715075423123216685509908153451278598090500673116809").unwrap()), ell_vw: Fq2::new(Fq::from_str("10584200644032922607821101570277578933623023426596682657073390742969073652895").unwrap(), Fq::from_str("12837953697025504547863231066503572653998301486985427806275987315834659812986").unwrap()), ell_vv: Fq2::new(Fq::from_str("8117107548208718171682922744254403704391859020172151233388364110731129041045").unwrap(), Fq::from_str("13049999463955485268550904857465746325902043660001708636685369259450564610454").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("8364961604780202742303921161081604160922838515374423608219736048018999860853").unwrap(), Fq::from_str("2351019119583633406807532653396626109486234484023597064733508979159842183382").unwrap()), ell_vw: Fq2::new(Fq::from_str("21028432880575926523752559333294402827459370060933557257525967464997077946986").unwrap(), Fq::from_str("8046872296156140750811104526560806223677002170113323621491280498583760192616").unwrap()), ell_vv: Fq2::new(Fq::from_str("12898024833491487058567514428449740829256150305937499595589780150928565842102").unwrap(), Fq::from_str("4285065448422006185829060466560686578745712460867295748911655362509144300948").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("8072232359889306936722303478835558056241116613774315547062190240728474923834").unwrap(), Fq::from_str("11177008054227634237663439894104913366867843529726469624208391026893438331544").unwrap()), ell_vw: Fq2::new(Fq::from_str("6809080759521318546573252139200749353665236974149238674422052122591837770337").unwrap(), Fq::from_str("17038522274461048605370004543286094433565235245947661777349160189041706387246").unwrap()), ell_vv: Fq2::new(Fq::from_str("14024834693814357899802394991122512208020612837103453479716792821268632856005").unwrap(), Fq::from_str("12162253519201182595066522680875457472054933146636012038079782236991342232964").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("7517621744472550238729147112582763290619696603888340562320661653468259785459").unwrap(), Fq::from_str("5775958628431902009928696648706206677051589697801230225207961079701351805991").unwrap()), ell_vw: Fq2::new(Fq::from_str("17101537508347637393046612945060881733341630882242907674489426818282232192765").unwrap(), Fq::from_str("11439924192045730134857855444951867487673481853061712676505751821713592417848").unwrap()), ell_vv: Fq2::new(Fq::from_str("16524831954193049226680237913906512849861747525312574297009556588583016003485").unwrap(), Fq::from_str("17961833510100001957206410577764797351002903290115563676683242640020206248673").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("16912575460065719860154998920476780598787708817030774617128864437198545819929").unwrap(), Fq::from_str("18523390726102623209642393197047309169880798198922216024709291880134883505164").unwrap()), ell_vw: Fq2::new(Fq::from_str("8711744315121016096444521009902336544777137803294883395587631628756743226403").unwrap(), Fq::from_str("384665796631718061484470820770858149414181395445222225443039170963661807931").unwrap()), ell_vv: Fq2::new(Fq::from_str("15774538322747462731385931838340979144243491164741595437453112329195248618714").unwrap(), Fq::from_str("11632275451601229685994821619837954138410931366643826111939350416539675376582").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("4193777370819241345875098285464776579390947200963349948619027946434227330979").unwrap(), Fq::from_str("509866391215626792725431227507051906680902316754753085884521978793646270092").unwrap()), ell_vw: Fq2::new(Fq::from_str("6037175194865804429042789680410737585033925490367655507316663717551980954540").unwrap(), Fq::from_str("18041663784750649636833746105782119188697038394822653258671508993118228507622").unwrap()), ell_vv: Fq2::new(Fq::from_str("5823434994217288412340362524966578492063846889148854293277163714066988854109").unwrap(), Fq::from_str("5672396008289604057966364197930061596003079754530820851713923211431895221182").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("4038166397773785982077422827894979131906966567246129821324959462204170431412").unwrap(), Fq::from_str("18437550501813119699199412531646289403507307693226141098858673836744295275291").unwrap()), ell_vw: Fq2::new(Fq::from_str("14284923653983275338299458156785711701885868498245540450770121796221847779083").unwrap(), Fq::from_str("12602897944901625024952695280664364105484904418233840369524699062991946110320").unwrap()), ell_vv: Fq2::new(Fq::from_str("5486498814295375008109094800710678785628460899485721363575459728835468025031").unwrap(), Fq::from_str("16200691130190012476008796365701719097160136918178762563260333528187553865070").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("12314546787967340174315157919495563500551846522430770700093005979883034971431").unwrap(), Fq::from_str("3526449836226747514409509674529965207620312182683788709671366453353950282439").unwrap()), ell_vw: Fq2::new(Fq::from_str("5093273668274685698542824313074411944254604209522407110762874543698544425100").unwrap(), Fq::from_str("14413934995969972893982204347979972725062293812449032149868616703831967954701").unwrap()), ell_vv: Fq2::new(Fq::from_str("5926816361020227072352913556957354039762675747717706660663155084741002986207").unwrap(), Fq::from_str("11598493987133778453993847174276604684683441101397418085992950965674684594056").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("14605016886275487293555440005295724730930746064070181829857167303747503806723").unwrap(), Fq::from_str("125488028578519699988648864181559879827394253209756997049515479109363959656").unwrap()), ell_vw: Fq2::new(Fq::from_str("13618798589781916917754240943924749450694915714095829562214241999969133303075").unwrap(), Fq::from_str("7919521676249934186275057369006254600339974917319129664598969307290994883444").unwrap()), ell_vv: Fq2::new(Fq::from_str("13618646464444173552701873510284352084626546810748938598145778229202447557390").unwrap(), Fq::from_str("9476574581262975213029808591741428214940842553330900266594833123385626035464").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("20330569420634759863618802200434379756882301920599312733161738334139846403423").unwrap(), Fq::from_str("16941424680188140264426604507849243567044731815188670894126312258606055646954").unwrap()), ell_vw: Fq2::new(Fq::from_str("4276245652046353558826181153633997046897244969830610112857888938076187872093").unwrap(), Fq::from_str("5218134511825218024578456539758678791071455837562669069391138178014247652283").unwrap()), ell_vv: Fq2::new(Fq::from_str("5074559277335429630740829563679793620243538875458585401599810469948320707914").unwrap(), Fq::from_str("8399181532479234777915440768158248538072911084393374644710065479068982458679").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("228456190514225889641604150077307010131630607092226273786181288134202785827").unwrap(), Fq::from_str("6829691501727142726895142554058769478742020003540473851444113314854839779896").unwrap()), ell_vw: Fq2::new(Fq::from_str("7239124420778199901449006013253700083679198792863456824897986348225520969318").unwrap(), Fq::from_str("9612815502366829000381428453387330590336692138445456863004983078118510731327").unwrap()), ell_vv: Fq2::new(Fq::from_str("3032095333965417803301146592958869325973325547652138997619203744855910423626").unwrap(), Fq::from_str("9140900126707960555227810341717821370276432039323462635381934661102939064717").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("15183121839690585327918066578070528787848569635478652970262274682157344344126").unwrap(), Fq::from_str("17032257654469593612895387885450704824985120918227803531218149228038247116993").unwrap()), ell_vw: Fq2::new(Fq::from_str("4449946720064992614521464039848704811573349967863291193814384877976473528180").unwrap(), Fq::from_str("16213473115750238158840972181958630393893069174331658126086552358850995045944").unwrap()), ell_vv: Fq2::new(Fq::from_str("14820364816206160787218030695278797219740217739100995460398841348849807343630").unwrap(), Fq::from_str("13391831978082769659226335015618034703517552702532093599836056107584855920091").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("19376350172596821938106071755601309991977917109030961898807003698880067029751").unwrap(), Fq::from_str("8432677549750770318577385075581043130355856390772751886981628666318547719067").unwrap()), ell_vw: Fq2::new(Fq::from_str("19907801542387995008974808935818119635245812040205527769417146511798917761987").unwrap(), Fq::from_str("5161102563798536355551540455475021097419898845343559944869423371507274303876").unwrap()), ell_vv: Fq2::new(Fq::from_str("8934093116311670296805867189709376457342708648536176030347455800639124491789").unwrap(), Fq::from_str("2703243001408501823489300103176986467031566088207792046257025791153593204110").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("14246576411085900561418149065736187542613160942023368670593924207716148262147").unwrap(), Fq::from_str("3894319734417282749436553474309402307762162253367446669496582090013706055610").unwrap()), ell_vw: Fq2::new(Fq::from_str("8375748526029577467798965722696242475050467009610481042928314348040697238273").unwrap(), Fq::from_str("5150547143250446746041060583106242553976637956855972603954741388536316657102").unwrap()), ell_vv: Fq2::new(Fq::from_str("15907790590424017252473639606980913201030067954178224590289142061847439584806").unwrap(), Fq::from_str("11159987404206640511954116561264499629862056981304725791275082386203982274572").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("8534750176510034607784585561677263909151653570221826711790070724427621935471").unwrap(), Fq::from_str("7817152482056867378682253092041139659293041868731996264204246731388352069129").unwrap()), ell_vw: Fq2::new(Fq::from_str("16961113237531514134058762308816541864350643273034658785830240248900850618232").unwrap(), Fq::from_str("1179078000730113993104369721112852911420394064162245727296751799866128074554").unwrap()), ell_vv: Fq2::new(Fq::from_str("20457374714339105220680185601375259716132143552962630535477097129353077443649").unwrap(), Fq::from_str("4422525926709336878773825319073188719230418589270084994878618648788354576429").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("7627847586963893527076371866753343482026765205726665650559945926181618692545").unwrap(), Fq::from_str("11854660093187915385325641424610609952992715573834872602540776698350931163771").unwrap()), ell_vw: Fq2::new(Fq::from_str("12061802008303461437760636908585485318885385024493098952040511326354153801540").unwrap(), Fq::from_str("2483328707802038203960641185873859652602098605879879702928670691083351376828").unwrap()), ell_vv: Fq2::new(Fq::from_str("3179938981834043637934998748622377818001614400201120097155888429286598162788").unwrap(), Fq::from_str("3387712903817602592588061982323804658595985412237324988504760476608661908669").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("18935668687475417752670528265758371887227072908926123795139093572891409166177").unwrap(), Fq::from_str("7137948322161355685471108734802366001522057971184080351365620151904917529667").unwrap()), ell_vw: Fq2::new(Fq::from_str("4339156032094542567838794153972725732752104151129467898887880650450609600794").unwrap(), Fq::from_str("12048671216064877807276853910375839173096391682288251978695349073893797560105").unwrap()), ell_vv: Fq2::new(Fq::from_str("1950772803812480571427275248471639178971170644950730072787717123474197730954").unwrap(), Fq::from_str("9181624035531055165333621510604231877029437698520448252367489746796255738472").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("11063654615400954717295218153285806614049542407631184476400061600226642090639").unwrap(), Fq::from_str("10322234663627758865726091502227262577963762851829427143908452375568424940506").unwrap()), ell_vw: Fq2::new(Fq::from_str("5088512375045637070831200431021199523774202712519349524231646029565850243137").unwrap(), Fq::from_str("13830164350727008423053294968318643994288249825064419288902919319598527803018").unwrap()), ell_vv: Fq2::new(Fq::from_str("5664114633031170805308124108402506179037402200165762737765348606356467981435").unwrap(), Fq::from_str("6420500289322096738785425211380588437700105870048213919996474830908364662080").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("12089017351045163205876906516764262596252750496306969979802152087250767442655").unwrap(), Fq::from_str("10210558808216547366164047615754858102665922901025755697993978799086515582216").unwrap()), ell_vw: Fq2::new(Fq::from_str("5805411494371431706302702029555168082350433630263841912576207114214821637903").unwrap(), Fq::from_str("13417649473832328299363309530328161775713579734460687495589323668315481789186").unwrap()), ell_vv: Fq2::new(Fq::from_str("7100161023424716303226449735068569399024504614273633218656928626719841360262").unwrap(), Fq::from_str("4608780226874321006375565633722017610093139802788639970961143211665122970425").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("7402656885283645637844399867437832414708002274719716809314550750025200002816").unwrap(), Fq::from_str("3359755418859898733774484791668766343777859472186478368952701397946858170923").unwrap()), ell_vw: Fq2::new(Fq::from_str("13972633490402551449646647488711159435500678729477306484155919447807727623768").unwrap(), Fq::from_str("7493072689132118410451346570771015894599566093753023472904593233187191189290").unwrap()), ell_vv: Fq2::new(Fq::from_str("10553632499535499938294059057022861371325444554911689769145218759723376859522").unwrap(), Fq::from_str("15031724096959569889887883166930429112400023614133986918505200177949730145621").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("19614666780105423133718014776333000423285968913902974939878172716992608453809").unwrap(), Fq::from_str("11603995082765755677014543451017101337084691386341657702599446172204728703168").unwrap()), ell_vw: Fq2::new(Fq::from_str("3640285626277218751703647194759635465446306209869384174139659856814138047153").unwrap(), Fq::from_str("19852827152161121673942769427296662914445783006570890136096053338872958018957").unwrap()), ell_vv: Fq2::new(Fq::from_str("15133277841235179795557173122463348349264603847086033764926472096701308210280").unwrap(), Fq::from_str("11477545358800195084620539652687203979589770655971679975137907269158992047904").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("14110112013826130627308137057440352762673176006101512326319632538666268095190").unwrap(), Fq::from_str("11939502323738533476551128890683182628994629048369580067812328846012456410750").unwrap()), ell_vw: Fq2::new(Fq::from_str("4903570904021379863440412383519324112841733022665060696302943089410552325193").unwrap(), Fq::from_str("14266803156893138980425811506240396635009593456307449490727186191181005152878").unwrap()), ell_vv: Fq2::new(Fq::from_str("14881072700587137091207448172180518390500891251946414292329369965686915854422").unwrap(), Fq::from_str("15926908436628594399583920130944357569706863129497064862365363758969852348774").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("12558706151877569296451844961986237406953575964353502214928234165049347203902").unwrap(), Fq::from_str("20986351119648698768747274174278231882270962658108929786748672594878002377020").unwrap()), ell_vw: Fq2::new(Fq::from_str("19049669232824485247328064933680299726722169816787419902393190195363675139966").unwrap(), Fq::from_str("9038064204016446065702766858179381867420605337773795708435102200352238225734").unwrap()), ell_vv: Fq2::new(Fq::from_str("20935836360010341957040867041308346409399314465762274744792154001186771274330").unwrap(), Fq::from_str("10203406114258970224185400536984459120400171034941205452745624846762632193745").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("5273485321155899320017031787320374433252287213043632779466667239152835002841").unwrap(), Fq::from_str("3053056741460414900892697870729818263336633719763669407853703091660209931803").unwrap()), ell_vw: Fq2::new(Fq::from_str("14340825980720147789753139378133767339849529468405152905403690972933175947523").unwrap(), Fq::from_str("10774533501875920615470512566140848784460524907884011188458749902565676606519").unwrap()), ell_vv: Fq2::new(Fq::from_str("3964278893337662256819152335242286127722672026565541866237058714197501521533").unwrap(), Fq::from_str("18022244127367262334687067361628890804799640514175862613915842388966614841200").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("16332502873726000727943290608233946425803707374189230473690022986938005644703").unwrap(), Fq::from_str("1567962947294031055122165158764344785171646063465144808831306664719326769256").unwrap()), ell_vw: Fq2::new(Fq::from_str("20252166217271346908492216031557669988560255907738941815783926955730505937199").unwrap(), Fq::from_str("17117796187826532609679851805915085724672631310980423533710032156677468384570").unwrap()), ell_vv: Fq2::new(Fq::from_str("17447467886258986370787137061379217598255780655216428459333948972673151012677").unwrap(), Fq::from_str("5126879980313476901131854669269919146638602687291477011061931179530480305128").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("6088970897125028663727223015788851329277811231049906821835210021445230531470").unwrap(), Fq::from_str("8269042454443245208032034843063473106607747629991910601236499943382389383909").unwrap()), ell_vw: Fq2::new(Fq::from_str("9390988514296458582044517134794855627538985266829573142624063333360969055582").unwrap(), Fq::from_str("14809192622791971819542419564616657251564340618078775692251491022692953050806").unwrap()), ell_vv: Fq2::new(Fq::from_str("12392403118986217155446218801253632414719473065843423820454574775060624572223").unwrap(), Fq::from_str("17315648550786131512187415565477299679361071890287921321719420341945848525344").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("8509704393616883918389461491307625357351608849562781376448108143686881397874").unwrap(), Fq::from_str("20273788844812669168701278525056990767412251703786041383351127813548056060175").unwrap()), ell_vw: Fq2::new(Fq::from_str("13587733380957601917144020513037833512715927912961534242274129611676082179670").unwrap(), Fq::from_str("18219007782547538891755131684432302965729501555658190323551992733355994409648").unwrap()), ell_vv: Fq2::new(Fq::from_str("15212235563035911012497734239429958536448576200220918241945992532007455729638").unwrap(), Fq::from_str("15868967159924843704138451124239001351258676891119605450618238262537266256000").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("566642944339341838872659990670775497018813579200118662738575067585620262683").unwrap(), Fq::from_str("80633561422247315402489148783136836267866433112423246780090434639448081501").unwrap()), ell_vw: Fq2::new(Fq::from_str("16776314135335488452712466347201660236874408701359052824185479604688604443386").unwrap(), Fq::from_str("9397286871750513195379912820117368990110423949006393550333044521067916303673").unwrap()), ell_vv: Fq2::new(Fq::from_str("11213241389292792759942014406027965770431479614859139543821417488303149467800").unwrap(), Fq::from_str("13968940593285273772185768016764839518710968380183853737025067101492115101555").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("6582524148524166281651049223972241418275829952930531957133089376081226599310").unwrap(), Fq::from_str("12912796945237135113112612780707994091197804899944634993937602107108817617069").unwrap()), ell_vw: Fq2::new(Fq::from_str("5915143763477343907746278070838860419057837978194254349309930636588762054019").unwrap(), Fq::from_str("16829580672314341143130263740133292836202318062773617917928807522786048686905").unwrap()), ell_vv: Fq2::new(Fq::from_str("6048262046277135430490078693604752857941808479983722633615831914037376701617").unwrap(), Fq::from_str("11339265196645696685132629437938382050130007247287625382019491127986659567508").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("16221040390720130701106818492787718054167551387591190581820867019772199748452").unwrap(), Fq::from_str("13816650124538024294513769170531808824213976236314395631188393480750356347832").unwrap()), ell_vw: Fq2::new(Fq::from_str("17884725815671017389514271919761437310258619463085199936932965564030171526767").unwrap(), Fq::from_str("18248973796984942533601080062907385035195180110558182226420098600644022965345").unwrap()), ell_vv: Fq2::new(Fq::from_str("3137416660392142651716799414483848222789075294982695366740419822188364577896").unwrap(), Fq::from_str("9694367271685190541749430052152010428005565838702388077250686840743397580564").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("6357610272973630149840225505274539351850867930627803566656124705077309138236").unwrap(), Fq::from_str("6532064867969562545750502093085673423036025281166230602025667552530979558092").unwrap()), ell_vw: Fq2::new(Fq::from_str("16786670094068652274677750111730144379312178408732829954083928254919389983983").unwrap(), Fq::from_str("6341529030893580474995202784328098673009525246309131907462803539077246776781").unwrap()), ell_vv: Fq2::new(Fq::from_str("1974498926273096258802796702119259182018579009158122203612017565775578915815").unwrap(), Fq::from_str("16202371700345289299248140495735793174324620209466889736695203996551421837181").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("7645596703863049949671625140189098533219782764649960681401839437305310089346").unwrap(), Fq::from_str("19662686877039093771787921401379726144714150453645437909640334762978238819386").unwrap()), ell_vw: Fq2::new(Fq::from_str("19675517233991845002128733409028273554556593365253876728982552859262257722391").unwrap(), Fq::from_str("17810553502824017140555114404262914434390762266054795883352268646885364159410").unwrap()), ell_vv: Fq2::new(Fq::from_str("18665072545928024509881282829034535118797458683193169885406548536617293330990").unwrap(), Fq::from_str("19411595053923298782146703680748498557982136082135126508584153240357549362078").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("14470514611309746071605255274946119202484366658141275541087020842175954067559").unwrap(), Fq::from_str("2594897518684998941907190849745518961756699273476562463311001408776339658780").unwrap()), ell_vw: Fq2::new(Fq::from_str("935873152179726234311073616163680460286310013909076125956427156452613025181").unwrap(), Fq::from_str("1833293684395168387765260233178732135932817355166697645828433012629591012612").unwrap()), ell_vv: Fq2::new(Fq::from_str("15506344864201849985196152671178663605696906693236430756412539843011256851197").unwrap(), Fq::from_str("21249626787247514663459670711353542343321937134901366527096486234949436840608").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("14844298589371811580114643149994926627935566419885127313767249060440840193460").unwrap(), Fq::from_str("21865492910219121705480174226763323783866193080747061232287801827306146286946").unwrap()), ell_vw: Fq2::new(Fq::from_str("17151100973274235243640792655189642956800303990602164052684336481398298289411").unwrap(), Fq::from_str("14392190221562986523826107157981466944764379735735544220135185056395768909930").unwrap()), ell_vv: Fq2::new(Fq::from_str("4350255565215767538405433120041486620343663937357318680211203103336875636557").unwrap(), Fq::from_str("9926081950735179732148156026738993553156731124330466866635273595752973758855").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("13218761129163938072615510267349114476398845757988353548581904005855607494827").unwrap(), Fq::from_str("5599121822325202341068442750237639791575579072233857524815548719976323285948").unwrap()), ell_vw: Fq2::new(Fq::from_str("15971120376538832228790474763995956330325241100841265267806125625593561200292").unwrap(), Fq::from_str("12523288625433720594713119420870235433458343202664663107083051275205911332632").unwrap()), ell_vv: Fq2::new(Fq::from_str("16334648365696506264242655102949242813244449705193389318653840069223698505218").unwrap(), Fq::from_str("5656148484258838154513173637616383061988381726003070282394323376702835298685").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("17813738829828434290335491956960173922453516626523775773406434062574864455163").unwrap(), Fq::from_str("20602369742562270094694203910247000789225384961205618501257582838246099644397").unwrap()), ell_vw: Fq2::new(Fq::from_str("19878549521928480892919729627745688905048802531306612605469432942519616386562").unwrap(), Fq::from_str("5017457857497268504225948371697060380418757800650743056471872023476645133890").unwrap()), ell_vv: Fq2::new(Fq::from_str("11253669115821132651967277438747600032147012598164512191026558281212405907059").unwrap(), Fq::from_str("17527101633359473938082679669763197242520773041737623688955597341205608093817").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("4215545173967763533303120784928042872427929334572670490629936371429427871071").unwrap(), Fq::from_str("12294795716682648009432653500746413304750644076589949908249345884220200424705").unwrap()), ell_vw: Fq2::new(Fq::from_str("18825720171660950734491494575282684507875213824417353005009465489701558195915").unwrap(), Fq::from_str("16705960611772713893529695679517111622754789317070946155202186944060021844706").unwrap()), ell_vv: Fq2::new(Fq::from_str("727717707174551758531261221179545916153780087179478433501944451338651779512").unwrap(), Fq::from_str("8112271980906996458859937804680221894935591766524365848766530295751121904748").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("10233024162783159141952008177866660607852287510058437429941412641351833243474").unwrap(), Fq::from_str("9594838646894429694778293318641489517645344438331800828882574864439372872574").unwrap()), ell_vw: Fq2::new(Fq::from_str("16645682756797429615195184700295338945089747965333104334845226609564841400023").unwrap(), Fq::from_str("4162352035394629024812623559649906911351409190415492302470361657064711946691").unwrap()), ell_vv: Fq2::new(Fq::from_str("6032119528092518381867266123093895943802103316999621961287699634448441322294").unwrap(), Fq::from_str("16978590002913514123461500015923670605519520926487110857783653291799516315656").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("8218196524386055233961282665657821183883422738571989211183480610420112668380").unwrap(), Fq::from_str("10415506866745570162620028320970394286106067432626512901854031549920339603997").unwrap()), ell_vw: Fq2::new(Fq::from_str("5210044150482547744429104832217210823770891140961861953610847155896342519957").unwrap(), Fq::from_str("7249576834325689075837879557249227252903732278842727723303356548006299693332").unwrap()), ell_vv: Fq2::new(Fq::from_str("10566250213397026997176223676662431109660830193211672199901829008963029461821").unwrap(), Fq::from_str("10109726446308488117042758961222375921166426648339199001670753736584698407831").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("14958658301977599538003236672511296104101995894564376511283101560012949683335").unwrap(), Fq::from_str("3048608469033125363311089547916169290720115268350540694596412955569118898494").unwrap()), ell_vw: Fq2::new(Fq::from_str("574489666015527527356241753990038643397861257676315289742492198353125347063").unwrap(), Fq::from_str("20135743140323312729475654806927894483883158119600963784497647679670289292384").unwrap()), ell_vv: Fq2::new(Fq::from_str("429037282750459454100296565403648201904295120071497517313998879109362521423").unwrap(), Fq::from_str("15853233349247701874380259741343021520320591334343875702627292292088300425145").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("18700030680590279184945769912899822403041358764323477965514340815228877402331").unwrap(), Fq::from_str("10683998632701796443517414309938660369292269662497518738662465524924019738910").unwrap()), ell_vw: Fq2::new(Fq::from_str("1310694906110452150099267122613851606459600617516350757492763828252192371857").unwrap(), Fq::from_str("4783076385012654125132398886450190041207227400254467276271580769333335151361").unwrap()), ell_vv: Fq2::new(Fq::from_str("19588709046778297839342456362815786766378666302331732601204786295010565114246").unwrap(), Fq::from_str("15768094978512846661820502043664035262466161785979259300748269497441797675546").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("13060364876687824737412369705619641318084803694180161837383771909738942399322").unwrap(), Fq::from_str("18217552047682533290011237332328856250732338018560605204821714348687071938690").unwrap()), ell_vw: Fq2::new(Fq::from_str("12869661643727171099080271622741099816758098286424339832028838960253051261267").unwrap(), Fq::from_str("9771005613420764951919925720405408601257471427085946188171574196917272408972").unwrap()), ell_vv: Fq2::new(Fq::from_str("50909581144623642750375885935872654509592339539711209410079054797052409354").unwrap(), Fq::from_str("18926528403466649263993967111065644144239906211315737390219036362821972145669").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("20549804669307980165961565026567668900175984309510667407934326603524068980568").unwrap(), Fq::from_str("2860891967505015044048360806986650576604923893849408592665407858213623973073").unwrap()), ell_vw: Fq2::new(Fq::from_str("15425380359872501526721873883293811604657689327716748473720120299487795555723").unwrap(), Fq::from_str("19476618728903347368164774764554105623361814261939990815870207909384053172994").unwrap()), ell_vv: Fq2::new(Fq::from_str("1391900673866549134877070137078952075347059393013450066057338968238992587955").unwrap(), Fq::from_str("4196316861965319041187484774757427170509633759722370813890904092023663955821").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("17275387403896977778123657295201685204905269368879700703438840824610915049467").unwrap(), Fq::from_str("19208898840130004899976179132541534722147305748950927170013798354273078356439").unwrap()), ell_vw: Fq2::new(Fq::from_str("2806476769978795611970161315323438770892084130100332531802098991493644141407").unwrap(), Fq::from_str("12027685137795878129532114449809592645364590672303716160923491816560315933093").unwrap()), ell_vv: Fq2::new(Fq::from_str("8236588507133996535585623995380369531792969451558421380293164596516832310024").unwrap(), Fq::from_str("7194816002224802847462171990762514812770222507636791853795854276216118633532").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("8988630595018387496849477785798469553390006328427704582362113137927154840355").unwrap(), Fq::from_str("2295777524525303845671711839319946449143831320846022233249545363172788331149").unwrap()), ell_vw: Fq2::new(Fq::from_str("19856063720990800344827304970869772811632532243020577718537161398914360273978").unwrap(), Fq::from_str("19822454791015240925652351664190735476954892251509376470361602440577337517739").unwrap()), ell_vv: Fq2::new(Fq::from_str("16437024240432855113321845721989731630224984534156413407578941758236945886728").unwrap(), Fq::from_str("10915573623321985107029933965796058231484519510670830600808943167010577314334").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("21286975582879196598650077668567283554136994229408922053641846042051746388886").unwrap(), Fq::from_str("5421529658065143080970166226046014500038127295988164061860891929414040927615").unwrap()), ell_vw: Fq2::new(Fq::from_str("5555485154537896512520569992616892769357525491868135294451583184259218193278").unwrap(), Fq::from_str("4831809465565304340310880868620257258305484211486308136368207144444325564275").unwrap()), ell_vv: Fq2::new(Fq::from_str("2583089702666220608189326901976548985583885001543445628170288810575076369034").unwrap(), Fq::from_str("7097788099664798415333498763272677734135137637990913349021964143797057593727").unwrap()) },
            EllCoeffs { ell_0: Fq2::new(Fq::from_str("21356455346146359274071336488740257199607341582407258426849457339859865454145").unwrap(), Fq::from_str("19512280220246743343068375192508909480135195681362473260894831078332945004959").unwrap()), ell_vw: Fq2::new(Fq::from_str("15435677241338788169971902812985017500898533325271813422525238374746874853456").unwrap(), Fq::from_str("6382067580057986766020919077216175344977116334133202456947134195629074647989").unwrap()), ell_vv: Fq2::new(Fq::from_str("1150360546643207487285835972411275438312694180812115547021545502055712764110").unwrap(), Fq::from_str("2194960903427600986795419731825236689333014397679688012527722471213485607323").unwrap()) }
        ]
    };

    assert!(expected_g2_p == g2_p);
    assert!(expected_g2_p.coeffs.len() == 102);
}

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
    use fields::Fq6;

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
