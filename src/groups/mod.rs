use std::ops::{Add,Sub,Neg,Mul};
use fields::{FieldElement, Fq, Fq2, Fr, const_fp};
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
}

pub struct G<P: GroupParams> {
    x: P::Base,
    y: P::Base,
    z: P::Base
}

pub struct AffineG<P: GroupParams> {
    x: P::Base,
    y: P::Base
}

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
    fn to_affine(&self) -> Option<AffineG<P>> {
        if self.z.is_zero() {
            None
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
    fn to_jacobian(&self) -> G<P> {
        G {
            x: self.x,
            y: self.y,
            z: P::Base::one()
        }
    }
}

impl<P: GroupParams> Encodable for G<P> {
    fn encode<S: Encoder>(&self, s: &mut S) -> Result<(), S::Error> {
        match self.to_affine() {
            None => {
                let l: u8 = 0;
                try!(l.encode(s));
            },
            Some(p) => {
                let l: u8 = 4;
                try!(l.encode(s));
                try!(p.x.encode(s));
                try!(p.y.encode(s));
            }
        }

        Ok(())
    }
}

impl<P: GroupParams> Decodable for G<P> {
    fn decode<S: Decoder>(s: &mut S) -> Result<G<P>, S::Error> {
        let l = try!(u8::decode(s));
        if l == 0 {
            Ok(G::zero())
        } else if l == 4 {
            let x = try!(P::Base::decode(s));
            let y = try!(P::Base::decode(s));

            // y^2 = x^3 + b
            if y.squared() == (x.squared() * x) + P::coeff_b() {
                Ok(G {
                    x: x,
                    y: y,
                    z: P::Base::one()
                })
            } else {
                Err(s.error("point is not on the curve"))
            }
        } else {
            Err(s.error("invalid leading byte for uncompressed group element"))
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
            y: const_fp([0xa6ba871b8b1e1b3a, 0x14f1d651eb8e167b, 0xccdd46def0f28c58, 0x1c14ef83340fbe5e]),
            z: Fq::one()
        }
    }

    fn coeff_b() -> Fq {
        const_fp([0x7a17caa950ad28d7, 0x1f6ac17ae15521b9, 0x334bea4e696bd284, 0x2a1f6744ce179d8e])
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
                const_fp([0x8e83b5d102bc2026, 0xdceb1935497b0172, 0xfbb8264797811adf, 0x19573841af96503b]),
                const_fp([0xafb4737da84c6140, 0x6043dd5a5802d8c4, 0x09e950fc52a02f86, 0x14fef0833aea7b6b])
            ),
            y: Fq2::new(
                const_fp([0x619dfa9d886be9f6, 0xfe7fd297f59e9b78, 0xff9e1a62231b7dfe, 0x28fd7eebae9e4206]),
                const_fp([0x64095b56c71856ee, 0xdc57f922327d3cbb, 0x55f935be33351076, 0x0da4a0e693fd6482])
            ),
            z: Fq2::one()
        }
    }

    fn coeff_b() -> Fq2 {
        Fq2::new(
            const_fp([0x3bf938e377b802a8, 0x020b1b273633535d, 0x26b7edf049755260, 0x2514c6324384a86d]),
            const_fp([0x38e7ecccd1dcff67, 0x65f0b37d93ce0d3e, 0xd749d0dd22ac00aa, 0x0141b9ce4a688d4d])
        )
    }
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
