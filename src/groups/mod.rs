use std::ops::{Add,Sub,Neg,Mul};
use fields::{FieldElement, Fq, Fq2, Fr, const_fp};
use arith::U256;
use std::fmt;
use rand::Rng;

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
    type Base: FieldElement;

    fn name() -> &'static str;
    fn one() -> Jacobian<Self>;
}

pub struct Jacobian<P: GroupParams> {
    x: P::Base,
    y: P::Base,
    z: P::Base
}

impl<P: GroupParams> fmt::Debug for Jacobian<P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}({:?}, {:?}, {:?})", P::name(), self.x, self.y, self.z)
    }
}

impl<P: GroupParams> Clone for Jacobian<P> {
    fn clone(&self) -> Self {
        Jacobian {
            x: self.x,
            y: self.y,
            z: self.z
        }
    }
}
impl<P: GroupParams> Copy for Jacobian<P> {}

impl<P: GroupParams> PartialEq for Jacobian<P> {
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
impl<P: GroupParams> Eq for Jacobian<P> { }

impl<P: GroupParams> GroupElement for Jacobian<P> {
    fn zero() -> Self {
        Jacobian {
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

        Jacobian {
            x: x3,
            y: e * (d - x3) - eight_c,
            z: y1z1 + y1z1
        }
    }
}

impl<P: GroupParams> Mul<Fr> for Jacobian<P> {
    type Output = Jacobian<P>;

    fn mul(self, other: Fr) -> Jacobian<P> {
        let mut res = Jacobian::zero();
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

impl<P: GroupParams> Add<Jacobian<P>> for Jacobian<P> {
    type Output = Jacobian<P>;

    fn add(self, other: Jacobian<P>) -> Jacobian<P> {
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

            Jacobian {
                x: x3,
                y: r * (v - x3) - (s1_j + s1_j),
                z: ((self.z + other.z).squared() - z1_squared - z2_squared) * h
            }
        }
    }
}

impl<P: GroupParams> Neg for Jacobian<P> {
    type Output = Jacobian<P>;

    fn neg(self) -> Jacobian<P> {
        Jacobian {
            x: self.x,
            y: -self.y,
            z: self.z
        }
    }
}

impl<P: GroupParams> Sub<Jacobian<P>> for Jacobian<P> {
    type Output = Jacobian<P>;

    fn sub(self, other: Jacobian<P>) -> Jacobian<P> {
        self + (-other)
    }
}

pub struct G1Params;

impl GroupParams for G1Params {
    type Base = Fq;

    fn name() -> &'static str { "G1" }

    fn one() -> Jacobian<Self> {
        Jacobian {
            x: Fq::one(),
            y: const_fp([0xa6ba871b8b1e1b3a, 0x14f1d651eb8e167b, 0xccdd46def0f28c58, 0x1c14ef83340fbe5e]),
            z: Fq::one()
        }
    }
}

pub type G1 = Jacobian<G1Params>;

pub struct G2Params;

impl GroupParams for G2Params {
    type Base = Fq2;

    fn name() -> &'static str { "G2" }

    fn one() -> Jacobian<Self> {
        Jacobian {
            x: Fq2::new(const_fp([0x8e83b5d102bc2026, 0xdceb1935497b0172, 0xfbb8264797811adf, 0x19573841af96503b]), const_fp([0xafb4737da84c6140, 0x6043dd5a5802d8c4, 0x09e950fc52a02f86, 0x14fef0833aea7b6b])),
            y: Fq2::new(const_fp([0x619dfa9d886be9f6, 0xfe7fd297f59e9b78, 0xff9e1a62231b7dfe, 0x28fd7eebae9e4206]), const_fp([0x64095b56c71856ee, 0xdc57f922327d3cbb, 0x55f935be33351076, 0x0da4a0e693fd6482])),
            z: Fq2::one()
        }
    }
}

pub type G2 = Jacobian<G2Params>;

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
