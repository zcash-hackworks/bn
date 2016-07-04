use fields::Field;
use fields::fp::{PrimeFieldParams, Fp};
use params::G2Params;
use super::{Fr,Fq,Fq2};
use rand::Rng;
use std::ops::{Add,Mul,Sub,Neg};
use std::fmt;

#[cfg(test)]
pub mod tests;

#[macro_use]
mod macros;

mod gt;
pub use self::gt::Gt;

pub trait GroupParams: Sized {
    type Base: Field;

    fn zero() -> Jacobian<Self>;
    fn one() -> Jacobian<Self>;
    fn coeff_b() -> Self::Base;
    fn name() -> &'static str;
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
            x: self.x.clone(),
            y: self.y.clone(),
            z: self.z.clone()
        }
    }
}
pub enum Affine<P: GroupParams> {
    Zero,
    Point{x: P::Base, y: P::Base}
}

impl<P: GroupParams> Clone for Affine<P> {
    fn clone(&self) -> Self {
        match *self {
            Affine::Zero => Affine::Zero,
            Affine::Point{ref x, ref y} => Affine::Point{x: x.clone(), y: y.clone()}
        }
    }
}

#[derive(PartialEq, Eq)]
pub struct EllCoeffs {
    pub ell_0: Fq2,
    pub ell_vw: Fq2,
    pub ell_vv: Fq2
}

impl<P: GroupParams> Affine<P> {
    pub fn get_x(&self) -> P::Base {
        match *self {
            Affine::Zero => P::Base::zero(),
            Affine::Point{ref x, ref y} => x.clone()
        }
    }

    pub fn get_y(&self) -> P::Base {
        match *self {
            Affine::Zero => P::Base::one(),
            Affine::Point{ref x, ref y} => y.clone()
        }
    }

    pub fn to_jacobian(&self) -> Jacobian<P> {
        match *self {
            Affine::Zero => P::zero(),
            Affine::Point{ref x, ref y} => Jacobian {
                x: x.clone(),
                y: y.clone(),
                z: P::Base::one()
            }
        }
    }

    pub fn neg(&self) -> Self {
        match *self {
            Affine::Zero => Affine::Zero,
            Affine::Point{ref x, ref y} => Affine::Point{x: x.clone(), y: y.neg()}
        }
    }
}

impl Jacobian<G2Params> {
    pub fn mul_by_q(&self) -> Jacobian<G2Params> {
        Jacobian {
            x: G2Params::twist_mul_by_q_x() * self.x.frobenius_map(1),
            y: G2Params::twist_mul_by_q_y() * self.y.frobenius_map(1),
            z: self.z.frobenius_map(1)
        }
    }
}

impl Affine<G2Params> {
    pub fn mul_by_q(&self) -> Affine<G2Params> {
        self.to_jacobian().mul_by_q().to_affine()
    }
}

impl<P: GroupParams> Jacobian<P> {
    pub fn new(x: P::Base, y: P::Base, z: P::Base) -> Option<Self> {
        let tmp = Jacobian {
            x: x,
            y: y,
            z: z
        };

        if tmp.is_well_formed() {
            Some(tmp)
        } else {
            None
        }
    }

    pub fn random<R: Rng>(rng: &mut R) -> Self {
        Self::one() * &Fr::random(rng)
    }

    pub fn is_well_formed(&self) -> bool {
        if self.is_zero() {
            true
        } else {
            let x2 = self.x.squared();
            let y2 = self.y.squared();
            let z2 = self.z.squared();

            let x3 = self.x.mul(&x2);
            let z3 = self.z.mul(&z2);
            let z6 = z3.squared();

            (y2.eq(&z6.mul(&P::coeff_b()).add(&x3)))
        }
    }

    pub fn to_affine(&self) -> Affine<P> {
        if self.is_zero() {
            Affine::Zero
        } else {
            let z_inv = self.z.inverse();
            let z2_inv = z_inv.squared();
            let z3_inv = z2_inv.mul(&z_inv);

            Affine::Point {
                x: self.x.mul(&z2_inv),
                y: self.y.mul(&z3_inv)
            }
        }
    }

    pub fn zero() -> Self {
        P::zero()
    }

    pub fn one() -> Self {
        P::one()
    }

    pub fn add(&self, other: &Self) -> Self {
        if self.is_zero() {
            return other.clone()
        }

        if other.is_zero() {
            return self.clone()
        }

        let z1_squared = self.z.squared();
        let z2_squared = other.z.squared();
        let u1 = self.x.mul(&z2_squared);
        let u2 = other.x.mul(&z1_squared);
        let z1_cubed = self.z.mul(&z1_squared);
        let z2_cubed = other.z.mul(&z2_squared);
        let s1 = self.y.mul(&z2_cubed);
        let s2 = other.y.mul(&z1_cubed);

        if u1.eq(&u2) && s1.eq(&s2) {
            self.double()
        } else {
            let h = u2.sub(&u1);
            let s2_minus_s1 = s2.sub(&s1);
            let i = h.add(&h).squared();
            let j = h.mul(&i);
            let r = s2_minus_s1.add(&s2_minus_s1);
            let v = u1.mul(&i);
            let s1_j = s1.mul(&j);
            let x3 = r.squared().sub(&j).sub(&v.add(&v));
            let y3 = r.mul(&v.sub(&x3)).sub(&s1_j.add(&s1_j));

            Jacobian {
                x: x3,
                y: y3,
                z: self.z.add(&other.z).squared().sub(&z1_squared).sub(&z2_squared).mul(&h)
            }
        }
    }

    pub fn double(&self) -> Self {
        let a = self.x.squared();
        let b = self.y.squared();
        let c = b.squared();
        let mut d = self.x.add(&b).squared().sub(&a).sub(&c);
        d = d.add(&d);
        let e = a.add(&a).add(&a);
        let f = e.squared();
        let x3 = f.sub(&d.add(&d));
        let mut eight_c = c.add(&c);
        eight_c = eight_c.add(&eight_c);
        eight_c = eight_c.add(&eight_c);
        let y3 = e.mul(&d.sub(&x3)).sub(&eight_c);
        let y1z1 = self.y.mul(&self.z);
        let z3 = y1z1.add(&y1z1);

        Jacobian {
            x: x3,
            y: y3,
            z: z3
        }
    }

    pub fn eq(&self, other: &Self) -> bool {
        if self.is_zero() {
            return other.is_zero()
        }

        if other.is_zero() {
            return false;
        }

        let z1_squared = self.z.squared();
        let z2_squared = other.z.squared();

        if self.x.mul(&z2_squared).ne(&other.x.mul(&z1_squared)) {
            return false;
        }

        let z1_cubed = self.z.mul(&z1_squared);
        let z2_cubed = other.z.mul(&z2_squared);

        if self.y.mul(&z2_cubed).ne(&other.y.mul(&z1_cubed)) {
            return false;
        }

        return true;
    }

    pub fn neg(&self) -> Self {
        Jacobian {
            x: self.x.clone(),
            y: self.y.neg(),
            z: self.z.clone()
        }
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.z.is_zero()
    }

    pub fn mul<S: PrimeFieldParams>(&self, other: &Fp<S>) -> Jacobian<P> {
        let mut result = Jacobian::<P>::zero();
        let mut found_one = false;
        for i in (0..S::bits()).rev() {
            if found_one {
                result = result.double();
            }

            if other.test_bit(i) {
                found_one = true;
                result = &result + self;
            }
        }

        result
    }

    #[inline]
    pub fn sub(&self, other: &Self) -> Jacobian<P> {
        self.add(&other.neg())
    }
}

impl Jacobian<G2Params> {
    pub fn mixed_addition_step_for_flipped_miller_loop(&mut self, base: &Affine<G2Params>) -> EllCoeffs {
        let d = &self.x - &self.z * &base.get_x();
        let e = &self.y - &self.z * &base.get_y();
        let f = d.squared();
        let g = e.squared();
        let h = &d * &f;
        let i = &self.x * &f;
        let j = &self.z * &g + &h - (&i + &i);

        self.x = &d * &j;
        self.y = &e * (&i - &j) - &h * &self.y;
        self.z = &self.z * &h;

        EllCoeffs {
            ell_0: G2Params::twist() * (&e * &base.get_x() - &d * &base.get_y()),
            ell_vv: e.neg(),
            ell_vw: d
        }
    }

    pub fn doubling_step_for_flipped_miller_loop(&mut self, two_inv: &Fq) -> EllCoeffs {
        let a = &(&self.x * &self.y) * two_inv;
        let b = self.y.squared();
        let c = self.z.squared();
        let d = &c + &c + &c;
        let e = G2Params::coeff_b() * &d;
        let f = &e + &e + &e;
        let g = &(&b + &f) * two_inv;
        let h = (&self.y + &self.z).squared() - (&b + &c);
        let i = &e - &b;
        let j = self.x.squared();
        let e_sq = e.squared();

        self.x = &a * (&b - &f);
        self.y = g.squared() - (&e_sq + &e_sq + &e_sq);
        self.z = &b * &h;

        EllCoeffs {
            ell_0: G2Params::twist() * &i,
            ell_vw: h.neg(),
            ell_vv: &j + &j + &j
        }
    }
}

forward_ops_to_group_ops!(impl(P: GroupParams) Jacobian<P>);

impl<P: GroupParams> PartialEq for Affine<P> {
    fn eq(&self, other: &Self) -> bool {
        self.to_jacobian() == other.to_jacobian()
    }
}

impl<P: GroupParams> Eq for Affine<P> { }
