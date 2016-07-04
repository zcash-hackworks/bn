use ::Fq12;
use ::Scalar;
use std::ops::{BitXor,Mul};
use fields::Field;
use std::cmp::{PartialEq, Eq};

#[derive(Debug,Eq,PartialEq)]
pub struct Gt(Fq12);

impl Gt {
    pub fn new(a: Fq12) -> Gt {
        return Gt(a);
    }
}

impl<'a, 'b> Mul<&'a Gt> for &'b Gt {
    type Output = Gt;

    fn mul(self, other: &Gt) -> Gt {
        Gt(&self.0 * &other.0)
    }
}

impl<'a, 'b> BitXor<&'a Scalar> for &'b Gt {
    type Output = Gt;

    fn bitxor(self, other: &Scalar) -> Gt {
        Gt(self.0.pow(other))
    }
}

forward_all_binop_to_ref_ref!(impl() Mul for Gt, mul, Gt);
forward_all_binop_to_ref_ref!(impl() BitXor for Gt, bitxor, Scalar);
