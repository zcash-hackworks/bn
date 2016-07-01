pub mod fp;

#[cfg(test)]
pub mod tests;

use rand::Rng;
use self::fp::{Fp, PrimeFieldParams};

pub trait Field: Sized + Clone {
	fn zero() -> Self;
	fn one() -> Self;
    fn random<R: Rng>(rng: &mut R) -> Self;
    fn is_zero(&self) -> bool;
    fn inverse(&self) -> Self;
    fn squared(&self) -> Self {
    	self.mul(self)
    }
    fn pow<P: PrimeFieldParams>(&self, exp: &Fp<P>) -> Self;
    fn eq(&self, other: &Self) -> bool;
    fn ne(&self, other: &Self) -> bool {
    	!self.eq(other)
    }

    fn neg(&self) -> Self;
    fn mul(&self, other: &Self) -> Self;
    fn sub(&self, other: &Self) -> Self;
    fn add(&self, other: &Self) -> Self;
}
