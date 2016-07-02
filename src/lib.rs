extern crate num;
extern crate rand;

mod fields;
mod params;

mod groups;

pub use fields::fp::Fp;
pub use params::{FrParams,FqParams,G1Params};
pub use fields::Field;
pub use groups::Jacobian;

pub type Fr = Fp<FrParams>;
pub type Fq = Fp<FqParams>;

pub type Scalar = Fr;
pub type G1 = Jacobian<G1Params>;
