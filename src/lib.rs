extern crate num;
extern crate rand;

mod fields;
mod params;

mod groups;

pub use fields::fp::Fp;
pub use fields::fp2::Fp2;
pub use params::{FrParams,FqParams,Fq2Params,G1Params,G2Params};
pub use fields::Field;
pub use groups::Jacobian;

pub type Fr = Fp<FrParams>;
pub type Fq = Fp<FqParams>;
pub type Fq2 = Fp2<Fq2Params>;

pub type Scalar = Fr;
pub type G1 = Jacobian<G1Params>;
pub type G2 = Jacobian<G2Params>;
