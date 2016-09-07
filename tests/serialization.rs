extern crate bincode;
extern crate rustc_serialize;

use self::bincode::SizeLimit::Infinite;
use self::bincode::rustc_serialize::{encode, decode};
use self::rustc_serialize::{Encodable, Decodable};
use self::rustc_serialize::hex::{FromHex, ToHex};

pub fn into_hex<S: Encodable>(obj: S) -> Option<String> {
    encode(&obj, Infinite).ok().map(|e| e.to_hex())
}

pub fn from_hex<S: Decodable>(s: &str) -> Option<S> {
    let s = s.from_hex().unwrap();

    decode(&s).ok()
}

pub fn reserialize<S: Encodable + Decodable>(obj: S) -> S {
    let s = into_hex(obj).unwrap();

    from_hex(&s).unwrap()
}
