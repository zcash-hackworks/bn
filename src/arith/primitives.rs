pub type Digit = u32;
pub const LIMBS: usize = 8;
type DoubleDigit = u64;
const BITS: usize = 32;
const BASE: DoubleDigit = 1 << BITS;
const LO_MASK: DoubleDigit = (-1i32 as DoubleDigit) >> BITS;

#[inline]
fn get_hi(n: DoubleDigit) -> Digit {
    (n >> BITS) as Digit
}

#[inline]
fn get_lo(n: DoubleDigit) -> Digit {
    (n & LO_MASK) as Digit
}

#[inline]
fn split(n: DoubleDigit) -> (Digit, Digit) {
    (get_hi(n), get_lo(n))
}

/// Divide by two
#[inline]
pub fn div2(a: &mut [Digit; LIMBS]) {
    let mut b = 0;
    for a in a.iter_mut().rev() {
        let t = *a << 31;
        *a = *a >> 1;
        *a = *a | b;
        b = t;
    }
}

#[inline]
pub fn add_nocarry(a: &mut [Digit; LIMBS], b: &[Digit; LIMBS]) {
    #[inline]
    fn adc(a: Digit, b: Digit, carry: &mut Digit) -> Digit {
        let (hi, lo) = split((a as DoubleDigit) + (b as DoubleDigit) +
                             (*carry as DoubleDigit));

        *carry = hi;
        lo
    }

    let mut carry = 0;

    for (a, b) in a.into_iter().zip(b.iter()) {
        *a = adc(*a, *b, &mut carry);
    }

    assert!(0 == carry);
}

#[inline]
pub fn sub_noborrow(a: &mut [Digit; LIMBS], b: &[Digit; LIMBS]) {
    #[inline]
    fn sbb(a: Digit, b: Digit, borrow: &mut Digit) -> Digit {
        let (hi, lo) = split(BASE + (a as DoubleDigit) -
                             (b as DoubleDigit) -
                             (*borrow as DoubleDigit));
        *borrow = (hi == 0) as Digit;
        lo
    }

    let mut borrow = 0;

    for (a, b) in a.into_iter().zip(b.iter()) {
        *a = sbb(*a, *b, &mut borrow);
    }

    assert!(0 == borrow);
}


pub fn mul_reduce(
    this: &mut [Digit; LIMBS],
    by: &[Digit; LIMBS],
    modulus: &[Digit; LIMBS],
    inv: Digit
)
{
    #[inline]
    fn mac_digit(acc: &mut [Digit], b: &[Digit], c: Digit)
    {
        #[inline]
        fn mac_with_carry(a: Digit, b: Digit, c: Digit, carry: &mut Digit) -> Digit {
            let (hi, lo) = split((a as DoubleDigit) +
                                 (b as DoubleDigit) * (c as DoubleDigit) +
                                 (*carry as DoubleDigit));
            *carry = hi;
            lo
        }

        if c == 0 {
            return;
        }

        let mut b_iter = b.iter();
        let mut carry = 0;

        for ai in acc.iter_mut() {
            if let Some(bi) = b_iter.next() {
                *ai = mac_with_carry(*ai, *bi, c, &mut carry);
            } else if carry != 0 {
                *ai = mac_with_carry(*ai, 0, c, &mut carry);
            } else {
                break;
            }
        }

        assert!(carry == 0);
    }

    // The Montgomery reduction here is based on Algorithm 14.32 in
    // Handbook of Applied Cryptography
    // <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.

    let mut res = [0; 2*LIMBS];
    for (i, xi) in this.iter().enumerate() {
        mac_digit(&mut res[i..], by, *xi);
    }

    for i in 0..LIMBS {
        let k = inv.wrapping_mul(res[i]);
        mac_digit(&mut res[i..], modulus, k);
    }

    this.copy_from_slice(&res[LIMBS..]);
}
