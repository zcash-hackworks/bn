macro_rules! forward_val_val_binop {
    (impl($($t:ident: $p:ident),*) $imp:ident for $res:ty, $method:ident, $rhs:ty) => {
        impl<$($t: $p),*> $imp<$rhs> for $res {
            type Output = $res;

            #[inline]
            fn $method(self, other: $rhs) -> $res {
                $imp::$method(&self, &other)
            }
        }
    }
}

macro_rules! forward_ref_val_binop {
    (impl($($t:ident: $p:ident),*) $imp:ident for $res:ty, $method:ident, $rhs:ty) => {
        impl<'a, $($t: $p),*> $imp<$rhs> for &'a $res {
            type Output = $res;

            #[inline]
            fn $method(self, other: $rhs) -> $res {
                $imp::$method(self, &other)
            }
        }
    }
}

macro_rules! forward_val_ref_binop {
    (impl($($t:ident: $p:ident),*) $imp:ident for $res:ty, $method:ident, $rhs:ty) => {
        impl<'a, $($t: $p),*> $imp<&'a $rhs> for $res {
            type Output = $res;

            #[inline]
            fn $method(self, other: &$rhs) -> $res {
                $imp::$method(&self, other)
            }
        }
    }
}

macro_rules! forward_all_binop_to_ref_ref {
    (impl($($t:ident: $p:ident),*) $imp:ident for $res:ty, $method:ident, $rhs:ty) => {
        forward_val_val_binop!(impl($($t: $p),*) $imp for $res, $method, $rhs);
        forward_ref_val_binop!(impl($($t: $p),*) $imp for $res, $method, $rhs);
        forward_val_ref_binop!(impl($($t: $p),*) $imp for $res, $method, $rhs);
    };
}

macro_rules! forward_ops_to_group_ops {
    (impl($($t:ident: $p:ident),*) $res:ty) => {
        impl<'a, 'b, $($t: $p),*> Add<&'a $res> for &'b $res {
            type Output = $res;

            #[inline]
            fn add(self, other: &'a $res) -> $res {
                Jacobian::add(self, other)
            }
        }

        impl<'a, 'b, $($t: $p),*> Mul<&'a Fr> for &'b $res {
            type Output = $res;

            #[inline]
            fn mul(self, other: &'a Fr) -> $res {
                Jacobian::mul(self, other)
            }
        }

        impl<'a, 'b, $($t: $p),*> Sub<&'a $res> for &'b $res {
            type Output = $res;

            #[inline]
            fn sub(self, other: &'a $res) -> $res {
                Jacobian::sub(self, other)
            }
        }

        impl<'a, $($t: $p),*> Neg for &'a $res {
            type Output = $res;

            #[inline]
            fn neg(self) -> $res {
                Jacobian::neg(self)
            }
        }

        impl<$($t: $p),*> Neg for $res {
            type Output = $res;

            #[inline]
            fn neg(self) -> $res {
                Jacobian::neg(&self)
            }
        }

        impl<$($t: $p),*> PartialEq for $res {
            fn eq(&self, other: &Self) -> bool {
                Jacobian::eq(self, other)
            }
        }

        impl<$($t: $p),*> Eq for $res {}

        forward_all_binop_to_ref_ref!(impl($($t: $p),*) Add for $res, add, $res);
        forward_all_binop_to_ref_ref!(impl($($t: $p),*) Sub for $res, sub, $res);
        forward_all_binop_to_ref_ref!(impl($($t: $p),*) Mul for $res, mul, Fr);
    }
}
