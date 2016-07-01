macro_rules! forward_val_val_binop {
    (impl($($t:ident: $p:ident),*) $imp:ident for $res:ty, $method:ident) => {
        impl<$($t: $p),*> $imp<$res> for $res {
            type Output = $res;

            #[inline]
            fn $method(self, other: $res) -> $res {
                $imp::$method(&self, &other)
            }
        }
    }
}

macro_rules! forward_ref_val_binop {
    (impl($($t:ident: $p:ident),*) $imp:ident for $res:ty, $method:ident) => {
        impl<'a, $($t: $p),*> $imp<$res> for &'a $res {
            type Output = $res;

            #[inline]
            fn $method(self, other: $res) -> $res {
                $imp::$method(self, &other)
            }
        }
    }
}

macro_rules! forward_val_ref_binop {
    (impl($($t:ident: $p:ident),*) $imp:ident for $res:ty, $method:ident) => {
        impl<'a, $($t: $p),*> $imp<&'a $res> for $res {
            type Output = $res;

            #[inline]
            fn $method(self, other: &$res) -> $res {
                $imp::$method(&self, other)
            }
        }
    }
}

macro_rules! forward_all_binop_to_ref_ref {
    (impl($($t:ident: $p:ident),*) $imp:ident for $res:ty, $method:ident) => {
        forward_val_val_binop!(impl($($t: $p),*) $imp for $res, $method);
        forward_ref_val_binop!(impl($($t: $p),*) $imp for $res, $method);
        forward_val_ref_binop!(impl($($t: $p),*) $imp for $res, $method);
    };
}

macro_rules! forward_ops_to_field_ops {
    (impl($($t:ident: $p:ident),*) $res:ty) => {
        impl<'a, 'b, $($t: $p),*> Add<&'a $res> for &'b $res {
            type Output = $res;

            #[inline]
            fn add(self, other: &'a $res) -> $res {
                Field::add(self, other)
            }
        }

        impl<'a, 'b, $($t: $p),*> Sub<&'a $res> for &'b $res {
            type Output = $res;

            #[inline]
            fn sub(self, other: &'a $res) -> $res {
                Field::sub(self, other)
            }
        }

        impl<'a, 'b, $($t: $p),*> Mul<&'a $res> for &'b $res {
            type Output = $res;

            #[inline]
            fn mul(self, other: &'a $res) -> $res {
                Field::mul(self, other)
            }
        }

        impl<'a, $($t: $p),*> Neg for &'a $res {
            type Output = $res;

            #[inline]
            fn neg(self) -> $res {
                Field::neg(self)
            }
        }

        impl<$($t: $p),*> Neg for $res {
            type Output = $res;

            #[inline]
            fn neg(self) -> $res {
                Field::neg(&self)
            }
        }

        impl<$($t: $p),*> PartialEq for $res {
            fn eq(&self, other: &Self) -> bool {
                Field::eq(self, other)
            }
        }

        impl<$($t: $p),*> Eq for $res {}

        forward_all_binop_to_ref_ref!(impl($($t: $p),*) Add for $res, add);
        forward_all_binop_to_ref_ref!(impl($($t: $p),*) Sub for $res, sub);
        forward_all_binop_to_ref_ref!(impl($($t: $p),*) Mul for $res, mul);
    }
}
