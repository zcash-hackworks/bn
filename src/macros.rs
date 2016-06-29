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
