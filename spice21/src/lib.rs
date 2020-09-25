//!
//! # Spice21
//!
//! SPICE for the 21st Century
//!
//! [https://github.com/HW21/Spice21](https://github.com/HW21/Spice21)
//!

///
/// # Spice21 Macros
///
/// Note: this module's unusual location, inline in lib.rs,
/// is the best way we've found to import it
/// to the rest of Spice21.
/// Note: this must be defined *before* any uses of it.
///


use serde::{Serialize, Deserialize};

#[macro_use]
pub(crate) mod macros {
    /// GetAttr-enabled struct builder.
    /// Creates structs from a list of field-definitions,
    /// adding a `getattr` method enabling by-string access.
    #[macro_export]
    macro_rules! attr {
    ( $src_name:ident, $struct_desc:literal, [
        $( ($attr_name:ident, $attr_type:ty, $default:literal, $desc:literal) ),* $(,)?
    ]) => {
        #[doc=$struct_desc]
        #[derive(Clone, Copy)]
        pub struct $src_name {
            $( #[doc=$desc]
                pub $attr_name : $attr_type ),*
        }
        impl $src_name {
            fn getattr<S: Into<String>>(&self, key: S) -> Option<f64> {
                let k: String = key.into();
                match &k as &str {
                    $( stringify!($attr_name) => Some(self.$attr_name)),*,
                    _ => None,
                }
            }
        }
        impl Default for $src_name {
            fn default() -> Self {
                Self {
                    $($attr_name : $default),*,
                }
            }
        }
    }
    }

    #[macro_export]
    macro_rules! paramstruct {
    ( $src_name:ident, $struct_desc:literal, [
        $( ($attr_name:ident, $attr_type:ty, $desc:literal) ),* $(,)?
    ]) => {
        #[doc=$struct_desc]
        #[derive(Clone, Copy, Default)]
        pub struct $src_name {
            $( #[doc=$desc]
                pub $attr_name : $attr_type ),*
        }
    }
    }

    #[macro_export]
    macro_rules! specgen {
    ( $specs_name:ident, $vals_name: ident, $struct_desc:literal, [
        $( ($attr_name:ident, $attr_type:ty, $desc:literal) ),* $(,)?
    ], [
        $( ($val_name:ident, $val_type:ty) ),* $(,)?
    ]$(,)? ) => {
        #[doc=$struct_desc]
        #[derive(Clone, Copy, Default)]
        pub(crate) struct $specs_name {
            $( #[doc=$desc]
                pub(crate) $attr_name : Option<$attr_type> ),*
        }
        #[doc=$struct_desc]
        #[derive(Clone, Copy, Default, Serialize, Deserialize, Debug)]
        pub(crate) struct $vals_name {
            $( pub(crate) $attr_name : $attr_type, )*
            $( pub(crate) $val_name : $val_type, )*
        }
    }
    }

    #[cfg(test)]
    mod tests {
        use crate::assert::*;
        use crate::spresult::TestResult;

        attr!(
            SampleModel,
            "Doc-string for Sample Model",
            [
                (param1, f64, 1.1, "First Param"),
                (param22, f64, 2.2, "Parameter #2"),
                (par33, f64, 3.3, "3param"),
            ]
        );

        #[test]
        fn test1() -> TestResult {
            // Instantiate a generated struct
            let s = SampleModel {
                param1: 1.1,
                param22: 22.22,
                par33: 333.333,
            };

            // Test fields
            assert(s.param1).eq(1.1)?;
            assert(s.param22).eq(22.22)?;
            assert(s.par33).eq(333.333)?;

            // Test getattr
            assert(s.getattr("param1")).eq(Some(s.param1))?;
            assert(s.getattr("param22")).eq(Some(s.param22))?;
            assert(s.getattr("par33")).eq(Some(s.par33))?;
            assert(s.getattr("fizzbuzz")).eq(None)?;

            Ok(())
        }
    }
}

// Modules
pub mod analysis;
pub mod circuit;
pub mod comps;
pub mod proto;

// Re-exports
pub use analysis::*;
pub use proto::*;

// Crate-wide public
pub(crate) use spnum::*;
pub(crate) use spresult::*;

// Private modules
mod assert;
mod sparse21;
mod spnum;
mod spresult;
mod tests;
