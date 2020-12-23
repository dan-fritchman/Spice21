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
/// is the best way we've found to import it to the rest of Spice21.
/// Note: this must be defined *before* any uses of it.
///

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
        #[allow(dead_code)]
        #[doc=$struct_desc]
        #[derive(Clone)]
        pub struct $src_name {
            $( #[doc=$desc]
                pub $attr_name : $attr_type ),*
        }
        impl $src_name {
            #[allow(dead_code)]
            fn getattr<S: Into<String>>(&self, key: S) -> Option<f64> {
                let k: String = key.into();
                match &k as &str {
                    $( stringify!($attr_name) => Some(self.$attr_name)),*,
                    _ => None,
                }
            }
        }
    }
    }

    /// The `from_opt_type` derives a type of the form `struct { x: X }` from `struct { x: Option<X> }`,
    /// given an existing type `from_type`, and default values for each field.
    #[macro_export]
    macro_rules! from_opt_type {
    ( $dest_type:ident, $from_type:ident, $struct_desc:literal, [
        $( ($attr_name:ident, $attr_type:ty, $default:literal, $desc:literal) ),* $(,)?
    ]) => {
        #[allow(dead_code)]
        #[doc=$struct_desc]
        #[derive(Clone)]
        pub struct $dest_type {
            $( #[doc=$desc]
                pub $attr_name : $attr_type ),*
        }
        impl From<$from_type> for $dest_type {
            pub fn from(specs: $from_type) -> Self {
                Self {
                    $($attr_name : if let Some(val) = specs.$attr_name { val } else { $default }; ),*,
                }
            }
        }
        impl Default for $dest_type {
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
        #[derive(Clone, Copy, Default, Serialize, Deserialize, Debug)]
        pub struct $specs_name {
            $(  #[doc=$desc]
                #[serde(default)]
                pub $attr_name : Option<$attr_type> ),*
        }
        #[doc=$struct_desc]
        #[derive(Clone, Copy, Default, Serialize, Deserialize, Debug)]
        pub struct $vals_name {
            $( pub $attr_name : $attr_type, )*
            $( pub $val_name : $val_type, )*
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
pub mod defs;
pub mod elab;
pub mod proto;
pub mod sparse21;
pub mod spresult;

// Re-exports
pub use analysis::*;
pub use proto::*;
pub use spresult::*;

// Crate-wide public
pub(crate) use spnum::*;

// Private modules
mod assert;
mod spnum;

#[cfg(test)]
mod tests;

// Use our internal friend-crate(s)
// Note macro-import *must* be here at crate-root, *not* down-hierarchy 
#[macro_use]
extern crate spice21procs;
