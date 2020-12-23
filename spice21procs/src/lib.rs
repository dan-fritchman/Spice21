//!
//! # Spice21 Procedural Macros
//!
//! A very simple crate, broken out solely for sake of rustc's requirement that `proc-macro` crates
//! *only* export their procedural macros. (These function as something like compiler plug-ins.)
//!

use proc_macro;
use quote::quote;
use syn::{parse_macro_input, DeriveInput};

///
/// Derive-Macro for Trait `SpProto`
/// Implementations are all empty, and defer to the default implementation.
///
#[proc_macro_derive(SpProto)]
pub fn derive_sp(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let name = input.ident;
    let expanded = quote! {
        impl spice21int::SpProto for #name { }
    };
    proc_macro::TokenStream::from(expanded)
}
