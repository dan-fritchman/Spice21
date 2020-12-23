//!
//! # Spice21 Internal Crate
//!
//! A very simple crate, broken out solely for sake of rustc's requirement that `proc-macro` crates
//! *only* export their procedural macros. (These function as something like compiler plug-ins.)
//!
//! This crate exposes the internals that `spiceproc` (the procedural-macros crate) depends on.
//! It is not designed to be particularly useful outside of Spice21.
//! All `pub` annotations should be read as "public for Spice21".
//!

use prost::Message;
use std::io::Cursor;

///
/// # SpProto Custom Message Trait
///
/// Largely an extension of `prost::Message`,
/// derived for each Spice21 `Message` type.
///
pub trait SpProto: Message + Sized + Default {
    /// Encode into Byte-Vector
    fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::<u8>::with_capacity(self.encoded_len());
        buf.reserve(self.encoded_len());
        self.encode(&mut buf).unwrap();
        buf
    }
    /// Decode from byte array/vector 
    fn from_bytes(bytes: &[u8]) -> Result<Self, prost::DecodeError> {
        Self::decode(&mut Cursor::new(bytes))
    }
}
