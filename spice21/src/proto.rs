/// Export everything from out protos-crate 
pub use spice21protos::*;

use crate::SpError;

/// Error Conversion
impl From<prost::DecodeError> for SpError {
    fn from(_e: prost::DecodeError) -> Self {
        SpError::new("Error Decoding Circuit Binary")
    }
}
