//!
//! Spice21 Result and Error Types
//!
use std::error::Error;
use std::fmt;

/// # Spice21 General Error Type 
#[derive(Debug)]
pub struct SpError {
    pub desc: String,
}
// Allow SpError in `dyn Error` contexts
impl Error for SpError {}
impl SpError {
    /// Spice Error Constructor, from anything String-convertible
    pub(crate) fn new<S: Into<String>>(s: S) -> SpError {
        SpError { desc: s.into() }
    }
    /// Create a Box'ed SpError 
    pub(crate) fn boxed<S: Into<String>>(s: S) -> Box<SpError> {
        Box::new(SpError { desc: s.into() })
    }
}
pub(crate) fn sperror<S: Into<String>>(s: S) -> SpError {
    SpError::new(s)
}

impl fmt::Display for SpError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.desc)
    }
}

// SpError Conversions 
impl From<&str> for SpError {
    fn from(s: &str) -> Self {
        SpError::new(s)
    }
}
impl From<& str> for Box<SpError> {
    fn from(s: &str) -> Self {
        SpError::boxed(s)
    }
}

/// # Spice21 General Result Type
pub type SpResult<T> = Result<T, SpError>;
/// # Spice21 Test Result Type 
pub type TestResult = SpResult<()>;
