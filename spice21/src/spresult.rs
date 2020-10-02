//! 
//! Spice21 Result and Error Types
//! 
use std::error::Error;
use std::fmt;

#[derive(Debug)]
pub struct SpError {
    pub desc: String,
}

pub(crate) fn sperror<S: Into<String>>(s: S) -> SpError {
    SpError { desc: s.into() }
}

impl Error for SpError {}

impl SpError {
    pub fn throw(s: &'static str) -> Box<SpError> {
        return Box::new(SpError { desc: String::from(s) });
    }
}

impl fmt::Display for SpError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.desc)
    }
}

impl From<&'static str> for SpError {
    fn from(s: &'static str) -> Self {
        SpError { desc: String::from(s) }
    }
}

impl From<&'static str> for Box<SpError> {
    fn from(s: &'static str) -> Self {
        Box::new(SpError { desc: String::from(s) })
    }
}

pub type SpResult<T> = Result<T, SpError>;
pub type TestResult = SpResult<()>;
