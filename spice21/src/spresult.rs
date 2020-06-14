use std::error::Error;
use std::fmt;

#[derive(Debug)]
pub struct SpError {
    desc: String
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

pub type SpResult<T> = Result<T, &'static str>;
pub type TestResult = SpResult<()>;

