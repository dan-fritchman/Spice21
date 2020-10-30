//!
//! Assertion-Based Debugging Utilities
//!
use super::spresult::{sperror, TestResult};
use std::fmt::Debug;

pub struct Assert<T> {
    val: T,
}

pub fn assert<T>(val: T) -> Assert<T> {
    Assert { val }
}

impl<T> Assert<T> {
    fn raise<S: Into<String>>(&self, s: S) -> TestResult {
        return Err(sperror(s.into())); // Breakpoint here catches all failures
    }
    /// Identity function
    /// For "grammatical" usage, such as
    /// `assert(3).is().gt(0);`
    pub fn is(self) -> Self {
        self
    }
}
impl<T: PartialEq + Debug> Assert<T> {
    pub fn eq(&self, other: T) -> TestResult {
        if self.val != other {
            self.raise(format!("Assert Eq Failed: {:?} != {:?}", self.val, other))
        } else {
            Ok(())
        }
    }
    pub fn ne(&self, other: T) -> TestResult {
        if self.val == other {
            self.raise(format!("Assert Neq Failed: {:?} == {:?}", self.val, other))
        } else {
            Ok(())
        }
    }
}
impl<T: PartialOrd + Debug> Assert<T> {
    pub fn gt(&self, other: T) -> TestResult {
        if self.val <= other {
            self.raise(format!("Assert Gt Failed: {:?} <= {:?}", self.val, other))
        } else {
            Ok(())
        }
    }
    pub fn lt(&self, other: T) -> TestResult {
        if self.val >= other {
            self.raise(format!("Assert Lt Failed: {:?} >= {:?}", self.val, other))
        } else {
            Ok(())
        }
    }
    pub fn ge(&self, other: T) -> TestResult {
        if self.val < other {
            self.raise(format!("Assert Geq Failed: {:?} < {:?}", self.val, other))
        } else {
            Ok(())
        }
    }
    pub fn le(&self, other: T) -> TestResult {
        if self.val > other {
            self.raise(format!("Assert Leq Failed: {:?} > {:?}", self.val, other))
        } else {
            Ok(())
        }
    }
}
impl Assert<f64> {
    pub fn abs(self) -> Assert<f64> {
        Assert { val: self.val.abs() }
    }
    pub fn isclose(&self, other: f64, tol: f64) -> TestResult {
        if (self.val - other).abs() > tol {
            self.raise(format!("Assert IsClose Failed: abs({:?} - {:?}) > {:?}", self.val, other, tol))
        } else {
            Ok(())
        }
    }
}
impl Assert<&Vec<f64>> {
    /// Last value in vector 
    pub fn last(&self) -> Assert<f64> {
        assert(self.val[self.val.len()-1])
    }
    /// Tests that a vector is strictly increasing,
    /// i.e. that each element k is > element k-1.
    pub fn increasing(&self) -> TestResult {
        for k in 0..self.val.len() - 1 {
            if self.val[k + 1] <= self.val[k] {
                return self.raise(format!(
                    "Non-increasing vector values {:?} and {:?} at index {}",
                    self.val[k],
                    self.val[k + 1],
                    k
                ));
            }
        }
        Ok(())
    }
    pub fn nondecreasing(&self) -> TestResult {
        for k in 0..self.val.len() - 1 {
            if self.val[k + 1] < self.val[k] {
                return self.raise(format!(
                    "Decreasing vector values {:?} and {:?} at index {}",
                    self.val[k],
                    self.val[k + 1],
                    k
                ));
            }
        }
        Ok(())
    }
    pub fn decreasing(&self) -> TestResult {
        for k in 0..self.val.len() - 1 {
            if self.val[k + 1] >= self.val[k] {
                return self.raise(format!(
                    "Non-decreasing vector values {:?} and {:?} at index {}",
                    self.val[k],
                    self.val[k + 1],
                    k
                ));
            }
        }
        Ok(())
    }
    /// Test for constant-valued vectors
    pub fn constant(&self, val: f64) -> TestResult {
        for k in 0..self.val.len() {
            if self.val[k] != val {
                return self.raise(format!("Non-constant vector with value {} at index {}", self.val[k], k));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eq() -> TestResult {
        assert(5).eq(5)
    }

    #[test]
    fn test_ne() -> TestResult {
        assert(1).ne(5)
    }

    #[test]
    fn test_gt() -> TestResult {
        assert(5).gt(1)
    }

    #[test]
    fn test_lt() -> TestResult {
        assert(1).lt(5)
    }

    #[test]
    fn test_ge() -> TestResult {
        assert(5).ge(5)
    }

    #[test]
    fn test_le() -> TestResult {
        assert(1).le(5)
    }
}
