//!
//! Assertion-Based Debugging Utilities
//!

#![allow(unused)] // Much of this is used here-and-there by tests

use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;

use super::spresult::{sperror, TestResult};

///
/// # Assertion Struct
///
/// Primarily used for semantic assertions, such as:
/// `assert(5).gt(3)?;`
/// `assert(vec).is().increasing()?;`
///
pub struct Assert<T> {
    val: T,
}
pub fn assert<T>(val: T) -> Assert<T> {
    Assert { val }
}
fn raise<S: Into<String>>(s: S) -> TestResult {
    return Err(sperror(s.into())); // Breakpoint here catches all failures
}

impl<T: PartialEq + Debug> Assert<T> {
    pub fn eq(&self, other: T) -> TestResult {
        if self.val != other {
            raise(format!("Assert Eq Failed: {:?} != {:?}", self.val, other))
        } else {
            Ok(())
        }
    }
    pub fn ne(&self, other: T) -> TestResult {
        if self.val == other {
            raise(format!("Assert Neq Failed: {:?} == {:?}", self.val, other))
        } else {
            Ok(())
        }
    }
}
impl<T: PartialOrd + Debug> Assert<T> {
    pub fn gt(&self, other: T) -> TestResult {
        if self.val <= other {
            raise(format!("Assert Gt Failed: {:?} <= {:?}", self.val, other))
        } else {
            Ok(())
        }
    }
}

impl<T> Assert<T> {
    /// Identity function
    /// For "grammatical" usage, such as
    /// `assert(3).is().gt(0);`
    pub fn is(self) -> Self {
        self
    }
}
impl Assert<&Vec<f64>> {
    /// Create an Assert regarding the last value in vector
    /// e.g. assert(vec![1,3,5]).last().gt(4)?;
    pub fn last(&self) -> Assert<f64> {
        assert(self.val[self.val.len() - 1])
    }
    /// Tests that a vector is strictly increasing,
    /// i.e. that each element k is > element k-1.
    pub fn increasing(&self) -> TestResult {
        for k in 0..self.val.len() - 1 {
            if self.val[k + 1] <= self.val[k] {
                return raise(format!(
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
                return raise(format!(
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
                return raise(format!(
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
                return raise(format!("Non-constant vector with value {} at index {}", self.val[k], k));
            }
        }
        Ok(())
    }
}
fn keys_match<T: Eq + Hash, U, V>(map1: &HashMap<T, U>, map2: &HashMap<T, V>) -> bool {
    map1.len() == map2.len() && map1.keys().all(|k| map2.contains_key(k))
}
impl Assert<&HashMap<String, Vec<f64>>> {
    pub fn isclose(&self, other: HashMap<String, Vec<f64>>, tol: f64) -> TestResult {
        if !keys_match(self.val, &other) {
            return raise(format!("HashMap keys do not match: {:?} vs {:?}", self.val.keys(), other.keys()));
        }
        for key in self.val.keys() {
            let s = self.val.get(key).unwrap();
            let o = other.get(key).unwrap();
            if s.len() != o.len() {
                return raise(format!("Mismatched Vector Lengths {} vs {} for Signal {}", s.len(), o.len(), key));
            }
            for i in 0..s.len() {
                if (s[i] - o[i]).abs() > tol {
                    return raise(format!("Assert IsClose Failed: abs({:?} - {:?}) > {:?}", s[i], o[i], tol));
                }
            }
        }
        Ok(())
    }
}
impl<T: PartialOrd + Debug> Assert<T> {
    pub fn lt(&self, other: T) -> TestResult {
        if self.val >= other {
            raise(format!("Assert Lt Failed: {:?} >= {:?}", self.val, other))
        } else {
            Ok(())
        }
    }
    pub fn ge(&self, other: T) -> TestResult {
        if self.val < other {
            raise(format!("Assert Geq Failed: {:?} < {:?}", self.val, other))
        } else {
            Ok(())
        }
    }
    pub fn le(&self, other: T) -> TestResult {
        if self.val > other {
            raise(format!("Assert Leq Failed: {:?} > {:?}", self.val, other))
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
            raise(format!("Assert IsClose Failed: abs({:?} - {:?}) > {:?}", self.val, other, tol))
        } else {
            Ok(())
        }
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
