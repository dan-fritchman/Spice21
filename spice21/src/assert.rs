///
/// Assertion-Based Debugging Utilities
///
use super::spresult::{TestResult, sperror};

pub struct Assert<T> {
    val: T,
}

pub fn assert<T>(val: T) -> Assert<T> {
    Assert { val }
}

impl<T> Assert<T> {
    fn raise(&self) -> TestResult {
        // Breakpoint here
        return Err(sperror("Assertion Failed"));
    }
}

impl<T: PartialEq> Assert<T> {
    pub fn eq(&self, other: T) -> TestResult {
        if self.val != other {
            self.raise()
        } else {
            Ok(())
        }
    }
    pub fn ne(&self, other: T) -> TestResult {
        if self.val == other {
            self.raise()
        } else {
            Ok(())
        }
    }
}

impl<T: PartialOrd> Assert<T> {
    pub fn gt(&self, other: T) -> TestResult {
        if self.val <= other {
            self.raise()
        } else {
            Ok(())
        }
    }
    pub fn lt(&self, other: T) -> TestResult {
        if self.val >= other {
            self.raise()
        } else {
            Ok(())
        }
    }
    pub fn ge(&self, other: T) -> TestResult {
        if self.val < other {
            self.raise()
        } else {
            Ok(())
        }
    }
    pub fn le(&self, other: T) -> TestResult {
        if self.val > other {
            self.raise()
        } else {
            Ok(())
        }
    }
}

impl Assert<f64> {
    pub fn abs(self) -> Assert<f64> {
        Assert {
            val: self.val.abs(),
        }
    }
    pub fn isclose(&self, other: f64, tol: f64) -> TestResult {
        if (self.val - other).abs() > tol {
            self.raise()
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
