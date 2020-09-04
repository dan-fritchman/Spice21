use num::traits::NumAssignOps;
use num::{Complex, Num, One, Zero};
use std::fmt;

// This long list of traits describes our required behavior for numeric types.
pub trait SpNum:
    Clone + Copy + NumAssignOps + Zero + Num + Abs + fmt::Display + fmt::Debug
{
}

impl<T> SpNum for T where
    T: Clone + Copy + NumAssignOps + Zero + One + Num + Abs + fmt::Display + fmt::Debug
{
}

/// Absolute Value Trait for numeric types
pub trait Abs {
    fn absv(&self) -> f64;
}

impl Abs for f64 {
    fn absv(&self) -> f64 {
        self.abs()
    }
}

impl Abs for Complex<f64> {
    fn absv(&self) -> f64 {
        self.norm()
    }
}
