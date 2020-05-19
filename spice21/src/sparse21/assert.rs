///
/// Assertion-Based Debugging Utilities
///

pub struct Assert<T> { val: T }

pub fn assert<T>(val: T) -> Assert<T> { Assert { val } }

impl<T> Assert<T> {
    fn raise(&self) -> Result<(), &'static str> {
        // Breakpoint here
        return Err("Assertion Failed");
    }
}

impl<T: PartialEq> Assert<T> {
    pub fn eq(&self, other: T) -> Result<(), &'static str> {
        if self.val != other { self.raise() } else { Ok(()) }
    }
    pub fn ne(&self, other: T) -> Result<(), &'static str> {
        if self.val == other { self.raise() } else { Ok(()) }
    }
}

impl<T: PartialOrd> Assert<T> {
    pub fn gt(&self, other: T) -> Result<(), &'static str> {
        if self.val <= other { self.raise() } else { Ok(()) }
    }
    pub fn lt(&self, other: T) -> Result<(), &'static str> {
        if self.val >= other { self.raise() } else { Ok(()) }
    }
    pub fn ge(&self, other: T) -> Result<(), &'static str> {
        if self.val < other { self.raise() } else { Ok(()) }
    }
    pub fn le(&self, other: T) -> Result<(), &'static str> {
        if self.val > other { self.raise() } else { Ok(()) }
    }
}
