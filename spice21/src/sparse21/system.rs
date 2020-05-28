/// Sparse Matrix System
///
/// Represents a linear system of the form `Ax=b`
///


use std::error::Error;

#[derive(Debug, Clone)]
struct NonRealNumError;

impl fmt::Display for NonRealNumError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid first item to double")
    }
}

impl Error for NonRealNumError {
    fn description(&self) -> &str {
        "invalid first item to double"
    }

    fn cause(&self) -> Option<&Error> {
        // Generic error, underlying cause isn't tracked.
        None
    }
}

pub struct System {
    mat: Matrix,
    rhs: Vec<f64>,
    title: Option<String>,
    size: usize,
}

use std::path::Path;

impl System {
    /// Splits a `System` into a two-tuple of `self.matrix` and `self.rhs`.
    /// Nothing is copied; `self` is consumed in the process.
    pub fn split(self) -> (Matrix, Vec<f64>) {
        (self.mat, self.rhs)
    }

    /// Solve the system `Ax=b`, where:
    /// * `A` is `self.matrix`
    /// * `b` is `self.rhs`
    /// * `x` is the return value.
    ///
    /// Returns a `Result` containing the `Vec<f64>` representing `x` if successful.
    /// Returns an `Err` if unsuccessful.
    ///
    /// Performs LU factorization, forward and backward substitution.
    pub fn solve(mut self) -> SpResult<Vec<f64>> {
        self.mat.solve(self.rhs)
    }

    /// Read a `System` from file
    pub fn from_file(filename: &Path) -> Result<System, Box<dyn Error>> {
        use std::fs::File;
        use std::io::prelude::*;
        use std::io::BufReader;

        let mut f = File::open(filename).unwrap();
        let mut f = BufReader::new(f);
        let mut buffer = String::new();
        let mut linesize = f.read_line(&mut buffer)?;

        // Convert the first line to a title
        let title = buffer.trim().to_string();

        // Read the size/ number-format line
        buffer.clear();
        f.read_line(&mut buffer).unwrap();
        let size_strs: Vec<String> = buffer
            .split_whitespace()
            .map(|s| String::from(s))
            .collect::<Vec<String>>();
        assert(size_strs.len()).eq(2);
        let size = size_strs[0].clone().parse::<usize>().unwrap();
        assert(size).gt(0);
        let num_type_str = size_strs[1].clone();
        if num_type_str != "real" {
            return Err(NonRealNumError.into());
        }

        // Header stuff checks out.  Create our Matrix.
        let mut m = Matrix::new();

        buffer.clear();
        linesize = f.read_line(&mut buffer).unwrap();
        while linesize != 0 {
            let line_split: Vec<String> = buffer
                .split_whitespace()
                .map(|s| String::from(s))
                .collect::<Vec<String>>();
            assert(line_split.len()).eq(3);

            let x = line_split[0].clone().parse::<usize>().unwrap();
            let y = line_split[1].clone().parse::<usize>().unwrap();
            let d = line_split[2].clone().parse::<f64>().unwrap();
            assert(x).le(size);
            assert(y).le(size);

            // Alternate "done" syntax: a line of three zeroes
            if (x == 0) && (y == 0) && (d == 0.0) {
                break;
            }
            // This is an Entry.  Add it!
            m.add_element(x - 1, y - 1, d);
            // Update for next iter
            buffer.clear();
            linesize = f.read_line(&mut buffer).unwrap();
        }

        // Read the RHS vector, if present
        let mut rhs: Vec<f64> = Vec::new();
        buffer.clear();
        linesize = f.read_line(&mut buffer).unwrap();
        while linesize != 0 {
            rhs.push(buffer.trim().parse::<f64>()?);
            buffer.clear();
            linesize = f.read_line(&mut buffer)?;
        }
        if rhs.len() > 0 {
            assert(rhs.len()).eq(size);
        }

        return Ok(System {
            mat: m,
            rhs: rhs,
            title: None,
            size: size,
        });
    }
}
