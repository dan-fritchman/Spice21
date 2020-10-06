//!
//! # Spice21 Python Bindings 
//! 
use std::collections::HashMap;

use pyo3::exceptions::RuntimeError;
use pyo3::prelude::*;
use pyo3::{PyErr, PyResult};

use spice21::circuit::Ckt;
use spice21::SpError;

// Note "spice21py" must be the name of the `.so` or `.pyd` file,
// i.e. it must be the `package` and/or `lib` name in Cargo.toml

// Error-Type Conversions
// Unfortunately requires this temporary middle-man defined internally,
// And a few calls to map_err(TempError::from)
struct TempError(SpError);
impl From<SpError> for TempError {
    fn from(err: SpError) -> TempError {
        TempError(err)
    }
}
impl From<TempError> for PyErr {
    fn from(err: TempError) -> PyErr {
        PyErr::new::<RuntimeError, String>(err.0.desc)
    }
}

///
/// # Spice21 Python Binding Module
///
/// Most methods perform some version of:
/// * Accept Python `bytes`/ Rust `&[u8]` as input
///   * (Serialization/ encoding is done on the Python side)
/// * Convert to Spice21 type via protobuf decoding
/// * Call a Spice21 core method
/// * Return core-generated values, either as Python primitive types
///   or as protobuf-encoded byte-string.
///
#[pymodule]
fn spice21py(_py: Python, m: &PyModule) -> PyResult<()> {
    /// "Health Check"
    #[pyfn(m, "health")]
    fn health_py(_py: Python) -> PyResult<String> {
        Ok("alive".to_string())
    }

    /// DC Operating Point
    #[pyfn(m, "_dcop")]
    fn dcop_py(_py: Python, bytes_: &[u8]) -> PyResult<HashMap<String, f64>> {
        use spice21::analysis::dcop;
        // Decode the proto-encoded circuit
        let ckt = Ckt::decode(bytes_).map_err(TempError::from)?;
        // Run DCOP
        let res = dcop(ckt).map_err(TempError::from)?;
        // And return the mapping signal <-> value
        Ok(res.map)
    }

    /// Transient
    #[pyfn(m, "_tran")]
    fn tran_py(_py: Python, bytes_: &[u8]) -> PyResult<HashMap<String, Vec<f64>>> {
        use spice21::analysis::{tran, TranOptions};
        // Decode the proto-encoded circuit
        let ckt = Ckt::decode(bytes_).map_err(TempError::from)?;
        let res = tran(ckt, TranOptions::default()).map_err(TempError::from)?;
        Ok(res.map)
    }

    /// AC Analysis 
    #[pyfn(m, "_ac")]
    fn ac_py(_py: Python, bytes_: &[u8]) -> PyResult<HashMap<String, Vec<(f64, f64)>>> {
        use spice21::analysis::{ac, AcOptions};
        // Decode the proto-encoded circuit
        let ckt = Ckt::decode(bytes_).map_err(TempError::from)?;
        let res = ac(ckt, AcOptions::default()).map_err(TempError::from)?;
        // PyO3 doesn't quite understand Rust's complex numbers, so we decode a bit 
        let mut map: HashMap<String, Vec<(f64, f64)>> = HashMap::new();
        for (name, vals) in &res.map {
            let mut v: Vec<(f64, f64)> = Vec::new();
            for k in 0..vals.len(){
                v.push((vals[k].re, vals[k].im));
            }
            map.insert(name.to_string(), v);
          }
        Ok(map)
    }

    Ok(())
}
