use pyo3::exceptions::RuntimeError;
use pyo3::prelude::*;
use pyo3::{PyErr, PyResult};

use spice21::circuit::Ckt;
use std::collections::HashMap;

// Note "spice21py" must be the name of the `.so` or `.pyd` file,
// i.e. it must be the `package` and/or `lib` name in Cargo.toml

///
/// Spice21 Python Binding Module
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
        let c = Ckt::decode(bytes_);
        // Unfortunately our Error types don't ?-convert, yet
        let ckt = match c {
            Ok(ckt) => ckt,
            Err(_e) => {
                return Err(PyErr::new::<RuntimeError, String>(
                    "Circuit Decode Failed".to_string(),
                ))
            }
        };
        let res = match dcop(ckt) {
            Ok(v) => v,
            Err(_e) => {
                return Err(PyErr::new::<RuntimeError, String>(
                    "Dcop Failed".to_string(),
                ))
            }
        };
        // Return the mapping signal <-> value
        Ok(res.map)
    }
    /// Transient
    #[pyfn(m, "_tran")]
    fn tran_py(_py: Python, bytes_: &[u8]) -> PyResult<Vec<Vec<f64>>> {
        // Decode the proto-encoded circuit
        let c = Ckt::decode(bytes_);
        // Unfortunately our Error types don't ?-convert, yet
        let ckt = match c {
            Ok(ckt) => ckt,
            Err(_e) => {
                return Err(PyErr::new::<RuntimeError, String>(
                    "Circuit Decode Failed".to_string(),
                ))
            }
        };

        use spice21::analysis::{tran, TranOptions};
        let res = tran(ckt, TranOptions::default());
        return match res {
            Ok(v) => Ok(v),
            Err(_e) => Err(PyErr::new::<RuntimeError, String>(
                "Tran Failed".to_string(),
            )),
        };
    }

    Ok(())
}
