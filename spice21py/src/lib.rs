//!
//! # Spice21 Python Bindings
//!

use pyo3::exceptions::RuntimeError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use pyo3::{PyErr, PyResult};

use spice21rs::proto::{Ac, CallableProto, Op, Tran};
use spice21rs::SpResult;

// Note "spice21py" must be the name of the `.so` or `.pyd` file,
// i.e. it must be the `package` and/or `lib` name in Cargo.toml

/// Convert results and errors
///
/// This might look nicer as a method. But defining one is a pain,
/// as neither `SpResult` nor `PyResult` are defined in this crate.
fn res(py: Python, sp: SpResult<Vec<u8>>) -> PyResult<PyObject> {
    match sp {
        Ok(bytes) => Ok(PyBytes::new(py, &bytes).into()),
        Err(e) => Err(PyErr::new::<RuntimeError, String>(e.desc)),
    }
}
///
/// # Spice21 Python Binding Module
///
/// Most methods perform some version of:
///
/// * Accept Python `bytes`/ Rust `&[u8]` as input
///   * Serialization/ encoding is done on the Python side
///   * PyO3 performs the `bytes` -> `&[u8]` conversion
/// * Convert to Spice21 type via protobuf decoding
/// * Call a Spice21 core method
/// * Convert the Rust-bytes (`Vec<u8>`) into `PyBites`
///   (PyO3 doesn't auto-generate this for the return-value)
/// * Return Python `bytes` encoding of a Spice21 result-type
/// * Decoding can then be performed in Python
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
    fn dcop_py(py: Python, bytes: &[u8]) -> PyResult<PyObject> {
        res(py, Op::call_bytes(bytes))
    }

    /// Transient
    #[pyfn(m, "_tran")]
    fn tran_py(py: Python, bytes: &[u8]) -> PyResult<PyObject> {
        res(py, Tran::call_bytes(bytes))
    }

    /// AC Analysis
    #[pyfn(m, "_ac")]
    fn ac_py(py: Python, bytes: &[u8]) -> PyResult<PyObject> {
        res(py, Ac::call_bytes(bytes))
    }

    Ok(())
}
