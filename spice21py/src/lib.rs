use prost::Message;
///
/// Spice21 Python Bindings
///
/// Note "spice21py" must be the name of the `.so` or `.pyd` file,
/// i.e. it must be the `package` and/or `lib` name in Cargo.toml
///
use pyo3::prelude::*;
use spice21::proto::proto::Circuit;
use std::io::Cursor;

///
/// Spice21 Python Binding Module
///
/// Most methods perform some version of:
/// * Accept Python `bytes`/ Rust `&[u8]` as input
/// * Convert to Spice21 type via protobuf decoding
/// * Call a Spice21 core method
/// * Encode returned values back to bytes
///
#[pymodule]
fn spice21py(py: Python, m: &PyModule) -> PyResult<()> {
    #[pyfn(m, "health")]
    fn health_py(_py: Python) -> PyResult<String> {
        Ok(health())
    }

    #[pyfn(m, "dcop")]
    fn dcop(_py: Python, bytes_: &[u8]) -> PyResult<()> {
        // use spice21::analysis::dcop;
        let c = Circuit::decode(&mut Cursor::new(bytes_));
        // Unfortunately these conversion errors don't convert to PyResult
        let r = match c {
            Ok(ckt) => Ok(()),
            Err(e) => Ok(()),
        };
        // dcop(c);
        r
    }

    Ok(())
}

/// Health Check
fn health() -> String {
    "alive".to_string()
}
