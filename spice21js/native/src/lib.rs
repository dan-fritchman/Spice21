//!
//! # Spice21 JavaScript Bindings
//!

use neon::prelude::*;
use neon::types::BinaryViewType;

///
/// Node-Buffer Conversion Trait
///
/// Courtesy the Neon team.
/// Converts binary array/vector-ish Rust-things (`Vec<u8>`, `&[u8]`) to Node's binary `Buffer`.
/// Hopefully this might some day be part of Neon itself.
///
trait ToBuffer {
    fn to_buffer<'a, C>(&self, cx: &mut C) -> JsResult<'a, JsBuffer>
    where
        C: Context<'a>;
}
impl<T: BinaryViewType + Copy> ToBuffer for &[T] {
    fn to_buffer<'a, C>(&self, cx: &mut C) -> JsResult<'a, JsBuffer>
    where
        C: Context<'a>,
    {
        let mut buf = cx.buffer(self.len() as u32)?;
        cx.borrow_mut(&mut buf, |buf| {
            let buf = buf.as_mut_slice();
            buf.copy_from_slice(self);
        });
        Ok(buf)
    }
}

use spice21::analysis;
use spice21::circuit::Ckt;

/// "Health Check"
fn health_js(mut cx: FunctionContext) -> JsResult<JsString> {
    Ok(cx.string("alive"))
}

/// DC Operating Point
fn dcop_js(mut cx: FunctionContext) -> JsResult<JsBuffer> {
    // Extract the Circuit and Options binary-encoded arguments
    let mut buffer = cx.argument::<JsBuffer>(0)?;
    let ckt_ = cx.borrow_mut(&mut buffer, |slice| slice.as_mut_slice::<u8>());

    // Decode the proto-encoded circuit
    let ckt = Ckt::decode(ckt_).unwrap();
    // Run DCOP
    let res = analysis::dcop(ckt).unwrap();
    // Convert the result back to bytes
    let rv = res.encode();
    // And return them as a Node Buffer/ byte-array
    rv.as_slice().to_buffer(&mut cx) 
}

/// Transient
fn tran_js(mut cx: FunctionContext) -> JsResult<JsBuffer> {
    // Extract the Circuit and Options binary-encoded arguments
    let mut buffer = cx.argument::<JsBuffer>(0)?;
    let ckt_ = cx.borrow_mut(&mut buffer, |slice| slice.as_mut_slice::<u8>());
    // let mut buffer = cx.argument::<JsBuffer>(1)?;
    // let opts_ = cx.borrow_mut(&mut buffer, |slice| slice.as_mut_slice::<u8>());

    // Decode the proto-encoded circuit
    let ckt = Ckt::decode(ckt_).unwrap();

    // Decode options, if any are provided
    use spice21::analysis::{tran, TranOptions};
    // let opts = if opts_.len() > 0 {
    //     TranOptions::decode(opts_).unwrap()
    // } else {
    //     TranOptions::default()
    // };
    let opts = TranOptions::default();
    // Run the transient analysis
    let res = tran(ckt, opts).unwrap(); //.map_err(TempError::from)?;

    // Convert the result back to bytes
    let rv = res.encode();
    // And return them as a Node Buffer/ byte-array
    rv.as_slice().to_buffer(&mut cx) 
}

/// AC Analysis
fn ac_js(mut cx: FunctionContext) -> JsResult<JsBuffer> {
    // Extract the Circuit and Options binary-encoded arguments
    let mut buffer = cx.argument::<JsBuffer>(0)?;
    let ckt_ = cx.borrow_mut(&mut buffer, |slice| slice.as_mut_slice::<u8>());

    // Decode the proto-encoded circuit
    let ckt = Ckt::decode(ckt_).unwrap();

    // TODO: include an AcOptions argument
    use spice21::analysis::{ac, AcOptions};
    let res = ac(ckt, AcOptions::default()).unwrap();

    // Convert the result back to bytes
    let rv = res.encode();
    // And return them as a Node Buffer/ byte-array
    rv.as_slice().to_buffer(&mut cx) 
}

// JavaScript Exports
register_module!(mut cx, {
    cx.export_function("_health", health_js)?;
    cx.export_function("_dcop", dcop_js)?;
    cx.export_function("_tran", tran_js)?;
    cx.export_function("_ac", ac_js)?;
    Ok(())
});
