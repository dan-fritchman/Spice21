//!
//! # Spice21 JavaScript Bindings
//!

use neon::prelude::*;
use neon::types::BinaryViewType;

/// "Health Check"
fn health_js(mut cx: FunctionContext) -> JsResult<JsString> {
    Ok(cx.string("alive"))
}

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

use spice21::proto::{Ac, CallableProto, Op, Tran};

trait JsCall: CallableProto {
    /// JavaScript call-wrapper for `CallableProto` types
    fn js_call(mut cx: FunctionContext) -> JsResult<JsBuffer> {
        // Extract the binary-encoded argument
        let mut buffer = cx.argument::<JsBuffer>(0)?;
        let bytes = cx.borrow_mut(&mut buffer, |slice| slice.as_mut_slice::<u8>());
        // Do our real work
        let rv = Self::call_bytes(bytes).unwrap();
        // And return the result as a Node Buffer/ byte-array
        rv.as_slice().to_buffer(&mut cx)
    }
}

// Apply `JsCall` to each of Spice21's proto-callable types
impl JsCall for Op {}
impl JsCall for Ac {}
impl JsCall for Tran {}

// JavaScript Exports
register_module!(mut cx, {
    cx.export_function("_health", health_js)?;
    cx.export_function("_dcop", Op::js_call)?;
    cx.export_function("_tran", Tran::js_call)?;
    cx.export_function("_ac", Ac::js_call)?;
    Ok(())
});
