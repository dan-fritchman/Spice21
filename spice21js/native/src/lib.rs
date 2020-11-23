//!
//! # Spice21 JavaScript Bindings
//!

use neon::prelude::*;

use spice21::circuit::Ckt;

/// "Health Check"
fn health(mut cx: FunctionContext) -> JsResult<JsString> {
    Ok(cx.string("alive"))
}

/// DC Operating Point
fn dcop_js(mut cx: FunctionContext) -> JsResult<JsObject> {
    // Extract the Circuit and Options binary-encoded arguments
    let mut buffer = cx.argument::<JsBuffer>(0)?;
    let ckt_ = cx.borrow_mut(&mut buffer, |slice| slice.as_mut_slice::<u8>());

    // Decode the proto-encoded circuit
    let ckt = Ckt::decode(ckt_).unwrap();
    // println!("{:?}", ckt.comps);
    // Run DCOP
    use spice21::analysis::dcop;
    let res = dcop(ckt).unwrap();
    println!("{:?}", res);

    // Convert the result-map into a JsObject
    let object = JsObject::new(&mut cx);
    for (k, v) in res.map {
        let val = cx.number(v);
        object.set(&mut cx, &*k, val)?;
    }
    Ok(object)
}

/// Transient
fn tran_js(mut cx: FunctionContext) -> JsResult<JsObject> {
    // Extract the Circuit and Options binary-encoded arguments
    let mut buffer = cx.argument::<JsBuffer>(0)?;
    let ckt_ = cx.borrow_mut(&mut buffer, |slice| slice.as_mut_slice::<u8>());
    let mut buffer = cx.argument::<JsBuffer>(1)?;
    let opts_ = cx.borrow_mut(&mut buffer, |slice| slice.as_mut_slice::<u8>());

    // Decode the proto-encoded circuit
    let ckt = Ckt::decode(ckt_).unwrap(); 

    // Decode options, if any are provided 
    use spice21::analysis::{tran, TranOptions};
    let opts = if opts_.len() > 0 {
        TranOptions::decode(opts_).unwrap() 
    } else {
        TranOptions::default()
    };
    // Run the transient analysis
    let res = tran(ckt, opts).unwrap(); //.map_err(TempError::from)?;

    // Convert the result-map into a JsObject
    let object = JsObject::new(&mut cx);
    for (k, vec) in res.map {
        let js_array = JsArray::new(&mut cx, vec.len() as u32);
        for (i, obj) in vec.iter().enumerate() {
            let val = cx.number(*obj);
            js_array.set(&mut cx, i as u32, val).unwrap();
        }
        object.set(&mut cx, &*k, js_array)?;
    }
    Ok(object)
}

/// AC Analysis
fn ac_js(mut cx: FunctionContext) -> JsResult<JsObject> {
    // Extract the Circuit and Options binary-encoded arguments
    let mut buffer = cx.argument::<JsBuffer>(0)?;
    let ckt_ = cx.borrow_mut(&mut buffer, |slice| slice.as_mut_slice::<u8>());

    // Decode the proto-encoded circuit
    let ckt = Ckt::decode(ckt_).unwrap(); 

    // TODO: include an AcOptions argument 
    use spice21::analysis::{ac, AcOptions};
    let res = ac(ckt, AcOptions::default()).unwrap(); 

    // There are no built-in JavaScript complex numbers, nor automatic conversion provided by neon
    // Results are instead returned in the form { signame : [ [re, im], [re, im], ...], ... }
    // I.e. as Object<String, Array<Array(2)>> (mixing notations).
    let object = JsObject::new(&mut cx);
    for (k, vec) in res.map {
        let js_array = JsArray::new(&mut cx, vec.len() as u32);
        for (i, obj) in vec.iter().enumerate() {
            let re = cx.number(obj.re);
            let im = cx.number(obj.im);
            let cmpl = JsArray::new(&mut cx, 2_u32); // Length-two array
            cmpl.set(&mut cx, 0_u32, re)?;
            cmpl.set(&mut cx, 1_u32, im)?;
            js_array.set(&mut cx, i as u32, cmpl).unwrap();
        }
        object.set(&mut cx, &*k, js_array)?;
    }
    Ok(object)
}

// JavaScript Exports
register_module!(mut cx, {
    cx.export_function("health", health)?;
    cx.export_function("dcop", dcop_js)?;
    cx.export_function("tran", tran_js)?;
    cx.export_function("ac", ac_js)?;
    Ok(())
});
