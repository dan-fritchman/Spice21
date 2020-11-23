// 
// Spice21 JavaScript Bindings 
// 

let native = require('../native'); // Native Rust Code 
let protos = require('./protos').spice21; // Protobuf-driven Type Definitions 

function dcop(ckt) {
    protos.Circuit.verify(ckt);
    let buffer = protos.Circuit.encode(ckt).finish();
    let res = native.dcop(buffer);
    return res;
}

module.exports.dcop = dcop;
module.exports.protos = protos;

