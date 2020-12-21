// 
// Spice21 JavaScript Bindings 
// 

let native = require('../native'); // Native Rust Code 
let protos = require('./protos').spice21; // Protobuf-driven Type Definitions 

function dcop(ckt) {
    protos.Circuit.verify(ckt);
    let buffer = protos.Circuit.encode(ckt).finish();
    let res = native._dcop(buffer);
    let rv = protos.OpResult.decode(res);
    return rv.vals;
}

module.exports._health = native._health;
module.exports._dcop = native._dcop;
module.exports.dcop = dcop;
module.exports.protos = protos;

