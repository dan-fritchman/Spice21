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
function tran(ckt) {
    protos.Circuit.verify(ckt);
    let buffer = protos.Circuit.encode(ckt).finish();
    let res = native._tran(buffer);
    let rv = protos.TranResult.decode(res);
    return rv.vals;
}
function ac(ckt) {
    protos.Circuit.verify(ckt);
    let buffer = protos.Circuit.encode(ckt).finish();
    let res = native._ac(buffer);
    let rv = protos.AcResult.decode(res);
    return rv.vals;
}
function sim(ckt) {
    protos.Circuit.verify(ckt);
    let buffer = protos.Sim.encode(ckt).finish();
    let res = native._sim(buffer);
    let rv = protos.SimResult.decode(res);
    return rv.vals;
}
// Public Exports 
module.exports.protos = protos;
module.exports._health = native._health;
module.exports.dcop = dcop;
module.exports._dcop = native._dcop;
module.exports.tran = tran;
module.exports._tran = native._tran;
module.exports.ac = ac;
module.exports._ac = native._ac;
module.exports.sim = sim;
module.exports._sim = native._sim;
