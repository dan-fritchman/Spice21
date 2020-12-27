// 
// Spice21 JavaScript Bindings 
// 

const native = require('../native'); // Native Rust Code 
const protos = require('./protos').spice21; // Protobuf-driven Type Definitions 

const { Op, OpResult, Tran, TranResult, Ac, AcResult, Sim, SimResult } = protos;

// Coming soon: multi-argument forms, 
// allowing provided the circuit, options, and/or arguments separately 

function dcop(arg) {
    const buffer = Op.encode(arg).finish();
    const res = native._dcop(buffer);
    const rv = OpResult.decode(res);
    return rv.vals;
}
function tran(arg) {
    const buffer = Tran.encode(arg).finish();
    const res = native._tran(buffer);
    const rv = TranResult.decode(res);
    return rv.vals;
}
function ac(arg) {
    const buffer = Ac.encode(arg).finish();
    const res = native._ac(buffer);
    const rv = AcResult.decode(res);
    return rv.vals;
}
function sim(arg) {
    const buffer = Sim.encode(arg).finish();
    const res = native._sim(buffer);
    const rv = SimResult.decode(res);
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
