// 
// Spice21Js Unit Tests 
// 

let assert = require('assert');
let spice21 = require('../lib');

describe('spice21js', function () {
    it('is alive', function () {
        let c = spice21.protos.Circuit.create({
            name: "ckt1",
            signals: ["a"],
            defs: [],
            comps: [
                { r: { name: "rr", p: "a", n: "", g: 1e-6 } },
                { c: { name: "cc", p: "a", n: "", c: 1e-6 } },
                { i: { name: "ii", p: "a", n: "", dc: 1e-6 } },
                { v: { name: "vv", p: "a", n: "", dc: 1.11 } },
            ]
        });
        
        let res = spice21.dcop(c);
        assert.strictEqual(res.a, 1.11);
        assert(res.vv < -1.1e-7);
        assert(res.vv > -1.11e-7);
    })
});

