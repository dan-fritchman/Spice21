// 
// Spice21Js Unit Tests 
// 

let assert = require('assert');
let axios = require('axios');
let spice21 = require('../../spice21js');


// HTTP POST Wrapper 
// Returns a Promise from axios 
const post = (path, data) => {
    const params = {
        headers: { 'Content-Type': 'application/protobuf' },
        responseType: 'arraybuffer',
    };
    return axios.post(`http://localhost:8080/${path}`, data, params);
}

describe('spice21js', function () {
    const { Circuit, OpResult, TranResult, AcResult } = spice21.protos;
    const c = Circuit.create({
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
    it('runs dcop', async () => {

        const buf = Circuit.encode(c).finish();
        const resp = await post('op', buf);
        const r = OpResult.decode(resp.data);

        const res = r.vals;
        assert.strictEqual(res.a, 1.11);
        assert(res.vv < -1.1e-7);
        assert(res.vv > -1.11e-7);
    });
    it('runs tran', async () => {

        const buf = Circuit.encode(c).finish();
        const resp = await post('tran', buf);
        const r = TranResult.decode(resp.data);

        const res = r.vals;
        assert.strictEqual(res.a.vals[0], 1.11);
        assert(res.vv.vals[0] < -1.1e-7);
        assert(res.vv.vals[0] > -1.11e-7);
    });
    it('runs ac', async () => {
        const c = Circuit.create({
            name: "ckt1",
            signals: ["a"],
            defs: [],
            comps: [
                { v: { name: "vv", p: "a", n: "", dc: 1.11, acm: 1.0} }, // A pretty simple circuit for now
            ]
        });

        const buf = Circuit.encode(c).finish();
        const resp = await post('ac', buf);
        const r = AcResult.decode(resp.data);

        const res = r.vals;
        assert.strictEqual(res.a.vals[0].re, 1.0);
    });
});

