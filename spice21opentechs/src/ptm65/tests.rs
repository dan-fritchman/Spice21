use spice21::analysis::{dcop, tran, OpResult, TranOptions, TranResult};
use spice21::circuit::{Ckt, NodeRef};
use spice21::{SpResult, TestResult};

fn cmos_ro3() -> Ckt {
    Ckt::from_yaml(
        r#"
            name: ro
            signals: ["1", "2", "3", vdd]
            defs: 
            - type: Module
              name: inv
              ports: [inp, out, vdd, vss] 
              params: {}
              signals: [] 
              comps: 
              - {type: M, name: p, ports: {g: inp, d: out, s: vdd, b: vdd}, params: default, model: pmos }
              - {type: M, name: n, ports: {g: inp, d: out, s: vss, b: vss}, params: default, model: nmos }
            comps:
              - {type: V, name: v1, p: vdd, n: "", dc: 1.0, acm: 0.0 }
              - {type: X, name: x1, module: inv, ports: {inp: "1",  out: "2", vdd: vdd, vss: "" }, params: {} }
              - {type: X, name: x2, module: inv, ports: {inp: "2",  out: "3", vdd: vdd, vss: "" }, params: {} }
              - {type: X, name: x3, module: inv, ports: {inp: "3",  out: "1", vdd: vdd, vss: "" }, params: {} }
        "#,
    )
    .unwrap()
}
/// Voltage-Biased, Diode-Connected PMOS
/// Model-name is `pmos`, and parameters are `default`.
fn pmos_diode_v() -> Ckt {
    Ckt::from_yaml(
        r#"
            name: nmos_diode
            signals: [gd]
            defs: []
            comps: 
            - {type: M, name: n, ports: {g: gd, d: gd, s: "", b: ""}, params: default, model: pmos }
            - {type: V, name: v1, p: gd, n: "", dc: -1.0, acm: 0.0 }
        "#,
    )
    .unwrap()
}
/// Voltage-Biased, Diode-Connected NMOS
/// Model-name is `nmos`, and parameters are `default`.
fn nmos_diode_v() -> Ckt {
    Ckt::from_yaml(
        r#"
            name: nmos_diode
            signals: [gd]
            defs: []
            comps: 
            - {type: M, name: n, ports: {g: gd, d: gd, s: "", b: ""}, params: default, model: nmos }
            - {type: V, name: v1, p: gd, n: "", dc: 1.0, acm: 0.0 }
        "#,
    )
    .unwrap()
}
#[test]
fn test1() -> TestResult {
    use super::defs;
    let mut d = defs();
    let mut ckt = cmos_ro3();
    let opts = TranOptions {
        tstep: 1e-10,
        tstop: 3e-7,
        ic: vec![(NodeRef::Num(1), 0.0)],
    };
    let soln = tran(ckt, None, Some(opts))?;
    Ok(())
}
