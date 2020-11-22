//
// # BSIM4 Module Tests 
// 

use super::{Bsim4, Bsim4InstSpecs, Bsim4ModelSpecs, Bsim4Ports, Bsim4Cache, };
// use super::bsim4ports::Bsim4Ports;
// use super::*;

// use super::cache::Bsim4Cache;
// use super::model::Bsim4ModelSpecs;
// use super::inst::Bsim4InstSpecs;
use crate::assert::assert;
use crate::{sperror, TestResult};

use crate::analysis::{AnalysisInfo,   VarIndex};
// use crate::comps::consts::*;
use crate::comps::mos::MosType;

/// "Direct" creation of the default Bsim4Solver,
/// without the simulation runtime or solver
/// Check that the operating-point function returns something sane
#[test]
fn test_bsim4_load() -> TestResult {
    // Create & retrieve model & instance params, the normal way
    let mut cache = Bsim4Cache::new();
    cache.add_model("default", Bsim4ModelSpecs::new(MosType::NMOS));
    cache.add_inst(Bsim4InstSpecs::default());
    let (model, inst) = cache.get(&"default".to_string(), &"".to_string()).ok_or(sperror("Model Not Found"))?;

    let ports = Bsim4Ports::<Option<VarIndex>>::default();
    let mut solver = Bsim4::new(ports, model, inst);

    let p = 1.0;
    let portvs: Bsim4Ports<f64> = Bsim4Ports {
        dNode: p,
        dNodePrime: p,
        sNode: 0.0,
        sNodePrime: 0.0,
        gNodeExt: p,
        gNodePrime: p,
        gNodeMid: p,
        bNode: 0.0,
        bNodePrime: 0.0,
        dbNode: 0.0,
        sbNode: 0.0,
        qNode: 0.0,
    };

    use crate::analysis::AnalysisInfo;
    let an = AnalysisInfo::OP;
    let op = solver.op(portvs, &an);
    println!("op.cd = {:?}", op.cd);

    Ok(())
}

#[test]
fn test_bsim4_nmos_dcop1() -> TestResult {
    use crate::analysis::dcop;
    use crate::circuit::*;
    use crate::circuit::{Mosi, Ckt, Comp, NodeRef};
    use crate::comps::mos::MosPorts;
    let inst = Bsim4InstSpecs::default();
    let ports = MosPorts {
        d: n("gd"),
        g: n("gd"),
        s: NodeRef::Gnd,
        b: NodeRef::Gnd,
    };
    let mut ckt = Ckt::new();
    ckt.defs.bsim4.add_model("default", Bsim4ModelSpecs::new(MosType::NMOS));
    ckt.defs.bsim4.add_inst(Bsim4InstSpecs::default());

    ckt.add(Mosi {
        name: "bsim4".to_string(),
        ports,
        model: "default".to_string(),
        params: "".to_string(),
    });
    let p = 1.0;
    ckt.add(Comp::vdc("v1", p, n("gd"), NodeRef::Gnd));
    ckt.add(Comp::r("r1", 1e-10, n("gd"), NodeRef::Gnd));
    let soln = dcop(ckt)?;
    let vgd = soln.get("gd")?;
    assert(vgd).eq(1.0)?;
    let id = soln.get("v1")?;
    assert(id).abs().isclose(150e-6, 1e-6)?;

    Ok(())
}
#[test]
fn test_bsim4_pmos_dcop1() -> TestResult {
    use crate::analysis::dcop;
    use crate::circuit::*;
    use crate::circuit::{Mosi, Ckt, Comp, NodeRef};
    use crate::comps::mos::MosPorts;
    use NodeRef::Gnd;

    let mut ckt = Ckt::new();
    ckt.defs.bsim4.add_model("pmos", Bsim4ModelSpecs::new(MosType::PMOS));
    ckt.defs.bsim4.add_inst(Bsim4InstSpecs::default());

    ckt.add(Mosi {
        name: "bsim4".to_string(),
        ports: [n("gd"), n("gd"), Gnd, Gnd].into(),
        model: "pmos".to_string(),
        params: "".into(),
    });
    let p = -1.0;
    ckt.add(Comp::vdc("v1", p, n("gd"), NodeRef::Gnd));
    ckt.add(Comp::r("r1", 1e-10, n("gd"), NodeRef::Gnd));
    let soln = dcop(ckt)?;
    let vgd = soln.get("gd")?;
    assert(vgd).eq(-1.0)?;
    let id = soln.get("v1")?;
    assert(id).abs().isclose(57e-6, 1e-6)?;

    Ok(())
}
#[test]
fn test_bsim4_inv_dcop() -> TestResult {
    use crate::analysis::dcop;
    use crate::circuit::*;
    use crate::circuit::{Mosi, Ckt, Comp, NodeRef};
    use crate::comps::mos::{MosPorts, MosType};
    use NodeRef::Gnd;

    let mut ckt = Ckt::new();
    ckt.defs.bsim4.add_model("nmos", Bsim4ModelSpecs::new(MosType::NMOS));
    ckt.defs.bsim4.add_model("pmos", Bsim4ModelSpecs::new(MosType::PMOS));
    ckt.defs.bsim4.add_inst(Bsim4InstSpecs::default());

    ckt.add(Mosi {
        name: "p".to_string(),
        ports: ("d", "inp", "vdd", "vdd").into(),
        model: "pmos".to_string(),
        params: "".into(),
    });
    ckt.add(Mosi {
        name: "n".to_string(),
        ports: ("d", "inp", Gnd, Gnd).into(),
        model: "nmos".to_string(),
        params: "".into(),
    });
    ckt.add(Comp::vdc("vinp", 0.0, n("inp"), NodeRef::Gnd));
    ckt.add(Comp::vdc("vvdd", 1.0, n("vdd"), NodeRef::Gnd));

    let soln = dcop(ckt)?;
    let vd = soln.get("vdd")?;
    assert(vd).eq(1.0)?;
    let vd = soln.get("inp")?;
    assert(vd).eq(0.0)?;
    let vd = soln.get("d")?;
    assert(vd).gt(0.95)?;
    let id = soln.get("vinp")?;
    assert(id).abs().lt(1e-6)?;
    let id = soln.get("vvdd")?;
    assert(id).abs().lt(1e-6)?;

    Ok(())
}
#[test]
fn test_bsim4_tran1() -> TestResult {
    use crate::analysis::{tran, TranOptions};
    use crate::circuit::*;
    use crate::circuit::{Mosi, Ckt, Comp, NodeRef};
    use crate::comps::mos::{MosPorts, MosType};
    use NodeRef::Gnd;

    let mut ckt = Ckt::new();
    ckt.defs.bsim4.add_model("nmos", Bsim4ModelSpecs::new(MosType::NMOS));
    ckt.defs.bsim4.add_model("pmos", Bsim4ModelSpecs::new(MosType::PMOS));
    ckt.defs.bsim4.add_inst(Bsim4InstSpecs::default());

    ckt.add(Mosi {
        name: "p".to_string(),
        ports: [n("d"), n("inp"), n("vdd"), n("vdd")].into(),
        model: "pmos".to_string(),
        params: "".into(),
    });
    ckt.add(Mosi {
        name: "n".to_string(),
        ports: [n("d"), n("inp"), Gnd, Gnd].into(),
        model: "nmos".to_string(),
        params: "".into(),
    });
    let p = 1.0;
    ckt.add(Comp::vdc("vinp", 0.0, n("inp"), NodeRef::Gnd));
    ckt.add(Comp::vdc("vvdd", 1.0, n("vdd"), NodeRef::Gnd));

    let opts = TranOptions {
        ic: vec![(n("inp"), 1.0)],
        ..Default::default()
    };
    let soln = tran(ckt, opts)?;
    println!("{:?}", soln.map);

    Ok(())
}
