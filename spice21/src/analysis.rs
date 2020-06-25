use std::cmp::PartialEq;
use std::convert::From;
use std::ops::{Index, IndexMut};


use std::fmt;
use std::cmp::PartialOrd;
use num::{Num,  Zero, One, Complex, Float};
use num::traits::NumAssignOps;

use crate::comps::{Component, ComponentSolver, Mos1InstanceParams, Mos1Model};
use crate::proto::{CktParse, CompParse, NodeRef};
use crate::sparse21::{Eindex, Matrix};
use crate::spresult::SpResult;
use crate::{SpNum, Abs};

/// `Stamps` are the interface between Components and Solvers.
/// Each Component returns `Stamps` from each call to `load`,
/// conveying its Matrix-contributions in `Stamps.j`
/// and its RHS contributions in `Stamps.b`.
#[derive(Debug)]
pub struct Stamps<NumT> {
    pub g: Vec<(Option<Eindex>, NumT)>,
    pub b: Vec<(Option<VarIndex>, NumT)>,
}

impl<NumT: SpNum> Stamps<NumT> {
    pub fn new() -> Stamps<NumT> {
        Stamps {
            g: vec![],
            b: vec![],
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum VarKind {
    V = 0,
    I,
}

#[derive(Debug, Clone, Copy)]
pub struct VarIndex(pub usize);

/// Conversions from Nodes and their references
impl From<NodeRef> for Option<VarIndex> {
    fn from(node: NodeRef) -> Self {
        match node {
            NodeRef::Gnd => None,
            NodeRef::Num(i) => Some(VarIndex(i)),
        }
    }
}

impl From<&NodeRef> for Option<VarIndex> {
    fn from(node: &NodeRef) -> Self {
        (*node).into()
    }
}

// This kinda thing doesn't quite work for arrays; maybe it can some day
//impl From<[NodeRef; 4]> for [Option<VarIndex>; 4] {
//    fn from(noderefs: [NodeRef; 4]) -> [Option<VarIndex>; 4] {
//        let mut opts: [Option<VarIndex>; 4] = [None, None, None, None];
//        for k in 0..noderefs.len() {
//            opts[k] = noderefs[k].into();
//        }
//        return opts;
//    }
//}

pub struct Variables<NumT> {
    kinds: Vec<VarKind>,
    values: Vec<NumT>,
}

impl<NumT: SpNum> Variables<NumT> {
    pub fn all_v(len: usize) -> Variables<NumT> {
        Variables {
            kinds: vec![VarKind::V; len],
            values: vec![NumT::zero(); len],
        }
    }
    fn add(&mut self, kind: VarKind) -> VarIndex {
        self.kinds.push(kind);
        self.values.push(NumT::zero());
        return VarIndex(self.kinds.len() - 1);
    }
    pub fn get(&self, i: Option<VarIndex>) -> NumT {
        match i {
            None => NumT::zero(),
            Some(ii) => self.values[ii.0],
        }
    }
    fn len(&self) -> usize {
        self.kinds.len()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum AnalysisMode {
    OP,
    DC,
    TRAN,
    AC,
}

struct Solver<NumT: SpNum> {
    comps: Vec<ComponentSolver>,
    vars: Variables<NumT>,
    mat: Matrix<NumT>,
    rhs: Vec<NumT>,
    history: Vec<Vec<NumT>>,
    an_mode: AnalysisMode,
}

impl Solver<f64> {
    fn new(ckt: CktParse) -> Solver<f64> {
        let mut op = Solver {
            comps: vec![],
            vars: Variables::all_v(ckt.nodes),
            mat: Matrix::new(),
            rhs: vec![],
            history: vec![],
            an_mode: AnalysisMode::OP,
        };

        // Convert each circuit-parser component into a corresponding component-solver
        for comp in ckt.comps.iter() {
            op.add_comp(comp);
        }
        // Create the corresponding matrix-elements
        for comp in op.comps.iter_mut() {
            comp.create_matrix_elems(&mut op.mat);
        }
        return op;
    }
    fn solve(&mut self, an: &AnalysisInfo) -> SpResult<Vec<f64>> {
        let mut dx = vec![0.0; self.vars.len()];

        for _k in 0..20 {
            // FIXME: number of iterations
            // Make a copy of state for tracking
            self.history.push(self.vars.values.clone());
            // Reset our matrix and RHS vector
            self.mat.reset();
            self.rhs = vec![0.0; self.vars.len()];

            // Load up component updates
            for comp in self.comps.iter_mut() {
                let updates = comp.load(&self.vars, an);
                // Make updates for G and b
                for upd in updates.g.iter() {
                    if let (Some(ei), val) = *upd {
                        self.mat.update(ei, val);
                    }
                }
                for upd in updates.b.iter() {
                    if let (Some(ei), val) = *upd {
                        self.rhs[ei.0] += val;
                    }
                }
            }
            // Calculate the residual error
            let res: Vec<f64> = self.mat.res(&self.vars.values, &self.rhs)?;
            // Check convergence
            if self.converged(&dx, &res) {
                return Ok(self.vars.values.clone());
            }
            // Solve for our update
            dx = self.mat.solve(res)?;
            let max_step = 1000e-3;
            let max_abs = dx
                .iter()
                .fold(0.0, |s, v| if v.abs() > s { v.abs() } else { s });
            if max_abs > max_step {
                for r in 0..dx.len() {
                    dx[r] = dx[r] * max_step / max_abs;
                }
            }
            // And update our guess
            for r in 0..self.vars.len() {
                self.vars.values[r] += dx[r];
            }
        }
        return Err("Convergence Failed");
    }
}

impl<NumT: SpNum> Solver<NumT> {
    fn add_comp(&mut self, comp: &CompParse) {
        match comp {
            CompParse::R(g, p, n) => {
                use crate::comps::Resistor;
                let comp = Resistor::new(*g, p.into(), n.into());
                self.comps.push(comp.into());
            }
            CompParse::C(c, p, n) => {
                use crate::comps::Capacitor;
                let comp = Capacitor::new(*c, p.into(), n.into());
                self.comps.push(comp.into());
            }
            CompParse::I(i, p, n) => {
                use crate::comps::Isrc;
                let comp = Isrc::new(*i, *p, *n);
                self.comps.push(comp.into());
            }
            CompParse::D(isat, vt, p, n) => {
                use crate::comps::Diode;
                let comp = Diode::new(*isat, *vt, *p, *n);
                self.comps.push(comp.into());
            }
            CompParse::V(v, p, n) => {
                use crate::comps::Vsrc;
                let ivar = self.vars.add(VarKind::I);
                let v = Vsrc::new(*v, p.into(), n.into(), ivar);
                self.comps.push(v.into());
            }
            CompParse::Mos(pol, g, d, s, b) => {
                use crate::comps::Mos;
                //let dp = self.vars.add(VarKind::V);
                //let sp = self.vars.add(VarKind::V);
                let comp = Mos::new(
                    [g, d, s, b].into(),
                    //                    dp.into(), ds.into(),
                    0.25,
                    50e-3,
                    3e-3,
                    *pol,
                );
                self.comps.push(comp.into());
            }
            CompParse::Mos1(model, params, g, d, s, b) => {
                use crate::comps::Mos1;
                let comp = Mos1::new(model.clone(), params.clone(), [g, d, s, b].into());
                self.comps.push(comp.into());
            }
        }
    }
    fn converged(&self, dx: &Vec<NumT>, res: &Vec<NumT>) -> bool {
        // Inter-step Newton convergence
        for e in dx.iter() {
            if e.absv() > 1e-3 {
                return false;
            }
        }
        // KCL convergence
        for e in res.iter() {
            if e.absv() > 1e-9 {
                return false;
            }
        }
        return true;
    }
}

pub fn dcop(ckt: CktParse) -> SpResult<Vec<f64>> {
    let mut s = Solver::<f64>::new(ckt);
    return s.solve(&AnalysisInfo::OP);
}

enum NumericalIntegration {
    BE,
    TRAP,
}

impl Default for NumericalIntegration {
    fn default() -> NumericalIntegration {
        NumericalIntegration::BE
    }
}

#[derive(Default)]
pub struct TranState {
    t: f64,
    dt: f64,
    vic: Vec<usize>,
    ric: Vec<usize>,
    ni: NumericalIntegration,
}

impl TranState {
    /// Numerical Integration
    pub fn integrate(&self, dq: f64, dq_dv: f64, vguess: f64, ip: f64) -> (f64, f64, f64) {
        let dt = self.dt;
        match self.ni {
            NumericalIntegration::BE => {
                let g = dq_dv / dt;
                let i = dq / dt;
                let rhs = i - g * vguess;
                (g, i, rhs)
            }
            NumericalIntegration::TRAP => {
                let g = 2.0 * dq_dv / dt;
                let i = 2.0 * dq / dt - ip;
                let rhs = i - g * vguess;
                (g, i, rhs)
            }
        }
    }
}

pub struct Tran {
    solver: Solver<f64>,
    state: TranState,
    pub opts: TranOptions,
}

pub enum AnalysisInfo<'a> {
    OP,
    TRAN(&'a TranOptions, &'a TranState),
}

impl Tran {
    pub fn new(ckt: CktParse, opts: TranOptions) -> Tran {
        return Tran {
            solver: Solver::new(ckt),
            opts,
            state: TranState {
                ..Default::default()
            },
        };
    }
    /// Create and set an initial condition on Node `n`, value `val`.
    pub fn ic(&mut self, n: NodeRef, val: f64) {
        use crate::comps::{Resistor, Vsrc};

        let fnode = self.solver.vars.add(VarKind::V);
        let ivar = self.solver.vars.add(VarKind::I);

        let mut r = Resistor::new(1.0, Some(fnode), n.into());
        r.create_matrix_elems(&mut self.solver.mat);
        self.solver.comps.push(r.into());
        self.state.ric.push(self.solver.comps.len() - 1);
        let mut v = Vsrc::new(val, Some(fnode), None, ivar);
        v.create_matrix_elems(&mut self.solver.mat);
        self.solver.comps.push(v.into());
        self.state.vic.push(self.solver.comps.len() - 1);
    }
    pub fn solve(&mut self) -> SpResult<Vec<Vec<f64>>> {
        use std::sync::mpsc;
        use std::thread;

        enum IoWriterMessage {
            STOP,
            DATA(Vec<f64>),
        }
        enum IoWriterResponse {
            OK,
            RESULT(Vec<Vec<f64>>),
        }
        let (tx, rx) = mpsc::channel::<IoWriterMessage>();
        let (tx2, rx2) = mpsc::channel::<IoWriterResponse>();

        let t = thread::spawn(move || {
            use serde::ser::{SerializeSeq, Serializer};
            use std::fs::File;

            let mut res: Vec<Vec<f64>> = vec![];

            let f = File::create("data.json").unwrap();
            let mut ser = serde_json::Serializer::new(f);
            let mut seq = ser.serialize_seq(None).unwrap();

            for msg in rx {
                match msg {
                    IoWriterMessage::DATA(d) => {
                        seq.serialize_element(&d).unwrap();
                        res.push(d);
                    }
                    IoWriterMessage::STOP => {
                        seq.end().unwrap();
                        tx2.send(IoWriterResponse::RESULT(res)).unwrap();
                        return;
                    }
                };
            }
        });
        let tsoln = self.solver.solve(&AnalysisInfo::OP);
        let tdata = match tsoln {
            Ok(x) => x,
            Err(e) => {
                println!("Failed to find initial solution");
                return Err(e);
            }
        };
        // Update initial-condition sources and resistances
        for c in self.state.vic.iter() {
            self.solver.comps[*c].update(0.0);
        }
        for c in self.state.ric.iter() {
            self.solver.comps[*c].update(1e-9);
        }
        for c in self.solver.comps.iter_mut() {
            c.tstep(&tdata);
        }
        tx.send(IoWriterMessage::DATA(tdata)).unwrap();

        self.solver.an_mode = AnalysisMode::TRAN;
        let mut tpoint: usize = 0;
        let max_tpoints: usize = 10000;
        self.state.dt = self.opts.tstep;
        while self.state.t < self.opts.tstop && tpoint < max_tpoints {
            let aninfo = AnalysisInfo::TRAN(&self.opts, &self.state);
            let tsoln = self.solver.solve(&aninfo);
            let tdata = match tsoln {
                Ok(x) => x,
                Err(e) => {
                    println!("Failed at t={}", self.state.t);
                    return Err(e);
                }
            };
            for c in self.solver.comps.iter_mut() {
                c.tstep(&tdata);
            }
            tx.send(IoWriterMessage::DATA(tdata)).unwrap();

            self.state.ni = NumericalIntegration::TRAP;
            tpoint += 1;
            self.state.t += self.opts.tstep;
        }
        tx.send(IoWriterMessage::STOP).unwrap();
        for msg in rx2 {
            match msg {
                IoWriterResponse::OK => {
                    continue;
                }
                IoWriterResponse::RESULT(res) => {
                    t.join().unwrap();
                    return Ok(res);
                }
            }
        }
        Err("Tran Results Failure")
    }
}

#[derive(Default)]
pub struct TranOptions {
    pub tstep: f64,
    pub tstop: f64,
}

pub fn tran(ckt: CktParse, opts: TranOptions) -> SpResult<Vec<Vec<f64>>> {
    return Tran::new(ckt, opts).solve();
}

// pub struct AcState {}

// pub struct AcOptions {}

// pub struct Ac {
//     solver: Solver<Complex<f64>>,
//     // state: AcState,
//     // pub opts: AcOptions,
// }

// impl Ac {
//     pub fn new(ckt: CktParse, opts: AcOptions) -> Ac {
//         return Ac {
//             solver: Solver::<Complex<f64>>::new(ckt),
//             // ..Default::default()
//         };
//     }
// }

// pub fn ac(ckt: CktParse, opts: AcOptions) -> SpResult<Vec<Vec<Complex<f64>>>> {
//     Err("????")
// }


#[cfg(test)]
mod tests {}
