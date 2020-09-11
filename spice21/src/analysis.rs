use std::cmp::PartialEq;
// use std::convert::From;

use num::{Complex, Float, Zero};

use crate::comps::{Component, ComponentSolver};
use crate::proto::{CktParse, CompParse, NodeRef};
use crate::sparse21::{Eindex, Matrix};
use crate::spresult::SpResult;
use crate::SpNum;

/// `Stamps` are the interface between Components and Solvers.
/// Each Component returns `Stamps` from each call to `load`,
/// conveying its Matrix-contributions in `Stamps.g`
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
pub enum VarKind {
    V = 0,
    I,
}

#[derive(Debug, Clone, Copy)]
pub struct VarIndex(pub usize);

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
    names: Vec<String>,
}

impl<NumT: SpNum> Variables<NumT> {
    pub fn new() -> Self {
        Variables {
            kinds: vec![],
            values: vec![],
            names: vec![],
        }
    }
    /// Convert Variables<OtherT> to Variables<NumT>
    /// Keeps all `kinds`, while resetting all values to zero.
    fn from<OtherT>(other: Variables<OtherT>) -> Self {
        Variables {
            kinds: other.kinds,
            names: other.names,
            values: vec![NumT::zero(); other.values.len()],
        }
    }
    pub fn add(&mut self, name: String, kind: VarKind) -> VarIndex {
        self.kinds.push(kind);
        self.names.push(name);
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

/// Solver Iteration Struct
/// Largely for debug of convergence and progress
struct Iteration<NumT: SpNum> {
    n: usize,
    x: Vec<NumT>,
    dx: Vec<NumT>,
    vtol: Vec<bool>,
    itol: Vec<bool>,
}

/// Newton-Style Iterative Solver
/// Owns each of its circuit's ComponentSolvers,
/// its SparseMatrix, and Variables.
pub struct Solver<NumT: SpNum> {
    pub comps: Vec<ComponentSolver>,
    pub vars: Variables<NumT>,
    pub mat: Matrix<NumT>,
    pub rhs: Vec<NumT>,
    pub history: Vec<Vec<NumT>>,
    pub opts: Options,
}

/// Real-valued Solver specifics
/// FIXME: nearly all of this *should* eventually be share-able with the Complex Solver
impl Solver<f64> {
    /// Collect and incorporate updates from all components
    fn update(&mut self, an: &AnalysisInfo) {
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
    }
    fn solve(&mut self, an: &AnalysisInfo) -> SpResult<Vec<f64>> {
        let mut dx = vec![0.0; self.vars.len()];

        for _k in 0..100 {
            // FIXME: number of iterations
            // Make a copy of state for tracking
            self.history.push(self.vars.values.clone());
            // Reset our matrix and RHS vector
            self.mat.reset();
            self.rhs = vec![0.0; self.vars.len()];

            // Load up component updates
            self.update(an);

            // Calculate the residual error
            let res: Vec<f64> = self.mat.res(&self.vars.values, &self.rhs)?;

            if self.converged(&dx, &res) {
                // Check convergence
                // Commit component results
                for c in self.comps.iter_mut() {
                    c.commit();
                }
                return Ok(self.vars.values.clone());
            }
            // Haven't Converged. Solve for our update.
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

/// Complex-Valued Solver Specifics
/// FIXME: nearly all of this *should* eventually be share-able with the Real Solver
impl Solver<Complex<f64>> {
    /// Create a Complex solver from a real-valued one.
    /// Commonly deployed when moving from DCOP to AC analysis.
    fn from(re: Solver<f64>) -> Self {
        let mut op = Solver::<Complex<f64>> {
            comps: re.comps,
            vars: Variables::<Complex<f64>>::from(re.vars),
            mat: Matrix::new(),
            rhs: vec![],
            history: vec![],
            opts: re.opts,
        };

        // Create matrix elements, over-writing each Component's pointers
        for comp in op.comps.iter_mut() {
            comp.create_matrix_elems(&mut op.mat);
        }
        return op;
    }

    /// Collect and incorporate updates from all components
    fn update(&mut self, an: &AnalysisInfo) {
        for comp in self.comps.iter_mut() {
            let updates = comp.load_ac(&self.vars, an);
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
    }
    fn solve(&mut self, an: &AnalysisInfo) -> SpResult<Vec<Complex<f64>>> {
        let mut dx = vec![Complex::zero(); self.vars.len()];
        let mut iters: Vec<Iteration<Complex<f64>>> = vec![];

        for _k in 0..20 {
            // FIXME: number of iterations
            // Make a copy of state for tracking
            self.history.push(self.vars.values.clone());
            // Reset our matrix and RHS vector
            self.mat.reset();
            self.rhs = vec![Complex::zero(); self.vars.len()];

            // Load up component updates
            self.update(an);

            // Calculate the residual error
            let res: Vec<Complex<f64>> = self.mat.res(&self.vars.values, &self.rhs)?;
            let vtol: Vec<bool> = dx.iter().map(|&v| v.norm() < 1e-3).collect();
            let itol: Vec<bool> = res.iter().map(|&v| v.norm() < 1e-9).collect();
            // Check convergence
            if vtol.iter().all(|v| *v) && itol.iter().all(|v| *v) {
                // Commit component results
                for c in self.comps.iter_mut() {
                    c.commit();
                }
                return Ok(self.vars.values.clone());
            }
            // Solve for our update
            dx = self.mat.solve(res)?;
            let max_step = 1.0;
            let max_abs = dx
                .iter()
                .fold(0.0, |s, v| if v.norm() > s { v.norm() } else { s });
            if max_abs > max_step {
                for r in 0..dx.len() {
                    dx[r] = dx[r] * max_step / max_abs;
                }
            }
            // And update our guess
            for r in 0..self.vars.len() {
                self.vars.values[r] += dx[r];
            }
            iters.push(Iteration {
                n: _k,
                x: self.vars.values.clone(),
                dx: dx.clone(),
                vtol: vtol,
                itol: itol,
            });
        }
        return Err("Convergence Failed");
    }
}

impl<NumT: SpNum> Solver<NumT> {
    /// Create a new Solver, translate `CktParse` Components into its `ComponentSolvers`.
    fn new(ckt: CktParse, opts:Options) -> Solver<NumT> {
        let mut op = Solver {
            comps: vec![],
            vars: Variables::new(),
            mat: Matrix::new(),
            rhs: vec![],
            history: vec![],
            opts,
        };

        // Convert each circuit-parser component into a corresponding component-solver
        for comp in ckt.comps.into_iter() {
            op.add_comp(comp);
        }
        // Create the corresponding matrix-elements
        for comp in op.comps.iter_mut() {
            comp.create_matrix_elems(&mut op.mat);
        }
        return op;
    }
    /// Retrieve the Variable corresponding to Node `node`,
    /// creating it if necessary.
    pub fn node_var(&mut self, node: NodeRef) -> Option<VarIndex> {
        match node {
            NodeRef::Gnd => None,
            NodeRef::Name(name) => {
                // FIXME: shouldn't have to clone all the names here
                match self.vars.names.iter().cloned().position(|x| x == name) {
                    Some(i) => Some(VarIndex(i)),
                    None => Some(self.vars.add(name.clone(), VarKind::V)),
                }
            }
            NodeRef::Num(num) => {
                let name = num.to_string();
                // FIXME: shouldn't have to clone all the names here
                match self.vars.names.iter().cloned().position(|x| x == name) {
                    Some(i) => Some(VarIndex(i)),
                    None => Some(self.vars.add(name.clone(), VarKind::V)),
                }
            }
        }
    }
    /// Add parser Component `comp`, and any related Variables
    fn add_comp(&mut self, comp: CompParse) {
        // Convert `comp` to a corresponding `ComponentSolver`
        let c: ComponentSolver = match comp {
            CompParse::R(g, p, n) => {
                use crate::comps::Resistor;
                let pvar = self.node_var(p.clone());
                let nvar = self.node_var(n.clone());
                Resistor::new(g, pvar, nvar).into()
            }
            CompParse::C(c, p, n) => {
                use crate::comps::Capacitor;
                let pvar = self.node_var(p.clone());
                let nvar = self.node_var(n.clone());
                Capacitor::new(c, pvar, nvar).into()
            }
            CompParse::I(i, p, n) => {
                use crate::comps::Isrc;
                let pvar = self.node_var(p.clone());
                let nvar = self.node_var(n.clone());
                Isrc::new(i, pvar, nvar).into()
            }
            // CompParse::D(isat, vt, p, n) => {
            //     // FIXME: incorporate new parameters
            //     // FIXME: add internal resistance detection/ node insertion
            //     use crate::comps::{Diode1, DiodeInstParams, DiodeModel, DiodePorts};
            //     let p = self.node_var(p.clone());
            //     let n = self.node_var(n.clone());
            //     let model = DiodeModel::default();
            //     // Internal resistance node addition
            //     let r = if model.has_rs() {
            //         Some(self.vars.add("diode_r".to_string(), VarKind::V)) // FIXME: name
            //     } else {
            //         p
            //     };
            //     let inst = DiodeInstParams::default();
            //     Diode1::new(DiodePorts { p, n, r }, model, inst).into()
            // }
            CompParse::D1(d) => {
                use crate::comps::Diode1;
                Diode1::from(d, self).into()
            }
            CompParse::Vb(vs) => {
                use crate::comps::Vsrc;
                let ivar = self.vars.add(vs.name.clone(), VarKind::I);
                let vc = vs;
                let p = self.node_var(vc.p.clone());
                let n = self.node_var(vc.n.clone());
                Vsrc::new(vc.vdc, vc.acm, p, n, ivar).into()
            }
            CompParse::Mos0(pol, g, d, s, b) => {
                use crate::comps::Mos0;
                //let dp = self.vars.add(VarKind::V);
                //let sp = self.vars.add(VarKind::V);
                let ports = [
                    self.node_var(g.clone()),
                    self.node_var(d.clone()),
                    self.node_var(s.clone()),
                    self.node_var(b.clone()),
                ];
                Mos0::new(ports.into(), pol).into()
            }
            CompParse::Mos1(model, params, g, d, s, b) => {
                use crate::comps::Mos1;
                let ports = [
                    self.node_var(g.clone()),
                    self.node_var(d.clone()),
                    self.node_var(s.clone()),
                    self.node_var(b.clone()),
                ];
                Mos1::new(model.clone(), params.clone(), ports.into()).into()
            }
        };

        // And add to our Component vector
        self.comps.push(c);
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
    let mut s = Solver::<f64>::new(ckt, Options::default());
    return s.solve(&AnalysisInfo::OP);
}

pub enum AnalysisInfo<'a> {
    OP,
    TRAN(&'a TranOptions, &'a TranState),
    AC(&'a AcOptions, &'a AcState),
}

pub enum NumericalIntegration {
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

#[derive(Default)]
pub struct TranOptions {
    pub tstep: f64,
    pub tstop: f64,
}

pub struct Tran {
    solver: Solver<f64>,
    state: TranState,
    pub opts: TranOptions,
}

impl Tran {
    pub fn new(ckt: CktParse, opts: TranOptions) -> Tran {
        return Tran {
            solver: Solver::new(ckt, Options::default()),
            opts,
            state: TranState::default(),
        };
    }
    /// Create and set an initial condition on Node `n`, value `val`.
    pub fn ic(&mut self, n: NodeRef, val: f64) {
        use crate::comps::{Resistor, Vsrc};

        let fnode = self.solver.vars.add("vfnode".to_string(), VarKind::V); // FIXME: names
        let ivar = self.solver.vars.add("ifsrc".to_string(), VarKind::I); // FIXME: names

        let mut r = Resistor::new(1.0, Some(fnode), self.solver.node_var(n));
        r.create_matrix_elems(&mut self.solver.mat);
        self.solver.comps.push(r.into());
        self.state.ric.push(self.solver.comps.len() - 1);
        let mut v = Vsrc::new(val, 0.0, Some(fnode), None, ivar);
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

            let f = File::create("data.json").unwrap(); // FIXME: name
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
        // Solve for our initial condition
        let tsoln = self.solver.solve(&AnalysisInfo::OP);
        let tdata = match tsoln {
            Ok(x) => x,
            Err(e) => {
                println!("Failed to find initial solution");
                return Err(e);
            }
        };
        tx.send(IoWriterMessage::DATA(tdata)).unwrap();
        // Update initial-condition sources and resistances
        for c in self.state.vic.iter() {
            self.solver.comps[*c].update(0.0);
        }
        for c in self.state.ric.iter() {
            self.solver.comps[*c].update(1e-9);
        }

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

pub fn tran(ckt: CktParse, opts: TranOptions) -> SpResult<Vec<Vec<f64>>> {
    return Tran::new(ckt, opts).solve();
}

/// Simulation Options
pub struct Options {
    pub temp: f64,
    pub nomTemp: f64,
    pub gmin: f64,
    pub abstol: f64,
    pub reltol: f64,
    pub chgtol: f64,
    pub voltTol: f64,
    pub trtol: usize,
    pub tranMaxIter: usize,
    pub dcMaxIter: usize,
    pub dcTrcvMaxIter: usize,
    pub integrateMethod: NumericalIntegration,
    pub order: usize,
    pub maxOrder: usize,
    pub pivotAbsTol: f64,
    pub pivotRelTol: f64,
    pub srcFactor: f64,
    pub diagGmin: f64,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            temp: 300.15,
            nomTemp: 300.15,
            gmin: 1e-12,
            abstol: 1e-12,
            reltol: 1e-3,
            chgtol: 1e-14,
            voltTol: 1e-6,
            trtol: 7,
            tranMaxIter: 10,
            dcMaxIter: 100,
            dcTrcvMaxIter: 50,
            integrateMethod: NumericalIntegration::TRAP,
            order: 1,
            maxOrder: 2,
            pivotAbsTol: 1e-13,
            pivotRelTol: 1e-3,
            srcFactor: 1.0,
            diagGmin: 0.0,
        }
    }
}

// struct IoWriter {

// }

// enum IoWriterMessage {
//     STOP,
//     DATA(Vec<f64>),
// }
// enum IoWriterResponse {
//     OK,
//     RESULT(Vec<Vec<f64>>),
// }

// use mpsc::channel;

// impl IoWriter {
//     fn new(rx: channel::<IoWriterMessage>, tx:channel::<IoWriterResponse>) -> IoWriter {
//         use std::sync::mpsc;
//         use std::thread;

//         let (tx, rx) = mpsc::channel::<IoWriterMessage>();
//         let (tx2, rx2) = mpsc::channel::<IoWriterResponse>();

//         let t = thread::spawn(move || {
//             use serde::ser::{SerializeSeq, Serializer};
//             use std::fs::File;

//             let mut res: Vec<Vec<f64>> = vec![];

//             let f = File::create("data.json").unwrap(); // FIXME: name
//             let mut ser = serde_json::Serializer::new(f);
//             let mut seq = ser.serialize_seq(None).unwrap();

//             for msg in rx {
//                 match msg {
//                     IoWriterMessage::DATA(d) => {
//                         seq.serialize_element(&d).unwrap();
//                         res.push(d);
//                     }
//                     IoWriterMessage::STOP => {
//                         seq.end().unwrap();
//                         tx2.send(IoWriterResponse::RESULT(res)).unwrap();
//                         return;
//                     }
//                 };
//             }
//         });

//         IoWriter { }
//     }
// }

#[derive(Default)]
pub struct AcState {
    pub omega: f64,
}

/// AC Analysis Options
pub struct AcOptions {
    pub fstart: usize,
    pub fstop: usize,
    pub npts: usize, // Total, not "per decade"
}

impl Default for AcOptions {
    fn default() -> Self {
        Self {
            fstart: 1,
            fstop: 1e15 as usize,
            npts: 1000,
        }
    }
}

use serde::{Deserialize, Serialize};

/// AcResult
/// In-Memory Store for Complex-Valued AC Data
#[derive(Default, Serialize, Deserialize)]
struct AcResult {
    signals: Vec<String>,
    freq: Vec<f64>,
    data: Vec<Vec<Complex<f64>>>,
}

impl AcResult {
    fn new() -> Self {
        Self::default()
    }
    fn signals<T>(&mut self, vars: &Variables<T>) {
        for name in vars.names.iter() {
            self.signals.push(name.to_string());
        }
    }
    fn push(&mut self, f: f64, vals: &Vec<Complex<f64>>) {
        self.freq.push(f);
        self.data.push(vals.clone());
        // FIXME: filter out un-saved and internal variables
    }
}

/// AC Analysis
/// FIXME: result saving is in flux, and essentially on three tracks:
/// * The in-memory format used by unit-tests returns vectors of complex numbers
/// * The first on-disk format, streaming JSON, falls down for nested data. It has complex numbers flattened, along with frequency.
/// * The AcResult struct holds all relevant data (in memory), and can serialize it to JSON (or any other serde format).
///
pub fn ac(ckt: CktParse, opts: AcOptions) -> SpResult<Vec<Vec<Complex<f64>>>> {
    use serde::ser::{SerializeSeq, Serializer};
    use std::fs::File;

    // Initial DCOP solver and solution
    let mut solver = Solver::<f64>::new(ckt, Options::default());
    let _dc_soln = solver.solve(&AnalysisInfo::OP)?;

    // Convert to an AC solver
    let mut solver = Solver::<Complex<f64>>::from(solver);
    let mut state = AcState::default();
    let mut soln = vec![];

    // Set up streaming writer
    let rf = File::create("stream.ac.json").unwrap(); // FIXME: name
    let mut ser = serde_json::Serializer::new(rf);
    let mut seq = ser.serialize_seq(None).unwrap();

    // Initialize results
    let mut results = AcResult::new();
    results.signals(&solver.vars);

    // Set up frequency sweep
    let mut f = opts.fstart as f64;
    let fstop = opts.fstop as f64;
    let fstep = (10.0).powf(f64::log10(fstop / f) / opts.npts as f64);

    // Main Frequency Loop
    while f <= fstop {
        use std::f64::consts::PI;
        state.omega = 2.0 * PI * f;
        let an = AnalysisInfo::AC(&opts, &state);
        let fsoln = solver.solve(&an)?;

        // Push to our in-mem data
        results.push(f, &fsoln);
        // AND push to the flattened, streaming data
        let mut flat: Vec<f64> = vec![f];
        for pt in fsoln.iter() {
            flat.push(pt.re);
            flat.push(pt.im);
        }
        seq.serialize_element(&flat).unwrap();
        // AND push to our simple vector-data
        soln.push(fsoln);
        // Last-iteration handling
        if f == fstop {
            break;
        }
        f = f64::min(f * fstep, fstop);
    }
    // Close up streaming results
    SerializeSeq::end(seq).unwrap();

    // Write the full-batch JSON result
    use std::io::prelude::*;
    let mut rfj = File::create("ac.json").unwrap(); // FIXME: name
    let s = serde_json::to_string(&results).unwrap();
    rfj.write_all(s.as_bytes()).unwrap();

    // And return the flattened results
    return Ok(soln);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::comps::MosType;
    use crate::proto::{n, s, CktParse, CompParse};
    use crate::spresult::TestResult;
    use CompParse::{Mos0, C, R};
    use NodeRef::{Gnd, Num};

    #[test]
    fn test_ac1() -> TestResult {
        let ckt = CktParse {
            nodes: 1,
            comps: vec![R(1.0, Num(0), Gnd)],
        };
        let soln = ac(ckt, AcOptions::default())?;

        Ok(())
    }

    #[test]
    fn test_ac2() -> TestResult {
        use crate::proto::Vs;
        use CompParse::{Vb, C, R};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                R(1e-3, Num(0), Num(1)),
                C(1e-9, Num(1), Gnd),
                Vb(Vs {
                    name: s("vi"),
                    vdc: 1.0,
                    acm: 1.0,
                    p: Num(0),
                    n: Gnd,
                }),
            ],
        };
        let soln = ac(ckt, AcOptions::default())?;
        Ok(())
    }

    #[test]
    fn test_ac3() -> TestResult {
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                R(1e-3, Num(0), Num(1)),
                C(1e-9, Num(1), Gnd),
                CompParse::V(1.0, Num(0), Gnd),
                Mos0(MosType::NMOS, Num(1), Num(0), Gnd, Gnd),
            ],
        };
        let soln = ac(ckt, AcOptions::default())?;

        Ok(())
    }

    // NMOS Common-Source Amp
    #[test]
    fn test_ac4() -> TestResult {
        use crate::comps::{Mos1InstanceParams, Mos1Model};
        use crate::proto::Vs;
        use CompParse::{Mos1, Vb, C, R};

        let ckt = CktParse {
            nodes: 3,
            comps: vec![
                R(1e-5, n("vdd"), n("d")),
                C(1e-9, n("d"), Gnd),
                Mos1(
                    Mos1Model::default(),
                    Mos1InstanceParams::default(),
                    n("g"),
                    n("d"),
                    Gnd,
                    Gnd,
                ),
                CompParse::V(1.0, n("vdd"), Gnd),
                Vb(Vs {
                    name: s("vg"),
                    vdc: 0.7,
                    acm: 1.0,
                    p: n("g"),
                    n: Gnd,
                }),
            ],
        };
        let soln = ac(ckt, AcOptions::default())?;

        Ok(())
    }

    /// Diode-Connected NMOS AC
    #[test]
    fn test_ac5() -> TestResult {
        use crate::comps::{Mos1InstanceParams, Mos1Model};
        use crate::proto::Vs;
        use CompParse::{Mos1, Vb};

        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                Vb(Vs {
                    name: s("vd"),
                    vdc: 0.5,
                    acm: 1.0,
                    p: Num(0),
                    n: Gnd,
                }),
                Mos1(
                    Mos1Model::default(),
                    Mos1InstanceParams::default(),
                    Num(0),
                    Num(0),
                    Gnd,
                    Gnd,
                ),
            ],
        };
        let soln = ac(ckt, AcOptions::default())?;

        Ok(())
    }
}
