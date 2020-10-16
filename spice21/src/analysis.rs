//! # Spice21 Analyses
//!
use serde::{Deserialize, Serialize};

use crate::circuit;
use crate::circuit::{Ckt, Comp, NodeRef};
use crate::comps::{Component, ComponentSolver};
use crate::sparse21::{Eindex, Matrix};
use crate::{sperror, SpNum, SpResult};
use num::{Complex, Float, Zero};

/// `Stamps` are the interface between Components and Solvers.
/// Each Component returns `Stamps` from each call to `load`,
/// conveying its Matrix-contributions in `Stamps.g`
/// and its RHS contributions in `Stamps.b`.
#[derive(Debug)]
pub(crate) struct Stamps<NumT> {
    pub(crate) g: Vec<(Option<Eindex>, NumT)>,
    pub(crate) b: Vec<(Option<VarIndex>, NumT)>,
}
impl<NumT: SpNum> Stamps<NumT> {
    pub fn new() -> Stamps<NumT> {
        Stamps { g: vec![], b: vec![] }
    }
}

#[derive(Debug, Clone, Copy, Deserialize, Serialize)]
pub(crate) enum VarKind {
    V = 0,
    I,
    Q,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct VarIndex(pub usize);

#[derive(Debug, Deserialize, Serialize)]
pub(crate) struct Variables<NumT> {
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
    /// Add a new Variable with attributes `name` and `kind`.
    pub fn add(&mut self, name: String, kind: VarKind) -> VarIndex {
        self.kinds.push(kind);
        self.names.push(name);
        self.values.push(NumT::zero());
        return VarIndex(self.kinds.len() - 1);
    }
    /// Retrieve a Variable value.
    /// `None` represents "ground" and always has value zero.
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
pub(crate) struct Solver<'a, NumT: SpNum> {
    pub(crate) comps: Vec<ComponentSolver<'a>>,
    pub(crate) vars: Variables<NumT>,
    pub(crate) mat: Matrix<NumT>,
    pub(crate) rhs: Vec<NumT>,
    pub(crate) history: Vec<Vec<NumT>>,
    pub(crate) models: circuit::ModelCache,
    pub(crate) opts: Options,
}

/// Real-valued Solver specifics
/// FIXME: nearly all of this *should* eventually be share-able with the Complex Solver
impl Solver<'_, f64> {
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
                return Ok(self.vars.values.clone()); // FIXME: stop cloning
            }
            // Haven't Converged. Solve for our update.
            dx = self.mat.solve(res)?;
            let max_step = 1000e-3;
            let max_abs = dx.iter().fold(0.0, |s, v| if v.abs() > s { v.abs() } else { s });
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
        return Err(sperror("Convergence Failed"));
    }
}

/// Complex-Valued Solver Specifics
/// FIXME: nearly all of this *should* eventually be share-able with the Real Solver
impl<'a> Solver<'a, Complex<f64>> {
    /// Create a Complex solver from a real-valued one.
    /// Commonly deployed when moving from DCOP to AC analysis.
    fn from(re: Solver<'a, f64>) -> Self {
        let mut op = Solver::<'a, Complex<f64>> {
            comps: re.comps,
            vars: Variables::<Complex<f64>>::from(re.vars),
            mat: Matrix::new(),
            rhs: vec![],
            history: vec![],
            models: re.models,
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
            let max_abs = dx.iter().fold(0.0, |s, v| if v.norm() > s { v.norm() } else { s });
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
        return Err(sperror("Convergence Failed"));
    }
}

impl<'a, NumT: SpNum> Solver<'a, NumT> {
    /// Create a new Solver, translate `Ckt` Components into its `ComponentSolvers`.
    pub(crate) fn new(ckt: Ckt, opts: Options) -> Solver<'a, NumT> {
        let Ckt { comps, models } = ckt;
        let mut op = Solver {
            comps: vec![],
            vars: Variables::new(),
            mat: Matrix::new(),
            rhs: vec![],
            history: vec![],
            models,
            opts,
        };

        // Convert each circuit-parser component into a corresponding component-solver
        for comp in comps.into_iter() {
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
    fn add_comp(&mut self, comp: Comp) {
        // Convert `comp` to a corresponding `ComponentSolver`
        let c: ComponentSolver = match comp {
            Comp::R(g, p, n) => {
                use crate::comps::Resistor;
                let pvar = self.node_var(p.clone());
                let nvar = self.node_var(n.clone());
                Resistor::new(g, pvar, nvar).into()
            }
            Comp::C(c, p, n) => {
                use crate::comps::Capacitor;
                let pvar = self.node_var(p.clone());
                let nvar = self.node_var(n.clone());
                Capacitor::new(c, pvar, nvar).into()
            }
            Comp::I(i, p, n) => {
                use crate::comps::Isrc;
                let pvar = self.node_var(p.clone());
                let nvar = self.node_var(n.clone());
                Isrc::new(i, pvar, nvar).into()
            }
            Comp::D(d) => {
                use crate::comps::Diode;
                Diode::from(d, self).into()
            }
            Comp::V(vs) => {
                use crate::comps::Vsrc;
                let ivar = self.vars.add(vs.name.clone(), VarKind::I);
                let vc = vs;
                let p = self.node_var(vc.p.clone());
                let n = self.node_var(vc.n.clone());
                Vsrc::new(vc.vdc, vc.acm, p, n, ivar).into()
            }
            Comp::Mos0(m) => {
                use crate::comps::mos::MosPorts;
                use crate::comps::Mos0;

                let circuit::Mos0i { name, mos_type, ports } = m;

                //let dp = self.vars.add(VarKind::V);
                //let sp = self.vars.add(VarKind::V);
                let ports: MosPorts<Option<VarIndex>> = [
                    self.node_var(ports.d.clone()),
                    self.node_var(ports.g.clone()),
                    self.node_var(ports.s.clone()),
                    self.node_var(ports.b.clone()),
                ]
                .into();
                Mos0::new(ports.into(), mos_type).into()
            }
            Comp::Mos1(m) => {
                use crate::comps::mos::MosPorts;
                use crate::comps::Mos1;

                let circuit::Mos1i { name, model, params, ports } = m;

                let ports: MosPorts<Option<VarIndex>> = [
                    self.node_var(ports.d.clone()),
                    self.node_var(ports.g.clone()),
                    self.node_var(ports.s.clone()),
                    self.node_var(ports.b.clone()),
                ]
                .into();
                Mos1::new(model.clone(), params.clone(), ports.into()).into()
            }
            Comp::Bsim4(b4i) => {
                use crate::comps::bsim4::bsim4ports::Bsim4Ports;
                use crate::comps::bsim4::Bsim4;
                use crate::comps::mos::MosPorts;

                let circuit::Bsim4i { name, ports, model, params } = b4i;
                
                // let iname = params.name.clone();// FIXME: elsewhere 
                // self.models.bsim4.add_inst(params); // FIXME: elsewhere 
                let (model, inst) = self.models.bsim4.get(&model, &params).unwrap();

                let ports: MosPorts<Option<VarIndex>> = [
                    self.node_var(ports.d.clone()),
                    self.node_var(ports.g.clone()),
                    self.node_var(ports.s.clone()),
                    self.node_var(ports.b.clone()),
                ]
                .into();
                let ports = Bsim4Ports::from(name, &ports, &model.vals, &inst.intp, self);
                Bsim4::new(ports, model, inst).into()
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
use std::collections::HashMap;
use std::ops::Index;

/// Operating Point Result
#[derive(Debug)]
pub struct OpResult {
    pub names: Vec<String>,
    pub values: Vec<f64>,
    pub map: HashMap<String, f64>,
}
impl OpResult {
    /// Create an OpResult from a (typically final) set of `Variables`.
    fn from(vars: Variables<f64>) -> Self {
        let mut map: HashMap<String, f64> = HashMap::new();
        for i in 0..vars.names.len() {
            map.insert(vars.names[i].clone(), vars.values[i]);
        }
        let Variables { names, values, .. } = vars;
        OpResult { names, values, map }
    }
    /// Get the value of signal `signame`, or an `SpError` if not present 
    pub(crate) fn get<S:Into<String>>(&self, signame: S) -> SpResult<f64> {
        match self.map.get(&signame.into()) {
            Some(v) => Ok(v.clone()), 
            None => Err(sperror("Signal Not Found")),
        }
    }
}
/// Maintain much (most?) of our original vector-result-format
/// via enabling integer indexing
impl Index<usize> for OpResult {
    type Output = f64;
    fn index(&self, t: usize) -> &f64 {
        &self.values[t]
    }
}

/// Dc Operating Point Analysis
pub fn dcop(ckt: Ckt) -> SpResult<OpResult> {
    let mut s = Solver::<f64>::new(ckt, Options::default());
    let _r = s.solve(&AnalysisInfo::OP)?;
    return Ok(OpResult::from(s.vars));
}

pub(crate) enum AnalysisInfo<'a> {
    OP,
    TRAN(&'a TranOptions, &'a TranState),
    AC(&'a AcOptions, &'a AcState),
}

pub(crate) enum NumericalIntegration {
    BE,
    TRAP,
}

impl Default for NumericalIntegration {
    fn default() -> NumericalIntegration {
        NumericalIntegration::BE
    }
}

/// # TranState
///
/// Internal state of transient analysis
/// Often passed to Components for timesteps etc
#[derive(Default)]
pub(crate) struct TranState {
    pub(crate) t: f64,
    pub(crate) dt: f64,
    pub(crate) vic: Vec<usize>,
    pub(crate) ric: Vec<usize>,
    pub(crate) ni: NumericalIntegration,
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

/// Transient Analysis Options
pub struct TranOptions {
    pub tstep: f64,
    pub tstop: f64,
    pub ic: Vec<(NodeRef, f64)>,
}

impl Default for TranOptions {
    fn default() -> TranOptions {
        TranOptions {
            tstep: 1e-6,
            tstop: 1e-3,
            ic: vec![],
        }
    }
}

pub(crate) struct Tran<'a> {
    solver: Solver<'a, f64>,
    state: TranState,
    pub(crate) opts: TranOptions,
}

impl<'a> Tran<'a> {
    pub fn new(ckt: Ckt, opts: TranOptions) -> Tran<'a> {
        let ics = opts.ic.clone();
        let mut t = Tran {
            solver: Solver::new(ckt, Options::default()),
            opts,
            state: TranState::default(),
        };
        for (node, val) in &ics {
            t.ic(node.clone(), *val);
        }
        t
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
    pub fn solve(&mut self) -> SpResult<TranResult> {
        // Initialize results
        let mut results = TranResult::new();
        results.signals(&self.solver.vars);

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
        results.push(self.state.t, &tdata);
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
            results.push(self.state.t, &tdata);
            tx.send(IoWriterMessage::DATA(tdata)).unwrap();

            self.state.ni = NumericalIntegration::TRAP;
            tpoint += 1;
            self.state.t += self.opts.tstep;
        }
        results.end();
        tx.send(IoWriterMessage::STOP).unwrap();
        for msg in rx2 {
            match msg {
                IoWriterResponse::OK => {
                    continue;
                }
                IoWriterResponse::RESULT(_) => {
                    t.join().unwrap();
                    return Ok(results);
                }
            }
        }
        Err(sperror("Tran Failure"))
    }
}

/// TranResult
/// In-Memory Store for transient data
#[derive(Default, Serialize, Deserialize)]
pub struct TranResult {
    pub signals: Vec<String>,
    pub time: Vec<f64>,
    pub data: Vec<Vec<f64>>,
    pub map: HashMap<String, Vec<f64>>,
}
impl TranResult {
    pub fn new() -> Self {
        Self {
            signals: vec![],
            time: vec![],
            data: vec![],
            map: HashMap::new(),
        }
    }
    fn signals(&mut self, vars: &Variables<f64>) {
        for name in vars.names.iter() {
            self.signals.push(name.to_string());
        }
    }
    fn push(&mut self, t: f64, vals: &Vec<f64>) {
        self.time.push(t);
        self.data.push(vals.clone());
        // FIXME: filter out un-saved and internal variables
    }
    /// Simulation complete, re-org data into hash-map of signals
    fn end(&mut self) {
        self.map.insert("time".to_string(), self.time.clone());
        for i in 0..self.signals.len() {
            let mut vals: Vec<f64> = vec![];
            for v in 0..self.time.len() {
                vals.push(self.data[v][i]);
            }
            self.map.insert(self.signals[i].clone(), vals);
        }
    }
    pub fn len(&self) -> usize {
        self.time.len()
    }
}
/// Maintain much (most?) of our original vector-result-format
/// via enabling integer indexing
impl Index<usize> for TranResult {
    type Output = Vec<f64>;
    fn index(&self, t: usize) -> &Vec<f64> {
        &self.data[t]
    }
}

/// Transient Analysis
pub fn tran(ckt: Ckt, opts: TranOptions) -> SpResult<TranResult> {
    return Tran::new(ckt, opts).solve();
}

/// Simulation Options
pub struct Options {
    pub temp: f64,
    pub nom_temp: f64,
    pub gmin: f64,
    pub abstol: f64,
    pub reltol: f64,
    pub chgtol: f64,
    pub volt_tol: f64,
    pub trtol: usize,
    pub tran_max_iter: usize,
    pub dc_max_iter: usize,
    pub dc_trcv_max_iter: usize,
    pub integrate_method: usize,
    pub order: usize,
    pub max_order: usize,
    pub pivot_abs_tol: f64,
    pub pivot_rel_tol: f64,
    pub src_factor: f64,
    pub diag_gmin: f64,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            temp: 300.15,
            nom_temp: 300.15,
            gmin: 1e-12,
            abstol: 1e-12,
            reltol: 1e-3,
            chgtol: 1e-14,
            volt_tol: 1e-6,
            trtol: 7,
            tran_max_iter: 10,
            dc_max_iter: 100,
            dc_trcv_max_iter: 50,
            integrate_method: 0,
            order: 1,
            max_order: 2,
            pivot_abs_tol: 1e-13,
            pivot_rel_tol: 1e-3,
            src_factor: 1.0,
            diag_gmin: 0.0,
        }
    }
}

#[derive(Default)]
pub(crate) struct AcState {
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

/// AcResult
/// In-Memory Store for Complex-Valued AC Data
#[derive(Default, Serialize, Deserialize)]
pub struct AcResult {
    pub signals: Vec<String>,
    pub freq: Vec<f64>,
    pub data: Vec<Vec<Complex<f64>>>,
    pub map: HashMap<String, Vec<Complex<f64>>>,
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
    /// Simulation complete, re-org data into hash-map of signals
    fn end(&mut self) {
        // self.map.insert("freq".to_string(), self.freq.clone()); // FIXME: will require multi datatypes per result
        for i in 0..self.signals.len() {
            let mut vals: Vec<Complex<f64>> = vec![];
            for v in 0..self.freq.len() {
                vals.push(self.data[v][i]);
            }
            self.map.insert(self.signals[i].clone(), vals);
        }
    }
    pub fn len(&self) -> usize {
        self.freq.len()
    }
}

/// AC Analysis
pub fn ac(ckt: Ckt, opts: AcOptions) -> SpResult<AcResult> {
    /// FIXME: result saving is in flux, and essentially on three tracks:
    /// * The in-memory format used by unit-tests returns vectors of complex numbers
    /// * The first on-disk format, streaming JSON, falls down for nested data. It has complex numbers flattened, along with frequency.
    /// * The AcResult struct holds all relevant data (in memory), and can serialize it to JSON (or any other serde format).
    ///
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

    // And return our results
    results.end();
    return Ok(results);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circuit::*;
    use crate::comps::mos::MosPorts;
    use crate::comps::MosType;
    use crate::spresult::TestResult;
    use NodeRef::{Gnd, Num};

    #[test]
    fn test_ac1() -> TestResult {
        let ckt = Ckt::from_comps(vec![Comp::R(1.0, Num(0), Gnd)]);
        let soln = ac(ckt, AcOptions::default())?;

        Ok(())
    }

    #[test]
    fn test_ac2() -> TestResult {
        use crate::circuit::Vs;
        use Comp::{C, R, V};
        let ckt = Ckt::from_comps(vec![
            R(1e-3, Num(0), Num(1)),
            C(1e-9, Num(1), Gnd),
            V(Vs {
                name: s("vi"),
                vdc: 1.0,
                acm: 1.0,
                p: Num(0),
                n: Gnd,
            }),
        ]);
        let soln = ac(ckt, AcOptions::default())?;
        Ok(())
    }

    #[test]
    #[ignore] // FIXME: aint no Mos0 AC! 
    fn test_ac3() -> TestResult {
        let ckt = Ckt::from_comps(vec![
            Comp::R(1e-3, Num(0), Num(1)),
            Comp::C(1e-9, Num(1), Gnd),
            Comp::vdc("v1", 1.0, Num(0), Gnd),
            Comp::Mos0(Mos0i {
                name: s("m"),
                mos_type: MosType::NMOS,
                ports: MosPorts {
                    g: Num(1),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
        ]);
        let soln = ac(ckt, AcOptions::default())?;

        Ok(())
    }

    // NMOS Common-Source Amp
    #[test]
    fn test_ac4() -> TestResult {
        use crate::comps::{Mos1InstanceParams, Mos1Model};

        let ckt = Ckt::from_comps(vec![
            Comp::C(1e-9, n("d"), Gnd),
            Comp::Mos1(Mos1i {
                name: s("m"),
                model: Mos1Model::default(),
                params: Mos1InstanceParams::default(),
                ports: MosPorts {
                    g: n("g"),
                    d: n("d"),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::V(Vs {
                name: s("vg"),
                vdc: 0.7,
                acm: 1.0,
                p: n("g"),
                n: Gnd,
            }),
        ]);
        let soln = ac(ckt, AcOptions::default())?;

        Ok(())
    }

    /// Diode-Connected NMOS AC
    #[test]
    fn test_ac5() -> TestResult {
        use crate::circuit::Vs;
        use crate::comps::{Mos1InstanceParams, Mos1Model};

        let ckt = Ckt::from_comps(vec![
            Comp::V(Vs {
                name: s("vd"),
                vdc: 0.5,
                acm: 1.0,
                p: Num(0),
                n: Gnd,
            }),
            Comp::Mos1(Mos1i {
                name: s("m"),
                model: Mos1Model::default(),
                params: Mos1InstanceParams::default(),
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
        ]);
        let soln = ac(ckt, AcOptions::default())?;

        Ok(())
    }
}
