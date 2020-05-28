use std::ops::{Index, IndexMut};
use std::cmp::PartialEq;

pub mod spresult;
pub mod sparse21;
pub mod assert;

use spresult::SpResult;
use sparse21::{Eindex, Matrix};


enum CompParse {
    R(f64, NodeRef, NodeRef),
    I(f64, NodeRef, NodeRef),
    V(f64, NodeRef, NodeRef),
    D(f64, f64, NodeRef, NodeRef),
    Mos(bool, NodeRef, NodeRef, NodeRef, NodeRef),
    C(f64, NodeRef, NodeRef),
}

pub struct CktParse {
    nodes: usize,
    comps: Vec<CompParse>,
}

/// Helper function to create matrix element at (row,col) if both are non-ground
fn make_matrix_elem(mat: &mut Matrix, row: NodeRef, col: NodeRef) {
    if let (NodeRef::Num(r), NodeRef::Num(c)) = (row, col) {
        mat.make(r, c);
    }
}

/// Helper function to get matrix element-pointer at (row, col) if present. or None if not 
fn get_matrix_elem(mat: &Matrix, row: NodeRef, col: NodeRef) -> Option<Eindex> {
    match (row, col) {
        (NodeRef::Num(r), NodeRef::Num(c)) => mat.get_elem(r, c),
        _ => None,
    }
}

fn get_v(x: &Vec<f64>, n: NodeRef) -> f64 {
    match n {
        NodeRef::Gnd => 0.0,
        NodeRef::Num(k) => x[k]
    }
}

trait Component {
    fn tstep(&mut self, _an: &Vec<f64>) {}
    fn update(&mut self, _val: f64) {}
    // FIXME: prob not for every Component
    fn load(&self, an: &DcOp) -> Stamps;
    fn create_matrix_elems(&self, mat: &mut Matrix);
    fn get_matrix_elems(&mut self, mat: &Matrix);
}

struct Vsrc {
    v: f64,
    p: NodeRef,
    n: NodeRef,
    ivar: NodeRef,
    pi: Option<Eindex>,
    ip: Option<Eindex>,
    ni: Option<Eindex>,
    in_: Option<Eindex>,
}

impl Vsrc {
    fn new(v: f64, p: NodeRef, n: NodeRef, ivar: NodeRef) -> Vsrc {
        Vsrc { v, p, n, ivar, pi: None, ip: None, ni: None, in_: None }
    }
}

impl Component for Vsrc {
    fn update(&mut self, val: f64) { self.v = val; }
    fn create_matrix_elems(&self, mat: &mut Matrix) {
        make_matrix_elem(mat, self.p, self.ivar);
        make_matrix_elem(mat, self.ivar, self.p);
        make_matrix_elem(mat, self.n, self.ivar);
        make_matrix_elem(mat, self.ivar, self.n);
    }
    fn get_matrix_elems(&mut self, mat: &Matrix) {
        self.pi = get_matrix_elem(mat, self.p, self.ivar);
        self.ip = get_matrix_elem(mat, self.ivar, self.p);
        self.ni = get_matrix_elem(mat, self.n, self.ivar);
        self.in_ = get_matrix_elem(mat, self.ivar, self.n);
    }
    fn load(&self, _an: &DcOp) -> Stamps {
        return Stamps {
            g: vec![
                (self.pi, 1.0),
                (self.ip, 1.0),
                (self.ni, -1.0),
                (self.in_, -1.0),
            ],
            j: vec![],
            b: vec![(self.ivar, self.v)],
        };
    }
}

struct Capacitor {
    c: f64,
    p: NodeRef,
    n: NodeRef,
    g: f64,
    i: f64,
    vp: f64,
    pp: Option<Eindex>,
    nn: Option<Eindex>,
    pn: Option<Eindex>,
    np: Option<Eindex>,
}

impl Capacitor {
    fn new(c: f64, p: NodeRef, n: NodeRef) -> Capacitor {
        Capacitor {
            c,
            p,
            n,
            g: 0.0,
            i: 0.0,
            vp: 0.0,
            pp: None,
            pn: None,
            np: None,
            nn: None,
        }
    }
}

const THE_TIMESTEP: f64 = 1e-9;

impl Component for Capacitor {
    fn create_matrix_elems(&self, mat: &mut Matrix) {
        make_matrix_elem(mat, self.p, self.p);
        make_matrix_elem(mat, self.p, self.n);
        make_matrix_elem(mat, self.n, self.p);
        make_matrix_elem(mat, self.n, self.n);
    }
    fn get_matrix_elems(&mut self, mat: &Matrix) {
        self.pp = get_matrix_elem(mat, self.p, self.p);
        self.pn = get_matrix_elem(mat, self.p, self.n);
        self.np = get_matrix_elem(mat, self.n, self.p);
        self.nn = get_matrix_elem(mat, self.n, self.n);
    }
    fn tstep(&mut self, x: &Vec<f64>) {
        let vp = get_v(x, self.p);
        let vn = get_v(x, self.n);
        let vd = vp - vn;
        self.vp = vd;
        self.i = vd * self.c / THE_TIMESTEP;
        self.g = self.c / THE_TIMESTEP;
    }
    fn load(&self, an: &DcOp) -> Stamps {
        if an.an_mode != AnalysisMode::TRAN {
            return Stamps::new();
        }

        let i = self.c * self.vp / THE_TIMESTEP;
        let g = self.c / THE_TIMESTEP;

        // FIXME: make a real index-attribute for this b-vector 
        let mut b: Vec<(NodeRef, f64)> = vec![];
        if let NodeRef::Num(_p) = self.p { b.push((self.p, i)) };
        if let NodeRef::Num(_n) = self.n { b.push((self.n, -i)) };

        return Stamps {
            j: vec![],
            g: vec![
                (self.pp, g),
                (self.nn, g),
                (self.pn, -g),
                (self.np, -g),
            ],
            b: b,
        };
    }
}

struct Resistor {
    g: f64,
    p: NodeRef,
    n: NodeRef,
    pp: Option<Eindex>,
    nn: Option<Eindex>,
    pn: Option<Eindex>,
    np: Option<Eindex>,
}

impl Resistor {
    fn new(g: f64, p: NodeRef, n: NodeRef) -> Resistor {
        Resistor { g, p, n, pp: None, pn: None, np: None, nn: None }
    }
}

impl Component for Resistor {
    fn update(&mut self, val: f64) { self.g = val; }
    fn create_matrix_elems(&self, mat: &mut Matrix) {
        make_matrix_elem(mat, self.p, self.p);
        make_matrix_elem(mat, self.p, self.n);
        make_matrix_elem(mat, self.n, self.p);
        make_matrix_elem(mat, self.n, self.n);
    }
    fn get_matrix_elems(&mut self, mat: &Matrix) {
        self.pp = get_matrix_elem(mat, self.p, self.p);
        self.pn = get_matrix_elem(mat, self.p, self.n);
        self.np = get_matrix_elem(mat, self.n, self.p);
        self.nn = get_matrix_elem(mat, self.n, self.n);
    }
    fn load(&self, _an: &DcOp) -> Stamps {
        return Stamps {
            g: vec![
                (self.pp, self.g),
                (self.nn, self.g),
                (self.pn, -self.g),
                (self.np, -self.g)
            ],
            j: vec![],
            b: vec![],
        };
    }
}


#[derive(Clone, Copy)]
enum MosTerm { G = 0, D = 1, S = 2, B = 3 }

// SPICE order: g, d, s, b
impl MosTerm {
    fn iterator() -> impl Iterator<Item=MosTerm> {
        use MosTerm::{G, D, S, B};
        [G, D, S, B].iter().copied()
    }
}

struct MosTerminals([NodeRef; 4]);

impl Index<MosTerm> for MosTerminals {
    type Output = NodeRef;
    fn index(&self, t: MosTerm) -> &NodeRef { &self.0[t as usize] }
}

struct MosMatrixPointers([[Option<Eindex>; 4]; 4]);

impl Index<(MosTerm, MosTerm)> for MosMatrixPointers {
    type Output = Option<Eindex>;
    fn index(&self, ts: (MosTerm, MosTerm)) -> &Option<Eindex> { &self.0[ts.0 as usize][ts.1 as usize] }
}

impl IndexMut<(MosTerm, MosTerm)> for MosMatrixPointers {
    fn index_mut(&mut self, ts: (MosTerm, MosTerm)) -> &mut Self::Output {
        &mut self.0[ts.0 as usize][ts.1 as usize]
    }
}

struct Mos {
    vth: f64,
    beta: f64,
    lam: f64,
    polarity: bool,
    ports: MosTerminals,
    matps: MosMatrixPointers,
}

impl Mos {
    fn new(ports: &[NodeRef; 4], vth: f64, beta: f64, lam: f64, polarity: bool) -> Mos {
        Mos {
            vth,
            beta,
            lam,
            polarity,
            ports: MosTerminals(ports.clone()),
            matps: MosMatrixPointers([[None; 4]; 4]),
        }
    }
}

impl Component for Mos {
    fn create_matrix_elems(&self, mat: &mut Matrix) {
        use MosTerm::{G, D, S};
        let matps = [(D, D), (S, S), (D, S), (S, D), (D, G), (S, G)];
        for (t1, t2) in matps.iter() {
            make_matrix_elem(mat, self.ports[*t1], self.ports[*t2]);
        }
    }
    fn get_matrix_elems(&mut self, mat: &Matrix) {
        use MosTerm::{G, D, S};
        let matps = [(D, D), (S, S), (D, S), (S, D), (D, G), (S, G)];
        for (t1, t2) in matps.iter() {
            self.matps[(*t1, *t2)] = get_matrix_elem(mat, self.ports[*t1], self.ports[*t2]);
        }
    }
    fn load(&self, an: &DcOp) -> Stamps {
        use MosTerm::{G, D, S, B};

        let vg = an.get_v(self.ports[G]);
        let vd = an.get_v(self.ports[D]);
        let vs = an.get_v(self.ports[S]);
        let vb = an.get_v(self.ports[B]);
        println!("vg={} vd={} vs={} vb={}", vg, vd, vs, vb);

        let p = if self.polarity { 1.0 } else { -1.0 };
        let vds1 = p * (vd - vs);
        let reversed = vds1 < 0.0;
        let vgs = if reversed { p * (vg - vd) } else { p * (vg - vs) };
        let vds = if reversed { -vds1 } else { vds1 };
        let vov = vgs - self.vth;

        let mut ids = 0.0;
        let mut gm = 0.0;
        let mut gds = 0.0;
        if vov <= 0.0 { // Cutoff
            // Already set
            println!("CUTOFF: vgs={} vds={}", vgs, vds);
        } else if vds >= vov { // Sat
            ids = self.beta / 2.0 * vov.powi(2) * (1.0 + self.lam * vds);
            gm = self.beta * vov * (1.0 + self.lam * vds);
            gds = self.lam * self.beta / 2.0 * vov.powi(2);
            println!("SAT: vgs={} vds={}, ids={}, gm={}, gds={}", vgs, vds, ids, gm, gds);
        } else { //Triode 
            ids = self.beta * (vov * vds - vds.powi(2) / 2.0) * (1.0 + self.lam * vds);
            gm = self.beta * vds * (1.0 + self.lam * vds);
            gds = self.beta * ((vov - vds) * (1.0 + self.lam * vds) + self.lam * ((vov * vds) - vds.powi(2) / 2.0));
            println!("LIN: vgs={} vds={}, ids={}, gm={}, gds={}", vgs, vds, ids, gm, gds);
        }
        // Sort out which are the "reported" drain and source terminals (sr, dr)
        let (sr, dr) = if !reversed { (S, D) } else { (D, S) };
        return Stamps {
            g: vec![],
            j: vec![
                (self.matps[(dr, dr)], gds),
                (self.matps[(sr, sr)], (gm + gds)),
                (self.matps[(dr, sr)], -(gm + gds)),
                (self.matps[(sr, dr)], -gds),
                (self.matps[(dr, G)], gm),
                (self.matps[(sr, G)], -gm),
            ],
            b: vec![
                (self.ports[dr], -p * ids),
                (self.ports[sr], p * ids),
            ],
        };
    }
}


struct Diode {
    isat: f64,
    vt: f64,
    p: NodeRef,
    n: NodeRef,
    pp: Option<Eindex>,
    nn: Option<Eindex>,
    pn: Option<Eindex>,
    np: Option<Eindex>,
}

impl Component for Diode {
    fn create_matrix_elems(&self, mat: &mut Matrix) {
        make_matrix_elem(mat, self.p, self.p);
        make_matrix_elem(mat, self.p, self.n);
        make_matrix_elem(mat, self.n, self.p);
        make_matrix_elem(mat, self.n, self.n);
    }
    fn get_matrix_elems(&mut self, mat: &Matrix) {
        self.pp = get_matrix_elem(mat, self.p, self.p);
        self.pn = get_matrix_elem(mat, self.p, self.n);
        self.np = get_matrix_elem(mat, self.n, self.p);
        self.nn = get_matrix_elem(mat, self.n, self.n);
    }
    fn load(&self, an: &DcOp) -> Stamps {
        let vp = an.get_v(self.p);
        let vn = an.get_v(self.n);
        let vd = (vp - vn).max(-1.5).min(1.5);
        let i = self.isat * ((vd / self.vt).exp() - 1.0);
        let di_dv = (self.isat / self.vt) * (vd / self.vt).exp();

        // FIXME: make a real index-attribute for this b-vector 
        let mut b: Vec<(NodeRef, f64)> = vec![];
        // FIXME: signs
        if let NodeRef::Num(_p) = self.p { b.push((self.p, -i)) };
        if let NodeRef::Num(_n) = self.n { b.push((self.n, -i)) };

        return Stamps {
            g: vec![],
            j: vec![
                (self.pp, di_dv),
                (self.nn, di_dv),
                (self.pn, -di_dv),
                (self.np, -di_dv)
            ],
            b: b,
        };
    }
}

struct Isrc {
    i: f64,
    p: NodeRef,
    n: NodeRef,
}

impl Component for Isrc {
    fn create_matrix_elems(&self, _mat: &mut Matrix) {}
    fn get_matrix_elems(&mut self, _mat: &Matrix) {}
    fn load(&self, _an: &DcOp) -> Stamps {
        return Stamps {
            g: vec![],
            j: vec![],
            b: vec![
                (self.p, self.i),
                (self.n, -self.i)
            ],
        };
    }
}

#[derive(Debug, Copy, Clone)]
enum NodeRef {
    Gnd,
    Num(usize),
}

#[derive(Debug)]
struct Stamps {
    g: Vec<(Option<Eindex>, f64)>,
    j: Vec<(Option<Eindex>, f64)>,
    b: Vec<(NodeRef, f64)>,
}

impl Stamps {
    fn new() -> Stamps {
        Stamps {
            g: vec![],
            j: vec![],
            b: vec![],
        }
    }
}

#[derive(Debug, Clone)]
enum Var {
    V,
    I,
}

struct Variables(Vec<Var>);

impl Variables {
    fn all_v(len: usize) -> Variables { Variables(vec![Var::V; len]) }
    fn add(&mut self, kind: Var) -> usize {
        self.0.push(kind);
        return self.0.len() - 1;
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum AnalysisMode { OP, DC, TRAN, AC }

struct DcOp {
    comps: Vec<Box<dyn Component>>,
    vars: Variables,
    mat: Matrix,
    rhs: Vec<f64>,
    x: Vec<f64>,
    history: Vec<Vec<f64>>,
    an_mode: AnalysisMode,
}

impl DcOp {
    fn add_comp(&mut self, comp: &CompParse) {
        match *comp {
            CompParse::R(g, p, n) => {
                let r = Resistor {
                    g,
                    p,
                    n,
                    pp: None,
                    pn: None,
                    np: None,
                    nn: None,
                };
                self.comps.push(Box::new(r));
            }
            CompParse::C(c, p, n) => {
                let c = Capacitor::new(c, p, n);
                self.comps.push(Box::new(c));
            }
            CompParse::I(i, p, n) => {
                let i = Isrc { i, p, n };
                self.comps.push(Box::new(i));
            }
            CompParse::D(isat, vt, p, n) => {
                let c = Diode {
                    isat,
                    vt,
                    p,
                    n,
                    pp: None,
                    pn: None,
                    np: None,
                    nn: None,
                };
                self.comps.push(Box::new(c));
            }
            CompParse::V(v, p, n) => {
                let ivar = self.vars.add(Var::I);
                let v = Vsrc {
                    v,
                    p,
                    n,
                    ivar: NodeRef::Num(ivar), // FIXME: eventually a variable-index thing
                    pi: None,
                    ip: None,
                    ni: None,
                    in_: None,
                };
                self.comps.push(Box::new(v));
            }
            CompParse::Mos(pol, g, d, s, b) => {
                let x = Mos::new(
                    &[g, d, s, b],
                    0.25, 50e-3, 3e-3, pol,
                );
                self.comps.push(Box::new(x));
            }
        }
    }
    fn get_v(&self, node: NodeRef) -> f64 {
        match node {
            NodeRef::Num(k) => self.x[k],
            NodeRef::Gnd => 0.0
        }
    }
    fn new(ckt: CktParse) -> DcOp {
        let mut op = DcOp {
            comps: vec![],
            vars: Variables::all_v(ckt.nodes),
            mat: Matrix::new(),
            x: vec![],
            rhs: vec![],
            history: vec![],
            an_mode: AnalysisMode::OP,
        };

        for comp in ckt.comps.iter() {
            op.add_comp(comp);
        }
        // Set up matrix elements per variable, and append them to Components
        // Sadly our borrow-check fighting requires two loop through the comp-list,
        // First to create matrix elements, and a second to append their references to Components.
        // I expect there's a way around this, although don't know one yet.
        for comp in op.comps.iter() {
            comp.create_matrix_elems(&mut op.mat);
        }
        for comp in op.comps.iter_mut() {
            comp.get_matrix_elems(&op.mat);
        }

        return op;
    }
    fn solve(&mut self) -> SpResult<Vec<f64>> {
        if self.x.len() == 0 { self.x = vec![0.0; self.vars.0.len()]; }
        let dx = vec![0.0; self.vars.0.len()];

        for _k in 0..20 {
            // FIXME: number of iterations
            // Make a copy of state for tracking
            self.history.push(self.x.clone());
            // Reset our matrix and RHS vector
            self.mat.reset();
            self.rhs = vec![0.0; self.vars.0.len()];

            // Load up component updates
            let mut jupdates: Vec<(Option<Eindex>, f64)> = vec![];
            for comp in self.comps.iter() {
                let updates = comp.load(&self);
                // Make updates for G and b
                for upd in updates.g.iter() {
                    if let (Some(ei), val) = *upd {
                        self.mat.update(ei, val);
                    }
                }
                for upd in updates.b.iter() {
                    if let (NodeRef::Num(ei), val) = *upd {
                        self.rhs[ei] += val;
                    }
                }
                // And save J-updates for later
                jupdates.extend(updates.j);
            }
            // Calculate the residual error
            let res: Vec<f64> = self.mat.res(&self.x, &self.rhs)?;
            // Check convergence
            if self.converged(&dx, &res) {
                return Ok(self.x.clone());
            }
            // Didn't converge, add in the Jacobian terms
            for upd in jupdates.iter() {
                if let (Some(ei), val) = *upd { self.mat.update(ei, val); }
            }
            // Solve for our update
            let mut dx = self.mat.solve(res)?;
            let max_step = 1000e-3;
            let max_abs = dx.iter().fold(0.0, |s, v| if v.abs() > s { v.abs() } else { s });
            if max_abs > max_step {
                for r in 0..dx.len() {
                    dx[r] = dx[r] * max_step / max_abs;
                }
            }
            // And update our guess
            for r in 0..self.x.len() {
                self.x[r] += dx[r];
            }
        }
        return Err("Convergence Failed");
    }
    fn converged(&self, dx: &Vec<f64>, res: &Vec<f64>) -> bool {
        // Inter-step Newton convergence
        for e in dx.iter() {
            if e.abs() > 1e-3 {
                return false;
            }
        }
        // KCL convergence
        for e in res.iter() {
            if e.abs() > 1e-9 {
                return false;
            }
        }
        return true;
    }
}

pub fn dcop(ckt: CktParse) -> SpResult<Vec<f64>> {
    let mut op = DcOp::new(ckt);
    return op.solve();
}

struct Tran {
    solver: DcOp,
    tstop: usize,
    vic: Vec<usize>,
    ric: Vec<usize>,
}

impl Tran {
    fn new(ckt: CktParse) -> Tran {
        let solver = DcOp::new(ckt);
        return Tran {
            solver,
            tstop: 20,
            vic: vec![],
            ric: vec![],
        };
    }
    fn solve(&mut self) -> SpResult<Vec<Vec<f64>>> {
        let mut res: Vec<Vec<f64>> = vec![];

        let tsoln = self.solver.solve();
        let tpoint = match tsoln {
            Ok(x) => x,
            Err(e) => {
                println!("Failed to find initial solution");
                return Err(e);
            }
        };
        for c in self.vic.iter() {
            self.solver.comps[*c].update(0.0);
        }
        for c in self.ric.iter() {
            self.solver.comps[*c].update(1e-9);
        }
        for c in self.solver.comps.iter_mut() {
            c.tstep(&tpoint);
        }
        res.push(tpoint);

        self.solver.an_mode = AnalysisMode::TRAN;
        for _t in 1..self.tstop {
            let tsoln = self.solver.solve();
            let tpoint = match tsoln {
                Ok(x) => x,
                Err(e) => {
                    println!("Failed at t={}", _t);
                    return Err(e);
                }
            };
            for c in self.solver.comps.iter_mut() {
                c.tstep(&tpoint);
            }
            res.push(tpoint);
        }
        return Ok(res);
    }
    fn ic(&mut self, n: NodeRef, val: f64) {
        let fnode = NodeRef::Num(self.solver.vars.add(Var::V));
        let ivar = NodeRef::Num(self.solver.vars.add(Var::I));

        let mut r = Resistor::new(1.0, fnode, n);
        r.create_matrix_elems(&mut self.solver.mat);
        r.get_matrix_elems(&self.solver.mat);
        self.solver.comps.push(Box::new(r));
        self.ric.push(self.solver.comps.len() - 1);
        let mut v = Vsrc::new(val, fnode, NodeRef::Gnd, ivar);
        v.create_matrix_elems(&mut self.solver.mat);
        v.get_matrix_elems(&self.solver.mat);
        self.solver.comps.push(Box::new(v));
        self.vic.push(self.solver.comps.len() - 1);
    }
}

pub fn tran(ckt: CktParse) -> SpResult<Vec<Vec<f64>>> {
    let mut tran = Tran::new(ckt);
    return tran.solve();
}


#[cfg(test)]
mod tests {
    use super::*;
    use super::assert::assert;
    use super::spresult::TestResult;

    /// Create a very basic Circuit
    fn parse_ckt() -> CktParse {
        CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(1e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Gnd),
            ],
        }
    }

    #[test]
    fn test_ckt_parse() -> TestResult {
        parse_ckt();
        Ok(())
    }

    #[test]
    fn test_dcop1() -> TestResult {
        let ckt = CktParse {
            nodes: 1,
            comps: vec![CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Gnd)],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert_eq!(soln, vec![0.0]);
        Ok(())
    }

    #[test]
    fn test_dcop2() -> TestResult {
        let ckt = parse_ckt();
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert_eq!(soln, vec![1.0, ]);
        Ok(())
    }

    #[test]
    fn test_dcop3() -> TestResult {
        // I - R - R divider
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::R(1e-3, NodeRef::Num(1), NodeRef::Num(0)),
                CompParse::I(1e-3, NodeRef::Num(1), NodeRef::Gnd),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert!((soln[0] - 1.0).abs() < 1e-4);
        assert!((soln[1] - 2.0).abs() < 1e-4);
        Ok(())
    }

    #[test]
    fn test_dcop4() -> TestResult {
        // V - R - R divider
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, NodeRef::Num(1), NodeRef::Gnd),
                CompParse::R(2e-3, NodeRef::Num(1), NodeRef::Num(0)),
                CompParse::R(2e-3, NodeRef::Num(0), NodeRef::Gnd),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert_eq!(soln, vec![0.5, 1.0, -1e-3]);
        Ok(())
    }

    #[test]
    fn test_dcop5() -> TestResult {
        // I - R - Diode
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::D(2e-16, 25e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::I(1e-3, NodeRef::Num(0), NodeRef::Gnd),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert!((soln[0] - 0.7).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_dcop5b() -> TestResult {
        // V - R - Diode
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::D(2e-16, 25e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Num(1)),
                CompParse::V(1.0, NodeRef::Num(1), NodeRef::Gnd),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert!((soln[0] - 0.7).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_dcop6() -> TestResult {
        // NMOS Char
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::V(1.0, NodeRef::Num(1), NodeRef::Gnd),
                CompParse::Mos(true, Num(0), Num(1), Gnd, Gnd),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert_eq!(soln[0], 1.0);
        assert_eq!(soln[1], 1.0);
        assert_eq!(soln[2], 0.0);
        assert!((soln[3] + 14.1e-3).abs() < 1e-4);
        Ok(())
    }

    #[test]
    fn test_dcop7() -> TestResult {
        // PMOS Char
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(-1.0, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::V(-1.0, NodeRef::Num(1), NodeRef::Gnd),
                CompParse::Mos(false, Num(0), Num(1), Gnd, Gnd),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert_eq!(soln[0], -1.0);
        assert_eq!(soln[1], -1.0);
        assert_eq!(soln[2], 0.0);
        assert!((soln[3] - 14.1e-3).abs() < 1e-4);
        Ok(())
    }

    #[test]
    fn test_dcop8() -> TestResult {
        // Diode NMOS
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(5e-3, Num(0), Gnd),
                CompParse::Mos(true, Num(0), Num(0), Gnd, Gnd),
                CompParse::R(1e-12, Num(0), Gnd), // "gmin"
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert!((soln[0] - 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_dcop8b() -> TestResult {
        // Diode NMOS, S/D Swapped
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(5e-3, Num(0), Gnd),
                CompParse::Mos(true, Num(0), Gnd, Num(0), Gnd),
                CompParse::R(1e-12, Num(0), Gnd), // "gmin"
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert!((soln[0] - 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_dcop8c() -> TestResult {
        // Diode PMOS
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(-5e-3, Num(0), Gnd),
                CompParse::Mos(false, Num(0), Num(0), Gnd, Gnd),
                CompParse::R(1e-12, Num(0), Gnd), // "gmin"
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert!((soln[0] + 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_dcop8d() -> TestResult {
        // Diode PMOS, S/D Swapped
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(-5e-3, Num(0), Gnd),
                CompParse::Mos(false, Num(0), Gnd, Num(0), Gnd),
                CompParse::R(1e-12, Num(0), Gnd), // "gmin"
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert!((soln[0] + 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_dcop9() -> TestResult {
        // NMOS-R, "Grounded"
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::R(1e-3, Num(0), Gnd),
                CompParse::Mos(true, Num(0), Num(0), Gnd, Gnd),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert_eq!(soln[0], 0.0);
        Ok(())
    }

    #[test]
    fn test_dcop9b() -> TestResult {
        // NMOS-R, "Grounded", S/D Swapped
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::R(1e-3, Num(0), Gnd),
                CompParse::Mos(true, Num(0), Gnd, Num(0), Gnd),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert_eq!(soln[0], 0.0);
        Ok(())
    }

    #[test]
    fn test_dcop9c() -> TestResult {
        // PMOS-R, "Grounded"
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::R(1e-3, Num(0), Gnd),
                CompParse::Mos(false, Num(0), Num(0), Gnd, Gnd),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert_eq!(soln[0], 0.0);
        Ok(())
    }

    #[test]
    fn test_dcop9d() -> TestResult {
        // PMOS-R, "Grounded", S/D Swapped
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::R(1e-3, Num(0), Gnd),
                CompParse::Mos(false, Num(0), Gnd, Num(0), Gnd),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert_eq!(soln[0], 0.0);
        Ok(())
    }

    #[test]
    fn test_dcop10() -> TestResult {
        // NMOS-R Inverter
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::R(1e-3, Num(1), Num(0)),
                CompParse::Mos(true, Num(0), Num(1), Gnd, Gnd),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert_eq!(soln[0], 1.0);
        assert!(soln[1] < 50e-3);
        assert!((soln[2] + 1e-3).abs() < 0.1e-3);
        Ok(())
    }

    #[test]
    fn test_dcop10b() -> TestResult {
        // PMOS-R Inverter
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(-1.0, Num(0), Gnd),
                CompParse::R(1e-3, Num(1), Num(0)),
                CompParse::Mos(false, Num(0), Num(1), Gnd, Gnd),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert_eq!(soln[0], -1.0);
        assert!(soln[1].abs() < 50e-3);
        assert!((soln[2] - 1e-3).abs() < 0.1e-3);
        Ok(())
    }

    #[test]
    fn test_dcop11() -> TestResult {
        // CMOS Inverter
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::Mos(false, Num(0), Num(1), Num(0), Num(0)),
                CompParse::Mos(true, Num(0), Num(1), Gnd, Gnd),
                CompParse::R(1e-9, Num(1), Gnd),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert_eq!(soln, vec![1.0, 0.0, 0.0]);
        Ok(())
    }

    #[test]
    fn test_dcop11b() -> TestResult {
        // CMOS Inverter
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::Mos(false, Gnd, Num(1), Num(0), Num(0)),
                CompParse::Mos(true, Gnd, Num(1), Gnd, Gnd),
                CompParse::R(1e-9, Num(1), Num(0)),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert_eq!(soln, vec![1.0, 1.0, 0.0]);
        Ok(())
    }

    #[test]
    fn test_dcop12() -> TestResult {
        // Several CMOS Inverters
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 5,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::Mos(false, Num(0), Num(1), Num(0), Num(0)),
                CompParse::Mos(true, Num(0), Num(1), Gnd, Gnd),
                CompParse::R(1e-9, Num(1), Gnd),
                CompParse::Mos(false, Num(1), Num(2), Num(0), Num(0)),
                CompParse::Mos(true, Num(1), Num(2), Gnd, Gnd),
                CompParse::R(1e-9, Num(2), Gnd),
                CompParse::Mos(false, Num(2), Num(3), Num(0), Num(0)),
                CompParse::Mos(true, Num(2), Num(3), Gnd, Gnd),
                CompParse::R(1e-9, Num(3), Gnd),
                CompParse::Mos(false, Num(3), Num(4), Num(0), Num(0)),
                CompParse::Mos(true, Num(3), Num(4), Gnd, Gnd),
                CompParse::R(1e-9, Num(4), Gnd),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert(soln[0]).eq(1.0)?;
        assert!(soln[1].abs() < 1e-3);
        assert!((soln[2] - 1.0).abs() < 1e-3);
        assert!(soln[3].abs() < 1e-3);
        assert!((soln[4] - 1.0).abs() < 1e-3);
        assert!(soln[5].abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn test_dcop13() -> TestResult {
        // RC Low-Pass Filter
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::R(1e-3, Num(1), Num(0)),
                CompParse::C(1e-9, Num(1), Gnd),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert_eq!(soln, vec![1.0, 1.0, 0.0]);
        Ok(())
    }

    #[test]
    fn test_dcop13b() -> TestResult {
        // RC High-Pass Filter
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::R(1e-3, Num(1), Gnd),
                CompParse::C(1e-9, Num(1), Num(0)),
            ],
        };
        let mut dcop = DcOp::new(ckt);
        let soln = dcop.solve()?;
        assert_eq!(soln, vec![1.0, 0.0, 0.0]);
        Ok(())
    }

    #[test]
    fn test_tran1() -> TestResult {
        // RC Low-Pass Filter
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::R(1e-3, Num(1), Num(0)),
                CompParse::C(1e-9, Num(1), Gnd),
            ],
        };
        let soln = tran(ckt)?;
        for point in soln {
            assert(point).eq(vec![1.0, 1.0, 0.0])?;
        }
        Ok(())
    }

    #[test]
    fn test_tran2() -> TestResult {
        // RC Low-Pass Filter, with Initial Condition
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                // CompParse::V(1.0, Num(1), Gnd),
                // CompParse::R(1e-3, Num(1), Num(0)),
                CompParse::I(5e-3, Num(0), Gnd),
                CompParse::C(1e-9, Num(0), Gnd),
            ],
        };
        let mut tran = Tran::new(ckt);
        tran.ic(Num(0), 0.0);
        let soln = tran.solve()?;

        assert(soln[0][0]).eq(5e-3)?;
        assert(soln[0][1]).eq(0.0)?;
        assert(soln[0][2]).eq(5e-3)?;
        for k in 1..soln.len() {
            assert(soln[k][0] - soln[k - 1][0] - 5e-3).lt(1e-6)?;
            assert(soln[k][1]).eq(0.0)?;
            assert(soln[k][2]).lt(1e-6)?;
        }
        Ok(())
    }
}

