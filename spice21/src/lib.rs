use std::ops::{Index, IndexMut};

mod sparse21;
use sparse21::{Eindex, Matrix};

enum CompParse {
    R(f64, NodeRef, NodeRef),
    I(f64, NodeRef, NodeRef),
    V(f64, NodeRef, NodeRef),
    D(f64, f64, NodeRef, NodeRef),
}

struct CktParse {
    nodes: usize,
    comps: Vec<CompParse>,
}

type SpResult<T> = Result<T, &'static str>;

struct Node {
    rf: NodeRef,
    solve: bool,
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

trait Component {
    fn load(&self, an: &DcOp) -> Stamps;
    fn terminals(&self) -> Vec<NodeRef>;
    fn create_matrix_elems(&self, mat: &mut Matrix);
    fn get_matrix_elems(&mut self, mat: &Matrix);
}

struct Vsrc {
    v: f64,
    p: NodeRef,
    n: NodeRef,
    ivar: usize,
    pi: Option<Eindex>,
    ip: Option<Eindex>,
    ni: Option<Eindex>,
    in_: Option<Eindex>,
}

impl Vsrc {
    fn new(v: f64, p: NodeRef, n: NodeRef) -> Vsrc {
        Vsrc {
            v,
            p,
            n,
            ivar: 0,
            pi: None,
            ip: None,
            ni: None,
            in_: None,
        }
    }
}

impl Component for Vsrc {
    fn terminals(&self) -> Vec<NodeRef> {
        return vec![self.p, self.n];
    }
    fn create_matrix_elems(&self, mat: &mut Matrix) {
        if let NodeRef::Num(p) = self.p {
            mat.make(p, self.ivar);
            mat.make(self.ivar, p);
        }
        if let NodeRef::Num(n) = self.n {
            mat.make(n, self.ivar);
            mat.make(self.ivar, n);
        }
    }
    fn get_matrix_elems(&mut self, mat: &Matrix) {
        if let NodeRef::Num(p) = self.p {
            self.pi = mat.get_elem(p, self.ivar);
            self.ip = mat.get_elem(self.ivar, p);
        }
        if let NodeRef::Num(n) = self.n {
            self.ni = mat.get_elem(n, self.ivar);
            self.in_ = mat.get_elem(self.ivar, n);
        }
    }
    fn load(&self, an: &DcOp) -> Stamps {
        return Stamps {
            G: vec![
                (self.pi, 1.0),
                (self.ip, 1.0),
                (self.ni, -1.0),
                (self.in_, -1.0),
            ],
            J: vec![],
            b: vec![(NodeRef::Num(self.ivar), self.v)], // FIXME: NodeRef here is pretty hacky
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

impl Component for Resistor {
    fn terminals(&self) -> Vec<NodeRef> {
        return vec![self.p, self.n];
    }
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
        return Stamps {
            G: vec![
                (self.pp, self.g),
                (self.nn, self.g),
                (self.pn, -self.g),
                (self.np, -self.g)
            ],
            J: vec![],
            b: vec![],
        };
    }
}


#[derive(Clone, Copy)]
enum MosTerm { g=0, d=1, s=2, b=3 }
impl MosTerm {
    pub fn iterator() -> impl Iterator<Item = MosTerm> { 
        use MosTerm::{g, d, s, b};
        [g, d, s, b].iter().copied() 
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

struct Mos {
    vth: f64, beta: f64, lam: f64, polarity: bool, 
    ports: Vec<NodeRef>, // SPICE order: g, d, s, b
    ei: Vec<Vec<Option<Eindex>>>,
}

impl Mos {
    fn new(ports:&Vec<NodeRef>, vth:f64, beta:f64, lam:f64, polarity:bool) -> Mos {
        Mos { vth, beta, lam, polarity, ports: ports.clone(), ei: vec![vec![None; 4]; 4] }
    }
}

impl Component for Mos {
    fn terminals(&self) -> Vec<NodeRef> {
        return self.ports.clone();
    }
    fn create_matrix_elems(&self, mat: &mut Matrix) {
        for t in MosTerm::iterator() { // FIXME: make this the real thing 
            println!("{}", t as usize);
        }
        for li in 0..self.ports.len() {
            for ri in 0..self.ports.len() {
                make_matrix_elem(mat, self.ports[li], self.ports[ri]);
            }
        }
    }
    fn get_matrix_elems(&mut self, mat: &Matrix) {
        for li in 0..self.ports.len() {
            for ri in 0..self.ports.len() {
                self.ei[li][ri] = get_matrix_elem(mat, self.ports[li], self.ports[ri]);
            }
        }
    }
    fn load(&self, an: &DcOp) -> Stamps {
        let vg = an.get_v(self.ports[0]);
        let vd = an.get_v(self.ports[1]);
        let vs = an.get_v(self.ports[2]);
        let vb = an.get_v(self.ports[3]);
        let p = if self.polarity { 1.0 } else { -1.0 };
        let vds1 = p * (vd - vs);
        let reversed = vds1 < 0.0;
        let vgs = if reversed { vg - vd } else { vg - vs };
        let vds = if reversed { -vds1 } else { vds1 };
        let vov = vgs - self.vth;

        let mut ids = 0.0;
        let mut gm = 0.0;
        let mut gds = 0.0;
        if vov <= 0.0 { // Cutoff
            // Already set
        } else if vds >= vov { // Sat
            ids = self.beta / 2.0 * vov.powi(2) * (1.0 + self.lam * vds);
            gm = self.beta * vov * (1.0 + self.lam * vds);
            gds = self.lam * self.beta / 2.0 * vov.powi(2);
        } else { //Triode 
            ids = self.beta * (vov * vds - vds.powi(2) / 2.0) * (1.0 + self.lam * vds);
            gm = self.beta * vds * (1.0 + self.lam * vds);
            gds = self.beta * ((vov - vds) * (1.0 + self.lam * vds) + self.lam * ((vov * vds) - vds.powi(2) / 2.0));
        }

        let sgn = if reversed ^ self.polarity { -1.0 } else { 1.0 };
        return Stamps {
            G: vec![],
            J: vec![
                (self.ei[1][1],  gds), // FIXME: sign?
                (self.ei[2][2], -sgn * (gm + gds)),
                (self.ei[1][2],  sgn * (gm + gds)),
                (self.ei[2][1],  sgn * gds),
                (self.ei[1][0], -sgn * gm),
                (self.ei[2][0],  sgn * gm),
            ],
            b: vec![
                (self.ports[1], -sgn * ids),
                (self.ports[2],  sgn * ids),
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
    fn terminals(&self) -> Vec<NodeRef> {
        return vec![self.p, self.n];
    }
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

        // FIXME: make a real index-attribute for b-vector 
        let mut b: Vec<(NodeRef, f64)> = vec![];
        if let NodeRef::Num(p) = self.p { b.push((self.p, -i)) };
        if let NodeRef::Num(n) = self.n { b.push((self.n, -i)) };

        return Stamps {
            G: vec![],
            J: vec![
                (self.pp, di_dv),
                (self.nn, di_dv),
                (self.pn, -di_dv),
                (self.np, -di_dv)
            ],
            b: b
        };
    }
}

struct Isrc {
    i: f64,
    p: NodeRef,
    n: NodeRef,
}

impl Component for Isrc {
    fn terminals(&self) -> Vec<NodeRef> {
        return vec![self.p, self.n];
    }
    fn create_matrix_elems(&self, mat: &mut Matrix) { }
    fn get_matrix_elems(&mut self, mat: &Matrix) { }
    fn load(&self, an: &DcOp) -> Stamps {
        return Stamps {
            G: vec![],
            J: vec![],
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
    G: Vec<(Option<Eindex>, f64)>,
    J: Vec<(Option<Eindex>, f64)>,
    b: Vec<(NodeRef, f64)>,
}

impl Stamps {
    fn new() -> Stamps {
        Stamps {
            G: vec![],
            J: vec![],
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
    fn new() -> Variables {
        Variables(vec![])
    }
    fn all_v(len: usize) -> Variables {
        Variables(vec![Var::V; len])
    }
    fn add_ivar(&mut self) -> usize {
        self.0.push(Var::I);
        return self.0.len() - 1;
    }
}

struct DcOp {
    comps: Vec<Box<dyn Component>>,
    vars: Variables,
    mat: Matrix,
    rhs: Vec<f64>,
    x: Vec<f64>,
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
            CompParse::I(i, p, n) => {
                let i = Isrc { i, p, n };
                self.comps.push(Box::new(i));
            }
            CompParse::D(isat, vt, p, n) => {
                let c = Diode {
                    isat, vt, 
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
                let ivar = self.vars.add_ivar();
                let v = Vsrc {
                    v,
                    p,
                    n,
                    ivar,
                    pi: None,
                    ip: None,
                    ni: None,
                    in_: None,
                };
                self.comps.push(Box::new(v));
            }
        }
    }
    fn get_v(&self, node:NodeRef) -> f64 {
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
        for mut comp in op.comps.iter_mut() {
            comp.get_matrix_elems(&op.mat);
        }

        op.x = vec![0.0; op.vars.0.len()];
        op.rhs = vec![0.0; op.vars.0.len()];
        return op;
    }
    fn solve(&mut self) -> SpResult<Vec<f64>> {
        self.x = vec![0.0; self.rhs.len()];
        let dx = vec![0.0; self.rhs.len()];

        for k in 0..100 {
            // FIXME: number of iterations
            // Reset our matrix and RHS vector
            self.mat.reset();
            self.rhs = vec![0.0; self.rhs.len()];

            // Load up component updates
            let mut jupdates: Vec<(Option<Eindex>, f64)> = vec![];
            for comp in self.comps.iter() {
                let updates = comp.load(&self);

                // Make updates for G and b
                for upd in updates.G.iter() {
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
                jupdates.extend(updates.J);
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
            let dx = self.mat.solve(res)?;
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

#[cfg(test)]
mod tests {
    use super::*;
    type TestResult = Result<(), &'static str>;

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
        let ckt = parse_ckt();
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
        assert_eq!(soln, vec![1.0,]);
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
        assert_eq!(soln, vec![1.0, 2.0]);
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
}
