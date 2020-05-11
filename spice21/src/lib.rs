use sparse21::{Eindex, Matrix};

enum CompParse {
    R(f64, NodeRef, NodeRef),
    I(f64, NodeRef, NodeRef),
    V(f64, NodeRef, NodeRef),
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

trait Component {
    fn load(&self) -> Stamps;
    fn terminals(&self) -> Vec<NodeRef>;
    fn create_matrix_elems(&self, mat: &mut Matrix);
    fn get_matrix_elems(&mut self, mat: &Matrix);
    fn add_vars(&mut self, vars: &mut Variables) {}
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
    fn load(&self) -> Stamps {
        let mut v: Vec<(Eindex, f64)> = vec![];
        if let Some(pi) = self.pi {
            v.push((pi, 1.0));
        }
        if let Some(ip) = self.ip {
            v.push((ip, 1.0));
        }
        if let Some(ni) = self.ni {
            v.push((ni, -1.0));
        }
        if let Some(in_) = self.in_ {
            v.push((in_, -1.0));
        }
        return Stamps {
            G: v,
            J: vec![],
            b: vec![(self.ivar, self.v)],
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

fn get_matrix_elem(mat: &Matrix, row: NodeRef, col: NodeRef) -> Option<Eindex> {
    match (row, col) {
        (NodeRef::Num(r), NodeRef::Num(c)) => mat.get_elem(r, c),
        _ => None,
    }
}
impl Component for Resistor {
    fn terminals(&self) -> Vec<NodeRef> {
        return vec![self.p, self.n];
    }
    fn create_matrix_elems(&self, mat: &mut Matrix) {
        if let (NodeRef::Num(l), NodeRef::Num(r)) = (self.p, self.n) {
            mat.make(l, r);
            mat.make(r, l);
        }
        if let NodeRef::Num(p) = self.p {
            mat.make(p, p);
        }
        if let NodeRef::Num(n) = self.n {
            mat.make(n, n);
        }
    }
    fn get_matrix_elems(&mut self, mat: &Matrix) {
        self.pp = get_matrix_elem(mat, self.p, self.p);
        self.pn = get_matrix_elem(mat, self.p, self.n);
        self.np = get_matrix_elem(mat, self.n, self.p);
        self.nn = get_matrix_elem(mat, self.n, self.n);
    }
    fn load(&self) -> Stamps {
        let mut v: Vec<(Eindex, f64)> = vec![];
        if let Some(pp) = self.pp {
            v.push((pp, self.g));
        }
        if let Some(nn) = self.nn {
            v.push((nn, self.g));
        }
        if let Some(pn) = self.pn {
            v.push((pn, -1.0 * self.g));
        }
        if let Some(np) = self.np {
            v.push((np, -1.0 * self.g));
        }
        return Stamps {
            G: v,
            J: vec![],
            b: vec![],
        };
    }
}

struct Isrc {
    i: f64,
    p: NodeRef,
    n: NodeRef,

    pp: Option<Eindex>,
    nn: Option<Eindex>,
    pn: Option<Eindex>,
    np: Option<Eindex>,
}

impl Component for Isrc {
    fn terminals(&self) -> Vec<NodeRef> {
        return vec![self.p, self.n];
    }
    fn create_matrix_elems(&self, mat: &mut Matrix) {
        if let (NodeRef::Num(l), NodeRef::Num(r)) = (self.p, self.n) {
            mat.make(l, r);
            mat.make(r, l);
        }
        if let NodeRef::Num(p) = self.p {
            mat.make(p, p);
        }
        if let NodeRef::Num(n) = self.n {
            mat.make(n, n);
        }
    }

    fn get_matrix_elems(&mut self, mat: &Matrix) {
        self.pp = get_matrix_elem(mat, self.p, self.p);
        self.pn = get_matrix_elem(mat, self.p, self.n);
        self.np = get_matrix_elem(mat, self.n, self.p);
        self.nn = get_matrix_elem(mat, self.n, self.n);
    }
    // fn setup(&mut self, mat: &mut Matrix) {}
    fn load(&self) -> Stamps {
        let mut b: Vec<(usize, f64)> = vec![];
        if let NodeRef::Num(pp) = self.p {
            b.push((pp, self.i))
        }
        if let NodeRef::Num(nn) = self.n {
            b.push((nn, -1.0 * self.i))
        }
        return Stamps {
            G: vec![],
            J: vec![],
            b: b,
        };
    }
}

#[derive(Copy, Clone)]
enum NodeRef {
    Gnd,
    Num(usize),
}

#[derive(Debug)]
struct Stamps {
    G: Vec<(Eindex, f64)>,
    J: Vec<(Eindex, f64)>,
    b: Vec<(usize, f64)>,
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
                let i = Isrc {
                    i,
                    p,
                    n,
                    pp: None,
                    pn: None,
                    np: None,
                    nn: None,
                };
                self.comps.push(Box::new(i));
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
    fn new(ckt: CktParse) -> DcOp {
        let mut op = DcOp {
            // ckt: Circuit::new(),
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

        for k in 0..100 { // FIXME: number of iterations
            // Reset our matrix and RHS vector
            self.mat.reset();
            self.rhs = vec![0.0; self.rhs.len()];

            // Load up component updates
            let mut jupdates: Vec<(Eindex, f64)> = vec![];
            for comp in self.comps.iter() {
                let updates = comp.load();

                // Make updates for G and b
                for upd in updates.G.iter() {
                    self.mat.update(upd.0, upd.1);
                }
                for upd in updates.b.iter() {
                    self.rhs[upd.0] += upd.1;
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
                self.mat.update(upd.0, upd.1);
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
            if *e > 1e-3 {
                return false;
            }
        }
        // KCL convergence
        for e in res.iter() {
            if *e > 1e-9 {
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
                // CompParse::V(1.00, NodeRef::Num(0), NodeRef::Gnd),
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
        // Getting crazy: I - R - R divider
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
        // Getting crazy: I - R - R divider
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
}
