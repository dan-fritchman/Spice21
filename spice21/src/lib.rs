use sparse21::{Eindex, Matrix};

type SpResult<T> = Result<T, &'static str>;

struct Node {
    rf: NodeRef,
    solve: bool,
}

trait Component {
    // fn setup(&mut self, mat: &mut Matrix);
    fn load(&self) -> Stamps;

    fn terminals(&self) -> Vec<NodeRef>;
    fn create_matrix_elems(&self, mat: &mut Matrix);
    
    fn get_matrix_elems(&mut self, mat: &Matrix);

    fn on_add(&self, ckt: &Circuit) {}
}

// struct Vsrc {
//     v: f64,
//     p: NodeRef,
//     n: NodeRef,
//     ivar: usize,
//     iindex: usize,
//     piv:Option<Eindex>,
//     ivp:Option<Eindex>,
//     niv:Option<Eindex>,
//     ivn:Option<Eindex>,
// }

// impl Vsrc {
//     fn new (v:f64, p: NodeRef, n:NodeRef) -> Vsrc {
//         Vsrc { v, p, n, ivar:0, iindex:0, piv: None, ivp: None, niv: None, ivn: None}
//     }
// }

// impl Component for Vsrc {
//     fn on_add(&mut self, ckt: &mut Circuit) {
//         self.ivar = ckt.add_var();
//     }
//     fn setup(&mut self, mat: &mut Matrix) {
//         self.iindex = self.ivar + mat.nodes.len();

//     }
//     fn load (&self) -> Stamps {
//         let mut s = Stamps::new();

//         return s;
//     }
// }

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
            mat.add_element(l, r, 0.0);
            mat.add_element(r, l, 0.0);
        }
        if let NodeRef::Num(p) = self.p { 
            mat.add_element(p, p, 0.0);
        }
        if let NodeRef::Num(n) = self.n { 
            mat.add_element(n, n, 0.0);
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
            mat.add_element(l, r, 0.0);
            mat.add_element(r, l, 0.0);
        }
        if let NodeRef::Num(p) = self.p { 
            mat.add_element(p, p, 0.0);
        }
        if let NodeRef::Num(n) = self.n { 
            mat.add_element(n, n, 0.0);
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


struct Circuit {
    comps: Vec<Box<dyn Component>>,
    nodes: Vec<Node>,
    vars: usize, // FIXME: more elaborate non-node variable handling
    node0: Node,
}

impl Circuit {
    fn new() -> Circuit {
        Circuit {
            comps: vec![],
            nodes: vec![],
            vars: 0,
            node0: Node {
                rf: NodeRef::Gnd,
                solve: false,
            },
        }
    }
    fn add_node(&mut self) -> NodeRef {
        let rf = NodeRef::Num(self.nodes.len());
        let node = Node {
            rf: rf.clone(),
            solve: true,
        };
        self.nodes.push(node);
        return rf;
    }
    fn add_var(&mut self) -> usize {
        let nv = self.vars;
        self.vars += 1;
        return nv;
    }
    fn add_comp<C: Component + 'static>(&mut self, comp: C) {
        comp.on_add(&self);
        self.comps.push(Box::new(comp));
    }
}

#[derive(Copy, Clone)]
enum NodeRef {
    Gnd,
    Num(usize),
}

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

struct DcOp {
    ckt: Circuit,
    mat: Matrix,
    rhs: Vec<f64>,
    x: Vec<f64>,
}

impl DcOp {
    fn new(mut ckt: Circuit) -> DcOp {
        let size = ckt.nodes.len() + ckt.vars;
        let mut mat = Matrix::new();
        
        // Sadly our borrow-check fighting requires two loop through the comp-list,
        // First to create matrix elements, and a second to append their references to Components. 
        // I expect there's a way around this, although don't know one yet. 
        for mut comp in ckt.comps.iter() {
            comp.create_matrix_elems(&mut mat);
        }
        for mut comp in ckt.comps.iter_mut() {
            comp.get_matrix_elems(&mat);
        }

        return DcOp {
            ckt,
            mat,
            x: vec![0.0; size],
            rhs: vec![0.0; size],
        };
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
            let mut jupdates: Vec<(Eindex, f64)> = vec![];
            for comp in self.ckt.comps.iter() {
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
            let mut res: Vec<f64> = self.mat.vecmul(&self.x)?;
            for k in 0..self.x.len() {
                res[k] = -1.0 * res[k];
                res[k] += self.rhs[k];
            }
            // Check convergence
            if self.converged(&dx, &res) {
                return Ok(self.x.clone());
            }
            // Didn't converge, add Jacobian terms
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
            println!("DX ENTRY:");
            println!("{}", *e);
            if *e > 1e-3 {
                return false;
            }
        }
        // KCL convergence
        for e in res.iter() {
            println!("RES ENTRY:");
            println!("{}", *e);
            if *e > 1e-9 {
                return false;
            }
        }
        return true;
    }
}

fn op(ckt: Circuit, x0: Option<Vec<f64>>) -> SpResult<Vec<f64>> {
    // Create the circuit - components, nodes, variables
    // Walk it to get matrix elements
    // Create a matrix of elements
    // Pointers to the matrix elements associate with the components, somehow
    //
    let mut dcop = DcOp::new(ckt);
    return dcop.solve();
}

#[cfg(test)]
mod tests {
    use super::*;
    type TestResult = Result<(), &'static str>;

    #[test]
    fn test_add_node() -> TestResult {
        let mut c = Circuit::new();
        let rf = c.add_node();
        match rf {
            NodeRef::Gnd => return Err("gnd?"),
            NodeRef::Num(n) => assert_eq!(n, 0),
        };
        assert_eq!(c.nodes.len(), 1);
        assert_eq!(c.comps.len(), 0);
        Ok(())
    }

    #[test]
    fn test_add_res() -> TestResult {
        let mut c = Circuit::new();
        let rf = c.add_node();
        let r = Resistor {
            g: 1e-3,
            p: rf,
            n: NodeRef::Gnd,
            pp: None,
            nn: None,
            np: None,
            pn: None,
        };
        c.add_comp(r);
        Ok(())
    }

    #[test]
    fn test_dcop1() -> TestResult {
        let mut c = Circuit::new();
        let rf = c.add_node();
        let r = Resistor {
            g: 1e-3,
            p: rf,
            n: NodeRef::Gnd,
            pp: None,
            nn: None,
            np: None,
            pn: None,
        };
        c.add_comp(r);

        let mut dcop = DcOp::new(c);
        let soln = dcop.solve()?;
        assert_eq!(soln, vec![0.0]);
        Ok(())
    }

    #[test]
    fn test_dcop2() -> TestResult {
        // Current source + resistor(!)
        let mut c = Circuit::new();
        let rf = c.add_node();
        let r = Resistor {
            g: 1e-3,
            p: rf,
            n: NodeRef::Gnd,
            pp: None,
            nn: None,
            np: None,
            pn: None,
        };
        c.add_comp(r);
        let i = Isrc {
            i: 1e-3,
            p: rf,
            n: NodeRef::Gnd,

            pp: None,
            nn: None,
            np: None,
            pn: None,
        };
        c.add_comp(i);

        let mut dcop = DcOp::new(c);
        let soln = dcop.solve()?;
        assert_eq!(soln, vec![1.0]);
        Ok(())
    }

    #[test]
    fn test_dcop3() -> TestResult {
        let mut c = Circuit::new();
        let n1 = c.add_node();
        let n2 = c.add_node();
        let r1 = Resistor {
            g: 1e-3,
            p: n1,
            n: NodeRef::Gnd,
            pp: None,
            nn: None,
            np: None,
            pn: None,
        };
        c.add_comp(r1);
        let r2 = Resistor {
            g: 1e-3,
            p: n2,
            n: n1,
            pp: None,
            nn: None,
            np: None,
            pn: None,
        };
        c.add_comp(r2);
        let i = Isrc {
            i: 1e-3,
            p: n2,
            n: NodeRef::Gnd,

            pp: None,
            nn: None,
            np: None,
            pn: None,
        };
        c.add_comp(i);

        let mut dcop = DcOp::new(c);
        let soln = dcop.solve()?;
        assert_eq!(soln, vec![1.0, 2.0]);
        Ok(())
    }
}
