use std::convert::From;
use std::ops::{Index, IndexMut};
use std::cmp::PartialEq;

pub mod spresult;
pub mod sparse21;
pub mod assert;

use spresult::SpResult;
use sparse21::{Eindex, Matrix};

use prost::Message;


// Include the `items` module, which is generated from items.proto.
pub mod proto {
    include!(concat!(env!("OUT_DIR"), "/spice21.proto.rs"));
    fn some() -> Circuit {
        return Circuit {
            name: String::from("tbd"),
            statements: vec![
//                Statement { statement: Some(Statement::Instance(Instance { name: String::from("???") })) }
            ],
        };
    }
}

enum CompParse {
    R(f64, NodeRef, NodeRef),
    I(f64, NodeRef, NodeRef),
    V(f64, NodeRef, NodeRef),
    D(f64, f64, NodeRef, NodeRef),
    Mos(bool, NodeRef, NodeRef, NodeRef, NodeRef),
    Mos1(Mos1Model, Mos1InstanceParams, NodeRef, NodeRef, NodeRef, NodeRef),
    C(f64, NodeRef, NodeRef),
}

pub struct CktParse {
    nodes: usize,
    comps: Vec<CompParse>,
}

/// Helper function to create matrix element at (row,col) if both are non-ground
fn make_matrix_elem(mat: &mut Matrix, row: Option<VarIndex>, col: Option<VarIndex>) -> Option<Eindex> {
    if let (Some(r), Some(c)) = (row, col) {
        return Some(mat.make(r.0, c.0));
    }
    return None;
}

/// Helper function to get matrix element-pointer at (row, col) if present. or None if not 
fn get_matrix_elem(mat: &Matrix, row: Option<VarIndex>, col: Option<VarIndex>) -> Option<Eindex> {
    match (row, col) {
        (Some(r), Some(c)) => mat.get_elem(r.0, c.0),
        _ => None,
    }
}

fn get_v(x: &Vec<f64>, n: Option<VarIndex>) -> f64 {
    match n {
        None => 0.0,
        Some(k) => x[k.0]
    }
}

trait Component {
    fn tstep(&mut self, _guess: &Vec<f64>) {}

    fn update(&mut self, _val: f64) {}
    // FIXME: prob not for every Component

    fn load(&mut self, guess: &Variables, an: &AnalysisInfo) -> Stamps;
    fn create_matrix_elems(&mut self, mat: &mut Matrix);
}

struct Vsrc {
    v: f64,
    p: Option<VarIndex>,
    n: Option<VarIndex>,
    ivar: VarIndex,
    pi: Option<Eindex>,
    ip: Option<Eindex>,
    ni: Option<Eindex>,
    in_: Option<Eindex>,
}

impl Vsrc {
    fn new(v: f64, p: Option<VarIndex>, n: Option<VarIndex>, ivar: VarIndex) -> Vsrc {
        Vsrc { v, p, n, ivar, pi: None, ip: None, ni: None, in_: None }
    }
}

impl Component for Vsrc {
    fn update(&mut self, val: f64) { self.v = val; }
    fn create_matrix_elems(&mut self, mat: &mut Matrix) {
        self.pi = make_matrix_elem(mat, self.p, Some(self.ivar));
        self.ip = make_matrix_elem(mat, Some(self.ivar), self.p);
        self.ni = make_matrix_elem(mat, self.n, Some(self.ivar));
        self.in_ = make_matrix_elem(mat, Some(self.ivar), self.n);
    }
    fn load(&mut self, _guess: &Variables, _an: &AnalysisInfo) -> Stamps {
        return Stamps {
            g: vec![
                (self.pi, 1.0),
                (self.ip, 1.0),
                (self.ni, -1.0),
                (self.in_, -1.0),
            ],
            b: vec![(Some(self.ivar), self.v)],
        };
    }
}

/// Numerical Integration
/// Calculated (G,I) from charge and time changes and derivative
fn num_integ(q: f64, qp: f64, dq_dv: f64, vguess: f64, ip: Option<f64>, dt: f64) -> (f64, f64, f64) {
    // You can have any method you want, as long as its Backward Euler
    // let g = dq_dv / dt;
    // let i = (q - qp) / dt;
    // let rhs = i - g * vguess;

    // You can have any method you want, as long as its TRAP
    let g = 2.0 * dq_dv / dt;
    let i = 2.0 * (q - qp) / dt - ip.unwrap();
    let rhs = i - g * vguess;
    return (g, i, rhs);
}

#[derive(Default)]
struct Capacitor {
    c: f64,
    p: Option<VarIndex>,
    n: Option<VarIndex>,
    pp: Option<Eindex>,
    nn: Option<Eindex>,
    pn: Option<Eindex>,
    np: Option<Eindex>,
    op: CapOpPoint,
    guess: CapOpPoint,
}


#[derive(Clone, Copy, Default)]
struct CapOpPoint {
    v: f64,
    q: f64,
    i: f64,
}

impl Capacitor {
    fn new(c: f64, p: Option<VarIndex>, n: Option<VarIndex>) -> Capacitor {
        Capacitor {
            c,
            p,
            n,
            ..Default::default()
        }
    }
    fn q(&self, v: f64) -> f64 {
        return self.c * v;
    }
    fn dq_dv(&self, _v: f64) -> f64 {
        return self.c;
    }
}

const THE_TIMESTEP: f64 = 1e-9;


impl Component for Capacitor {
    fn create_matrix_elems(&mut self, mat: &mut Matrix) {
        self.pp = make_matrix_elem(mat, self.p, self.p);
        self.pn = make_matrix_elem(mat, self.p, self.n);
        self.np = make_matrix_elem(mat, self.n, self.p);
        self.nn = make_matrix_elem(mat, self.n, self.n);
    }
    /// Load our last guess as the new operating point
    fn tstep(&mut self, _x: &Vec<f64>) {
        self.op = self.guess;
    }
    fn load(&mut self, guess: &Variables, an: &AnalysisInfo) -> Stamps {
        let vd = guess.get(self.p) - guess.get(self.n);
        let q = self.q(vd);

        let _data = match *an {
            AnalysisInfo::TRAN(d) => d,
            _ => {
                // FIXME: calculating this during DCOP, so we copy cleanly afterward
                // Should probably just find a way to calculate it then
                self.guess = CapOpPoint { v: vd, q: q, i: 0.0 };
                return Stamps::new();
            }
        };
        let (g, i, rhs) = num_integ(q, self.op.q, self.dq_dv(vd), vd, Some(self.op.i), THE_TIMESTEP);
        self.guess = CapOpPoint { v: vd, q: q, i: i };

        return Stamps {
            g: vec![
                (self.pp, g),
                (self.nn, g),
                (self.pn, -g),
                (self.np, -g),
            ],
            b: vec![
                (self.p, -rhs),
                (self.n, rhs)
            ],
        };
    }
}


#[derive(Clone, Copy)]
enum TwoTerm { P = 0, N = 1 }

struct TwoTerminals([Option<VarIndex>; 2]);

impl Index<TwoTerm> for TwoTerminals {
    type Output = Option<VarIndex>;
    fn index(&self, t: TwoTerm) -> &Option<VarIndex> { &self.0[t as usize] }
}

struct TwoTermMatrixPointers([[Option<Eindex>; 2]; 2]);

impl Index<(TwoTerm, TwoTerm)> for TwoTermMatrixPointers {
    type Output = Option<Eindex>;
    fn index(&self, ts: (TwoTerm, TwoTerm)) -> &Option<Eindex> {
        &self.0[ts.0 as usize][ts.1 as usize]
    }
}

impl IndexMut<(TwoTerm, TwoTerm)> for TwoTermMatrixPointers {
    fn index_mut(&mut self, ts: (TwoTerm, TwoTerm)) -> &mut Self::Output {
        &mut self.0[ts.0 as usize][ts.1 as usize]
    }
}

struct Resistor {
    g: f64,
    terms: TwoTerminals,
    matps: TwoTermMatrixPointers,
}

impl Resistor {
    fn new(g: f64, p: Option<VarIndex>, n: Option<VarIndex>) -> Resistor {
        Resistor { g, terms: TwoTerminals([p, n]), matps: TwoTermMatrixPointers([[None; 2]; 2]) }
    }
}

impl Component for Resistor {
    fn update(&mut self, val: f64) { self.g = val; }
    fn create_matrix_elems(&mut self, mat: &mut Matrix) {
        use TwoTerm::{P, N};
        for l in [P, N].into_iter() {
            for r in [P, N].into_iter() {
                self.matps[(*l, *r)] = make_matrix_elem(mat, self.terms[*l], self.terms[*r]);
            }
        }
    }
    fn load(&mut self, _guess: &Variables, _an: &AnalysisInfo) -> Stamps {
        use TwoTerm::{P, N};
        return Stamps {
            g: vec![
                (self.matps[(P, P)], self.g),
                (self.matps[(N, N)], self.g),
                (self.matps[(P, N)], -self.g),
                (self.matps[(N, P)], -self.g)
            ],
            b: vec![],
        };
    }
}


/// Mos Terminals, in SPICE order: g, d, s, b
#[derive(Clone, Copy)]
enum MosTerm { G = 0, D = 1, S = 2, B = 3 }

struct MosTerminals([Option<VarIndex>; 4]);

impl Index<MosTerm> for MosTerminals {
    type Output = Option<VarIndex>;
    fn index(&self, t: MosTerm) -> &Option<VarIndex> { &self.0[t as usize] }
}

impl From<[&NodeRef; 4]> for MosTerminals {
    fn from(n: [&NodeRef; 4]) -> Self {
        return MosTerminals(
            [n[0].into(), n[1].into(), n[2].into(), n[3].into()]
        );
    }
}

#[derive(Default)]
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

#[derive(Clone)]
enum MosType { NMOS, PMOS }

// FIXME: should reference instead of cloning, when we can get it to play nicely with `Box<dyn Comp>`
#[derive(Clone)]
struct Mos1Model {
    mos_type: MosType,
    vt0: f64,
    kp: f64,
    gamma: f64,
    phi: f64,
    lambda: f64,
    rd: f64,
    rs: f64,
    cbd: f64,
    cbs: f64,
    is: f64,
    pb: f64,
    cgso: f64,
    cgdo: f64,
    cgbo: f64,
    rsh: f64,
    cj: f64,
    mj: f64,
    cjsw: f64,
    mjsw: f64,
    js: f64,
    tox: f64,
    ld: f64,
    u0: f64,
    fc: f64,
    nsub: f64,
    tpg: bool,
    nss: f64,
    tnom: f64,
    kf: f64,
    af: f64,
}

impl Default for Mos1Model {
    fn default() -> Self {
        Mos1Model {
            mos_type: MosType::NMOS,
            vt0: 0.0,
            kp: 2.0e-5,
            gamma: 0.0,
            phi: 0.6,
            lambda: 0.0,
            rd: 0.0,
            rs: 0.0,
            cbd: 0.0,
            cbs: 0.0,
            is: 1.0e-14,
            pb: 0.8,
            cgso: 0.0,
            cgdo: 0.0,
            cgbo: 0.0,
            rsh: 0.0,
            cj: 0.0,
            mj: 0.5,
            cjsw: 0.0,
            mjsw: 0.5,
            js: 1.0e-8, // FIXME
            tox: 1.0e-7,
            nsub: 0.0,
            nss: 0.0,
            tpg: false,
            ld: 0.0,
            u0: 600.0,
            fc: 0.5,
            kf: 0.0,
            af: 1.0,
            tnom: 27.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Mos1InstanceParams {
    m: f64,
    l: f64,
    w: f64,
    a_d: f64,
    a_s: f64,
    pd: f64,
    ps: f64,
    nrd: f64,
    nrs: f64,
    temp: f64,
    // FIXME: explicitly ignore these
    dtemp: f64,
    off: bool,
    icvds: f64,
    icvgs: f64,
    icvbs: f64,
    ic: f64,
}


impl Default for Mos1InstanceParams {
    fn default() -> Self {
        Mos1InstanceParams {
            m: 0.0,
            l: 1e-6,
            w: 1e-6,
            a_d: 1e-12,
            a_s: 1e-12,
            pd: 1e-6,
            ps: 1e-6,
            nrd: 1.0,
            nrs: 1.0,
            temp: TEMP_REF,
            // FIXME: explicitly ignore these
            off: false,
            dtemp: 0.0,
            icvds: 0.0,
            icvgs: 0.0,
            icvbs: 0.0,
            ic: 0.0,
        }
    }
}

/// Mos1 Internal "Parameters", derived at instance-construction
/// and updated only on changes in temperature
#[derive(Default)]
struct Mos1InternalParams {
    temp: f64,
    vt_t: f64,
    kp_t: f64,
    phi_t: f64,
    beta: f64,
    cox: f64,
    cgs_ov: f64,
    cgd_ov: f64,
    cgb_ov: f64,
    leff: f64,
    isat_bd: f64,
    isat_bs: f64,
}

/// Mos1 DC & Transient Operating Point
#[derive(Default, Clone)]
struct Mos1OpPoint {
    //    id: f64,
//    is: f64,
//    ig: f64,
//    ib: f64,
//    ibd: f64,
//    ibs: f64,
    vgs: f64,
    vgd: f64,
    vgb: f64,
//    vds: f64,
//    vbs: f64,
//    vbd: f64,

//    von: f64,
//    vdsat: f64,
//    sourcevcrit: f64,
//    drainvcrit: f64,
//    rs: f64,
//    sourceconductance: f64,
//    rd: f64,
//    drainconductance: f64,
//
//    gm: f64,
//    gds: f64,
//    gmb: f64,
//    gmbs: f64,
//    gbd: f64,
//    gbs: f64,

    //    cbd: f64,
//    cbs: f64,
    cgs: f64,
    cgd: f64,
    cgb: f64,

//    cqgs: f64,
//    cqgd: f64,
//    cqgb: f64,
//    cqbd: f64,
//    cqbs: f64,
//
//    cbd0: f64,
//    cbdsw0: f64,
//    cbs0: f64,
//    cbssw0: f64,

    qgs: f64,
    qgd: f64,
    qgb: f64,
//    qbd: f64,
//    qbs: f64,
//    pwr: f64,

    icgs: f64,
    icgd: f64,
    icgb: f64,
}

struct Mos1 {
    model: Mos1Model,
    params: Mos1InstanceParams,
    intparams: Mos1InternalParams,
    op: Mos1OpPoint,
    guess: Mos1OpPoint,
    ports: MosTerminals,
    matps: MosMatrixPointers,
}

const KB: f64 = 1.3806226e-23;
const Q: f64 = 1.6021918e-19;
const KB_OVER_Q: f64 = KB / Q;
const TEMP_REF: f64 = 300.15;
const KELVIN_TO_C: f64 = 273.15;
const SIO2_PERMITTIVITY: f64 = 3.9 * 8.854214871e-12;

/// Mosfet Level 1 Instance
impl Mos1 {
    fn new(model: Mos1Model, params: Mos1InstanceParams, ports: MosTerminals) -> Mos1 {
        let intparams = Mos1::derive(&model, &params);
        Mos1 { model, params, intparams, ports, op: Mos1OpPoint::default(), guess: Mos1OpPoint::default(), matps: MosMatrixPointers::default() }
    }
    /// Calculate derived parameters from instance parameters
    fn derive(model: &Mos1Model, inst: &Mos1InstanceParams) -> Mos1InternalParams {
        let temp = inst.temp;
        let vtherm = temp * KB_OVER_Q;
        // FIXME: all temperature dependences
        let phi_t = model.phi;

        let leff = inst.l - model.ld * 2.;

        let cox_per_area = SIO2_PERMITTIVITY / model.tox;
        let cox = cox_per_area * leff * inst.w;

        let kp_t = model.u0 * cox_per_area * 1e-4;
        let beta = kp_t * inst.w / leff;

        let isat_bd = 0.0;
        let isat_bs = 0.0;

        Mos1InternalParams { temp, leff, cox, beta, phi_t, isat_bd, isat_bs, ..Default::default() }
    }
}

impl Component for Mos1 {
    fn create_matrix_elems(&mut self, mat: &mut Matrix) {
        use MosTerm::{G, D, S, B};
        for t1 in [G, D, S, B].iter() {
            for t2 in [G, D, S, B].iter() {
                self.matps[(*t1, *t2)] = make_matrix_elem(mat, self.ports[*t1], self.ports[*t2]);
            }
        }
    }
    /// Load our last guess as the new operating point
    fn tstep(&mut self, _x: &Vec<f64>) {
        self.op = self.guess.clone();
    }
    fn load(&mut self, guess: &Variables, an: &AnalysisInfo) -> Stamps {
        use MosTerm::{G, D, S, B};

        let vg = guess.get(self.ports[G]);
        let vd = guess.get(self.ports[D]);
        let vs = guess.get(self.ports[S]);
        let vb = guess.get(self.ports[B]);

        // Initially factor out polarity of NMOS/PMOS and source/drain swapping
        // All math after this block uses increasing vgs,vds <=> increasing ids,
        // e.g. the polarities typically expressed for NMOS
        let p = match self.model.mos_type { MosType::NMOS => 1.0, MosType::PMOS => -1.0 };
        let vds1 = p * (vd - vs);
        let reversed = vds1 < 0.0;
        let vgs = if reversed { p * (vg - vd) } else { p * (vg - vs) };
        let vds = if reversed { -vds1 } else { vds1 };
        let vsb = if reversed { vd - vb } else { vs - vb };
        let vdb = if reversed { vs - vb } else { vd - vb };

        let von = if vsb > 0.0 {
            self.intparams.vt_t + self.model.gamma * ((self.intparams.phi_t + vsb).sqrt() - self.intparams.phi_t.sqrt())
        } else {
            self.intparams.vt_t // FIXME: body effect for Vsb < 0
        };

        let vov = vgs - von;
        let vdsat = vov.max(0.0);

        let mut ids = 0.0;
        let mut gm = 0.0;
        let mut gds = 0.0;
        let mut gmbs = 0.0;
        if vov <= 0.0 { // Cutoff
            // Already set
        } else {
            if vds >= vov { // Sat
                ids = self.intparams.beta / 2.0 * vov.powi(2) * (1.0 + self.model.lambda * vds);
                gm = self.intparams.beta * vov * (1.0 + self.model.lambda * vds);
                gds = self.model.lambda * self.intparams.beta / 2.0 * vov.powi(2);
            } else { // Triode
                ids = self.intparams.beta * (vov * vds - vds.powi(2) / 2.0) * (1.0 + self.model.lambda * vds);
                gm = self.intparams.beta * vds * (1.0 + self.model.lambda * vds);
                gds = self.intparams.beta * ((vov - vds) * (1.0 + self.model.lambda * vds) + self.model.lambda * ((vov * vds) - vds.powi(2) / 2.0));
            }
            gmbs = if self.intparams.phi_t + vsb > 0.0 {
                gm * self.model.gamma / 2.0 / (self.intparams.phi_t + vsb).sqrt()
            } else { 0.0 };
        }

        let mut cgs1: f64;
        let mut cgd1: f64;
        let mut cgb1: f64;
        let cox = self.intparams.cox;

        if vov <= -self.intparams.phi_t {
            cgb1 = cox / 2.0;
            cgs1 = 0.0;
            cgd1 = 0.0;
        } else if vov <= -self.intparams.phi_t / 2.0 {
            cgb1 = -vov * cox / (2.0 * self.intparams.phi_t);
            cgs1 = 0.0;
            cgd1 = 0.0;
        } else if vov <= 0.0 {
            cgb1 = -vov * cox / (2.0 * self.intparams.phi_t);
            cgs1 = vov * cox / (1.5 * self.intparams.phi_t) + cox / 3.0;
            cgd1 = 0.0;
        } else if vdsat <= vds {
            cgs1 = cox / 3.0;
            cgd1 = 0.0;
            cgb1 = 0.0;
        } else {
            let vddif = 2.0 * vdsat - vds;
            let vddif1 = vdsat - vds;
            let vddif2 = vddif * vddif;
            cgd1 = cox * (1.0 - vdsat * vdsat / vddif2) / 3.0;
            cgs1 = cox * (1.0 - vddif1 * vddif1 / vddif2) / 3.0;
            cgb1 = 0.0;
        }

        // Now start incorporating past history
        // FIXME: gotta sort out swaps in polarity between time-points
        // FIXME: this isnt quite right as we move from OP into first TRAN point. hacking that for now
        let cgs2 = if self.op.cgs == 0.0 {
            cgs1
        } else {
            self.op.cgs
        };
        let cgs = cgs1 + cgs2 + self.intparams.cgs_ov;
        let cgd = cgd1 + self.op.cgd + self.intparams.cgd_ov;
        let cgb = cgb1 + self.op.cgb + self.intparams.cgb_ov;

        let vgd = vgs - vds;
        let vgb = vgs + vsb;

        let vgs1 = self.op.vgs;
        let vgd1 = self.op.vgd;
        let vgb1 = self.op.vgb;
        let qgs = (vgs - vgs1) * cgs + self.op.qgs;
        let qgd = (vgd - vgd1) * cgd + self.op.qgd;
        let qgb = (vgb - vgb1) * cgb + self.op.qgb;


//        let (mut icgs, mut icgd, icgb) = (0.0, 0.0, 0.0);
//        let (mut gcgs, mut gcgd, gcgb) = (0.0, 0.0, 0.0);
//        let (mut _rhsgs, mut _rhsgd, _rhsgb) = (0.0, 0.0, 0.0);

        let (gcgs, icgs, rhsgs) = if let AnalysisInfo::TRAN(_d) = an {
            num_integ(qgs, self.op.qgs, cgs, vgs, Some(self.op.icgs), THE_TIMESTEP)
        } else {
            (0.0, 0.0, 0.0)
        };
        let (gcgd, icgd, rhsgd) = if let AnalysisInfo::TRAN(_d) = an {
            num_integ(qgd, self.op.qgd, cgd, vgd, Some(self.op.icgd), THE_TIMESTEP)
        } else {
            (0.0, 0.0, 0.0)
        };
        let (gcgb, icgb, rhsgb) = if let AnalysisInfo::TRAN(_d) = an {
            num_integ(qgd, self.op.qgd, cgd, vgd, Some(self.op.icgd), THE_TIMESTEP)
        } else {
            (0.0, 0.0, 0.0)
        };
//            icgs += -gcgs * vgs - 2.0 * qgs / THE_TIMESTEP;
//            icgd += -gcgd * vgd - 2.0 * qgd / THE_TIMESTEP;
//            icgb += -gcgd * vgb - 2.0 * qgb / THE_TIMESTEP;

        // FIXME: bulk junction diodes
        let cbs = 0.0;
        let cbd = 0.0;
        let gbd = 1e-9;
        let gbs = 1e-9;
        let gcbd = 0.0;
        let gcbs = 0.0;
        let gcbg = 0.0;
        // FIXME: terminal resistances
        let grd = 0.0;
        let grs = 0.0;

        // Store as our op point for next time
        self.guess = Mos1OpPoint { vgs, vgd, vgb, cgs: cgs1, cgd: cgd1, cgb: cgb1, qgs, qgd, qgb, icgs, icgd, icgb };

        // Bulk junction caps RHS adjustment
        let ceqbs = 0.0; //p * (cbs + gbs * vsb);
        let ceqbd = 0.0; //p * (cbd + gbd * vdb);
        let irhs = ids - gm * vgs - gds * vds;

        // Sort out which are the "reported" drain and source terminals (sr, dr)
        let (sr, sx, dr, dx) = if !reversed { (S, S, D, D) } else { (D, D, S, S) };
        return Stamps {
            g: vec![
                (self.matps[(dr, dr)], gds + grd + gbd + gcgd),
                (self.matps[(sr, sr)], gm + gds + grs + gbs + gmbs + gcgs),
                (self.matps[(dr, sr)], -gm - gds - gmbs),
                (self.matps[(sr, dr)], -gds),
                (self.matps[(dr, G)], gm - gcgd),
                (self.matps[(sr, G)], -gm - gcgs),
                (self.matps[(G, G)], (gcgd + gcgs + gcgb)),
                (self.matps[(B, B)], (gbd + gbs + gcgb)),
                (self.matps[(G, B)], -gcgb),
                (self.matps[(G, dr)], -gcgd),
                (self.matps[(G, sr)], -gcgs),
                (self.matps[(B, G)], -gcbg),
                (self.matps[(G, dr)], -gcgd),
                (self.matps[(B, dr)], -gbd),
                (self.matps[(B, sr)], -gbs),
                (self.matps[(dr, B)], -gbd + gmbs),
                (self.matps[(sr, B)], -gbs - gmbs),
                (self.matps[(dx, dr)], -grd),
                (self.matps[(dr, dx)], -grd),
                (self.matps[(dx, dx)], grd),
                (self.matps[(sx, sr)], -grs),
                (self.matps[(sr, sx)], -grs),
                (self.matps[(sx, sx)], grs),
            ],
            b: vec![
                (self.ports[dr], -p * irhs + ceqbd + p * rhsgd),
                (self.ports[sr], p * irhs + ceqbs + p * rhsgs),
                (self.ports[G], -p * (rhsgs + rhsgb + rhsgd)),
                (self.ports[B], -(ceqbs + ceqbd - p * rhsgb)),
            ],
        };
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
    fn new(ports: MosTerminals, vth: f64, beta: f64, lam: f64, polarity: bool) -> Mos {
        Mos {
            vth,
            beta,
            lam,
            polarity,
            ports,
            matps: MosMatrixPointers([[None; 4]; 4]),
        }
    }
}

impl Component for Mos {
    fn create_matrix_elems(&mut self, mat: &mut Matrix) {
        use MosTerm::{G, D, S};
        let matps = [(D, D), (S, S), (D, S), (S, D), (D, G), (S, G)];
        for (t1, t2) in matps.iter() {
            self.matps[(*t1, *t2)] = make_matrix_elem(mat, self.ports[*t1], self.ports[*t2]);
        }
    }
    fn load(&mut self, guess: &Variables, an: &AnalysisInfo) -> Stamps {
        use MosTerm::{G, D, S, B};

        let vg = guess.get(self.ports[G]);
        let vd = guess.get(self.ports[D]);
        let vs = guess.get(self.ports[S]);
        let vb = guess.get(self.ports[B]);
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
        let irhs = ids - gm * vgs - gds * vds;
        return Stamps {
            g: vec![
                (self.matps[(dr, dr)], gds),
                (self.matps[(sr, sr)], (gm + gds)),
                (self.matps[(dr, sr)], -(gm + gds)),
                (self.matps[(sr, dr)], -gds),
                (self.matps[(dr, G)], gm),
                (self.matps[(sr, G)], -gm),
            ],
            b: vec![
                (self.ports[dr], -p * irhs),
                (self.ports[sr], p * irhs),
            ],
        };
    }
}

#[derive(Default)]
struct Diode {
    isat: f64,
    vt: f64,
    p: Option<VarIndex>,
    n: Option<VarIndex>,
    pp: Option<Eindex>,
    nn: Option<Eindex>,
    pn: Option<Eindex>,
    np: Option<Eindex>,
}

impl Component for Diode {
    fn create_matrix_elems(&mut self, mat: &mut Matrix) {
        self.pp = make_matrix_elem(mat, self.p, self.p);
        self.pn = make_matrix_elem(mat, self.p, self.n);
        self.np = make_matrix_elem(mat, self.n, self.p);
        self.nn = make_matrix_elem(mat, self.n, self.n);
    }
    fn load(&mut self, guess: &Variables, an: &AnalysisInfo) -> Stamps {
        let vp = guess.get(self.p);
        let vn = guess.get(self.n);
        let vd = (vp - vn).max(-1.5).min(1.5);
        let i = self.isat * ((vd / self.vt).exp() - 1.0);
        let di_dv = (self.isat / self.vt) * (vd / self.vt).exp();
        let irhs = i - vd * di_dv;

        return Stamps {
            g: vec![
                (self.pp, di_dv),
                (self.nn, di_dv),
                (self.pn, -di_dv),
                (self.np, -di_dv)
            ],
            b: vec![
                (self.p, -irhs),
                (self.n, irhs)
            ],
        };
    }
}

#[derive(Default)]
struct Isrc {
    i: f64,
    p: Option<VarIndex>,
    n: Option<VarIndex>,
}

impl Component for Isrc {
    fn create_matrix_elems(&mut self, _mat: &mut Matrix) {}
    fn load(&mut self, _guess: &Variables, _an: &AnalysisInfo) -> Stamps {
        return Stamps {
            g: vec![],
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
    b: Vec<(Option<VarIndex>, f64)>,
}

impl Stamps {
    fn new() -> Stamps {
        Stamps { g: vec![], b: vec![] }
    }
}

#[derive(Debug, Clone, Copy)]
enum VarKind { V = 0, I }

#[derive(Debug, Clone, Copy)]
struct VarIndex(usize);


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
    fn from(node: &NodeRef) -> Self { (*node).into() }
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

struct Variables {
    kinds: Vec<VarKind>,
    values: Vec<f64>,
}

impl Variables {
    fn all_v(len: usize) -> Variables {
        Variables { kinds: vec![VarKind::V; len], values: vec![0.0; len] }
    }
    fn add(&mut self, kind: VarKind) -> VarIndex {
        self.kinds.push(kind);
        self.values.push(0.0);
        return VarIndex(self.kinds.len() - 1);
    }
    fn get(&self, i: Option<VarIndex>) -> f64 {
        match i {
            None => 0.0,
            Some(ii) => self.values[ii.0],
        }
    }
    fn len(&self) -> usize { self.kinds.len() }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum AnalysisMode { OP, DC, TRAN, AC }

struct Solver {
    comps: Vec<Box<dyn Component>>,
    vars: Variables,
    mat: Matrix,
    rhs: Vec<f64>,
    history: Vec<Vec<f64>>,
    an_mode: AnalysisMode,
}

impl Solver {
    fn add_comp(&mut self, comp: &CompParse) {
        match comp {
            CompParse::R(g, p, n) => {
                let r = Resistor::new(*g, p.into(), n.into());
                self.comps.push(Box::new(r));
            }
            CompParse::C(c, p, n) => {
                let c = Capacitor::new(*c, p.into(), n.into());
                self.comps.push(Box::new(c));
            }
            CompParse::I(i, p, n) => {
                let i = Isrc { i: *i, p: p.into(), n: n.into() };
                self.comps.push(Box::new(i));
            }
            CompParse::D(isat, vt, p, n) => {
                let c = Diode {
                    isat: *isat,
                    vt: *vt,
                    p: p.into(),
                    n: n.into(),
                    ..Default::default()
                };
                self.comps.push(Box::new(c));
            }
            CompParse::V(v, p, n) => {
                let ivar = self.vars.add(VarKind::I);
                let v = Vsrc::new(*v, p.into(), n.into(), ivar);
                self.comps.push(Box::new(v));
            }
            CompParse::Mos(pol, g, d, s, b) => {
//                let dp = self.vars.add(VarKind::V);
//                let sp = self.vars.add(VarKind::V);
                let x = Mos::new(
                    [g, d, s, b].into(),
//                    dp.into(), ds.into(),
                    0.25, 50e-3, 3e-3, *pol,
                );
                self.comps.push(Box::new(x));
            }
            CompParse::Mos1(model, params, g, d, s, b) => {
                let x = Mos1::new(model.clone(), params.clone(), [g, d, s, b].into());
                self.comps.push(Box::new(x));
            }
        }
    }
    fn new(ckt: CktParse) -> Solver {
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
    let mut s = Solver::new(ckt);
    return s.solve(&AnalysisInfo::OP);
}

struct Tran {
    solver: Solver,
    tstop: usize,
    vic: Vec<usize>,
    ric: Vec<usize>,
}

enum AnalysisInfo<'a> {
    OP,
    TRAN(&'a Vec<f64>),
}

impl Tran {
    fn new(ckt: CktParse) -> Tran {
        let solver = Solver::new(ckt);
        return Tran {
            solver,
            tstop: 5000,
            vic: vec![],
            ric: vec![],
        };
    }
    fn ic(&mut self, n: NodeRef, val: f64) {
        let fnode = self.solver.vars.add(VarKind::V);
        let ivar = self.solver.vars.add(VarKind::I);

        let mut r = Resistor::new(1.0, Some(fnode), n.into());
        r.create_matrix_elems(&mut self.solver.mat);
        self.solver.comps.push(Box::new(r));
        self.ric.push(self.solver.comps.len() - 1);
        let mut v = Vsrc::new(val, Some(fnode), None, ivar);
        v.create_matrix_elems(&mut self.solver.mat);
        self.solver.comps.push(Box::new(v));
        self.vic.push(self.solver.comps.len() - 1);
    }
    fn solve(&mut self) -> SpResult<Vec<Vec<f64>>> {
        use std::thread;
        use std::sync::mpsc;

        enum IoWriterMessage { STOP, DATA(Vec<f64>) }
        enum IoWriterResponse { OK, RESULT(Vec<Vec<f64>>) }
        let (tx, rx) = mpsc::channel::<IoWriterMessage>();
        let (tx2, rx2) = mpsc::channel::<IoWriterResponse>();

        let t = thread::spawn(move || {
            use std::fs::File;
            use serde::ser::{SerializeSeq, Serializer};

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
        let tpoint = match tsoln {
            Ok(x) => x,
            Err(e) => {
                println!("Failed to find initial solution");
                return Err(e);
            }
        };
        // Update initial-condition sources and resistances
        for c in self.vic.iter() {
            self.solver.comps[*c].update(0.0);
        }
        for c in self.ric.iter() {
            self.solver.comps[*c].update(1e-9);
        }
        for c in self.solver.comps.iter_mut() {
            c.tstep(&tpoint);
        }
        let mut state = tpoint.clone();
        tx.send(IoWriterMessage::DATA(tpoint)).unwrap();

        self.solver.an_mode = AnalysisMode::TRAN;
        for _t in 1..self.tstop {
            let tsoln = self.solver.solve(&AnalysisInfo::TRAN(&state));
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
            tx.send(IoWriterMessage::DATA(tpoint)).unwrap();
        }
        tx.send(IoWriterMessage::STOP).unwrap();
        for msg in rx2 {
            match msg {
                IoWriterResponse::OK => { continue; }
                IoWriterResponse::RESULT(res) => {
                    t.join().unwrap();
                    return Ok(res);
                }
            }
        }
        Err("Tran Results Failure")
    }
}

pub struct TranOptions {
    tstep: f64,
    tstop: f64,
}

pub fn tran(ckt: CktParse) -> SpResult<Vec<Vec<f64>>> {
    return Tran::new(ckt).solve();
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
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(1e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Gnd),
            ],
        };
        Ok(())
    }

    #[test]
    fn test_dcop1() -> TestResult {
        let ckt = CktParse {
            nodes: 1,
            comps: vec![CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Gnd)],
        };

        let soln = dcop(ckt)?;
        assert_eq!(soln, vec![0.0]);
        Ok(())
    }

    #[test]
    fn test_dcop2() -> TestResult {
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(1e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Gnd),
            ],
        };

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
        assert!((soln[0] - 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_diode_nmos_tran() -> TestResult {
        // Diode NMOS Tran
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(5e-3, Num(0), Gnd),
                CompParse::Mos(true, Num(0), Num(0), Gnd, Gnd),
                CompParse::R(1e-12, Num(0), Gnd), // "gmin"
            ],
        };

        let soln = tran(ckt)?;
        for point in soln.iter() {
            assert!((point[0] - 0.697).abs() < 1e-3);
        }
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

        let soln = dcop(ckt)?;
        assert!((soln[0] - 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_diode_pmos_dcop() -> TestResult {
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

        let soln = dcop(ckt)?;
        assert!((soln[0] + 0.697).abs() < 1e-3);
        Ok(())
    }


    #[test]
    fn test_diode_pmos_tran() -> TestResult {
        // Diode PMOS Tran
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(-5e-3, Num(0), Gnd),
                CompParse::Mos(false, Num(0), Num(0), Gnd, Gnd),
                CompParse::R(1e-12, Num(0), Gnd), // "gmin"
            ],
        };

        let soln = tran(ckt)?;
        for point in soln.iter() {
            assert!((point[0] + 0.697).abs() < 1e-3);
        }
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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

        let soln = dcop(ckt)?;
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
        for point in soln.into_iter() {
            assert(point).eq(vec![1.0, 1.0, 0.0])?;
        }
        Ok(())
    }

    #[test]
    fn test_tran2() -> TestResult {
        // I-C Integrator, with Initial Condition
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
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

    #[test]
    fn test_tran3() -> TestResult {
        // Ring Oscillator
        use NodeRef::{Gnd, Num};
        let c = 1e-10;
        let ckt = CktParse {
            nodes: 4,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::R(1e-3, Num(0), Gnd),
                CompParse::Mos(false, Num(3), Num(1), Num(0), Num(0)),
                CompParse::Mos(true, Num(3), Num(1), Gnd, Gnd),
                CompParse::R(1e-5, Num(1), Gnd),
                CompParse::C(c, Num(1), Gnd),
                CompParse::Mos(false, Num(1), Num(2), Num(0), Num(0)),
                CompParse::Mos(true, Num(1), Num(2), Gnd, Gnd),
                CompParse::R(1e-5, Num(2), Gnd),
                CompParse::C(c, Num(2), Gnd),
                CompParse::Mos(false, Num(2), Num(3), Num(0), Num(0)),
                CompParse::Mos(true, Num(2), Num(3), Gnd, Gnd),
                CompParse::R(1e-5, Num(3), Gnd),
                CompParse::C(c, Num(3), Gnd),
            ],
        };

        let mut tran = Tran::new(ckt);
        tran.ic(Num(1), 0.0);
        let soln = tran.solve()?;
        // FIXME: dream up some checks
        Ok(())
    }

    #[test]
    fn test_mos1_op() -> TestResult {
        let model = Mos1Model::default();
        let params = Mos1InstanceParams::default();
        use NodeRef::{Num, Gnd};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::Mos1(model, params, Num(0), Num(0), Gnd, Gnd),
                CompParse::V(1.0, Num(0), Gnd),
            ],
        };
        let soln = dcop(ckt)?;
        assert(soln[0]).eq(1.0)?;
        assert(soln[1]).lt(0.0)?;
        assert(soln[1]).gt(-1e-3)?;
        Ok(())
    }

    #[test]
    fn test_mos1_tran() -> TestResult {
        let model = Mos1Model::default();
        let params = Mos1InstanceParams::default();
        use NodeRef::{Num, Gnd};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::Mos1(model, params, Num(0), Num(0), Gnd, Gnd),
                CompParse::V(1.0, Num(0), Gnd),
            ],
        };
        let soln = tran(ckt)?;
        for k in 1..soln.len() {
            assert(soln[k][0]).eq(1.0)?;
            assert(soln[k][1]).lt(0.0)?;
            assert(soln[k][1]).gt(-1e-3)?;
        }
        Ok(())
    }

    #[test]
    fn test_mos1_inv_dcop() -> TestResult {
        // Mos1 Inverter DCOP

        let nmos = Mos1Model::default();
        let pmos = Mos1Model { mos_type: MosType::PMOS, ..Mos1Model::default() };
        let params = Mos1InstanceParams::default();
        use NodeRef::{Gnd, Num};
        let c = 1e-10;
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::Mos1(pmos.clone(), params, Gnd, Num(1), Num(0), Num(0)),
                CompParse::Mos1(nmos.clone(), params, Gnd, Num(1), Gnd, Gnd),
                CompParse::R(1e-4, Num(1), Gnd),
            ],
        };
        let soln = dcop(ckt)?;
        Ok(())
    }

    #[test]
    fn test_mos1_ro_dcop() -> TestResult {
        // Mos1 Ring Oscillator Dc Op

        let nmos = Mos1Model::default();
        let pmos = Mos1Model { mos_type: MosType::PMOS, ..Mos1Model::default() };
        let params = Mos1InstanceParams::default();
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 4,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::Mos1(pmos.clone(), params, Num(3), Num(1), Num(0), Num(0)),
                CompParse::Mos1(nmos.clone(), params, Num(3), Num(1), Gnd, Gnd),
                CompParse::Mos1(pmos.clone(), params, Num(1), Num(2), Num(0), Num(0)),
                CompParse::Mos1(nmos.clone(), params, Num(1), Num(2), Gnd, Gnd),
                CompParse::Mos1(pmos.clone(), params, Num(2), Num(3), Num(0), Num(0)),
                CompParse::Mos1(nmos.clone(), params, Num(2), Num(3), Gnd, Gnd),
            ],
        };

        let soln = dcop(ckt)?;
        assert(soln[0]).eq(1.0)?;
        for k in 1..3 {
            assert(soln[k]).gt(0.45)?;
            assert(soln[k]).lt(0.55)?;
        }
        Ok(())
    }

    #[test]
    fn test_mos1_ro_tran() -> TestResult {
        // Mos1 Ring Oscillator Tran

        let nmos = Mos1Model::default();
        let pmos = Mos1Model { mos_type: MosType::PMOS, ..Mos1Model::default() };
        let params = Mos1InstanceParams::default();
        use NodeRef::{Gnd, Num};
        let c = 1e-12;
        let ckt = CktParse {
            nodes: 4,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::Mos1(pmos.clone(), params, Num(3), Num(1), Num(0), Num(0)),
                CompParse::Mos1(nmos.clone(), params, Num(3), Num(1), Gnd, Gnd),
                CompParse::C(c, Num(1), Gnd),
                CompParse::R(1e-9, Num(1), Gnd),
                CompParse::Mos1(pmos.clone(), params, Num(1), Num(2), Num(0), Num(0)),
                CompParse::Mos1(nmos.clone(), params, Num(1), Num(2), Gnd, Gnd),
                CompParse::C(c, Num(2), Gnd),
                CompParse::R(1e-9, Num(2), Gnd),
                CompParse::Mos1(pmos.clone(), params, Num(2), Num(3), Num(0), Num(0)),
                CompParse::Mos1(nmos.clone(), params, Num(2), Num(3), Gnd, Gnd),
                CompParse::C(c, Num(3), Gnd),
                CompParse::R(1e-9, Num(3), Gnd),
            ],
        };

        let mut tran = Tran::new(ckt);
        tran.ic(Num(1), 0.0);
        let soln = tran.solve()?;
        // FIXME: dream up some checks
        Ok(())
    }
}

