use enum_dispatch::enum_dispatch;
use std::convert::From;
use std::ops::{Index, IndexMut};
use num::Complex;

use super::SpNum;
use super::analysis::{AnalysisInfo, Stamps, VarIndex, Variables};
use super::proto::NodeRef;
use super::sparse21::{Eindex, Matrix};

#[enum_dispatch]
pub enum ComponentSolver {
    Vsrc,
    Isrc,
    Capacitor,
    Resistor,
    Diode,
    Mos0,
    Mos1,
}

#[enum_dispatch(ComponentSolver)]
pub trait Component {
    fn tstep(&mut self, _guess: &Vec<f64>) {}

    fn update(&mut self, _val: f64) {}
    // FIXME: prob not for every Component

    fn load_ac(&mut self, guess: &Variables<Complex<f64>>, an: &AnalysisInfo) -> Stamps<Complex<f64>> { Stamps::<Complex<f64>>::new()}
    fn load(&mut self, guess: &Variables<f64>, an: &AnalysisInfo) -> Stamps<f64>;
    fn create_matrix_elems<T: SpNum>(&mut self, mat: &mut Matrix<T>);
}

pub struct Vsrc {
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
    pub fn new(v: f64, p: Option<VarIndex>, n: Option<VarIndex>, ivar: VarIndex) -> Vsrc {
        Vsrc {
            v,
            p,
            n,
            ivar,
            pi: None,
            ip: None,
            ni: None,
            in_: None,
        }
    }
}

impl Component for Vsrc {
    fn update(&mut self, val: f64) {
        self.v = val;
    }
    fn create_matrix_elems<T: SpNum>(&mut self, mat: &mut Matrix<T>) {
        self.pi = make_matrix_elem(mat, self.p, Some(self.ivar));
        self.ip = make_matrix_elem(mat, Some(self.ivar), self.p);
        self.ni = make_matrix_elem(mat, self.n, Some(self.ivar));
        self.in_ = make_matrix_elem(mat, Some(self.ivar), self.n);
    }
    fn load(&mut self, _guess: &Variables<f64>, _an: &AnalysisInfo) -> Stamps<f64> {
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
    fn load_ac(&mut self, _guess: &Variables<Complex<f64>>, _an: &AnalysisInfo) -> Stamps<Complex<f64>> {
        // FIXME: always has ACM=1, for now 
        return Stamps {
            g: vec![
                (self.pi, Complex::new(1.0, 0.0)),
                (self.ip, Complex::new(1.0, 0.0)),
                (self.ni, Complex::new(-1.0, 0.0)),
                (self.in_, Complex::new(-1.0, 0.0)),
            ],
            b: vec![
                (Some(self.ivar), Complex::new(self.v, 0.0))
            ],
        };
    }
}

#[derive(Default)]
pub struct Capacitor {
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
    pub fn new(c: f64, p: Option<VarIndex>, n: Option<VarIndex>) -> Capacitor {
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

impl Component for Capacitor {
    fn create_matrix_elems<T: SpNum>(&mut self, mat: &mut Matrix<T>) {
        self.pp = make_matrix_elem(mat, self.p, self.p);
        self.pn = make_matrix_elem(mat, self.p, self.n);
        self.np = make_matrix_elem(mat, self.n, self.p);
        self.nn = make_matrix_elem(mat, self.n, self.n);
    }
    /// Load our last guess as the new operating point
    fn tstep(&mut self, _x: &Vec<f64>) {
        self.op = self.guess;
    }
    fn load(&mut self, guess: &Variables<f64>, an: &AnalysisInfo) -> Stamps<f64> {
        let vd = guess.get(self.p) - guess.get(self.n);
        let q = self.q(vd);

        match *an {
            AnalysisInfo::OP => {
                // FIXME: calculating this during DCOP, so we copy cleanly afterward
                // Should probably just find a way to calculate it then
                self.guess = CapOpPoint {
                    v: vd,
                    q: q,
                    i: 0.0,
                };
                return Stamps::new();
            }
            AnalysisInfo::TRAN(opts, state) => {
                let (g, i, rhs) = state.integrate(q - self.op.q, self.dq_dv(vd), vd, self.op.i);
                self.guess = CapOpPoint { v: vd, q: q, i: i };

                return Stamps {
                    g: vec![(self.pp, g), (self.nn, g), (self.pn, -g), (self.np, -g)],
                    b: vec![(self.p, -rhs), (self.n, rhs)],
                };
            }
            AnalysisInfo::AC(_o, _s) => panic!("HOW WE GET HERE?!?")
        }
    }
    fn load_ac(&mut self, _guess: &Variables<Complex<f64>>, an: &AnalysisInfo) -> Stamps<Complex<f64>> {
        let an_st = match an {
            AnalysisInfo::AC(opts, state) => state,
            _ => panic!("Invalid AC AnalysisInfo")
        };
        let c = self.dq_dv(0.0);
        return Stamps {
            g: vec![
                (self.pp, Complex::new(0.0, an_st.omega * c)),
                (self.nn, Complex::new(0.0, an_st.omega * c)),
                (self.pn, Complex::new(0.0, -an_st.omega * c)),
                (self.np, Complex::new(0.0, -an_st.omega * c)),
            ],
            b: vec![],
        };
    }
}

#[derive(Clone, Copy)]
enum TwoTerm {
    P = 0,
    N = 1,
}

struct TwoTerminals([Option<VarIndex>; 2]);

impl Index<TwoTerm> for TwoTerminals {
    type Output = Option<VarIndex>;
    fn index(&self, t: TwoTerm) -> &Option<VarIndex> {
        &self.0[t as usize]
    }
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

pub struct Resistor {
    g: f64,
    terms: TwoTerminals,
    matps: TwoTermMatrixPointers,
}

impl Resistor {
    pub fn new(g: f64, p: Option<VarIndex>, n: Option<VarIndex>) -> Resistor {
        Resistor {
            g,
            terms: TwoTerminals([p, n]),
            matps: TwoTermMatrixPointers([[None; 2]; 2]),
        }
    }
}

impl Component for Resistor {
    fn update(&mut self, val: f64) {
        self.g = val;
    }
    fn create_matrix_elems<T: SpNum>(&mut self, mat: &mut Matrix<T>) {
        use TwoTerm::{N, P};
        for l in [P, N].into_iter() {
            for r in [P, N].into_iter() {
                self.matps[(*l, *r)] = make_matrix_elem(mat, self.terms[*l], self.terms[*r]);
            }
        }
    }
    fn load(&mut self, _guess: &Variables<f64>, _an: &AnalysisInfo) -> Stamps<f64> {
        use TwoTerm::{N, P};
        return Stamps {
            g: vec![
                (self.matps[(P, P)], self.g),
                (self.matps[(N, N)], self.g),
                (self.matps[(P, N)], -self.g),
                (self.matps[(N, P)], -self.g),
            ],
            b: vec![],
        };
    }
    fn load_ac(&mut self, _guess: &Variables<Complex<f64>>, _an: &AnalysisInfo) -> Stamps<Complex<f64>> {
        use TwoTerm::{N, P};
        return Stamps {
            g: vec![
                (self.matps[(P, P)], Complex::new(self.g, 0.0)),
                (self.matps[(N, N)], Complex::new(self.g, 0.0)),
                (self.matps[(P, N)], Complex::new(-self.g, 0.0)),
                (self.matps[(N, P)], Complex::new(-self.g, 0.0)),
            ],
            b: vec![],
        };
    }
}

/// Mos Terminals, in SPICE order: g, d, s, b
#[derive(Clone, Copy)]
enum MosTerm {
    G = 0,
    D = 1,
    S = 2,
    B = 3,
}

pub struct MosTerminals([Option<VarIndex>; 4]);

impl Index<MosTerm> for MosTerminals {
    type Output = Option<VarIndex>;
    fn index(&self, t: MosTerm) -> &Option<VarIndex> {
        &self.0[t as usize]
    }
}

impl From<[&NodeRef; 4]> for MosTerminals {
    fn from(n: [&NodeRef; 4]) -> Self {
        return MosTerminals([n[0].into(), n[1].into(), n[2].into(), n[3].into()]);
    }
}

#[derive(Default)]
struct MosMatrixPointers([[Option<Eindex>; 4]; 4]);

impl Index<(MosTerm, MosTerm)> for MosMatrixPointers {
    type Output = Option<Eindex>;
    fn index(&self, ts: (MosTerm, MosTerm)) -> &Option<Eindex> {
        &self.0[ts.0 as usize][ts.1 as usize]
    }
}

impl IndexMut<(MosTerm, MosTerm)> for MosMatrixPointers {
    fn index_mut(&mut self, ts: (MosTerm, MosTerm)) -> &mut Self::Output {
        &mut self.0[ts.0 as usize][ts.1 as usize]
    }
}

#[derive(Clone, Copy)]
pub enum MosType { NMOS, PMOS }

// FIXME: should reference instead of cloning, when we can get it to play nicely with `Box<dyn Comp>`
#[derive(Clone)]
pub struct Mos1Model {
    pub mos_type: MosType,
    pub vt0: f64,
    pub kp: f64,
    pub gamma: f64,
    pub phi: f64,
    pub lambda: f64,
    pub rd: f64,
    pub rs: f64,
    pub cbd: f64,
    pub cbs: f64,
    pub is: f64,
    pub pb: f64,
    pub cgso: f64,
    pub cgdo: f64,
    pub cgbo: f64,
    pub rsh: f64,
    pub cj: f64,
    pub mj: f64,
    pub cjsw: f64,
    pub mjsw: f64,
    pub js: f64,
    pub tox: f64,
    pub ld: f64,
    pub u0: f64,
    pub fc: f64,
    pub nsub: f64,
    pub tpg: bool,
    pub nss: f64,
    pub tnom: f64,
    pub kf: f64,
    pub af: f64,
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
pub struct Mos1InstanceParams {
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
    grd: f64, 
    grs: f64, 
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
       gm: f64,
       gds: f64,
    //    gmb: f64,
       gmbs: f64,
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
    reversed: bool,
}

pub struct Mos1 {
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
    pub fn new(model: Mos1Model, params: Mos1InstanceParams, ports: MosTerminals) -> Mos1 {
        let intparams = Mos1::derive(&model, &params);
        Mos1 {
            model,
            params,
            intparams,
            ports,
            op: Mos1OpPoint::default(),
            guess: Mos1OpPoint::default(),
            matps: MosMatrixPointers::default(),
        }
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

        Mos1InternalParams {
            temp,
            leff,
            cox,
            beta,
            phi_t,
            isat_bd,
            isat_bs,
            ..Default::default()
        }
    }
}

impl Component for Mos1 {
    fn create_matrix_elems<T: SpNum>(&mut self, mat: &mut Matrix<T>) {
        use MosTerm::{B, D, G, S};
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
    fn load(&mut self, guess: &Variables<f64>, an: &AnalysisInfo) -> Stamps<f64> {
        use MosTerm::{B, D, G, S};

        let vg = guess.get(self.ports[G]);
        let vd = guess.get(self.ports[D]);
        let vs = guess.get(self.ports[S]);
        let vb = guess.get(self.ports[B]);

        // Initially factor out polarity of NMOS/PMOS and source/drain swapping
        // All math after this block uses increasing vgs,vds <=> increasing ids,
        // i.e. the polarities typically expressed for NMOS
        let p = match self.model.mos_type {
            MosType::NMOS =>  1.0,
            MosType::PMOS => -1.0,
        };
        let vds1 = p * (vd - vs);
        let reversed = vds1 < 0.0;
        let vgs = if reversed {
            p * (vg - vd)
        } else {
            p * (vg - vs)
        };
        let vds = if reversed { -vds1 } else { vds1 };
        let vsb = if reversed { vd - vb } else { vs - vb };
        let vdb = if reversed { vs - vb } else { vd - vb };

        let von = if vsb > 0.0 {
            self.intparams.vt_t
                + self.model.gamma
                    * ((self.intparams.phi_t + vsb).sqrt() - self.intparams.phi_t.sqrt())
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
            if vds >= vov {
                // Sat
                ids = self.intparams.beta / 2.0 * vov.powi(2) * (1.0 + self.model.lambda * vds);
                gm = self.intparams.beta * vov * (1.0 + self.model.lambda * vds);
                gds = self.model.lambda * self.intparams.beta / 2.0 * vov.powi(2);
            } else {
                // Triode
                ids = self.intparams.beta
                    * (vov * vds - vds.powi(2) / 2.0)
                    * (1.0 + self.model.lambda * vds);
                gm = self.intparams.beta * vds * (1.0 + self.model.lambda * vds);
                gds = self.intparams.beta
                    * ((vov - vds) * (1.0 + self.model.lambda * vds)
                        + self.model.lambda * ((vov * vds) - vds.powi(2) / 2.0));
            }
            gmbs = if self.intparams.phi_t + vsb > 0.0 {
                gm * self.model.gamma / 2.0 / (self.intparams.phi_t + vsb).sqrt()
            } else {
                0.0
            };
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

        let (gcgs, icgs, rhsgs) = if let AnalysisInfo::TRAN(_, state) = an {
            state.integrate(qgs - self.op.qgs, cgs, vgs, self.op.icgs)
        } else {
            (0.0, 0.0, 0.0)
        };
        let (gcgd, icgd, rhsgd) = if let AnalysisInfo::TRAN(_, state) = an {
            state.integrate(qgd - self.op.qgd, cgd, vgd, self.op.icgd)
        } else {
            (0.0, 0.0, 0.0)
        };
        let (gcgb, icgb, rhsgb) = if let AnalysisInfo::TRAN(_, state) = an {
            state.integrate(qgd - self.op.qgd, cgd, vgd, self.op.icgd)
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

        // Store as our op point for next time
        self.guess = Mos1OpPoint {
            vgs,
            vgd,
            vgb,
            cgs: cgs1,
            cgd: cgd1,
            cgb: cgb1,
            qgs,
            qgd,
            qgb,
            icgs,
            icgd,
            icgb,
            gm,
            gds,
            gmbs,
            reversed
        };

        // Bulk junction caps RHS adjustment
        let ceqbs = 0.0; //p * (cbs + gbs * vsb);
        let ceqbd = 0.0; //p * (cbd + gbd * vdb);
        let irhs = ids - gm * vgs - gds * vds;

        // Sort out which are the "reported" drain and source terminals (sr, dr)
        // FIXME: this also needs the "prime" vs "external" source & drains 
        let (sr, sx, dr, dx) = if !reversed {
            (S, S, D, D)
        } else {
            (D, D, S, S)
        };
        // Include our terminal resistances
        let grd = self.intparams.grd;
        let grs = self.intparams.grs;
        // And finally send back our matrix contributions 
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
    fn load_ac(&mut self, _guess: &Variables<Complex<f64>>, an: &AnalysisInfo) -> Stamps<Complex<f64>> {
        // Grab the frequency-variable from our analysis
        let omega = match an {
            AnalysisInfo::AC(opts, state) => state.omega,
            _ => panic!("Invalid AC AnalysisInfo")
        };
        use MosTerm::{B, D, G, S};

        // Short-hand the conductances from our op-point. 
        // (Rustc should be smart enough not to copy these.)
        let (gm, gds, gmbs) = (self.op.gm, self.op.gds, self.op.gmbs);

        // FIXME: bulk junction diodes
        let cbs = 0.0;
        let cbd = 0.0;
        let gbd = 1e-9;
        let gbs = 1e-9;
        let gcbd = 0.0;
        let gcbs = 0.0;
        let gcbg = 0.0;
        
        // FIXME: cap impedances
        let gcgs = 0.0;
        let gcgd = 0.0;
        let gcgb = 0.0;

        // Sort out which are the "reported" drain and source terminals (sr, dr)
        // FIXME: this also needs the "prime" vs "external" source & drains 
        let (sr, sx, dr, dx) = if !self.op.reversed {
            (S, S, D, D)
        } else {
            (D, D, S, S)
        };
        // Include our terminal resistances
        let grd = self.intparams.grd;
        let grs = self.intparams.grs;
        // And finally, send back our AC-matrix contributions
        return Stamps {
            g: vec![
                (self.matps[(dr, dr)], Complex::new(0.0, gds + grd + gbd + gcgd)),
                (self.matps[(sr, sr)], Complex::new(0.0, gm + gds + grs + gbs + gmbs + gcgs)),
                (self.matps[(dr, sr)], Complex::new(0.0, -gm - gds - gmbs)),
                (self.matps[(sr, dr)], Complex::new(0.0, -gds)),
                (self.matps[(dr, G)], Complex::new(0.0, gm - gcgd)),
                (self.matps[(sr, G)], Complex::new(0.0, -gm - gcgs)),
                (self.matps[(G, G)], Complex::new(0.0, (gcgd + gcgs + gcgb))),
                (self.matps[(B, B)], Complex::new(0.0, (gbd + gbs + gcgb))),
                (self.matps[(G, B)], Complex::new(0.0, -gcgb)),
                (self.matps[(G, dr)], Complex::new(0.0, -gcgd)),
                (self.matps[(G, sr)], Complex::new(0.0, -gcgs)),
                (self.matps[(B, G)], Complex::new(0.0, -gcbg)),
                (self.matps[(G, dr)], Complex::new(0.0, -gcgd)),
                (self.matps[(B, dr)], Complex::new(0.0, -gbd)),
                (self.matps[(B, sr)], Complex::new(0.0, -gbs)),
                (self.matps[(dr, B)], Complex::new(0.0, -gbd + gmbs)),
                (self.matps[(sr, B)], Complex::new(0.0, -gbs - gmbs)),
                (self.matps[(dx, dr)], Complex::new(0.0, -grd)),
                (self.matps[(dr, dx)], Complex::new(0.0, -grd)),
                (self.matps[(dx, dx)], Complex::new(0.0, grd)),
                (self.matps[(sx, sr)], Complex::new(0.0, -grs)),
                (self.matps[(sr, sx)], Complex::new(0.0, -grs)),
                (self.matps[(sx, sx)], Complex::new(0.0, grs)),
            ],
            b: vec![],
        };
    }
}

pub struct Mos0Params {
    mos_type: MosType,
    vth: f64,
    beta: f64,
    lam: f64,
}

impl Default for Mos0Params {
    fn default() -> Self {
        Mos0Params {
            mos_type: MosType::NMOS,
            vth: 0.25,
            beta: 50e-3,
            lam: 3e-3,
        }
    }
}

pub struct Mos0 {
    params: Mos0Params,
    ports: MosTerminals,
    matps: MosMatrixPointers,
}

impl Mos0 {
    pub fn new(ports: MosTerminals, mos_type: MosType) -> Self {
        Mos0 {
            params: Mos0Params {
                mos_type: mos_type,
                ..Mos0Params::default()
            },
            ports,
            matps: MosMatrixPointers([[None; 4]; 4]),
        }
    }
}

impl Component for Mos0 {
    fn create_matrix_elems<T: SpNum>(&mut self, mat: &mut Matrix<T>) {
        use MosTerm::{D, G, S};
        let matps = [(D, D), (S, S), (D, S), (S, D), (D, G), (S, G)];
        for (t1, t2) in matps.iter() {
            self.matps[(*t1, *t2)] = make_matrix_elem(mat, self.ports[*t1], self.ports[*t2]);
        }
    }
    fn load(&mut self, guess: &Variables<f64>, an: &AnalysisInfo) -> Stamps<f64> {
        use MosTerm::{B, D, G, S};

        let vg = guess.get(self.ports[G]);
        let vd = guess.get(self.ports[D]);
        let vs = guess.get(self.ports[S]);
        let vb = guess.get(self.ports[B]);
        println!("vg={} vd={} vs={} vb={}", vg, vd, vs, vb);

        let p = match self.params.mos_type {
            MosType::NMOS => 1.0,
            MosType::PMOS => -1.0,
        };    
        let vds1 = p * (vd - vs);
        let reversed = vds1 < 0.0;
        let vgs = if reversed {
            p * (vg - vd)
        } else {
            p * (vg - vs)
        };
        let vds = if reversed { -vds1 } else { vds1 };
        let vov = vgs - self.params.vth;

        let mut ids = 0.0;
        let mut gm = 0.0;
        let mut gds = 0.0;
        let lam = self.params.lam;
        let beta = self.params.beta;
        if vov <= 0.0 {
            // Cutoff
            // Already set
            println!("CUTOFF: vgs={} vds={}", vgs, vds);
        } else if vds >= vov {
            // Sat
            ids = beta / 2.0 * vov.powi(2) * (1.0 + lam * vds);
            gm = beta * vov * (1.0 + lam * vds);
            gds = lam * beta / 2.0 * vov.powi(2);
            println!(
                "SAT: vgs={} vds={}, ids={}, gm={}, gds={}",
                vgs, vds, ids, gm, gds
            );
        } else {
            //Triode
            ids = beta * (vov * vds - vds.powi(2) / 2.0) * (1.0 + lam * vds);
            gm = beta * vds * (1.0 + lam * vds);
            gds = beta
                * ((vov - vds) * (1.0 + lam * vds)
                    + lam * ((vov * vds) - vds.powi(2) / 2.0));
            println!(
                "LIN: vgs={} vds={}, ids={}, gm={}, gds={}",
                vgs, vds, ids, gm, gds
            );
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
            b: vec![(self.ports[dr], -p * irhs), (self.ports[sr], p * irhs)],
        };
    }
}

#[derive(Default)]
pub struct Diode {
    isat: f64,
    vt: f64,
    p: Option<VarIndex>,
    n: Option<VarIndex>,
    pp: Option<Eindex>,
    nn: Option<Eindex>,
    pn: Option<Eindex>,
    np: Option<Eindex>,
}

impl Diode {
    pub fn new(isat: f64, vt: f64, p: NodeRef, n: NodeRef) -> Diode {
        Diode {
            isat: isat,
            vt: vt,
            p: p.into(),
            n: n.into(),
            ..Default::default()
        }
    }
}

impl Component for Diode {
    fn create_matrix_elems<T: SpNum>(&mut self, mat: &mut Matrix<T>) {
        self.pp = make_matrix_elem(mat, self.p, self.p);
        self.pn = make_matrix_elem(mat, self.p, self.n);
        self.np = make_matrix_elem(mat, self.n, self.p);
        self.nn = make_matrix_elem(mat, self.n, self.n);
    }
    fn load(&mut self, guess: &Variables<f64>, an: &AnalysisInfo) -> Stamps<f64> {
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
                (self.np, -di_dv),
            ],
            b: vec![(self.p, -irhs), (self.n, irhs)],
        };
    }
}

#[derive(Default)]
pub struct Isrc {
    i: f64,
    p: Option<VarIndex>,
    n: Option<VarIndex>,
}

impl Isrc {
    pub fn new(i: f64, p: NodeRef, n: NodeRef) -> Isrc {
        Isrc {
            i: i,
            p: p.into(),
            n: n.into(),
        }
    }
}

impl Component for Isrc {
    fn create_matrix_elems<T:SpNum>(&mut self, _mat: &mut Matrix<T>) {}
    fn load(&mut self, _guess: &Variables<f64>, _an: &AnalysisInfo) -> Stamps<f64> {
        return Stamps {
            g: vec![],
            b: vec![(self.p, self.i), (self.n, -self.i)],
        };
    }
}

/// Helper function to create matrix element at (row,col) if both are non-ground
fn make_matrix_elem<T: SpNum>(
    mat: &mut Matrix<T>,
    row: Option<VarIndex>,
    col: Option<VarIndex>,
) -> Option<Eindex> {
    if let (Some(r), Some(c)) = (row, col) {
        return Some(mat.make(r.0, c.0));
    }
    return None;
}
