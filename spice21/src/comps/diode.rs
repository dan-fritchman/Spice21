///
/// Diode Solver(s)
///
use super::consts;
use super::{make_matrix_elem, Component};
use crate::analysis::{AnalysisInfo, Options, Solver, Stamps, VarIndex, VarKind, Variables};
use crate::proto::D1;
use crate::sparse21::{Eindex, Matrix};
use crate::{attr, SpNum, SpResult};

/// Diode Model Parameters
attr!(
    DiodeModel,
    [
        (tnom, f64, 300.15, "Parameter measurement temperature"), // FIXME: defaults have to be literals in our macro
        (is, f64, 1e-14, "Saturation current"),
        (n, f64, 1.0, "Emission Coefficient"), //
        (tt, f64, 0.0, "Transit Time"),        //
        (vj, f64, 1.0, "Junction potential"),  //
        (m, f64, 0.5, "Grading coefficient"),  //
        (eg, f64, 1.11, "Activation energy"),  //
        (xti, f64, 3.0, "Saturation current temperature exp."), //
        (kf, f64, 0.0, "flicker noise coefficient"), //
        (af, f64, 1.0, "flicker noise exponent"), //
        (fc, f64, 0.5, "Forward bias junction fit parameter"), //
        (bv, f64, 0.0, "Reverse breakdown voltage"), // FIXME: Optional, default effectively -inf, default val
        (ibv, f64, 1e-3, "Current at reverse breakdown voltage"), //
        (rs, f64, 0.0, "Ohmic resistance"),          
        (cj0, f64, 0.0, "Junction capacitance"),     
        // Removed: 
        // (cjo, f64, 0.0, "Junction capacitance"), // FIXME: disallow                                             
        // (cond, f64, 0.0, "Ohmic conductance"),       // probably move this to internal
    ]
);

impl DiodeModel {
    /// Boolean indication of non-zero terminal resistance
    /// Used to determine whether to add an internal node
    pub fn has_rs(&self) -> bool {
        self.rs != 0.0
    }
    pub fn has_bv(&self) -> bool {
        self.bv != 0.0
    }
}

#[derive(Default)]
pub struct DiodeInstParams {
    pub temp: Option<f64>, // Instance temperature
    pub area: Option<f64>, // Area factor
}

#[derive(Clone, Copy, Default)]
pub struct DiodeOpPoint {
    pub vd: f64,     // "Diode voltage"),
    pub id: f64,     // "Diode current"),
    pub gd: f64,     // "Diode conductance"),
    pub cd: f64,     // "Diode capacitance"),
    pub charge: f64, // "Diode capacitor charge"),
    pub capcur: f64, // "Diode capacitor current"),
    pub p: f64,      // "Diode power"),
}

/// Diode Ports (or Variables)
/// Includes internal "r" node for terminal resistance,
/// on the "p" (cathode) side.
#[derive(Default)]
pub struct DiodePorts {
    pub p: Option<VarIndex>,
    pub n: Option<VarIndex>,
    pub r: Option<VarIndex>,
}

/// Diode Matrix-Pointers
/// Includes internal "r" node-pointers on "p" (cathode) side.
/// Total of seven matrix elements for the D+R series combo.
#[derive(Default)]
pub struct DiodeMatps {
    pp: Option<Eindex>,
    pr: Option<Eindex>,
    rp: Option<Eindex>,
    rr: Option<Eindex>,
    nr: Option<Eindex>,
    rn: Option<Eindex>,
    nn: Option<Eindex>,
}

/// Diode Internal Params
/// Derived from model and instance params at creation
#[derive(Default)]
pub struct DiodeIntParams {
    pub vt: f64,
    pub vte: f64,
    pub vcrit: f64,
    pub isat: f64,
    pub gspr: f64,
    pub cz: f64,
    pub cz2: f64,
    pub dep_threshold: f64,
    pub f1: f64,
    pub f2: f64,
    pub f3: f64,
    pub bv: f64, // Breakdown Voltage
}

impl DiodeIntParams {
    /// Derive Diode internal parameters from model, instance, and circuit options.
    fn derive(model: &DiodeModel, inst: &DiodeInstParams, opts: &Options) -> Self {
        let tnom = model.tnom;
        let temp = if let Some(t) = inst.temp {
            t
        } else {
            opts.temp
        };
        let area = if let Some(a) = inst.area { a } else { 1.0 };
        let gs = if model.has_rs() { 1.0 / model.rs } else { 0.0 };

        // Thermal voltage(s)
        let vt = consts::KB_OVER_Q * temp;
        let vtnom = consts::KB_OVER_Q * tnom;

        // This part gets really ugly - I won't even try to explain these equations.
        // (That's a SPICE joke.)
        let fact2 = temp / consts::TEMP_REF;
        let egfet = 1.16 - (7.02e-4 * temp * temp) / (temp + 1108.0);
        let arg =
            -egfet / (2.0 * consts::KB * temp) + 1.1150877 / (2.0 * consts::KB * consts::TEMP_REF);
        let pbfact = -2.0 * vt * (1.5 * fact2.ln() + consts::Q * arg);
        let egfet1 = 1.16 - (7.02e-4 * tnom) / (tnom + 1108.0);
        let arg1 =
            -egfet1 / (consts::KB * 2.0 * tnom) + 1.1150877 / (2.0 * consts::KB * consts::TEMP_REF);
        let fact1 = tnom / consts::TEMP_REF;
        let pbfact1 = -2.0 * vtnom * (1.5 * fact1.ln() + consts::Q * arg1);
        let pbo = (model.vj - pbfact1) / fact1;
        let gmaold = (model.vj - pbo) / pbo;
        let mut cjunc = model.cj0 / (1.0 + model.m * (400e-6 * (tnom - consts::TEMP_REF) - gmaold));
        let vjunc = pbfact + fact2 * pbo;
        let gmanew = (vjunc - pbo) / pbo;
        cjunc *= 1.0 + model.m * (400e-6 * (temp - consts::TEMP_REF) - gmanew);

        // Temperature-dependent saturation current
        let isat = model.is
            * (((temp / tnom) - 1.0) * model.eg / model.n * vt
                + model.xti / model.n * (temp / tnom).ln())
            .exp();
        let xfc = 1.0 - model.fc.ln();
        let f1 = vjunc * (1.0 - (1.0 - model.m * xfc).exp()) / (1.0 - model.m);
        let dep_threshold = model.fc * model.vj;
        let vte = model.n * vt;
        let vcrit = vte * (vte / (2.0 as f64).sqrt() / isat);

        let mut bv = model.bv;
        if model.has_bv() {
            // Fun part: iteratively update breakdown (voltage,current) for temperature
            // SPICE models this as `ibv` being constant across temperature, and `bv` changing to meet it.
            // We skip any convergence check here, and just get as close as we can in N iterations. 
            let ibv = model.ibv;
            for i in 0..25 {
                bv = model.bv - vt * (ibv / isat + 1.0 - bv / vt).ln();
            }
        }

        // Forward-bias depletion-cap fitting
        let f2 = (xfc * (1.0 + model.m)).exp();
        let f3 = 1.0 - model.fc * (1.0 + model.m);
        let gspr = gs * area;
        let cz = model.cj0 * area;
        let cz2 = cz / f2;

        return DiodeIntParams {
            vt,
            vte,
            vcrit,
            isat,
            gspr,
            cz,
            cz2,
            dep_threshold,
            f1,
            f2,
            f3,
            bv,
        };
    }
}

/// Diode Solver
#[derive(Default)]
pub struct Diode1 {
    pub ports: DiodePorts,
    pub model: DiodeModel,
    pub inst: DiodeInstParams,
    pub intp: DiodeIntParams,
    pub matps: DiodeMatps,
    pub op: DiodeOpPoint,
    pub guess: DiodeOpPoint,
}

impl Diode1 {
    /// Create a new Diode solver from a Circuit (parser) Diode
    pub fn from<T: SpNum>(d: D1, solver: &mut Solver<T>) -> Diode1 {
        // Destruct the key parser-diode attributes
        let D1 {
            mut name,
            model,
            inst,
            p,
            n,
        } = d;
        // Create or retrive the solver node-variables
        let p = solver.node_var(p);
        let n = solver.node_var(n);
        // Internal resistance node addition
        let r = if d.model.has_rs() {
            name.push_str("_r");
            Some(solver.vars.add(name, VarKind::V))
        } else {
            p
        };
        // Derive internal params
        let intp = DiodeIntParams::derive(&model, &inst, &solver.opts);
        // And create our solver
        return Diode1 {
            ports: DiodePorts { p, n, r },
            model,
            inst,
            intp,
            ..Default::default()
        };
    }
    /// Create a new Diode
    pub fn new(ports: DiodePorts, model: DiodeModel, inst: DiodeInstParams) -> Diode1 {
        Diode1 {
            ports,
            model,
            inst,
            ..Default::default()
        }
    }
    /// Voltage limiting
    fn limit(&self, vd: f64, past: Option<f64>) -> f64 {
        let vnew = vd;
        let vold = if let Some(v) = past { v } else { self.guess.vd };
        let DiodeIntParams { vte, vcrit, .. } = self.intp;

        // Typical case - unchanged
        if vnew <= vcrit || (vnew - vold).abs() <= 2.0 * vte {
            return vnew;
        }
        // Limiting cases
        if vold > 0.0 {
            let arg = 1.0 + (vnew - vold) / vte;
            if arg > 0.0 {
                return vold + vte * arg.ln();
            }
            return vcrit;
        }
        return vte * (vnew / vte).ln();
    }
}
impl Component for Diode1 {
    fn create_matrix_elems<T: SpNum>(&mut self, mat: &mut Matrix<T>) {
        self.matps.pp = make_matrix_elem(mat, self.ports.p, self.ports.p);
        self.matps.pr = make_matrix_elem(mat, self.ports.p, self.ports.r);
        self.matps.rp = make_matrix_elem(mat, self.ports.r, self.ports.p);
        self.matps.rr = make_matrix_elem(mat, self.ports.r, self.ports.r);
        self.matps.nr = make_matrix_elem(mat, self.ports.n, self.ports.r);
        self.matps.rn = make_matrix_elem(mat, self.ports.r, self.ports.n);
        self.matps.nn = make_matrix_elem(mat, self.ports.n, self.ports.n);
    }
    /// Parameter Validation
    fn validate(&self) -> SpResult<()> {
        // FIXME: convert these to Error types
        if self.model.m > 0.9 {
            panic!("diode grading coefficient too big!");
        }
        if self.model.eg < 0.1 {
            panic!("Diode activation energy too small!");
        }
        if self.model.fc > 0.95 {
            panic!("Diode fc too big!");
        }
        if self.model.bv < 0.0 {
            panic!("Diode breakdown voltages must be positive values");
        }
        Ok(())
    }
    /// Load our last guess as the new operating point
    fn commit(&mut self) {
        self.op = self.guess;
    }
    /// DC & Transient Stamp Loading
    fn load(&mut self, guess: &Variables<f64>, an: &AnalysisInfo) -> Stamps<f64> {
        let gmin_temp = 1e-15; // FIXME: from ckt
        let gmin = gmin_temp;
        // Destructure most key parameters
        let DiodeIntParams {
            vt,
            vte,
            isat,
            gspr,
            f1,
            f2,
            f3,
            cz,
            cz2,
            dep_threshold,
            bv,
            ..
        } = self.intp;

        // Extract our differential voltage
        let mut vd = guess.get(self.ports.r) - guess.get(self.ports.n);
        // Apply inter-estimate limits
        if self.model.has_bv() && vd < (10.0 * vte - bv).min(0.0) {
            let vtemp = self.limit(-bv, Some(bv - self.guess.vd));
            vd = vtemp - bv;
        } else {
            vd = self.limit(vd, None);
        }
        // Calculate diode current and its derivative, conductance
        let (mut id, mut gd) = if !self.model.has_bv() || vd >= -bv {
            // Regular (non-breakdown) operation
            let e = (vd / vte).exp();
            (isat * (e - 1.0) + gmin * vd, isat * e / vte + gmin)
        } else {
            // Breakdown - vd < BV
            let e = ((vd - bv) / vte).exp();
            (-isat * e + gmin * vd, isat * e / vte + gmin)
        };

        // Charge Storage Calculations
        let (qd, cd) = if vd < dep_threshold {
            let a = 1.0 - vd / self.model.vj;
            let s = -self.model.m * a.ln();
            let qd = self.model.tt * self.model.vj * cz * (1.0 - a * s) / (1.0 - self.model.m);
            let cd = self.model.tt * gd + cz * s;
            (qd, cd)
        } else {
            // Forward-bias model, adapted from Spice's polynomial approach
            let qd = self.model.tt * id
                + cz * f1
                + cz2
                    * (f3 * (vd - dep_threshold)
                        + self.model.m / 2.0 / self.model.vj
                            * (vd * vd - dep_threshold * dep_threshold));
            let cd = self.model.tt + cz2 * f3 + self.model.m * vd / self.model.vj;
            (qd, cd)
        };
        // If in transient, add the cap current and conductance
        let (gc, ic, _) = if let AnalysisInfo::TRAN(_, state) = an {
            state.integrate(qd - self.op.charge, cd, vd, self.op.capcur)
        } else {
            (0.0, 0.0, 0.0)
        };
        id += ic;
        gd += gc;

        // Update our op-point guess
        self.guess = DiodeOpPoint {
            vd,
            id,
            gd,
            cd,
            charge: qd,
            capcur: ic,
            p: vd * id,
        };
        // And finally return our matrix stamps
        let irhs = id - vd * gd;
        return Stamps {
            g: vec![
                (self.matps.nn, gd),
                (self.matps.rn, -gd),
                (self.matps.nr, -gd),
                (self.matps.rr, gd + gspr),
                (self.matps.pp, gspr),
                (self.matps.pr, -gspr),
                (self.matps.rp, -gspr),
            ],
            b: vec![(self.ports.r, -irhs), (self.ports.n, irhs)],
        };
    }
}

/// Simplified Diode Model, Level "Zero"
#[derive(Default)]
pub struct Diode0 {
    isat: f64,
    vt: f64,
    p: Option<VarIndex>,
    n: Option<VarIndex>,
    pp: Option<Eindex>,
    nn: Option<Eindex>,
    pn: Option<Eindex>,
    np: Option<Eindex>,
}

impl Diode0 {
    pub fn new(isat: f64, vt: f64, p: Option<VarIndex>, n: Option<VarIndex>) -> Diode0 {
        Diode0 {
            isat,
            vt,
            p,
            n,
            ..Default::default()
        }
    }
}

impl Component for Diode0 {
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
        let gd = (self.isat / self.vt) * (vd / self.vt).exp();
        let irhs = i - vd * gd;

        return Stamps {
            g: vec![(self.pp, gd), (self.nn, gd), (self.pn, -gd), (self.np, -gd)],
            b: vec![(self.p, -irhs), (self.n, irhs)],
        };
    }
}
