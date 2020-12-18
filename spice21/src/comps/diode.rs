//!
//! # Diode Solver(s)
//!
use std::collections::HashMap;

use super::consts;
use super::{make_matrix_elem, Component};
use crate::analysis::{AnalysisInfo, Options, Stamps, VarIndex, VarKind, Variables};
use crate::defs::DefPtr;
use crate::proto;
use crate::sparse21::{Eindex, Matrix};
use crate::{attr, from_opt_type, sperror, SpNum, SpResult};

pub(crate) use crate::proto::DiodeInstParams;

// Diode Model Parameters
attr!(
    DiodeModel,
    "Diode Model Parameters",
    [
        (tnom, f64, 300.15, "Parameter measurement temperature"), // FIXME: defaults have to be literals in our macro
        (is, f64, 1e-14, "Saturation current"),
        (n, f64, 1.0, "Emission Coefficient"),                    //
        (tt, f64, 0.0, "Transit Time"),                           //
        (vj, f64, 1.0, "Junction potential"),                     //
        (m, f64, 0.5, "Grading coefficient"),                     //
        (eg, f64, 1.11, "Activation energy"),                     //
        (xti, f64, 3.0, "Saturation current temperature exp."),   //
        (kf, f64, 0.0, "flicker noise coefficient"),              //
        (af, f64, 1.0, "flicker noise exponent"),                 //
        (fc, f64, 0.5, "Forward bias junction fit parameter"),    //
        (bv, f64, 0.0, "Reverse breakdown voltage"),              // FIXME: Optional, default effectively -inf, default val
        (ibv, f64, 1e-3, "Current at reverse breakdown voltage"), //
        (rs, f64, 0.0, "Ohmic resistance"),
        (cj0, f64, 0.0, "Junction capacitance"),
        // Removed, redudant params:
        // (cjo, f64, 0.0, "Junction capacitance"),
        // (cond, f64, 0.0, "Ohmic conductance"),
    ]
);
impl DiodeModel {
    /// Boolean indication of non-zero terminal resistance
    /// Used to determine whether to add an internal node
    pub(crate) fn has_rs(&self) -> bool {
        self.rs != 0.0
    }
    pub(crate) fn has_bv(&self) -> bool {
        self.bv != 0.0
    }
    /// Derive a `DiodeModel` from (`Option`-based) `proto::DiodeModel`
    /// Apply defaults for all unspecified fields
    pub(crate) fn from(specs: proto::DiodeModel) -> Self {
        Self {
            tnom: if let Some(val) = specs.tnom { val } else { 300.15 },
            is: if let Some(val) = specs.is { val } else { 1e-14 },
            n: if let Some(val) = specs.n { val } else { 1.0 },
            tt: if let Some(val) = specs.tt { val } else { 0.0 },
            vj: if let Some(val) = specs.vj { val } else { 1.0 },
            m: if let Some(val) = specs.m { val } else { 0.5 },
            eg: if let Some(val) = specs.eg { val } else { 1.11 },
            xti: if let Some(val) = specs.xti { val } else { 3.0 },
            kf: if let Some(val) = specs.kf { val } else { 0.0 },
            af: if let Some(val) = specs.af { val } else { 1.0 },
            fc: if let Some(val) = specs.fc { val } else { 0.5 },
            bv: if let Some(val) = specs.bv { val } else { 0.0 },
            ibv: if let Some(val) = specs.ibv { val } else { 1e-3 },
            rs: if let Some(val) = specs.rs { val } else { 0.0 },
            cj0: if let Some(val) = specs.cj0 { val } else { 0.0 },
        }
    }
}
impl Default for DiodeModel {
    /// Default DiodeModel, derived from the all-default-value proto-
    fn default() -> Self {
        Self::from(proto::DiodeModel::default())
    }
}

/// Diode Operating Point
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
impl DiodePorts {
    pub(crate) fn from<P: Clone + Into<Option<VarIndex>>, T: SpNum>(path: String, model: &DiodeModel, p: P, n: P, vars: &mut Variables<T>) -> Self {
        // Internal resistance node addition
        let r = if model.has_rs() {
            let name = format!("{}.{}", path, "r");
            Some(vars.add(name, VarKind::V))
        } else {
            p.clone().into()
        };
        Self { p: p.into(), n: n.into(), r }
    }
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
    pub(crate) fn derive(model: &DiodeModel, inst: &DiodeInstParams, opts: &Options) -> Self {
        let tnom = model.tnom;
        let temp = if let Some(t) = inst.temp { t } else { opts.temp };
        let area = if let Some(a) = inst.area { a } else { 1.0 };
        let gs = if model.has_rs() { 1.0 / model.rs } else { 0.0 };

        // Thermal voltage(s)
        let vt = consts::KB_OVER_Q * temp;
        let vtnom = consts::KB_OVER_Q * tnom;

        // This part gets really ugly - I won't even try to explain these equations.
        // (That's a SPICE joke.)
        let fact2 = temp / consts::TEMP_REF;
        let egfet = 1.16 - (7.02e-4 * temp * temp) / (temp + 1108.0);
        let arg = -egfet / (2.0 * consts::KB * temp) + 1.1150877 / (2.0 * consts::KB * consts::TEMP_REF);
        let pbfact = -2.0 * vt * (1.5 * fact2.ln() + consts::Q * arg);
        let egfet1 = 1.16 - (7.02e-4 * tnom) / (tnom + 1108.0);
        let arg1 = -egfet1 / (consts::KB * 2.0 * tnom) + 1.1150877 / (2.0 * consts::KB * consts::TEMP_REF);
        let fact1 = tnom / consts::TEMP_REF;
        let pbfact1 = -2.0 * vtnom * (1.5 * fact1.ln() + consts::Q * arg1);
        let pbo = (model.vj - pbfact1) / fact1;
        let gmaold = (model.vj - pbo) / pbo;
        let mut cjunc = model.cj0 / (1.0 + model.m * (400e-6 * (tnom - consts::TEMP_REF) - gmaold));
        let vjunc = pbfact + fact2 * pbo;
        let gmanew = (vjunc - pbo) / pbo;
        cjunc *= 1.0 + model.m * (400e-6 * (temp - consts::TEMP_REF) - gmanew);

        // Temperature-dependent saturation current
        let isat = model.is * (((temp / tnom) - 1.0) * model.eg / model.n * vt + model.xti / model.n * (temp / tnom).ln()).exp();
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
            for _i in 0..25 {
                bv = model.bv - vt * (ibv / isat + 1.0 - bv / vt).ln();
            }
        }
        // Forward-bias depletion-cap fitting
        let f2 = (xfc * (1.0 + model.m)).exp();
        let f3 = 1.0 - model.fc * (1.0 + model.m);
        let gspr = gs * area;
        let cz = model.cj0 * area;
        let cz2 = cz / f2;

        DiodeIntParams {
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
        }
    }
}

/// Diode Solver
#[derive(Default)]
pub struct Diode {
    pub ports: DiodePorts,
    pub model: DefPtr<DiodeModel>,
    // pub inst: DiodeInstParams,
    pub intp: DefPtr<DiodeIntParams>,
    pub matps: DiodeMatps,
    pub op: DiodeOpPoint,
    pub guess: DiodeOpPoint,
}
impl Diode {
    /// Voltage limiting
    fn limit(&self, vd: f64, past: Option<f64>) -> f64 {
        let vnew = vd;
        let vold = if let Some(v) = past { v } else { self.guess.vd };
        let intp = &*self.intp.read();
        // Typical case - unchanged
        if vnew <= intp.vcrit || (vnew - vold).abs() <= 2.0 * intp.vte {
            return vnew;
        }
        // Limiting cases
        if vold > 0.0 {
            let arg = 1.0 + (vnew - vold) / intp.vte;
            if arg > 0.0 {
                return vold + intp.vte * arg.ln();
            }
            return intp.vcrit;
        }
        return intp.vte * (vnew / intp.vte).ln();
    }
}
impl Component for Diode {
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
        let model = &*self.model.read();
        if model.m > 0.9 {
            return Err(sperror("diode grading coefficient too big!"));
        }
        if model.eg < 0.1 {
            return Err(sperror("Diode activation energy too small!"));
        }
        if model.fc > 0.95 {
            return Err(sperror("Diode fc too big!"));
        }
        if model.bv < 0.0 {
            return Err(sperror("Diode breakdown voltages must be positive values"));
        }
        Ok(())
    }
    /// Load our last guess as the new operating point
    fn commit(&mut self) {
        self.op = self.guess;
    }
    /// DC & Transient Stamp Loading
    fn load(&mut self, guess: &Variables<f64>, an: &AnalysisInfo, opts: &Options) -> Stamps<f64> {
        // Grab the data from our shared attributes
        let model = &*self.model.read();
        let intp = &*self.intp.read();
        // Grab all relevant options
        let gmin = opts.gmin;

        // Extract our differential voltage
        let mut vd = guess.get(self.ports.r) - guess.get(self.ports.n);
        // Apply inter-estimate limits
        if model.has_bv() && vd < (10.0 * intp.vte - intp.bv).min(0.0) {
            let vtemp = self.limit(-intp.bv, Some(intp.bv - self.guess.vd));
            vd = vtemp - intp.bv;
        } else {
            vd = self.limit(vd, None);
        }
        // Calculate diode current and its derivative, conductance
        let (mut id, mut gd) = if !model.has_bv() || vd >= -intp.bv {
            // Regular (non-breakdown) operation
            let e = (vd / intp.vte).exp();
            (intp.isat * (e - 1.0) + gmin * vd, intp.isat * e / intp.vte + gmin)
        } else {
            // Breakdown - vd < BV
            let e = ((vd - intp.bv) / intp.vte).exp();
            (-intp.isat * e + gmin * vd, intp.isat * e / intp.vte + gmin)
        };

        // Charge Storage Calculations
        let (qd, cd) = if vd < intp.dep_threshold {
            let a = 1.0 - vd / model.vj;
            let s = -model.m * a.ln();
            let qd = model.tt * model.vj * intp.cz * (1.0 - a * s) / (1.0 - model.m);
            let cd = model.tt * gd + intp.cz * s;
            (qd, cd)
        } else {
            // Forward-bias model, adapted from Spice's polynomial approach
            let qd = model.tt * id
                + intp.cz * intp.f1
                + intp.cz2 * (intp.f3 * (vd - intp.dep_threshold) + model.m / 2.0 / model.vj * (vd * vd - intp.dep_threshold * intp.dep_threshold));
            let cd = model.tt + intp.cz2 * intp.f3 + model.m * vd / model.vj;
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
                (self.matps.rr, gd + intp.gspr),
                (self.matps.pp, intp.gspr),
                (self.matps.pr, -intp.gspr),
                (self.matps.rp, -intp.gspr),
            ],
            b: vec![(self.ports.r, -irhs), (self.ports.n, irhs)],
        };
    }
}

/// Simplified Diode Model, Level "Zero"
#[derive(Default)]
pub(crate) struct Diode0 {
    isat: f64,
    vt: f64,
    p: Option<VarIndex>,
    n: Option<VarIndex>,
    pp: Option<Eindex>,
    nn: Option<Eindex>,
    pn: Option<Eindex>,
    np: Option<Eindex>,
}
impl Component for Diode0 {
    fn create_matrix_elems<T: SpNum>(&mut self, mat: &mut Matrix<T>) {
        self.pp = make_matrix_elem(mat, self.p, self.p);
        self.pn = make_matrix_elem(mat, self.p, self.n);
        self.np = make_matrix_elem(mat, self.n, self.p);
        self.nn = make_matrix_elem(mat, self.n, self.n);
    }
    fn load(&mut self, guess: &Variables<f64>, _an: &AnalysisInfo, _opts: &Options) -> Stamps<f64> {
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

///
/// # Diode Model and Instance-Param Definitions
///
/// Definitions of DiodeModels and DiodeInstanceParams are stored in HashMaps
/// `models` and `insts` as they are defined.
/// When instances are requested via the `get` method,
/// internal parameters (`DiodeInternalParams`) are derived,
/// and stored alongside pointers to their models in the `cache` Hash.
///
#[derive(Default)]
pub struct DiodeDefs {
    models: HashMap<String, DefPtr<DiodeModel>>,
    insts: HashMap<String, DefPtr<DiodeInstParams>>,
    cache: HashMap<String, DiodeCacheEntry>,
}
/// Diode Cache Entry
/// Includes the internal/ derived, instance, and model parameters
/// that fully characterize a Diode instance
pub(crate) struct DiodeCacheEntry {
    pub(crate) model: DefPtr<DiodeModel>,
    pub(crate) inst: DefPtr<DiodeInstParams>,
    pub(crate) intp: DefPtr<DiodeIntParams>,
}
impl DiodeCacheEntry {
    fn clone(&self) -> Self {
        Self {
            model: DefPtr::clone(&self.model),
            inst: DefPtr::clone(&self.inst),
            intp: DefPtr::clone(&self.intp),
        }
    }
}
impl DiodeDefs {
    pub(crate) fn add_model(&mut self, name: &str, specs: DiodeModel) {
        self.models.insert(name.to_string(), DefPtr::new(specs));
    }
    pub(crate) fn add_inst(&mut self, name: &str, inst: DiodeInstParams) {
        self.insts.insert(name.to_string(), DefPtr::new(inst));
    }
    pub(crate) fn get(&mut self, inst_name: &str, opts: &Options) -> Option<DiodeCacheEntry> {
        // If we've already derived these parameters, clone a new pointer to them 
        if let Some(e) = self.cache.get(inst_name) {
            return Some(e.clone());
        }

        // Not in cache, check whether we have definitions.
        let instptr = self.insts.get(inst_name)?;
        let inst = &*instptr.read();
        let modelptr = self.models.get(&inst.model)?;

        // If we get here, we found definitions of both instance and model params.
        // Now derive the internal ones, including any circuit options.
        let intp = DiodeIntParams::derive(&*modelptr.read(), inst, opts);

        // Create new pointers and a new cache entry for the new combo.
        let intp = DefPtr::new(intp);
        let entry = DiodeCacheEntry {
            intp,
            inst: DefPtr::clone(&instptr),
            model: DefPtr::clone(&modelptr),
        };
        self.cache.insert(inst_name.to_string(), entry.clone());
        Some(entry)
    }
}
