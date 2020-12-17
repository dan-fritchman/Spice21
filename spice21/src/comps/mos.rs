//!
//! # MOS Solvers
//!

use num::Complex;
use serde::{Deserialize, Serialize};
use std::convert::From;
use std::ops::{Index, IndexMut};

use super::consts;
use super::{make_matrix_elem, Component};
use crate::analysis::{AnalysisInfo, Options, Stamps, VarIndex, Variables};
use crate::sparse21::{Eindex, Matrix};
use crate::{proto, SpNum};

/// Mos Terminals, in SPICE order: d, g, s, b
#[derive(Clone, Copy)]
pub enum MosTerm {
    D = 0,
    G = 1,
    S = 2,
    B = 3,
}
use MosTerm::{B, D, G, S};

#[derive(Default)]
pub struct MosPorts<T> {
    pub d: T,
    pub g: T,
    pub s: T,
    pub b: T,
}
/// Index MosPorts by the `MosTerm` enum
impl<T> Index<MosTerm> for MosPorts<T> {
    type Output = T;
    fn index(&self, t: MosTerm) -> &T {
        match t {
            D => &self.d,
            G => &self.g,
            S => &self.s,
            B => &self.b,
        }
    }
}
/// Very fun conversion from four-element arrays into MosPorts of `From`-able types.
impl<S, T: Clone + Into<S>> From<[T; 4]> for MosPorts<S> {
    fn from(n: [T; 4]) -> MosPorts<S> {
        return MosPorts {
            d: n[0].clone().into(),
            g: n[1].clone().into(),
            s: n[2].clone().into(),
            b: n[3].clone().into(),
        };
    }
}
/// Even more fun conversion from four-element tuples into MosPorts of `From`-able types.
/// Note in this case, each of the four elements can be of distinct types.
impl<S, T: Clone + Into<S>, U: Clone + Into<S>, V: Clone + Into<S>, W: Clone + Into<S>> From<(T, U, V, W)> for MosPorts<S> {
    fn from(n: (T, U, V, W)) -> MosPorts<S> {
        return MosPorts {
            d: n.0.clone().into(),
            g: n.1.clone().into(),
            s: n.2.clone().into(),
            b: n.3.clone().into(),
        };
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

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub enum MosType {
    NMOS,
    PMOS,
}
impl Default for MosType {
    fn default() -> MosType {
        MosType::NMOS
    }
}
impl MosType {
    /// Polarity Function
    /// The very common need to negate values for PMOS, and leave NMOS unchanged.
    pub fn p(&self) -> f64 {
        match self {
            MosType::PMOS => -1.0,
            MosType::NMOS => 1.0,
        }
    }
}

/// Mos Level 1 Model Parameters
#[derive(Clone)]
pub struct Mos1Model {
    pub mos_type: MosType,
    pub vt0: f64,
    pub kp: f64,
    pub gamma: f64,
    pub phi: f64,
    pub lambda: f64,
    pub cbd: f64,
    pub cbs: f64,
    pub is: f64,
    pub pb: f64,
    pub cgso: f64,
    pub cgdo: f64,
    pub cgbo: f64,
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
    pub nss: f64,
    pub tnom: f64,
    pub kf: f64,
    pub af: f64,
    pub rd: Option<f64>,
    pub rs: Option<f64>,
    pub rsh: Option<f64>,
    pub tpg: bool,
}
impl Mos1Model {
    pub(crate) fn resolve(specs: &proto::Mos1Model) -> Self {
        Self {
            mos_type: if specs.mos_type == 1 { MosType::PMOS } else { MosType::NMOS },
            vt0: if let Some(val) = specs.vt0 { val } else { 0.0 },
            kp: if let Some(val) = specs.kp { val } else { 2.0e-5 },
            gamma: if let Some(val) = specs.gamma { val } else { 0.0 },
            phi: if let Some(val) = specs.phi { val } else { 0.6 },
            lambda: if let Some(val) = specs.lambda { val } else { 0.0 },
            cbd: if let Some(val) = specs.cbd { val } else { 0.0 },
            cbs: if let Some(val) = specs.cbs { val } else { 0.0 },
            is: if let Some(val) = specs.is { val } else { 1.0e-14 },
            pb: if let Some(val) = specs.pb { val } else { 0.8 },
            cgso: if let Some(val) = specs.cgso { val } else { 0.0 },
            cgdo: if let Some(val) = specs.cgdo { val } else { 0.0 },
            cgbo: if let Some(val) = specs.cgbo { val } else { 0.0 },
            cj: if let Some(val) = specs.cj { val } else { 0.0 },
            mj: if let Some(val) = specs.mj { val } else { 0.5 },
            cjsw: if let Some(val) = specs.cjsw { val } else { 0.0 },
            mjsw: if let Some(val) = specs.mjsw { val } else { 0.5 },
            js: if let Some(val) = specs.js { val } else { 1.0e-8 }, // FIXME
            tox: if let Some(val) = specs.tox { val } else { 1.0e-7 },
            nsub: if let Some(val) = specs.nsub { val } else { 0.0 },
            nss: if let Some(val) = specs.nss { val } else { 0.0 },
            ld: if let Some(val) = specs.ld { val } else { 0.0 },
            u0: if let Some(val) = specs.u0 { val } else { 600.0 },
            fc: if let Some(val) = specs.fc { val } else { 0.5 },
            kf: if let Some(val) = specs.kf { val } else { 0.0 },
            af: if let Some(val) = specs.af { val } else { 1.0 },
            tnom: if let Some(val) = specs.tnom {
                val + consts::KELVIN_TO_C
            } else {
                consts::TEMP_REF
            }, // C to Kelvin conversion, right here
            rd: specs.rd, // Options
            rs: specs.rs,
            rsh: specs.rsh,
            tpg: specs.tpg,
        }
    }
    pub(crate) fn p(&self) -> f64 {
        self.mos_type.p()
    }
}
impl Default for Mos1Model {
    fn default() -> Self {
        Self::resolve(&proto::Mos1Model::default())
    }
}

/// Mos Level 1 Instance Parameters
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
    temp: Option<f64>,
    // FIXME: maybe even more explicitly ignore these
    // dtemp: Option<f64>,
    // off: bool,
    // icvds: f64,
    // icvgs: f64,
    // icvbs: f64,
    // ic: f64,
}
impl Mos1InstanceParams {
    pub(crate) fn resolve(specs: &proto::Mos1InstParams) -> Self {
        Mos1InstanceParams {
            m: if let Some(val) = specs.m { val } else { 0.0 },
            l: if let Some(val) = specs.l { val } else { 1e-6 },
            w: if let Some(val) = specs.w { val } else { 1e-6 },
            a_d: if let Some(val) = specs.a_d { val } else { 1e-12 },
            a_s: if let Some(val) = specs.a_s { val } else { 1e-12 },
            pd: if let Some(val) = specs.pd { val } else { 1e-6 },
            ps: if let Some(val) = specs.ps { val } else { 1e-6 },
            nrd: if let Some(val) = specs.nrd { val } else { 1.0 },
            nrs: if let Some(val) = specs.nrs { val } else { 1.0 },
            temp: specs.temp
            // dtemp: if let Some(val) = specs.dtemp { val } else { 0.0 },
            // icvds: if let Some(val) = specs.icvds { val } else { 0.0 },
            // icvgs: if let Some(val) = specs.icvgs { val } else { 0.0 },
            // icvbs: if let Some(val) = specs.icvbs { val } else { 0.0 },
            // ic: if let Some(val) = specs.ic { val } else { 0.0 },
            // off: specs.off,
        }
    }
}
impl Default for Mos1InstanceParams {
    fn default() -> Self {
        Self::resolve(&proto::Mos1InstParams::default())
    }
}

/// Mos1 Internal "Parameters", derived at instance-construction
/// and updated only on changes in temperature
struct Mos1InternalParams {
    temp: f64,
    vtherm: f64,
    vt0_t: f64,
    kp_t: f64,
    phi_t: f64,
    beta: f64,
    cox: f64,
    cgs_ov: f64,
    cgd_ov: f64,
    cgb_ov: f64,
    leff: f64,
    drain_junc: MosJunction,
    source_junc: MosJunction,
    grd: f64,
    grs: f64,
}

enum SourceDrain {
    S,
    D,
}
struct MosJunction {
    area: f64,
    isat: f64,
    depletion_threshold: f64,
    vcrit: f64,
    czb: f64,
    czbsw: f64,
    f2: f64,
    f3: f64,
    f4: f64,
    sd: SourceDrain,
}

/// Mos1 DC & Transient Operating Point
#[derive(Default, Clone)]
struct Mos1OpPoint {
    ids: f64,
    vgs: f64,
    vds: f64,
    vgd: f64,
    vgb: f64,
    gm: f64,
    gds: f64,
    gmbs: f64,
    cgs: f64,
    cgd: f64,
    cgb: f64,
    qgs: f64,
    qgd: f64,
    qgb: f64,
    icgs: f64,
    icgd: f64,
    icgb: f64,
    reversed: bool,
}

///
/// # Mos Level 1 Solver
///
pub struct Mos1 {
    model: Mos1Model,
    _params: Mos1InstanceParams,
    intparams: Mos1InternalParams,
    op: Mos1OpPoint,
    guess: Mos1OpPoint,
    ports: MosPorts<Option<VarIndex>>,
    matps: MosMatrixPointers,
}
impl Mos1 {
    pub(crate) fn new(model: Mos1Model, params: Mos1InstanceParams, ports: MosPorts<Option<VarIndex>>, opts: &Options) -> Mos1 {
        let intparams = Mos1::derive(&model, &params, opts);
        Mos1 {
            model,
            _params: params,
            intparams,
            ports,
            op: Mos1OpPoint::default(),
            guess: Mos1OpPoint::default(),
            matps: MosMatrixPointers::default(),
        }
    }
    /// Calculate derived parameters from instance parameters
    fn derive(model: &Mos1Model, inst: &Mos1InstanceParams, opts: &Options) -> Mos1InternalParams {
        if let Some(t) = inst.temp {
            panic!("Mos1 Instance Temperatures Are Not Supported");
        }
        let temp = opts.temp; // Note: in Kelvin

        use consts::{KB, KB_OVER_Q, Q, TEMP_REF};

        // FIXME: offload these parts to derived Models
        let cox_per_area = consts::SIO2_PERMITTIVITY / model.tox;
        // Nominal temperature params
        let fact1 = model.tnom / TEMP_REF;
        let vtnom = model.tnom * KB_OVER_Q;
        let kt1 = KB * model.tnom;
        let egfet1 = 1.16 - (7.02e-4 * model.tnom.powi(2)) / (model.tnom + 1108.0);
        let arg1 = -egfet1 / 2.0 / kt1 + 1.1150877 / (KB * 2.0 * TEMP_REF);
        let pbfact1 = -2.0 * vtnom * (1.5 * fact1.ln() + Q * arg1);

        // FIXME: more model derivations to come

        // Instance temperature params
        let kt = temp * KB;
        let vtherm = temp * KB_OVER_Q;
        let temp_ratio = temp / model.tnom;
        let fact2 = temp / TEMP_REF;
        let egfet = 1.16 - (7.02e-4 * temp.powi(2)) / (temp + 1108.0);
        let arg = -egfet / 2.0 / kt + 1.1150877 / (KB * 2.0 * TEMP_REF);
        let pbfact = -2.0 * vtherm * (1.5 * fact2.ln() + Q * arg);

        // Effective Length
        let leff = inst.l - 2.0 * model.ld;
        if leff < 0.0 {
            panic!("Mos1 Effective Length < 0");
        }

        let phio = (model.phi - pbfact1) / fact1;
        let phi_t = fact2 * phio + pbfact;
        let vbi_t = model.vt0 - model.p() * (model.gamma * model.phi.sqrt()) + 0.5 * (egfet1 - egfet) + model.p() * 0.5 * (phi_t - model.phi);
        let vt0_t = vbi_t + model.p() * model.gamma * phi_t.sqrt();
        let isat_t = model.is * (-egfet / vtherm + egfet1 / vtnom).exp();
        let jsat_t = model.js * (-egfet / vtherm + egfet1 / vtnom).exp();

        let pbo = (model.pb - pbfact1) / fact1;
        let gmaold = (model.pb - pbo) / pbo;
        let capfact = 1.0 / (1.0 + model.mj * (4e-4 * (model.tnom - TEMP_REF) - gmaold));
        let mut cbd_t = model.cbd * capfact;
        let mut cbs_t = model.cbs * capfact;
        let mut cj_t = model.cj * capfact;
        let capfact = 1.0 / (1.0 + model.mjsw * (4e-4 * (model.tnom - TEMP_REF) - gmaold));
        let mut cjsw_t = model.cjsw * capfact;
        let bulkpot_t = fact2 * pbo + pbfact;
        let gmanew = (bulkpot_t - pbo) / pbo;
        let capfact = 1.0 / (1.0 + model.mj * (4e-4 * (temp - TEMP_REF) - gmanew));
        cbd_t *= capfact;
        cbs_t *= capfact;
        cj_t *= capfact;
        let capfact = 1.0 / (1.0 + model.mjsw * (4e-4 * (temp - TEMP_REF) - gmanew));
        cjsw_t *= capfact;

        // S/D Junction Params
        let depletion_threshold = model.fc * bulkpot_t;
        let arg = 1.0 - model.fc;
        let sarg = ((-model.mj) * arg.ln()).exp();
        let sargsw = ((-model.mjsw) * arg.ln()).exp();
        let use_default_isat: bool = jsat_t == 0.0 || inst.a_d == 0.0 || inst.a_s == 0.0;

        // MosJunction construction-closure
        let junc_new = |area: f64, perim: f64, sd: SourceDrain| {
            let isat = if use_default_isat { isat_t } else { jsat_t * area };
            let vcrit = vtherm * (vtherm / (consts::SQRT2 * isat)).ln();
            let czb = match sd {
                SourceDrain::D => {
                    if model.cbd != 0.0 {
                        cbd_t
                    } else {
                        cj_t * area
                    }
                }
                SourceDrain::S => {
                    if model.cbs != 0.0 {
                        cbs_t
                    } else {
                        cj_t * area
                    }
                }
            };
            let czbsw = cjsw_t * perim;
            let f2 = czb * (1.0 - model.fc * (1.0 + model.mj)) * sarg / arg + czbsw * (1.0 - model.fc * (1.0 + model.mjsw)) * sargsw / arg;
            let f3 = czb * model.mj * sarg / arg / bulkpot_t + czbsw * model.mjsw * sargsw / arg / bulkpot_t;
            let f4 = czb * bulkpot_t * (1.0 - arg * sarg) / (1.0 - model.mj) + czbsw * bulkpot_t * (1.0 - arg * sargsw) / (1.0 - model.mjsw)
                - f3 / 2.0 * (depletion_threshold * depletion_threshold)
                - depletion_threshold * f2;

            MosJunction {
                area,
                isat,
                depletion_threshold,
                vcrit,
                czb,
                czbsw,
                f2,
                f3,
                f4,
                sd,
            }
        };

        // Create the source & drain junction params
        let drain_junc = junc_new(inst.a_d, inst.pd, SourceDrain::D);
        let source_junc = junc_new(inst.a_s, inst.ps, SourceDrain::S);

        // Terminal Ohmic Resistances
        let grs = if let Some(r) = model.rs {
            if r <= 0.0 {
                println!("Warning: Mos1 Model with rs <= 0");
                0.0
            } else {
                1.0 / r
            }
        } else if let Some(rsh) = model.rsh {
            if rsh <= 0.0 {
                println!("Warning: Mos1 Model with rsh <= 0");
                0.0
            } else {
                1.0 / rsh / inst.nrs
            }
        } else {
            0.0
        };
        let grd = if let Some(r) = model.rd {
            if r <= 0.0 {
                println!("Warning: Mos1 Model with rd <= 0");
                0.0
            } else {
                1.0 / r
            }
        } else if let Some(rsh) = model.rsh {
            if rsh <= 0.0 {
                println!("Warning: Mos1 Model with rsh <= 0");
                0.0
            } else {
                1.0 / rsh / inst.nrd
            }
        } else {
            0.0
        };

        // Temperature-adjusted transconductance 
        let kp_t = model.kp / temp_ratio * temp_ratio.sqrt();
        Mos1InternalParams {
            vt0_t,
            kp_t,
            temp,
            vtherm,
            leff,
            cox: cox_per_area * leff * inst.w,
            beta: kp_t * inst.w / leff,
            phi_t,
            drain_junc,
            source_junc,
            cgs_ov: inst.w * model.cgso,
            cgd_ov: inst.w * model.cgdo,
            cgb_ov: leff * model.cgbo,
            grs,
            grd,
        }
    }
    /// Gather the voltages on each of our node-variables from `Variables` `guess`.
    fn vs(&self, vars: &Variables<f64>) -> MosPorts<f64> {
        MosPorts {
            d: vars.get(self.ports[D]),
            g: vars.get(self.ports[G]),
            s: vars.get(self.ports[S]),
            b: vars.get(self.ports[B]),
        }
    }
    /// Primary action behind dc & transient loading.
    /// Returns calculated "guess" operating point, plus matrix stamps
    fn op_stamp(&self, v: MosPorts<f64>, an: &AnalysisInfo, opts: &Options) -> (Mos1OpPoint, Stamps<f64>) {
        let gmin = opts.gmin;
        // Initially factor out polarity of NMOS/PMOS and source/drain swapping
        // All math after this block uses increasing vgs,vds <=> increasing ids,
        // i.e. the polarities typically expressed for NMOS
        let p = self.model.mos_type.p();
        let reversed = p * (v.d - v.s) < 0.0;
        // FIXME: add inter-step limiting
        let (vd, vs) = if reversed { (v.s, v.d) } else { (v.d, v.s) };
        let vgs = p * (v.g - vs);
        let vgd = p * (v.g - vd);
        let vds = p * (vd - vs);
        let vgb = p * (v.g - v.b);
        // Same for bulk junction diodes - polarities such that more `vsb`, `vdb` = more *reverse* bias.
        let vsb = p * (vs - v.b);
        let vdb = p * (vd - v.b);

        // Threshold & body effect calcs
        let von = if vsb > 0.0 {
            self.intparams.vt0_t + self.model.gamma * ((self.intparams.phi_t + vsb).sqrt() - self.intparams.phi_t.sqrt())
        } else {
            self.intparams.vt0_t // FIXME: body effect for Vsb < 0
        };
        let vov = vgs - von;
        let vdsat = vov.max(0.0);

        // Drain current & its g-derivatives
        // Default to cutoff values
        let mut ids = 0.0;
        let mut gm = 0.0;
        let mut gds = 0.0;
        let mut gmbs = 0.0;
        if vov > 0.0 {
            if vds >= vov {
                // Sat
                ids = self.intparams.beta / 2.0 * vov.powi(2) * (1.0 + self.model.lambda * vds);
                gm = self.intparams.beta * vov * (1.0 + self.model.lambda * vds);
                gds = self.model.lambda * self.intparams.beta / 2.0 * vov.powi(2);
            } else {
                // Triode
                ids = self.intparams.beta * (vov * vds - vds.powi(2) / 2.0) * (1.0 + self.model.lambda * vds);
                gm = self.intparams.beta * vds * (1.0 + self.model.lambda * vds);
                gds = self.intparams.beta * ((vov - vds) * (1.0 + self.model.lambda * vds) + self.model.lambda * ((vov * vds) - vds.powi(2) / 2.0));
            }
            gmbs = if self.intparams.phi_t + vsb > 0.0 {
                gm * self.model.gamma / 2.0 / (self.intparams.phi_t + vsb).sqrt()
            } else {
                0.0
            };
        }

        // Bulk Junction Diodes
        let vtherm = self.intparams.vtherm;
        let isat_bs = self.intparams.source_junc.isat;
        let isat_bd = self.intparams.drain_junc.isat;
        // Source-Bulk
        let ibs = isat_bs * ((-vsb / vtherm).exp() - 1.0);
        let gbs = (isat_bs / vtherm) * (-vsb / vtherm).exp() + gmin;
        let ibs_rhs = ibs + vsb * gbs;
        // Drain-Bulk
        let ibd = isat_bd * ((-vdb / vtherm).exp() - 1.0);
        let gbd = (isat_bd / vtherm) * (-vdb / vtherm).exp() + gmin;
        let ibd_rhs = ibd + vdb * gbd;

        // Capacitance Calculations
        let cox = self.intparams.cox;
        let cgs1: f64;
        let cgd1: f64;
        let cgb1: f64;
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
            // This is the fake initial-time check to be cleaned
            cgs1
        } else if reversed == self.op.reversed {
            self.op.cgs
        } else {
            self.op.cgd
        };
        let cgs = cgs1 + cgs2 + self.intparams.cgs_ov;
        let cgd = cgd1 + self.intparams.cgd_ov + if reversed == self.op.reversed { self.op.cgd } else { self.op.cgs };
        let cgb = cgb1 + self.intparams.cgb_ov + self.op.cgb;

        // Numerical integrations for cap currents and impedances
        let (gcgs, icgs, rhsgs) = if let AnalysisInfo::TRAN(_, state) = an {
            let dqgs = if reversed == self.op.reversed {
                (vgs - self.op.vgs) * cgs
            } else {
                (vgs - self.op.vgd) * cgs
            };
            let ip = if reversed == self.op.reversed { self.op.icgs } else { self.op.icgd };
            state.integrate(dqgs, cgs, vgs, ip)
        } else {
            (0.0, 0.0, 0.0)
        };
        let (gcgd, icgd, rhsgd) = if let AnalysisInfo::TRAN(_, state) = an {
            let dqgd = if reversed == self.op.reversed {
                (vgd - self.op.vgd) * cgd
            } else {
                (vgd - self.op.vgs) * cgd
            };
            let ip = if reversed == self.op.reversed { self.op.icgd } else { self.op.icgs };
            state.integrate(dqgd, cgd, vgd, ip)
        } else {
            (0.0, 0.0, 0.0)
        };
        let (gcgb, icgb, rhsgb) = if let AnalysisInfo::TRAN(_, state) = an {
            let dqgb = (vgb - self.op.vgb) * cgb;
            state.integrate(dqgb, cgb, vgb, self.op.icgb)
        } else {
            (0.0, 0.0, 0.0)
        };

        // FIXME: bulk junction diode caps
        let _cbs = 0.0;
        let _cbd = 0.0;
        let _gcbd = 0.0;
        let _gcbs = 0.0;

        // Store as our op point for next time
        let guess = Mos1OpPoint {
            ids,
            vgs,
            vds,
            vgd,
            vgb,
            cgs: cgs1,
            cgd: cgd1,
            cgb: cgb1,
            qgs: 0.0, // FIXME: unused
            qgd: 0.0,
            qgb: 0.0,
            icgs,
            icgd,
            icgb,
            gm,
            gds,
            gmbs,
            reversed,
        };

        // FIXME: Bulk junction caps RHS adjustment
        // let _ceqbs = p * (cbs + gbs * vsb);
        // let _ceqbd = 0.0; //p * (cbd + gbd * vdb);
        let irhs = ids - gm * vgs - gds * vds;

        // Sort out which are the "reported" drain and source terminals (sr, dr)
        // FIXME: this also needs the "prime" vs "external" source & drains
        let (sr, sx, dr, dx) = if !reversed { (S, S, D, D) } else { (D, D, S, S) };
        // Include our terminal resistances
        let grd = self.intparams.grd;
        let grs = self.intparams.grs;
        // And finally send back our matrix contributions
        let stamps = Stamps {
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
                (self.matps[(B, G)], -gcgb),
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
                (self.ports[dr], p * (-irhs + ibd_rhs + rhsgd)),
                (self.ports[sr], p * (irhs + ibs_rhs + rhsgs)),
                (self.ports[G], -p * (rhsgs + rhsgb + rhsgd)),
                (self.ports[B], -p * (ibd_rhs + ibs_rhs - rhsgb)),
            ],
        };
        (guess, stamps)
    }
}

impl Component for Mos1 {
    fn create_matrix_elems<T: SpNum>(&mut self, mat: &mut Matrix<T>) {
        for t1 in [G, D, S, B].iter() {
            for t2 in [G, D, S, B].iter() {
                self.matps[(*t1, *t2)] = make_matrix_elem(mat, self.ports[*t1], self.ports[*t2]);
            }
        }
    }
    fn commit(&mut self) {
        // Load our last guess as the new operating point
        self.op = self.guess.clone();
    }
    fn load(&mut self, vars: &Variables<f64>, an: &AnalysisInfo, opts: &Options) -> Stamps<f64> {
        let v = self.vs(vars); // Collect terminal voltages
        let (op, stamps) = self.op_stamp(v, an, opts); // Do most of our work here
        self.guess = op; // Save the calculated operating point
        stamps // And return our matrix stamps
    }
    fn load_ac(&mut self, _guess: &Variables<Complex<f64>>, an: &AnalysisInfo, _opts: &Options) -> Stamps<Complex<f64>> {
        // Grab the frequency-variable from our analysis
        let omega = match an {
            AnalysisInfo::AC(_opts, state) => state.omega,
            _ => panic!("Invalid AC AnalysisInfo"),
        };

        // Short-hand the conductances from our op-point.
        // (Rustc should be smart enough not to copy these.)
        let (gm, gds, gmbs) = (self.op.gm, self.op.gds, self.op.gmbs);

        // Cap admittances
        let gcgs = omega * self.op.cgs;
        let gcgd = omega * self.op.cgd;
        let gcgb = omega * self.op.cgb;

        // FIXME: bulk junction diodes
        let _cbs = 0.0;
        let _cbd = 0.0;
        let gbd = 1e-9;
        let gbs = 1e-9;
        let _gcbd = 0.0;
        let _gcbs = 0.0;
        let gcbg = 0.0;

        // Sort out which are the "reported" drain and source terminals (sr, dr)
        // FIXME: this also needs the "prime" vs "external" source & drains
        let (sr, sx, dr, dx) = if !self.op.reversed { (S, S, D, D) } else { (D, D, S, S) };

        // Include our terminal resistances
        let grd = self.intparams.grd;
        let grs = self.intparams.grs;

        // And finally, send back our AC-matrix contributions
        return Stamps {
            g: vec![
                (self.matps[(dr, dr)], Complex::new(gds + grd + gbd, gcgd)),
                (self.matps[(sr, sr)], Complex::new(gm + gds + grs + gbs + gmbs, gcgs)),
                (self.matps[(dr, sr)], Complex::new(-gm - gds - gmbs, 0.0)),
                (self.matps[(sr, dr)], Complex::new(-gds, 0.0)),
                (self.matps[(dr, G)], Complex::new(gm, -gcgd)),
                (self.matps[(sr, G)], Complex::new(-gm, -gcgs)),
                (self.matps[(G, G)], Complex::new(0.0, gcgd + gcgs + gcgb)),
                (self.matps[(B, B)], Complex::new(gbd + gbs, gcgb)),
                (self.matps[(G, B)], Complex::new(0.0, -gcgb)),
                (self.matps[(G, dr)], Complex::new(0.0, -gcgd)),
                (self.matps[(G, sr)], Complex::new(0.0, -gcgs)),
                (self.matps[(B, G)], Complex::new(0.0, -gcbg)),
                (self.matps[(G, dr)], Complex::new(0.0, -gcgd)),
                (self.matps[(B, dr)], Complex::new(-gbd, 0.0)),
                (self.matps[(B, sr)], Complex::new(-gbs, 0.0)),
                (self.matps[(dr, B)], Complex::new(-gbd + gmbs, 0.0)),
                (self.matps[(sr, B)], Complex::new(-gbs - gmbs, 0.0)),
                (self.matps[(dx, dr)], Complex::new(-grd, 0.0)),
                (self.matps[(dr, dx)], Complex::new(-grd, 0.0)),
                (self.matps[(dx, dx)], Complex::new(grd, 0.0)),
                (self.matps[(sx, sr)], Complex::new(-grs, 0.0)),
                (self.matps[(sr, sx)], Complex::new(-grs, 0.0)),
                (self.matps[(sx, sx)], Complex::new(grs, 0.0)),
            ],
            b: vec![],
        };
    }
}

use std::collections::HashMap;
///
/// # Mos1 Model and Instance-Param Definitions
///
#[derive(Default)]
pub struct Mos1Defs {
    pub(crate) models: HashMap<String, proto::Mos1Model>,
    pub(crate) insts: HashMap<String, proto::Mos1InstParams>,
}

/// Mos Level-Zero Instance Parameters
pub(crate) struct Mos0Params {
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

/// Mos "Level Zero" Simplified Solver
pub struct Mos0 {
    params: Mos0Params,
    ports: MosPorts<Option<VarIndex>>,
    matps: MosMatrixPointers,
}
impl Mos0 {
    pub(crate) fn new(ports: MosPorts<Option<VarIndex>>, mos_type: MosType) -> Self {
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
        let matps = [(D, D), (S, S), (D, S), (S, D), (D, G), (S, G)];
        for (t1, t2) in matps.iter() {
            self.matps[(*t1, *t2)] = make_matrix_elem(mat, self.ports[*t1], self.ports[*t2]);
        }
    }
    fn load(&mut self, guess: &Variables<f64>, _an: &AnalysisInfo, _opts: &Options) -> Stamps<f64> {
        let vg = guess.get(self.ports[G]);
        let vd = guess.get(self.ports[D]);
        let vs = guess.get(self.ports[S]);

        let p = self.params.mos_type.p();
        let vds1 = p * (vd - vs);
        let reversed = vds1 < 0.0;
        let vgs = if reversed { p * (vg - vd) } else { p * (vg - vs) };
        let vds = if reversed { -vds1 } else { vds1 };
        let vov = vgs - self.params.vth;

        // Cutoff conditions
        let mut ids = 0.0;
        let mut gm = 0.0;
        let mut gds = 0.0;
        if vov > 0.0 {
            let lam = self.params.lam;
            let beta = self.params.beta;
            if vds >= vov {
                // Saturation
                ids = beta / 2.0 * vov.powi(2) * (1.0 + lam * vds);
                gm = beta * vov * (1.0 + lam * vds);
                gds = lam * beta / 2.0 * vov.powi(2);
            } else {
                // Triode
                ids = beta * (vov * vds - vds.powi(2) / 2.0) * (1.0 + lam * vds);
                gm = beta * vds * (1.0 + lam * vds);
                gds = beta * ((vov - vds) * (1.0 + lam * vds) + lam * ((vov * vds) - vds.powi(2) / 2.0));
            }
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
