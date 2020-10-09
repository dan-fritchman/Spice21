//!
//! BSIM4 MOSFET Implementation
//!

#![allow(nonstandard_style,	warnings, unused)] // FIXME: work through these 

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

pub mod bsim4defs;
pub mod bsim4derive;
pub mod bsim4inst;
pub mod bsim4modelvals;
pub mod bsim4ports;
pub mod bsim4solver;

pub use bsim4defs::*;
pub use bsim4solver::*;

use super::consts::*;
use crate::sparse21::{Eindex, Matrix};



// FIXME: get from circuit/ analysis
pub(crate) const gmin: f64 = 1e-12;

pub(crate) const EXPL_THRESHOLD: f64 = 100.0;
pub(crate) const EXP_THRESHOLD: f64 = 34.0;
pub(crate) const MAX_EXP: f64 = 5.834617425e14;
pub(crate) const MIN_EXP: f64 = 1.713908431e-15;
pub(crate) const MAX_EXPL: f64 = 2.688117142e+43;
pub(crate) const MIN_EXPL: f64 = 3.720075976e-44;
pub(crate) const DELTA_1: f64 = 0.02;
pub(crate) const DELTA_2: f64 = 0.02;
pub(crate) const DELTA_3: f64 = 0.02;
pub(crate) const DELTA_4: f64 = 0.02;

#[derive(Clone)]
pub(crate) struct Bsim4ModelEntry {
    pub(crate) specs: Bsim4ModelSpecs,
    pub(crate) vals: Bsim4ModelVals,
    pub(crate) derived: Bsim4ModelDerivedParams,
    pub(crate) insts: Vec<Bsim4InstEntry>,
}

impl Bsim4ModelEntry {
    fn new(specs: Bsim4ModelSpecs) -> Self {
        use bsim4derive::derive;
        use bsim4modelvals::resolve;

        let vals = resolve(&specs);
        let derived = derive(&vals);
        Self {
            specs,
            vals,
            derived,
            insts: vec![],
        }
    }
}

#[derive(Clone)]
pub(crate) struct Bsim4InstEntry {
    pub(crate) specs: Bsim4InstSpecs,
    pub(crate) intp: Bsim4InternalParams,
    pub(crate) size_params: Bsim4SizeDepParams,
}

impl Bsim4InstEntry {
    fn new(specs: Bsim4InstSpecs, model: &Bsim4ModelEntry) -> Self {
        use bsim4inst::from;
        let (intp, size_params) = from(&model.vals, &model.derived, &specs);
        Self { specs, intp, size_params }
    }
}

pub(crate) struct Bsim4ModelCache(HashMap<String, Bsim4ModelEntry>);

impl Bsim4ModelCache {
    pub(crate) fn new() -> Self {
        Self(HashMap::new())
    }
    pub(crate) fn add<S: Into<String>>(&mut self, name: S, specs: Bsim4ModelSpecs) {
        let entry = Bsim4ModelEntry::new(specs);
        self.0.insert(name.into(), entry);
    }
    pub(crate) fn model(&mut self, model_name: &String) -> Option<&mut Bsim4ModelEntry> {
        self.0.get_mut(model_name)
    }
    pub(crate) fn inst(&mut self, model_name: &String, specs: Bsim4InstSpecs) -> Option<(Bsim4ModelEntry, Bsim4InstEntry)> {
        // FIXME: actually check whether these things are already in the cache!
        let model: &mut Bsim4ModelEntry = self.0.get_mut(model_name)?;
        let inst = Bsim4InstEntry::new(specs, &model);
        model.insts.push(inst.clone());
        // FIXME: stop cloning, return references or pointers 
        // Some((&*model, &model.insts[model.insts.len() - 1]))
        Some((model.clone(), inst))
    }
}

// Some helper math
// C-style call syntax, e.g. `min(a,b)` instead of `a.min(b)`
pub(crate) fn sqrt(a: f64) -> f64 {
    a.sqrt()
}
pub(crate) fn log(a: f64) -> f64 {
    a.ln()
}
pub(crate) fn exp(a: f64) -> f64 {
    a.exp()
}
pub(crate) fn MAX(a: f64, b: f64) -> f64 {
    a.max(b)
}
pub(crate) fn MIN(a: f64, b: f64) -> f64 {
    a.min(b)
}
pub(crate) fn abs(a: f64) -> f64 {
    a.abs()
}
pub(crate) fn atan(a: f64) -> f64 {
    a.atan()
}
pub(crate) fn pow(a: f64, b: f64) -> f64 {
    a.powf(b)
}
pub(crate) fn max(a: f64, b: f64) -> f64 {
    a.max(b)
}
pub(crate) fn dexpb(A: f64) -> f64 {
    if A > EXP_THRESHOLD {
        MAX_EXP * (1.0 + (A) - EXP_THRESHOLD)
    } else if A < -EXP_THRESHOLD {
        MIN_EXP
    } else {
        exp(A)
    }
}
pub(crate) fn dexpc(A: f64) -> f64 {
    if A > EXP_THRESHOLD {
        MAX_EXP
    } else if A < -EXP_THRESHOLD {
        0.0
    } else {
        exp(A)
    }
}

/// Bsim4 Internal, Derived Parameters
/// These are the numbers calculated offline and used during sim-time,
/// i.e. during `load`, `load_ac`, etc.
#[derive(Clone, Default, Debug, Deserialize, Serialize)]
pub(crate) struct Bsim4InternalParams {
    pub(crate) l: f64,
    pub(crate) w: f64,
    pub(crate) ad: f64,
    pub(crate) r#as: f64,
    pub(crate) nrd: f64,
    pub(crate) nrs: f64,
    pub(crate) pd: f64,
    pub(crate) ps: f64,
    pub(crate) sourceConductance: f64,
    pub(crate) drainConductance: f64,

    /* stress effect instance param */
    pub(crate) sa: f64,
    pub(crate) sb: f64,
    pub(crate) sd: f64,
    pub(crate) sca: f64,
    pub(crate) scb: f64,
    pub(crate) scc: f64,
    pub(crate) sc: f64,

    pub(crate) delvto: f64,
    pub(crate) xgw: f64,
    pub(crate) ngcon: f64,

    pub(crate) vjsmFwd: f64,
    pub(crate) vjsmRev: f64,
    pub(crate) vjdmFwd: f64,
    pub(crate) vjdmRev: f64,
    pub(crate) XExpBVS: f64,
    pub(crate) XExpBVD: f64,
    pub(crate) SslpFwd: f64,
    pub(crate) SslpRev: f64,
    pub(crate) DslpFwd: f64,
    pub(crate) DslpRev: f64,
    pub(crate) IVjsmFwd: f64,
    pub(crate) IVjsmRev: f64,
    pub(crate) IVjdmFwd: f64,
    pub(crate) IVjdmRev: f64,

    pub(crate) grgeltd: f64,
    pub(crate) Pseff: f64,
    pub(crate) Pdeff: f64,
    pub(crate) Aseff: f64,
    pub(crate) Adeff: f64,

    /* added here to account stress effect instance dependence */
    pub(crate) u0temp: f64,
    pub(crate) vsattemp: f64,
    pub(crate) vth0: f64,
    pub(crate) vfb: f64,
    pub(crate) vfbzb: f64,
    pub(crate) vtfbphi1: f64,
    pub(crate) vtfbphi2: f64,
    pub(crate) k2: f64,
    pub(crate) vbsc: f64,
    pub(crate) k2ox: f64,
    pub(crate) eta0: f64,
    pub(crate) nf: f64,
    pub(crate) grbsb: f64,
    pub(crate) grbdb: f64,
    pub(crate) grbpb: f64,
    pub(crate) grbps: f64,
    pub(crate) grbpd: f64,
    pub(crate) SjctTempRevSatCur: f64,
    pub(crate) DjctTempRevSatCur: f64,
    pub(crate) SswTempRevSatCur: f64,
    pub(crate) DswTempRevSatCur: f64,
    pub(crate) SswgTempRevSatCur: f64,
    pub(crate) DswgTempRevSatCur: f64,
    pub(crate) SourceSatCurrent: f64,
    pub(crate) DrainSatCurrent: f64,

    // Modes
    pub(crate) mode: usize,
    pub(crate) min: usize,
    pub(crate) rgeomod: usize,

    // Note these are *also* the names of model parameters
    pub(crate) toxp: f64,
    pub(crate) rbdb: f64,
    pub(crate) rbsb: f64,
    pub(crate) rbpb: f64,
    pub(crate) rbps: f64,
    pub(crate) rbpd: f64,

    // Instance-specific mode-selections not (yet?) supported
    pub(crate) trnqsmod: usize,
    pub(crate) acnqsmod: usize,
    pub(crate) rbodymod: usize,
    pub(crate) rgatemod: usize,
    geomod: usize,

    // Note these are *also* the names of (derived) model parameters
    pub(crate) coxp: f64,
}

/// Derived Bsim4 Model Parameters
/// Primarily params which depend on temperature, etc.
#[derive(Clone, Default, Debug, Deserialize, Serialize)]
pub(crate) struct Bsim4ModelDerivedParams {
    // Note these are *also* the names of internal instance parameters
    pub(crate) coxp: f64,

    pub(crate) Eg0: f64,
    pub(crate) vtm: f64,
    pub(crate) vtm0: f64,
    pub(crate) coxe: f64,
    pub(crate) vcrit: f64,
    pub(crate) factor1: f64,
    pub(crate) PhiBS: f64,
    pub(crate) PhiBSWS: f64,
    pub(crate) PhiBSWGS: f64,
    pub(crate) SjctTempSatCurDensity: f64,
    pub(crate) SjctSidewallTempSatCurDensity: f64,
    pub(crate) SjctGateSidewallTempSatCurDensity: f64,
    pub(crate) PhiBD: f64,
    pub(crate) PhiBSWD: f64,
    pub(crate) PhiBSWGD: f64,
    pub(crate) DjctTempSatCurDensity: f64,
    pub(crate) DjctSidewallTempSatCurDensity: f64,
    pub(crate) DjctGateSidewallTempSatCurDensity: f64,
    pub(crate) SunitAreaTempJctCap: f64,
    pub(crate) DunitAreaTempJctCap: f64,
    pub(crate) SunitLengthSidewallTempJctCap: f64,
    pub(crate) DunitLengthSidewallTempJctCap: f64,
    pub(crate) SunitLengthGateSidewallTempJctCap: f64,
    pub(crate) DunitLengthGateSidewallTempJctCap: f64,
    pub(crate) njtsstemp: f64,
    pub(crate) njtsswstemp: f64,
    pub(crate) njtsswgstemp: f64,
    pub(crate) njtsdtemp: f64,
    pub(crate) njtsswdtemp: f64,
    pub(crate) njtsswgdtemp: f64,
    pub(crate) TempRatio: f64,
    pub(crate) epssub: f64,
    pub(crate) ni: f64,
    pub(crate) Nvtms: f64,
    pub(crate) Nvtmd: f64,
    pub(crate) Nvtmrss: f64,
    pub(crate) Nvtmrssws: f64,
    pub(crate) Nvtmrsswgs: f64,
    pub(crate) Nvtmrsd: f64,
    pub(crate) Nvtmrsswd: f64,
    pub(crate) Nvtmrsswgd: f64,
}

#[derive(Clone, Default, Debug, Deserialize, Serialize)]
pub(crate) struct Bsim4OpPoint {
    pub(crate) mode: isize,

    pub(crate) vbd: f64,
    pub(crate) vbs: f64,
    pub(crate) vgs: f64,
    pub(crate) vds: f64,
    pub(crate) vdbs: f64,
    pub(crate) vdbd: f64,
    pub(crate) vsbs: f64,
    pub(crate) vges: f64,
    pub(crate) vgms: f64,
    pub(crate) vses: f64,
    pub(crate) vdes: f64,

    pub(crate) qb: f64,
    pub(crate) cqb: f64,
    pub(crate) qg: f64,
    pub(crate) cqg: f64,
    pub(crate) qd: f64,
    pub(crate) cqd: f64,
    pub(crate) qgmid: f64,
    pub(crate) cqgmid: f64,

    pub(crate) qbs: f64,
    pub(crate) cqbs: f64,
    pub(crate) qbd: f64,
    pub(crate) cqbd: f64,

    pub(crate) qcheq: f64,
    pub(crate) cqcheq: f64,
    pub(crate) qcdump: f64,
    pub(crate) cqcdump: f64,
    pub(crate) qdef: f64,
    pub(crate) qs: f64,

    pub(crate) von: f64,
    pub(crate) vdsat: f64,
    pub(crate) cgdo: f64,
    pub(crate) qgdo: f64,
    pub(crate) cgso: f64,
    pub(crate) qgso: f64,

    pub(crate) Vgsteff: f64,
    pub(crate) vgs_eff: f64,
    pub(crate) vgd_eff: f64,
    pub(crate) dvgs_eff_dvg: f64,
    pub(crate) dvgd_eff_dvg: f64,
    pub(crate) Vdseff: f64,
    pub(crate) nstar: f64,
    pub(crate) qinv: f64,
    pub(crate) cd: f64,
    pub(crate) cbs: f64,
    pub(crate) cbd: f64,
    pub(crate) csub: f64,
    pub(crate) Igidl: f64,
    pub(crate) Igisl: f64,
    pub(crate) gm: f64,
    pub(crate) gds: f64,
    pub(crate) gmbs: f64,
    pub(crate) gbd: f64,
    pub(crate) gbs: f64,
    pub(crate) noiGd0: f64,
    pub(crate) Coxeff: f64,

    pub(crate) gbbs: f64,
    pub(crate) gbgs: f64,
    pub(crate) gbds: f64,
    pub(crate) ggidld: f64,
    pub(crate) ggidlg: f64,
    pub(crate) ggidls: f64,
    pub(crate) ggidlb: f64,
    pub(crate) ggisld: f64,
    pub(crate) ggislg: f64,
    pub(crate) ggisls: f64,
    pub(crate) ggislb: f64,

    pub(crate) Igcs: f64,
    pub(crate) gIgcsg: f64,
    pub(crate) gIgcsd: f64,
    pub(crate) gIgcss: f64,
    pub(crate) gIgcsb: f64,
    pub(crate) Igcd: f64,
    pub(crate) gIgcdg: f64,
    pub(crate) gIgcdd: f64,
    pub(crate) gIgcds: f64,
    pub(crate) gIgcdb: f64,

    pub(crate) Igs: f64,
    pub(crate) gIgsg: f64,
    pub(crate) gIgss: f64,
    pub(crate) Igd: f64,
    pub(crate) gIgdg: f64,
    pub(crate) gIgdd: f64,

    pub(crate) Igb: f64,
    pub(crate) gIgbg: f64,
    pub(crate) gIgbd: f64,
    pub(crate) gIgbs: f64,
    pub(crate) gIgbb: f64,

    pub(crate) grdsw: f64,
    pub(crate) IdovVds: f64,
    pub(crate) gcrg: f64,
    pub(crate) gcrgd: f64,
    pub(crate) gcrgg: f64,
    pub(crate) gcrgs: f64,
    pub(crate) gcrgb: f64,

    pub(crate) gstot: f64,
    pub(crate) gstotd: f64,
    pub(crate) gstotg: f64,
    pub(crate) gstots: f64,
    pub(crate) gstotb: f64,

    pub(crate) gdtot: f64,
    pub(crate) gdtotd: f64,
    pub(crate) gdtotg: f64,
    pub(crate) gdtots: f64,
    pub(crate) gdtotb: f64,

    pub(crate) cggb: f64,
    pub(crate) cgdb: f64,
    pub(crate) cgsb: f64,
    pub(crate) cbgb: f64,
    pub(crate) cbdb: f64,
    pub(crate) cbsb: f64,
    pub(crate) cdgb: f64,
    pub(crate) cddb: f64,
    pub(crate) cdsb: f64,
    pub(crate) csgb: f64,
    pub(crate) csdb: f64,
    pub(crate) cssb: f64,
    pub(crate) cgbb: f64,
    pub(crate) cdbb: f64,
    pub(crate) csbb: f64,
    pub(crate) cbbb: f64,
    pub(crate) capbd: f64,
    pub(crate) capbs: f64,

    pub(crate) cqgb: f64,
    pub(crate) cqdb: f64,
    pub(crate) cqsb: f64,
    pub(crate) cqbb: f64,

    pub(crate) qgate: f64,
    pub(crate) qbulk: f64,
    pub(crate) qdrn: f64,
    pub(crate) qsrc: f64,

    pub(crate) qchqs: f64,
    pub(crate) taunet: f64,
    pub(crate) gtau: f64,
    pub(crate) gtg: f64,
    pub(crate) gtd: f64,
    pub(crate) gts: f64,
    pub(crate) gtb: f64,

    pub(crate) thetavth: f64,
    pub(crate) ueff: f64,
    pub(crate) Abulk: f64,
    pub(crate) EsatL: f64,
    pub(crate) AbovVgst2Vtm: f64,

    pub(crate) vgd: f64,
    pub(crate) sxpart: f64,
    pub(crate) dxpart: f64,
    pub(crate) vbs_jct: f64,
    pub(crate) vbd_jct: f64,
    pub(crate) gqdef: f64,
    pub(crate) ggtg: f64,
    pub(crate) ggtd: f64,
    pub(crate) ggts: f64,
    pub(crate) ggtb: f64,
    pub(crate) gcsbsb: f64,
    pub(crate) gcdbdb: f64,
    pub(crate) gcqsb: f64,
    pub(crate) gcqdb: f64,
    pub(crate) gcqgb: f64,

    pub(crate) cqdef: f64,
    pub(crate) gcqbb: f64,
    pub(crate) gcgbb: f64,
    pub(crate) ceqqg: f64,
    pub(crate) ceqqd: f64,
    pub(crate) ceqqb: f64,
    pub(crate) ceqqjs: f64,
    pub(crate) ceqqjd: f64,
    pub(crate) ceqqgmid: f64,
    pub(crate) gcsbb: f64,
    pub(crate) gcssb: f64,
    pub(crate) gcsgb: f64,
    pub(crate) gcsdb: f64,
    pub(crate) gcdbb: f64,
    pub(crate) gcggb: f64,
    pub(crate) gcgdb: f64,
    pub(crate) gcgsb: f64,
    pub(crate) gcgmgmb: f64,
    pub(crate) gcgmsb: f64,
    pub(crate) gcgmdb: f64,
    pub(crate) gcgmbb: f64,
    pub(crate) gcdgmb: f64,
    pub(crate) gcsgmb: f64,
    pub(crate) ddxpart_dVd: f64,
    pub(crate) ddxpart_dVs: f64,
    pub(crate) ddxpart_dVb: f64,
    pub(crate) dsxpart_dVb: f64,
    pub(crate) dsxpart_dVs: f64,
    pub(crate) dsxpart_dVg: f64,
    pub(crate) dsxpart_dVd: f64,
    pub(crate) ddxpart_dVg: f64,
    pub(crate) gcdgb: f64,
    pub(crate) gcdsb: f64,
    pub(crate) gcbdb: f64,
    pub(crate) gcbgb: f64,
    pub(crate) gcbsb: f64,
    pub(crate) gcbbb: f64,
    pub(crate) gcbgmb: f64,
    pub(crate) gcddb: f64,
}

#[derive(Clone, Default, Debug, Deserialize, Serialize)]
pub(crate) struct Bsim4SizeDepParams {
    pub(crate) Width: f64,
    pub(crate) Length: f64,
    pub(crate) NFinger: f64,

    pub(crate) cdsc: f64,
    pub(crate) cdscb: f64,
    pub(crate) cdscd: f64,
    pub(crate) cit: f64,
    pub(crate) nfactor: f64,
    pub(crate) xj: f64,
    pub(crate) vsat: f64,
    pub(crate) at: f64,
    pub(crate) a0: f64,
    pub(crate) ags: f64,
    pub(crate) a1: f64,
    pub(crate) a2: f64,
    pub(crate) keta: f64,
    pub(crate) nsub: f64,
    pub(crate) ndep: f64,
    pub(crate) nsd: f64,
    pub(crate) phin: f64,
    pub(crate) ngate: f64,
    pub(crate) gamma1: f64,
    pub(crate) gamma2: f64,
    pub(crate) vbx: f64,
    pub(crate) vbi: f64,
    pub(crate) vbm: f64,
    pub(crate) xt: f64,
    pub(crate) phi: f64,
    pub(crate) litl: f64,
    pub(crate) k1: f64,
    pub(crate) kt1: f64,
    pub(crate) kt1l: f64,
    pub(crate) kt2: f64,
    pub(crate) k2: f64,
    pub(crate) k3: f64,
    pub(crate) k3b: f64,
    pub(crate) w0: f64,
    pub(crate) dvtp0: f64,
    pub(crate) dvtp1: f64,
    pub(crate) dvtp2: f64,
    pub(crate) dvtp3: f64,
    pub(crate) dvtp4: f64,
    pub(crate) dvtp5: f64,
    pub(crate) lpe0: f64,
    pub(crate) lpeb: f64,
    pub(crate) dvt0: f64,
    pub(crate) dvt1: f64,
    pub(crate) dvt2: f64,
    pub(crate) dvt0w: f64,
    pub(crate) dvt1w: f64,
    pub(crate) dvt2w: f64,
    pub(crate) drout: f64,
    pub(crate) dsub: f64,
    pub(crate) vth0: f64,
    pub(crate) ua: f64,
    pub(crate) ua1: f64,
    pub(crate) ub: f64,
    pub(crate) ub1: f64,
    pub(crate) uc: f64,
    pub(crate) uc1: f64,
    pub(crate) ud: f64,
    pub(crate) ud1: f64,
    pub(crate) up: f64,
    pub(crate) lp: f64,
    pub(crate) u0: f64,
    pub(crate) eu: f64,
    pub(crate) ucs: f64,
    pub(crate) ute: f64,
    pub(crate) ucste: f64,
    pub(crate) voff: f64,
    pub(crate) tvoff: f64,
    pub(crate) tnfactor: f64,
    pub(crate) teta0: f64,
    pub(crate) tvoffcv: f64,
    pub(crate) minv: f64,
    pub(crate) minvcv: f64,
    pub(crate) vfb: f64,
    pub(crate) delta: f64,
    pub(crate) rdsw: f64,
    pub(crate) rds0: f64,
    pub(crate) rs0: f64,
    pub(crate) rd0: f64,
    pub(crate) rsw: f64,
    pub(crate) rdw: f64,
    pub(crate) prwg: f64,
    pub(crate) prwb: f64,
    pub(crate) prt: f64,
    pub(crate) eta0: f64,
    pub(crate) etab: f64,
    pub(crate) pclm: f64,
    pub(crate) pdibl1: f64,
    pub(crate) pdibl2: f64,
    pub(crate) pdiblb: f64,
    pub(crate) fprout: f64,
    pub(crate) pdits: f64,
    pub(crate) pditsd: f64,
    pub(crate) pscbe1: f64,
    pub(crate) pscbe2: f64,
    pub(crate) pvag: f64,
    pub(crate) wr: f64,
    pub(crate) dwg: f64,
    pub(crate) dwb: f64,
    pub(crate) b0: f64,
    pub(crate) b1: f64,
    pub(crate) alpha0: f64,
    pub(crate) alpha1: f64,
    pub(crate) beta0: f64,
    pub(crate) agidl: f64,
    pub(crate) bgidl: f64,
    pub(crate) cgidl: f64,
    pub(crate) egidl: f64,
    pub(crate) fgidl: f64,
    pub(crate) kgidl: f64,
    pub(crate) rgidl: f64,
    pub(crate) agisl: f64,
    pub(crate) bgisl: f64,
    pub(crate) cgisl: f64,
    pub(crate) egisl: f64,
    pub(crate) fgisl: f64,
    pub(crate) kgisl: f64,
    pub(crate) rgisl: f64,
    pub(crate) aigc: f64,
    pub(crate) bigc: f64,
    pub(crate) cigc: f64,
    pub(crate) aigsd: f64,
    pub(crate) bigsd: f64,
    pub(crate) cigsd: f64,
    pub(crate) aigs: f64,
    pub(crate) bigs: f64,
    pub(crate) cigs: f64,
    pub(crate) aigd: f64,
    pub(crate) bigd: f64,
    pub(crate) cigd: f64,
    pub(crate) aigbacc: f64,
    pub(crate) bigbacc: f64,
    pub(crate) cigbacc: f64,
    pub(crate) aigbinv: f64,
    pub(crate) bigbinv: f64,
    pub(crate) cigbinv: f64,
    pub(crate) nigc: f64,
    pub(crate) nigbacc: f64,
    pub(crate) nigbinv: f64,
    pub(crate) ntox: f64,
    pub(crate) eigbinv: f64,
    pub(crate) pigcd: f64,
    pub(crate) poxedge: f64,
    pub(crate) xrcrg1: f64,
    pub(crate) xrcrg2: f64,
    pub(crate) lambda: f64,   /* overshoot */
    pub(crate) vtl: f64,      /* thermal velocity limit */
    pub(crate) xn: f64,       /* back scattering parameter */
    pub(crate) lc: f64,       /* back scattering parameter */
    pub(crate) tfactor: f64,  /* ballistic transportation factor  */
    pub(crate) vfbsdoff: f64, /* S/D flatband offset voltage  */
    pub(crate) tvfbsdoff: f64,

    /* added for stress effect */
    pub(crate) ku0: f64,
    pub(crate) kvth0: f64,
    pub(crate) ku0temp: f64,
    pub(crate) rho_ref: f64,
    pub(crate) inv_od_ref: f64,
    /* added for well proximity effect */
    pub(crate) kvth0we: f64,
    pub(crate) k2we: f64,
    pub(crate) ku0we: f64,

    /* CV model */
    pub(crate) cgsl: f64,
    pub(crate) cgdl: f64,
    pub(crate) ckappas: f64,
    pub(crate) ckappad: f64,
    pub(crate) cf: f64,
    pub(crate) clc: f64,
    pub(crate) cle: f64,
    pub(crate) vfbcv: f64,
    pub(crate) noff: f64,
    pub(crate) voffcv: f64,
    pub(crate) acde: f64,
    pub(crate) moin: f64,

    /* Pre-calculated constants */
    pub(crate) dw: f64,
    pub(crate) dl: f64,
    pub(crate) leff: f64,
    pub(crate) weff: f64,

    pub(crate) dwc: f64,
    pub(crate) dlc: f64,
    pub(crate) dwj: f64,
    pub(crate) leffCV: f64,
    pub(crate) weffCV: f64,
    pub(crate) weffCJ: f64,
    pub(crate) abulkCVfactor: f64,
    pub(crate) cgso: f64,
    pub(crate) cgdo: f64,
    pub(crate) cgbo: f64,

    pub(crate) u0temp: f64,
    pub(crate) vsattemp: f64,
    pub(crate) sqrtPhi: f64,
    pub(crate) phis3: f64,
    pub(crate) Xdep0: f64,
    pub(crate) sqrtXdep0: f64,
    pub(crate) theta0vb0: f64,
    pub(crate) thetaRout: f64,
    pub(crate) mstar: f64,
    pub(crate) VgsteffVth: f64,
    pub(crate) mstarcv: f64,
    pub(crate) voffcbn: f64,
    pub(crate) voffcbncv: f64,
    pub(crate) rdswmin: f64,
    pub(crate) rdwmin: f64,
    pub(crate) rswmin: f64,
    pub(crate) vfbsd: f64,

    pub(crate) cof1: f64,
    pub(crate) cof2: f64,
    pub(crate) cof3: f64,
    pub(crate) cof4: f64,
    pub(crate) cdep0: f64,
    pub(crate) ToxRatio: f64,
    pub(crate) Aechvb: f64,
    pub(crate) Bechvb: f64,
    pub(crate) ToxRatioEdge: f64,
    pub(crate) AechvbEdgeS: f64,
    pub(crate) AechvbEdgeD: f64,
    pub(crate) BechvbEdge: f64,
    pub(crate) ldeb: f64,
    pub(crate) k1ox: f64,
    pub(crate) k2ox: f64,
    pub(crate) vfbzbfactor: f64,
    pub(crate) dvtp2factor: f64, /* v4.7 */
}

/// # BSIM4 Matrix Pointers
#[derive(Default)]
pub(crate) struct Bsim4MatrixPointers {
    GEgePtr: Option<Eindex>,
    GPgePtr: Option<Eindex>,
    GEgpPtr: Option<Eindex>,
    GPgpPtr: Option<Eindex>,
    GPdpPtr: Option<Eindex>,
    GPspPtr: Option<Eindex>,
    GPbpPtr: Option<Eindex>,
    GEdpPtr: Option<Eindex>,
    GEspPtr: Option<Eindex>,
    GEbpPtr: Option<Eindex>,

    GEgmPtr: Option<Eindex>,
    GMgePtr: Option<Eindex>,
    GMgmPtr: Option<Eindex>,

    GMdpPtr: Option<Eindex>,
    GMgpPtr: Option<Eindex>,
    GMspPtr: Option<Eindex>,
    GMbpPtr: Option<Eindex>,

    DPgmPtr: Option<Eindex>,
    GPgmPtr: Option<Eindex>,
    SPgmPtr: Option<Eindex>,
    BPgmPtr: Option<Eindex>,

    DgpPtr: Option<Eindex>,
    DspPtr: Option<Eindex>,
    DbpPtr: Option<Eindex>,
    SdpPtr: Option<Eindex>,
    SgpPtr: Option<Eindex>,
    SbpPtr: Option<Eindex>,

    DPdpPtr: Option<Eindex>,
    DPdPtr: Option<Eindex>,
    DPgpPtr: Option<Eindex>,
    DPspPtr: Option<Eindex>,
    DPbpPtr: Option<Eindex>,

    DdpPtr: Option<Eindex>,
    DdPtr: Option<Eindex>,

    SPdpPtr: Option<Eindex>,
    SPgpPtr: Option<Eindex>,
    SPspPtr: Option<Eindex>,
    SPsPtr: Option<Eindex>,
    SPbpPtr: Option<Eindex>,

    SspPtr: Option<Eindex>,
    SsPtr: Option<Eindex>,

    BPdpPtr: Option<Eindex>,
    BPgpPtr: Option<Eindex>,
    BPspPtr: Option<Eindex>,
    BPbpPtr: Option<Eindex>,

    DPdbPtr: Option<Eindex>,
    SPsbPtr: Option<Eindex>,

    DBdpPtr: Option<Eindex>,
    DBdbPtr: Option<Eindex>,
    DBbpPtr: Option<Eindex>,
    DBbPtr: Option<Eindex>,

    BPdbPtr: Option<Eindex>,
    BPbPtr: Option<Eindex>,
    BPsbPtr: Option<Eindex>,

    SBspPtr: Option<Eindex>,
    SBbpPtr: Option<Eindex>,
    SBbPtr: Option<Eindex>,
    SBsbPtr: Option<Eindex>,

    BdbPtr: Option<Eindex>,
    BbpPtr: Option<Eindex>,
    BsbPtr: Option<Eindex>,
    BbPtr: Option<Eindex>,

    QqPtr: Option<Eindex>,
    QgpPtr: Option<Eindex>,
    QdpPtr: Option<Eindex>,
    QspPtr: Option<Eindex>,
    QbpPtr: Option<Eindex>,

    DPqPtr: Option<Eindex>,
    SPqPtr: Option<Eindex>,
    GPqPtr: Option<Eindex>,
}
