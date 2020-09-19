//!
//! BSIM4 MOSFET Implementation
//!

use super::consts::{Q, VT_REF};
use super::Component;
use crate::analysis::{AnalysisInfo, Stamps, VarIndex, Variables};
use crate::sparse21::{Eindex, Matrix};
use crate::SpNum;
use std::f64::{MAX_EXP as MAX_EXPI, MIN_EXP as MIN_EXPI};

use super::bsim4defs::*;

// FIXME: get from circuit/ analysis
const gmin: f64 = 1e-9;

const MAX_EXP: f64 = MAX_EXPI as f64;
const MIN_EXP: f64 = MIN_EXPI as f64;
const MAX_EXPL: f64 = 2.688117142e+43;
const MIN_EXPL: f64 = 3.720075976e-44;
const EXPL_THRESHOLD: f64 = 100.0;
const EXP_THRESHOLD: f64 = 34.0;
const EPS0: f64 = 8.85418e-12;
const EPSSI: f64 = 1.03594e-10;
const DELTA_1: f64 = 0.02;
const DELTA_2: f64 = 0.02;
const DELTA_3: f64 = 0.02;
const DELTA_4: f64 = 0.02;
const MM: f64 = 3.0;

fn DEXPb(A: f64) -> f64 {
    if A > EXP_THRESHOLD {
        MAX_EXP * (1.0 + (A) - EXP_THRESHOLD)
    } else if A < -EXP_THRESHOLD {
        MIN_EXP
    } else {
        exp(A)
    }
}
fn DEXPc(A: f64) -> f64 {
    if A > EXP_THRESHOLD {
        MAX_EXP
    } else if A < -EXP_THRESHOLD {
        0.0
    } else {
        exp(A)
    }
}

impl Bsim4Model {
    /// Polarity function
    fn p(&self) -> f64 {
        // FIXME: use mos_type enum
        if self.nmos {
            1.0
        } else {
            -1.0
        }
    }
}

pub struct Bsim4Ports {
    dNode: Option<VarIndex>,
    gNodeExt: Option<VarIndex>,
    sNode: Option<VarIndex>,
    bNode: Option<VarIndex>,
    dNodePrime: Option<VarIndex>,
    gNodePrime: Option<VarIndex>,
    gNodeMid: Option<VarIndex>,
    sNodePrime: Option<VarIndex>,
    bNodePrime: Option<VarIndex>,
    dbNode: Option<VarIndex>,
    sbNode: Option<VarIndex>,
    qNode: Option<VarIndex>,
}

pub struct Bsim4InternalParams {
    vjsmFwd: f64,
    vjsmRev: f64,
    vjdmFwd: f64,
    vjdmRev: f64,
    XExpBVS: f64,
    XExpBVD: f64,
    SslpFwd: f64,
    SslpRev: f64,
    DslpFwd: f64,
    DslpRev: f64,
    IVjsmFwd: f64,
    IVjsmRev: f64,
    IVjdmFwd: f64,
    IVjdmRev: f64,

    grgeltd: f64,
    Pseff: f64,
    Pdeff: f64,
    Aseff: f64,
    Adeff: f64,

    l: f64,
    w: f64,
    drainArea: f64,
    sourceArea: f64,
    drainSquares: f64,
    sourceSquares: f64,
    drainPerimeter: f64,
    sourcePerimeter: f64,
    sourceConductance: f64,
    drainConductance: f64,
    /* stress effect instance param */
    sa: f64,
    sb: f64,
    sd: f64,
    sca: f64,
    scb: f64,
    scc: f64,
    sc: f64,

    rbdb: f64,
    rbsb: f64,
    rbpb: f64,
    rbps: f64,
    rbpd: f64,
    delvto: f64,
    xgw: f64,
    ngcon: f64,

    /* added here to account stress effect instance dependence */
    u0temp: f64,
    vsattemp: f64,
    vth0: f64,
    vfb: f64,
    vfbzb: f64,
    vtfbphi1: f64,
    vtfbphi2: f64,
    k2: f64,
    vbsc: f64,
    k2ox: f64,
    eta0: f64,
    toxp: f64,
    coxp: f64,

    icVDS: f64,
    icVGS: f64,
    icVBS: f64,
    nf: f64,

    off: isize,
    mode: isize,
    trnqsMod: isize,
    acnqsMod: isize,
    rbodyMod: isize,
    rgateMod: isize,
    geoMod: isize,
    rgeoMod: isize,
    min: isize,

    SjctTempRevSatCur: f64,
    DjctTempRevSatCur: f64,
    SswTempRevSatCur: f64,
    DswTempRevSatCur: f64,
    SswgTempRevSatCur: f64,
    DswgTempRevSatCur: f64,

    // Model-derived; perhaps break out
    Eg0: f64,
    vtm: f64,
    vtm0: f64,
    coxe: f64,
    cof1: f64,
    cof2: f64,
    cof3: f64,
    cof4: f64,
    vcrit: f64,
    factor1: f64,
    PhiBS: f64,
    PhiBSWS: f64,
    PhiBSWGS: f64,
    SjctTempSatCurDensity: f64,
    SjctSidewallTempSatCurDensity: f64,
    SjctGateSidewallTempSatCurDensity: f64,
    PhiBD: f64,
    PhiBSWD: f64,
    PhiBSWGD: f64,
    DjctTempSatCurDensity: f64,
    DjctSidewallTempSatCurDensity: f64,
    DjctGateSidewallTempSatCurDensity: f64,
    SunitAreaTempJctCap: f64,
    DunitAreaTempJctCap: f64,
    SunitLengthSidewallTempJctCap: f64,
    DunitLengthSidewallTempJctCap: f64,
    SunitLengthGateSidewallTempJctCap: f64,
    DunitLengthGateSidewallTempJctCap: f64,

    oxideTrapDensityA: f64,
    oxideTrapDensityB: f64,
    oxideTrapDensityC: f64,
    em: f64,
    ef: f64,
    af: f64,
    kf: f64,
    sheetResistance: f64,
    SjctSatCurDensity: f64,
    DjctSatCurDensity: f64,
    SjctSidewallSatCurDensity: f64,
    DjctSidewallSatCurDensity: f64,
    SjctGateSidewallSatCurDensity: f64,
    DjctGateSidewallSatCurDensity: f64,
    SbulkJctPotential: f64,
    DbulkJctPotential: f64,
    SbulkJctBotGradingCoeff: f64,
    DbulkJctBotGradingCoeff: f64,
    SbulkJctSideGradingCoeff: f64,
    DbulkJctSideGradingCoeff: f64,
    SbulkJctGateSideGradingCoeff: f64,
    DbulkJctGateSideGradingCoeff: f64,
    SsidewallJctPotential: f64,
    DsidewallJctPotential: f64,
    SGatesidewallJctPotential: f64,
    DGatesidewallJctPotential: f64,
    SunitAreaJctCap: f64,
    DunitAreaJctCap: f64,
    SunitLengthSidewallJctCap: f64,
    DunitLengthSidewallJctCap: f64,
    SunitLengthGateSidewallJctCap: f64,
    DunitLengthGateSidewallJctCap: f64,
    SjctEmissionCoeff: f64,
    DjctEmissionCoeff: f64,
    SjctTempExponent: f64,
    DjctTempExponent: f64,
    njtsstemp: f64,
    njtsswstemp: f64,
    njtsswgstemp: f64,
    njtsdtemp: f64,
    njtsswdtemp: f64,
    njtsswgdtemp: f64,

    TempRatio: f64,
}

#[derive(Default)]
struct Bsim4OpPoint {
    mode: isize,

    vbd: f64,
    vbs: f64,
    vgs: f64,
    vds: f64,
    vdbs: f64,
    vdbd: f64,
    vsbs: f64,
    vges: f64,
    vgms: f64,
    vses: f64,
    vdes: f64,

    qb: f64,
    cqb: f64,
    qg: f64,
    cqg: f64,
    qd: f64,
    cqd: f64,
    qgmid: f64,
    cqgmid: f64,

    qbs: f64,
    cqbs: f64,
    qbd: f64,
    cqbd: f64,

    qcheq: f64,
    cqcheq: f64,
    qcdump: f64,
    cqcdump: f64,
    qdef: f64,
    qs: f64,

    von: f64,
    vdsat: f64,
    cgdo: f64,
    qgdo: f64,
    cgso: f64,
    qgso: f64,
    grbsb: f64,
    grbdb: f64,
    grbpb: f64,
    grbps: f64,
    grbpd: f64,

    Vgsteff: f64,
    vgs_eff: f64,
    vgd_eff: f64,
    dvgs_eff_dvg: f64,
    dvgd_eff_dvg: f64,
    Vdseff: f64,
    nstar: f64,
    qinv: f64,
    cd: f64,
    cbs: f64,
    cbd: f64,
    csub: f64,
    Igidl: f64,
    Igisl: f64,
    gm: f64,
    gds: f64,
    gmbs: f64,
    gbd: f64,
    gbs: f64,
    noiGd0: f64, /* tnoiMod=2 (v4.7) */
    Coxeff: f64,

    gbbs: f64,
    gbgs: f64,
    gbds: f64,
    ggidld: f64,
    ggidlg: f64,
    ggidls: f64,
    ggidlb: f64,
    ggisld: f64,
    ggislg: f64,
    ggisls: f64,
    ggislb: f64,

    Igcs: f64,
    gIgcsg: f64,
    gIgcsd: f64,
    gIgcss: f64,
    gIgcsb: f64,
    Igcd: f64,
    gIgcdg: f64,
    gIgcdd: f64,
    gIgcds: f64,
    gIgcdb: f64,

    Igs: f64,
    gIgsg: f64,
    gIgss: f64,
    Igd: f64,
    gIgdg: f64,
    gIgdd: f64,

    Igb: f64,
    gIgbg: f64,
    gIgbd: f64,
    gIgbs: f64,
    gIgbb: f64,

    grdsw: f64,
    IdovVds: f64,
    gcrg: f64,
    gcrgd: f64,
    gcrgg: f64,
    gcrgs: f64,
    gcrgb: f64,

    gstot: f64,
    gstotd: f64,
    gstotg: f64,
    gstots: f64,
    gstotb: f64,

    gdtot: f64,
    gdtotd: f64,
    gdtotg: f64,
    gdtots: f64,
    gdtotb: f64,

    cggb: f64,
    cgdb: f64,
    cgsb: f64,
    cbgb: f64,
    cbdb: f64,
    cbsb: f64,
    cdgb: f64,
    cddb: f64,
    cdsb: f64,
    csgb: f64,
    csdb: f64,
    cssb: f64,
    cgbb: f64,
    cdbb: f64,
    csbb: f64,
    cbbb: f64,
    capbd: f64,
    capbs: f64,

    cqgb: f64,
    cqdb: f64,
    cqsb: f64,
    cqbb: f64,

    qgate: f64,
    qbulk: f64,
    qdrn: f64,
    qsrc: f64,

    qchqs: f64,
    taunet: f64,
    gtau: f64,
    gtg: f64,
    gtd: f64,
    gts: f64,
    gtb: f64,

    thetavth: f64,
    ueff: f64,
    Abulk: f64,
    EsatL: f64,
    AbovVgst2Vtm: f64,
}

pub struct Bsim4SizeDepParams {
    Width: f64,
    Length: f64,
    NFinger: f64,

    cdsc: f64,
    cdscb: f64,
    cdscd: f64,
    cit: f64,
    nfactor: f64,
    xj: f64,
    vsat: f64,
    at: f64,
    a0: f64,
    ags: f64,
    a1: f64,
    a2: f64,
    keta: f64,
    nsub: f64,
    ndep: f64,
    nsd: f64,
    phin: f64,
    ngate: f64,
    gamma1: f64,
    gamma2: f64,
    vbx: f64,
    vbi: f64,
    vbm: f64,
    xt: f64,
    phi: f64,
    litl: f64,
    k1: f64,
    kt1: f64,
    kt1l: f64,
    kt2: f64,
    k2: f64,
    k3: f64,
    k3b: f64,
    w0: f64,
    dvtp0: f64,
    dvtp1: f64,
    dvtp2: f64, /* New DIBL/Rout */
    dvtp3: f64,
    dvtp4: f64,
    dvtp5: f64,
    lpe0: f64,
    lpeb: f64,
    dvt0: f64,
    dvt1: f64,
    dvt2: f64,
    dvt0w: f64,
    dvt1w: f64,
    dvt2w: f64,
    drout: f64,
    dsub: f64,
    vth0: f64,
    ua: f64,
    ua1: f64,
    ub: f64,
    ub1: f64,
    uc: f64,
    uc1: f64,
    ud: f64,
    ud1: f64,
    up: f64,
    lp: f64,
    u0: f64,
    eu: f64,
    ucs: f64,
    ute: f64,
    ucste: f64,
    voff: f64,
    tvoff: f64,
    tnfactor: f64, /* v4.7 Temp dep of leakage current */
    teta0: f64,    /* v4.7 temp dep of leakage current */
    tvoffcv: f64,  /* v4.7 temp dep of leakage current */
    minv: f64,
    minvcv: f64,
    vfb: f64,
    delta: f64,
    rdsw: f64,
    rds0: f64,
    rs0: f64,
    rd0: f64,
    rsw: f64,
    rdw: f64,
    prwg: f64,
    prwb: f64,
    prt: f64,
    eta0: f64,
    etab: f64,
    pclm: f64,
    pdibl1: f64,
    pdibl2: f64,
    pdiblb: f64,
    fprout: f64,
    pdits: f64,
    pditsd: f64,
    pscbe1: f64,
    pscbe2: f64,
    pvag: f64,
    wr: f64,
    dwg: f64,
    dwb: f64,
    b0: f64,
    b1: f64,
    alpha0: f64,
    alpha1: f64,
    beta0: f64,
    agidl: f64,
    bgidl: f64,
    cgidl: f64,
    egidl: f64,
    fgidl: f64, /* v4.7 New GIDL/GISL */
    kgidl: f64, /* v4.7 New GIDL/GISL */
    rgidl: f64, /* v4.7 New GIDL/GISL */
    agisl: f64,
    bgisl: f64,
    cgisl: f64,
    egisl: f64,
    fgisl: f64, /* v4.7 New GIDL/GISL */
    kgisl: f64, /* v4.7 New GIDL/GISL */
    rgisl: f64, /* v4.7 New GIDL/GISL */
    aigc: f64,
    bigc: f64,
    cigc: f64,
    aigsd: f64,
    bigsd: f64,
    cigsd: f64,
    aigs: f64,
    bigs: f64,
    cigs: f64,
    aigd: f64,
    bigd: f64,
    cigd: f64,
    aigbacc: f64,
    bigbacc: f64,
    cigbacc: f64,
    aigbinv: f64,
    bigbinv: f64,
    cigbinv: f64,
    nigc: f64,
    nigbacc: f64,
    nigbinv: f64,
    ntox: f64,
    eigbinv: f64,
    pigcd: f64,
    poxedge: f64,
    xrcrg1: f64,
    xrcrg2: f64,
    lambda: f64,   /* overshoot */
    vtl: f64,      /* thermal velocity limit */
    xn: f64,       /* back scattering parameter */
    lc: f64,       /* back scattering parameter */
    tfactor: f64,  /* ballistic transportation factor  */
    vfbsdoff: f64, /* S/D flatband offset voltage  */
    tvfbsdoff: f64,

    /* added for stress effect */
    ku0: f64,
    kvth0: f64,
    ku0temp: f64,
    rho_ref: f64,
    inv_od_ref: f64,
    /* added for well proximity effect */
    kvth0we: f64,
    k2we: f64,
    ku0we: f64,

    /* CV model */
    cgsl: f64,
    cgdl: f64,
    ckappas: f64,
    ckappad: f64,
    cf: f64,
    clc: f64,
    cle: f64,
    vfbcv: f64,
    noff: f64,
    voffcv: f64,
    acde: f64,
    moin: f64,

    /* Pre-calculated constants */
    dw: f64,
    dl: f64,
    leff: f64,
    weff: f64,

    dwc: f64,
    dlc: f64,
    dwj: f64,
    leffCV: f64,
    weffCV: f64,
    weffCJ: f64,
    abulkCVfactor: f64,
    cgso: f64,
    cgdo: f64,
    cgbo: f64,

    u0temp: f64,
    vsattemp: f64,
    sqrtPhi: f64,
    phis3: f64,
    Xdep0: f64,
    sqrtXdep0: f64,
    theta0vb0: f64,
    thetaRout: f64,
    mstar: f64,
    VgsteffVth: f64,
    mstarcv: f64,
    voffcbn: f64,
    voffcbncv: f64,
    rdswmin: f64,
    rdwmin: f64,
    rswmin: f64,
    vfbsd: f64,

    cof1: f64,
    cof2: f64,
    cof3: f64,
    cof4: f64,
    cdep0: f64,
    ToxRatio: f64,
    Aechvb: f64,
    Bechvb: f64,
    ToxRatioEdge: f64,
    AechvbEdgeS: f64,
    AechvbEdgeD: f64,
    BechvbEdge: f64,
    ldeb: f64,
    k1ox: f64,
    k2ox: f64,
    vfbzbfactor: f64,
    dvtp2factor: f64, /* v4.7 */
}

/// # BSIM4 Matrix Pointers
struct Bsim4MatrixPointers {
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

/// BSIM4 MOSFET Solver
pub struct Bsim4 {
    ports: Bsim4Ports,
    inst: Bsim4Inst,
    model: Bsim4Model,
    size_params: Bsim4SizeDepParams,
    intp: Bsim4InternalParams,
    guess: Bsim4OpPoint,
    op: Bsim4OpPoint,
    matps: Bsim4MatrixPointers,
}

impl Component for Bsim4 {
    fn create_matrix_elems<T: SpNum>(&mut self, mat: &mut Matrix<T>) {} // FIXME!
    fn load(&mut self, guess: &Variables<f64>, an: &AnalysisInfo) -> Stamps<f64> {
        // Start by declaring about 700 local float variables!
        let mut ceqgstot: f64;
        let mut dgstot_dvd: f64;
        let mut dgstot_dvg: f64;
        let mut dgstot_dvs: f64;
        let mut dgstot_dvb: f64;
        let mut ceqgdtot: f64;
        let mut dgdtot_dvd: f64;
        let mut dgdtot_dvg: f64;
        let mut dgdtot_dvs: f64;
        let mut dgdtot_dvb: f64;
        let mut gstot: f64;
        let mut gstotd: f64;
        let mut gstotg: f64;
        let mut gstots: f64;
        let mut gstotb: f64;
        let mut gspr: f64;
        let mut Rs: f64;
        let mut Rd: f64;
        let mut gdtot: f64;
        let mut gdtotd: f64;
        let mut gdtotg: f64;
        let mut gdtots: f64;
        let mut gdtotb: f64;
        let mut gdpr: f64;
        let mut vgs_eff: f64;
        let mut vgd_eff: f64;
        let mut dvgs_eff_dvg: f64;
        let mut dvgd_eff_dvg: f64;
        let mut dRs_dvg: f64;
        let mut dRd_dvg: f64;
        let mut dRs_dvb: f64;
        let mut dRd_dvb: f64;
        let mut dT0_dvg: f64;
        let mut dT1_dvb: f64;
        let mut dT3_dvg: f64;
        let mut dT3_dvb: f64;
        let mut vses: f64;
        let mut vdes: f64;
        let mut vdedo: f64;
        let mut delvses: f64;
        let mut delvded: f64;
        let mut delvdes: f64;
        let mut Isestot: f64;
        let mut cseshat: f64;
        let mut Idedtot: f64;
        let mut cdedhat: f64;
        let mut tol0: f64;
        let mut tol1: f64;
        let mut tol2: f64;
        let mut tol3: f64;
        let mut tol4: f64;
        let mut tol5: f64;
        let mut tol6: f64;

        let mut geltd: f64;
        let mut gcrg: f64;
        let mut gcrgg: f64;
        let mut gcrgd: f64;
        let mut gcrgs: f64;
        let mut gcrgb: f64;
        let mut ceqgcrg: f64;
        let mut vges: f64;
        let mut vgms: f64;
        let mut vgedo: f64;
        let mut vgmdo: f64;
        let mut vged: f64;
        let mut vgmd: f64;
        let mut delvged: f64;
        let mut delvgmd: f64;
        let mut delvges: f64;
        let mut delvgms: f64;
        let mut vgmb: f64;
        let mut gcgmgmb: f64;
        let mut gcgmdb: f64;
        let mut gcgmsb: f64;
        let mut gcdgmb: f64;
        let mut gcsgmb: f64;
        let mut gcgmbb: f64;
        let mut gcbgmb: f64;
        let mut qgmb: f64;

        let mut vbd: f64;
        let mut vbs: f64;
        let mut vds: f64;
        let mut vgb: f64;
        let mut vgd: f64;
        let mut vgs: f64;
        let mut vgdo: f64;
        let mut xfact: f64;
        let mut vdbs: f64;
        let mut vdbd: f64;
        let mut vsbs: f64;
        let mut vsbdo: f64;
        let mut vsbd: f64;
        let mut delvdbs: f64;
        let mut delvdbd: f64;
        let mut delvsbs: f64;
        let mut delvbd_jct: f64;
        let mut delvbs_jct: f64;
        let mut vbs_jct: f64;
        let mut vbd_jct: f64;

        let mut SourceSatCurrent: f64;
        let mut DrainSatCurrent: f64;
        let mut ag0: f64;
        let mut qgd: f64;
        let mut qgs: f64;
        let mut qgb: f64;
        let mut von: f64;
        let mut cbhat: f64;
        let mut VgstNVt: f64;
        let mut ExpVgst: f64;
        let mut ceqqb: f64;
        let mut ceqqd: f64;
        let mut ceqqg: f64;
        let mut ceqqjd: f64;
        let mut ceqqjs: f64;
        let mut ceq: f64;
        let mut geq: f64;
        let mut cdrain: f64;
        let mut cdhat: f64;
        let mut ceqdrn: f64;
        let mut ceqbd: f64;
        let mut ceqbs: f64;
        let mut ceqjd: f64;
        let mut ceqjs: f64;
        let mut gjbd: f64;
        let mut gjbs: f64;
        let mut czbd: f64;
        let mut czbdsw: f64;
        let mut czbdswg: f64;
        let mut czbs: f64;
        let mut czbssw: f64;
        let mut czbsswg: f64;
        let mut evbd: f64;
        let mut evbs: f64;
        let mut arg: f64;
        let mut sarg: f64;
        let mut delvbd: f64;
        let mut delvbs: f64;
        let mut delvds: f64;
        let mut delvgd: f64;
        let mut delvgs: f64;
        let mut Vfbeff: f64;
        let mut dVfbeff_dVg: f64;
        let mut dVfbeff_dVb: f64;
        let mut V3: f64;
        let mut V4: f64;
        let mut gcbdb: f64;
        let mut gcbgb: f64;
        let mut gcbsb: f64;
        let mut gcddb: f64;
        let mut gcdgb: f64;
        let mut gcdsb: f64;
        let mut gcgdb: f64;
        let mut gcggb: f64;
        let mut gcgsb: f64;
        let mut gcsdb: f64;
        let mut gcgbb: f64;
        let mut gcdbb: f64;
        let mut gcsbb: f64;
        let mut gcbbb: f64;
        let mut gcdbdb: f64;
        let mut gcsbsb: f64;
        let mut gcsgb: f64;
        let mut gcssb: f64;
        let mut MJD: f64;
        let mut MJSWD: f64;
        let mut MJSWGD: f64;
        let mut MJS: f64;
        let mut MJSWS: f64;
        let mut MJSWGS: f64;
        let mut cqgate: f64;
        let mut cqbody: f64;
        let mut cqdrn: f64;
        let mut Vdb: f64;
        let mut Vds: f64;
        let mut Vgs: f64;
        let mut Vbs: f64;
        let mut Gmbs: f64;
        let mut FwdSum: f64;
        let mut RevSum: f64;
        let mut Igidl: f64;
        let mut Ggidld: f64;
        let mut Ggidlg: f64;
        let mut Ggidlb: f64;
        let mut VxNVt: f64;
        let mut ExpVxNVt: f64;
        let mut Igc: f64;
        let mut dIgc_dVg: f64;
        let mut dIgc_dVd: f64;
        let mut dIgc_dVb: f64;
        let mut Igcs: f64;
        let mut dIgcs_dVg: f64;
        let mut dIgcs_dVd: f64;
        let mut dIgcs_dVb: f64;
        let mut Igcd: f64;
        let mut dIgcd_dVg: f64;
        let mut dIgcd_dVd: f64;
        let mut dIgcd_dVb: f64;
        let mut Igs: f64;
        let mut dIgs_dVg: f64;
        let mut dIgs_dVs: f64;
        let mut Igd: f64;
        let mut dIgd_dVg: f64;
        let mut dIgd_dVd: f64;
        let mut Igbacc: f64;
        let mut dIgbacc_dVg: f64;
        let mut dIgbacc_dVd: f64;
        let mut dIgbacc_dVb: f64;
        let mut Igbinv: f64;
        let mut dIgbinv_dVg: f64;
        let mut dIgbinv_dVd: f64;
        let mut dIgbinv_dVb: f64;
        let mut Igb: f64;
        let mut dIgb_dVg: f64;
        let mut dIgb_dVd: f64;
        let mut dIgb_dVb: f64;
        let mut Pigcd: f64;
        let mut dPigcd_dVg: f64;
        let mut dPigcd_dVd: f64;
        let mut dPigcd_dVb: f64;
        let mut Istoteq: f64;
        let mut gIstotg: f64;
        let mut gIstotd: f64;
        let mut gIstots: f64;
        let mut gIstotb: f64;
        let mut Idtoteq: f64;
        let mut gIdtotg: f64;
        let mut gIdtotd: f64;
        let mut gIdtots: f64;
        let mut gIdtotb: f64;
        let mut Ibtoteq: f64;
        let mut gIbtotg: f64;
        let mut gIbtotd: f64;
        let mut gIbtots: f64;
        let mut gIbtotb: f64;
        let mut Igtoteq: f64;
        let mut gIgtotg: f64;
        let mut gIgtotd: f64;
        let mut gIgtots: f64;
        let mut gIgtotb: f64;
        let mut Igstot: f64;
        let mut cgshat: f64;
        let mut Igdtot: f64;
        let mut cgdhat: f64;
        let mut Igbtot: f64;
        let mut cgbhat: f64;
        let mut Vgs_eff: f64;
        let mut dVbs_dVb: f64;
        let mut Vth_NarrowW: f64;
        let mut Phis: f64;
        let mut dPhis_dVb: f64;
        let mut sqrtPhis: f64;
        let mut dsqrtPhis_dVb: f64;
        let mut Vth: f64;
        let mut dVth_dVb: f64;
        let mut dVth_dVd: f64;
        let mut Vgst: f64;
        let mut dVgst_dVg: f64;
        let mut dVgst_dVb: f64;
        let mut dVgs_eff_dVg: f64;
        let mut Nvtms: f64;
        let mut Nvtmd: f64;
        let mut Vgdt: f64;
        let mut Vgsaddvth: f64;
        let mut Vgsaddvth2: f64;
        let mut Vgsaddvth1o3: f64;
        let mut Vtm: f64;
        let mut Vtm0: f64;
        let mut n: f64;
        let mut dn_dVb: f64;
        let mut dn_dVd: f64;
        let mut voffcv: f64;
        let mut noff: f64;
        let mut dnoff_dVd: f64;
        let mut dnoff_dVb: f64;
        let mut ExpArg: f64;
        let mut ExpArg1: f64;
        let mut V0: f64;
        let mut CoxWLcen: f64;
        let mut QovCox: f64;
        let mut LINK: f64;
        let mut DeltaPhi: f64;
        let mut dDeltaPhi_dVg: f64;
        let mut VgDP: f64;
        let mut dVgDP_dVg: f64;
        let mut Cox: f64;
        let mut Tox: f64;
        let mut Tcen: f64;
        let mut dTcen_dVg: f64;
        let mut dTcen_dVd: f64;
        let mut dTcen_dVb: f64;
        let mut Ccen: f64;
        let mut Coxeff: f64;
        let mut dCoxeff_dVd: f64;
        let mut dCoxeff_dVg: f64;
        let mut dCoxeff_dVb: f64;
        let mut Denomi: f64;
        let mut dDenomi_dVg: f64;
        let mut dDenomi_dVd: f64;
        let mut dDenomi_dVb: f64;
        let mut ueff: f64;
        let mut dueff_dVg: f64;
        let mut dueff_dVd: f64;
        let mut dueff_dVb: f64;
        let mut Esat: f64;
        let mut dEsat_dVg: f64;
        let mut dEsat_dVd: f64;
        let mut dEsat_dVb: f64;
        let mut Vdsat: f64;
        let mut Vdsat0: f64;
        let mut EsatL: f64;
        let mut dEsatL_dVg: f64;
        let mut dEsatL_dVd: f64;
        let mut dEsatL_dVb: f64;
        let mut Ilimit: f64;
        let mut Iexp: f64;
        let mut dIexp_dVg: f64;
        let mut dIexp_dVd: f64;
        let mut dIexp_dVb: f64;
        let mut dVdsat_dVg: f64;
        let mut dVdsat_dVb: f64;
        let mut dVdsat_dVd: f64;
        let mut Vasat: f64;
        let mut dAlphaz_dVg: f64;
        let mut dAlphaz_dVb: f64;
        let mut dVasat_dVg: f64;
        let mut dVasat_dVb: f64;
        let mut dVasat_dVd: f64;
        let mut Va: f64;
        let mut Va2: f64;
        let mut dVa_dVd: f64;
        let mut dVa_dVg: f64;
        let mut dVa_dVb: f64;
        let mut Vbseff: f64;
        let mut dVbseff_dVb: f64;
        let mut VbseffCV: f64;
        let mut dVbseffCV_dVb: f64;
        let mut VgsteffVth: f64;
        let mut dT11_dVg: f64;
        let mut Arg1: f64;
        let mut Arg2: f64;
        let mut One_Third_CoxWL: f64;
        let mut Two_Third_CoxWL: f64;
        let mut Alphaz: f64;
        let mut CoxWL: f64;
        let mut T0: f64;
        let mut dT0_dVg: f64;
        let mut dT0_dVd: f64;
        let mut dT0_dVb: f64;
        let mut T1: f64;
        let mut dT1_dVg: f64;
        let mut dT1_dVd: f64;
        let mut dT1_dVb: f64;
        let mut T2: f64;
        let mut dT2_dVg: f64;
        let mut dT2_dVd: f64;
        let mut dT2_dVb: f64;
        let mut T3: f64;
        let mut dT3_dVg: f64;
        let mut dT3_dVd: f64;
        let mut dT3_dVb: f64;
        let mut T4: f64;
        let mut dT4_dVg: f64;
        let mut dT4_dVd: f64;
        let mut dT4_dVb: f64;
        let mut T5: f64;
        let mut dT5_dVg: f64;
        let mut dT5_dVd: f64;
        let mut dT5_dVb: f64;
        let mut T6: f64;
        let mut dT6_dVg: f64;
        let mut dT6_dVd: f64;
        let mut dT6_dVb: f64;
        let mut T7: f64;
        let mut dT7_dVg: f64;
        let mut dT7_dVd: f64;
        let mut dT7_dVb: f64;
        let mut T8: f64;
        let mut dT8_dVg: f64;
        let mut dT8_dVd: f64;
        let mut dT8_dVb: f64;
        let mut T9: f64;
        let mut dT9_dVg: f64;
        let mut dT9_dVd: f64;
        let mut dT9_dVb: f64;
        let mut T10: f64;
        let mut dT10_dVg: f64;
        let mut dT10_dVb: f64;
        let mut dT10_dVd: f64;
        let mut T11: f64;
        let mut T12: f64;
        let mut T13: f64;
        let mut T14: f64;
        let mut tmp: f64;
        let mut Abulk: f64;
        let mut dAbulk_dVb: f64;
        let mut Abulk0: f64;
        let mut dAbulk0_dVb: f64;
        let mut Cclm: f64;
        let mut dCclm_dVg: f64;
        let mut dCclm_dVd: f64;
        let mut dCclm_dVb: f64;
        let mut FP: f64;
        let mut dFP_dVg: f64;
        let mut PvagTerm: f64;
        let mut dPvagTerm_dVg: f64;
        let mut dPvagTerm_dVd: f64;
        let mut dPvagTerm_dVb: f64;
        let mut VADITS: f64;
        let mut dVADITS_dVg: f64;
        let mut dVADITS_dVd: f64;
        let mut Lpe_Vb: f64;
        let mut DITS_Sft: f64;
        let mut dDITS_Sft_dVb: f64;
        let mut dDITS_Sft_dVd: f64;
        let mut DITS_Sft2: f64;
        let mut dDITS_Sft2_dVd: f64; /* v4.7 New DITS */
        let mut VACLM: f64;
        let mut dVACLM_dVg: f64;
        let mut dVACLM_dVd: f64;
        let mut dVACLM_dVb: f64;
        let mut VADIBL: f64;
        let mut dVADIBL_dVg: f64;
        let mut dVADIBL_dVd: f64;
        let mut dVADIBL_dVb: f64;
        let mut Xdep: f64;
        let mut dXdep_dVb: f64;
        let mut lt1: f64;
        let mut dlt1_dVb: f64;
        let mut ltw: f64;
        let mut dltw_dVb: f64;
        let mut Delt_vth: f64;
        let mut dDelt_vth_dVb: f64;
        let mut Theta0: f64;
        let mut dTheta0_dVb: f64;
        let mut Theta1: f64;
        let mut dTheta1_dVb: f64;
        let mut Thetarout: f64;
        let mut dThetarout_dVb: f64;
        let mut TempRatio: f64;
        let mut tmp1: f64;
        let mut tmp2: f64;
        let mut tmp3: f64;
        let mut tmp4: f64;
        let mut DIBL_Sft: f64;
        let mut dDIBL_Sft_dVd: f64;
        let mut DIBL_fact: f64;
        let mut Lambda: f64;
        let mut dLambda_dVg: f64;
        let mut Idtot: f64;
        let mut Ibtot: f64;
        let mut a1: f64;
        let mut ScalingFactor: f64;

        let mut Vgsteff: f64;
        let mut dVgsteff_dVg: f64;
        let mut dVgsteff_dVd: f64;
        let mut dVgsteff_dVb: f64;
        let mut Vdseff: f64;
        let mut dVdseff_dVg: f64;
        let mut dVdseff_dVd: f64;
        let mut dVdseff_dVb: f64;
        let mut VdseffCV: f64;
        let mut dVdseffCV_dVg: f64;
        let mut dVdseffCV_dVd: f64;
        let mut dVdseffCV_dVb: f64;
        let mut diffVds: f64;
        let mut diffVdsCV: f64;
        let mut dAbulk_dVg: f64;
        let mut beta: f64;
        let mut dbeta_dVg: f64;
        let mut dbeta_dVd: f64;
        let mut dbeta_dVb: f64;
        let mut gche: f64;
        let mut dgche_dVg: f64;
        let mut dgche_dVd: f64;
        let mut dgche_dVb: f64;
        let mut fgche1: f64;
        let mut dfgche1_dVg: f64;
        let mut dfgche1_dVd: f64;
        let mut dfgche1_dVb: f64;
        let mut fgche2: f64;
        let mut dfgche2_dVg: f64;
        let mut dfgche2_dVd: f64;
        let mut dfgche2_dVb: f64;
        let mut Idl: f64;
        let mut dIdl_dVg: f64;
        let mut dIdl_dVd: f64;
        let mut dIdl_dVb: f64;
        let mut Idsa: f64;
        let mut dIdsa_dVg: f64;
        let mut dIdsa_dVd: f64;
        let mut dIdsa_dVb: f64;
        let mut Ids: f64;
        let mut Gm: f64;
        let mut Gds: f64;
        let mut Gmb: f64;
        let mut devbs_dvb: f64;
        let mut devbd_dvb: f64;
        let mut Isub: f64;
        let mut Gbd: f64;
        let mut Gbg: f64;
        let mut Gbb: f64;
        let mut VASCBE: f64;
        let mut dVASCBE_dVg: f64;
        let mut dVASCBE_dVd: f64;
        let mut dVASCBE_dVb: f64;
        let mut CoxeffWovL: f64;
        let mut Rds: f64;
        let mut dRds_dVg: f64;
        let mut dRds_dVb: f64;
        let mut WVCox: f64;
        let mut WVCoxRds: f64;
        let mut Vgst2Vtm: f64;
        let mut VdsatCV: f64;
        let mut dVdsatCV_dVd: f64;
        let mut dVdsatCV_dVg: f64;
        let mut dVdsatCV_dVb: f64;
        let mut Leff: f64;
        let mut Weff: f64;
        let mut dWeff_dVg: f64;
        let mut dWeff_dVb: f64;
        let mut AbulkCV: f64;
        let mut dAbulkCV_dVb: f64;
        let mut qcheq: f64;
        let mut qdef: f64;
        let mut gqdef: f64;
        let mut cqdef: f64;
        let mut cqcheq: f64;
        // let mut gcqdb: f64;
        // let mut gcqsb: f64;
        // let mut gcqgb: f64;
        // let mut gcqbb: f64;
        // let mut dxpart: f64;
        // let mut sxpart: f64;
        // let mut ggtg: f64;
        // let mut ggtd: f64;
        // let mut ggts: f64;
        // let mut ggtb: f64;
        // let mut ddxpart_dVd: f64;
        // let mut ddxpart_dVg: f64;
        // let mut ddxpart_dVb: f64;
        // let mut ddxpart_dVs: f64;
        // let mut dsxpart_dVd: f64;
        // let mut dsxpart_dVg: f64;
        // let mut dsxpart_dVb: f64;
        // let mut dsxpart_dVs: f64;
        let mut gbspsp: f64;
        let mut gbbdp: f64;
        let mut gbbsp: f64;
        let mut gbspg: f64;
        let mut gbspb: f64;
        let mut gbspdp: f64;
        let mut gbdpdp: f64;
        let mut gbdpg: f64;
        let mut gbdpb: f64;
        let mut gbdpsp: f64;
        let mut qgdo: f64;
        let mut qgso: f64;
        let mut cgdo: f64;
        let mut cgso: f64;
        let mut cqbs: f64;
        let mut cqbd: f64;
        let mut Cgg: f64;
        let mut Cgd: f64;
        let mut Cgs: f64;
        let mut Cgb: f64;
        let mut Cdg: f64;
        let mut Cdd: f64;
        let mut Cds: f64;
        let mut Cdb: f64;
        let mut Qg: f64;
        let mut Qd: f64;
        let mut Csg: f64;
        let mut Csd: f64;
        let mut Css: f64;
        let mut Csb: f64;
        let mut Cbg: f64;
        let mut Cbd: f64;
        let mut Cbs: f64;
        let mut Cbb: f64;
        let mut Qs: f64;
        let mut Qb: f64;
        let mut Cgg1: f64;
        let mut Cgb1: f64;
        let mut Cgd1: f64;
        let mut Cbg1: f64;
        let mut Cbb1: f64;
        let mut Cbd1: f64;
        let mut Csg1: f64;
        let mut Csd1: f64;
        let mut Csb1: f64;
        let mut Qac0: f64;
        let mut Qsub0: f64;
        let mut dQac0_dVg: f64;
        let mut dQac0_dVb: f64;
        let mut dQsub0_dVg: f64;
        let mut dQsub0_dVd: f64;
        let mut dQsub0_dVb: f64;
        let mut ggidld: f64;
        let mut ggidlg: f64;
        let mut ggidlb: f64;
        let mut ggisld: f64;
        let mut ggislg: f64;
        let mut ggislb: f64;
        let mut ggisls: f64;
        let mut Igisl: f64;
        let mut Ggisld: f64;
        let mut Ggislg: f64;
        let mut Ggislb: f64;
        let mut Ggisls: f64;
        let mut Nvtmrss: f64;
        let mut Nvtmrssws: f64;
        let mut Nvtmrsswgs: f64;
        let mut Nvtmrsd: f64;
        let mut Nvtmrsswd: f64;
        let mut Nvtmrsswgd: f64;

        let mut vs: f64;
        let mut Fsevl: f64;
        let mut dvs_dVg: f64;
        let mut dvs_dVd: f64;
        let mut dvs_dVb: f64;
        let mut dFsevl_dVg: f64;
        let mut dFsevl_dVd: f64;
        let mut dFsevl_dVb: f64;
        let mut vgdx: f64;
        let mut vgsx: f64;
        let mut epssub: f64;
        let mut toxe: f64;
        let mut epsrox: f64;

        // Initialized locals. Complicated code-paths do not otherwise ensure these are ever set.
        let mut qgmid = 0.0;
        let mut ceqqgmid = 0.0;
        let mut qgate = 0.0;
        let mut qbulk = 0.0;
        let mut qdrn = 0.0;
        let mut qsrc = 0.0;

        let mut Voxdepinv = 0.0;
        let mut dVoxdepinv_dVg = 0.0;
        let mut dVoxdepinv_dVd = 0.0;
        let mut dVoxdepinv_dVb = 0.0;
        let mut dVoxacc_dVg = 0.0;
        let mut dVoxacc_dVb = 0.0;

        let mut Vaux = 0.0;
        let mut dVaux_dVg = 0.0;
        let mut dVaux_dVd = 0.0;
        let mut dVaux_dVb = 0.0;
        let mut Voxacc = 0.0;
        let mut Vfb = 0.0;

        let ScalingFactor = 1.0e-9;
        let ChargeComputationNeeded = if let AnalysisInfo::TRAN(_a, _b) = an {
            true
        } else {
            false
        };

        // Create a new operating point, which we'll fill in along the way
        let mut newop = Bsim4OpPoint::default();

        let mut vds =
            self.model.p() * (guess.get(self.ports.dNodePrime) - guess.get(self.ports.sNodePrime));
        let mut vgs =
            self.model.p() * (guess.get(self.ports.gNodePrime) - guess.get(self.ports.sNodePrime));
        let mut vbs =
            self.model.p() * (guess.get(self.ports.bNodePrime) - guess.get(self.ports.sNodePrime));
        let mut vges =
            self.model.p() * (guess.get(self.ports.gNodeExt) - guess.get(self.ports.sNodePrime));
        let mut vgms =
            self.model.p() * (guess.get(self.ports.gNodeMid) - guess.get(self.ports.sNodePrime));
        let mut vdbs =
            self.model.p() * (guess.get(self.ports.dbNode) - guess.get(self.ports.sNodePrime));
        let mut vsbs =
            self.model.p() * (guess.get(self.ports.sbNode) - guess.get(self.ports.sNodePrime));
        let mut vses =
            self.model.p() * (guess.get(self.ports.sNode) - guess.get(self.ports.sNodePrime));
        let mut vdes =
            self.model.p() * (guess.get(self.ports.dNode) - guess.get(self.ports.sNodePrime));
        let mut qdef = self.model.p() * (guess.get(self.ports.qNode));

        vgdo = self.guess.vgs - self.guess.vds;
        vgedo = self.guess.vges - self.guess.vds;
        vgmdo = self.guess.vgms - self.guess.vds;

        vbd = vbs - vds;
        vdbd = vdbs - vds;
        vgd = vgs - vds;
        vged = vges - vds;
        vgmd = vgms - vds;

        delvbd = vbd - self.guess.vbd;
        delvdbd = vdbd - self.guess.vdbd;
        delvgd = vgd - vgdo;
        delvged = vged - vgedo;
        delvgmd = vgmd - vgmdo;

        delvds = vds - self.guess.vds;
        delvgs = vgs - self.guess.vgs;
        delvges = vges - self.guess.vges;
        delvgms = vgms - self.guess.vgms;
        delvbs = vbs - self.guess.vbs;
        delvdbs = vdbs - self.guess.vdbs;
        delvsbs = vsbs - self.guess.vsbs;

        delvses = vses - (self.guess.vses);
        vdedo = self.guess.vdes - self.guess.vds;
        delvdes = vdes - self.guess.vdes;
        delvded = vdes - vds - vdedo;

        delvbd_jct = if self.intp.rbodyMod != 0 {
            delvbd
        } else {
            delvdbd
        };
        delvbs_jct = if self.intp.rbodyMod != 0 {
            delvbs
        } else {
            delvsbs
        };

        // Convergence/ LTE criteria - to be moved elsewhere
        // if self.guess.mode >= 0 {
        //     Idtot = self.guess.cd + self.guess.csub - self.guess.cbd + self.guess.Igidl;
        //     cdhat = Idtot - self.guess.gbd * delvbd_jct
        //         + (self.guess.gmbs + self.guess.gbbs + self.guess.ggidlb) * delvbs
        //         + (self.guess.gm + self.guess.gbgs + self.guess.ggidlg) * delvgs
        //         + (self.guess.gds + self.guess.gbds + self.guess.ggidld) * delvds;
        //     Ibtot = self.guess.cbs + self.guess.cbd
        //         - self.guess.Igidl
        //         - self.guess.Igisl
        //         - self.guess.csub;
        //     cbhat = Ibtot + self.guess.gbd * delvbd_jct + self.guess.gbs * delvbs_jct
        //         - (self.guess.gbbs + self.guess.ggidlb) * delvbs
        //         - (self.guess.gbgs + self.guess.ggidlg) * delvgs
        //         - (self.guess.gbds + self.guess.ggidld - self.guess.ggisls) * delvds
        //         - self.guess.ggislg * delvgd
        //         - self.guess.ggislb * delvbd;

        //     Igstot = self.guess.Igs + self.guess.Igcs;
        //     cgshat = Igstot
        //         + (self.guess.gIgsg + self.guess.gIgcsg) * delvgs
        //         + self.guess.gIgcsd * delvds
        //         + self.guess.gIgcsb * delvbs;

        //     Igdtot = self.guess.Igd + self.guess.Igcd;
        //     cgdhat = Igdtot
        //         + self.guess.gIgdg * delvgd
        //         + self.guess.gIgcdg * delvgs
        //         + self.guess.gIgcdd * delvds
        //         + self.guess.gIgcdb * delvbs;

        //     Igbtot = self.guess.Igb;
        //     cgbhat = self.guess.Igb
        //         + self.guess.gIgbg * delvgs
        //         + self.guess.gIgbd * delvds
        //         + self.guess.gIgbb * delvbs;
        // } else {
        //     Idtot = self.guess.cd + self.guess.cbd - self.guess.Igidl;
        //     cdhat = Idtot
        //         + self.guess.gbd * delvbd_jct
        //         + self.guess.gmbs * delvbd
        //         + self.guess.gm * delvgd
        //         - (self.guess.gds + self.guess.ggidls) * delvds
        //         - self.guess.ggidlg * delvgs
        //         - self.guess.ggidlb * delvbs;
        //     Ibtot = self.guess.cbs + self.guess.cbd
        //         - self.guess.Igidl
        //         - self.guess.Igisl
        //         - self.guess.csub;
        //     cbhat = Ibtot + self.guess.gbs * delvbs_jct + self.guess.gbd * delvbd_jct
        //         - (self.guess.gbbs + self.guess.ggislb) * delvbd
        //         - (self.guess.gbgs + self.guess.ggislg) * delvgd
        //         + (self.guess.gbds + self.guess.ggisld - self.guess.ggidls) * delvds
        //         - self.guess.ggidlg * delvgs
        //         - self.guess.ggidlb * delvbs;

        //     Igstot = self.guess.Igs + self.guess.Igcd;
        //     cgshat = Igstot + self.guess.gIgsg * delvgs + self.guess.gIgcdg * delvgd
        //         - self.guess.gIgcdd * delvds
        //         + self.guess.gIgcdb * delvbd;

        //     Igdtot = self.guess.Igd + self.guess.Igcs;
        //     cgdhat = Igdtot + (self.guess.gIgdg + self.guess.gIgcsg) * delvgd
        //         - self.guess.gIgcsd * delvds
        //         + self.guess.gIgcsb * delvbd;

        //     Igbtot = self.guess.Igb;
        //     cgbhat = self.guess.Igb + self.guess.gIgbg * delvgd - self.guess.gIgbd * delvds
        //         + self.guess.gIgbb * delvbd;
        // }

        Isestot = self.guess.gstot * self.guess.vses;
        cseshat = Isestot
            + self.guess.gstot * delvses
            + self.guess.gstotd * delvds
            + self.guess.gstotg * delvgs
            + self.guess.gstotb * delvbs;

        Idedtot = self.guess.gdtot * vdedo;
        cdedhat = Idedtot
            + self.guess.gdtot * delvded
            + self.guess.gdtotd * delvds
            + self.guess.gdtotg * delvgs
            + self.guess.gdtotb * delvbs;

        von = self.guess.von;
        if self.guess.vds >= 0.0 {
            vgs = DEVfetlim(vgs, self.guess.vgs, von);
            vds = vgs - vgd;
            vds = DEVlimvds(vds, self.guess.vds);
            vgd = vgs - vds;
            if self.intp.rgateMod == 3 {
                vges = DEVfetlim(vges, self.guess.vges, von);
                vgms = DEVfetlim(vgms, self.guess.vgms, von);
                vged = vges - vds;
                vgmd = vgms - vds;
            } else if (self.intp.rgateMod == 1) || (self.intp.rgateMod == 2) {
                vges = DEVfetlim(vges, self.guess.vges, von);
                vged = vges - vds;
            }
            if self.model.rdsmod != 0 {
                vdes = DEVlimvds(vdes, self.guess.vdes);
                vses = -DEVlimvds(-vses, -self.guess.vses);
            }
        } else {
            vgd = DEVfetlim(vgd, vgdo, von);
            vds = vgs - vgd;
            vds = -DEVlimvds(-vds, -self.guess.vds);
            vgs = vgd + vds;

            if self.intp.rgateMod == 3 {
                vged = DEVfetlim(vged, vgedo, von);
                vges = vged + vds;
                vgmd = DEVfetlim(vgmd, vgmdo, von);
                vgms = vgmd + vds;
            }
            if (self.intp.rgateMod == 1) || (self.intp.rgateMod == 2) {
                vged = DEVfetlim(vged, vgedo, von);
                vges = vged + vds;
            }

            if self.model.rdsmod != 0 {
                vdes = -DEVlimvds(-vdes, -self.guess.vdes);
                vses = DEVlimvds(vses, self.guess.vses);
            }
        }

        if vds >= 0.0 {
            vbs = DEVpnjlim(vbs, self.guess.vbs, VT_REF, self.intp.vcrit);
            vbd = vbs - vds;
            if self.intp.rbodyMod != 0 {
                vdbs = DEVpnjlim(vdbs, self.guess.vdbs, VT_REF, self.intp.vcrit);
                vdbd = vdbs - vds;
                vsbs = DEVpnjlim(vsbs, self.guess.vsbs, VT_REF, self.intp.vcrit);
            }
        } else {
            vbd = DEVpnjlim(vbd, self.guess.vbd, VT_REF, self.intp.vcrit);
            vbs = vbd + vds;
            if self.intp.rbodyMod != 0 {
                vdbd = DEVpnjlim(vdbd, self.guess.vdbd, VT_REF, self.intp.vcrit);
                vdbs = vdbd + vds;
                vsbdo = self.guess.vsbs - self.guess.vds;
                vsbd = vsbs - vds;
                vsbd = DEVpnjlim(vsbd, vsbdo, VT_REF, self.intp.vcrit);
                vsbs = vsbd + vds;
            }
        }

        /* Calculate DC currents and their derivatives */
        vbd = vbs - vds;
        vgd = vgs - vds;
        vgb = vgs - vbs;
        vged = vges - vds;
        vgmd = vgms - vds;
        vgmb = vgms - vbs;
        vdbd = vdbs - vds;

        vbs_jct = if self.intp.rbodyMod != 0 { vbs } else { vsbs };
        vbd_jct = if self.intp.rbodyMod != 0 { vbd } else { vdbd };

        // Source/drain junction diode DC model begins
        Nvtms = self.intp.vtm * self.intp.SjctEmissionCoeff;
        if (self.intp.Aseff <= 0.0) && (self.intp.Pseff <= 0.0) {
            SourceSatCurrent = 0.0;
        } else {
            SourceSatCurrent = self.intp.Aseff * self.intp.SjctTempSatCurDensity
                + self.intp.Pseff * self.intp.SjctSidewallTempSatCurDensity
                + self.size_params.weffCJ
                    * self.intp.nf
                    * self.intp.SjctGateSidewallTempSatCurDensity;
        }

        if SourceSatCurrent <= 0.0 {
            newop.gbs = gmin;
            newop.cbs = newop.gbs * vbs_jct;
        } else {
            match self.model.diomod {
                0 => {
                    evbs = exp(vbs_jct / Nvtms);
                    T1 = self.model.xjbvs * exp(-(self.model.bvs + vbs_jct) / Nvtms);

                    newop.gbs = SourceSatCurrent * (evbs + T1) / Nvtms + gmin;
                    newop.cbs =
                        SourceSatCurrent * (evbs + self.intp.XExpBVS - T1 - 1.0) + gmin * vbs_jct;
                }

                1 => {
                    T2 = vbs_jct / Nvtms;
                    if T2 < -EXP_THRESHOLD {
                        newop.gbs = gmin;
                        newop.cbs = SourceSatCurrent * (MIN_EXP - 1.0) + gmin * vbs_jct;
                    } else if vbs_jct <= self.intp.vjsmFwd {
                        evbs = exp(T2);
                        newop.gbs = SourceSatCurrent * evbs / Nvtms + gmin;
                        newop.cbs = SourceSatCurrent * (evbs - 1.0) + gmin * vbs_jct;
                    } else {
                        T0 = self.intp.IVjsmFwd / Nvtms;
                        newop.gbs = T0 + gmin;
                        newop.cbs = self.intp.IVjsmFwd - SourceSatCurrent
                            + T0 * (vbs_jct - self.intp.vjsmFwd)
                            + gmin * vbs_jct;
                    }
                }
                2 => {
                    if vbs_jct < self.intp.vjsmRev {
                        T0 = vbs_jct / Nvtms;
                        if T0 < -EXP_THRESHOLD {
                            evbs = MIN_EXP;
                            devbs_dvb = 0.0;
                        } else {
                            evbs = exp(T0);
                            devbs_dvb = evbs / Nvtms;
                        }

                        T1 = evbs - 1.0;
                        T2 = self.intp.IVjsmRev + self.intp.SslpRev * (vbs_jct - self.intp.vjsmRev);
                        newop.gbs = devbs_dvb * T2 + T1 * self.intp.SslpRev + gmin;
                        newop.cbs = T1 * T2 + gmin * vbs_jct;
                    } else if vbs_jct <= self.intp.vjsmFwd {
                        T0 = vbs_jct / Nvtms;
                        if T0 < -EXP_THRESHOLD {
                            evbs = MIN_EXP;
                            devbs_dvb = 0.0;
                        } else {
                            evbs = exp(T0);
                            devbs_dvb = evbs / Nvtms;
                        }

                        T1 = (self.model.bvs + vbs_jct) / Nvtms;
                        if T1 > EXP_THRESHOLD {
                            T2 = MIN_EXP;
                            T3 = 0.0;
                        } else {
                            T2 = exp(-T1);
                            T3 = -T2 / Nvtms;
                        }
                        newop.gbs = SourceSatCurrent * (devbs_dvb - self.model.xjbvs * T3) + gmin;
                        newop.cbs = SourceSatCurrent
                            * (evbs + self.intp.XExpBVS - 1.0 - self.model.xjbvs * T2)
                            + gmin * vbs_jct;
                    } else {
                        newop.gbs = self.intp.SslpFwd + gmin;
                        newop.cbs = self.intp.IVjsmFwd
                            + self.intp.SslpFwd * (vbs_jct - self.intp.vjsmFwd)
                            + gmin * vbs_jct;
                    }
                }
                _ => panic!("Invalid dioMod!"),
            }
        }

        Nvtmd = self.intp.vtm * self.intp.DjctEmissionCoeff;

        if (self.intp.Adeff <= 0.0) && (self.intp.Pdeff <= 0.0) {
            DrainSatCurrent = 0.0;
        } else {
            DrainSatCurrent = self.intp.Adeff * self.intp.DjctTempSatCurDensity
                + self.intp.Pdeff * self.intp.DjctSidewallTempSatCurDensity
                + self.size_params.weffCJ
                    * self.intp.nf
                    * self.intp.DjctGateSidewallTempSatCurDensity;
        }

        if DrainSatCurrent <= 0.0 {
            newop.gbd = gmin;
            newop.cbd = newop.gbd * vbd_jct;
        } else {
            match self.model.diomod {
                0 => {
                    evbd = exp(vbd_jct / Nvtmd);
                    T1 = self.model.xjbvd * exp(-(self.model.bvd + vbd_jct) / Nvtmd);
                    newop.gbd = DrainSatCurrent * (evbd + T1) / Nvtmd + gmin;
                    newop.cbd =
                        DrainSatCurrent * (evbd + self.intp.XExpBVD - T1 - 1.0) + gmin * vbd_jct;
                }
                1 => {
                    T2 = vbd_jct / Nvtmd;
                    if T2 < -EXP_THRESHOLD {
                        newop.gbd = gmin;
                        newop.cbd = DrainSatCurrent * (MIN_EXP - 1.0) + gmin * vbd_jct;
                    } else if vbd_jct <= self.intp.vjdmFwd {
                        evbd = exp(T2);
                        newop.gbd = DrainSatCurrent * evbd / Nvtmd + gmin;
                        newop.cbd = DrainSatCurrent * (evbd - 1.0) + gmin * vbd_jct;
                    } else {
                        T0 = self.intp.IVjdmFwd / Nvtmd;
                        newop.gbd = T0 + gmin;
                        newop.cbd = self.intp.IVjdmFwd - DrainSatCurrent
                            + T0 * (vbd_jct - self.intp.vjdmFwd)
                            + gmin * vbd_jct;
                    }
                }
                2 => {
                    if vbd_jct < self.intp.vjdmRev {
                        T0 = vbd_jct / Nvtmd;
                        if T0 < -EXP_THRESHOLD {
                            evbd = MIN_EXP;
                            devbd_dvb = 0.0;
                        } else {
                            evbd = exp(T0);
                            devbd_dvb = evbd / Nvtmd;
                        }

                        T1 = evbd - 1.0;
                        T2 = self.intp.IVjdmRev + self.intp.DslpRev * (vbd_jct - self.intp.vjdmRev);
                        newop.gbd = devbd_dvb * T2 + T1 * self.intp.DslpRev + gmin;
                        newop.cbd = T1 * T2 + gmin * vbd_jct;
                    } else if vbd_jct <= self.intp.vjdmFwd {
                        T0 = vbd_jct / Nvtmd;
                        if T0 < -EXP_THRESHOLD {
                            evbd = MIN_EXP;
                            devbd_dvb = 0.0;
                        } else {
                            evbd = exp(T0);
                            devbd_dvb = evbd / Nvtmd;
                        }

                        T1 = (self.model.bvd + vbd_jct) / Nvtmd;
                        if T1 > EXP_THRESHOLD {
                            T2 = MIN_EXP;
                            T3 = 0.0;
                        } else {
                            T2 = exp(-T1);
                            T3 = -T2 / Nvtmd;
                        }
                        newop.gbd = DrainSatCurrent * (devbd_dvb - self.model.xjbvd * T3) + gmin;
                        newop.cbd = DrainSatCurrent
                            * (evbd + self.intp.XExpBVD - 1.0 - self.model.xjbvd * T2)
                            + gmin * vbd_jct;
                    } else {
                        newop.gbd = self.intp.DslpFwd + gmin;
                        newop.cbd = self.intp.IVjdmFwd
                            + self.intp.DslpFwd * (vbd_jct - self.intp.vjdmFwd)
                            + gmin * vbd_jct;
                    }
                }
                _ => panic!("Invalid dioMod!"),
            }
        }

        /* trap-assisted tunneling and recombination current for reverse bias  */
        Nvtmrssws = self.intp.vtm0 * self.intp.njtsswstemp;
        Nvtmrsswgs = self.intp.vtm0 * self.intp.njtsswgstemp;
        Nvtmrss = self.intp.vtm0 * self.intp.njtsstemp;
        Nvtmrsswd = self.intp.vtm0 * self.intp.njtsswdtemp;
        Nvtmrsswgd = self.intp.vtm0 * self.intp.njtsswgdtemp;
        Nvtmrsd = self.intp.vtm0 * self.intp.njtsdtemp;

        if (self.model.vtss - vbs_jct) < (self.model.vtss * 1e-3) {
            T9 = 1.0e3;
            T0 = -vbs_jct / Nvtmrss * T9;
            T1 = DEXPb(T0);
            T10 = DEXPc(T0);
            dT1_dVb = T10 / Nvtmrss * T9;
        } else {
            T9 = 1.0 / (self.model.vtss - vbs_jct);
            T0 = -vbs_jct / Nvtmrss * self.model.vtss * T9;
            dT0_dVb = self.model.vtss / Nvtmrss * (T9 + vbs_jct * T9 * T9);
            T1 = DEXPb(T0);
            T10 = DEXPc(T0);
            dT1_dVb = T10 * dT0_dVb;
        }

        if (self.model.vtsd - vbd_jct) < (self.model.vtsd * 1e-3) {
            T9 = 1.0e3;
            T0 = -vbd_jct / Nvtmrsd * T9;
            T2 = DEXPb(T0);
            T10 = DEXPc(T0);
            dT2_dVb = T10 / Nvtmrsd * T9;
        } else {
            T9 = 1.0 / (self.model.vtsd - vbd_jct);
            T0 = -vbd_jct / Nvtmrsd * self.model.vtsd * T9;
            dT0_dVb = self.model.vtsd / Nvtmrsd * (T9 + vbd_jct * T9 * T9);
            T2 = DEXPb(T0);
            T10 = DEXPc(T0);
            dT2_dVb = T10 * dT0_dVb;
        }

        if (self.model.vtssws - vbs_jct) < (self.model.vtssws * 1e-3) {
            T9 = 1.0e3;
            T0 = -vbs_jct / Nvtmrssws * T9;
            T3 = DEXPb(T0);
            T10 = DEXPc(T0);
            dT3_dVb = T10 / Nvtmrssws * T9;
        } else {
            T9 = 1.0 / (self.model.vtssws - vbs_jct);
            T0 = -vbs_jct / Nvtmrssws * self.model.vtssws * T9;
            dT0_dVb = self.model.vtssws / Nvtmrssws * (T9 + vbs_jct * T9 * T9);
            T3 = DEXPb(T0);
            T10 = DEXPc(T0);
            dT3_dVb = T10 * dT0_dVb;
        }

        if (self.model.vtsswd - vbd_jct) < (self.model.vtsswd * 1e-3) {
            T9 = 1.0e3;
            T0 = -vbd_jct / Nvtmrsswd * T9;
            T4 = DEXPb(T0);
            T10 = DEXPc(T0);
            dT4_dVb = T10 / Nvtmrsswd * T9;
        } else {
            T9 = 1.0 / (self.model.vtsswd - vbd_jct);
            T0 = -vbd_jct / Nvtmrsswd * self.model.vtsswd * T9;
            dT0_dVb = self.model.vtsswd / Nvtmrsswd * (T9 + vbd_jct * T9 * T9);
            T4 = DEXPb(T0);
            T10 = DEXPc(T0);
            dT4_dVb = T10 * dT0_dVb;
        }

        if (self.model.vtsswgs - vbs_jct) < (self.model.vtsswgs * 1e-3) {
            T9 = 1.0e3;
            T0 = -vbs_jct / Nvtmrsswgs * T9;
            T5 = DEXPb(T0);
            T10 = DEXPc(T0);
            dT5_dVb = T10 / Nvtmrsswgs * T9;
        } else {
            T9 = 1.0 / (self.model.vtsswgs - vbs_jct);
            T0 = -vbs_jct / Nvtmrsswgs * self.model.vtsswgs * T9;
            dT0_dVb = self.model.vtsswgs / Nvtmrsswgs * (T9 + vbs_jct * T9 * T9);
            T5 = DEXPb(T0);
            T10 = DEXPc(T0);
            dT5_dVb = T10 * dT0_dVb;
        }

        if (self.model.vtsswgd - vbd_jct) < (self.model.vtsswgd * 1e-3) {
            T9 = 1.0e3;
            T0 = -vbd_jct / Nvtmrsswgd * T9;
            T6 = DEXPb(T0);
            T10 = DEXPc(T0);
            dT6_dVb = T10 / Nvtmrsswgd * T9;
        } else {
            T9 = 1.0 / (self.model.vtsswgd - vbd_jct);
            T0 = -vbd_jct / Nvtmrsswgd * self.model.vtsswgd * T9;
            dT0_dVb = self.model.vtsswgd / Nvtmrsswgd * (T9 + vbd_jct * T9 * T9);
            T6 = DEXPb(T0);
            T10 = DEXPc(T0);
            dT6_dVb = T10 * dT0_dVb;
        }

        newop.gbs += self.intp.SjctTempRevSatCur * dT1_dVb
            + self.intp.SswTempRevSatCur * dT3_dVb
            + self.intp.SswgTempRevSatCur * dT5_dVb;
        newop.cbs -= self.intp.SjctTempRevSatCur * (T1 - 1.0)
            + self.intp.SswTempRevSatCur * (T3 - 1.0)
            + self.intp.SswgTempRevSatCur * (T5 - 1.0);
        newop.gbd += self.intp.DjctTempRevSatCur * dT2_dVb
            + self.intp.DswTempRevSatCur * dT4_dVb
            + self.intp.DswgTempRevSatCur * dT6_dVb;
        newop.cbd -= self.intp.DjctTempRevSatCur * (T2 - 1.0)
            + self.intp.DswTempRevSatCur * (T4 - 1.0)
            + self.intp.DswgTempRevSatCur * (T6 - 1.0);

        /* End of diode DC model */

        if vds >= 0.0 {
            newop.mode = 1;
            Vds = vds;
            Vgs = vgs;
            Vbs = vbs;
            Vdb = vds - vbs;
        } else {
            newop.mode = -1;
            Vds = -vds;
            Vgs = vgd;
            Vbs = vbd;
            Vdb = -vbs;
        }

        if self.model.mtrlmod != 0 {
            epsrox = 3.9;
            toxe = self.model.eot;
            epssub = EPS0 * self.model.epsrsub;
        } else {
            epsrox = self.model.epsrox;
            toxe = self.model.toxe;
            epssub = EPSSI;
        }

        T0 = Vbs - self.intp.vbsc - 0.001;
        T1 = sqrt(T0 * T0 - 0.004 * self.intp.vbsc);
        if T0 >= 0.0 {
            Vbseff = self.intp.vbsc + 0.5 * (T0 + T1);
            dVbseff_dVb = 0.5 * (1.0 + T0 / T1);
        } else {
            T2 = -0.002 / (T1 - T0);
            Vbseff = self.intp.vbsc * (1.0 + T2);
            dVbseff_dVb = T2 * self.intp.vbsc / T1;
        }

        /* JX: Correction to forward body bias  */
        T9 = 0.95 * self.size_params.phi;
        T0 = T9 - Vbseff - 0.001;
        T1 = sqrt(T0 * T0 + 0.004 * T9);
        Vbseff = T9 - 0.5 * (T0 + T1);
        dVbseff_dVb *= 0.5 * (1.0 + T0 / T1);
        Phis = self.size_params.phi - Vbseff;
        dPhis_dVb = -1.0;
        sqrtPhis = sqrt(Phis);
        dsqrtPhis_dVb = -0.5 / sqrtPhis;

        Xdep = self.size_params.Xdep0 * sqrtPhis / self.size_params.sqrtPhi;
        dXdep_dVb = (self.size_params.Xdep0 / self.size_params.sqrtPhi) * dsqrtPhis_dVb;

        Leff = self.size_params.leff;
        Vtm = self.intp.vtm;
        Vtm0 = self.intp.vtm0;

        /* Vth Calculation */
        T3 = sqrt(Xdep);
        V0 = self.size_params.vbi - self.size_params.phi;

        T0 = self.size_params.dvt2 * Vbseff;
        if T0 >= -0.5 {
            T1 = 1.0 + T0;
            T2 = self.size_params.dvt2;
        } else {
            T4 = 1.0 / (3.0 + 8.0 * T0);
            T1 = (1.0 + 3.0 * T0) * T4;
            T2 = self.size_params.dvt2 * T4 * T4;
        }
        lt1 = self.intp.factor1 * T3 * T1;
        dlt1_dVb = self.intp.factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

        T0 = self.size_params.dvt2w * Vbseff;
        if T0 >= -0.5 {
            T1 = 1.0 + T0;
            T2 = self.size_params.dvt2w;
        } else {
            T4 = 1.0 / (3.0 + 8.0 * T0);
            T1 = (1.0 + 3.0 * T0) * T4;
            T2 = self.size_params.dvt2w * T4 * T4;
        }
        ltw = self.intp.factor1 * T3 * T1;
        dltw_dVb = self.intp.factor1 * (0.5 / T3 * T1 * dXdep_dVb + T3 * T2);

        T0 = self.size_params.dvt1 * Leff / lt1;
        if T0 < EXP_THRESHOLD {
            T1 = exp(T0);
            T2 = T1 - 1.0;
            T3 = T2 * T2;
            T4 = T3 + 2.0 * T1 * MIN_EXP;
            Theta0 = T1 / T4;
            dT1_dVb = -T0 * T1 * dlt1_dVb / lt1;
            dTheta0_dVb = dT1_dVb * (T4 - 2.0 * T1 * (T2 + MIN_EXP)) / T4 / T4;
        } else {
            Theta0 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
            dTheta0_dVb = 0.0;
        }
        newop.thetavth = self.size_params.dvt0 * Theta0;
        Delt_vth = newop.thetavth * V0;
        dDelt_vth_dVb = self.size_params.dvt0 * dTheta0_dVb * V0;

        T0 = self.size_params.dvt1w * self.size_params.weff * Leff / ltw;
        if T0 < EXP_THRESHOLD {
            T1 = exp(T0);
            T2 = T1 - 1.0;
            T3 = T2 * T2;
            T4 = T3 + 2.0 * T1 * MIN_EXP;
            T5 = T1 / T4;
            dT1_dVb = -T0 * T1 * dltw_dVb / ltw;
            dT5_dVb = dT1_dVb * (T4 - 2.0 * T1 * (T2 + MIN_EXP)) / T4 / T4;
        } else {
            T5 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
            dT5_dVb = 0.0;
        }
        T0 = self.size_params.dvt0w * T5;
        T2 = T0 * V0;
        dT2_dVb = self.size_params.dvt0w * dT5_dVb * V0;

        //   TempRatio =  ckt->CKTtemp / self.model.tnom - 1.0;
        T0 = sqrt(1.0 + self.size_params.lpe0 / Leff);
        T1 = self.size_params.k1ox * (T0 - 1.0) * self.size_params.sqrtPhi
            + (self.size_params.kt1 + self.size_params.kt1l / Leff + self.size_params.kt2 * Vbseff)
                * self.intp.TempRatio;
        Vth_NarrowW = toxe * self.size_params.phi / (self.size_params.weff + self.size_params.w0);

        T3 = self.intp.eta0 + self.size_params.etab * Vbseff;
        if T3 < 1.0e-4 {
            T9 = 1.0 / (3.0 - 2.0e4 * T3);
            T3 = (2.0e-4 - T3) * T9;
            T4 = T9 * T9;
        } else {
            T4 = 1.0;
        }
        dDIBL_Sft_dVd = T3 * self.size_params.theta0vb0;
        DIBL_Sft = dDIBL_Sft_dVd * Vds;

        Lpe_Vb = sqrt(1.0 + self.size_params.lpeb / Leff);

        Vth = self.model.p() * self.intp.vth0
            + (self.size_params.k1ox * sqrtPhis - self.size_params.k1 * self.size_params.sqrtPhi)
                * Lpe_Vb
            - self.intp.k2ox * Vbseff
            - Delt_vth
            - T2
            + (self.size_params.k3 + self.size_params.k3b * Vbseff) * Vth_NarrowW
            + T1
            - DIBL_Sft;

        dVth_dVb = Lpe_Vb * self.size_params.k1ox * dsqrtPhis_dVb
            - self.intp.k2ox
            - dDelt_vth_dVb
            - dT2_dVb
            + self.size_params.k3b * Vth_NarrowW
            - self.size_params.etab * Vds * self.size_params.theta0vb0 * T4
            + self.size_params.kt2 * self.intp.TempRatio;
        dVth_dVd = -dDIBL_Sft_dVd;

        /* Calculate n */
        tmp1 = epssub / Xdep;
        newop.nstar = self.intp.vtm / Q * (self.intp.coxe + tmp1 + self.size_params.cit);
        tmp2 = self.size_params.nfactor * tmp1;
        tmp3 =
            self.size_params.cdsc + self.size_params.cdscb * Vbseff + self.size_params.cdscd * Vds;
        tmp4 = (tmp2 + tmp3 * Theta0 + self.size_params.cit) / self.intp.coxe;
        if tmp4 >= -0.5 {
            n = 1.0 + tmp4;
            dn_dVb =
                (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb + self.size_params.cdscb * Theta0)
                    / self.intp.coxe;
            dn_dVd = self.size_params.cdscd * Theta0 / self.intp.coxe;
        } else {
            T0 = 1.0 / (3.0 + 8.0 * tmp4);
            n = (1.0 + 3.0 * tmp4) * T0;
            T0 *= T0;
            dn_dVb =
                (-tmp2 / Xdep * dXdep_dVb + tmp3 * dTheta0_dVb + self.size_params.cdscb * Theta0)
                    / self.intp.coxe
                    * T0;
            dn_dVd = self.size_params.cdscd * Theta0 / self.intp.coxe * T0;
        }

        /* Vth correction for Pocket implant */
        if self.size_params.dvtp0 > 0.0 {
            T0 = -self.size_params.dvtp1 * Vds;
            if T0 < -EXP_THRESHOLD {
                T2 = MIN_EXP;
                dT2_dVd = 0.0;
            } else {
                T2 = exp(T0);
                dT2_dVd = -self.size_params.dvtp1 * T2;
            }

            T3 = Leff + self.size_params.dvtp0 * (1.0 + T2);
            dT3_dVd = self.size_params.dvtp0 * dT2_dVd;
            if self.model.tempmod < 2 {
                T4 = Vtm * log(Leff / T3);
                dT4_dVd = -Vtm * dT3_dVd / T3;
            } else {
                T4 = self.intp.vtm0 * log(Leff / T3);
                dT4_dVd = -self.intp.vtm0 * dT3_dVd / T3;
            }
            dDITS_Sft_dVd = dn_dVd * T4 + n * dT4_dVd;
            dDITS_Sft_dVb = T4 * dn_dVb;

            Vth -= n * T4;
            dVth_dVd -= dDITS_Sft_dVd;
            dVth_dVb -= dDITS_Sft_dVb;
        }

        /* v4.7 DITS_SFT2  */
        if ((self.size_params.dvtp4 == 0.0) || (self.size_params.dvtp2factor == 0.0)) {
            T0 = 0.0;
            DITS_Sft2 = 0.0;
        } else {
            //T0 = exp(2.0 * self.size_params.dvtp4 * Vds);   /* beta code */
            T1 = 2.0 * self.size_params.dvtp4 * Vds;
            T0 = DEXPb(T1);
            T10 = DEXPc(T1);
            DITS_Sft2 = self.size_params.dvtp2factor * (T0 - 1.0) / (T0 + 1.0);
            //dDITS_Sft2_dVd = self.size_params.dvtp2factor * self.size_params.dvtp4 * 4.0 * T0 / ((T0+1) * (T0+1));   /* beta code */
            dDITS_Sft2_dVd = self.size_params.dvtp2factor * self.size_params.dvtp4 * 4.0 * T10
                / ((T0 + 1.0) * (T0 + 1.0));
            Vth -= DITS_Sft2;
            dVth_dVd -= dDITS_Sft2_dVd;
        }

        newop.von = Vth;

        /* Poly Gate Si Depletion Effect */
        T0 = self.intp.vfb + self.size_params.phi;
        if self.model.mtrlmod == 0 {
            T1 = EPSSI;
        } else {
            T1 = self.model.epsrgate * EPS0;
        }
        // Sad destructuring
        let (_v, _dv) = polyDepletion(T0, self.size_params.ngate, T1, self.intp.coxe, vgs);
        vgs_eff = _v;
        dvgs_eff_dvg = _dv;
        let (_v, _dv) = polyDepletion(T0, self.size_params.ngate, T1, self.intp.coxe, vgd);
        vgd_eff = _v;
        dvgd_eff_dvg = _dv;

        if self.guess.mode > 0 {
            Vgs_eff = vgs_eff;
            dVgs_eff_dVg = dvgs_eff_dvg;
        } else {
            Vgs_eff = vgd_eff;
            dVgs_eff_dVg = dvgd_eff_dvg;
        }
        newop.vgs_eff = vgs_eff;
        newop.vgd_eff = vgd_eff;
        newop.dvgs_eff_dvg = dvgs_eff_dvg;
        newop.dvgd_eff_dvg = dvgd_eff_dvg;

        Vgst = Vgs_eff - Vth;

        /* Calculate Vgsteff */
        T0 = n * Vtm;
        T1 = self.size_params.mstar * Vgst;
        T2 = T1 / T0;
        if T2 > EXP_THRESHOLD {
            T10 = T1;
            dT10_dVg = self.size_params.mstar * dVgs_eff_dVg;
            dT10_dVd = -dVth_dVd * self.size_params.mstar;
            dT10_dVb = -dVth_dVb * self.size_params.mstar;
        } else if T2 < -EXP_THRESHOLD {
            T10 = Vtm * log(1.0 + MIN_EXP);
            dT10_dVg = 0.0;
            dT10_dVd = T10 * dn_dVd;
            dT10_dVb = T10 * dn_dVb;
            T10 *= n;
        } else {
            ExpVgst = exp(T2);
            T3 = Vtm * log(1.0 + ExpVgst);
            T10 = n * T3;
            dT10_dVg = self.size_params.mstar * ExpVgst / (1.0 + ExpVgst);
            dT10_dVb = T3 * dn_dVb - dT10_dVg * (dVth_dVb + Vgst * dn_dVb / n);
            dT10_dVd = T3 * dn_dVd - dT10_dVg * (dVth_dVd + Vgst * dn_dVd / n);
            dT10_dVg *= dVgs_eff_dVg;
        }

        T1 = self.size_params.voffcbn - (1.0 - self.size_params.mstar) * Vgst;
        T2 = T1 / T0;
        if T2 < -EXP_THRESHOLD {
            T3 = self.intp.coxe * MIN_EXP / self.size_params.cdep0;
            T9 = self.size_params.mstar + T3 * n;
            dT9_dVg = 0.0;
            dT9_dVd = dn_dVd * T3;
            dT9_dVb = dn_dVb * T3;
        } else if T2 > EXP_THRESHOLD {
            T3 = self.intp.coxe * MAX_EXP / self.size_params.cdep0;
            T9 = self.size_params.mstar + T3 * n;
            dT9_dVg = 0.0;
            dT9_dVd = dn_dVd * T3;
            dT9_dVb = dn_dVb * T3;
        } else {
            ExpVgst = exp(T2);
            T3 = self.intp.coxe / self.size_params.cdep0;
            T4 = T3 * ExpVgst;
            T5 = T1 * T4 / T0;
            T9 = self.size_params.mstar + n * T4;
            dT9_dVg = T3 * (self.size_params.mstar - 1.0) * ExpVgst / Vtm;
            dT9_dVb = T4 * dn_dVb - dT9_dVg * dVth_dVb - T5 * dn_dVb;
            dT9_dVd = T4 * dn_dVd - dT9_dVg * dVth_dVd - T5 * dn_dVd;
            dT9_dVg *= dVgs_eff_dVg;
        }
        newop.Vgsteff = T10 / T9;
        Vgsteff = newop.Vgsteff;
        T11 = T9 * T9;
        dVgsteff_dVg = (T9 * dT10_dVg - T10 * dT9_dVg) / T11;
        dVgsteff_dVd = (T9 * dT10_dVd - T10 * dT9_dVd) / T11;
        dVgsteff_dVb = (T9 * dT10_dVb - T10 * dT9_dVb) / T11;

        /* Calculate Effective Channel Geometry */
        T9 = sqrtPhis - self.size_params.sqrtPhi;
        Weff = self.size_params.weff
            - 2.0 * (self.size_params.dwg * Vgsteff + self.size_params.dwb * T9);
        dWeff_dVg = -2.0 * self.size_params.dwg;
        dWeff_dVb = -2.0 * self.size_params.dwb * dsqrtPhis_dVb;

        if (Weff < 2.0e-8)
        /* to avoid the discontinuity problem due to Weff*/
        {
            T0 = 1.0 / (6.0e-8 - 2.0 * Weff);
            Weff = 2.0e-8 * (4.0e-8 - Weff) * T0;
            T0 *= T0 * 4.0e-16;
            dWeff_dVg *= T0;
            dWeff_dVb *= T0;
        }

        if (self.model.rdsmod == 1) {
            Rds = 0.0;
            dRds_dVg = 0.0;
            dRds_dVb = 0.0;
        } else {
            T0 = 1.0 + self.size_params.prwg * Vgsteff;
            dT0_dVg = -self.size_params.prwg / T0 / T0;
            T1 = self.size_params.prwb * T9;
            dT1_dVb = self.size_params.prwb * dsqrtPhis_dVb;

            T2 = 1.0 / T0 + T1;
            T3 = T2 + sqrt(T2 * T2 + 0.01); /* 0.01 = 4.0 * 0.05 * 0.05 */
            dT3_dVg = 1.0 + T2 / (T3 - T2);
            dT3_dVb = dT3_dVg * dT1_dVb;
            dT3_dVg *= dT0_dVg;

            T4 = self.size_params.rds0 * 0.5;
            Rds = self.size_params.rdswmin + T3 * T4;
            dRds_dVg = T4 * dT3_dVg;
            dRds_dVb = T4 * dT3_dVb;

            if (Rds > 0.0) {
                newop.grdsw = 1.0 / Rds * self.intp.nf; /*4.6.2*/
            } else {
                newop.grdsw = 0.0;
            }
        }

        /* Calculate Abulk */
        T9 = 0.5 * self.size_params.k1ox * Lpe_Vb / sqrtPhis;
        T1 = T9 + self.intp.k2ox - self.size_params.k3b * Vth_NarrowW;
        dT1_dVb = -T9 / sqrtPhis * dsqrtPhis_dVb;

        T9 = sqrt(self.size_params.xj * Xdep);
        tmp1 = Leff + 2.0 * T9;
        T5 = Leff / tmp1;
        tmp2 = self.size_params.a0 * T5;
        tmp3 = self.size_params.weff + self.size_params.b1;
        tmp4 = self.size_params.b0 / tmp3;
        T2 = tmp2 + tmp4;
        dT2_dVb = -T9 / tmp1 / Xdep * dXdep_dVb;
        T6 = T5 * T5;
        T7 = T5 * T6;

        Abulk0 = 1.0 + T1 * T2;
        dAbulk0_dVb = T1 * tmp2 * dT2_dVb + T2 * dT1_dVb;

        T8 = self.size_params.ags * self.size_params.a0 * T7;
        dAbulk_dVg = -T1 * T8;
        Abulk = Abulk0 + dAbulk_dVg * Vgsteff;
        dAbulk_dVb = dAbulk0_dVb - T8 * Vgsteff * (dT1_dVb + 3.0 * T1 * dT2_dVb);

        if (Abulk0 < 0.1)
        /* added to avoid the problems caused by Abulk0 */
        {
            T9 = 1.0 / (3.0 - 20.0 * Abulk0);
            Abulk0 = (0.2 - Abulk0) * T9;
            dAbulk0_dVb *= T9 * T9;
        }

        if Abulk < 0.1 {
            T9 = 1.0 / (3.0 - 20.0 * Abulk);
            Abulk = (0.2 - Abulk) * T9;
            T10 = T9 * T9;
            dAbulk_dVb *= T10;
            dAbulk_dVg *= T10;
        }
        newop.Abulk = Abulk;

        T2 = self.size_params.keta * Vbseff;
        if T2 >= -0.9 {
            T0 = 1.0 / (1.0 + T2);
            dT0_dVb = -self.size_params.keta * T0 * T0;
        } else {
            T1 = 1.0 / (0.8 + T2);
            T0 = (17.0 + 20.0 * T2) * T1;
            dT0_dVb = -self.size_params.keta * T1 * T1;
        }
        dAbulk_dVg *= T0;
        dAbulk_dVb = dAbulk_dVb * T0 + Abulk * dT0_dVb;
        dAbulk0_dVb = dAbulk0_dVb * T0 + Abulk0 * dT0_dVb;
        Abulk *= T0;
        Abulk0 *= T0;

        /* Mobility calculation */
        if self.model.mtrlmod != 0 && self.model.mtrlcompatmod == 0 {
            T14 = 2.0
                * self.model.p()
                * (self.model.phig - self.model.easub - 0.5 * self.intp.Eg0 + 0.45);
        } else {
            T14 = 0.0;
        }

        if self.model.mobmod == 0 {
            T0 = Vgsteff + Vth + Vth - T14;
            T2 = self.size_params.ua + self.size_params.uc * Vbseff;
            T3 = T0 / toxe;
            T12 = sqrt(Vth * Vth + 0.0001);
            T9 = 1.0 / (Vgsteff + 2.0 * T12);
            T10 = T9 * toxe;
            T8 = self.size_params.ud * T10 * T10 * Vth;
            T6 = T8 * Vth;
            T5 = T3 * (T2 + self.size_params.ub * T3) + T6;
            T7 = -2.0 * T6 * T9;
            T11 = T7 * Vth / T12;
            dDenomi_dVg = (T2 + 2.0 * self.size_params.ub * T3) / toxe;
            T13 = 2.0 * (dDenomi_dVg + T11 + T8);
            dDenomi_dVd = T13 * dVth_dVd;
            dDenomi_dVb = T13 * dVth_dVb + self.size_params.uc * T3;
            dDenomi_dVg += T7;
        } else if self.model.mobmod == 1 {
            T0 = Vgsteff + Vth + Vth - T14;
            T2 = 1.0 + self.size_params.uc * Vbseff;
            T3 = T0 / toxe;
            T4 = T3 * (self.size_params.ua + self.size_params.ub * T3);
            T12 = sqrt(Vth * Vth + 0.0001);
            T9 = 1.0 / (Vgsteff + 2.0 * T12);
            T10 = T9 * toxe;
            T8 = self.size_params.ud * T10 * T10 * Vth;
            T6 = T8 * Vth;
            T5 = T4 * T2 + T6;
            T7 = -2.0 * T6 * T9;
            T11 = T7 * Vth / T12;
            dDenomi_dVg = (self.size_params.ua + 2.0 * self.size_params.ub * T3) * T2 / toxe;
            T13 = 2.0 * (dDenomi_dVg + T11 + T8);
            dDenomi_dVd = T13 * dVth_dVd;
            dDenomi_dVb = T13 * dVth_dVb + self.size_params.uc * T4;
            dDenomi_dVg += T7;
        } else if self.model.mobmod == 2 {
            T0 = (Vgsteff + self.intp.vtfbphi1) / toxe;
            T1 = exp(self.size_params.eu * log(T0));
            dT1_dVg = T1 * self.size_params.eu / T0 / toxe;
            T2 = self.size_params.ua + self.size_params.uc * Vbseff;

            T12 = sqrt(Vth * Vth + 0.0001);
            T9 = 1.0 / (Vgsteff + 2.0 * T12);
            T10 = T9 * toxe;
            T8 = self.size_params.ud * T10 * T10 * Vth;
            T6 = T8 * Vth;
            T5 = T1 * T2 + T6;
            T7 = -2.0 * T6 * T9;
            T11 = T7 * Vth / T12;
            dDenomi_dVg = T2 * dT1_dVg + T7;
            T13 = 2.0 * (T11 + T8);
            dDenomi_dVd = T13 * dVth_dVd;
            dDenomi_dVb = T13 * dVth_dVb + T1 * self.size_params.uc;
        } else if (self.model.mobmod == 4)
        /* Synopsys 08/30/2013 add */
        {
            T0 = Vgsteff + self.intp.vtfbphi1 - T14;
            T2 = self.size_params.ua + self.size_params.uc * Vbseff;
            T3 = T0 / toxe;
            T12 = sqrt(self.intp.vtfbphi1 * self.intp.vtfbphi1 + 0.0001);
            T9 = 1.0 / (Vgsteff + 2.0 * T12);
            T10 = T9 * toxe;
            T8 = self.size_params.ud * T10 * T10 * self.intp.vtfbphi1;
            T6 = T8 * self.intp.vtfbphi1;
            T5 = T3 * (T2 + self.size_params.ub * T3) + T6;
            T7 = -2.0 * T6 * T9;
            dDenomi_dVg = (T2 + 2.0 * self.size_params.ub * T3) / toxe;
            dDenomi_dVd = 0.0;
            dDenomi_dVb = self.size_params.uc * T3;
            dDenomi_dVg += T7;
        } else if (self.model.mobmod == 5)
        /* Synopsys 08/30/2013 add */
        {
            T0 = Vgsteff + self.intp.vtfbphi1 - T14;
            T2 = 1.0 + self.size_params.uc * Vbseff;
            T3 = T0 / toxe;
            T4 = T3 * (self.size_params.ua + self.size_params.ub * T3);
            T12 = sqrt(self.intp.vtfbphi1 * self.intp.vtfbphi1 + 0.0001);
            T9 = 1.0 / (Vgsteff + 2.0 * T12);
            T10 = T9 * toxe;
            T8 = self.size_params.ud * T10 * T10 * self.intp.vtfbphi1;
            T6 = T8 * self.intp.vtfbphi1;
            T5 = T4 * T2 + T6;
            T7 = -2.0 * T6 * T9;
            dDenomi_dVg = (self.size_params.ua + 2.0 * self.size_params.ub * T3) * T2 / toxe;
            dDenomi_dVd = 0.0;
            dDenomi_dVb = self.size_params.uc * T4;
            dDenomi_dVg += T7;
        } else if (self.model.mobmod == 6)
        /* Synopsys 08/30/2013 modify */
        {
            T0 = (Vgsteff + self.intp.vtfbphi1) / toxe;
            T1 = exp(self.size_params.eu * log(T0));
            dT1_dVg = T1 * self.size_params.eu / T0 / toxe;
            T2 = self.size_params.ua + self.size_params.uc * Vbseff;

            T12 = sqrt(self.intp.vtfbphi1 * self.intp.vtfbphi1 + 0.0001);
            T9 = 1.0 / (Vgsteff + 2.0 * T12);
            T10 = T9 * toxe;
            T8 = self.size_params.ud * T10 * T10 * self.intp.vtfbphi1;
            T6 = T8 * self.intp.vtfbphi1;
            T5 = T1 * T2 + T6;
            T7 = -2.0 * T6 * T9;
            dDenomi_dVg = T2 * dT1_dVg + T7;
            dDenomi_dVd = 0.0;
            dDenomi_dVb = T1 * self.size_params.uc;
        }
        /*high K mobility*/
        else {
            /*univsersal mobility*/
            T0 = (Vgsteff + self.intp.vtfbphi1) * 1.0e-8 / toxe / 6.0;
            T1 = exp(self.size_params.eu * log(T0));
            dT1_dVg = T1 * self.size_params.eu * 1.0e-8 / T0 / toxe / 6.0;
            T2 = self.size_params.ua + self.size_params.uc * Vbseff;

            /*Coulombic*/
            VgsteffVth = self.size_params.VgsteffVth;

            T10 = exp(self.size_params.ucs * log(0.5 + 0.5 * Vgsteff / VgsteffVth));
            T11 = self.size_params.ud / T10;
            dT11_dVg =
                -0.5 * self.size_params.ucs * T11 / (0.5 + 0.5 * Vgsteff / VgsteffVth) / VgsteffVth;

            dDenomi_dVg = T2 * dT1_dVg + dT11_dVg;
            dDenomi_dVd = 0.0;
            dDenomi_dVb = T1 * self.size_params.uc;

            T5 = T1 * T2 + T11;
        }

        if T5 >= -0.8 {
            Denomi = 1.0 + T5;
        } else {
            T9 = 1.0 / (7.0 + 10.0 * T5);
            Denomi = (0.6 + T5) * T9;
            T9 *= T9;
            dDenomi_dVg *= T9;
            dDenomi_dVd *= T9;
            dDenomi_dVb *= T9;
        }

        newop.ueff = self.intp.u0temp / Denomi;
        ueff = newop.ueff;
        T9 = -ueff / Denomi;
        dueff_dVg = T9 * dDenomi_dVg;
        dueff_dVd = T9 * dDenomi_dVd;
        dueff_dVb = T9 * dDenomi_dVb;

        /* Saturation Drain Voltage  Vdsat */
        WVCox = Weff * self.intp.vsattemp * self.intp.coxe;
        WVCoxRds = WVCox * Rds;

        Esat = 2.0 * self.intp.vsattemp / ueff;
        newop.EsatL = Esat * Leff;
        EsatL = newop.EsatL;
        T0 = -EsatL / ueff;
        dEsatL_dVg = T0 * dueff_dVg;
        dEsatL_dVd = T0 * dueff_dVd;
        dEsatL_dVb = T0 * dueff_dVb;

        /* Sqrt() */
        a1 = self.size_params.a1;
        if a1 == 0.0 {
            Lambda = self.size_params.a2;
            dLambda_dVg = 0.0;
        } else if a1 > 0.0 {
            T0 = 1.0 - self.size_params.a2;
            T1 = T0 - self.size_params.a1 * Vgsteff - 0.0001;
            T2 = sqrt(T1 * T1 + 0.0004 * T0);
            Lambda = self.size_params.a2 + T0 - 0.5 * (T1 + T2);
            dLambda_dVg = 0.5 * self.size_params.a1 * (1.0 + T1 / T2);
        } else {
            T1 = self.size_params.a2 + self.size_params.a1 * Vgsteff - 0.0001;
            T2 = sqrt(T1 * T1 + 0.0004 * self.size_params.a2);
            Lambda = 0.5 * (T1 + T2);
            dLambda_dVg = 0.5 * self.size_params.a1 * (1.0 + T1 / T2);
        }

        Vgst2Vtm = Vgsteff + 2.0 * Vtm;
        if Rds > 0.0 {
            tmp2 = dRds_dVg / Rds + dWeff_dVg / Weff;
            tmp3 = dRds_dVb / Rds + dWeff_dVb / Weff;
        } else {
            tmp2 = dWeff_dVg / Weff;
            tmp3 = dWeff_dVb / Weff;
        }
        if (Rds == 0.0) && (Lambda == 1.0) {
            T0 = 1.0 / (Abulk * EsatL + Vgst2Vtm);
            tmp1 = 0.0;
            T1 = T0 * T0;
            T2 = Vgst2Vtm * T0;
            T3 = EsatL * Vgst2Vtm;
            Vdsat = T3 * T0;

            dT0_dVg = -(Abulk * dEsatL_dVg + EsatL * dAbulk_dVg + 1.0) * T1;
            dT0_dVd = -(Abulk * dEsatL_dVd) * T1;
            dT0_dVb = -(Abulk * dEsatL_dVb + dAbulk_dVb * EsatL) * T1;

            dVdsat_dVg = T3 * dT0_dVg + T2 * dEsatL_dVg + EsatL * T0;
            dVdsat_dVd = T3 * dT0_dVd + T2 * dEsatL_dVd;
            dVdsat_dVb = T3 * dT0_dVb + T2 * dEsatL_dVb;
        } else {
            tmp1 = dLambda_dVg / (Lambda * Lambda);
            T9 = Abulk * WVCoxRds;
            T8 = Abulk * T9;
            T7 = Vgst2Vtm * T9;
            T6 = Vgst2Vtm * WVCoxRds;
            T0 = 2.0 * Abulk * (T9 - 1.0 + 1.0 / Lambda);
            dT0_dVg =
                2.0 * (T8 * tmp2 - Abulk * tmp1 + (2.0 * T9 + 1.0 / Lambda - 1.0) * dAbulk_dVg);

            dT0_dVb =
                2.0 * (T8 * (2.0 / Abulk * dAbulk_dVb + tmp3) + (1.0 / Lambda - 1.0) * dAbulk_dVb);
            dT0_dVd = 0.0;
            T1 = Vgst2Vtm * (2.0 / Lambda - 1.0) + Abulk * EsatL + 3.0 * T7;

            dT1_dVg = (2.0 / Lambda - 1.0) - 2.0 * Vgst2Vtm * tmp1
                + Abulk * dEsatL_dVg
                + EsatL * dAbulk_dVg
                + 3.0 * (T9 + T7 * tmp2 + T6 * dAbulk_dVg);
            dT1_dVb = Abulk * dEsatL_dVb + EsatL * dAbulk_dVb + 3.0 * (T6 * dAbulk_dVb + T7 * tmp3);
            dT1_dVd = Abulk * dEsatL_dVd;

            T2 = Vgst2Vtm * (EsatL + 2.0 * T6);
            dT2_dVg = EsatL + Vgst2Vtm * dEsatL_dVg + T6 * (4.0 + 2.0 * Vgst2Vtm * tmp2);
            dT2_dVb = Vgst2Vtm * (dEsatL_dVb + 2.0 * T6 * tmp3);
            dT2_dVd = Vgst2Vtm * dEsatL_dVd;

            T3 = sqrt(T1 * T1 - 2.0 * T0 * T2);
            Vdsat = (T1 - T3) / T0;

            dT3_dVg = (T1 * dT1_dVg - 2.0 * (T0 * dT2_dVg + T2 * dT0_dVg)) / T3;
            dT3_dVd = (T1 * dT1_dVd - 2.0 * (T0 * dT2_dVd + T2 * dT0_dVd)) / T3;
            dT3_dVb = (T1 * dT1_dVb - 2.0 * (T0 * dT2_dVb + T2 * dT0_dVb)) / T3;

            dVdsat_dVg =
                (dT1_dVg - (T1 * dT1_dVg - dT0_dVg * T2 - T0 * dT2_dVg) / T3 - Vdsat * dT0_dVg)
                    / T0;
            dVdsat_dVb =
                (dT1_dVb - (T1 * dT1_dVb - dT0_dVb * T2 - T0 * dT2_dVb) / T3 - Vdsat * dT0_dVb)
                    / T0;
            dVdsat_dVd = (dT1_dVd - (T1 * dT1_dVd - T0 * dT2_dVd) / T3) / T0;
        }
        newop.vdsat = Vdsat;

        /* Calculate Vdseff */
        T1 = Vdsat - Vds - self.size_params.delta;
        dT1_dVg = dVdsat_dVg;
        dT1_dVd = dVdsat_dVd - 1.0;
        dT1_dVb = dVdsat_dVb;

        T2 = sqrt(T1 * T1 + 4.0 * self.size_params.delta * Vdsat);
        T0 = T1 / T2;
        T9 = 2.0 * self.size_params.delta;
        T3 = T9 / T2;
        dT2_dVg = T0 * dT1_dVg + T3 * dVdsat_dVg;
        dT2_dVd = T0 * dT1_dVd + T3 * dVdsat_dVd;
        dT2_dVb = T0 * dT1_dVb + T3 * dVdsat_dVb;

        if T1 >= 0.0 {
            Vdseff = Vdsat - 0.5 * (T1 + T2);
            dVdseff_dVg = dVdsat_dVg - 0.5 * (dT1_dVg + dT2_dVg);
            dVdseff_dVd = dVdsat_dVd - 0.5 * (dT1_dVd + dT2_dVd);
            dVdseff_dVb = dVdsat_dVb - 0.5 * (dT1_dVb + dT2_dVb);
        } else {
            T4 = T9 / (T2 - T1);
            T5 = 1.0 - T4;
            T6 = Vdsat * T4 / (T2 - T1);
            Vdseff = Vdsat * T5;
            dVdseff_dVg = dVdsat_dVg * T5 + T6 * (dT2_dVg - dT1_dVg);
            dVdseff_dVd = dVdsat_dVd * T5 + T6 * (dT2_dVd - dT1_dVd);
            dVdseff_dVb = dVdsat_dVb * T5 + T6 * (dT2_dVb - dT1_dVb);
        }

        if Vds == 0.0 {
            Vdseff = 0.0;
            dVdseff_dVg = 0.0;
            dVdseff_dVb = 0.0;
        }

        if (Vdseff > Vds) {
            Vdseff = Vds;
        }
        diffVds = Vds - Vdseff;
        newop.Vdseff = Vdseff;

        /* Velocity Overshoot */
        if self.model.lambda > 0.0 {
            T1 = Leff * ueff;
            T2 = self.size_params.lambda / T1;
            T3 = -T2 / T1 * Leff;
            dT2_dVd = T3 * dueff_dVd;
            dT2_dVg = T3 * dueff_dVg;
            dT2_dVb = T3 * dueff_dVb;
            T5 = 1.0 / (Esat * self.size_params.litl);
            T4 = -T5 / EsatL;
            dT5_dVg = dEsatL_dVg * T4;
            dT5_dVd = dEsatL_dVd * T4;
            dT5_dVb = dEsatL_dVb * T4;
            T6 = 1.0 + diffVds * T5;
            dT6_dVg = dT5_dVg * diffVds - dVdseff_dVg * T5;
            dT6_dVd = dT5_dVd * diffVds + (1.0 - dVdseff_dVd) * T5;
            dT6_dVb = dT5_dVb * diffVds - dVdseff_dVb * T5;
            T7 = 2.0 / (T6 * T6 + 1.0);
            T8 = 1.0 - T7;
            T9 = T6 * T7 * T7;
            dT8_dVg = T9 * dT6_dVg;
            dT8_dVd = T9 * dT6_dVd;
            dT8_dVb = T9 * dT6_dVb;
            T10 = 1.0 + T2 * T8;
            dT10_dVg = dT2_dVg * T8 + T2 * dT8_dVg;
            dT10_dVd = dT2_dVd * T8 + T2 * dT8_dVd;
            dT10_dVb = dT2_dVb * T8 + T2 * dT8_dVb;
            if T10 == 1.0 {
                dT10_dVg = 0.0;
                dT10_dVd = 0.0;
                dT10_dVb = 0.0;
            }

            dEsatL_dVg *= T10;
            dEsatL_dVg += EsatL * dT10_dVg;
            dEsatL_dVd *= T10;
            dEsatL_dVd += EsatL * dT10_dVd;
            dEsatL_dVb *= T10;
            dEsatL_dVb += EsatL * dT10_dVb;
            EsatL *= T10;
            Esat = EsatL / Leff; /* bugfix by Wenwei Yang (4.6.4) */
            newop.EsatL = EsatL;
        }

        /* Calculate Vasat */
        tmp4 = 1.0 - 0.5 * Abulk * Vdsat / Vgst2Vtm;
        T9 = WVCoxRds * Vgsteff;
        T8 = T9 / Vgst2Vtm;
        T0 = EsatL + Vdsat + 2.0 * T9 * tmp4;

        T7 = 2.0 * WVCoxRds * tmp4;
        dT0_dVg = dEsatL_dVg + dVdsat_dVg + T7 * (1.0 + tmp2 * Vgsteff)
            - T8 * (Abulk * dVdsat_dVg - Abulk * Vdsat / Vgst2Vtm + Vdsat * dAbulk_dVg);

        dT0_dVb = dEsatL_dVb + dVdsat_dVb + T7 * tmp3 * Vgsteff
            - T8 * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
        dT0_dVd = dEsatL_dVd + dVdsat_dVd - T8 * Abulk * dVdsat_dVd;

        T9 = WVCoxRds * Abulk;
        T1 = 2.0 / Lambda - 1.0 + T9;
        dT1_dVg = -2.0 * tmp1 + WVCoxRds * (Abulk * tmp2 + dAbulk_dVg);
        dT1_dVb = dAbulk_dVb * WVCoxRds + T9 * tmp3;

        Vasat = T0 / T1;
        dVasat_dVg = (dT0_dVg - Vasat * dT1_dVg) / T1;
        dVasat_dVb = (dT0_dVb - Vasat * dT1_dVb) / T1;
        dVasat_dVd = dT0_dVd / T1;

        /* Calculate Idl first */

        tmp1 = self.intp.vtfbphi2;
        tmp2 = 2.0e8 * self.intp.toxp;
        dT0_dVg = 1.0 / tmp2;
        T0 = (Vgsteff + tmp1) * dT0_dVg;

        tmp3 = exp(self.model.bdos * 0.7 * log(T0));
        T1 = 1.0 + tmp3;
        T2 = self.model.bdos * 0.7 * tmp3 / T0;
        Tcen = self.model.ados * 1.9e-9 / T1;
        dTcen_dVg = -Tcen * T2 * dT0_dVg / T1;

        Coxeff = epssub * self.intp.coxp / (epssub + self.intp.coxp * Tcen);
        newop.Coxeff = Coxeff;
        dCoxeff_dVg = -Coxeff * Coxeff * dTcen_dVg / epssub;

        CoxeffWovL = Coxeff * Weff / Leff;
        beta = ueff * CoxeffWovL;
        T3 = ueff / Leff;
        dbeta_dVg = CoxeffWovL * dueff_dVg + T3 * (Weff * dCoxeff_dVg + Coxeff * dWeff_dVg);
        dbeta_dVd = CoxeffWovL * dueff_dVd;
        dbeta_dVb = CoxeffWovL * dueff_dVb + T3 * Coxeff * dWeff_dVb;

        newop.AbovVgst2Vtm = Abulk / Vgst2Vtm;
        T0 = 1.0 - 0.5 * Vdseff * newop.AbovVgst2Vtm;
        dT0_dVg = -0.5 * (Abulk * dVdseff_dVg - Abulk * Vdseff / Vgst2Vtm + Vdseff * dAbulk_dVg)
            / Vgst2Vtm;
        dT0_dVd = -0.5 * Abulk * dVdseff_dVd / Vgst2Vtm;
        dT0_dVb = -0.5 * (Abulk * dVdseff_dVb + dAbulk_dVb * Vdseff) / Vgst2Vtm;

        fgche1 = Vgsteff * T0;
        dfgche1_dVg = Vgsteff * dT0_dVg + T0;
        dfgche1_dVd = Vgsteff * dT0_dVd;
        dfgche1_dVb = Vgsteff * dT0_dVb;

        T9 = Vdseff / EsatL;
        fgche2 = 1.0 + T9;
        dfgche2_dVg = (dVdseff_dVg - T9 * dEsatL_dVg) / EsatL;
        dfgche2_dVd = (dVdseff_dVd - T9 * dEsatL_dVd) / EsatL;
        dfgche2_dVb = (dVdseff_dVb - T9 * dEsatL_dVb) / EsatL;

        gche = beta * fgche1 / fgche2;
        dgche_dVg = (beta * dfgche1_dVg + fgche1 * dbeta_dVg - gche * dfgche2_dVg) / fgche2;
        dgche_dVd = (beta * dfgche1_dVd + fgche1 * dbeta_dVd - gche * dfgche2_dVd) / fgche2;
        dgche_dVb = (beta * dfgche1_dVb + fgche1 * dbeta_dVb - gche * dfgche2_dVb) / fgche2;

        T0 = 1.0 + gche * Rds;
        Idl = gche / T0;
        T1 = (1.0 - Idl * Rds) / T0;
        T2 = Idl * Idl;
        dIdl_dVg = T1 * dgche_dVg - T2 * dRds_dVg;
        dIdl_dVd = T1 * dgche_dVd;
        dIdl_dVb = T1 * dgche_dVb - T2 * dRds_dVb;

        /* Calculate degradation factor due to pocket implant */

        if self.size_params.fprout <= 0.0 {
            FP = 1.0;
            dFP_dVg = 0.0;
        } else {
            T9 = self.size_params.fprout * sqrt(Leff) / Vgst2Vtm;
            FP = 1.0 / (1.0 + T9);
            dFP_dVg = FP * FP * T9 / Vgst2Vtm;
        }

        /* Calculate VACLM */
        T8 = self.size_params.pvag / EsatL;
        T9 = T8 * Vgsteff;
        if T9 > -0.9 {
            PvagTerm = 1.0 + T9;
            dPvagTerm_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL);
            dPvagTerm_dVb = -T9 * dEsatL_dVb / EsatL;
            dPvagTerm_dVd = -T9 * dEsatL_dVd / EsatL;
        } else {
            T4 = 1.0 / (17.0 + 20.0 * T9);
            PvagTerm = (0.8 + T9) * T4;
            T4 *= T4;
            dPvagTerm_dVg = T8 * (1.0 - Vgsteff * dEsatL_dVg / EsatL) * T4;
            T9 *= T4 / EsatL;
            dPvagTerm_dVb = -T9 * dEsatL_dVb;
            dPvagTerm_dVd = -T9 * dEsatL_dVd;
        }

        if (self.size_params.pclm > MIN_EXP) && (diffVds > 1.0e-10) {
            T0 = 1.0 + Rds * Idl;
            dT0_dVg = dRds_dVg * Idl + Rds * dIdl_dVg;
            dT0_dVd = Rds * dIdl_dVd;
            dT0_dVb = dRds_dVb * Idl + Rds * dIdl_dVb;

            T2 = Vdsat / Esat;
            T1 = Leff + T2;
            dT1_dVg = (dVdsat_dVg - T2 * dEsatL_dVg / Leff) / Esat;
            dT1_dVd = (dVdsat_dVd - T2 * dEsatL_dVd / Leff) / Esat;
            dT1_dVb = (dVdsat_dVb - T2 * dEsatL_dVb / Leff) / Esat;

            Cclm = FP * PvagTerm * T0 * T1 / (self.size_params.pclm * self.size_params.litl);
            dCclm_dVg =
                Cclm * (dFP_dVg / FP + dPvagTerm_dVg / PvagTerm + dT0_dVg / T0 + dT1_dVg / T1);
            dCclm_dVb = Cclm * (dPvagTerm_dVb / PvagTerm + dT0_dVb / T0 + dT1_dVb / T1);
            dCclm_dVd = Cclm * (dPvagTerm_dVd / PvagTerm + dT0_dVd / T0 + dT1_dVd / T1);
            VACLM = Cclm * diffVds;

            dVACLM_dVg = dCclm_dVg * diffVds - dVdseff_dVg * Cclm;
            dVACLM_dVb = dCclm_dVb * diffVds - dVdseff_dVb * Cclm;
            dVACLM_dVd = dCclm_dVd * diffVds + (1.0 - dVdseff_dVd) * Cclm;
        } else {
            VACLM = MAX_EXP;
            Cclm = MAX_EXP;
            dVACLM_dVd = 0.0;
            dVACLM_dVg = 0.0;
            dVACLM_dVb = 0.0;
            dCclm_dVd = 0.0;
            dCclm_dVg = 0.0;
            dCclm_dVb = 0.0;
        }

        /* Calculate VADIBL */
        if self.size_params.thetaRout > MIN_EXP {
            T8 = Abulk * Vdsat;
            T0 = Vgst2Vtm * T8;
            dT0_dVg = Vgst2Vtm * Abulk * dVdsat_dVg + T8 + Vgst2Vtm * Vdsat * dAbulk_dVg;
            dT0_dVb = Vgst2Vtm * (dAbulk_dVb * Vdsat + Abulk * dVdsat_dVb);
            dT0_dVd = Vgst2Vtm * Abulk * dVdsat_dVd;

            T1 = Vgst2Vtm + T8;
            dT1_dVg = 1.0 + Abulk * dVdsat_dVg + Vdsat * dAbulk_dVg;
            dT1_dVb = Abulk * dVdsat_dVb + dAbulk_dVb * Vdsat;
            dT1_dVd = Abulk * dVdsat_dVd;

            T9 = T1 * T1;
            T2 = self.size_params.thetaRout;
            VADIBL = (Vgst2Vtm - T0 / T1) / T2;
            dVADIBL_dVg = (1.0 - dT0_dVg / T1 + T0 * dT1_dVg / T9) / T2;
            dVADIBL_dVb = (-dT0_dVb / T1 + T0 * dT1_dVb / T9) / T2;
            dVADIBL_dVd = (-dT0_dVd / T1 + T0 * dT1_dVd / T9) / T2;

            T7 = self.size_params.pdiblb * Vbseff;
            if T7 >= -0.9 {
                T3 = 1.0 / (1.0 + T7);
                VADIBL *= T3;
                dVADIBL_dVg *= T3;
                dVADIBL_dVb = (dVADIBL_dVb - VADIBL * self.size_params.pdiblb) * T3;
                dVADIBL_dVd *= T3;
            } else {
                T4 = 1.0 / (0.8 + T7);
                T3 = (17.0 + 20.0 * T7) * T4;
                dVADIBL_dVg *= T3;
                dVADIBL_dVb = dVADIBL_dVb * T3 - VADIBL * self.size_params.pdiblb * T4 * T4;
                dVADIBL_dVd *= T3;
                VADIBL *= T3;
            }

            dVADIBL_dVg = dVADIBL_dVg * PvagTerm + VADIBL * dPvagTerm_dVg;
            dVADIBL_dVb = dVADIBL_dVb * PvagTerm + VADIBL * dPvagTerm_dVb;
            dVADIBL_dVd = dVADIBL_dVd * PvagTerm + VADIBL * dPvagTerm_dVd;
            VADIBL *= PvagTerm;
        } else {
            VADIBL = MAX_EXP;
            dVADIBL_dVd = 0.0;
            dVADIBL_dVg = 0.0;
            dVADIBL_dVb = 0.0;
        }

        /* Calculate Va */
        Va = Vasat + VACLM;
        dVa_dVg = dVasat_dVg + dVACLM_dVg;
        dVa_dVb = dVasat_dVb + dVACLM_dVb;
        dVa_dVd = dVasat_dVd + dVACLM_dVd;

        /* Calculate VADITS */
        T0 = self.size_params.pditsd * Vds;
        if T0 > EXP_THRESHOLD {
            T1 = MAX_EXP;
            dT1_dVd = 0.0;
        } else {
            T1 = exp(T0);
            dT1_dVd = T1 * self.size_params.pditsd;
        }

        if self.size_params.pdits > MIN_EXP {
            T2 = 1.0 + self.model.pditsl * Leff;
            VADITS = (1.0 + T2 * T1) / self.size_params.pdits;
            dVADITS_dVg = VADITS * dFP_dVg;
            dVADITS_dVd = FP * T2 * dT1_dVd / self.size_params.pdits;
            VADITS *= FP;
        } else {
            VADITS = MAX_EXP;
            dVADITS_dVg = 0.0;
            dVADITS_dVd = 0.0;
        }

        /* Calculate VASCBE */
        if ((self.size_params.pscbe2 > 0.0) && (self.size_params.pscbe1 >= 0.0))
        /*4.6.2*/
        {
            if (diffVds > self.size_params.pscbe1 * self.size_params.litl / EXP_THRESHOLD) {
                T0 = self.size_params.pscbe1 * self.size_params.litl / diffVds;
                VASCBE = Leff * exp(T0) / self.size_params.pscbe2;
                T1 = T0 * VASCBE / diffVds;
                dVASCBE_dVg = T1 * dVdseff_dVg;
                dVASCBE_dVd = -T1 * (1.0 - dVdseff_dVd);
                dVASCBE_dVb = T1 * dVdseff_dVb;
            } else {
                VASCBE = MAX_EXP * Leff / self.size_params.pscbe2;
                dVASCBE_dVg = 0.0;
                dVASCBE_dVd = 0.0;
                dVASCBE_dVb = 0.0;
            }
        } else {
            VASCBE = MAX_EXP;
            dVASCBE_dVg = 0.0;
            dVASCBE_dVd = 0.0;
            dVASCBE_dVb = 0.0;
        }

        /* Add DIBL to Ids */
        T9 = diffVds / VADIBL;
        T0 = 1.0 + T9;
        Idsa = Idl * T0;
        dIdsa_dVg = T0 * dIdl_dVg - Idl * (dVdseff_dVg + T9 * dVADIBL_dVg) / VADIBL;
        dIdsa_dVd = T0 * dIdl_dVd + Idl * (1.0 - dVdseff_dVd - T9 * dVADIBL_dVd) / VADIBL;
        dIdsa_dVb = T0 * dIdl_dVb - Idl * (dVdseff_dVb + T9 * dVADIBL_dVb) / VADIBL;

        /* Add DITS to Ids */
        T9 = diffVds / VADITS;
        T0 = 1.0 + T9;
        dIdsa_dVg = T0 * dIdsa_dVg - Idsa * (dVdseff_dVg + T9 * dVADITS_dVg) / VADITS;
        dIdsa_dVd = T0 * dIdsa_dVd + Idsa * (1.0 - dVdseff_dVd - T9 * dVADITS_dVd) / VADITS;
        dIdsa_dVb = T0 * dIdsa_dVb - Idsa * dVdseff_dVb / VADITS;
        Idsa *= T0;

        /* Add CLM to Ids */
        T0 = log(Va / Vasat);
        dT0_dVg = dVa_dVg / Va - dVasat_dVg / Vasat;
        dT0_dVb = dVa_dVb / Va - dVasat_dVb / Vasat;
        dT0_dVd = dVa_dVd / Va - dVasat_dVd / Vasat;
        T1 = T0 / Cclm;
        T9 = 1.0 + T1;
        dT9_dVg = (dT0_dVg - T1 * dCclm_dVg) / Cclm;
        dT9_dVb = (dT0_dVb - T1 * dCclm_dVb) / Cclm;
        dT9_dVd = (dT0_dVd - T1 * dCclm_dVd) / Cclm;

        dIdsa_dVg = dIdsa_dVg * T9 + Idsa * dT9_dVg;
        dIdsa_dVb = dIdsa_dVb * T9 + Idsa * dT9_dVb;
        dIdsa_dVd = dIdsa_dVd * T9 + Idsa * dT9_dVd;
        Idsa *= T9;

        /* Substrate current begins */
        tmp = self.size_params.alpha0 + self.size_params.alpha1 * Leff;
        if (tmp <= 0.0) || (self.size_params.beta0 <= 0.0) {
            Isub = 0.0;
            Gbd = 0.0;
            Gbb = 0.0;
            Gbg = 0.0;
        } else {
            T2 = tmp / Leff;
            if diffVds > self.size_params.beta0 / EXP_THRESHOLD {
                T0 = -self.size_params.beta0 / diffVds;
                T1 = T2 * diffVds * exp(T0);
                T3 = T1 / diffVds * (T0 - 1.0);
                dT1_dVg = T3 * dVdseff_dVg;
                dT1_dVd = T3 * (dVdseff_dVd - 1.0);
                dT1_dVb = T3 * dVdseff_dVb;
            } else {
                T3 = T2 * MIN_EXP;
                T1 = T3 * diffVds;
                dT1_dVg = -T3 * dVdseff_dVg;
                dT1_dVd = T3 * (1.0 - dVdseff_dVd);
                dT1_dVb = -T3 * dVdseff_dVb;
            }
            T4 = Idsa * Vdseff;
            Isub = T1 * T4;
            Gbg = T1 * (dIdsa_dVg * Vdseff + Idsa * dVdseff_dVg) + T4 * dT1_dVg;
            Gbd = T1 * (dIdsa_dVd * Vdseff + Idsa * dVdseff_dVd) + T4 * dT1_dVd;
            Gbb = T1 * (dIdsa_dVb * Vdseff + Idsa * dVdseff_dVb) + T4 * dT1_dVb;

            Gbd += Gbg * dVgsteff_dVd;
            Gbb += Gbg * dVgsteff_dVb;
            Gbg *= dVgsteff_dVg;
            Gbb *= dVbseff_dVb;
        }
        newop.csub = Isub;
        newop.gbbs = Gbb;
        newop.gbgs = Gbg;
        newop.gbds = Gbd;

        /* Add SCBE to Ids */
        T9 = diffVds / VASCBE;
        T0 = 1.0 + T9;
        Ids = Idsa * T0;

        Gm = T0 * dIdsa_dVg - Idsa * (dVdseff_dVg + T9 * dVASCBE_dVg) / VASCBE;
        Gds = T0 * dIdsa_dVd + Idsa * (1.0 - dVdseff_dVd - T9 * dVASCBE_dVd) / VASCBE;
        Gmb = T0 * dIdsa_dVb - Idsa * (dVdseff_dVb + T9 * dVASCBE_dVb) / VASCBE;

        tmp1 = Gds + Gm * dVgsteff_dVd;
        tmp2 = Gmb + Gm * dVgsteff_dVb;
        tmp3 = Gm;

        Gm = (Ids * dVdseff_dVg + Vdseff * tmp3) * dVgsteff_dVg;
        Gds = Ids * (dVdseff_dVd + dVdseff_dVg * dVgsteff_dVd) + Vdseff * tmp1;
        Gmb = (Ids * (dVdseff_dVb + dVdseff_dVg * dVgsteff_dVb) + Vdseff * tmp2) * dVbseff_dVb;

        cdrain = Ids * Vdseff;

        /* Source End Velocity Limit  */
        // if ((self.model.vtlGiven) && (self.model.vtl > 0.0)) {
        if self.model.vtl > 0.0 {
            // FIXME: the reference implementation's default condition here is "not given",
            // (although with a default value of 2e5)
            // So far we default to zero, in which case this block is not executed.

            T12 = 1.0 / Leff / CoxeffWovL;
            T11 = T12 / Vgsteff;
            T10 = -T11 / Vgsteff;
            vs = cdrain * T11; /* vs */
            dvs_dVg = Gm * T11 + cdrain * T10 * dVgsteff_dVg;
            dvs_dVd = Gds * T11 + cdrain * T10 * dVgsteff_dVd;
            dvs_dVb = Gmb * T11 + cdrain * T10 * dVgsteff_dVb;
            T0 = 2.0 * MM;
            T1 = vs / (self.size_params.vtl * self.size_params.tfactor);
            if T1 > 0.0 {
                T2 = 1.0 + exp(T0 * log(T1));
                T3 = (T2 - 1.0) * T0 / vs;
                Fsevl = 1.0 / exp(log(T2) / T0);
                dT2_dVg = T3 * dvs_dVg;
                dT2_dVd = T3 * dvs_dVd;
                dT2_dVb = T3 * dvs_dVb;
                T4 = -1.0 / T0 * Fsevl / T2;
                dFsevl_dVg = T4 * dT2_dVg;
                dFsevl_dVd = T4 * dT2_dVd;
                dFsevl_dVb = T4 * dT2_dVb;
            } else {
                Fsevl = 1.0;
                dFsevl_dVg = 0.0;
                dFsevl_dVd = 0.0;
                dFsevl_dVb = 0.0;
            }
            Gm *= Fsevl;
            Gm += cdrain * dFsevl_dVg;
            Gmb *= Fsevl;
            Gmb += cdrain * dFsevl_dVb;
            Gds *= Fsevl;
            Gds += cdrain * dFsevl_dVd;

            cdrain *= Fsevl;
        }

        newop.gds = Gds;
        newop.gm = Gm;
        newop.gmbs = Gmb;
        newop.IdovVds = Ids;
        if (newop.IdovVds <= 1.0e-9) {
            newop.IdovVds = 1.0e-9;
        }

        /* Calculate Rg */
        if ((self.intp.rgateMod > 1) || (self.model.trnqsmod != 0) || (self.model.acnqsmod != 0)) {
            T9 = self.size_params.xrcrg2 * self.intp.vtm;
            T0 = T9 * beta;
            dT0_dVd = (dbeta_dVd + dbeta_dVg * dVgsteff_dVd) * T9;
            dT0_dVb = (dbeta_dVb + dbeta_dVg * dVgsteff_dVb) * T9;
            dT0_dVg = dbeta_dVg * T9;

            newop.gcrg = self.size_params.xrcrg1 * (T0 + Ids);
            newop.gcrgd = self.size_params.xrcrg1 * (dT0_dVd + tmp1);
            newop.gcrgb = self.size_params.xrcrg1 * (dT0_dVb + tmp2) * dVbseff_dVb;
            newop.gcrgg = self.size_params.xrcrg1 * (dT0_dVg + tmp3) * dVgsteff_dVg;

            if self.intp.nf != 1.0 {
                newop.gcrg *= self.intp.nf;
                newop.gcrgg *= self.intp.nf;
                newop.gcrgd *= self.intp.nf;
                newop.gcrgb *= self.intp.nf;
            }

            if self.intp.rgateMod == 2 {
                T10 = self.intp.grgeltd * self.intp.grgeltd;
                T11 = self.intp.grgeltd + newop.gcrg;
                newop.gcrg = self.intp.grgeltd * newop.gcrg / T11;
                T12 = T10 / T11 / T11;
                newop.gcrgg *= T12;
                newop.gcrgd *= T12;
                newop.gcrgb *= T12;
            }
            newop.gcrgs = -(newop.gcrgg + newop.gcrgd + newop.gcrgb);
        }

        /* Calculate bias-dependent external S/D resistance */
        if self.model.rdsmod != 0 {
            /* Rs(V) */
            T0 = vgs - self.size_params.vfbsd;
            T1 = sqrt(T0 * T0 + 1.0e-4);
            vgs_eff = 0.5 * (T0 + T1);
            dvgs_eff_dvg = vgs_eff / T1;

            T0 = 1.0 + self.size_params.prwg * vgs_eff;
            dT0_dvg = -self.size_params.prwg / T0 / T0 * dvgs_eff_dvg;
            T1 = -self.size_params.prwb * vbs;
            dT1_dvb = -self.size_params.prwb;

            T2 = 1.0 / T0 + T1;
            T3 = T2 + sqrt(T2 * T2 + 0.01);
            dT3_dvg = T3 / (T3 - T2);
            dT3_dvb = dT3_dvg * dT1_dvb;
            dT3_dvg *= dT0_dvg;

            T4 = self.size_params.rs0 * 0.5;
            Rs = self.size_params.rswmin + T3 * T4;
            dRs_dvg = T4 * dT3_dvg;
            dRs_dvb = T4 * dT3_dvb;

            T0 = 1.0 + self.intp.sourceConductance * Rs;
            newop.gstot = self.intp.sourceConductance / T0;
            T0 = -newop.gstot * newop.gstot;
            dgstot_dvd = 0.0;
            dgstot_dvg = T0 * dRs_dvg;
            dgstot_dvb = T0 * dRs_dvb;
            dgstot_dvs = -(dgstot_dvg + dgstot_dvb + dgstot_dvd);

            /* Rd(V) */
            T0 = vgd - self.size_params.vfbsd;
            T1 = sqrt(T0 * T0 + 1.0e-4);
            vgd_eff = 0.5 * (T0 + T1);
            dvgd_eff_dvg = vgd_eff / T1;

            T0 = 1.0 + self.size_params.prwg * vgd_eff;
            dT0_dvg = -self.size_params.prwg / T0 / T0 * dvgd_eff_dvg;
            T1 = -self.size_params.prwb * vbd;
            dT1_dvb = -self.size_params.prwb;

            T2 = 1.0 / T0 + T1;
            T3 = T2 + sqrt(T2 * T2 + 0.01);
            dT3_dvg = T3 / (T3 - T2);
            dT3_dvb = dT3_dvg * dT1_dvb;
            dT3_dvg *= dT0_dvg;

            T4 = self.size_params.rd0 * 0.5;
            Rd = self.size_params.rdwmin + T3 * T4;
            dRd_dvg = T4 * dT3_dvg;
            dRd_dvb = T4 * dT3_dvb;

            T0 = 1.0 + self.intp.drainConductance * Rd;
            newop.gdtot = self.intp.drainConductance / T0;
            T0 = -newop.gdtot * newop.gdtot;
            dgdtot_dvs = 0.0;
            dgdtot_dvg = T0 * dRd_dvg;
            dgdtot_dvb = T0 * dRd_dvb;
            dgdtot_dvd = -(dgdtot_dvg + dgdtot_dvb + dgdtot_dvs);

            newop.gstotd = vses * dgstot_dvd;
            newop.gstotg = vses * dgstot_dvg;
            newop.gstots = vses * dgstot_dvs;
            newop.gstotb = vses * dgstot_dvb;

            T2 = vdes - vds;
            newop.gdtotd = T2 * dgdtot_dvd;
            newop.gdtotg = T2 * dgdtot_dvg;
            newop.gdtots = T2 * dgdtot_dvs;
            newop.gdtotb = T2 * dgdtot_dvb;
        } else
        /* WDLiu: for bypass */
        {
            newop.gstot = 0.0;
            newop.gstotd = 0.0;
            newop.gstotg = 0.0;
            newop.gstots = 0.0;
            newop.gstotb = 0.0;
            newop.gdtot = 0.0;
            newop.gdtotd = 0.0;
            newop.gdtotg = 0.0;
            newop.gdtots = 0.0;
            newop.gdtotb = 0.0;
        }

        /* GIDL/GISL Models */

        if (self.model.mtrlmod == 0) {
            T0 = 3.0 * toxe;
        } else {
            T0 = self.model.epsrsub * toxe / epsrox;
        }

        /* Calculate GIDL current */

        vgs_eff = newop.vgs_eff;
        dvgs_eff_dvg = newop.dvgs_eff_dvg;
        vgd_eff = newop.vgd_eff;
        dvgd_eff_dvg = newop.dvgd_eff_dvg;

        if (self.model.gidlmod == 0) {
            if (self.model.mtrlmod == 0) {
                T1 = (vds - vgs_eff - self.size_params.egidl) / T0;
            } else {
                T1 = (vds - vgs_eff - self.size_params.egidl + self.size_params.vfbsd) / T0;
            }

            if ((self.size_params.agidl <= 0.0)
                || (self.size_params.bgidl <= 0.0)
                || (T1 <= 0.0)
                || (self.size_params.cgidl <= 0.0)
                || (vbd > 0.0))
            {
                Igidl = 0.0;
                Ggidld = 0.0;
                Ggidlg = 0.0;
                Ggidlb = 0.0;
            } else {
                dT1_dVd = 1.0 / T0;
                dT1_dVg = -dvgs_eff_dvg * dT1_dVd;
                T2 = self.size_params.bgidl / T1;
                if T2 < 100.0 {
                    Igidl = self.size_params.agidl * self.size_params.weffCJ * T1 * exp(-T2);
                    T3 = Igidl * (1.0 + T2) / T1;
                    Ggidld = T3 * dT1_dVd;
                    Ggidlg = T3 * dT1_dVg;
                } else {
                    Igidl = self.size_params.agidl * self.size_params.weffCJ * 3.720075976e-44;
                    Ggidld = Igidl * dT1_dVd;
                    Ggidlg = Igidl * dT1_dVg;
                    Igidl *= T1;
                }

                T4 = vbd * vbd;
                T5 = -vbd * T4;
                T6 = self.size_params.cgidl + T5;
                T7 = T5 / T6;
                T8 = 3.0 * self.size_params.cgidl * T4 / T6 / T6;
                Ggidld = Ggidld * T7 + Igidl * T8;
                Ggidlg = Ggidlg * T7;
                Ggidlb = -Igidl * T8;
                Igidl *= T7;
            }
            newop.Igidl = Igidl;
            newop.ggidld = Ggidld;
            newop.ggidlg = Ggidlg;
            newop.ggidlb = Ggidlb;
            /* Calculate GISL current  */

            if (self.model.mtrlmod == 0) {
                T1 = (-vds - vgd_eff - self.size_params.egisl) / T0;
            } else {
                T1 = (-vds - vgd_eff - self.size_params.egisl + self.size_params.vfbsd) / T0;
            }

            if ((self.size_params.agisl <= 0.0)
                || (self.size_params.bgisl <= 0.0)
                || (T1 <= 0.0)
                || (self.size_params.cgisl <= 0.0)
                || (vbs > 0.0))
            {
                Igisl = 0.0;
                Ggisls = 0.0;
                Ggislg = 0.0;
                Ggislb = 0.0;
            } else {
                dT1_dVd = 1.0 / T0;
                dT1_dVg = -dvgd_eff_dvg * dT1_dVd;
                T2 = self.size_params.bgisl / T1;
                if T2 < 100.0 {
                    Igisl = self.size_params.agisl * self.size_params.weffCJ * T1 * exp(-T2);
                    T3 = Igisl * (1.0 + T2) / T1;
                    Ggisls = T3 * dT1_dVd;
                    Ggislg = T3 * dT1_dVg;
                } else {
                    Igisl = self.size_params.agisl * self.size_params.weffCJ * 3.720075976e-44;
                    Ggisls = Igisl * dT1_dVd;
                    Ggislg = Igisl * dT1_dVg;
                    Igisl *= T1;
                }

                T4 = vbs * vbs;
                T5 = -vbs * T4;
                T6 = self.size_params.cgisl + T5;
                T7 = T5 / T6;
                T8 = 3.0 * self.size_params.cgisl * T4 / T6 / T6;
                Ggisls = Ggisls * T7 + Igisl * T8;
                Ggislg = Ggislg * T7;
                Ggislb = -Igisl * T8;
                Igisl *= T7;
            }
            newop.Igisl = Igisl;
            newop.ggisls = Ggisls;
            newop.ggislg = Ggislg;
            newop.ggislb = Ggislb;
        } else {
            /* v4.7 New Gidl/GISL model */

            /* GISL */
            if (self.model.mtrlmod == 0) {
                T1 = (-vds - self.size_params.rgisl * vgd_eff - self.size_params.egisl) / T0;
            } else {
                T1 = (-vds - self.size_params.rgisl * vgd_eff - self.size_params.egisl
                    + self.size_params.vfbsd)
                    / T0;
            }

            if ((self.size_params.agisl <= 0.0)
                || (self.size_params.bgisl <= 0.0)
                || (T1 <= 0.0)
                || (self.size_params.cgisl < 0.0))
            {
                Igisl = 0.0;
                Ggisls = 0.0;
                Ggislg = 0.0;
                Ggislb = 0.0;
            } else {
                dT1_dVd = 1.0 / T0;
                dT1_dVg = -self.size_params.rgisl * dT1_dVd * dvgd_eff_dvg;
                T2 = self.size_params.bgisl / T1;
                if T2 < EXPL_THRESHOLD {
                    Igisl = self.size_params.weffCJ * self.size_params.agisl * T1 * exp(-T2);
                    T3 = Igisl / T1 * (T2 + 1.0);
                    Ggisls = T3 * dT1_dVd;
                    Ggislg = T3 * dT1_dVg;
                } else {
                    T3 = self.size_params.weffCJ * self.size_params.agisl * MIN_EXPL;
                    Igisl = T3 * T1;
                    Ggisls = T3 * dT1_dVd;
                    Ggislg = T3 * dT1_dVg;
                }
                T4 = vbs - self.size_params.fgisl;

                if (T4 == 0.0) {
                    T5 = EXPL_THRESHOLD;
                } else {
                    T5 = self.size_params.kgisl / T4;
                }
                if T5 < EXPL_THRESHOLD {
                    T6 = exp(T5);
                    Ggislb = -Igisl * T6 * T5 / T4;
                } else {
                    T6 = MAX_EXPL;
                    Ggislb = 0.0;
                }
                Ggisls *= T6;
                Ggislg *= T6;
                Igisl *= T6;
            }
            newop.Igisl = Igisl;
            newop.ggisls = Ggisls;
            newop.ggislg = Ggislg;
            newop.ggislb = Ggislb;
            /* End of GISL */

            /* GIDL */
            if (self.model.mtrlmod == 0) {
                T1 = (vds - self.size_params.rgidl * vgs_eff - self.size_params.egidl) / T0;
            } else {
                T1 = (vds - self.size_params.rgidl * vgs_eff - self.size_params.egidl
                    + self.size_params.vfbsd)
                    / T0;
            }

            if ((self.size_params.agidl <= 0.0)
                || (self.size_params.bgidl <= 0.0)
                || (T1 <= 0.0)
                || (self.size_params.cgidl < 0.0))
            {
                Igidl = 0.0;
                Ggidld = 0.0;
                Ggidlg = 0.0;
                Ggidlb = 0.0;
            } else {
                dT1_dVd = 1.0 / T0;
                dT1_dVg = -self.size_params.rgidl * dT1_dVd * dvgs_eff_dvg;
                T2 = self.size_params.bgidl / T1;
                if T2 < EXPL_THRESHOLD {
                    Igidl = self.size_params.weffCJ * self.size_params.agidl * T1 * exp(-T2);
                    T3 = Igidl / T1 * (T2 + 1.0);
                    Ggidld = T3 * dT1_dVd;
                    Ggidlg = T3 * dT1_dVg;
                } else {
                    T3 = self.size_params.weffCJ * self.size_params.agidl * MIN_EXPL;
                    Igidl = T3 * T1;
                    Ggidld = T3 * dT1_dVd;
                    Ggidlg = T3 * dT1_dVg;
                }
                T4 = vbd - self.size_params.fgidl;
                if (T4 == 0.0) {
                    T5 = EXPL_THRESHOLD;
                } else {
                    T5 = self.size_params.kgidl / T4;
                }
                if T5 < EXPL_THRESHOLD {
                    T6 = exp(T5);
                    Ggidlb = -Igidl * T6 * T5 / T4;
                } else {
                    T6 = MAX_EXPL;
                    Ggidlb = 0.0;
                }
                Ggidld *= T6;
                Ggidlg *= T6;
                Igidl *= T6;
            }
            newop.Igidl = Igidl;
            newop.ggidld = Ggidld;
            newop.ggidlg = Ggidlg;
            newop.ggidlb = Ggidlb;
            /* End of New GIDL */
        }
        /*End of Gidl*/

        /* Calculate gate tunneling current */
        if (self.model.igcmod != 0) || (self.model.igbmod != 0) {
            Vfb = self.intp.vfbzb;
            V3 = Vfb - Vgs_eff + Vbseff - DELTA_3;
            if (Vfb <= 0.0) {
                T0 = sqrt(V3 * V3 - 4.0 * DELTA_3 * Vfb);
            } else {
                T0 = sqrt(V3 * V3 + 4.0 * DELTA_3 * Vfb);
            }
            T1 = 0.5 * (1.0 + V3 / T0);
            Vfbeff = Vfb - 0.5 * (V3 + T0);
            dVfbeff_dVg = T1 * dVgs_eff_dVg;
            dVfbeff_dVb = -T1;

            Voxacc = Vfb - Vfbeff;
            dVoxacc_dVg = -dVfbeff_dVg;
            dVoxacc_dVb = -dVfbeff_dVb;
            if (Voxacc < 0.0)
            /* WDLiu: Avoiding numerical instability. */
            {
                Voxacc = 0.0;
                dVoxacc_dVg = 0.0;
                dVoxacc_dVb = 0.0;
            }

            T0 = 0.5 * self.size_params.k1ox;
            T3 = Vgs_eff - Vfbeff - Vbseff - Vgsteff;
            if (self.size_params.k1ox == 0.0) {
                Voxdepinv = 0.0;
                dVoxdepinv_dVg = 0.0;
                dVoxdepinv_dVd = 0.0;
                dVoxdepinv_dVb = 0.0;
            } else if T3 < 0.0 {
                Voxdepinv = -T3;
                dVoxdepinv_dVg = -dVgs_eff_dVg + dVfbeff_dVg + dVgsteff_dVg;
                dVoxdepinv_dVd = dVgsteff_dVd;
                dVoxdepinv_dVb = dVfbeff_dVb + 1.0 + dVgsteff_dVb;
            } else {
                T1 = sqrt(T0 * T0 + T3);
                T2 = T0 / T1;
                Voxdepinv = self.size_params.k1ox * (T1 - T0);
                dVoxdepinv_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg);
                dVoxdepinv_dVd = -T2 * dVgsteff_dVd;
                dVoxdepinv_dVb = -T2 * (dVfbeff_dVb + 1.0 + dVgsteff_dVb);
            }

            Voxdepinv += Vgsteff;
            dVoxdepinv_dVg += dVgsteff_dVg;
            dVoxdepinv_dVd += dVgsteff_dVd;
            dVoxdepinv_dVb += dVgsteff_dVb;
        }

        if (self.model.tempmod < 2) {
            tmp = Vtm;
        } else
        /* self.model.tempmod = 2, 3*/
        {
            tmp = Vtm0;
        }
        if self.model.igcmod != 0 {
            T0 = tmp * self.size_params.nigc;
            // FIXME: enum-ize
            if self.model.igcmod == 1 {
                VxNVt = (Vgs_eff - self.model.p() * self.intp.vth0) / T0;
                if VxNVt > EXP_THRESHOLD {
                    Vaux = Vgs_eff - self.model.p() * self.intp.vth0;
                    dVaux_dVg = dVgs_eff_dVg;
                    dVaux_dVd = 0.0;
                    dVaux_dVb = 0.0;
                }
            } else {
                // if (self.model.igcmod == 2) {
                VxNVt = (Vgs_eff - newop.von) / T0;
                if VxNVt > EXP_THRESHOLD {
                    Vaux = Vgs_eff - newop.von;
                    dVaux_dVg = dVgs_eff_dVg;
                    dVaux_dVd = -dVth_dVd;
                    dVaux_dVb = -dVth_dVb;
                }
            }
            if VxNVt < -EXP_THRESHOLD {
                Vaux = T0 * log(1.0 + MIN_EXP);
                dVaux_dVg = 0.0;
                dVaux_dVd = 0.0;
                dVaux_dVb = 0.0;
            } else if (VxNVt >= -EXP_THRESHOLD) && (VxNVt <= EXP_THRESHOLD) {
                ExpVxNVt = exp(VxNVt);
                Vaux = T0 * log(1.0 + ExpVxNVt);
                dVaux_dVg = ExpVxNVt / (1.0 + ExpVxNVt);
                // FIXME: enum-ize
                if self.model.igcmod == 1 {
                    dVaux_dVd = 0.0;
                    dVaux_dVb = 0.0;
                } else {
                    // if (self.model.igcmod == 2) {
                    dVaux_dVd = -dVaux_dVg * dVth_dVd; /* Synopsys 08/30/2013 modify */
                    dVaux_dVb = -dVaux_dVg * dVth_dVb; /* Synopsys 08/30/2013 modify */
                }
                dVaux_dVg *= dVgs_eff_dVg;
            }

            T2 = Vgs_eff * Vaux;
            dT2_dVg = dVgs_eff_dVg * Vaux + Vgs_eff * dVaux_dVg;
            dT2_dVd = Vgs_eff * dVaux_dVd;
            dT2_dVb = Vgs_eff * dVaux_dVb;

            T11 = self.size_params.Aechvb;
            T12 = self.size_params.Bechvb;
            T3 = self.size_params.aigc * self.size_params.cigc - self.size_params.bigc;
            T4 = self.size_params.bigc * self.size_params.cigc;
            T5 = T12 * (self.size_params.aigc + T3 * Voxdepinv - T4 * Voxdepinv * Voxdepinv);

            if T5 > EXP_THRESHOLD {
                T6 = MAX_EXP;
                dT6_dVg = 0.0;
                dT6_dVd = 0.0;
                dT6_dVb = 0.0;
            } else if T5 < -EXP_THRESHOLD {
                T6 = MIN_EXP;
                dT6_dVg = 0.0;
                dT6_dVd = 0.0;
                dT6_dVb = 0.0;
            } else {
                T6 = exp(T5);
                dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxdepinv);
                dT6_dVd = dT6_dVg * dVoxdepinv_dVd;
                dT6_dVb = dT6_dVg * dVoxdepinv_dVb;
                dT6_dVg *= dVoxdepinv_dVg;
            }

            Igc = T11 * T2 * T6;
            dIgc_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
            dIgc_dVd = T11 * (T2 * dT6_dVd + T6 * dT2_dVd);
            dIgc_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);

            if self.model.pigcd == 0.0 {
                // FIXME: reference implementation condition is "if pigcd given",
                // maybe make this an Option.
                // Here we are using the 0.0 default value.
                Pigcd = self.size_params.pigcd;
                dPigcd_dVg = 0.0;
                dPigcd_dVd = 0.0;
                dPigcd_dVb = 0.0;
            } else {
                T11 = -self.size_params.Bechvb;
                T12 = Vgsteff + 1.0e-20;
                T13 = T11 / T12 / T12;
                T14 = -T13 / T12;
                Pigcd = T13 * (1.0 - 0.5 * Vdseff / T12);
                dPigcd_dVg = T14 * (2.0 + 0.5 * (dVdseff_dVg - 3.0 * Vdseff / T12));
                dPigcd_dVd = 0.5 * T14 * dVdseff_dVd;
                dPigcd_dVb = 0.5 * T14 * dVdseff_dVb;
            }

            T7 = -Pigcd * Vdseff;
            dT7_dVg = -Vdseff * dPigcd_dVg - Pigcd * dVdseff_dVg;
            dT7_dVd = -Vdseff * dPigcd_dVd - Pigcd * dVdseff_dVd + dT7_dVg * dVgsteff_dVd;
            dT7_dVb = -Vdseff * dPigcd_dVb - Pigcd * dVdseff_dVb + dT7_dVg * dVgsteff_dVb;
            dT7_dVg *= dVgsteff_dVg;
            T8 = T7 * T7 + 2.0e-4;
            dT8_dVg = 2.0 * T7;
            dT8_dVd = dT8_dVg * dT7_dVd;
            dT8_dVb = dT8_dVg * dT7_dVb;
            dT8_dVg *= dT7_dVg;

            if T7 > EXP_THRESHOLD {
                T9 = MAX_EXP;
                dT9_dVg = 0.0;
                dT9_dVd = 0.0;
                dT9_dVb = 0.0;
            } else if T7 < -EXP_THRESHOLD {
                T9 = MIN_EXP;
                dT9_dVg = 0.0;
                dT9_dVd = 0.0;
                dT9_dVb = 0.0;
            } else {
                T9 = exp(T7);
                dT9_dVg = T9 * dT7_dVg;
                dT9_dVd = T9 * dT7_dVd;
                dT9_dVb = T9 * dT7_dVb;
            }

            T0 = T8 * T8;
            T1 = T9 - 1.0 + 1.0e-4;
            T10 = (T1 - T7) / T8;
            dT10_dVg = (dT9_dVg - dT7_dVg - T10 * dT8_dVg) / T8;
            dT10_dVd = (dT9_dVd - dT7_dVd - T10 * dT8_dVd) / T8;
            dT10_dVb = (dT9_dVb - dT7_dVb - T10 * dT8_dVb) / T8;

            Igcs = Igc * T10;
            dIgcs_dVg = dIgc_dVg * T10 + Igc * dT10_dVg;
            dIgcs_dVd = dIgc_dVd * T10 + Igc * dT10_dVd;
            dIgcs_dVb = dIgc_dVb * T10 + Igc * dT10_dVb;

            T1 = T9 - 1.0 - 1.0e-4;
            T10 = (T7 * T9 - T1) / T8;
            dT10_dVg = (dT7_dVg * T9 + (T7 - 1.0) * dT9_dVg - T10 * dT8_dVg) / T8;
            dT10_dVd = (dT7_dVd * T9 + (T7 - 1.0) * dT9_dVd - T10 * dT8_dVd) / T8;
            dT10_dVb = (dT7_dVb * T9 + (T7 - 1.0) * dT9_dVb - T10 * dT8_dVb) / T8;
            Igcd = Igc * T10;
            dIgcd_dVg = dIgc_dVg * T10 + Igc * dT10_dVg;
            dIgcd_dVd = dIgc_dVd * T10 + Igc * dT10_dVd;
            dIgcd_dVb = dIgc_dVb * T10 + Igc * dT10_dVb;

            newop.Igcs = Igcs;
            newop.gIgcsg = dIgcs_dVg;
            newop.gIgcsd = dIgcs_dVd;
            newop.gIgcsb = dIgcs_dVb * dVbseff_dVb;
            newop.Igcd = Igcd;
            newop.gIgcdg = dIgcd_dVg;
            newop.gIgcdd = dIgcd_dVd;
            newop.gIgcdb = dIgcd_dVb * dVbseff_dVb;

            T0 = vgs - (self.size_params.vfbsd + self.size_params.vfbsdoff);
            vgs_eff = sqrt(T0 * T0 + 1.0e-4);
            dvgs_eff_dvg = T0 / vgs_eff;

            T2 = vgs * vgs_eff;
            dT2_dVg = vgs * dvgs_eff_dvg + vgs_eff;
            T11 = self.size_params.AechvbEdgeS;
            T12 = self.size_params.BechvbEdge;
            T3 = self.size_params.aigs * self.size_params.cigs - self.size_params.bigs;
            T4 = self.size_params.bigs * self.size_params.cigs;
            T5 = T12 * (self.size_params.aigs + T3 * vgs_eff - T4 * vgs_eff * vgs_eff);
            if T5 > EXP_THRESHOLD {
                T6 = MAX_EXP;
                dT6_dVg = 0.0;
            } else if T5 < -EXP_THRESHOLD {
                T6 = MIN_EXP;
                dT6_dVg = 0.0;
            } else {
                T6 = exp(T5);
                dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * vgs_eff) * dvgs_eff_dvg;
            }
            Igs = T11 * T2 * T6;
            dIgs_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
            dIgs_dVs = -dIgs_dVg;

            T0 = vgd - (self.size_params.vfbsd + self.size_params.vfbsdoff);
            vgd_eff = sqrt(T0 * T0 + 1.0e-4);
            dvgd_eff_dvg = T0 / vgd_eff;

            T2 = vgd * vgd_eff;
            dT2_dVg = vgd * dvgd_eff_dvg + vgd_eff;
            T11 = self.size_params.AechvbEdgeD;
            T3 = self.size_params.aigd * self.size_params.cigd - self.size_params.bigd;
            T4 = self.size_params.bigd * self.size_params.cigd;
            T5 = T12 * (self.size_params.aigd + T3 * vgd_eff - T4 * vgd_eff * vgd_eff);
            if T5 > EXP_THRESHOLD {
                T6 = MAX_EXP;
                dT6_dVg = 0.0;
            } else if T5 < -EXP_THRESHOLD {
                T6 = MIN_EXP;
                dT6_dVg = 0.0;
            } else {
                T6 = exp(T5);
                dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * vgd_eff) * dvgd_eff_dvg;
            }
            Igd = T11 * T2 * T6;
            dIgd_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
            dIgd_dVd = -dIgd_dVg;

            newop.Igs = Igs;
            newop.gIgsg = dIgs_dVg;
            newop.gIgss = dIgs_dVs;
            newop.Igd = Igd;
            newop.gIgdg = dIgd_dVg;
            newop.gIgdd = dIgd_dVd;
        } else {
            newop.Igcs = 0.0;
            newop.gIgcsg = 0.0;
            newop.gIgcsd = 0.0;
            newop.gIgcsb = 0.0;
            newop.Igcd = 0.0;
            newop.gIgcdg = 0.0;
            newop.gIgcdd = 0.0;
            newop.gIgcdb = 0.0;
            newop.Igs = 0.0;
            newop.gIgsg = 0.0;
            newop.gIgss = 0.0;
            newop.Igd = 0.0;
            newop.gIgdg = 0.0;
            newop.gIgdd = 0.0;
        }

        if self.model.igbmod != 0 {
            T0 = tmp * self.size_params.nigbacc;
            T1 = -Vgs_eff + Vbseff + Vfb;
            VxNVt = T1 / T0;
            if VxNVt > EXP_THRESHOLD {
                Vaux = T1;
                dVaux_dVg = -dVgs_eff_dVg;
                dVaux_dVb = 1.0;
            } else if VxNVt < -EXP_THRESHOLD {
                Vaux = T0 * log(1.0 + MIN_EXP);
                dVaux_dVg = 0.0;
                dVaux_dVb = 0.0;
            } else {
                ExpVxNVt = exp(VxNVt);
                Vaux = T0 * log(1.0 + ExpVxNVt);
                dVaux_dVb = ExpVxNVt / (1.0 + ExpVxNVt);
                dVaux_dVg = -dVaux_dVb * dVgs_eff_dVg;
            }

            T2 = (Vgs_eff - Vbseff) * Vaux;
            dT2_dVg = dVgs_eff_dVg * Vaux + (Vgs_eff - Vbseff) * dVaux_dVg;
            dT2_dVb = -Vaux + (Vgs_eff - Vbseff) * dVaux_dVb;

            T11 = 4.97232e-7
                * self.size_params.weff
                * self.size_params.leff
                * self.size_params.ToxRatio;
            T12 = -7.45669e11 * toxe;
            T3 = self.size_params.aigbacc * self.size_params.cigbacc - self.size_params.bigbacc;
            T4 = self.size_params.bigbacc * self.size_params.cigbacc;
            T5 = T12 * (self.size_params.aigbacc + T3 * Voxacc - T4 * Voxacc * Voxacc);

            if T5 > EXP_THRESHOLD {
                T6 = MAX_EXP;
                dT6_dVg = 0.0;
                dT6_dVb = 0.0;
            } else if T5 < -EXP_THRESHOLD {
                T6 = MIN_EXP;
                dT6_dVg = 0.0;
                dT6_dVb = 0.0;
            } else {
                T6 = exp(T5);
                dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxacc);
                dT6_dVb = dT6_dVg * dVoxacc_dVb;
                dT6_dVg *= dVoxacc_dVg;
            }

            Igbacc = T11 * T2 * T6;
            dIgbacc_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
            dIgbacc_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);

            T0 = tmp * self.size_params.nigbinv;
            T1 = Voxdepinv - self.size_params.eigbinv;
            VxNVt = T1 / T0;
            if VxNVt > EXP_THRESHOLD {
                Vaux = T1;
                dVaux_dVg = dVoxdepinv_dVg;
                dVaux_dVd = dVoxdepinv_dVd;
                dVaux_dVb = dVoxdepinv_dVb;
            } else if VxNVt < -EXP_THRESHOLD {
                Vaux = T0 * log(1.0 + MIN_EXP);
                dVaux_dVg = 0.0;
                dVaux_dVd = 0.0;
                dVaux_dVb = 0.0;
            } else {
                ExpVxNVt = exp(VxNVt);
                Vaux = T0 * log(1.0 + ExpVxNVt);
                dVaux_dVg = ExpVxNVt / (1.0 + ExpVxNVt);
                dVaux_dVd = dVaux_dVg * dVoxdepinv_dVd;
                dVaux_dVb = dVaux_dVg * dVoxdepinv_dVb;
                dVaux_dVg *= dVoxdepinv_dVg;
            }

            T2 = (Vgs_eff - Vbseff) * Vaux;
            dT2_dVg = dVgs_eff_dVg * Vaux + (Vgs_eff - Vbseff) * dVaux_dVg;
            dT2_dVd = (Vgs_eff - Vbseff) * dVaux_dVd;
            dT2_dVb = -Vaux + (Vgs_eff - Vbseff) * dVaux_dVb;

            T11 *= 0.75610;
            T12 *= 1.31724;
            T3 = self.size_params.aigbinv * self.size_params.cigbinv - self.size_params.bigbinv;
            T4 = self.size_params.bigbinv * self.size_params.cigbinv;
            T5 = T12 * (self.size_params.aigbinv + T3 * Voxdepinv - T4 * Voxdepinv * Voxdepinv);

            if T5 > EXP_THRESHOLD {
                T6 = MAX_EXP;
                dT6_dVg = 0.0;
                dT6_dVd = 0.0;
                dT6_dVb = 0.0;
            } else if T5 < -EXP_THRESHOLD {
                T6 = MIN_EXP;
                dT6_dVg = 0.0;
                dT6_dVd = 0.0;
                dT6_dVb = 0.0;
            } else {
                T6 = exp(T5);
                dT6_dVg = T6 * T12 * (T3 - 2.0 * T4 * Voxdepinv);
                dT6_dVd = dT6_dVg * dVoxdepinv_dVd;
                dT6_dVb = dT6_dVg * dVoxdepinv_dVb;
                dT6_dVg *= dVoxdepinv_dVg;
            }

            Igbinv = T11 * T2 * T6;
            dIgbinv_dVg = T11 * (T2 * dT6_dVg + T6 * dT2_dVg);
            dIgbinv_dVd = T11 * (T2 * dT6_dVd + T6 * dT2_dVd);
            dIgbinv_dVb = T11 * (T2 * dT6_dVb + T6 * dT2_dVb);

            newop.Igb = Igbinv + Igbacc;
            newop.gIgbg = dIgbinv_dVg + dIgbacc_dVg;
            newop.gIgbd = dIgbinv_dVd;
            newop.gIgbb = (dIgbinv_dVb + dIgbacc_dVb) * dVbseff_dVb;
        } else {
            newop.Igb = 0.0;
            newop.gIgbg = 0.0;
            newop.gIgbd = 0.0;
            newop.gIgbs = 0.0;
            newop.gIgbb = 0.0;
        } /* End of Gate current */

        if self.intp.nf != 1.0 {
            cdrain *= self.intp.nf;
            newop.gds *= self.intp.nf;
            newop.gm *= self.intp.nf;
            newop.gmbs *= self.intp.nf;
            newop.IdovVds *= self.intp.nf;

            newop.gbbs *= self.intp.nf;
            newop.gbgs *= self.intp.nf;
            newop.gbds *= self.intp.nf;
            newop.csub *= self.intp.nf;

            newop.Igidl *= self.intp.nf;
            newop.ggidld *= self.intp.nf;
            newop.ggidlg *= self.intp.nf;
            newop.ggidlb *= self.intp.nf;

            newop.Igisl *= self.intp.nf;
            newop.ggisls *= self.intp.nf;
            newop.ggislg *= self.intp.nf;
            newop.ggislb *= self.intp.nf;

            newop.Igcs *= self.intp.nf;
            newop.gIgcsg *= self.intp.nf;
            newop.gIgcsd *= self.intp.nf;
            newop.gIgcsb *= self.intp.nf;
            newop.Igcd *= self.intp.nf;
            newop.gIgcdg *= self.intp.nf;
            newop.gIgcdd *= self.intp.nf;
            newop.gIgcdb *= self.intp.nf;

            newop.Igs *= self.intp.nf;
            newop.gIgsg *= self.intp.nf;
            newop.gIgss *= self.intp.nf;
            newop.Igd *= self.intp.nf;
            newop.gIgdg *= self.intp.nf;
            newop.gIgdd *= self.intp.nf;

            newop.Igb *= self.intp.nf;
            newop.gIgbg *= self.intp.nf;
            newop.gIgbd *= self.intp.nf;
            newop.gIgbb *= self.intp.nf;
        }

        newop.ggidls = -(newop.ggidld + newop.ggidlg + newop.ggidlb);
        newop.ggisld = -(newop.ggisls + newop.ggislg + newop.ggislb);
        newop.gIgbs = -(newop.gIgbg + newop.gIgbd + newop.gIgbb);
        newop.gIgcss = -(newop.gIgcsg + newop.gIgcsd + newop.gIgcsb);
        newop.gIgcds = -(newop.gIgcdg + newop.gIgcdd + newop.gIgcdb);
        newop.cd = cdrain;

        /* Calculations for noise analysis */

        if self.model.tnoimod == 0 {
            Abulk = Abulk0 * self.size_params.abulkCVfactor;
            Vdsat = Vgsteff / Abulk;
            T0 = Vdsat - Vds - DELTA_4;
            T1 = sqrt(T0 * T0 + 4.0 * DELTA_4 * Vdsat);
            if T0 >= 0.0 {
                Vdseff = Vdsat - 0.5 * (T0 + T1);
            } else {
                T3 = (DELTA_4 + DELTA_4) / (T1 - T0);
                T4 = 1.0 - T3;
                T5 = Vdsat * T3 / (T1 - T0);
                Vdseff = Vdsat * T4;
            }
            if Vds == 0.0 {
                Vdseff = 0.0;
            }

            T0 = Abulk * Vdseff;
            T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.0e-20);
            T2 = Vdseff / T1;
            T3 = T0 * T2;
            newop.qinv = Coxeff
                * self.size_params.weffCV
                * self.intp.nf
                * self.size_params.leffCV
                * (Vgsteff - 0.5 * T0 + Abulk * T3);
        } else if self.model.tnoimod == 2 {
            newop.noiGd0 = self.intp.nf * beta * Vgsteff / (1.0 + gche * Rds);
        }

        /*
         *  BSIM4 C-V begins
         */

        if (self.model.xpart < 0.0) || (!ChargeComputationNeeded) {
            qgate = 0.0;
            qdrn = 0.0;
            qsrc = 0.0;
            qbulk = 0.0;
            newop.cggb = 0.0;
            newop.cgsb = 0.0;
            newop.cgdb = 0.0;
            newop.cdgb = 0.0;
            newop.cdsb = 0.0;
            newop.cddb = 0.0;
            newop.cbgb = 0.0;
            newop.cbsb = 0.0;
            newop.cbdb = 0.0;
            newop.csgb = 0.0;
            newop.cssb = 0.0;
            newop.csdb = 0.0;
            newop.cgbb = 0.0;
            newop.csbb = 0.0;
            newop.cdbb = 0.0;
            newop.cbbb = 0.0;
            newop.cqdb = 0.0;
            newop.cqsb = 0.0;
            newop.cqgb = 0.0;
            newop.cqbb = 0.0;
            newop.gtau = 0.0;
        } else if let AnalysisInfo::TRAN(_, state) = an {
            //FIXME: always true, or we should generate an error
            if self.model.capmod == 0 {
                if Vbseff < 0.0 {
                    VbseffCV = Vbs; /*4.6.2*/
                    dVbseffCV_dVb = 1.0;
                } else {
                    VbseffCV = self.size_params.phi - Phis;
                    dVbseffCV_dVb = -dPhis_dVb * dVbseff_dVb; /*4.6.2*/
                }

                Vfb = self.size_params.vfbcv;
                Vth = Vfb + self.size_params.phi + self.size_params.k1ox * sqrtPhis;
                Vgst = Vgs_eff - Vth;
                dVth_dVb = self.size_params.k1ox * dsqrtPhis_dVb * dVbseff_dVb; /*4.6.2*/
                dVgst_dVb = -dVth_dVb;
                dVgst_dVg = dVgs_eff_dVg;

                CoxWL = self.intp.coxe
                    * self.size_params.weffCV
                    * self.size_params.leffCV
                    * self.intp.nf;
                Arg1 = Vgs_eff - VbseffCV - Vfb;

                if Arg1 <= 0.0 {
                    qgate = CoxWL * Arg1;
                    qbulk = -qgate;
                    qdrn = 0.0;

                    newop.cggb = CoxWL * dVgs_eff_dVg;
                    newop.cgdb = 0.0;
                    newop.cgsb = CoxWL * (dVbseffCV_dVb - dVgs_eff_dVg);

                    newop.cdgb = 0.0;
                    newop.cddb = 0.0;
                    newop.cdsb = 0.0;

                    newop.cbgb = -CoxWL * dVgs_eff_dVg;
                    newop.cbdb = 0.0;
                    newop.cbsb = -newop.cgsb;
                }
                /* Arg1 <= 0.0, end of accumulation */
                else if Vgst <= 0.0 {
                    T1 = 0.5 * self.size_params.k1ox;
                    T2 = sqrt(T1 * T1 + Arg1);
                    qgate = CoxWL * self.size_params.k1ox * (T2 - T1);
                    qbulk = -qgate;
                    qdrn = 0.0;

                    T0 = CoxWL * T1 / T2;
                    newop.cggb = T0 * dVgs_eff_dVg;
                    newop.cgdb = 0.0;
                    newop.cgsb = T0 * (dVbseffCV_dVb - dVgs_eff_dVg);

                    newop.cdgb = 0.0;
                    newop.cddb = 0.0;
                    newop.cdsb = 0.0;

                    newop.cbgb = -newop.cggb;
                    newop.cbdb = 0.0;
                    newop.cbsb = -newop.cgsb;
                }
                /* Vgst <= 0.0, end of depletion */
                else {
                    One_Third_CoxWL = CoxWL / 3.0;
                    Two_Third_CoxWL = 2.0 * One_Third_CoxWL;

                    AbulkCV = Abulk0 * self.size_params.abulkCVfactor;
                    dAbulkCV_dVb = self.size_params.abulkCVfactor * dAbulk0_dVb * dVbseff_dVb;

                    dVdsat_dVg = 1.0 / AbulkCV; /*4.6.2*/
                    Vdsat = Vgst * dVdsat_dVg;
                    dVdsat_dVb = -(Vdsat * dAbulkCV_dVb + dVth_dVb) * dVdsat_dVg;

                    if self.model.xpart > 0.5 {
                        /* 0/100 Charge partition model */
                        if Vdsat <= Vds {
                            /* saturation region */
                            T1 = Vdsat / 3.0;
                            qgate = CoxWL * (Vgs_eff - Vfb - self.size_params.phi - T1);
                            T2 = -Two_Third_CoxWL * Vgst;
                            qbulk = -(qgate + T2);
                            qdrn = 0.0;

                            newop.cggb = One_Third_CoxWL * (3.0 - dVdsat_dVg) * dVgs_eff_dVg;
                            T2 = -One_Third_CoxWL * dVdsat_dVb;
                            newop.cgsb = -(newop.cggb + T2);
                            newop.cgdb = 0.0;

                            newop.cdgb = 0.0;
                            newop.cddb = 0.0;
                            newop.cdsb = 0.0;

                            newop.cbgb = -(newop.cggb - Two_Third_CoxWL * dVgs_eff_dVg);
                            T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                            newop.cbsb = -(newop.cbgb + T3);
                            newop.cbdb = 0.0;
                        } else {
                            /* linear region */
                            Alphaz = Vgst / Vdsat;
                            T1 = 2.0 * Vdsat - Vds;
                            T2 = Vds / (3.0 * T1);
                            T3 = T2 * Vds;
                            T9 = 0.25 * CoxWL;
                            T4 = T9 * Alphaz;
                            T7 = 2.0 * Vds - T1 - 3.0 * T3;
                            T8 = T3 - T1 - 2.0 * Vds;
                            qgate =
                                CoxWL * (Vgs_eff - Vfb - self.size_params.phi - 0.5 * (Vds - T3));
                            T10 = T4 * T8;
                            qdrn = T4 * T7;
                            qbulk = -(qgate + qdrn + T10);

                            T5 = T3 / T1;
                            newop.cggb = CoxWL * (1.0 - T5 * dVdsat_dVg) * dVgs_eff_dVg;
                            T11 = -CoxWL * T5 * dVdsat_dVb;
                            newop.cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                            newop.cgsb = -(newop.cggb + T11 + newop.cgdb);
                            T6 = 1.0 / Vdsat;
                            dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                            dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);
                            T7 = T9 * T7;
                            T8 = T9 * T8;
                            T9 = 2.0 * T4 * (1.0 - 3.0 * T5);
                            newop.cdgb = (T7 * dAlphaz_dVg - T9 * dVdsat_dVg) * dVgs_eff_dVg;
                            T12 = T7 * dAlphaz_dVb - T9 * dVdsat_dVb;
                            newop.cddb = T4 * (3.0 - 6.0 * T2 - 3.0 * T5);
                            newop.cdsb = -(newop.cdgb + T12 + newop.cddb);

                            T9 = 2.0 * T4 * (1.0 + T5);
                            T10 = (T8 * dAlphaz_dVg - T9 * dVdsat_dVg) * dVgs_eff_dVg;
                            T11 = T8 * dAlphaz_dVb - T9 * dVdsat_dVb;
                            T12 = T4 * (2.0 * T2 + T5 - 1.0);
                            T0 = -(T10 + T11 + T12);

                            newop.cbgb = -(newop.cggb + newop.cdgb + T10);
                            newop.cbdb = -(newop.cgdb + newop.cddb + T12);
                            newop.cbsb = -(newop.cgsb + newop.cdsb + T0);
                        }
                    } else if self.model.xpart < 0.5 {
                        /* 40/60 Charge partition model */
                        if Vds >= Vdsat {
                            /* saturation region */
                            T1 = Vdsat / 3.0;
                            qgate = CoxWL * (Vgs_eff - Vfb - self.size_params.phi - T1);
                            T2 = -Two_Third_CoxWL * Vgst;
                            qbulk = -(qgate + T2);
                            qdrn = 0.4 * T2;

                            newop.cggb = One_Third_CoxWL * (3.0 - dVdsat_dVg) * dVgs_eff_dVg;
                            T2 = -One_Third_CoxWL * dVdsat_dVb;
                            newop.cgsb = -(newop.cggb + T2);
                            newop.cgdb = 0.0;

                            T3 = 0.4 * Two_Third_CoxWL;
                            newop.cdgb = -T3 * dVgs_eff_dVg;
                            newop.cddb = 0.0;
                            T4 = T3 * dVth_dVb;
                            newop.cdsb = -(T4 + newop.cdgb);

                            newop.cbgb = -(newop.cggb - Two_Third_CoxWL * dVgs_eff_dVg);
                            T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                            newop.cbsb = -(newop.cbgb + T3);
                            newop.cbdb = 0.0;
                        } else {
                            /* linear region  */
                            Alphaz = Vgst / Vdsat;
                            T1 = 2.0 * Vdsat - Vds;
                            T2 = Vds / (3.0 * T1);
                            T3 = T2 * Vds;
                            T9 = 0.25 * CoxWL;
                            T4 = T9 * Alphaz;
                            qgate =
                                CoxWL * (Vgs_eff - Vfb - self.size_params.phi - 0.5 * (Vds - T3));

                            T5 = T3 / T1;
                            newop.cggb = CoxWL * (1.0 - T5 * dVdsat_dVg) * dVgs_eff_dVg;
                            tmp = -CoxWL * T5 * dVdsat_dVb;
                            newop.cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                            newop.cgsb = -(newop.cggb + newop.cgdb + tmp);

                            T6 = 1.0 / Vdsat;
                            dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                            dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

                            T6 = 8.0 * Vdsat * Vdsat - 6.0 * Vdsat * Vds + 1.2 * Vds * Vds;
                            T8 = T2 / T1;
                            T7 = Vds - T1 - T8 * T6;
                            qdrn = T4 * T7;
                            T7 *= T9;
                            tmp = T8 / T1;
                            tmp1 = T4 * (2.0 - 4.0 * tmp * T6 + T8 * (16.0 * Vdsat - 6.0 * Vds));

                            newop.cdgb = (T7 * dAlphaz_dVg - tmp1 * dVdsat_dVg) * dVgs_eff_dVg;
                            T10 = T7 * dAlphaz_dVb - tmp1 * dVdsat_dVb;
                            newop.cddb = T4
                                * (2.0 - (1.0 / (3.0 * T1 * T1) + 2.0 * tmp) * T6
                                    + T8 * (6.0 * Vdsat - 2.4 * Vds));
                            newop.cdsb = -(newop.cdgb + T10 + newop.cddb);

                            T7 = 2.0 * (T1 + T3);
                            qbulk = -(qgate - T4 * T7);
                            T7 *= T9;
                            T0 = 4.0 * T4 * (1.0 - T5);
                            T12 = (-T7 * dAlphaz_dVg - T0 * dVdsat_dVg) * dVgs_eff_dVg - newop.cdgb; /*4.6.2*/
                            T11 = -T7 * dAlphaz_dVb - T10 - T0 * dVdsat_dVb;
                            T10 = -4.0 * T4 * (T2 - 0.5 + 0.5 * T5) - newop.cddb;
                            tmp = -(T10 + T11 + T12);

                            newop.cbgb = -(newop.cggb + newop.cdgb + T12);
                            newop.cbdb = -(newop.cgdb + newop.cddb + T10);
                            newop.cbsb = -(newop.cgsb + newop.cdsb + tmp);
                        }
                    } else {
                        /* 50/50 partitioning */
                        if Vds >= Vdsat {
                            /* saturation region */
                            T1 = Vdsat / 3.0;
                            qgate = CoxWL * (Vgs_eff - Vfb - self.size_params.phi - T1);
                            T2 = -Two_Third_CoxWL * Vgst;
                            qbulk = -(qgate + T2);
                            qdrn = 0.5 * T2;

                            newop.cggb = One_Third_CoxWL * (3.0 - dVdsat_dVg) * dVgs_eff_dVg;
                            T2 = -One_Third_CoxWL * dVdsat_dVb;
                            newop.cgsb = -(newop.cggb + T2);
                            newop.cgdb = 0.0;

                            newop.cdgb = -One_Third_CoxWL * dVgs_eff_dVg;
                            newop.cddb = 0.0;
                            T4 = One_Third_CoxWL * dVth_dVb;
                            newop.cdsb = -(T4 + newop.cdgb);

                            newop.cbgb = -(newop.cggb - Two_Third_CoxWL * dVgs_eff_dVg);
                            T3 = -(T2 + Two_Third_CoxWL * dVth_dVb);
                            newop.cbsb = -(newop.cbgb + T3);
                            newop.cbdb = 0.0;
                        } else {
                            /* linear region */
                            Alphaz = Vgst / Vdsat;
                            T1 = 2.0 * Vdsat - Vds;
                            T2 = Vds / (3.0 * T1);
                            T3 = T2 * Vds;
                            T9 = 0.25 * CoxWL;
                            T4 = T9 * Alphaz;
                            qgate =
                                CoxWL * (Vgs_eff - Vfb - self.size_params.phi - 0.5 * (Vds - T3));

                            T5 = T3 / T1;
                            newop.cggb = CoxWL * (1.0 - T5 * dVdsat_dVg) * dVgs_eff_dVg;
                            tmp = -CoxWL * T5 * dVdsat_dVb;
                            newop.cgdb = CoxWL * (T2 - 0.5 + 0.5 * T5);
                            newop.cgsb = -(newop.cggb + newop.cgdb + tmp);

                            T6 = 1.0 / Vdsat;
                            dAlphaz_dVg = T6 * (1.0 - Alphaz * dVdsat_dVg);
                            dAlphaz_dVb = -T6 * (dVth_dVb + Alphaz * dVdsat_dVb);

                            T7 = T1 + T3;
                            qdrn = -T4 * T7;
                            qbulk = -(qgate + qdrn + qdrn);
                            T7 *= T9;
                            T0 = T4 * (2.0 * T5 - 2.0);

                            newop.cdgb = (T0 * dVdsat_dVg - T7 * dAlphaz_dVg) * dVgs_eff_dVg;
                            T12 = T0 * dVdsat_dVb - T7 * dAlphaz_dVb;
                            newop.cddb = T4 * (1.0 - 2.0 * T2 - T5);
                            newop.cdsb = -(newop.cdgb + T12 + newop.cddb);
                            newop.cbgb = -(newop.cggb + 2.0 * newop.cdgb);
                            newop.cbdb = -(newop.cgdb + 2.0 * newop.cddb);
                            newop.cbsb = -(newop.cgsb + 2.0 * newop.cdsb);
                        } /* end of linear region */
                    } /* end of 50/50 partition */
                } /* end of inversion */
            }
            /* end of capmod=0 */
            else {
                if Vbseff < 0.0 {
                    VbseffCV = Vbseff;
                    dVbseffCV_dVb = 1.0;
                } else {
                    VbseffCV = self.size_params.phi - Phis;
                    dVbseffCV_dVb = -dPhis_dVb;
                }

                CoxWL = self.intp.coxe
                    * self.size_params.weffCV
                    * self.size_params.leffCV
                    * self.intp.nf;

                if (self.model.cvchargemod == 0) {
                    /* Seperate VgsteffCV with noff and voffcv */
                    noff = n * self.size_params.noff;
                    dnoff_dVd = self.size_params.noff * dn_dVd;
                    dnoff_dVb = self.size_params.noff * dn_dVb;
                    T0 = Vtm * noff;
                    voffcv = self.size_params.voffcv;
                    VgstNVt = (Vgst - voffcv) / T0;

                    if VgstNVt > EXP_THRESHOLD {
                        Vgsteff = Vgst - voffcv;
                        dVgsteff_dVg = dVgs_eff_dVg;
                        dVgsteff_dVd = -dVth_dVd;
                        dVgsteff_dVb = -dVth_dVb;
                    } else if VgstNVt < -EXP_THRESHOLD {
                        Vgsteff = T0 * log(1.0 + MIN_EXP);
                        dVgsteff_dVg = 0.0;
                        dVgsteff_dVd = Vgsteff / noff;
                        dVgsteff_dVb = dVgsteff_dVd * dnoff_dVb;
                        dVgsteff_dVd *= dnoff_dVd;
                    } else {
                        ExpVgst = exp(VgstNVt);
                        Vgsteff = T0 * log(1.0 + ExpVgst);
                        dVgsteff_dVg = ExpVgst / (1.0 + ExpVgst);
                        dVgsteff_dVd = -dVgsteff_dVg
                            * (dVth_dVd + (Vgst - voffcv) / noff * dnoff_dVd)
                            + Vgsteff / noff * dnoff_dVd;
                        dVgsteff_dVb = -dVgsteff_dVg
                            * (dVth_dVb + (Vgst - voffcv) / noff * dnoff_dVb)
                            + Vgsteff / noff * dnoff_dVb;
                        dVgsteff_dVg *= dVgs_eff_dVg;
                    }
                /* End of VgsteffCV for cvchargeMod = 0 */
                } else {
                    T0 = n * Vtm;
                    T1 = self.size_params.mstarcv * Vgst;
                    T2 = T1 / T0;
                    if T2 > EXP_THRESHOLD {
                        T10 = T1;
                        dT10_dVg = self.size_params.mstarcv * dVgs_eff_dVg;
                        dT10_dVd = -dVth_dVd * self.size_params.mstarcv;
                        dT10_dVb = -dVth_dVb * self.size_params.mstarcv;
                    } else if T2 < -EXP_THRESHOLD {
                        T10 = Vtm * log(1.0 + MIN_EXP);
                        dT10_dVg = 0.0;
                        dT10_dVd = T10 * dn_dVd;
                        dT10_dVb = T10 * dn_dVb;
                        T10 *= n;
                    } else {
                        ExpVgst = exp(T2);
                        T3 = Vtm * log(1.0 + ExpVgst);
                        T10 = n * T3;
                        dT10_dVg = self.size_params.mstarcv * ExpVgst / (1.0 + ExpVgst);
                        dT10_dVb = T3 * dn_dVb - dT10_dVg * (dVth_dVb + Vgst * dn_dVb / n);
                        dT10_dVd = T3 * dn_dVd - dT10_dVg * (dVth_dVd + Vgst * dn_dVd / n);
                        dT10_dVg *= dVgs_eff_dVg;
                    }

                    T1 = self.size_params.voffcbncv - (1.0 - self.size_params.mstarcv) * Vgst;
                    T2 = T1 / T0;
                    if T2 < -EXP_THRESHOLD {
                        T3 = self.intp.coxe * MIN_EXP / self.size_params.cdep0;
                        T9 = self.size_params.mstarcv + T3 * n;
                        dT9_dVg = 0.0;
                        dT9_dVd = dn_dVd * T3;
                        dT9_dVb = dn_dVb * T3;
                    } else if T2 > EXP_THRESHOLD {
                        T3 = self.intp.coxe * MAX_EXP / self.size_params.cdep0;
                        T9 = self.size_params.mstarcv + T3 * n;
                        dT9_dVg = 0.0;
                        dT9_dVd = dn_dVd * T3;
                        dT9_dVb = dn_dVb * T3;
                    } else {
                        ExpVgst = exp(T2);
                        T3 = self.intp.coxe / self.size_params.cdep0;
                        T4 = T3 * ExpVgst;
                        T5 = T1 * T4 / T0;
                        T9 = self.size_params.mstarcv + n * T4;
                        dT9_dVg = T3 * (self.size_params.mstarcv - 1.0) * ExpVgst / Vtm;
                        dT9_dVb = T4 * dn_dVb - dT9_dVg * dVth_dVb - T5 * dn_dVb;
                        dT9_dVd = T4 * dn_dVd - dT9_dVg * dVth_dVd - T5 * dn_dVd;
                        dT9_dVg *= dVgs_eff_dVg;
                    }

                    Vgsteff = T10 / T9;
                    T11 = T9 * T9;
                    dVgsteff_dVg = (T9 * dT10_dVg - T10 * dT9_dVg) / T11;
                    dVgsteff_dVd = (T9 * dT10_dVd - T10 * dT9_dVd) / T11;
                    dVgsteff_dVb = (T9 * dT10_dVb - T10 * dT9_dVb) / T11;
                    /* End of VgsteffCV for cvchargeMod = 1 */
                }

                if self.model.capmod == 1 {
                    Vfb = self.intp.vfbzb;
                    V3 = Vfb - Vgs_eff + VbseffCV - DELTA_3;
                    if (Vfb <= 0.0) {
                        T0 = sqrt(V3 * V3 - 4.0 * DELTA_3 * Vfb);
                    } else {
                        T0 = sqrt(V3 * V3 + 4.0 * DELTA_3 * Vfb);
                    }

                    T1 = 0.5 * (1.0 + V3 / T0);
                    Vfbeff = Vfb - 0.5 * (V3 + T0);
                    dVfbeff_dVg = T1 * dVgs_eff_dVg;
                    dVfbeff_dVb = -T1 * dVbseffCV_dVb;
                    Qac0 = CoxWL * (Vfbeff - Vfb);
                    dQac0_dVg = CoxWL * dVfbeff_dVg;
                    dQac0_dVb = CoxWL * dVfbeff_dVb;

                    T0 = 0.5 * self.size_params.k1ox;
                    T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
                    if self.size_params.k1ox == 0.0 {
                        T1 = 0.0;
                        T2 = 0.0;
                    } else if T3 < 0.0 {
                        T1 = T0 + T3 / self.size_params.k1ox;
                        T2 = CoxWL;
                    } else {
                        T1 = sqrt(T0 * T0 + T3);
                        T2 = CoxWL * T0 / T1;
                    }

                    Qsub0 = CoxWL * self.size_params.k1ox * (T1 - T0);

                    dQsub0_dVg = T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg);
                    dQsub0_dVd = -T2 * dVgsteff_dVd;
                    dQsub0_dVb = -T2 * (dVfbeff_dVb + dVbseffCV_dVb + dVgsteff_dVb);

                    AbulkCV = Abulk0 * self.size_params.abulkCVfactor;
                    dAbulkCV_dVb = self.size_params.abulkCVfactor * dAbulk0_dVb;
                    VdsatCV = Vgsteff / AbulkCV;

                    T0 = VdsatCV - Vds - DELTA_4;
                    dT0_dVg = 1.0 / AbulkCV;
                    dT0_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV;
                    T1 = sqrt(T0 * T0 + 4.0 * DELTA_4 * VdsatCV);
                    dT1_dVg = (T0 + DELTA_4 + DELTA_4) / T1;
                    dT1_dVd = -T0 / T1;
                    dT1_dVb = dT1_dVg * dT0_dVb;
                    dT1_dVg *= dT0_dVg;
                    if T0 >= 0.0 {
                        VdseffCV = VdsatCV - 0.5 * (T0 + T1);
                        dVdseffCV_dVg = 0.5 * (dT0_dVg - dT1_dVg);
                        dVdseffCV_dVd = 0.5 * (1.0 - dT1_dVd);
                        dVdseffCV_dVb = 0.5 * (dT0_dVb - dT1_dVb);
                    } else {
                        T3 = (DELTA_4 + DELTA_4) / (T1 - T0);
                        T4 = 1.0 - T3;
                        T5 = VdsatCV * T3 / (T1 - T0);
                        VdseffCV = VdsatCV * T4;
                        dVdseffCV_dVg = dT0_dVg * T4 + T5 * (dT1_dVg - dT0_dVg);
                        dVdseffCV_dVd = T5 * (dT1_dVd + 1.0);
                        dVdseffCV_dVb = dT0_dVb * (T4 - T5) + T5 * dT1_dVb;
                    }

                    if Vds == 0.0 {
                        VdseffCV = 0.0;
                        dVdseffCV_dVg = 0.0;
                        dVdseffCV_dVb = 0.0;
                    }

                    T0 = AbulkCV * VdseffCV;
                    T1 = 12.0 * (Vgsteff - 0.5 * T0 + 1.0e-20);
                    T2 = VdseffCV / T1;
                    T3 = T0 * T2;

                    T4 = (1.0 - 12.0 * T2 * T2 * AbulkCV);
                    T5 = (6.0 * T0 * (4.0 * Vgsteff - T0) / (T1 * T1) - 0.5);
                    T6 = 12.0 * T2 * T2 * Vgsteff;

                    qgate = CoxWL * (Vgsteff - 0.5 * VdseffCV + T3);
                    Cgg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
                    Cgd1 = CoxWL * T5 * dVdseffCV_dVd + Cgg1 * dVgsteff_dVd;
                    Cgb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb) + Cgg1 * dVgsteff_dVb;
                    Cgg1 *= dVgsteff_dVg;

                    T7 = 1.0 - AbulkCV;
                    qbulk = CoxWL * T7 * (0.5 * VdseffCV - T3);
                    T4 = -T7 * (T4 - 1.0);
                    T5 = -T7 * T5;
                    T6 = -(T7 * T6 + (0.5 * VdseffCV - T3));
                    Cbg1 = CoxWL * (T4 + T5 * dVdseffCV_dVg);
                    Cbd1 = CoxWL * T5 * dVdseffCV_dVd + Cbg1 * dVgsteff_dVd;
                    Cbb1 = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb) + Cbg1 * dVgsteff_dVb;
                    Cbg1 *= dVgsteff_dVg;

                    if self.model.xpart > 0.5 {
                        /* 0/100 Charge petition model */
                        T1 = T1 + T1;
                        qsrc = -CoxWL * (0.5 * Vgsteff + 0.25 * T0 - T0 * T0 / T1);
                        T7 = (4.0 * Vgsteff - T0) / (T1 * T1);
                        T4 = -(0.5 + 24.0 * T0 * T0 / (T1 * T1));
                        T5 = -(0.25 * AbulkCV - 12.0 * AbulkCV * T0 * T7);
                        T6 = -(0.25 * VdseffCV - 12.0 * T0 * VdseffCV * T7);
                        Csg = CoxWL * (T4 + T5 * dVdseffCV_dVg);
                        Csd = CoxWL * T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
                        Csb = CoxWL * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb) + Csg * dVgsteff_dVb;
                        Csg *= dVgsteff_dVg;
                    } else if self.model.xpart < 0.5 {
                        /* 40/60 Charge petition model */
                        T1 = T1 / 12.0;
                        T2 = 0.5 * CoxWL / (T1 * T1);
                        T3 = Vgsteff * (2.0 * T0 * T0 / 3.0 + Vgsteff * (Vgsteff - 4.0 * T0 / 3.0))
                            - 2.0 * T0 * T0 * T0 / 15.0;
                        qsrc = -T2 * T3;
                        T7 = 4.0 / 3.0 * Vgsteff * (Vgsteff - T0) + 0.4 * T0 * T0;
                        T4 = -2.0 * qsrc / T1
                            - T2 * (Vgsteff * (3.0 * Vgsteff - 8.0 * T0 / 3.0)
                                + 2.0 * T0 * T0 / 3.0);
                        T5 = (qsrc / T1 + T2 * T7) * AbulkCV;
                        T6 = (qsrc / T1 * VdseffCV + T2 * T7 * VdseffCV);
                        Csg = (T4 + T5 * dVdseffCV_dVg);
                        Csd = T5 * dVdseffCV_dVd + Csg * dVgsteff_dVd;
                        Csb = (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb) + Csg * dVgsteff_dVb;
                        Csg *= dVgsteff_dVg;
                    } else {
                        /* 50/50 Charge petition model */
                        qsrc = -0.5 * (qgate + qbulk);
                        Csg = -0.5 * (Cgg1 + Cbg1);
                        Csb = -0.5 * (Cgb1 + Cbb1);
                        Csd = -0.5 * (Cgd1 + Cbd1);
                    }

                    qgate += Qac0 + Qsub0;
                    qbulk -= (Qac0 + Qsub0);
                    qdrn = -(qgate + qbulk + qsrc);

                    Cgg = dQac0_dVg + dQsub0_dVg + Cgg1;
                    Cgd = dQsub0_dVd + Cgd1;
                    Cgb = dQac0_dVb + dQsub0_dVb + Cgb1;

                    Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
                    Cbd = Cbd1 - dQsub0_dVd;
                    Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

                    Cgb *= dVbseff_dVb;
                    Cbb *= dVbseff_dVb;
                    Csb *= dVbseff_dVb;

                    newop.cggb = Cgg;
                    newop.cgsb = -(Cgg + Cgd + Cgb);
                    newop.cgdb = Cgd;
                    newop.cdgb = -(Cgg + Cbg + Csg);
                    newop.cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb + Csg + Csd + Csb);
                    newop.cddb = -(Cgd + Cbd + Csd);
                    newop.cbgb = Cbg;
                    newop.cbsb = -(Cbg + Cbd + Cbb);
                    newop.cbdb = Cbd;
                }
                /* Charge-Thickness capmod (CTM) begins */
                else if self.model.capmod == 2 {
                    V3 = self.intp.vfbzb - Vgs_eff + VbseffCV - DELTA_3;
                    if (self.intp.vfbzb <= 0.0) {
                        T0 = sqrt(V3 * V3 - 4.0 * DELTA_3 * self.intp.vfbzb);
                    } else {
                        T0 = sqrt(V3 * V3 + 4.0 * DELTA_3 * self.intp.vfbzb);
                    }

                    T1 = 0.5 * (1.0 + V3 / T0);
                    Vfbeff = self.intp.vfbzb - 0.5 * (V3 + T0);
                    dVfbeff_dVg = T1 * dVgs_eff_dVg;
                    dVfbeff_dVb = -T1 * dVbseffCV_dVb;

                    Cox = self.intp.coxp;
                    Tox = 1.0e8 * self.intp.toxp;
                    T0 = (Vgs_eff - VbseffCV - self.intp.vfbzb) / Tox;
                    dT0_dVg = dVgs_eff_dVg / Tox;
                    dT0_dVb = -dVbseffCV_dVb / Tox;

                    tmp = T0 * self.size_params.acde;
                    if (-EXP_THRESHOLD < tmp) && (tmp < EXP_THRESHOLD) {
                        Tcen = self.size_params.ldeb * exp(tmp);
                        dTcen_dVg = self.size_params.acde * Tcen;
                        dTcen_dVb = dTcen_dVg * dT0_dVb;
                        dTcen_dVg *= dT0_dVg;
                    } else if tmp <= -EXP_THRESHOLD {
                        Tcen = self.size_params.ldeb * MIN_EXP;
                        dTcen_dVg = 0.0;
                        dTcen_dVb = 0.0;
                    } else {
                        Tcen = self.size_params.ldeb * MAX_EXP;
                        dTcen_dVg = 0.0;
                        dTcen_dVb = 0.0;
                    }

                    LINK = 1.0e-3 * self.intp.toxp;
                    V3 = self.size_params.ldeb - Tcen - LINK;
                    V4 = sqrt(V3 * V3 + 4.0 * LINK * self.size_params.ldeb);
                    Tcen = self.size_params.ldeb - 0.5 * (V3 + V4);
                    T1 = 0.5 * (1.0 + V3 / V4);
                    dTcen_dVg *= T1;
                    dTcen_dVb *= T1;

                    Ccen = epssub / Tcen;
                    T2 = Cox / (Cox + Ccen);
                    Coxeff = T2 * Ccen;
                    T3 = -Ccen / Tcen;
                    dCoxeff_dVg = T2 * T2 * T3;
                    dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
                    dCoxeff_dVg *= dTcen_dVg;
                    CoxWLcen = CoxWL * Coxeff / self.intp.coxe;

                    Qac0 = CoxWLcen * (Vfbeff - self.intp.vfbzb);
                    QovCox = Qac0 / Coxeff;
                    dQac0_dVg = CoxWLcen * dVfbeff_dVg + QovCox * dCoxeff_dVg;
                    dQac0_dVb = CoxWLcen * dVfbeff_dVb + QovCox * dCoxeff_dVb;

                    T0 = 0.5 * self.size_params.k1ox;
                    T3 = Vgs_eff - Vfbeff - VbseffCV - Vgsteff;
                    if self.size_params.k1ox == 0.0 {
                        T1 = 0.0;
                        T2 = 0.0;
                    } else if T3 < 0.0 {
                        T1 = T0 + T3 / self.size_params.k1ox;
                        T2 = CoxWLcen;
                    } else {
                        T1 = sqrt(T0 * T0 + T3);
                        T2 = CoxWLcen * T0 / T1;
                    }

                    Qsub0 = CoxWLcen * self.size_params.k1ox * (T1 - T0);
                    QovCox = Qsub0 / Coxeff;
                    dQsub0_dVg =
                        T2 * (dVgs_eff_dVg - dVfbeff_dVg - dVgsteff_dVg) + QovCox * dCoxeff_dVg;
                    dQsub0_dVd = -T2 * dVgsteff_dVd;
                    dQsub0_dVb =
                        -T2 * (dVfbeff_dVb + dVbseffCV_dVb + dVgsteff_dVb) + QovCox * dCoxeff_dVb;

                    /* Gate-bias dependent delta Phis begins */
                    if self.size_params.k1ox <= 0.0 {
                        Denomi = 0.25 * self.size_params.moin * Vtm;
                        T0 = 0.5 * self.size_params.sqrtPhi;
                    } else {
                        Denomi = self.size_params.moin
                            * Vtm
                            * self.size_params.k1ox
                            * self.size_params.k1ox;
                        T0 = self.size_params.k1ox * self.size_params.sqrtPhi;
                    }
                    T1 = 2.0 * T0 + Vgsteff;

                    DeltaPhi = Vtm * log(1.0 + T1 * Vgsteff / Denomi);
                    dDeltaPhi_dVg = 2.0 * Vtm * (T1 - T0) / (Denomi + T1 * Vgsteff);
                    /* End of delta Phis */

                    /* VgDP = Vgsteff - DeltaPhi */
                    T0 = Vgsteff - DeltaPhi - 0.001;
                    dT0_dVg = 1.0 - dDeltaPhi_dVg;
                    T1 = sqrt(T0 * T0 + Vgsteff * 0.004);
                    VgDP = 0.5 * (T0 + T1);
                    dVgDP_dVg = 0.5 * (dT0_dVg + (T0 * dT0_dVg + 0.002) / T1);

                    Tox += Tox;
                    T0 = (Vgsteff + self.intp.vtfbphi2) / Tox;
                    tmp = exp(self.model.bdos * 0.7 * log(T0));
                    T1 = 1.0 + tmp;
                    T2 = self.model.bdos * 0.7 * tmp / (T0 * Tox);
                    Tcen = self.model.ados * 1.9e-9 / T1;
                    dTcen_dVg = -Tcen * T2 / T1;
                    dTcen_dVd = dTcen_dVg * dVgsteff_dVd;
                    dTcen_dVb = dTcen_dVg * dVgsteff_dVb;
                    dTcen_dVg *= dVgsteff_dVg;

                    Ccen = epssub / Tcen;
                    T0 = Cox / (Cox + Ccen);
                    Coxeff = T0 * Ccen;
                    T1 = -Ccen / Tcen;
                    dCoxeff_dVg = T0 * T0 * T1;
                    dCoxeff_dVd = dCoxeff_dVg * dTcen_dVd;
                    dCoxeff_dVb = dCoxeff_dVg * dTcen_dVb;
                    dCoxeff_dVg *= dTcen_dVg;
                    CoxWLcen = CoxWL * Coxeff / self.intp.coxe;

                    AbulkCV = Abulk0 * self.size_params.abulkCVfactor;
                    dAbulkCV_dVb = self.size_params.abulkCVfactor * dAbulk0_dVb;
                    VdsatCV = VgDP / AbulkCV;

                    T0 = VdsatCV - Vds - DELTA_4;
                    dT0_dVg = dVgDP_dVg / AbulkCV;
                    dT0_dVb = -VdsatCV * dAbulkCV_dVb / AbulkCV;
                    T1 = sqrt(T0 * T0 + 4.0 * DELTA_4 * VdsatCV);
                    dT1_dVg = (T0 + DELTA_4 + DELTA_4) / T1;
                    dT1_dVd = -T0 / T1;
                    dT1_dVb = dT1_dVg * dT0_dVb;
                    dT1_dVg *= dT0_dVg;
                    if T0 >= 0.0 {
                        VdseffCV = VdsatCV - 0.5 * (T0 + T1);
                        dVdseffCV_dVg = 0.5 * (dT0_dVg - dT1_dVg);
                        dVdseffCV_dVd = 0.5 * (1.0 - dT1_dVd);
                        dVdseffCV_dVb = 0.5 * (dT0_dVb - dT1_dVb);
                    } else {
                        T3 = (DELTA_4 + DELTA_4) / (T1 - T0);
                        T4 = 1.0 - T3;
                        T5 = VdsatCV * T3 / (T1 - T0);
                        VdseffCV = VdsatCV * T4;
                        dVdseffCV_dVg = dT0_dVg * T4 + T5 * (dT1_dVg - dT0_dVg);
                        dVdseffCV_dVd = T5 * (dT1_dVd + 1.0);
                        dVdseffCV_dVb = dT0_dVb * (T4 - T5) + T5 * dT1_dVb;
                    }

                    if Vds == 0.0 {
                        VdseffCV = 0.0;
                        dVdseffCV_dVg = 0.0;
                        dVdseffCV_dVb = 0.0;
                    }

                    T0 = AbulkCV * VdseffCV;
                    T1 = VgDP;
                    T2 = 12.0 * (T1 - 0.5 * T0 + 1.0e-20);
                    T3 = T0 / T2;
                    T4 = 1.0 - 12.0 * T3 * T3;
                    T5 = AbulkCV * (6.0 * T0 * (4.0 * T1 - T0) / (T2 * T2) - 0.5);
                    T6 = T5 * VdseffCV / AbulkCV;

                    qgate = CoxWLcen * (T1 - T0 * (0.5 - T3));
                    QovCox = qgate / Coxeff;
                    Cgg1 = CoxWLcen * (T4 * dVgDP_dVg + T5 * dVdseffCV_dVg);
                    Cgd1 =
                        CoxWLcen * T5 * dVdseffCV_dVd + Cgg1 * dVgsteff_dVd + QovCox * dCoxeff_dVd;
                    Cgb1 = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                        + Cgg1 * dVgsteff_dVb
                        + QovCox * dCoxeff_dVb;
                    Cgg1 = Cgg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;

                    T7 = 1.0 - AbulkCV;
                    T8 = T2 * T2;
                    T9 = 12.0 * T7 * T0 * T0 / (T8 * AbulkCV);
                    T10 = T9 * dVgDP_dVg;
                    T11 = -T7 * T5 / AbulkCV;
                    T12 = -(T9 * T1 / AbulkCV + VdseffCV * (0.5 - T0 / T2));

                    qbulk = CoxWLcen * T7 * (0.5 * VdseffCV - T0 * VdseffCV / T2);
                    QovCox = qbulk / Coxeff;
                    Cbg1 = CoxWLcen * (T10 + T11 * dVdseffCV_dVg);
                    Cbd1 =
                        CoxWLcen * T11 * dVdseffCV_dVd + Cbg1 * dVgsteff_dVd + QovCox * dCoxeff_dVd;
                    Cbb1 = CoxWLcen * (T11 * dVdseffCV_dVb + T12 * dAbulkCV_dVb)
                        + Cbg1 * dVgsteff_dVb
                        + QovCox * dCoxeff_dVb;
                    Cbg1 = Cbg1 * dVgsteff_dVg + QovCox * dCoxeff_dVg;

                    if self.model.xpart > 0.5 {
                        /* 0/100 partition */
                        qsrc = -CoxWLcen * (T1 / 2.0 + T0 / 4.0 - 0.5 * T0 * T0 / T2);
                        QovCox = qsrc / Coxeff;
                        T2 += T2;
                        T3 = T2 * T2;
                        T7 = -(0.25 - 12.0 * T0 * (4.0 * T1 - T0) / T3);
                        T4 = -(0.5 + 24.0 * T0 * T0 / T3) * dVgDP_dVg;
                        T5 = T7 * AbulkCV;
                        T6 = T7 * VdseffCV;

                        Csg = CoxWLcen * (T4 + T5 * dVdseffCV_dVg);
                        Csd = CoxWLcen * T5 * dVdseffCV_dVd
                            + Csg * dVgsteff_dVd
                            + QovCox * dCoxeff_dVd;
                        Csb = CoxWLcen * (T5 * dVdseffCV_dVb + T6 * dAbulkCV_dVb)
                            + Csg * dVgsteff_dVb
                            + QovCox * dCoxeff_dVb;
                        Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
                    } else if self.model.xpart < 0.5 {
                        /* 40/60 partition */
                        T2 = T2 / 12.0;
                        T3 = 0.5 * CoxWLcen / (T2 * T2);
                        T4 = T1 * (2.0 * T0 * T0 / 3.0 + T1 * (T1 - 4.0 * T0 / 3.0))
                            - 2.0 * T0 * T0 * T0 / 15.0;
                        qsrc = -T3 * T4;
                        QovCox = qsrc / Coxeff;
                        T8 = 4.0 / 3.0 * T1 * (T1 - T0) + 0.4 * T0 * T0;
                        T5 = -2.0 * qsrc / T2
                            - T3 * (T1 * (3.0 * T1 - 8.0 * T0 / 3.0) + 2.0 * T0 * T0 / 3.0);
                        T6 = AbulkCV * (qsrc / T2 + T3 * T8);
                        T7 = T6 * VdseffCV / AbulkCV;

                        Csg = T5 * dVgDP_dVg + T6 * dVdseffCV_dVg;
                        Csd = Csg * dVgsteff_dVd + T6 * dVdseffCV_dVd + QovCox * dCoxeff_dVd;
                        Csb = Csg * dVgsteff_dVb
                            + T6 * dVdseffCV_dVb
                            + T7 * dAbulkCV_dVb
                            + QovCox * dCoxeff_dVb;
                        Csg = Csg * dVgsteff_dVg + QovCox * dCoxeff_dVg;
                    } else {
                        /* 50/50 partition */
                        qsrc = -0.5 * qgate;
                        Csg = -0.5 * Cgg1;
                        Csd = -0.5 * Cgd1;
                        Csb = -0.5 * Cgb1;
                    }

                    qgate += Qac0 + Qsub0 - qbulk;
                    qbulk -= (Qac0 + Qsub0);
                    qdrn = -(qgate + qbulk + qsrc);

                    Cbg = Cbg1 - dQac0_dVg - dQsub0_dVg;
                    Cbd = Cbd1 - dQsub0_dVd;
                    Cbb = Cbb1 - dQac0_dVb - dQsub0_dVb;

                    Cgg = Cgg1 - Cbg;
                    Cgd = Cgd1 - Cbd;
                    Cgb = Cgb1 - Cbb;

                    Cgb *= dVbseff_dVb;
                    Cbb *= dVbseff_dVb;
                    Csb *= dVbseff_dVb;

                    newop.cggb = Cgg;
                    newop.cgsb = -(Cgg + Cgd + Cgb);
                    newop.cgdb = Cgd;
                    newop.cdgb = -(Cgg + Cbg + Csg);
                    newop.cdsb = (Cgg + Cgd + Cgb + Cbg + Cbd + Cbb + Csg + Csd + Csb);
                    newop.cddb = -(Cgd + Cbd + Csd);
                    newop.cbgb = Cbg;
                    newop.cbsb = -(Cbg + Cbd + Cbb);
                    newop.cbdb = Cbd;
                } /* End of CTM */
            }

            newop.csgb = -newop.cggb - newop.cdgb - newop.cbgb;
            newop.csdb = -newop.cgdb - newop.cddb - newop.cbdb;
            newop.cssb = -newop.cgsb - newop.cdsb - newop.cbsb;
            newop.cgbb = -newop.cgdb - newop.cggb - newop.cgsb;
            newop.cdbb = -newop.cddb - newop.cdgb - newop.cdsb;
            newop.cbbb = -newop.cbgb - newop.cbdb - newop.cbsb;
            newop.csbb = -newop.cgbb - newop.cdbb - newop.cbbb;
            newop.qgate = qgate;
            newop.qbulk = qbulk;
            newop.qdrn = qdrn;
            newop.qsrc = -(qgate + qbulk + qdrn);

            /* NQS begins */
            if (self.model.trnqsmod != 0) || (self.model.acnqsmod != 0) {
                qcheq = -(qbulk + qgate);
                newop.qchqs = qcheq;
                newop.cqgb = -(newop.cggb + newop.cbgb);
                newop.cqdb = -(newop.cgdb + newop.cbdb);
                newop.cqsb = -(newop.cgsb + newop.cbsb);
                newop.cqbb = -(newop.cqgb + newop.cqdb + newop.cqsb);

                CoxWL = self.intp.coxe
                    * self.size_params.weffCV
                    * self.intp.nf
                    * self.size_params.leffCV;
                T1 = newop.gcrg / CoxWL; /* 1 / tau */
                newop.gtau = T1 * ScalingFactor;

                if self.model.acnqsmod != 0 {
                    newop.taunet = 1.0 / T1;
                }

                newop.qcheq = qcheq;

                if self.model.trnqsmod != 0 {
                    let mut _tmp1: f64;
                    let mut _tmp2: f64;
                    let (_g, i, _r) =
                        state.integrate(newop.qcheq - self.op.qcheq, 0.0, 0.0, self.op.cqcheq);
                    newop.cqcheq = i;
                    // error = NIintegrate(ckt, &geq, &ceq, 0.0, newop.qcheq);
                }
            }
        }

        // Charge computations
        if let AnalysisInfo::TRAN(_, state) = an {
            czbd = self.intp.DunitAreaTempJctCap * self.intp.Adeff;
            czbs = self.intp.SunitAreaTempJctCap * self.intp.Aseff;
            czbdsw = self.intp.DunitLengthSidewallTempJctCap * self.intp.Pdeff;
            czbdswg = self.intp.DunitLengthGateSidewallTempJctCap
                * self.size_params.weffCJ
                * self.intp.nf;
            czbssw = self.intp.SunitLengthSidewallTempJctCap * self.intp.Pseff;
            czbsswg = self.intp.SunitLengthGateSidewallTempJctCap
                * self.size_params.weffCJ
                * self.intp.nf;

            MJS = self.intp.SbulkJctBotGradingCoeff;
            MJSWS = self.intp.SbulkJctSideGradingCoeff;
            MJSWGS = self.intp.SbulkJctGateSideGradingCoeff;

            MJD = self.intp.DbulkJctBotGradingCoeff;
            MJSWD = self.intp.DbulkJctSideGradingCoeff;
            MJSWGD = self.intp.DbulkJctGateSideGradingCoeff;

            /* Source Bulk Junction */
            if vbs_jct == 0.0 {
                newop.qbs = 0.0;
                newop.capbs = czbs + czbssw + czbsswg;
            } else if vbs_jct < 0.0 {
                if czbs > 0.0 {
                    arg = 1.0 - vbs_jct / self.intp.PhiBS;
                    if (MJS == 0.5) {
                        sarg = 1.0 / sqrt(arg);
                    } else {
                        sarg = exp(-MJS * log(arg));
                    }
                    newop.qbs = self.intp.PhiBS * czbs * (1.0 - arg * sarg) / (1.0 - MJS);
                    newop.capbs = czbs * sarg;
                } else {
                    newop.qbs = 0.0;
                    newop.capbs = 0.0;
                }
                if czbssw > 0.0 {
                    arg = 1.0 - vbs_jct / self.intp.PhiBSWS;
                    if (MJSWS == 0.5) {
                        sarg = 1.0 / sqrt(arg);
                    } else {
                        sarg = exp(-MJSWS * log(arg));
                    }
                    newop.qbs += self.intp.PhiBSWS * czbssw * (1.0 - arg * sarg) / (1.0 - MJSWS);
                    newop.capbs += czbssw * sarg;
                }
                if czbsswg > 0.0 {
                    arg = 1.0 - vbs_jct / self.intp.PhiBSWGS;
                    if (MJSWGS == 0.5) {
                        sarg = 1.0 / sqrt(arg);
                    } else {
                        sarg = exp(-MJSWGS * log(arg));
                    }
                    newop.qbs += self.intp.PhiBSWGS * czbsswg * (1.0 - arg * sarg) / (1.0 - MJSWGS);
                    newop.capbs += czbsswg * sarg;
                }
            } else {
                T0 = czbs + czbssw + czbsswg;
                T1 = vbs_jct
                    * (czbs * MJS / self.intp.PhiBS
                        + czbssw * MJSWS / self.intp.PhiBSWS
                        + czbsswg * MJSWGS / self.intp.PhiBSWGS);
                newop.qbs = vbs_jct * (T0 + 0.5 * T1);
                newop.capbs = T0 + T1;
            }

            /* Drain Bulk Junction */
            if vbd_jct == 0.0 {
                newop.qbd = 0.0;
                newop.capbd = czbd + czbdsw + czbdswg;
            } else if vbd_jct < 0.0 {
                if czbd > 0.0 {
                    arg = 1.0 - vbd_jct / self.intp.PhiBD;
                    if (MJD == 0.5) {
                        sarg = 1.0 / sqrt(arg);
                    } else {
                        sarg = exp(-MJD * log(arg));
                    }
                    newop.qbd = self.intp.PhiBD * czbd * (1.0 - arg * sarg) / (1.0 - MJD);
                    newop.capbd = czbd * sarg;
                } else {
                    newop.qbd = 0.0;
                    newop.capbd = 0.0;
                }
                if czbdsw > 0.0 {
                    arg = 1.0 - vbd_jct / self.intp.PhiBSWD;
                    if (MJSWD == 0.5) {
                        sarg = 1.0 / sqrt(arg);
                    } else {
                        sarg = exp(-MJSWD * log(arg));
                    }
                    newop.qbd += self.intp.PhiBSWD * czbdsw * (1.0 - arg * sarg) / (1.0 - MJSWD);
                    newop.capbd += czbdsw * sarg;
                }
                if czbdswg > 0.0 {
                    arg = 1.0 - vbd_jct / self.intp.PhiBSWGD;
                    if (MJSWGD == 0.5) {
                        sarg = 1.0 / sqrt(arg);
                    } else {
                        sarg = exp(-MJSWGD * log(arg));
                    }
                    newop.qbd += self.intp.PhiBSWGD * czbdswg * (1.0 - arg * sarg) / (1.0 - MJSWGD);
                    newop.capbd += czbdswg * sarg;
                }
            } else {
                T0 = czbd + czbdsw + czbdswg;
                T1 = vbd_jct
                    * (czbd * MJD / self.intp.PhiBD
                        + czbdsw * MJSWD / self.intp.PhiBSWD
                        + czbdswg * MJSWGD / self.intp.PhiBSWGD);
                newop.qbd = vbd_jct * (T0 + 0.5 * T1);
                newop.capbd = T0 + T1;
            }
        }

        newop.vds = vds;
        newop.vgs = vgs;
        newop.vbs = vbs;
        newop.vbd = vbd;
        newop.vges = vges;
        newop.vgms = vgms;
        newop.vdbs = vdbs;
        newop.vdbd = vdbd;
        newop.vsbs = vsbs;
        newop.vses = vses;
        newop.vdes = vdes;
        newop.qdef = qdef;

        // Initially zero all capacitances and their impedances
        // Many complicated paths through the code below do not assure they are otherwise initialized.
        let mut ceqqg = 0.0;
        let mut ceqqb = 0.0;
        let mut ceqqd = 0.0;
        let mut ceqqjd = 0.0;
        let mut ceqqjs = 0.0;
        let mut cqcheq = 0.0;
        let mut cqdef = 0.0;

        let mut gcdgb = 0.0;
        let mut gcddb = 0.0;
        let mut gcdsb = 0.0;
        let mut gcdbb = 0.0;
        let mut gcsgb = 0.0;
        let mut gcsdb = 0.0;
        let mut gcssb = 0.0;
        let mut gcsbb = 0.0;
        let mut gcggb = 0.0;
        let mut gcgdb = 0.0;
        let mut gcgsb = 0.0;
        let mut gcgbb = 0.0;
        let mut gcbdb = 0.0;
        let mut gcbgb = 0.0;
        let mut gcbsb = 0.0;
        let mut gcbbb = 0.0;

        let mut gcgmgmb = 0.0;
        let mut gcgmdb = 0.0;
        let mut gcgmsb = 0.0;
        let mut gcgmbb = 0.0;
        let mut gcdgmb = 0.0;
        let mut gcsgmb = 0.0;
        let mut gcbgmb = 0.0;
        let mut ceqqgmid = 0.0;
        let mut gcdbdb = 0.0;
        let mut gcsbsb = 0.0;

        let mut gqdef = 0.0;
        let mut gcqgb = 0.0;
        let mut gcqdb = 0.0;
        let mut gcqsb = 0.0;
        let mut gcqbb = 0.0;
        let mut ggtg = 0.0;
        let mut ggtd = 0.0;
        let mut ggtb = 0.0;
        let mut ggts = 0.0;
        let mut dxpart = if self.guess.mode > 0 { 0.4 } else { 0.6 };
        let mut sxpart = 1.0 - dxpart;
        let mut ddxpart_dVd = 0.0;
        let mut ddxpart_dVg = 0.0;
        let mut ddxpart_dVb = 0.0;
        let mut ddxpart_dVs = 0.0;
        let mut dsxpart_dVd = 0.0;
        let mut dsxpart_dVg = 0.0;
        let mut dsxpart_dVb = 0.0;
        let mut dsxpart_dVs = 0.0;

        if self.model.trnqsmod != 0 {
            CoxWL =
                self.intp.coxe * self.size_params.weffCV * self.intp.nf * self.size_params.leffCV;
            T1 = newop.gcrg / CoxWL;
            newop.gtau = T1 * ScalingFactor;
        } else {
            newop.gtau = 0.0;
        }

        // Charge computations
        if let AnalysisInfo::TRAN(_, state) = an {
            if self.intp.rgateMod == 3 {
                vgdx = vgmd;
                vgsx = vgms;
            } else {
                vgdx = vgd;
                vgsx = vgs;
            }
            if self.model.capmod == 0 {
                cgdo = self.size_params.cgdo;
                qgdo = self.size_params.cgdo * vgdx;
                cgso = self.size_params.cgso;
                qgso = self.size_params.cgso * vgsx;
            } else {
                T0 = vgdx + DELTA_1;
                T1 = sqrt(T0 * T0 + 4.0 * DELTA_1);
                T2 = 0.5 * (T0 - T1);

                T3 = self.size_params.weffCV * self.size_params.cgdl;
                T4 = sqrt(1.0 - 4.0 * T2 / self.size_params.ckappad);
                cgdo = self.size_params.cgdo + T3 - T3 * (1.0 - 1.0 / T4) * (0.5 - 0.5 * T0 / T1);
                qgdo = (self.size_params.cgdo + T3) * vgdx
                    - T3 * (T2 + 0.5 * self.size_params.ckappad * (T4 - 1.0));

                T0 = vgsx + DELTA_1;
                T1 = sqrt(T0 * T0 + 4.0 * DELTA_1);
                T2 = 0.5 * (T0 - T1);
                T3 = self.size_params.weffCV * self.size_params.cgsl;
                T4 = sqrt(1.0 - 4.0 * T2 / self.size_params.ckappas);
                cgso = self.size_params.cgso + T3 - T3 * (1.0 - 1.0 / T4) * (0.5 - 0.5 * T0 / T1);
                qgso = (self.size_params.cgso + T3) * vgsx
                    - T3 * (T2 + 0.5 * self.size_params.ckappas * (T4 - 1.0));
            }

            if self.intp.nf != 1.0 {
                cgdo *= self.intp.nf;
                cgso *= self.intp.nf;
                qgdo *= self.intp.nf;
                qgso *= self.intp.nf;
            }
            newop.cgdo = cgdo;
            newop.qgdo = qgdo;
            newop.cgso = cgso;
            newop.qgso = qgso;

            // TODO: the BSIM4 reference implementation essentially bakes numerical integration in here,
            // ignoring the circuit/ analysis integration method.
            // All of these impedances are calculated as g = C/dt, e.g. using Backward Euler.
            // Should figure out whether this is the implementation intent, or just for reference.
            ag0 = 1.0 / state.dt;
            if self.guess.mode > 0 {
                if self.model.trnqsmod == 0 {
                    qdrn -= qgdo;
                    if self.intp.rgateMod == 3 {
                        gcgmgmb = (cgdo + cgso + self.size_params.cgbo) * ag0;
                        gcgmdb = -cgdo * ag0;
                        gcgmsb = -cgso * ag0;
                        gcgmbb = -self.size_params.cgbo * ag0;

                        gcdgmb = gcgmdb;
                        gcsgmb = gcgmsb;
                        gcbgmb = gcgmbb;

                        gcggb = newop.cggb * ag0;
                        gcgdb = newop.cgdb * ag0;
                        gcgsb = newop.cgsb * ag0;
                        gcgbb = -(gcggb + gcgdb + gcgsb);

                        gcdgb = newop.cdgb * ag0;
                        gcsgb = -(newop.cggb + newop.cbgb + newop.cdgb) * ag0;
                        gcbgb = newop.cbgb * ag0;

                        qgmb = self.size_params.cgbo * vgmb;
                        qgmid = qgdo + qgso + qgmb;
                        qbulk -= qgmb;
                        qsrc = -(qgate + qgmid + qbulk + qdrn);
                    } else {
                        gcggb = (newop.cggb + cgdo + cgso + self.size_params.cgbo) * ag0;
                        gcgdb = (newop.cgdb - cgdo) * ag0;
                        gcgsb = (newop.cgsb - cgso) * ag0;
                        gcgbb = -(gcggb + gcgdb + gcgsb);

                        gcdgb = (newop.cdgb - cgdo) * ag0;
                        gcsgb = -(newop.cggb + newop.cbgb + newop.cdgb + cgso) * ag0;
                        gcbgb = (newop.cbgb - self.size_params.cgbo) * ag0;

                        gcdgmb = 0.0;
                        gcsgmb = 0.0;
                        gcbgmb = 0.0;

                        qgb = self.size_params.cgbo * vgb;
                        qgate += qgdo + qgso + qgb;
                        qbulk -= qgb;
                        qsrc = -(qgate + qbulk + qdrn);
                    }
                    gcddb = (newop.cddb + newop.capbd + cgdo) * ag0;
                    gcdsb = newop.cdsb * ag0;

                    gcsdb = -(newop.cgdb + newop.cbdb + newop.cddb) * ag0;
                    gcssb = (newop.capbs + cgso - (newop.cgsb + newop.cbsb + newop.cdsb)) * ag0;

                    if self.intp.rbodyMod == 0 {
                        gcdbb = -(gcdgb + gcddb + gcdsb + gcdgmb);
                        gcsbb = -(gcsgb + gcsdb + gcssb + gcsgmb);
                        gcbdb = (newop.cbdb - newop.capbd) * ag0;
                        gcbsb = (newop.cbsb - newop.capbs) * ag0;
                        gcdbdb = 0.0;
                        gcsbsb = 0.0;
                    } else {
                        gcdbb = -(newop.cddb + newop.cdgb + newop.cdsb) * ag0;
                        gcsbb = -(gcsgb + gcsdb + gcssb + gcsgmb) + newop.capbs * ag0;
                        gcbdb = newop.cbdb * ag0;
                        gcbsb = newop.cbsb * ag0;

                        gcdbdb = -newop.capbd * ag0;
                        gcsbsb = -newop.capbs * ag0;
                    }
                    gcbbb = -(gcbdb + gcbgb + gcbsb + gcbgmb);

                    ggtg = 0.0;
                    ggtd = 0.0;
                    ggtb = 0.0;
                    ggts = 0.0;
                    sxpart = 0.6;
                    dxpart = 0.4;
                    ddxpart_dVd = 0.0;
                    ddxpart_dVg = 0.0;
                    ddxpart_dVb = 0.0;
                    ddxpart_dVs = 0.0;
                    dsxpart_dVd = 0.0;
                    dsxpart_dVg = 0.0;
                    dsxpart_dVb = 0.0;
                    dsxpart_dVs = 0.0;
                } else {
                    qcheq = newop.qchqs;
                    CoxWL = self.intp.coxe
                        * self.size_params.weffCV
                        * self.intp.nf
                        * self.size_params.leffCV;
                    T0 = qdef * ScalingFactor / CoxWL;

                    ggtg = T0 * newop.gcrgg;
                    newop.gtg = ggtg;
                    ggtd = T0 * newop.gcrgd;
                    newop.gtd = ggtd;
                    ggts = T0 * newop.gcrgs;
                    newop.gts = ggts;
                    ggtb = T0 * newop.gcrgb;
                    newop.gtb = ggtb;
                    gqdef = ScalingFactor * ag0;

                    gcqgb = newop.cqgb * ag0;
                    gcqdb = newop.cqdb * ag0;
                    gcqsb = newop.cqsb * ag0;
                    gcqbb = newop.cqbb * ag0;

                    if qcheq.abs() <= 1.0e-5 * CoxWL {
                        if self.model.xpart < 0.5 {
                            dxpart = 0.4;
                        } else if self.model.xpart > 0.5 {
                            dxpart = 0.0;
                        } else {
                            dxpart = 0.5;
                        }
                        ddxpart_dVd = 0.0;
                        ddxpart_dVg = 0.0;
                        ddxpart_dVb = 0.0;
                        ddxpart_dVs = 0.0;
                    } else {
                        dxpart = qdrn / qcheq;
                        Cdd = newop.cddb;
                        Csd = -(newop.cgdb + newop.cddb + newop.cbdb);
                        ddxpart_dVd = (Cdd - dxpart * (Cdd + Csd)) / qcheq;
                        Cdg = newop.cdgb;
                        Csg = -(newop.cggb + newop.cdgb + newop.cbgb);
                        ddxpart_dVg = (Cdg - dxpart * (Cdg + Csg)) / qcheq;

                        Cds = newop.cdsb;
                        Css = -(newop.cgsb + newop.cdsb + newop.cbsb);
                        ddxpart_dVs = (Cds - dxpart * (Cds + Css)) / qcheq;

                        ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);
                    }
                    sxpart = 1.0 - dxpart;
                    dsxpart_dVd = -ddxpart_dVd;
                    dsxpart_dVg = -ddxpart_dVg;
                    dsxpart_dVs = -ddxpart_dVs;
                    dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);

                    if self.intp.rgateMod == 3 {
                        gcgmgmb = (cgdo + cgso + self.size_params.cgbo) * ag0;
                        gcgmdb = -cgdo * ag0;
                        gcgmsb = -cgso * ag0;
                        gcgmbb = -self.size_params.cgbo * ag0;

                        gcdgmb = gcgmdb;
                        gcsgmb = gcgmsb;
                        gcbgmb = gcgmbb;

                        gcdgb = 0.0;
                        gcsgb = 0.0;
                        gcbgb = 0.0;
                        gcggb = 0.0;
                        gcgdb = 0.0;
                        gcgsb = 0.0;
                        gcgbb = 0.0;

                        qgmb = self.size_params.cgbo * vgmb;
                        qgmid = qgdo + qgso + qgmb;
                        qgate = 0.0;
                        qbulk = -qgmb;
                        qdrn = -qgdo;
                        qsrc = -(qgmid + qbulk + qdrn);
                    } else {
                        gcggb = (cgdo + cgso + self.size_params.cgbo) * ag0;
                        gcgdb = -cgdo * ag0;
                        gcgsb = -cgso * ag0;
                        gcgbb = -self.size_params.cgbo * ag0;

                        gcdgb = gcgdb;
                        gcsgb = gcgsb;
                        gcbgb = gcgbb;
                        gcdgmb = 0.0;
                        gcsgmb = 0.0;
                        gcbgmb = 0.0;

                        qgb = self.size_params.cgbo * vgb;
                        qgate = qgdo + qgso + qgb;
                        qbulk = -qgb;
                        qdrn = -qgdo;
                        qsrc = -(qgate + qbulk + qdrn);
                    }

                    gcddb = (newop.capbd + cgdo) * ag0;
                    gcdsb = 0.0;
                    gcsdb = 0.0;
                    gcssb = (newop.capbs + cgso) * ag0;

                    if self.intp.rbodyMod == 0 {
                        gcdbb = -(gcdgb + gcddb + gcdgmb);
                        gcsbb = -(gcsgb + gcssb + gcsgmb);
                        gcbdb = -newop.capbd * ag0;
                        gcbsb = -newop.capbs * ag0;
                        gcdbdb = 0.0;
                        gcsbsb = 0.0;
                    } else {
                        gcdbb = 0.0;
                        gcsbb = 0.0;
                        gcbdb = 0.0;
                        gcbsb = 0.0;
                        gcdbdb = -newop.capbd * ag0;
                        gcsbsb = -newop.capbs * ag0;
                    }
                    gcbbb = -(gcbdb + gcbgb + gcbsb + gcbgmb);
                }
            } else {
                if self.model.trnqsmod == 0 {
                    qsrc = qdrn - qgso;
                    if self.intp.rgateMod == 3 {
                        gcgmgmb = (cgdo + cgso + self.size_params.cgbo) * ag0;
                        gcgmdb = -cgdo * ag0;
                        gcgmsb = -cgso * ag0;
                        gcgmbb = -self.size_params.cgbo * ag0;

                        gcdgmb = gcgmdb;
                        gcsgmb = gcgmsb;
                        gcbgmb = gcgmbb;

                        gcggb = newop.cggb * ag0;
                        gcgdb = newop.cgsb * ag0;
                        gcgsb = newop.cgdb * ag0;
                        gcgbb = -(gcggb + gcgdb + gcgsb);

                        gcdgb = -(newop.cggb + newop.cbgb + newop.cdgb) * ag0;
                        gcsgb = newop.cdgb * ag0;
                        gcbgb = newop.cbgb * ag0;

                        qgmb = self.size_params.cgbo * vgmb;
                        qgmid = qgdo + qgso + qgmb;
                        qbulk -= qgmb;
                        qdrn = -(qgate + qgmid + qbulk + qsrc);
                    } else {
                        gcggb = (newop.cggb + cgdo + cgso + self.size_params.cgbo) * ag0;
                        gcgdb = (newop.cgsb - cgdo) * ag0;
                        gcgsb = (newop.cgdb - cgso) * ag0;
                        gcgbb = -(gcggb + gcgdb + gcgsb);

                        gcdgb = -(newop.cggb + newop.cbgb + newop.cdgb + cgdo) * ag0;
                        gcsgb = (newop.cdgb - cgso) * ag0;
                        gcbgb = (newop.cbgb - self.size_params.cgbo) * ag0;

                        gcdgmb = 0.0;
                        gcsgmb = 0.0;
                        gcbgmb = 0.0;

                        qgb = self.size_params.cgbo * vgb;
                        qgate += qgdo + qgso + qgb;
                        qbulk -= qgb;
                        qdrn = -(qgate + qbulk + qsrc);
                    }
                    gcddb = (newop.capbd + cgdo - (newop.cgsb + newop.cbsb + newop.cdsb)) * ag0;
                    gcdsb = -(newop.cgdb + newop.cbdb + newop.cddb) * ag0;

                    gcsdb = newop.cdsb * ag0;
                    gcssb = (newop.cddb + newop.capbs + cgso) * ag0;

                    if self.intp.rbodyMod == 0 {
                        gcdbb = -(gcdgb + gcddb + gcdsb + gcdgmb);
                        gcsbb = -(gcsgb + gcsdb + gcssb + gcsgmb);
                        gcbdb = (newop.cbsb - newop.capbd) * ag0;
                        gcbsb = (newop.cbdb - newop.capbs) * ag0;
                        gcdbdb = 0.0;
                        gcsbsb = 0.0;
                    } else {
                        gcdbb = -(gcdgb + gcddb + gcdsb + gcdgmb) + newop.capbd * ag0;
                        gcsbb = -(newop.cddb + newop.cdgb + newop.cdsb) * ag0;
                        gcbdb = newop.cbsb * ag0;
                        gcbsb = newop.cbdb * ag0;
                        gcdbdb = -newop.capbd * ag0;
                        gcsbsb = -newop.capbs * ag0;
                    }
                    gcbbb = -(gcbgb + gcbdb + gcbsb + gcbgmb);

                    ggtg = 0.0;
                    ggtd = 0.0;
                    ggtb = 0.0;
                    ggts = 0.0;
                    sxpart = 0.4;
                    dxpart = 0.6;
                    ddxpart_dVd = 0.0;
                    ddxpart_dVg = 0.0;
                    ddxpart_dVb = 0.0;
                    ddxpart_dVs = 0.0;
                    dsxpart_dVd = 0.0;
                    dsxpart_dVg = 0.0;
                    dsxpart_dVb = 0.0;
                    dsxpart_dVs = 0.0;
                } else {
                    qcheq = newop.qchqs;
                    CoxWL = self.intp.coxe
                        * self.size_params.weffCV
                        * self.intp.nf
                        * self.size_params.leffCV;
                    T0 = qdef * ScalingFactor / CoxWL;
                    ggtg = T0 * newop.gcrgg;
                    newop.gtg = ggtg;
                    ggts = T0 * newop.gcrgd;
                    newop.gts = ggts;
                    ggtd = T0 * newop.gcrgs;
                    newop.gtd = ggtd;
                    ggtb = T0 * newop.gcrgb;
                    newop.gtb = ggtb;
                    gqdef = ScalingFactor * ag0;

                    gcqgb = newop.cqgb * ag0;
                    gcqdb = newop.cqsb * ag0;
                    gcqsb = newop.cqdb * ag0;
                    gcqbb = newop.cqbb * ag0;

                    if qcheq.abs() <= 1.0e-5 * CoxWL {
                        if self.model.xpart < 0.5 {
                            sxpart = 0.4;
                        } else if self.model.xpart > 0.5 {
                            sxpart = 0.0;
                        } else {
                            sxpart = 0.5;
                        }
                        dsxpart_dVd = 0.0;
                        dsxpart_dVg = 0.0;
                        dsxpart_dVb = 0.0;
                        dsxpart_dVs = 0.0;
                    } else {
                        sxpart = qdrn / qcheq;
                        Css = newop.cddb;
                        Cds = -(newop.cgdb + newop.cddb + newop.cbdb);
                        dsxpart_dVs = (Css - sxpart * (Css + Cds)) / qcheq;
                        Csg = newop.cdgb;
                        Cdg = -(newop.cggb + newop.cdgb + newop.cbgb);
                        dsxpart_dVg = (Csg - sxpart * (Csg + Cdg)) / qcheq;

                        Csd = newop.cdsb;
                        Cdd = -(newop.cgsb + newop.cdsb + newop.cbs);
                        dsxpart_dVd = (Csd - sxpart * (Csd + Cdd)) / qcheq;

                        dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);
                    }
                    dxpart = 1.0 - sxpart;
                    ddxpart_dVd = -dsxpart_dVd;
                    ddxpart_dVg = -dsxpart_dVg;
                    ddxpart_dVs = -dsxpart_dVs;
                    ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);

                    if self.intp.rgateMod == 3 {
                        gcgmgmb = (cgdo + cgso + self.size_params.cgbo) * ag0;
                        gcgmdb = -cgdo * ag0;
                        gcgmsb = -cgso * ag0;
                        gcgmbb = -self.size_params.cgbo * ag0;

                        gcdgmb = gcgmdb;
                        gcsgmb = gcgmsb;
                        gcbgmb = gcgmbb;

                        gcdgb = 0.0;
                        gcsgb = 0.0;
                        gcbgb = 0.0;
                        gcggb = 0.0;
                        gcgdb = 0.0;
                        gcgsb = 0.0;
                        gcgbb = 0.0;

                        qgmb = self.size_params.cgbo * vgmb;
                        qgmid = qgdo + qgso + qgmb;
                        qgate = 0.0;
                        qbulk = -qgmb;
                        qdrn = -qgdo;
                        qsrc = -qgso;
                    } else {
                        gcggb = (cgdo + cgso + self.size_params.cgbo) * ag0;
                        gcgdb = -cgdo * ag0;
                        gcgsb = -cgso * ag0;
                        gcgbb = -self.size_params.cgbo * ag0;

                        gcdgb = gcgdb;
                        gcsgb = gcgsb;
                        gcbgb = gcgbb;
                        gcdgmb = 0.0;
                        gcsgmb = 0.0;
                        gcbgmb = 0.0;

                        qgb = self.size_params.cgbo * vgb;
                        qgate = qgdo + qgso + qgb;
                        qbulk = -qgb;
                        qdrn = -qgdo;
                        qsrc = -qgso;
                    }

                    gcddb = (newop.capbd + cgdo) * ag0;
                    gcdsb = 0.0;
                    gcsdb = 0.0;
                    gcssb = (newop.capbs + cgso) * ag0;
                    if self.intp.rbodyMod == 0 {
                        gcdbb = -(gcdgb + gcddb + gcdgmb);
                        gcsbb = -(gcsgb + gcssb + gcsgmb);
                        gcbdb = -newop.capbd * ag0;
                        gcbsb = -newop.capbs * ag0;
                        gcdbdb = 0.0;
                        gcsbsb = 0.0;
                    } else {
                        gcdbb = 0.0;
                        gcsbb = 0.0;
                        gcbdb = 0.0;
                        gcbsb = 0.0;
                        gcdbdb = -newop.capbd * ag0;
                        gcsbsb = -newop.capbs * ag0;
                    }
                    gcbbb = -(gcbdb + gcbgb + gcbsb + gcbgmb);
                }
            }
            let mut _tmp1: f64;
            let mut _tmp2: f64;
            if self.model.trnqsmod != 0 {
                newop.qcdump = qdef * ScalingFactor;
                let (_g, i, _r) =
                    state.integrate(newop.qcdump - self.op.qcdump, 0.0, 0.0, self.op.cqcdump);
                newop.cqcdump = i;
            }

            newop.qg = qgate;
            newop.qd = qdrn - newop.qbd;
            newop.qs = qsrc - newop.qbs;
            if self.intp.rgateMod == 3 {
                newop.qgmid = qgmid;
            }

            if self.intp.rbodyMod != 0 {
                newop.qb = qbulk + newop.qbd + newop.qbs;
            } else {
                newop.qb = qbulk;
            }

            // error = NIintegrate(ckt, &geq, &ceq, 0.0, newop.qb);
            let (_g, i, _r) = state.integrate(newop.qb - self.op.qb, 0.0, 0.0, self.op.cqb);
            newop.cqb = i;
            // error = NIintegrate(ckt, &geq, &ceq, 0.0, newop.qg);
            let (_g, i, _r) = state.integrate(newop.qg - self.op.qg, 0.0, 0.0, self.op.cqg);
            newop.cqg = i;
            // error = NIintegrate(ckt, &geq, &ceq, 0.0, newop.qd);
            let (_g, i, _r) = state.integrate(newop.qd - self.op.qd, 0.0, 0.0, self.op.cqd);
            newop.cqd = i;

            if self.intp.rgateMod == 3 {
                // error = NIintegrate(ckt, &geq, &ceq, 0.0, newop.qgmid);
                let (_g, i, _r) =
                    state.integrate(newop.qgmid - self.op.qgmid, 0.0, 0.0, self.op.cqgmid);
                newop.cqgmid = i;
            }

            if self.intp.rbodyMod != 0 {
                // error = NIintegrate(ckt, &geq, &ceq, 0.0, newop.qbs);
                let (_g, i, _r) = state.integrate(newop.qbs - self.op.qbs, 0.0, 0.0, self.op.cqbs);
                newop.cqbs = i;
                // error = NIintegrate(ckt, &geq, &ceq, 0.0, newop.qbd);
                let (_g, i, _r) = state.integrate(newop.qbd - self.op.qbd, 0.0, 0.0, self.op.cqbd);
                newop.cqbd = i;
            }

            /* Calculate equivalent charge current */

            cqgate = newop.cqg;
            cqbody = newop.cqb;
            cqdrn = newop.cqd;

            ceqqg = cqgate - gcggb * vgb + gcgdb * vbd + gcgsb * vbs;
            ceqqd = cqdrn - gcdgb * vgb - gcdgmb * vgmb + (gcddb + gcdbdb) * vbd - gcdbdb * vbd_jct
                + gcdsb * vbs;
            ceqqb = cqbody - gcbgb * vgb - gcbgmb * vgmb + gcbdb * vbd + gcbsb * vbs;

            if self.intp.rgateMod == 3 {
                ceqqgmid = newop.cqgmid + gcgmdb * vbd + gcgmsb * vbs - gcgmgmb * vgmb;
            } else {
                ceqqgmid = 0.0;
            }

            if self.intp.rbodyMod != 0 {
                ceqqjs = newop.cqbs + gcsbsb * vbs_jct;
                ceqqjd = newop.cqbd + gcdbdb * vbd_jct;
            }

            if self.model.trnqsmod != 0 {
                T0 = ggtg * vgb - ggtd * vbd - ggts * vbs;
                ceqqg += T0;
                T1 = qdef * newop.gtau;
                ceqqd -=
                    dxpart * T0 + T1 * (ddxpart_dVg * vgb - ddxpart_dVd * vbd - ddxpart_dVs * vbs);
                cqdef = newop.cqcdump - gqdef * qdef;
                cqcheq = newop.cqcheq - (gcqgb * vgb - gcqdb * vbd - gcqsb * vbs) + T0;
            }
        }

        /*
         *  Load current vector
         */

        if self.guess.mode >= 0 {
            Gm = newop.gm;
            Gmbs = newop.gmbs;
            FwdSum = Gm + Gmbs;
            RevSum = 0.0;

            ceqdrn = self.model.p() * (cdrain - newop.gds * vds - Gm * vgs - Gmbs * vbs);
            ceqbd = self.model.p()
                * (newop.csub + newop.Igidl
                    - (newop.gbds + newop.ggidld) * vds
                    - (newop.gbgs + newop.ggidlg) * vgs
                    - (newop.gbbs + newop.ggidlb) * vbs);
            ceqbs = self.model.p()
                * (newop.Igisl + newop.ggisls * vds - newop.ggislg * vgd - newop.ggislb * vbd);

            gbbdp = -(newop.gbds);
            gbbsp = newop.gbds + newop.gbgs + newop.gbbs;

            gbdpg = newop.gbgs;
            gbdpdp = newop.gbds;
            gbdpb = newop.gbbs;
            gbdpsp = -(gbdpg + gbdpdp + gbdpb);

            gbspg = 0.0;
            gbspdp = 0.0;
            gbspb = 0.0;
            gbspsp = 0.0;

            if self.model.igcmod != 0 {
                gIstotg = newop.gIgsg + newop.gIgcsg;
                gIstotd = newop.gIgcsd;
                gIstots = newop.gIgss + newop.gIgcss;
                gIstotb = newop.gIgcsb;
                Istoteq = self.model.p()
                    * (newop.Igs + newop.Igcs
                        - gIstotg * vgs
                        - newop.gIgcsd * vds
                        - newop.gIgcsb * vbs);

                gIdtotg = newop.gIgdg + newop.gIgcdg;
                gIdtotd = newop.gIgdd + newop.gIgcdd;
                gIdtots = newop.gIgcds;
                gIdtotb = newop.gIgcdb;
                Idtoteq = self.model.p()
                    * (newop.Igd + newop.Igcd
                        - newop.gIgdg * vgd
                        - newop.gIgcdg * vgs
                        - newop.gIgcdd * vds
                        - newop.gIgcdb * vbs);
            } else {
                gIstotg = 0.0;
                gIstotd = 0.0;
                gIstots = 0.0;
                gIstotb = 0.0;
                Istoteq = 0.0;
                gIdtotg = 0.0;
                gIdtotd = 0.0;
                gIdtots = 0.0;
                gIdtotb = 0.0;
                Idtoteq = 0.0;
            }

            if self.model.igbmod != 0 {
                gIbtotg = newop.gIgbg;
                gIbtotd = newop.gIgbd;
                gIbtots = newop.gIgbs;
                gIbtotb = newop.gIgbb;
                Ibtoteq = self.model.p()
                    * (newop.Igb - newop.gIgbg * vgs - newop.gIgbd * vds - newop.gIgbb * vbs);
            } else {
                gIbtotg = 0.0;
                gIbtotd = 0.0;
                gIbtots = 0.0;
                gIbtotb = 0.0;
                Ibtoteq = 0.0;
            }

            if (self.model.igcmod != 0) || (self.model.igbmod != 0) {
                gIgtotg = gIstotg + gIdtotg + gIbtotg;
                gIgtotd = gIstotd + gIdtotd + gIbtotd;
                gIgtots = gIstots + gIdtots + gIbtots;
                gIgtotb = gIstotb + gIdtotb + gIbtotb;
                Igtoteq = Istoteq + Idtoteq + Ibtoteq;
            } else {
                gIgtotg = 0.0;
                gIgtotd = 0.0;
                gIgtots = 0.0;
                gIgtotb = 0.0;
                Igtoteq = 0.0;
            }

            if (self.intp.rgateMod == 2) {
                T0 = vges - vgs;
            } else if (self.intp.rgateMod == 3) {
                T0 = vgms - vgs;
            }
            if self.intp.rgateMod > 1 {
                gcrgd = newop.gcrgd * T0;
                gcrgg = newop.gcrgg * T0;
                gcrgs = newop.gcrgs * T0;
                gcrgb = newop.gcrgb * T0;
                ceqgcrg = -(gcrgd * vds + gcrgg * vgs + gcrgb * vbs);
                gcrgg -= newop.gcrg;
                gcrg = newop.gcrg;
            } else {
                ceqgcrg = 0.0;
                gcrg = 0.0;
                gcrgd = 0.0;
                gcrgg = 0.0;
                gcrgs = 0.0;
                gcrgb = 0.0;
            }
        } else {
            Gm = -newop.gm;
            Gmbs = -newop.gmbs;
            FwdSum = 0.0;
            RevSum = -(Gm + Gmbs);

            ceqdrn = -self.model.p() * (cdrain + newop.gds * vds + Gm * vgd + Gmbs * vbd);

            ceqbs = self.model.p()
                * (newop.csub + newop.Igisl + (newop.gbds + newop.ggisls) * vds
                    - (newop.gbgs + newop.ggislg) * vgd
                    - (newop.gbbs + newop.ggislb) * vbd);
            ceqbd = self.model.p()
                * (newop.Igidl - newop.ggidld * vds - newop.ggidlg * vgs - newop.ggidlb * vbs);

            gbbsp = -(newop.gbds);
            gbbdp = newop.gbds + newop.gbgs + newop.gbbs;

            gbdpg = 0.0;
            gbdpsp = 0.0;
            gbdpb = 0.0;
            gbdpdp = 0.0;

            gbspg = newop.gbgs;
            gbspsp = newop.gbds;
            gbspb = newop.gbbs;
            gbspdp = -(gbspg + gbspsp + gbspb);

            if self.model.igcmod != 0 {
                gIstotg = newop.gIgsg + newop.gIgcdg;
                gIstotd = newop.gIgcds;
                gIstots = newop.gIgss + newop.gIgcdd;
                gIstotb = newop.gIgcdb;
                Istoteq = self.model.p()
                    * (newop.Igs + newop.Igcd - newop.gIgsg * vgs - newop.gIgcdg * vgd
                        + newop.gIgcdd * vds
                        - newop.gIgcdb * vbd);

                gIdtotg = newop.gIgdg + newop.gIgcsg;
                gIdtotd = newop.gIgdd + newop.gIgcss;
                gIdtots = newop.gIgcsd;
                gIdtotb = newop.gIgcsb;
                Idtoteq = self.model.p()
                    * (newop.Igd + newop.Igcs - (newop.gIgdg + newop.gIgcsg) * vgd
                        + newop.gIgcsd * vds
                        - newop.gIgcsb * vbd);
            } else {
                gIstotg = 0.0;
                gIstotd = 0.0;
                gIstots = 0.0;
                gIstotb = 0.0;
                Istoteq = 0.0;
                gIdtotg = 0.0;
                gIdtotd = 0.0;
                gIdtots = 0.0;
                gIdtotb = 0.0;
                Idtoteq = 0.0;
            }

            if self.model.igbmod != 0 {
                gIbtotg = newop.gIgbg;
                gIbtotd = newop.gIgbs;
                gIbtots = newop.gIgbd;
                gIbtotb = newop.gIgbb;
                Ibtoteq = self.model.p()
                    * (newop.Igb - newop.gIgbg * vgd + newop.gIgbd * vds - newop.gIgbb * vbd);
            } else {
                gIbtotg = 0.0;
                gIbtotd = 0.0;
                gIbtots = 0.0;
                gIbtotb = 0.0;
                Ibtoteq = 0.0;
            }

            if (self.model.igcmod != 0) || (self.model.igbmod != 0) {
                gIgtotg = gIstotg + gIdtotg + gIbtotg;
                gIgtotd = gIstotd + gIdtotd + gIbtotd;
                gIgtots = gIstots + gIdtots + gIbtots;
                gIgtotb = gIstotb + gIdtotb + gIbtotb;
                Igtoteq = Istoteq + Idtoteq + Ibtoteq;
            } else {
                gIgtotg = 0.0;
                gIgtotd = 0.0;
                gIgtots = 0.0;
                gIgtotb = 0.0;
                Igtoteq = 0.0;
            }

            if (self.intp.rgateMod == 2) {
                T0 = vges - vgs;
            } else if (self.intp.rgateMod == 3) {
                T0 = vgms - vgs;
            }
            if self.intp.rgateMod > 1 {
                gcrgd = newop.gcrgs * T0;
                gcrgg = newop.gcrgg * T0;
                gcrgs = newop.gcrgd * T0;
                gcrgb = newop.gcrgb * T0;
                ceqgcrg = -(gcrgg * vgd - gcrgs * vds + gcrgb * vbd);
                gcrgg -= newop.gcrg;
                gcrg = newop.gcrg;
            } else {
                ceqgcrg = 0.0;
                gcrg = 0.0;
                gcrgd = 0.0;
                gcrgg = 0.0;
                gcrgs = 0.0;
                gcrgb = 0.0;
            }
        }

        if self.model.rdsmod == 1 {
            ceqgstot =
                self.model.p() * (newop.gstotd * vds + newop.gstotg * vgs + newop.gstotb * vbs);
            gstot = newop.gstot;
            gstotd = newop.gstotd;
            gstotg = newop.gstotg;
            gstots = newop.gstots - gstot;
            gstotb = newop.gstotb;

            ceqgdtot =
                -self.model.p() * (newop.gdtotd * vds + newop.gdtotg * vgs + newop.gdtotb * vbs);
            gdtot = newop.gdtot;
            gdtotd = newop.gdtotd - gdtot;
            gdtotg = newop.gdtotg;
            gdtots = newop.gdtots;
            gdtotb = newop.gdtotb;
        } else {
            gstot = 0.0;
            gstotd = 0.0;
            gstotg = 0.0;
            gstots = 0.0;
            gstotb = 0.0;
            ceqgstot = 0.0;
            gdtot = 0.0;
            gdtotd = 0.0;
            gdtotg = 0.0;
            gdtots = 0.0;
            gdtotb = 0.0;
            ceqgdtot = 0.0;
        }

        if self.model.p() > 0.0 {
            ceqjs = (newop.cbs - newop.gbs * vbs_jct);
            ceqjd = (newop.cbd - newop.gbd * vbd_jct);
        } else {
            ceqjs = -(newop.cbs - newop.gbs * vbs_jct);
            ceqjd = -(newop.cbd - newop.gbd * vbd_jct);
            ceqqg = -ceqqg;
            ceqqd = -ceqqd;
            ceqqb = -ceqqb;
            ceqgcrg = -ceqgcrg;

            if self.model.trnqsmod != 0 {
                cqdef = -cqdef;
                cqcheq = -cqcheq;
            }

            if self.intp.rbodyMod != 0 {
                ceqqjs = -ceqqjs;
                ceqqjd = -ceqqjd;
            }

            if self.intp.rgateMod == 3 {
                ceqqgmid = -ceqqgmid;
            }
        }

        // Gather up RHS current-vector terms

        let mut b: Vec<(Option<VarIndex>, f64)> = vec![];

        b.push((
            self.ports.dNodePrime,
            (ceqjd - ceqbd + ceqgdtot - ceqdrn - ceqqd + Idtoteq),
        ));
        b.push((self.ports.gNodePrime, -(ceqqg - ceqgcrg + Igtoteq)));

        if (self.intp.rgateMod == 2) {
            b.push((self.ports.gNodeExt, -ceqgcrg));
        } else if (self.intp.rgateMod == 3) {
            b.push((self.ports.gNodeMid, -(ceqqgmid + ceqgcrg)));
        }

        if self.intp.rbodyMod == 0 {
            b.push((
                self.ports.bNodePrime,
                (ceqbd + ceqbs - ceqjd - ceqjs - ceqqb + Ibtoteq),
            ));
            b.push((
                self.ports.sNodePrime,
                (ceqdrn - ceqbs + ceqjs + ceqqg + ceqqb + ceqqd + ceqqgmid - ceqgstot + Istoteq),
            ));
        } else {
            b.push((self.ports.dbNode, -(ceqjd + ceqqjd)));
            b.push((self.ports.bNodePrime, (ceqbd + ceqbs - ceqqb + Ibtoteq)));
            b.push((self.ports.sbNode, -(ceqjs + ceqqjs)));
            b.push((
                self.ports.sNodePrime,
                (ceqdrn - ceqbs + ceqjs + ceqqd + ceqqg + ceqqb + ceqqjd + ceqqjs + ceqqgmid
                    - ceqgstot
                    + Istoteq),
            ));
        }

        if self.model.rdsmod != 0 {
            b.push((self.ports.dNode, -ceqgdtot));
            b.push((self.ports.sNode, ceqgstot));
        }
        if self.model.trnqsmod != 0 {
            b.push((self.ports.qNode, cqcheq - cqdef));
        }

        // Gather up matrix Jacobian terms
        let mut j: Vec<(Option<Eindex>, f64)> = vec![];

        if self.intp.rbodyMod != 0 {
            gjbd = newop.gbd;
            gjbs = newop.gbs;
        } else {
            gjbd = 0.0;
            gjbs = 0.0;
        }

        if self.model.rdsmod != 0 {
            gdpr = self.intp.drainConductance;
            gspr = self.intp.sourceConductance;
        } else {
            gdpr = 0.0;
            gspr = 0.0;
        }

        geltd = self.intp.grgeltd;
        T1 = qdef * newop.gtau;

        if self.intp.rgateMod == 1 {
            j.push((self.matps.GEgePtr, geltd));
            j.push((self.matps.GPgePtr, -(geltd)));
            j.push((self.matps.GEgpPtr, -(geltd)));
            j.push((self.matps.GPgpPtr, gcggb + geltd - ggtg + gIgtotg));
            j.push((self.matps.GPdpPtr, gcgdb - ggtd + gIgtotd));
            j.push((self.matps.GPspPtr, gcgsb - ggts + gIgtots));
            j.push((self.matps.GPbpPtr, gcgbb - ggtb + gIgtotb));
        } else if self.intp.rgateMod == 2 {
            j.push((self.matps.GEgePtr, gcrg));
            j.push((self.matps.GEgpPtr, gcrgg));
            j.push((self.matps.GEdpPtr, gcrgd));
            j.push((self.matps.GEspPtr, gcrgs));
            j.push((self.matps.GEbpPtr, gcrgb));

            j.push((self.matps.GPgePtr, -(gcrg)));
            j.push((self.matps.GPgpPtr, gcggb - gcrgg - ggtg + gIgtotg));
            j.push((self.matps.GPdpPtr, gcgdb - gcrgd - ggtd + gIgtotd));
            j.push((self.matps.GPspPtr, gcgsb - gcrgs - ggts + gIgtots));
            j.push((self.matps.GPbpPtr, gcgbb - gcrgb - ggtb + gIgtotb));
        } else if self.intp.rgateMod == 3 {
            j.push((self.matps.GEgePtr, geltd));
            j.push((self.matps.GEgmPtr, -(geltd)));
            j.push((self.matps.GMgePtr, -(geltd)));
            j.push((self.matps.GMgmPtr, geltd + gcrg + gcgmgmb));

            j.push((self.matps.GMdpPtr, gcrgd + gcgmdb));
            j.push((self.matps.GMgpPtr, gcrgg));
            j.push((self.matps.GMspPtr, gcrgs + gcgmsb));
            j.push((self.matps.GMbpPtr, gcrgb + gcgmbb));

            j.push((self.matps.DPgmPtr, gcdgmb));
            j.push((self.matps.GPgmPtr, -(gcrg)));
            j.push((self.matps.SPgmPtr, gcsgmb));
            j.push((self.matps.BPgmPtr, gcbgmb));

            j.push((self.matps.GPgpPtr, gcggb - gcrgg - ggtg + gIgtotg));
            j.push((self.matps.GPdpPtr, gcgdb - gcrgd - ggtd + gIgtotd));
            j.push((self.matps.GPspPtr, gcgsb - gcrgs - ggts + gIgtots));
            j.push((self.matps.GPbpPtr, gcgbb - gcrgb - ggtb + gIgtotb));
        } else {
            j.push((self.matps.GPgpPtr, gcggb - ggtg + gIgtotg));
            j.push((self.matps.GPdpPtr, gcgdb - ggtd + gIgtotd));
            j.push((self.matps.GPspPtr, gcgsb - ggts + gIgtots));
            j.push((self.matps.GPbpPtr, gcgbb - ggtb + gIgtotb));
        }

        if self.model.rdsmod != 0 {
            j.push((self.matps.DgpPtr, gdtotg));
            j.push((self.matps.DspPtr, gdtots));
            j.push((self.matps.DbpPtr, gdtotb));
            j.push((self.matps.SdpPtr, gstotd));
            j.push((self.matps.SgpPtr, gstotg));
            j.push((self.matps.SbpPtr, gstotb));
        }

        j.push((
            self.matps.DPdpPtr,
            gdpr + newop.gds + newop.gbd + T1 * ddxpart_dVd - gdtotd
                + RevSum
                + gcddb
                + gbdpdp
                + dxpart * ggtd
                - gIdtotd,
        ));
        j.push((self.matps.DPdPtr, -(gdpr + gdtot)));
        j.push((
            self.matps.DPgpPtr,
            Gm + gcdgb - gdtotg + gbdpg - gIdtotg + dxpart * ggtg + T1 * ddxpart_dVg,
        ));
        j.push((
            self.matps.DPspPtr,
            -(newop.gds + gdtots - dxpart * ggts + gIdtots - T1 * ddxpart_dVs + FwdSum
                - gcdsb
                - gbdpsp),
        ));
        j.push((
            self.matps.DPbpPtr,
            -(gjbd + gdtotb - Gmbs - gcdbb - gbdpb + gIdtotb - T1 * ddxpart_dVb - dxpart * ggtb),
        ));

        j.push((self.matps.DdpPtr, -(gdpr - gdtotd)));
        j.push((self.matps.DdPtr, gdpr + gdtot));

        j.push((
            self.matps.SPdpPtr,
            -(newop.gds + gstotd + RevSum - gcsdb - gbspdp - T1 * dsxpart_dVd - sxpart * ggtd
                + gIstotd),
        ));
        j.push((
            self.matps.SPgpPtr,
            gcsgb - Gm - gstotg + gbspg + sxpart * ggtg + T1 * dsxpart_dVg - gIstotg,
        ));
        j.push((
            self.matps.SPspPtr,
            gspr + newop.gds + newop.gbs + T1 * dsxpart_dVs - gstots
                + FwdSum
                + gcssb
                + gbspsp
                + sxpart * ggts
                - gIstots,
        ));
        j.push((self.matps.SPsPtr, -(gspr + gstot)));
        j.push((
            self.matps.SPbpPtr,
            -(gjbs + gstotb + Gmbs - gcsbb - gbspb - sxpart * ggtb - T1 * dsxpart_dVb + gIstotb),
        ));

        j.push((self.matps.SspPtr, -(gspr - gstots)));
        j.push((self.matps.SsPtr, gspr + gstot));

        j.push((self.matps.BPdpPtr, gcbdb - gjbd + gbbdp - gIbtotd));
        j.push((self.matps.BPgpPtr, gcbgb - newop.gbgs - gIbtotg));
        j.push((self.matps.BPspPtr, gcbsb - gjbs + gbbsp - gIbtots));
        j.push((
            self.matps.BPbpPtr,
            gjbd + gjbs + gcbbb - newop.gbbs - gIbtotb,
        ));

        ggidld = newop.ggidld;
        ggidlg = newop.ggidlg;
        ggidlb = newop.ggidlb;
        ggislg = newop.ggislg;
        ggisls = newop.ggisls;
        ggislb = newop.ggislb;

        /* stamp gidl */
        j.push((self.matps.DPdpPtr, ggidld));
        j.push((self.matps.DPgpPtr, ggidlg));
        j.push((self.matps.DPspPtr, -(ggidlg + ggidld + ggidlb)));
        j.push((self.matps.DPbpPtr, ggidlb));
        j.push((self.matps.BPdpPtr, -(ggidld)));
        j.push((self.matps.BPgpPtr, -(ggidlg)));
        j.push((self.matps.BPspPtr, (ggidlg + ggidld + ggidlb)));
        j.push((self.matps.BPbpPtr, -(ggidlb)));
        /* stamp gisl */
        j.push((self.matps.SPdpPtr, -(ggisls + ggislg + ggislb)));
        j.push((self.matps.SPgpPtr, ggislg));
        j.push((self.matps.SPspPtr, ggisls));
        j.push((self.matps.SPbpPtr, ggislb));
        j.push((self.matps.BPdpPtr, (ggislg + ggisls + ggislb)));
        j.push((self.matps.BPgpPtr, -(ggislg)));
        j.push((self.matps.BPspPtr, -(ggisls)));
        j.push((self.matps.BPbpPtr, -(ggislb)));

        if self.intp.rbodyMod != 0 {
            j.push((self.matps.DPdbPtr, gcdbdb - newop.gbd));
            j.push((self.matps.SPsbPtr, -(newop.gbs - gcsbsb)));

            j.push((self.matps.DBdpPtr, gcdbdb - newop.gbd));
            j.push((
                self.matps.DBdbPtr,
                newop.gbd - gcdbdb + newop.grbpd + newop.grbdb,
            ));
            j.push((self.matps.DBbpPtr, -(newop.grbpd)));
            j.push((self.matps.DBbPtr, -(newop.grbdb)));

            j.push((self.matps.BPdbPtr, -(newop.grbpd)));
            j.push((self.matps.BPbPtr, -(newop.grbpb)));
            j.push((self.matps.BPsbPtr, -(newop.grbps)));
            j.push((self.matps.BPbpPtr, newop.grbpd + newop.grbps + newop.grbpb));

            j.push((self.matps.SBspPtr, gcsbsb - newop.gbs));
            j.push((self.matps.SBbpPtr, -(newop.grbps)));
            j.push((self.matps.SBbPtr, -(newop.grbsb)));
            j.push((
                self.matps.SBsbPtr,
                newop.gbs - gcsbsb + newop.grbps + newop.grbsb,
            ));

            j.push((self.matps.BdbPtr, -(newop.grbdb)));
            j.push((self.matps.BbpPtr, -(newop.grbpb)));
            j.push((self.matps.BsbPtr, -(newop.grbsb)));
            j.push((self.matps.BbPtr, newop.grbsb + newop.grbdb + newop.grbpb));
        }

        if self.model.trnqsmod != 0 {
            j.push((self.matps.QqPtr, gqdef + newop.gtau));
            j.push((self.matps.QgpPtr, ggtg - gcqgb));
            j.push((self.matps.QdpPtr, ggtd - gcqdb));
            j.push((self.matps.QspPtr, ggts - gcqsb));
            j.push((self.matps.QbpPtr, ggtb - gcqbb));

            j.push((self.matps.DPqPtr, dxpart * newop.gtau));
            j.push((self.matps.SPqPtr, sxpart * newop.gtau));
            j.push((self.matps.GPqPtr, -(newop.gtau)));
        }
        // Update our best-guess operating point
        self.guess = newop;
        // And return our matrix stamps
        return Stamps { g: j, b };
    }
}

/// compute poly depletion effect
fn polyDepletion(phi: f64, ngate: f64, epsgate: f64, coxe: f64, Vgs: f64) -> (f64, f64) {
    if (ngate > 1.0e18) && (ngate < 1.0e25) && (Vgs > phi) && (epsgate != 0.0) {
        let T1 = 1.0e6 * Q * epsgate * ngate / (coxe * coxe);
        let T8 = Vgs - phi;
        let T4 = sqrt(1.0 + 2.0 * T8 / T1);
        let T2 = 2.0 * T8 / (T4 + 1.0);
        let T3 = 0.5 * T2 * T2 / T1; /* T3 = Vpoly */
        let T7 = 1.12 - T3 - 0.05;
        let T6 = sqrt(T7 * T7 + 0.224);
        let T5 = 1.12 - 0.5 * (T7 + T6);
        let Vgs_eff = Vgs - T5;
        let dVgs_eff_dVg = 1.0 - (0.5 - 0.5 / T4) * (1.0 + T7 / T6);
        (Vgs_eff, dVgs_eff_dVg)
    } else {
        (Vgs, 1.0)
    }
}

// FIXME - implement these here limiting methods!
fn DEVlimvds(_a: f64, _b: f64) -> f64 {
    _a
}
fn DEVfetlim(a: f64, b: f64, c: f64) -> f64 {
    a
}
fn DEVpnjlim(a: f64, b: f64, c: f64, d: f64) -> f64 {
    a
}

// Some helper math
fn sqrt(a: f64) -> f64 {
    a.sqrt()
}
fn log(a: f64) -> f64 {
    a.ln()
}
fn exp(a: f64) -> f64 {
    a.exp()
}
