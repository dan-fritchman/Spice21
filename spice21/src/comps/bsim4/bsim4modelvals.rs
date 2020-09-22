use super::bsim4defs::*;
use super::log;
use crate::comps::consts::*;
use crate::comps::MosType;

/// Resolve input-provided model-specs into their values, 
/// incorporating defaults and limiting constraints. 
fn resolve(specs: &Bsim4ModelSpecs) -> Bsim4ModelVals {
    let mut vals = Bsim4ModelVals::default();
    use MosType::{NMOS, PMOS};

    // vals.type = if let Some(val) = specs.type { if val > 0 { NMOS } else { PMOS} } else { NMOS };
    vals.mos_type = NMOS; // FIXME!
    let tnom = 300.15; // FIXME: from ckt->CKTnomTemp

    vals.mobmod = if let Some(val) = specs.mobmod { val } else { 0 };
    if vals.mobmod > 6 {
        vals.mobmod = 0;
        println!("Warning: mobmod has been set to its default value: 0.\n");
    }
    vals.diomod = if let Some(val) = specs.diomod { val } else { 1 };
    if vals.diomod > 2 {
        vals.diomod = 1;
        println!("Warning: diomod has been set to its default value: 1.\n");
    }
    vals.capmod = if let Some(val) = specs.capmod { val } else { 2 };
    if vals.capmod > 2 {
        vals.capmod = 2;
        println!("Warning: capmod has been set to its default value: 2.\n");
    }
    vals.rdsmod = if let Some(val) = specs.rdsmod { val } else { 0 };
    if vals.rdsmod > 1 {
        vals.rdsmod = 0;
        println!("Warning: rdsmod has been set to its default value: 0.\n");
    }
    vals.rbodymod = if let Some(val) = specs.rbodymod { val } else { 0 };
    if vals.rbodymod > 2 {
        vals.rbodymod = 0;
        println!("Warning: rbodymod has been set to its default value: 0.\n");
    }
    vals.rgatemod = if let Some(val) = specs.rgatemod { val } else { 0 };
    if vals.rgatemod > 3 {
        vals.rgatemod = 0;
        println!("Warning: rgatemod has been set to its default value: 0.\n");
    }
    vals.permod = if let Some(val) = specs.permod { val } else { 1 };
    if vals.permod > 1 {
        vals.permod = 1;
        println!("Warning: permod has been set to its default value: 1.\n");
    }
    vals.fnoimod = if let Some(val) = specs.fnoimod { val } else { 1 };
    if vals.fnoimod > 1 {
        vals.fnoimod = 1;
        println!("Warning: fnoimod has been set to its default value: 1.\n");
    }
    vals.tnoimod = if let Some(val) = specs.tnoimod { val } else { 0 };
    if vals.tnoimod > 2 {
        vals.tnoimod = 0;
        println!("Warning: tnoimod has been set to its default value: 0.\n");
    }
    vals.trnqsmod = if let Some(val) = specs.trnqsmod { val } else { 0 };
    if vals.trnqsmod > 1 {
        vals.trnqsmod = 0;
        println!("Warning: trnqsmod has been set to its default value: 0.\n");
    }
    vals.acnqsmod = if let Some(val) = specs.acnqsmod { val } else { 0 };
    if vals.acnqsmod > 1 {
        vals.acnqsmod = 0;
        println!("Warning: acnqsmod has been set to its default value: 0.\n");
    }
    vals.mtrlmod = if let Some(val) = specs.mtrlmod { val } else { 0 };
    if vals.mtrlmod > 1 {
        vals.mtrlmod = 0;
        println!("Warning: mtrlmod has been set to its default value: 0.\n");
    }
    vals.mtrlcompatmod = if let Some(val) = specs.mtrlcompatmod { val } else { 0 };
    if vals.mtrlcompatmod > 1 {
        vals.mtrlcompatmod = 0;
        println!("Warning: mtrlcompatmod has been set to its default value: 0.\n");
    }
    vals.igcmod = if let Some(val) = specs.igcmod { val } else { 0 };
    if vals.igcmod > 2 {
        vals.igcmod = 0;
        println!("Warning: igcmod has been set to its default value: 0.\n");
    }
    vals.igbmod = if let Some(val) = specs.igbmod { val } else { 0 };
    if vals.igbmod > 1 {
        vals.igbmod = 0;
        println!("Warning: igbmod has been set to its default value: 0.\n");
    }
    vals.tempmod = if let Some(val) = specs.tempmod { val } else { 0 };
    if vals.tempmod > 3 {
        vals.tempmod = 0;
        println!("Warning: tempmod has been set to its default value: 0.\n");
    }
    vals.wpemod = if let Some(val) = specs.wpemod { val } else { 0 };
    if vals.wpemod > 1 {
        vals.wpemod = 0;
        println!("Warning: wpemod has been set to its default value: 0.\n");
    }

    // FIXME: range check these
    vals.gidlmod = if let Some(val) = specs.gidlmod { val } else { 0 };
    vals.geomod = if let Some(val) = specs.geomod { val } else { 0 };
    vals.cvchargemod = if let Some(val) = specs.cvchargemod { val } else { 0 };
    vals.binunit = if let Some(val) = specs.binunit { val } else { 1 };
    vals.paramchk = if let Some(val) = specs.paramchk { val } else { 1 };

    // Beginning primary double-valued params
    vals.version = if let Some(val) = specs.version { val } else { 4.80 };
    vals.toxref = if let Some(val) = specs.toxref { val } else { 30.0e-10 };
    vals.eot = if let Some(val) = specs.eot { val } else { 15.0e-10 };
    vals.vddeot = if let Some(val) = specs.vddeot {
        val
    } else {
        match vals.mos_type {
            NMOS => 1.5,
            PMOS => -1.5,
        }
    };
    vals.tempeot = if let Some(val) = specs.tempeot { val } else { 300.15 };
    vals.leffeot = if let Some(val) = specs.leffeot { val } else { 1.0 };
    vals.weffeot = if let Some(val) = specs.weffeot { val } else { 10.0 };
    vals.ados = if let Some(val) = specs.ados { val } else { 1.0 };
    vals.bdos = if let Some(val) = specs.bdos { val } else { 1.0 };
    vals.toxe = if let Some(val) = specs.toxe { val } else { 30.0e-10 };
    vals.toxp = if let Some(val) = specs.toxp { val } else { vals.toxe };
    vals.toxm = if let Some(val) = specs.toxm { val } else { vals.toxe };
    vals.dtox = if let Some(val) = specs.dtox { val } else { 0.0 };
    vals.epsrox = if let Some(val) = specs.epsrox { val } else { 3.9 };
    vals.cdsc = if let Some(val) = specs.cdsc { val } else { 2.4e-4 };
    vals.cdscb = if let Some(val) = specs.cdscb { val } else { 0.0 };
    vals.cdscd = if let Some(val) = specs.cdscd { val } else { 0.0 };
    vals.cit = if let Some(val) = specs.cit { val } else { 0.0 };
    vals.nfactor = if let Some(val) = specs.nfactor { val } else { 1.0 };
    vals.xj = if let Some(val) = specs.xj { val } else { 0.15e-6 };
    vals.vsat = if let Some(val) = specs.vsat { val } else { 8.0e4 };
    vals.at = if let Some(val) = specs.at { val } else { 3.3e4 };
    vals.a0 = if let Some(val) = specs.a0 { val } else { 1.0 };
    vals.ags = if let Some(val) = specs.ags { val } else { 0.0 };
    vals.a1 = if let Some(val) = specs.a1 { val } else { 0.0 };
    vals.a2 = if let Some(val) = specs.a2 { val } else { 1.0 };
    vals.keta = if let Some(val) = specs.keta { val } else { -0.047 };
    vals.nsub = if let Some(val) = specs.nsub { val } else { 6.0e16 };
    vals.phig = if let Some(val) = specs.phig { val } else { 4.05 };
    vals.epsrgate = if let Some(val) = specs.epsrgate { val } else { 11.7 };
    vals.easub = if let Some(val) = specs.easub { val } else { 4.05 };
    vals.epsrsub = if let Some(val) = specs.epsrsub { val } else { 11.7 };
    vals.ni0sub = if let Some(val) = specs.ni0sub { val } else { 1.45e10 };
    vals.bg0sub = if let Some(val) = specs.bg0sub { val } else { 1.16 };
    vals.tbgasub = if let Some(val) = specs.tbgasub { val } else { 7.02e-4 };
    vals.tbgbsub = if let Some(val) = specs.tbgbsub { val } else { 1108.0 };
    vals.ndep = if let Some(val) = specs.ndep { val } else { 1.7e17 };
    vals.nsd = if let Some(val) = specs.nsd { val } else { 1.0e20 };
    vals.phin = if let Some(val) = specs.phin { val } else { 0.0 };
    vals.ngate = if let Some(val) = specs.ngate { val } else { 0.0 };
    vals.vbm = if let Some(val) = specs.vbm { val } else { -3.0 };
    vals.xt = if let Some(val) = specs.xt { val } else { 1.55e-7 };
    vals.kt1 = if let Some(val) = specs.kt1 { val } else { -0.11 };
    vals.kt1l = if let Some(val) = specs.kt1l { val } else { 0.0 };
    vals.kt2 = if let Some(val) = specs.kt2 { val } else { 0.022 };
    vals.k3 = if let Some(val) = specs.k3 { val } else { 80.0 };
    vals.k3b = if let Some(val) = specs.k3b { val } else { 0.0 };
    vals.w0 = if let Some(val) = specs.w0 { val } else { 2.5e-6 };
    vals.lpe0 = if let Some(val) = specs.lpe0 { val } else { 1.74e-7 };
    vals.lpeb = if let Some(val) = specs.lpeb { val } else { 0.0 };
    vals.dvtp0 = if let Some(val) = specs.dvtp0 { val } else { 0.0 };
    vals.dvtp1 = if let Some(val) = specs.dvtp1 { val } else { 0.0 };
    vals.dvtp2 = if let Some(val) = specs.dvtp2 { val } else { 0.0 };
    vals.dvtp3 = if let Some(val) = specs.dvtp3 { val } else { 0.0 };
    vals.dvtp4 = if let Some(val) = specs.dvtp4 { val } else { 0.0 };
    vals.dvtp5 = if let Some(val) = specs.dvtp5 { val } else { 0.0 };
    vals.dvt0 = if let Some(val) = specs.dvt0 { val } else { 2.2 };
    vals.dvt1 = if let Some(val) = specs.dvt1 { val } else { 0.53 };
    vals.dvt2 = if let Some(val) = specs.dvt2 { val } else { -0.032 };

    vals.dvt0w = if let Some(val) = specs.dvt0w { val } else { 0.0 };
    vals.dvt1w = if let Some(val) = specs.dvt1w { val } else { 5.3e6 };
    vals.dvt2w = if let Some(val) = specs.dvt2w { val } else { -0.032 };

    vals.drout = if let Some(val) = specs.drout { val } else { 0.56 };
    vals.dsub = if let Some(val) = specs.dsub { val } else { vals.drout };
    vals.vth0 = if let Some(val) = specs.vth0 {
        val
    } else {
        match vals.mos_type {
            NMOS => 0.7,
            PMOS => -0.7,
        }
    };
    vals.vfb = if let Some(val) = specs.vfb { val } else { -1.0 };
    vals.eu = if let Some(val) = specs.eu {
        val
    } else {
        match vals.mos_type {
            NMOS => 1.67,
            PMOS => 1.0,
        }
    };
    vals.ucs = if let Some(val) = specs.ucs {
        val
    } else {
        match vals.mos_type {
            NMOS => 1.67,
            PMOS => 1.0,
        }
    };
    vals.ua = if let Some(val) = specs.ua {
        val
    } else if vals.mobmod == 2 {
        1.0e-15
    } else {
        1.0e-9
    };
    vals.ua1 = if let Some(val) = specs.ua1 { val } else { 1.0e-9 };
    vals.ub = if let Some(val) = specs.ub { val } else { 1.0e-19 };
    vals.ub1 = if let Some(val) = specs.ub1 { val } else { -1.0e-18 };
    vals.uc = if let Some(val) = specs.uc {
        val
    } else if vals.mobmod == 1 {
        -0.0465
    } else {
        -0.0465e-9
    };
    vals.uc1 = if let Some(val) = specs.uc1 {
        val
    } else if vals.mobmod == 1 {
        -0.056
    } else {
        -0.056e-9
    };
    vals.ud = if let Some(val) = specs.ud { val } else { 0.0 };
    vals.ud1 = if let Some(val) = specs.ud1 { val } else { 0.0 };
    vals.up = if let Some(val) = specs.up { val } else { 0.0 };
    vals.lp = if let Some(val) = specs.lp { val } else { 1.0e-8 };
    vals.u0 = if let Some(val) = specs.u0 {
        val
    } else {
        match vals.mos_type {
            NMOS => 0.067,
            PMOS => 0.025,
        }
    };
    vals.ute = if let Some(val) = specs.ute { val } else { -1.5 };
    vals.ucste = if let Some(val) = specs.ucste { val } else { -4.775e-3 };
    vals.voff = if let Some(val) = specs.voff { val } else { -0.08 };
    vals.voffl = if let Some(val) = specs.voffl { val } else { 0.0 };
    vals.voffcvl = if let Some(val) = specs.voffcvl { val } else { 0.0 };
    vals.minv = if let Some(val) = specs.minv { val } else { 0.0 };
    vals.minvcv = if let Some(val) = specs.minvcv { val } else { 0.0 };
    vals.fprout = if let Some(val) = specs.fprout { val } else { 0.0 };
    vals.pdits = if let Some(val) = specs.pdits { val } else { 0.0 };
    vals.pditsd = if let Some(val) = specs.pditsd { val } else { 0.0 };
    vals.pditsl = if let Some(val) = specs.pditsl { val } else { 0.0 };
    vals.delta = if let Some(val) = specs.delta { val } else { 0.01 };
    vals.rdswmin = if let Some(val) = specs.rdswmin { val } else { 0.0 };
    vals.rdwmin = if let Some(val) = specs.rdwmin { val } else { 0.0 };
    vals.rswmin = if let Some(val) = specs.rswmin { val } else { 0.0 };
    vals.rdsw = if let Some(val) = specs.rdsw { val } else { 200.0 }; /* in ohm*um */
    vals.rdw = if let Some(val) = specs.rdw { val } else { 100.0 };
    vals.rsw = if let Some(val) = specs.rsw { val } else { 100.0 };
    vals.prwg = if let Some(val) = specs.prwg { val } else { 1.0 }; /* in 1/V */
    vals.prwb = if let Some(val) = specs.prwb { val } else { 0.0 };
    vals.prt = if let Some(val) = specs.prt { val } else { 0.0 };
    vals.eta0 = if let Some(val) = specs.eta0 { val } else { 0.08 }; /* no unit  */
    vals.etab = if let Some(val) = specs.etab { val } else { -0.07 };
    vals.pclm = if let Some(val) = specs.pclm { val } else { 1.3 }; /* no unit  */
    vals.pdiblc1 = if let Some(val) = specs.pdiblc1 { val } else { 0.39 }; /* no unit  */
    vals.pdiblc2 = if let Some(val) = specs.pdiblc2 { val } else { 0.0086 }; /* no unit  */
    vals.pdiblcb = if let Some(val) = specs.pdiblcb { val } else { 0.0 }; /* 1/V  */
    vals.pscbe1 = if let Some(val) = specs.pscbe1 { val } else { 4.24e8 };
    vals.pscbe2 = if let Some(val) = specs.pscbe2 { val } else { 1.0e-5 };
    vals.pvag = if let Some(val) = specs.pvag { val } else { 0.0 };
    vals.wr = if let Some(val) = specs.wr { val } else { 1.0 };
    vals.dwg = if let Some(val) = specs.dwg { val } else { 0.0 };
    vals.dwb = if let Some(val) = specs.dwb { val } else { 0.0 };
    vals.b0 = if let Some(val) = specs.b0 { val } else { 0.0 };
    vals.b1 = if let Some(val) = specs.b1 { val } else { 0.0 };
    vals.alpha0 = if let Some(val) = specs.alpha0 { val } else { 0.0 };
    vals.alpha1 = if let Some(val) = specs.alpha1 { val } else { 0.0 };
    vals.beta0 = if let Some(val) = specs.beta0 { val } else { 0.0 };
    vals.agidl = if let Some(val) = specs.agidl { val } else { 0.0 };
    vals.bgidl = if let Some(val) = specs.bgidl { val } else { 2.3e9 }; /* V/m */
    vals.cgidl = if let Some(val) = specs.cgidl { val } else { 0.5 }; /* V^3 */
    vals.egidl = if let Some(val) = specs.egidl { val } else { 0.8 }; /* V */
    vals.rgidl = if let Some(val) = specs.rgidl { val } else { 1.0 };
    vals.kgidl = if let Some(val) = specs.kgidl { val } else { 0.0 };
    vals.fgidl = if let Some(val) = specs.fgidl { val } else { 1.0 };

    /*Default value of agisl, bgisl, cgisl, egisl, rgisl, kgisl, and fgisl are set as follows */
    vals.agisl = if let Some(val) = specs.agisl { val } else { vals.agidl };
    vals.bgisl = if let Some(val) = specs.bgisl { val } else { vals.bgidl };
    vals.cgisl = if let Some(val) = specs.cgisl { val } else { vals.cgidl };
    vals.egisl = if let Some(val) = specs.egisl { val } else { vals.egidl };
    vals.rgisl = if let Some(val) = specs.rgisl { val } else { vals.rgidl };
    vals.kgisl = if let Some(val) = specs.kgisl { val } else { vals.kgidl };
    vals.fgisl = if let Some(val) = specs.fgisl { val } else { vals.fgidl };

    vals.aigc = if let Some(val) = specs.aigc {
        val
    } else {
        match vals.mos_type {
            NMOS => 1.36e-2,
            PMOS => 9.80e-3,
        }
    };
    vals.bigc = if let Some(val) = specs.bigc {
        val
    } else {
        match vals.mos_type {
            NMOS => 1.71e-3,
            PMOS => 7.59e-4,
        }
    };
    vals.cigc = if let Some(val) = specs.cigc {
        val
    } else {
        match vals.mos_type {
            NMOS => 0.075,
            PMOS => 0.03,
        }
    };
    if let Some(aigsd) = specs.aigsd {
        vals.aigs = aigsd;
        vals.aigd = aigsd;
    } else {
        vals.aigsd = match vals.mos_type {
            NMOS => 1.36e-2,
            PMOS => 9.80e-3,
        };
        vals.aigs = if let Some(val) = specs.aigs {
            val
        } else {
            match vals.mos_type {
                NMOS => 1.36e-2,
                PMOS => 9.80e-3,
            }
        };
        vals.aigd = if let Some(val) = specs.aigd {
            val
        } else {
            match vals.mos_type {
                NMOS => 1.36e-2,
                PMOS => 9.80e-3,
            }
        };
    }
    if let Some(bigsd) = specs.bigsd {
        vals.bigs = bigsd;
        vals.bigd = bigsd;
    } else {
        vals.bigsd = match vals.mos_type {
            NMOS => 1.71e-3,
            PMOS => 7.59e-4,
        };
        vals.bigs = if let Some(val) = specs.bigs {
            val
        } else {
            match vals.mos_type {
                NMOS => 1.71e-3,
                PMOS => 7.59e-4,
            }
        };
        vals.bigd = if let Some(val) = specs.bigd {
            val
        } else {
            match vals.mos_type {
                NMOS => 1.71e-3,
                PMOS => 7.59e-4,
            }
        };
    }
    if let Some(cigsd) = specs.cigsd {
        vals.cigs = cigsd;
        vals.cigd = cigsd;
    } else {
        vals.cigsd = match vals.mos_type {
            NMOS => 0.075,
            PMOS => 0.03,
        };
        vals.cigs = if let Some(val) = specs.cigs {
            val
        } else {
            match vals.mos_type {
                NMOS => 0.075,
                PMOS => 0.03,
            }
        };
        vals.cigd = if let Some(val) = specs.cigd {
            val
        } else {
            match vals.mos_type {
                NMOS => 0.075,
                PMOS => 0.03,
            }
        };
    }
    vals.aigbacc = if let Some(val) = specs.aigbacc { val } else { 1.36e-2 };
    vals.bigbacc = if let Some(val) = specs.bigbacc { val } else { 1.71e-3 };
    vals.cigbacc = if let Some(val) = specs.cigbacc { val } else { 0.075 };
    vals.aigbinv = if let Some(val) = specs.aigbinv { val } else { 1.11e-2 };
    vals.bigbinv = if let Some(val) = specs.bigbinv { val } else { 9.49e-4 };
    vals.cigbinv = if let Some(val) = specs.cigbinv { val } else { 0.006 };
    vals.nigc = if let Some(val) = specs.nigc { val } else { 1.0 };
    vals.nigbinv = if let Some(val) = specs.nigbinv { val } else { 3.0 };
    vals.nigbacc = if let Some(val) = specs.nigbacc { val } else { 1.0 };
    vals.ntox = if let Some(val) = specs.ntox { val } else { 1.0 };
    vals.eigbinv = if let Some(val) = specs.eigbinv { val } else { 1.1 };
    vals.pigcd = if let Some(val) = specs.pigcd { val } else { 1.0 };
    vals.poxedge = if let Some(val) = specs.poxedge { val } else { 1.0 };
    vals.xrcrg1 = if let Some(val) = specs.xrcrg1 { val } else { 12.0 };
    vals.xrcrg2 = if let Some(val) = specs.xrcrg2 { val } else { 1.0 };
    vals.ijthsfwd = if let Some(val) = specs.ijthsfwd { val } else { 0.1 };
    vals.ijthdfwd = if let Some(val) = specs.ijthdfwd { val } else { vals.ijthsfwd };
    vals.ijthsrev = if let Some(val) = specs.ijthsrev { val } else { 0.1 };
    vals.ijthdrev = if let Some(val) = specs.ijthdrev { val } else { vals.ijthsrev };
    vals.tnoia = if let Some(val) = specs.tnoia { val } else { 1.5 };
    vals.tnoib = if let Some(val) = specs.tnoib { val } else { 3.5 };
    vals.tnoic = if let Some(val) = specs.tnoic { val } else { 0.0 };
    vals.rnoia = if let Some(val) = specs.rnoia { val } else { 0.577 };
    vals.rnoib = if let Some(val) = specs.rnoib { val } else { 0.5164 };
    vals.rnoic = if let Some(val) = specs.rnoic { val } else { 0.395 };
    vals.ntnoi = if let Some(val) = specs.ntnoi { val } else { 1.0 };
    vals.lambda = if let Some(val) = specs.lambda { val } else { 0.0 };
    vals.vtl = if let Some(val) = specs.vtl { val } else { 2.0e5 };
    vals.xn = if let Some(val) = specs.xn { val } else { 3.0 };
    vals.lc = if let Some(val) = specs.lc { val } else { 5.0e-9 };
    vals.vfbsdoff = if let Some(val) = specs.vfbsdoff { val } else { 0.0 };
    vals.tvfbsdoff = if let Some(val) = specs.tvfbsdoff { val } else { 0.0 };
    vals.tvoff = if let Some(val) = specs.tvoff { val } else { 0.0 };
    vals.tnfactor = if let Some(val) = specs.tnfactor { val } else { 0.0 };
    vals.teta0 = if let Some(val) = specs.teta0 { val } else { 0.0 };
    vals.tvoffcv = if let Some(val) = specs.tvoffcv { val } else { 0.0 };

    vals.lintnoi = if let Some(val) = specs.lintnoi { val } else { 0.0 };

    vals.xjbvs = if let Some(val) = specs.xjbvs { val } else { 1.0 }; /* no unit */
    vals.xjbvd = if let Some(val) = specs.xjbvd { val } else { vals.xjbvs };
    vals.bvs = if let Some(val) = specs.bvs { val } else { 10.0 }; /* V */
    vals.bvd = if let Some(val) = specs.bvd { val } else { vals.bvs };

    vals.gbmin = if let Some(val) = specs.gbmin { val } else { 1.0e-12 }; /* in mho */
    vals.rbdb = if let Some(val) = specs.rbdb { val } else { 50.0 }; /* in ohm */
    vals.rbpb = if let Some(val) = specs.rbpb { val } else { 50.0 };
    vals.rbsb = if let Some(val) = specs.rbsb { val } else { 50.0 };
    vals.rbps = if let Some(val) = specs.rbps { val } else { 50.0 };
    vals.rbpd = if let Some(val) = specs.rbpd { val } else { 50.0 };

    vals.rbps0 = if let Some(val) = specs.rbps0 { val } else { 50.0 };
    vals.rbpsl = if let Some(val) = specs.rbpsl { val } else { 0.0 };
    vals.rbpsw = if let Some(val) = specs.rbpsw { val } else { 0.0 };
    vals.rbpsnf = if let Some(val) = specs.rbpsnf { val } else { 0.0 };

    vals.rbpd0 = if let Some(val) = specs.rbpd0 { val } else { 50.0 };
    vals.rbpdl = if let Some(val) = specs.rbpdl { val } else { 0.0 };
    vals.rbpdw = if let Some(val) = specs.rbpdw { val } else { 0.0 };
    vals.rbpdnf = if let Some(val) = specs.rbpdnf { val } else { 0.0 };

    vals.rbpbx0 = if let Some(val) = specs.rbpbx0 { val } else { 100.0 };
    vals.rbpbxl = if let Some(val) = specs.rbpbxl { val } else { 0.0 };
    vals.rbpbxw = if let Some(val) = specs.rbpbxw { val } else { 0.0 };
    vals.rbpbxnf = if let Some(val) = specs.rbpbxnf { val } else { 0.0 };
    vals.rbpby0 = if let Some(val) = specs.rbpby0 { val } else { 100.0 };
    vals.rbpbyl = if let Some(val) = specs.rbpbyl { val } else { 0.0 };
    vals.rbpbyw = if let Some(val) = specs.rbpbyw { val } else { 0.0 };
    vals.rbpbynf = if let Some(val) = specs.rbpbynf { val } else { 0.0 };

    vals.rbsbx0 = if let Some(val) = specs.rbsbx0 { val } else { 100.0 };
    vals.rbsby0 = if let Some(val) = specs.rbsby0 { val } else { 100.0 };
    vals.rbdbx0 = if let Some(val) = specs.rbdbx0 { val } else { 100.0 };
    vals.rbdby0 = if let Some(val) = specs.rbdby0 { val } else { 100.0 };

    vals.rbsdbxl = if let Some(val) = specs.rbsdbxl { val } else { 0.0 };
    vals.rbsdbxw = if let Some(val) = specs.rbsdbxw { val } else { 0.0 };
    vals.rbsdbxnf = if let Some(val) = specs.rbsdbxnf { val } else { 0.0 };
    vals.rbsdbyl = if let Some(val) = specs.rbsdbyl { val } else { 0.0 };
    vals.rbsdbyw = if let Some(val) = specs.rbsdbyw { val } else { 0.0 };
    vals.rbsdbynf = if let Some(val) = specs.rbsdbynf { val } else { 0.0 };

    vals.cgsl = if let Some(val) = specs.cgsl { val } else { 0.0 };
    vals.cgdl = if let Some(val) = specs.cgdl { val } else { 0.0 };
    vals.ckappas = if let Some(val) = specs.ckappas { val } else { 0.6 };
    vals.ckappad = if let Some(val) = specs.ckappad { val } else { vals.ckappas };
    vals.clc = if let Some(val) = specs.clc { val } else { 0.1e-6 };
    vals.cle = if let Some(val) = specs.cle { val } else { 0.6 };
    vals.vfbcv = if let Some(val) = specs.vfbcv { val } else { -1.0 };
    vals.acde = if let Some(val) = specs.acde { val } else { 1.0 };
    vals.moin = if let Some(val) = specs.moin { val } else { 15.0 };
    vals.noff = if let Some(val) = specs.noff { val } else { 1.0 };
    vals.voffcv = if let Some(val) = specs.voffcv { val } else { 0.0 };
    vals.dmcg = if let Some(val) = specs.dmcg { val } else { 0.0 };
    vals.dmci = if let Some(val) = specs.dmci { val } else { vals.dmcg };
    vals.dmdg = if let Some(val) = specs.dmdg { val } else { 0.0 };
    vals.dmcgt = if let Some(val) = specs.dmcgt { val } else { 0.0 };
    vals.xgw = if let Some(val) = specs.xgw { val } else { 0.0 };
    vals.xgl = if let Some(val) = specs.xgl { val } else { 0.0 };
    vals.rshg = if let Some(val) = specs.rshg { val } else { 0.1 };
    vals.ngcon = if let Some(val) = specs.ngcon { val } else { 1.0 };
    vals.tcj = if let Some(val) = specs.tcj { val } else { 0.0 };
    vals.tpb = if let Some(val) = specs.tpb { val } else { 0.0 };
    vals.tcjsw = if let Some(val) = specs.tcjsw { val } else { 0.0 };
    vals.tpbsw = if let Some(val) = specs.tpbsw { val } else { 0.0 };
    vals.tcjswg = if let Some(val) = specs.tcjswg { val } else { 0.0 };
    vals.tpbswg = if let Some(val) = specs.tpbswg { val } else { 0.0 };

    /* Length dependence */
    vals.lcdsc = if let Some(val) = specs.lcdsc { val } else { 0.0 };
    vals.lcdscb = if let Some(val) = specs.lcdscb { val } else { 0.0 };
    vals.lcdscd = if let Some(val) = specs.lcdscd { val } else { 0.0 };
    vals.lcit = if let Some(val) = specs.lcit { val } else { 0.0 };
    vals.lnfactor = if let Some(val) = specs.lnfactor { val } else { 0.0 };
    vals.lxj = if let Some(val) = specs.lxj { val } else { 0.0 };
    vals.lvsat = if let Some(val) = specs.lvsat { val } else { 0.0 };
    vals.lat = if let Some(val) = specs.lat { val } else { 0.0 };
    vals.la0 = if let Some(val) = specs.la0 { val } else { 0.0 };
    vals.lags = if let Some(val) = specs.lags { val } else { 0.0 };
    vals.la1 = if let Some(val) = specs.la1 { val } else { 0.0 };
    vals.la2 = if let Some(val) = specs.la2 { val } else { 0.0 };
    vals.lketa = if let Some(val) = specs.lketa { val } else { 0.0 };
    vals.lnsub = if let Some(val) = specs.lnsub { val } else { 0.0 };
    vals.lndep = if let Some(val) = specs.lndep { val } else { 0.0 };
    vals.lnsd = if let Some(val) = specs.lnsd { val } else { 0.0 };
    vals.lphin = if let Some(val) = specs.lphin { val } else { 0.0 };
    vals.lngate = if let Some(val) = specs.lngate { val } else { 0.0 };
    vals.lvbm = if let Some(val) = specs.lvbm { val } else { 0.0 };
    vals.lxt = if let Some(val) = specs.lxt { val } else { 0.0 };
    vals.lk1 = if let Some(val) = specs.lk1 { val } else { 0.0 };
    vals.lkt1 = if let Some(val) = specs.lkt1 { val } else { 0.0 };
    vals.lkt1l = if let Some(val) = specs.lkt1l { val } else { 0.0 };
    vals.lkt2 = if let Some(val) = specs.lkt2 { val } else { 0.0 };
    vals.lk2 = if let Some(val) = specs.lk2 { val } else { 0.0 };
    vals.lk3 = if let Some(val) = specs.lk3 { val } else { 0.0 };
    vals.lk3b = if let Some(val) = specs.lk3b { val } else { 0.0 };
    vals.lw0 = if let Some(val) = specs.lw0 { val } else { 0.0 };
    vals.llpe0 = if let Some(val) = specs.llpe0 { val } else { 0.0 };
    vals.llpeb = if let Some(val) = specs.llpeb { val } else { 0.0 };
    vals.ldvtp0 = if let Some(val) = specs.ldvtp0 { val } else { 0.0 };
    vals.ldvtp1 = if let Some(val) = specs.ldvtp1 { val } else { 0.0 };
    vals.ldvtp2 = if let Some(val) = specs.ldvtp2 { val } else { 0.0 };
    vals.ldvtp3 = if let Some(val) = specs.ldvtp3 { val } else { 0.0 };
    vals.ldvtp4 = if let Some(val) = specs.ldvtp4 { val } else { 0.0 };
    vals.ldvtp5 = if let Some(val) = specs.ldvtp5 { val } else { 0.0 };
    vals.ldvt0 = if let Some(val) = specs.ldvt0 { val } else { 0.0 };
    vals.ldvt1 = if let Some(val) = specs.ldvt1 { val } else { 0.0 };
    vals.ldvt2 = if let Some(val) = specs.ldvt2 { val } else { 0.0 };
    vals.ldvt0w = if let Some(val) = specs.ldvt0w { val } else { 0.0 };
    vals.ldvt1w = if let Some(val) = specs.ldvt1w { val } else { 0.0 };
    vals.ldvt2w = if let Some(val) = specs.ldvt2w { val } else { 0.0 };
    vals.ldrout = if let Some(val) = specs.ldrout { val } else { 0.0 };
    vals.ldsub = if let Some(val) = specs.ldsub { val } else { 0.0 };
    vals.lvth0 = if let Some(val) = specs.lvth0 { val } else { 0.0 };
    vals.lua = if let Some(val) = specs.lua { val } else { 0.0 };
    vals.lua1 = if let Some(val) = specs.lua1 { val } else { 0.0 };
    vals.lub = if let Some(val) = specs.lub { val } else { 0.0 };
    vals.lub1 = if let Some(val) = specs.lub1 { val } else { 0.0 };
    vals.luc = if let Some(val) = specs.luc { val } else { 0.0 };
    vals.luc1 = if let Some(val) = specs.luc1 { val } else { 0.0 };
    vals.lud = if let Some(val) = specs.lud { val } else { 0.0 };
    vals.lud1 = if let Some(val) = specs.lud1 { val } else { 0.0 };
    vals.lup = if let Some(val) = specs.lup { val } else { 0.0 };
    vals.llp = if let Some(val) = specs.llp { val } else { 0.0 };
    vals.lu0 = if let Some(val) = specs.lu0 { val } else { 0.0 };
    vals.lute = if let Some(val) = specs.lute { val } else { 0.0 };
    vals.lucste = if let Some(val) = specs.lucste { val } else { 0.0 };
    vals.lvoff = if let Some(val) = specs.lvoff { val } else { 0.0 };
    vals.lminv = if let Some(val) = specs.lminv { val } else { 0.0 };
    vals.lminvcv = if let Some(val) = specs.lminvcv { val } else { 0.0 };
    vals.lfprout = if let Some(val) = specs.lfprout { val } else { 0.0 };
    vals.lpdits = if let Some(val) = specs.lpdits { val } else { 0.0 };
    vals.lpditsd = if let Some(val) = specs.lpditsd { val } else { 0.0 };
    vals.ldelta = if let Some(val) = specs.ldelta { val } else { 0.0 };
    vals.lrdsw = if let Some(val) = specs.lrdsw { val } else { 0.0 };
    vals.lrdw = if let Some(val) = specs.lrdw { val } else { 0.0 };
    vals.lrsw = if let Some(val) = specs.lrsw { val } else { 0.0 };
    vals.lprwb = if let Some(val) = specs.lprwb { val } else { 0.0 };
    vals.lprwg = if let Some(val) = specs.lprwg { val } else { 0.0 };
    vals.lprt = if let Some(val) = specs.lprt { val } else { 0.0 };
    vals.leta0 = if let Some(val) = specs.leta0 { val } else { 0.0 };
    vals.letab = if let Some(val) = specs.letab { val } else { -0.0 };
    vals.lpclm = if let Some(val) = specs.lpclm { val } else { 0.0 };
    vals.lpdiblc1 = if let Some(val) = specs.lpdiblc1 { val } else { 0.0 };
    vals.lpdiblc2 = if let Some(val) = specs.lpdiblc2 { val } else { 0.0 };
    vals.lpdiblcb = if let Some(val) = specs.lpdiblcb { val } else { 0.0 };
    vals.lpscbe1 = if let Some(val) = specs.lpscbe1 { val } else { 0.0 };
    vals.lpscbe2 = if let Some(val) = specs.lpscbe2 { val } else { 0.0 };
    vals.lpvag = if let Some(val) = specs.lpvag { val } else { 0.0 };
    vals.lwr = if let Some(val) = specs.lwr { val } else { 0.0 };
    vals.ldwg = if let Some(val) = specs.ldwg { val } else { 0.0 };
    vals.ldwb = if let Some(val) = specs.ldwb { val } else { 0.0 };
    vals.lb0 = if let Some(val) = specs.lb0 { val } else { 0.0 };
    vals.lb1 = if let Some(val) = specs.lb1 { val } else { 0.0 };
    vals.lalpha0 = if let Some(val) = specs.lalpha0 { val } else { 0.0 };
    vals.lalpha1 = if let Some(val) = specs.lalpha1 { val } else { 0.0 };
    vals.lbeta0 = if let Some(val) = specs.lbeta0 { val } else { 0.0 };
    vals.lagidl = if let Some(val) = specs.lagidl { val } else { 0.0 };
    vals.lbgidl = if let Some(val) = specs.lbgidl { val } else { 0.0 };
    vals.lcgidl = if let Some(val) = specs.lcgidl { val } else { 0.0 };
    vals.legidl = if let Some(val) = specs.legidl { val } else { 0.0 };
    vals.lrgidl = if let Some(val) = specs.lrgidl { val } else { 0.0 };
    vals.lkgidl = if let Some(val) = specs.lkgidl { val } else { 0.0 };
    vals.lfgidl = if let Some(val) = specs.lfgidl { val } else { 0.0 };

    /*Default value of lagisl, lbgisl, lcgisl, legisl, lrgisl, lkgisl, and lfgisl are set as follows */
    vals.lagisl = if let Some(val) = specs.lagisl { val } else { vals.lagidl };
    vals.lbgisl = if let Some(val) = specs.lbgisl { val } else { vals.lbgidl };
    vals.lcgisl = if let Some(val) = specs.lcgisl { val } else { vals.lcgidl };
    vals.legisl = if let Some(val) = specs.legisl { val } else { vals.legidl };
    vals.lrgisl = if let Some(val) = specs.lrgisl { val } else { vals.lrgidl };
    vals.lkgisl = if let Some(val) = specs.lkgisl { val } else { vals.lkgidl };
    vals.lfgisl = if let Some(val) = specs.lfgisl { val } else { vals.lfgidl };

    vals.laigc = if let Some(val) = specs.laigc { val } else { 0.0 };
    vals.lbigc = if let Some(val) = specs.lbigc { val } else { 0.0 };
    vals.lcigc = if let Some(val) = specs.lcigc { val } else { 0.0 };
    if !specs.aigsd.is_some() && (specs.aigs.is_some() || specs.aigd.is_some()) {
        vals.laigs = if let Some(val) = specs.laigs { val } else { 0.0 };
        vals.laigd = if let Some(val) = specs.laigd { val } else { 0.0 };
    } else {
        vals.laigsd = if let Some(val) = specs.laigsd { val } else { 0.0 };
        vals.laigs = vals.laigsd;
        vals.laigd = vals.laigsd;
    }
    if !specs.bigsd.is_some() && (specs.bigs.is_some() || specs.bigd.is_some()) {
        vals.lbigs = if let Some(val) = specs.lbigs { val } else { 0.0 };
        vals.lbigd = if let Some(val) = specs.lbigd { val } else { 0.0 };
    } else {
        vals.lbigsd = if let Some(val) = specs.lbigsd { val } else { 0.0 };
        vals.lbigs = vals.lbigsd;
        vals.lbigd = vals.lbigsd;
    }
    if !specs.cigsd.is_some() && (specs.cigs.is_some() || specs.cigd.is_some()) {
        vals.lcigs = if let Some(val) = specs.lcigs { val } else { 0.0 };
        vals.lcigd = if let Some(val) = specs.lcigd { val } else { 0.0 };
    } else {
        vals.lcigsd = if let Some(val) = specs.lcigsd { val } else { 0.0 };
        vals.lcigs = vals.lcigsd;
        vals.lcigd = vals.lcigsd;
    }
    vals.laigbacc = if let Some(val) = specs.laigbacc { val } else { 0.0 };
    vals.lbigbacc = if let Some(val) = specs.lbigbacc { val } else { 0.0 };
    vals.lcigbacc = if let Some(val) = specs.lcigbacc { val } else { 0.0 };
    vals.laigbinv = if let Some(val) = specs.laigbinv { val } else { 0.0 };
    vals.lbigbinv = if let Some(val) = specs.lbigbinv { val } else { 0.0 };
    vals.lcigbinv = if let Some(val) = specs.lcigbinv { val } else { 0.0 };
    vals.lnigc = if let Some(val) = specs.lnigc { val } else { 0.0 };
    vals.lnigbinv = if let Some(val) = specs.lnigbinv { val } else { 0.0 };
    vals.lnigbacc = if let Some(val) = specs.lnigbacc { val } else { 0.0 };
    vals.lntox = if let Some(val) = specs.lntox { val } else { 0.0 };
    vals.leigbinv = if let Some(val) = specs.leigbinv { val } else { 0.0 };
    vals.lpigcd = if let Some(val) = specs.lpigcd { val } else { 0.0 };
    vals.lpoxedge = if let Some(val) = specs.lpoxedge { val } else { 0.0 };
    vals.lxrcrg1 = if let Some(val) = specs.lxrcrg1 { val } else { 0.0 };
    vals.lxrcrg2 = if let Some(val) = specs.lxrcrg2 { val } else { 0.0 };
    vals.leu = if let Some(val) = specs.leu { val } else { 0.0 };
    vals.lucs = if let Some(val) = specs.lucs { val } else { 0.0 };
    vals.lvfb = if let Some(val) = specs.lvfb { val } else { 0.0 };
    vals.llambda = if let Some(val) = specs.llambda { val } else { 0.0 };
    vals.lvtl = if let Some(val) = specs.lvtl { val } else { 0.0 };
    vals.lxn = if let Some(val) = specs.lxn { val } else { 0.0 };
    vals.lvfbsdoff = if let Some(val) = specs.lvfbsdoff { val } else { 0.0 };
    vals.ltvfbsdoff = if let Some(val) = specs.ltvfbsdoff { val } else { 0.0 };
    vals.ltvoff = if let Some(val) = specs.ltvoff { val } else { 0.0 };
    vals.ltnfactor = if let Some(val) = specs.ltnfactor { val } else { 0.0 };
    vals.lteta0 = if let Some(val) = specs.lteta0 { val } else { 0.0 };
    vals.ltvoffcv = if let Some(val) = specs.ltvoffcv { val } else { 0.0 };

    vals.lcgsl = if let Some(val) = specs.lcgsl { val } else { 0.0 };
    vals.lcgdl = if let Some(val) = specs.lcgdl { val } else { 0.0 };
    vals.lckappas = if let Some(val) = specs.lckappas { val } else { 0.0 };
    vals.lckappad = if let Some(val) = specs.lckappad { val } else { 0.0 };
    vals.lclc = if let Some(val) = specs.lclc { val } else { 0.0 };
    vals.lcle = if let Some(val) = specs.lcle { val } else { 0.0 };
    vals.lcf = if let Some(val) = specs.lcf { val } else { 0.0 };
    vals.lvfbcv = if let Some(val) = specs.lvfbcv { val } else { 0.0 };
    vals.lacde = if let Some(val) = specs.lacde { val } else { 0.0 };
    vals.lmoin = if let Some(val) = specs.lmoin { val } else { 0.0 };
    vals.lnoff = if let Some(val) = specs.lnoff { val } else { 0.0 };
    vals.lvoffcv = if let Some(val) = specs.lvoffcv { val } else { 0.0 };

    /* Width dependence */
    vals.wcdsc = if let Some(val) = specs.wcdsc { val } else { 0.0 };
    vals.wcdscb = if let Some(val) = specs.wcdscb { val } else { 0.0 };
    vals.wcdscd = if let Some(val) = specs.wcdscd { val } else { 0.0 };
    vals.wcit = if let Some(val) = specs.wcit { val } else { 0.0 };
    vals.wnfactor = if let Some(val) = specs.wnfactor { val } else { 0.0 };
    vals.wxj = if let Some(val) = specs.wxj { val } else { 0.0 };
    vals.wvsat = if let Some(val) = specs.wvsat { val } else { 0.0 };
    vals.wat = if let Some(val) = specs.wat { val } else { 0.0 };
    vals.wa0 = if let Some(val) = specs.wa0 { val } else { 0.0 };
    vals.wags = if let Some(val) = specs.wags { val } else { 0.0 };
    vals.wa1 = if let Some(val) = specs.wa1 { val } else { 0.0 };
    vals.wa2 = if let Some(val) = specs.wa2 { val } else { 0.0 };
    vals.wketa = if let Some(val) = specs.wketa { val } else { 0.0 };
    vals.wnsub = if let Some(val) = specs.wnsub { val } else { 0.0 };
    vals.wndep = if let Some(val) = specs.wndep { val } else { 0.0 };
    vals.wnsd = if let Some(val) = specs.wnsd { val } else { 0.0 };
    vals.wphin = if let Some(val) = specs.wphin { val } else { 0.0 };
    vals.wngate = if let Some(val) = specs.wngate { val } else { 0.0 };
    vals.wvbm = if let Some(val) = specs.wvbm { val } else { 0.0 };
    vals.wxt = if let Some(val) = specs.wxt { val } else { 0.0 };
    vals.wk1 = if let Some(val) = specs.wk1 { val } else { 0.0 };
    vals.wkt1 = if let Some(val) = specs.wkt1 { val } else { 0.0 };
    vals.wkt1l = if let Some(val) = specs.wkt1l { val } else { 0.0 };
    vals.wkt2 = if let Some(val) = specs.wkt2 { val } else { 0.0 };
    vals.wk2 = if let Some(val) = specs.wk2 { val } else { 0.0 };
    vals.wk3 = if let Some(val) = specs.wk3 { val } else { 0.0 };
    vals.wk3b = if let Some(val) = specs.wk3b { val } else { 0.0 };
    vals.ww0 = if let Some(val) = specs.ww0 { val } else { 0.0 };
    vals.wlpe0 = if let Some(val) = specs.wlpe0 { val } else { 0.0 };
    vals.wlpeb = if let Some(val) = specs.wlpeb { val } else { 0.0 };
    vals.wdvtp0 = if let Some(val) = specs.wdvtp0 { val } else { 0.0 };
    vals.wdvtp1 = if let Some(val) = specs.wdvtp1 { val } else { 0.0 };
    vals.wdvtp2 = if let Some(val) = specs.wdvtp2 { val } else { 0.0 };
    vals.wdvtp3 = if let Some(val) = specs.wdvtp3 { val } else { 0.0 };
    vals.wdvtp4 = if let Some(val) = specs.wdvtp4 { val } else { 0.0 };
    vals.wdvtp5 = if let Some(val) = specs.wdvtp5 { val } else { 0.0 };
    vals.wdvt0 = if let Some(val) = specs.wdvt0 { val } else { 0.0 };
    vals.wdvt1 = if let Some(val) = specs.wdvt1 { val } else { 0.0 };
    vals.wdvt2 = if let Some(val) = specs.wdvt2 { val } else { 0.0 };
    vals.wdvt0w = if let Some(val) = specs.wdvt0w { val } else { 0.0 };
    vals.wdvt1w = if let Some(val) = specs.wdvt1w { val } else { 0.0 };
    vals.wdvt2w = if let Some(val) = specs.wdvt2w { val } else { 0.0 };
    vals.wdrout = if let Some(val) = specs.wdrout { val } else { 0.0 };
    vals.wdsub = if let Some(val) = specs.wdsub { val } else { 0.0 };
    vals.wvth0 = if let Some(val) = specs.wvth0 { val } else { 0.0 };
    vals.wua = if let Some(val) = specs.wua { val } else { 0.0 };
    vals.wua1 = if let Some(val) = specs.wua1 { val } else { 0.0 };
    vals.wub = if let Some(val) = specs.wub { val } else { 0.0 };
    vals.wub1 = if let Some(val) = specs.wub1 { val } else { 0.0 };
    vals.wuc = if let Some(val) = specs.wuc { val } else { 0.0 };
    vals.wuc1 = if let Some(val) = specs.wuc1 { val } else { 0.0 };
    vals.wud = if let Some(val) = specs.wud { val } else { 0.0 };
    vals.wud1 = if let Some(val) = specs.wud1 { val } else { 0.0 };
    vals.wup = if let Some(val) = specs.wup { val } else { 0.0 };
    vals.wlp = if let Some(val) = specs.wlp { val } else { 0.0 };
    vals.wu0 = if let Some(val) = specs.wu0 { val } else { 0.0 };
    vals.wute = if let Some(val) = specs.wute { val } else { 0.0 };
    vals.wucste = if let Some(val) = specs.wucste { val } else { 0.0 };
    vals.wvoff = if let Some(val) = specs.wvoff { val } else { 0.0 };
    vals.wminv = if let Some(val) = specs.wminv { val } else { 0.0 };
    vals.wminvcv = if let Some(val) = specs.wminvcv { val } else { 0.0 };
    vals.wfprout = if let Some(val) = specs.wfprout { val } else { 0.0 };
    vals.wpdits = if let Some(val) = specs.wpdits { val } else { 0.0 };
    vals.wpditsd = if let Some(val) = specs.wpditsd { val } else { 0.0 };
    vals.wdelta = if let Some(val) = specs.wdelta { val } else { 0.0 };
    vals.wrdsw = if let Some(val) = specs.wrdsw { val } else { 0.0 };
    vals.wrdw = if let Some(val) = specs.wrdw { val } else { 0.0 };
    vals.wrsw = if let Some(val) = specs.wrsw { val } else { 0.0 };
    vals.wprwb = if let Some(val) = specs.wprwb { val } else { 0.0 };
    vals.wprwg = if let Some(val) = specs.wprwg { val } else { 0.0 };
    vals.wprt = if let Some(val) = specs.wprt { val } else { 0.0 };
    vals.weta0 = if let Some(val) = specs.weta0 { val } else { 0.0 };
    vals.wetab = if let Some(val) = specs.wetab { val } else { 0.0 };
    vals.wpclm = if let Some(val) = specs.wpclm { val } else { 0.0 };
    vals.wpdiblc1 = if let Some(val) = specs.wpdiblc1 { val } else { 0.0 };
    vals.wpdiblc2 = if let Some(val) = specs.wpdiblc2 { val } else { 0.0 };
    vals.wpdiblcb = if let Some(val) = specs.wpdiblcb { val } else { 0.0 };
    vals.wpscbe1 = if let Some(val) = specs.wpscbe1 { val } else { 0.0 };
    vals.wpscbe2 = if let Some(val) = specs.wpscbe2 { val } else { 0.0 };
    vals.wpvag = if let Some(val) = specs.wpvag { val } else { 0.0 };
    vals.wwr = if let Some(val) = specs.wwr { val } else { 0.0 };
    vals.wdwg = if let Some(val) = specs.wdwg { val } else { 0.0 };
    vals.wdwb = if let Some(val) = specs.wdwb { val } else { 0.0 };
    vals.wb0 = if let Some(val) = specs.wb0 { val } else { 0.0 };
    vals.wb1 = if let Some(val) = specs.wb1 { val } else { 0.0 };
    vals.walpha0 = if let Some(val) = specs.walpha0 { val } else { 0.0 };
    vals.walpha1 = if let Some(val) = specs.walpha1 { val } else { 0.0 };
    vals.wbeta0 = if let Some(val) = specs.wbeta0 { val } else { 0.0 };
    vals.wagidl = if let Some(val) = specs.wagidl { val } else { 0.0 };
    vals.wbgidl = if let Some(val) = specs.wbgidl { val } else { 0.0 };
    vals.wcgidl = if let Some(val) = specs.wcgidl { val } else { 0.0 };
    vals.wegidl = if let Some(val) = specs.wegidl { val } else { 0.0 };
    vals.wrgidl = if let Some(val) = specs.wrgidl { val } else { 0.0 };
    vals.wkgidl = if let Some(val) = specs.wkgidl { val } else { 0.0 };
    vals.wfgidl = if let Some(val) = specs.wfgidl { val } else { 0.0 };

    /*Default value of wagisl, wbgisl, wcgisl, wegisl, wrgisl, wkgisl, and wfgisl are set as follows */
    vals.wagisl = if let Some(val) = specs.wagisl { val } else { vals.wagidl };
    vals.wbgisl = if let Some(val) = specs.wbgisl { val } else { vals.wbgidl };
    vals.wcgisl = if let Some(val) = specs.wcgisl { val } else { vals.wcgidl };
    vals.wegisl = if let Some(val) = specs.wegisl { val } else { vals.wegidl };
    vals.wrgisl = if let Some(val) = specs.wrgisl { val } else { vals.wrgidl };
    vals.wkgisl = if let Some(val) = specs.wkgisl { val } else { vals.wkgidl };
    vals.wfgisl = if let Some(val) = specs.wfgisl { val } else { vals.wfgidl };

    vals.waigc = if let Some(val) = specs.waigc { val } else { 0.0 };
    vals.wbigc = if let Some(val) = specs.wbigc { val } else { 0.0 };
    vals.wcigc = if let Some(val) = specs.wcigc { val } else { 0.0 };
    if !specs.aigsd.is_some() && (specs.aigs.is_some() || specs.aigd.is_some()) {
        vals.waigs = if let Some(val) = specs.waigs { val } else { 0.0 };
        vals.waigd = if let Some(val) = specs.waigd { val } else { 0.0 };
    } else {
        vals.waigsd = if let Some(val) = specs.waigsd { val } else { 0.0 };
        vals.waigs = vals.waigsd;
        vals.waigd = vals.waigsd;
    }
    if !specs.bigsd.is_some() && (specs.bigs.is_some() || specs.bigd.is_some()) {
        vals.wbigs = if let Some(val) = specs.wbigs { val } else { 0.0 };
        vals.wbigd = if let Some(val) = specs.wbigd { val } else { 0.0 };
    } else {
        vals.wbigsd = if let Some(val) = specs.wbigsd { val } else { 0.0 };
        vals.wbigs = vals.wbigsd;
        vals.wbigd = vals.wbigsd;
    }
    if !specs.cigsd.is_some() && (specs.cigs.is_some() || specs.cigd.is_some()) {
        vals.wcigs = if let Some(val) = specs.wcigs { val } else { 0.0 };
        vals.wcigd = if let Some(val) = specs.wcigd { val } else { 0.0 };
    } else {
        vals.wcigsd = if let Some(val) = specs.wcigsd { val } else { 0.0 };
        vals.wcigs = vals.wcigsd;
        vals.wcigd = vals.wcigsd;
    }
    vals.waigbacc = if let Some(val) = specs.waigbacc { val } else { 0.0 };
    vals.wbigbacc = if let Some(val) = specs.wbigbacc { val } else { 0.0 };
    vals.wcigbacc = if let Some(val) = specs.wcigbacc { val } else { 0.0 };
    vals.waigbinv = if let Some(val) = specs.waigbinv { val } else { 0.0 };
    vals.wbigbinv = if let Some(val) = specs.wbigbinv { val } else { 0.0 };
    vals.wcigbinv = if let Some(val) = specs.wcigbinv { val } else { 0.0 };
    vals.wnigc = if let Some(val) = specs.wnigc { val } else { 0.0 };
    vals.wnigbinv = if let Some(val) = specs.wnigbinv { val } else { 0.0 };
    vals.wnigbacc = if let Some(val) = specs.wnigbacc { val } else { 0.0 };
    vals.wntox = if let Some(val) = specs.wntox { val } else { 0.0 };
    vals.weigbinv = if let Some(val) = specs.weigbinv { val } else { 0.0 };
    vals.wpigcd = if let Some(val) = specs.wpigcd { val } else { 0.0 };
    vals.wpoxedge = if let Some(val) = specs.wpoxedge { val } else { 0.0 };
    vals.wxrcrg1 = if let Some(val) = specs.wxrcrg1 { val } else { 0.0 };
    vals.wxrcrg2 = if let Some(val) = specs.wxrcrg2 { val } else { 0.0 };
    vals.weu = if let Some(val) = specs.weu { val } else { 0.0 };
    vals.wucs = if let Some(val) = specs.wucs { val } else { 0.0 };
    vals.wvfb = if let Some(val) = specs.wvfb { val } else { 0.0 };
    vals.wlambda = if let Some(val) = specs.wlambda { val } else { 0.0 };
    vals.wvtl = if let Some(val) = specs.wvtl { val } else { 0.0 };
    vals.wxn = if let Some(val) = specs.wxn { val } else { 0.0 };
    vals.wvfbsdoff = if let Some(val) = specs.wvfbsdoff { val } else { 0.0 };
    vals.wtvfbsdoff = if let Some(val) = specs.wtvfbsdoff { val } else { 0.0 };
    vals.wtvoff = if let Some(val) = specs.wtvoff { val } else { 0.0 };
    vals.wtnfactor = if let Some(val) = specs.wtnfactor { val } else { 0.0 };
    vals.wteta0 = if let Some(val) = specs.wteta0 { val } else { 0.0 };
    vals.wtvoffcv = if let Some(val) = specs.wtvoffcv { val } else { 0.0 };

    vals.wcgsl = if let Some(val) = specs.wcgsl { val } else { 0.0 };
    vals.wcgdl = if let Some(val) = specs.wcgdl { val } else { 0.0 };
    vals.wckappas = if let Some(val) = specs.wckappas { val } else { 0.0 };
    vals.wckappad = if let Some(val) = specs.wckappad { val } else { 0.0 };
    vals.wcf = if let Some(val) = specs.wcf { val } else { 0.0 };
    vals.wclc = if let Some(val) = specs.wclc { val } else { 0.0 };
    vals.wcle = if let Some(val) = specs.wcle { val } else { 0.0 };
    vals.wvfbcv = if let Some(val) = specs.wvfbcv { val } else { 0.0 };
    vals.wacde = if let Some(val) = specs.wacde { val } else { 0.0 };
    vals.wmoin = if let Some(val) = specs.wmoin { val } else { 0.0 };
    vals.wnoff = if let Some(val) = specs.wnoff { val } else { 0.0 };
    vals.wvoffcv = if let Some(val) = specs.wvoffcv { val } else { 0.0 };

    /* Cross-term dependence */
    vals.pcdsc = if let Some(val) = specs.pcdsc { val } else { 0.0 };
    vals.pcdscb = if let Some(val) = specs.pcdscb { val } else { 0.0 };
    vals.pcdscd = if let Some(val) = specs.pcdscd { val } else { 0.0 };
    vals.pcit = if let Some(val) = specs.pcit { val } else { 0.0 };
    vals.pnfactor = if let Some(val) = specs.pnfactor { val } else { 0.0 };
    vals.pxj = if let Some(val) = specs.pxj { val } else { 0.0 };
    vals.pvsat = if let Some(val) = specs.pvsat { val } else { 0.0 };
    vals.pat = if let Some(val) = specs.pat { val } else { 0.0 };
    vals.pa0 = if let Some(val) = specs.pa0 { val } else { 0.0 };
    vals.pags = if let Some(val) = specs.pags { val } else { 0.0 };
    vals.pa1 = if let Some(val) = specs.pa1 { val } else { 0.0 };
    vals.pa2 = if let Some(val) = specs.pa2 { val } else { 0.0 };
    vals.pketa = if let Some(val) = specs.pketa { val } else { 0.0 };
    vals.pnsub = if let Some(val) = specs.pnsub { val } else { 0.0 };
    vals.pndep = if let Some(val) = specs.pndep { val } else { 0.0 };
    vals.pnsd = if let Some(val) = specs.pnsd { val } else { 0.0 };
    vals.pphin = if let Some(val) = specs.pphin { val } else { 0.0 };
    vals.pngate = if let Some(val) = specs.pngate { val } else { 0.0 };
    vals.pvbm = if let Some(val) = specs.pvbm { val } else { 0.0 };
    vals.pxt = if let Some(val) = specs.pxt { val } else { 0.0 };
    vals.pk1 = if let Some(val) = specs.pk1 { val } else { 0.0 };
    vals.pkt1 = if let Some(val) = specs.pkt1 { val } else { 0.0 };
    vals.pkt1l = if let Some(val) = specs.pkt1l { val } else { 0.0 };
    vals.pk2 = if let Some(val) = specs.pk2 { val } else { 0.0 };
    vals.pkt2 = if let Some(val) = specs.pkt2 { val } else { 0.0 };
    vals.pk3 = if let Some(val) = specs.pk3 { val } else { 0.0 };
    vals.pk3b = if let Some(val) = specs.pk3b { val } else { 0.0 };
    vals.pw0 = if let Some(val) = specs.pw0 { val } else { 0.0 };
    vals.plpe0 = if let Some(val) = specs.plpe0 { val } else { 0.0 };
    vals.plpeb = if let Some(val) = specs.plpeb { val } else { 0.0 };
    vals.pdvtp0 = if let Some(val) = specs.pdvtp0 { val } else { 0.0 };
    vals.pdvtp1 = if let Some(val) = specs.pdvtp1 { val } else { 0.0 };
    vals.pdvtp2 = if let Some(val) = specs.pdvtp2 { val } else { 0.0 };
    vals.pdvtp3 = if let Some(val) = specs.pdvtp3 { val } else { 0.0 };
    vals.pdvtp4 = if let Some(val) = specs.pdvtp4 { val } else { 0.0 };
    vals.pdvtp5 = if let Some(val) = specs.pdvtp5 { val } else { 0.0 };
    vals.pdvt0 = if let Some(val) = specs.pdvt0 { val } else { 0.0 };
    vals.pdvt1 = if let Some(val) = specs.pdvt1 { val } else { 0.0 };
    vals.pdvt2 = if let Some(val) = specs.pdvt2 { val } else { 0.0 };
    vals.pdvt0w = if let Some(val) = specs.pdvt0w { val } else { 0.0 };
    vals.pdvt1w = if let Some(val) = specs.pdvt1w { val } else { 0.0 };
    vals.pdvt2w = if let Some(val) = specs.pdvt2w { val } else { 0.0 };
    vals.pdrout = if let Some(val) = specs.pdrout { val } else { 0.0 };
    vals.pdsub = if let Some(val) = specs.pdsub { val } else { 0.0 };
    vals.pvth0 = if let Some(val) = specs.pvth0 { val } else { 0.0 };
    vals.pua = if let Some(val) = specs.pua { val } else { 0.0 };
    vals.pua1 = if let Some(val) = specs.pua1 { val } else { 0.0 };
    vals.r#pub = if let Some(val) = specs.r#pub { val } else { 0.0 };
    vals.pub1 = if let Some(val) = specs.pub1 { val } else { 0.0 };
    vals.puc = if let Some(val) = specs.puc { val } else { 0.0 };
    vals.puc1 = if let Some(val) = specs.puc1 { val } else { 0.0 };
    vals.pud = if let Some(val) = specs.pud { val } else { 0.0 };
    vals.pud1 = if let Some(val) = specs.pud1 { val } else { 0.0 };
    vals.pup = if let Some(val) = specs.pup { val } else { 0.0 };
    vals.plp = if let Some(val) = specs.plp { val } else { 0.0 };
    vals.pu0 = if let Some(val) = specs.pu0 { val } else { 0.0 };
    vals.pute = if let Some(val) = specs.pute { val } else { 0.0 };
    vals.pucste = if let Some(val) = specs.pucste { val } else { 0.0 };
    vals.pvoff = if let Some(val) = specs.pvoff { val } else { 0.0 };
    vals.pminv = if let Some(val) = specs.pminv { val } else { 0.0 };
    vals.pminvcv = if let Some(val) = specs.pminvcv { val } else { 0.0 };
    vals.pfprout = if let Some(val) = specs.pfprout { val } else { 0.0 };
    vals.ppdits = if let Some(val) = specs.ppdits { val } else { 0.0 };
    vals.ppditsd = if let Some(val) = specs.ppditsd { val } else { 0.0 };
    vals.pdelta = if let Some(val) = specs.pdelta { val } else { 0.0 };
    vals.prdsw = if let Some(val) = specs.prdsw { val } else { 0.0 };
    vals.prdw = if let Some(val) = specs.prdw { val } else { 0.0 };
    vals.prsw = if let Some(val) = specs.prsw { val } else { 0.0 };
    vals.pprwb = if let Some(val) = specs.pprwb { val } else { 0.0 };
    vals.pprwg = if let Some(val) = specs.pprwg { val } else { 0.0 };
    vals.pprt = if let Some(val) = specs.pprt { val } else { 0.0 };
    vals.peta0 = if let Some(val) = specs.peta0 { val } else { 0.0 };
    vals.petab = if let Some(val) = specs.petab { val } else { 0.0 };
    vals.ppclm = if let Some(val) = specs.ppclm { val } else { 0.0 };
    vals.ppdiblc1 = if let Some(val) = specs.ppdiblc1 { val } else { 0.0 };
    vals.ppdiblc2 = if let Some(val) = specs.ppdiblc2 { val } else { 0.0 };
    vals.ppdiblcb = if let Some(val) = specs.ppdiblcb { val } else { 0.0 };
    vals.ppscbe1 = if let Some(val) = specs.ppscbe1 { val } else { 0.0 };
    vals.ppscbe2 = if let Some(val) = specs.ppscbe2 { val } else { 0.0 };
    vals.ppvag = if let Some(val) = specs.ppvag { val } else { 0.0 };
    vals.pwr = if let Some(val) = specs.pwr { val } else { 0.0 };
    vals.pdwg = if let Some(val) = specs.pdwg { val } else { 0.0 };
    vals.pdwb = if let Some(val) = specs.pdwb { val } else { 0.0 };
    vals.pb0 = if let Some(val) = specs.pb0 { val } else { 0.0 };
    vals.pb1 = if let Some(val) = specs.pb1 { val } else { 0.0 };
    vals.palpha0 = if let Some(val) = specs.palpha0 { val } else { 0.0 };
    vals.palpha1 = if let Some(val) = specs.palpha1 { val } else { 0.0 };
    vals.pbeta0 = if let Some(val) = specs.pbeta0 { val } else { 0.0 };
    vals.pagidl = if let Some(val) = specs.pagidl { val } else { 0.0 };
    vals.pbgidl = if let Some(val) = specs.pbgidl { val } else { 0.0 };
    vals.pcgidl = if let Some(val) = specs.pcgidl { val } else { 0.0 };
    vals.pegidl = if let Some(val) = specs.pegidl { val } else { 0.0 };
    vals.prgidl = if let Some(val) = specs.prgidl { val } else { 0.0 };
    vals.pkgidl = if let Some(val) = specs.pkgidl { val } else { 0.0 };
    vals.pfgidl = if let Some(val) = specs.pfgidl { val } else { 0.0 };

    /*Default value of pagisl, pbgisl, pcgisl, pegisl, prgisl, pkgisl, and pfgisl are set as follows */
    vals.pagisl = if let Some(val) = specs.pagisl { val } else { vals.pagidl };
    vals.pbgisl = if let Some(val) = specs.pbgisl { val } else { vals.pbgidl };
    vals.pcgisl = if let Some(val) = specs.pcgisl { val } else { vals.pcgidl };
    vals.pegisl = if let Some(val) = specs.pegisl { val } else { vals.pegidl };
    vals.prgisl = if let Some(val) = specs.prgisl { val } else { vals.prgidl };
    vals.pkgisl = if let Some(val) = specs.pkgisl { val } else { vals.pkgidl };
    vals.pfgisl = if let Some(val) = specs.pfgisl { val } else { vals.pfgidl };

    vals.paigc = if let Some(val) = specs.paigc { val } else { 0.0 };
    vals.pbigc = if let Some(val) = specs.pbigc { val } else { 0.0 };
    vals.pcigc = if let Some(val) = specs.pcigc { val } else { 0.0 };
    if !specs.aigsd.is_some() && (specs.aigs.is_some() || specs.aigd.is_some()) {
        vals.paigs = if let Some(val) = specs.paigs { val } else { 0.0 };
        vals.paigd = if let Some(val) = specs.paigd { val } else { 0.0 };
    } else {
        vals.paigsd = if let Some(val) = specs.paigsd { val } else { 0.0 };
        vals.paigs = vals.paigsd;
        vals.paigd = vals.paigsd;
    }
    if !specs.bigsd.is_some() && (specs.bigs.is_some() || specs.bigd.is_some()) {
        vals.pbigs = if let Some(val) = specs.pbigs { val } else { 0.0 };
        vals.pbigd = if let Some(val) = specs.pbigd { val } else { 0.0 };
    } else {
        vals.pbigsd = if let Some(val) = specs.pbigsd { val } else { 0.0 };
        vals.pbigs = vals.pbigsd;
        vals.pbigd = vals.pbigsd;
    }
    if !specs.cigsd.is_some() && (specs.cigs.is_some() || specs.cigd.is_some()) {
        vals.pcigs = if let Some(val) = specs.pcigs { val } else { 0.0 };
        vals.pcigd = if let Some(val) = specs.pcigd { val } else { 0.0 };
    } else {
        vals.pcigsd = if let Some(val) = specs.pcigsd { val } else { 0.0 };
        vals.pcigs = vals.pcigsd;
        vals.pcigd = vals.pcigsd;
    }
    vals.paigbacc = if let Some(val) = specs.paigbacc { val } else { 0.0 };
    vals.pbigbacc = if let Some(val) = specs.pbigbacc { val } else { 0.0 };
    vals.pcigbacc = if let Some(val) = specs.pcigbacc { val } else { 0.0 };
    vals.paigbinv = if let Some(val) = specs.paigbinv { val } else { 0.0 };
    vals.pbigbinv = if let Some(val) = specs.pbigbinv { val } else { 0.0 };
    vals.pcigbinv = if let Some(val) = specs.pcigbinv { val } else { 0.0 };
    vals.pnigc = if let Some(val) = specs.pnigc { val } else { 0.0 };
    vals.pnigbinv = if let Some(val) = specs.pnigbinv { val } else { 0.0 };
    vals.pnigbacc = if let Some(val) = specs.pnigbacc { val } else { 0.0 };
    vals.pntox = if let Some(val) = specs.pntox { val } else { 0.0 };
    vals.peigbinv = if let Some(val) = specs.peigbinv { val } else { 0.0 };
    vals.ppigcd = if let Some(val) = specs.ppigcd { val } else { 0.0 };
    vals.ppoxedge = if let Some(val) = specs.ppoxedge { val } else { 0.0 };
    vals.pxrcrg1 = if let Some(val) = specs.pxrcrg1 { val } else { 0.0 };
    vals.pxrcrg2 = if let Some(val) = specs.pxrcrg2 { val } else { 0.0 };
    vals.peu = if let Some(val) = specs.peu { val } else { 0.0 };
    vals.pucs = if let Some(val) = specs.pucs { val } else { 0.0 };
    vals.pvfb = if let Some(val) = specs.pvfb { val } else { 0.0 };
    vals.plambda = if let Some(val) = specs.plambda { val } else { 0.0 };
    vals.pvtl = if let Some(val) = specs.pvtl { val } else { 0.0 };
    vals.pxn = if let Some(val) = specs.pxn { val } else { 0.0 };
    vals.pvfbsdoff = if let Some(val) = specs.pvfbsdoff { val } else { 0.0 };
    vals.ptvfbsdoff = if let Some(val) = specs.ptvfbsdoff { val } else { 0.0 };
    vals.ptvoff = if let Some(val) = specs.ptvoff { val } else { 0.0 };
    vals.ptnfactor = if let Some(val) = specs.ptnfactor { val } else { 0.0 };
    vals.pteta0 = if let Some(val) = specs.pteta0 { val } else { 0.0 };
    vals.ptvoffcv = if let Some(val) = specs.ptvoffcv { val } else { 0.0 };
    vals.pcgsl = if let Some(val) = specs.pcgsl { val } else { 0.0 };
    vals.pcgdl = if let Some(val) = specs.pcgdl { val } else { 0.0 };
    vals.pckappas = if let Some(val) = specs.pckappas { val } else { 0.0 };
    vals.pckappad = if let Some(val) = specs.pckappad { val } else { 0.0 };
    vals.pcf = if let Some(val) = specs.pcf { val } else { 0.0 };
    vals.pclc = if let Some(val) = specs.pclc { val } else { 0.0 };
    vals.pcle = if let Some(val) = specs.pcle { val } else { 0.0 };
    vals.pvfbcv = if let Some(val) = specs.pvfbcv { val } else { 0.0 };
    vals.pacde = if let Some(val) = specs.pacde { val } else { 0.0 };
    vals.pmoin = if let Some(val) = specs.pmoin { val } else { 0.0 };
    vals.pnoff = if let Some(val) = specs.pnoff { val } else { 0.0 };
    vals.pvoffcv = if let Some(val) = specs.pvoffcv { val } else { 0.0 };
    vals.gamma1 = if let Some(val) = specs.gamma1 { val } else { 0.0 };
    vals.lgamma1 = if let Some(val) = specs.lgamma1 { val } else { 0.0 };
    vals.wgamma1 = if let Some(val) = specs.wgamma1 { val } else { 0.0 };
    vals.pgamma1 = if let Some(val) = specs.pgamma1 { val } else { 0.0 };
    vals.gamma2 = if let Some(val) = specs.gamma2 { val } else { 0.0 };
    vals.lgamma2 = if let Some(val) = specs.lgamma2 { val } else { 0.0 };
    vals.wgamma2 = if let Some(val) = specs.wgamma2 { val } else { 0.0 };
    vals.pgamma2 = if let Some(val) = specs.pgamma2 { val } else { 0.0 };
    vals.vbx = if let Some(val) = specs.vbx { val } else { 0.0 };
    vals.lvbx = if let Some(val) = specs.lvbx { val } else { 0.0 };
    vals.wvbx = if let Some(val) = specs.wvbx { val } else { 0.0 };
    vals.pvbx = if let Some(val) = specs.pvbx { val } else { 0.0 };

    vals.tnom = if let Some(val) = specs.tnom { val } else { tnom };
    vals.lint = if let Some(val) = specs.lint { val } else { 0.0 };
    vals.ll = if let Some(val) = specs.ll { val } else { 0.0 };
    vals.llc = if let Some(val) = specs.llc { val } else { vals.ll };
    vals.lln = if let Some(val) = specs.lln { val } else { 1.0 };
    vals.lw = if let Some(val) = specs.lw { val } else { 0.0 };
    vals.lwc = if let Some(val) = specs.lwc { val } else { vals.lw };
    vals.lwn = if let Some(val) = specs.lwn { val } else { 1.0 };
    vals.lwl = if let Some(val) = specs.lwl { val } else { 0.0 };
    vals.lwlc = if let Some(val) = specs.lwlc { val } else { vals.lwl };
    vals.lmin = if let Some(val) = specs.lmin { val } else { 0.0 };
    vals.lmax = if let Some(val) = specs.lmax { val } else { 1.0 };
    vals.wint = if let Some(val) = specs.wint { val } else { 0.0 };
    vals.wl = if let Some(val) = specs.wl { val } else { 0.0 };
    vals.wlc = if let Some(val) = specs.wlc { val } else { vals.wl };
    vals.wln = if let Some(val) = specs.wln { val } else { 1.0 };
    vals.ww = if let Some(val) = specs.ww { val } else { 0.0 };
    vals.wwc = if let Some(val) = specs.wwc { val } else { vals.ww };
    vals.wwn = if let Some(val) = specs.wwn { val } else { 1.0 };
    vals.wwl = if let Some(val) = specs.wwl { val } else { 0.0 };
    vals.wwlc = if let Some(val) = specs.wwlc { val } else { vals.wwl };
    vals.wmin = if let Some(val) = specs.wmin { val } else { 0.0 };
    vals.wmax = if let Some(val) = specs.wmax { val } else { 1.0 };
    vals.dwc = if let Some(val) = specs.dwc { val } else { vals.wint };
    vals.dlc = if let Some(val) = specs.dlc { val } else { vals.lint };
    vals.xl = if let Some(val) = specs.xl { val } else { 0.0 };
    vals.xw = if let Some(val) = specs.xw { val } else { 0.0 };
    vals.dlcig = if let Some(val) = specs.dlcig { val } else { vals.lint };
    vals.dlcigd = if let Some(val) = specs.dwj { val } else { vals.lint };
    vals.dwj = if let Some(val) = specs.dwj { val } else { vals.dwc };

    vals.cf = if let Some(val) = specs.cf {
        val
    } else {
        2.0 * vals.epsrox * EPS0 / PI * log(1.0 + 0.4e-6 / vals.toxe)
    };

    vals.xpart = if let Some(val) = specs.xpart { val } else { 0.0 };
    vals.rsh = if let Some(val) = specs.rsh { val } else { 0.0 };
    vals.cjs = if let Some(val) = specs.cjs { val } else { 5.0E-4 };
    vals.cjd = if let Some(val) = specs.cjd { val } else { vals.cjs };
    vals.cjsws = if let Some(val) = specs.cjsws { val } else { 5.0E-10 };
    vals.cjswd = if let Some(val) = specs.cjswd { val } else { vals.cjsws };
    vals.cjswgs = if let Some(val) = specs.cjswgs { val } else { vals.cjsws };
    vals.cjswgd = if let Some(val) = specs.cjswgd { val } else { vals.cjswgs };
    vals.jss = if let Some(val) = specs.jss { val } else { 1.0E-4 };
    vals.jsd = if let Some(val) = specs.jsd { val } else { vals.jss };
    vals.jsws = if let Some(val) = specs.jsws { val } else { 0.0 };
    vals.jswd = if let Some(val) = specs.jswd { val } else { vals.jsws };
    vals.jswgs = if let Some(val) = specs.jswgs { val } else { 0.0 };
    vals.jswgd = if let Some(val) = specs.jswgd { val } else { vals.jswgs };
    vals.pbs = if let Some(val) = specs.pbs { val } else { 1.0 };
    vals.pbd = if let Some(val) = specs.pbd { val } else { vals.pbs };
    vals.pbsws = if let Some(val) = specs.pbsws { val } else { 1.0 };
    vals.pbswd = if let Some(val) = specs.pbswd { val } else { vals.pbsws };
    vals.pbswgs = if let Some(val) = specs.pbswgs { val } else { vals.pbsws };
    vals.pbswgd = if let Some(val) = specs.pbswgd { val } else { vals.pbswgs };
    vals.mjs = if let Some(val) = specs.mjs { val } else { 0.5 };
    vals.mjd = if let Some(val) = specs.mjd { val } else { vals.mjs };
    vals.mjsws = if let Some(val) = specs.mjsws { val } else { 0.33 };
    vals.mjswd = if let Some(val) = specs.mjswd { val } else { vals.mjsws };
    vals.mjswgs = if let Some(val) = specs.mjswgs { val } else { vals.mjsws };
    vals.mjswgd = if let Some(val) = specs.mjswgd { val } else { vals.mjswgs };
    vals.njs = if let Some(val) = specs.njs { val } else { 1.0 };
    vals.njd = if let Some(val) = specs.njd { val } else { vals.njs };
    vals.xtis = if let Some(val) = specs.xtis { val } else { 3.0 };
    vals.xtid = if let Some(val) = specs.xtid { val } else { vals.xtis };

    vals.jtss = if let Some(val) = specs.jtss { val } else { 0.0 };
    vals.jtsd = if let Some(val) = specs.jtsd { val } else { vals.jtss };
    vals.jtssws = if let Some(val) = specs.jtssws { val } else { 0.0 };
    vals.jtsswd = if let Some(val) = specs.jtsswd { val } else { vals.jtssws };
    vals.jtsswgs = if let Some(val) = specs.jtsswgs { val } else { 0.0 };
    vals.jtsswgd = if let Some(val) = specs.jtsswgd { val } else { vals.jtsswgs };
    vals.jtweff = if let Some(val) = specs.jtweff { val } else { 0.0 };
    vals.njts = if let Some(val) = specs.njts { val } else { 20.0 };
    vals.njtssw = if let Some(val) = specs.njtssw { val } else { 20.0 };
    vals.njtsswg = if let Some(val) = specs.njtsswg { val } else { 20.0 };
    vals.njtsd = if let Some(val) = specs.njtsd { val } else { vals.njts };
    vals.njtsswd = if let Some(val) = specs.njtsswd { val } else { vals.njtssw };
    vals.njtsswgd = if let Some(val) = specs.njtsswgd { val } else { vals.njtsswg };

    vals.xtss = if let Some(val) = specs.xtss { val } else { 0.02 };
    vals.xtsd = if let Some(val) = specs.xtsd { val } else { vals.xtss };
    vals.xtssws = if let Some(val) = specs.xtssws { val } else { 0.02 };
    vals.xtsswd = if let Some(val) = specs.xtsswd { val } else { vals.xtssws };
    vals.xtsswgs = if let Some(val) = specs.xtsswgs { val } else { 0.02 };
    vals.xtsswgd = if let Some(val) = specs.xtsswgd { val } else { vals.xtsswgs };
    vals.tnjts = if let Some(val) = specs.tnjts { val } else { 0.0 };
    vals.tnjtssw = if let Some(val) = specs.tnjtssw { val } else { 0.0 };
    vals.tnjtsswg = if let Some(val) = specs.tnjtsswg { val } else { 0.0 };
    vals.tnjtsd = if let Some(val) = specs.tnjtsd { val } else { vals.tnjts };
    vals.tnjtsswd = if let Some(val) = specs.tnjtsswd { val } else { vals.tnjtssw };
    vals.tnjtsswgd = if let Some(val) = specs.tnjtsswgd { val } else { vals.tnjtsswg };
    vals.vtss = if let Some(val) = specs.vtss { val } else { 10.0 };
    vals.vtsd = if let Some(val) = specs.vtsd { val } else { vals.vtss };
    vals.vtssws = if let Some(val) = specs.vtssws { val } else { 10.0 };
    vals.vtsswd = if let Some(val) = specs.vtsswd { val } else { vals.vtssws };
    vals.vtsswgs = if let Some(val) = specs.vtsswgs { val } else { 10.0 };
    vals.vtsswgd = if let Some(val) = specs.vtsswgd { val } else { vals.vtsswgs };

    vals.noia = if let Some(val) = specs.noia {
        val
    } else {
        match vals.mos_type {
            NMOS => 6.25e41,
            PMOS => 6.188e40,
        }
    };
    vals.noib = if let Some(val) = specs.noib {
        val
    } else {
        match vals.mos_type {
            NMOS => 3.125e26,
            PMOS => 1.5e25,
        }
    };
    vals.noic = if let Some(val) = specs.noic { val } else { 8.75e9 };
    vals.em = if let Some(val) = specs.em { val } else { 4.1e7 }; /* V/m */
    vals.ef = if let Some(val) = specs.ef { val } else { 1.0 };
    vals.af = if let Some(val) = specs.af { val } else { 1.0 };
    vals.kf = if let Some(val) = specs.kf { val } else { 0.0 };

    /* stress effect */
    vals.saref = if let Some(val) = specs.saref { val } else { 1e-6 }; /* m */
    vals.sbref = if let Some(val) = specs.sbref { val } else { 1e-6 }; /* m */
    vals.wlod = if let Some(val) = specs.wlod { val } else { 0.0 }; /* m */
    vals.ku0 = if let Some(val) = specs.ku0 { val } else { 0.0 }; /* 1/m */
    vals.kvsat = if let Some(val) = specs.kvsat { val } else { 0.0 };
    vals.kvth0 = if let Some(val) = specs.kvth0 { val } else { 0.0 };
    vals.tku0 = if let Some(val) = specs.tku0 { val } else { 0.0 };
    vals.llodku0 = if let Some(val) = specs.llodku0 { val } else { 0.0 };
    vals.wlodku0 = if let Some(val) = specs.wlodku0 { val } else { 0.0 };
    vals.llodvth = if let Some(val) = specs.llodvth { val } else { 0.0 };
    vals.wlodvth = if let Some(val) = specs.wlodvth { val } else { 0.0 };
    vals.lku0 = if let Some(val) = specs.lku0 { val } else { 0.0 };
    vals.wku0 = if let Some(val) = specs.wku0 { val } else { 0.0 };
    vals.pku0 = if let Some(val) = specs.pku0 { val } else { 0.0 };
    vals.lkvth0 = if let Some(val) = specs.lkvth0 { val } else { 0.0 };
    vals.wkvth0 = if let Some(val) = specs.wkvth0 { val } else { 0.0 };
    vals.pkvth0 = if let Some(val) = specs.pkvth0 { val } else { 0.0 };
    vals.stk2 = if let Some(val) = specs.stk2 { val } else { 0.0 };
    vals.lodk2 = if let Some(val) = specs.lodk2 { val } else { 1.0 };
    vals.steta0 = if let Some(val) = specs.steta0 { val } else { 0.0 };
    vals.lodeta0 = if let Some(val) = specs.lodeta0 { val } else { 1.0 };

    /* Well Proximity Effect  */
    vals.web = if let Some(val) = specs.web { val } else { 0.0 };
    vals.wec = if let Some(val) = specs.wec { val } else { 0.0 };
    vals.kvth0we = if let Some(val) = specs.kvth0we { val } else { 0.0 };
    vals.k2we = if let Some(val) = specs.k2we { val } else { 0.0 };
    vals.ku0we = if let Some(val) = specs.ku0we { val } else { 0.0 };
    vals.scref = if let Some(val) = specs.scref { val } else { 1.0E-6 }; /* m */

    vals.lkvth0we = if let Some(val) = specs.lkvth0we { val } else { 0.0 };
    vals.lk2we = if let Some(val) = specs.lk2we { val } else { 0.0 };
    vals.lku0we = if let Some(val) = specs.lku0we { val } else { 0.0 };
    vals.wkvth0we = if let Some(val) = specs.wkvth0we { val } else { 0.0 };
    vals.wk2we = if let Some(val) = specs.wk2we { val } else { 0.0 };
    vals.wku0we = if let Some(val) = specs.wku0we { val } else { 0.0 };
    vals.pkvth0we = if let Some(val) = specs.pkvth0we { val } else { 0.0 };
    vals.pk2we = if let Some(val) = specs.pk2we { val } else { 0.0 };
    vals.pku0we = if let Some(val) = specs.pku0we { val } else { 0.0 };

    // FIXME: this whole mtrlmod params setup
    // For now just bail if they don't add up
    if vals.toxe != vals.toxp + vals.dtox {
        panic!("Invalid toxe, toxp and dtox params");
    }
    // if model.mtrlmod {
    //     epsrox = 3.9;
    //     toxe = model.eot;
    //     epssub = EPS0 * model.epsrsub;
    // } else {
    //     epsrox = model.epsrox;
    //     toxe = model.toxe;
    //     epssub = EPSSI;
    // }
    // if vals.mtrlmod == 0 {
    //     if (model.toxeGiven) && (model.toxpGiven) && (model.dtoxGiven) && (model.toxe != (model.toxp + model.dtox)) {
    //         println!("Warning: toxe, toxp and dtox all given and toxe != toxp + dtox; dtox ignored.\n",);
    //     } else if (model.toxeGiven) && (!model.toxpGiven) {
    //         vals.toxp = model.toxe - model.dtox;
    //     } else if (!model.toxeGiven) && (model.toxpGiven) {
    //         {
    //             vals.toxe = model.toxp + model.dtox;
    //         }
    //         if !model.toxmGiven {
    //             vals.toxm = model.toxe;
    //         }
    //     }
    // } else if vals.mtrlcompatmod != 0 {
    //     let T0 = vals.epsrox / 3.9;
    //     if (model.eotGiven) && (model.toxpGiven) && (model.dtoxGiven) && (abs(model.eot * T0 - (model.toxp + model.dtox)) > 1.0e-20) {
    //         println!("Warning: eot, toxp and dtox all given and eot * EPSROX / 3.9 != toxp + dtox; dtox ignored.\n");
    //     } else if (model.eotGiven) && (!model.toxpGiven) {
    //         vals.toxp = T0 * model.eot - model.dtox;
    //     } else if (!model.eotGiven) && (model.toxpGiven) {
    //         vals.eot = (model.toxp + model.dtox) / T0;
    //         if !model.toxmGiven {
    //             vals.toxm = model.eot;
    //         }
    //     }
    // }

    // FIXME: `coxe` is also derived elsewhere, merge when possible
    let coxe = vals.epsrox * EPS0 / vals.toxe;
    vals.cgso = if let Some(val) =specs.cgso { val } else {
        if specs.dlc.is_some() && vals.dlc > 0.0 {
            vals.dlc * coxe - vals.cgsl
        } else {
            0.6 * vals.xj * coxe
        }
    };
    vals.cgdo = if let Some(val) =specs.cgdo { val } else {
        if specs.dlc.is_some() && vals.dlc > 0.0 {
            vals.dlc * coxe - vals.cgdl
        } else {
            0.6 * vals.xj * coxe
        }
    };
    vals.cgbo = if let Some(val) = specs.cgbo { val } else {
        2.0 * vals.dwc * coxe
    };

    // Value Range-Limiting and Related Stern Warnings
    if vals.pbs < 0.1 {
        vals.pbs = 0.1;
        println!("Given pbs is less than 0.1. Pbs is set to 0.1.\n");
    }
    if vals.pbsws < 0.1 {
        vals.pbsws = 0.1;
        println!("Given pbsws is less than 0.1. Pbsws is set to 0.1.\n",);
    }
    if vals.pbswgs < 0.1 {
        vals.pbswgs = 0.1;
        println!("Given pbswgs is less than 0.1. Pbswgs is set to 0.1.\n",);
    }
    if vals.pbd < 0.1 {
        vals.pbd = 0.1;
        println!("Given pbd is less than 0.1. Pbd is set to 0.1.\n");
    }
    if vals.pbswd < 0.1 {
        vals.pbswd = 0.1;
        println!("Given pbswd is less than 0.1. Pbswd is set to 0.1.\n",);
    }
    if vals.pbswgd < 0.1 {
        vals.pbswgd = 0.1;
        println!("Given pbswgd is less than 0.1. Pbswgd is set to 0.1.\n",);
    }
    if vals.ijthdfwd <= 0.0 {
        vals.ijthdfwd = 0.0;
        println!("Ijthdfwd reset to %g.\n"); //vals.ijthdfwd);
    }
    if vals.ijthsfwd <= 0.0 {
        vals.ijthsfwd = 0.0;
        println!("Ijthsfwd reset to %g.\n"); //vals.ijthsfwd);
    }
    if vals.ijthdrev <= 0.0 {
        vals.ijthdrev = 0.0;
        println!("Ijthdrev reset to %g.\n"); //vals.ijthdrev);
    }
    if vals.ijthsrev <= 0.0 {
        vals.ijthsrev = 0.0;
        println!("Ijthsrev reset to %g.\n"); //vals.ijthsrev);
    }
    if vals.xjbvd <= 0.0 && (vals.diomod == 2 || vals.diomod == 0) {
        vals.xjbvd = 0.0;
        println!("Xjbvd reset to %g.\n"); //vals.xjbvd);
    }
    if vals.xjbvs <= 0.0 && (vals.diomod == 2 || vals.diomod == 0) {
        vals.xjbvs = 0.0;
        println!("Xjbvs reset to %g.\n"); //vals.xjbvs);
    }
    if vals.bvd <= 0.0 {
        vals.bvd = 0.0;
        println!("BVD reset to %g.\n"); //vals.bvd);
    }
    if vals.bvs <= 0.0 {
        vals.bvs = 0.0;
        println!("BVS reset to %g.\n"); //vals.bvs);
    }
    if vals.jtweff < 0.0 {
        vals.jtweff = 0.0;
        println!("TAT width dependence effect is negative. Jtweff is clamped to zero.\n",);
    }
    if vals.cjsws < 0.0 {
        vals.cjsws = 0.0;
        println!("CJSWS is negative. Cjsws is clamped to zero.\n");
    }
    if vals.cjswd < 0.0 {
        vals.cjswd = 0.0;
        println!("CJSWD is negative. Cjswd is clamped to zero.\n");
    }

    return vals;
}
