use super::consts::*;
use super::bsim4::*;
use super::bsim4defs::{Bsim4Model};

/// BSIM4 Model
/// Derive internal parameters from specified params
fn derive(model: &Bsim4Model) -> Bsim4ModelDerivedParams {
    let mut tmp: f64;
    let mut tmp1: f64;
    let mut tmp2: f64;
    let mut tmp3: f64;
    let mut Eg: f64;
    let mut Eg0: f64;
    let mut ni: f64;
    let mut epssub: f64;
    let mut T0: f64;
    let mut T1: f64;
    let mut T2: f64;
    let mut T3: f64;
    let mut T4: f64;
    let mut T5: f64;
    let mut T6: f64;
    let mut T7: f64;
    let mut T8: f64;
    let mut T9: f64;
    let mut Lnew: f64;
    let mut Wnew: f64;
    let mut delTemp: f64;
    let mut Temp: f64;
    let mut TRatio: f64;
    let mut Inv_L: f64;
    let mut Inv_W: f64;
    let mut Inv_LW: f64;
    let mut Dw: f64;
    let mut Dl: f64;
    let mut Vtm0: f64;
    let mut Tnom: f64;
    let mut dumPs: f64;
    let mut dumPd: f64;
    let mut dumAs: f64;
    let mut dumAd: f64;
    let mut PowWeffWr: f64;
    let mut DMCGeff: f64;
    let mut DMCIeff: f64;
    let mut DMDGeff: f64;
    let mut Nvtms: f64;
    let mut Nvtmd: f64;
    let mut SourceSatCurrent: f64;
    let mut DrainSatCurrent: f64;
    let mut T10: f64;
    let mut T11: f64;
    let mut Inv_saref: f64;
    let mut Inv_sbref: f64;
    let mut Inv_sa: f64;
    let mut Inv_sb: f64;
    let mut rho: f64;
    let mut Ldrn: f64;
    let mut dvth0_lod: f64;
    let mut W_tmp: f64;
    let mut Inv_ODeff: f64;
    let mut OD_offset: f64;
    let mut dk2_lod: f64;
    let mut deta0_lod: f64;
    let mut lnl: f64;
    let mut lnw: f64;
    let mut lnnf: f64;
    let mut rbpbx: f64;
    let mut rbpby: f64;
    let mut rbsbx: f64;
    let mut rbsby: f64;
    let mut rbdbx: f64;
    let mut rbdby: f64;
    let mut kvsat: f64;
    let mut wlod: f64;
    let mut sceff: f64;
    let mut Wdrn: f64;
    let mut V0: f64;
    let mut lt1: f64;
    let mut ltw: f64;
    let mut Theta0: f64;
    let mut Delt_vth: f64;
    let mut TempRatio: f64;
    let mut Vth_NarrowW: f64;
    let mut Lpe_Vb: f64;
    let mut Vth: f64;
    let mut n: f64;
    let mut n0: f64;
    let mut Vgsteff: f64;
    let mut Vgs_eff: f64;
    let mut niter: f64;
    let mut toxpf: f64;
    let mut toxpi: f64;
    let mut Tcen: f64;
    let mut toxe: f64;
    let mut epsrox: f64;
    let mut vddeot: f64;
    let mut vtfbphi2eot: f64;
    let mut phieot: f64;
    let mut TempRatioeot: f64;
    let mut Vtm0eot: f64;
    let mut Vtmeot: f64;
    let mut vbieot: f64;

    let mut Size_Not_Found: bool;
    let mut i: usize;

    let model_derived = Bsim4ModelDerivedParams::default();
    let Temp = 300.15; // FIXME !ckt->CKTtemp;

    // FIXME: initial references to `derived` here arent the right idea
    if (model_derived.SbulkJctPotential < 0.1) {
        model_derived.SbulkJctPotential = 0.1;
        panic!("Given pbs is less than 0.1. Pbs is set to 0.1.\n");
    }
    if (model_derived.SsidewallJctPotential < 0.1) {
        model_derived.SsidewallJctPotential = 0.1;
        panic!("Given pbsws is less than 0.1. Pbsws is set to 0.1.\n",);
    }
    if (model_derived.SGatesidewallJctPotential < 0.1) {
        model_derived.SGatesidewallJctPotential = 0.1;
        panic!("Given pbswgs is less than 0.1. Pbswgs is set to 0.1.\n",);
    }

    if (model_derived.DbulkJctPotential < 0.1) {
        model_derived.DbulkJctPotential = 0.1;
        panic!("Given pbd is less than 0.1. Pbd is set to 0.1.\n");
    }
    if (model_derived.DsidewallJctPotential < 0.1) {
        model_derived.DsidewallJctPotential = 0.1;
        panic!("Given pbswd is less than 0.1. Pbswd is set to 0.1.\n",);
    }
    if (model_derived.DGatesidewallJctPotential < 0.1) {
        model_derived.DGatesidewallJctPotential = 0.1;
        panic!("Given pbswgd is less than 0.1. Pbswgd is set to 0.1.\n",);
    }

    if (model.mtrlmod == 0) {
        if ((model.toxeGiven)
            && (model.toxpGiven)
            && (model.dtoxGiven)
            && (model.toxe != (model.toxp + model.dtox)))
        {
            panic!(
                "Warning: toxe, toxp and dtox all given and toxe != toxp + dtox; dtox ignored.\n",
            );
        } else if ((model.toxeGiven) && (!model.toxpGiven)) {
            model_derived.toxp = model.toxe - model.dtox;
        } else if ((!model.toxeGiven) && (model.toxpGiven)) {
            {
                model_derived.toxe = model.toxp + model.dtox;
            }
            if (!model.toxmGiven)
            /* v4.7 */
            {
                model_derived.toxm = model.toxe;
            }
        }
    } else if (model.mtrlcompatmod != 0)
    /* v4.7 */
    {
        T0 = model.epsrox / 3.9;
        if ((model.eotGiven)
            && (model.toxpGiven)
            && (model.dtoxGiven)
            && (abs(model.eot * T0 - (model.toxp + model.dtox)) > 1.0e-20))
        {
            panic!("Warning: eot, toxp and dtox all given and eot * EPSROX / 3.9 != toxp + dtox; dtox ignored.\n");
        } else if ((model.eotGiven) && (!model.toxpGiven)) {
            model_derived.toxp = T0 * model.eot - model.dtox;
        } else if ((!model.eotGiven) && (model.toxpGiven)) {
            model_derived.eot = (model.toxp + model.dtox) / T0;
            if (!model.toxmGiven) {
                model_derived.toxm = model.eot;
            }
        }
    }

    if (model.mtrlmod) {
        epsrox = 3.9;
        toxe = model.eot;
        epssub = EPS0 * model.epsrsub;
    } else {
        epsrox = model.epsrox;
        toxe = model.toxe;
        epssub = EPSSI;
    }

    model_derived.coxe = epsrox * EPS0 / toxe;
    if (model.mtrlmod == 0 || model.mtrlcompatmod != 0) {
        model_derived.coxp = model.epsrox * EPS0 / model.toxp;
    }

    if (!model.cgdoGiven) {
        if (model.dlcGiven && (model.dlc > 0.0)) {
            model_derived.cgdo = model.dlc * model_derived.coxe - model.cgdl;
        } else {
            model_derived.cgdo = 0.6 * model.xj * model_derived.coxe;
        }
    }
    if (!model.cgsoGiven) {
        if (model.dlcGiven && (model.dlc > 0.0)) {
            model_derived.cgso = model.dlc * model_derived.coxe - model.cgsl;
        } else {
            model_derived.cgso = 0.6 * model.xj * model_derived.coxe;
        }
    }
    if (!model.cgboGiven) {
        model_derived.cgbo = 2.0 * model.dwc * model_derived.coxe;
    }

    Tnom = model.tnom;
    TRatio = Temp / Tnom;

    model_derived.vcrit = VT_REF * log(VT_REF / (SQRT2 * 1.0e-14));
    model_derived.factor1 = sqrt(epssub / (epsrox * EPS0) * toxe);

    Vtm0 = KB_OVER_Q * Tnom;
    model_derived.vtm0 = Vtm0;

    if (model.mtrlmod == 0) {
        Eg0 = 1.16 - 7.02e-4 * Tnom * Tnom / (Tnom + 1108.0);
        ni = 1.45e10 * (Tnom / 300.15) * sqrt(Tnom / 300.15) * exp(21.5565981 - Eg0 / (2.0 * Vtm0));
    } else {
        Eg0 = model.bg0sub - model.tbgasub * Tnom * Tnom / (Tnom + model.tbgbsub);
        T0 = model.bg0sub - model.tbgasub * 90090.0225 / (300.15 + model.tbgbsub);
        ni = model.ni0sub * (Tnom / 300.15) * sqrt(Tnom / 300.15) * exp((T0 - Eg0) / (2.0 * Vtm0));
    }

    model_derived.Eg0 = Eg0;
    model_derived.vtm = KB_OVER_Q * Temp;
    if (model.mtrlmod == 0) {
        Eg = 1.16 - 7.02e-4 * Temp * Temp / (Temp + 1108.0);
    } else {
        Eg = model.bg0sub - model.tbgasub * Temp * Temp / (Temp + model.tbgbsub);
    }
    if (Temp != Tnom) {
        T0 = Eg0 / Vtm0 - Eg / model_derived.vtm;
        T1 = log(Temp / Tnom);
        T2 = T0 + model_derived.SjctTempExponent * T1;
        T3 = exp(T2 / model_derived.SjctEmissionCoeff);
        model_derived.SjctTempSatCurDensity = model_derived.SjctSatCurDensity * T3;
        model_derived.SjctSidewallTempSatCurDensity = model_derived.SjctSidewallSatCurDensity * T3;
        model_derived.SjctGateSidewallTempSatCurDensity =
            model_derived.SjctGateSidewallSatCurDensity * T3;

        T2 = T0 + model_derived.DjctTempExponent * T1;
        T3 = exp(T2 / model_derived.DjctEmissionCoeff);
        model_derived.DjctTempSatCurDensity = model_derived.DjctSatCurDensity * T3;
        model_derived.DjctSidewallTempSatCurDensity = model_derived.DjctSidewallSatCurDensity * T3;
        model_derived.DjctGateSidewallTempSatCurDensity =
            model_derived.DjctGateSidewallSatCurDensity * T3;
    } else {
        model_derived.SjctTempSatCurDensity = model_derived.SjctSatCurDensity;
        model_derived.SjctSidewallTempSatCurDensity = model_derived.SjctSidewallSatCurDensity;
        model_derived.SjctGateSidewallTempSatCurDensity =
            model_derived.SjctGateSidewallSatCurDensity;
        model_derived.DjctTempSatCurDensity = model_derived.DjctSatCurDensity;
        model_derived.DjctSidewallTempSatCurDensity = model_derived.DjctSidewallSatCurDensity;
        model_derived.DjctGateSidewallTempSatCurDensity =
            model_derived.DjctGateSidewallSatCurDensity;
    }

    if (model_derived.SjctTempSatCurDensity < 0.0) {
        model_derived.SjctTempSatCurDensity = 0.0;
    }
    if (model_derived.SjctSidewallTempSatCurDensity < 0.0) {
        model_derived.SjctSidewallTempSatCurDensity = 0.0;
    }
    if (model_derived.SjctGateSidewallTempSatCurDensity < 0.0) {
        model_derived.SjctGateSidewallTempSatCurDensity = 0.0;
    }
    if (model_derived.DjctTempSatCurDensity < 0.0) {
        model_derived.DjctTempSatCurDensity = 0.0;
    }
    if (model_derived.DjctSidewallTempSatCurDensity < 0.0) {
        model_derived.DjctSidewallTempSatCurDensity = 0.0;
    }
    if (model_derived.DjctGateSidewallTempSatCurDensity < 0.0) {
        model_derived.DjctGateSidewallTempSatCurDensity = 0.0;
    }

    /* Temperature dependence of D/B and S/B diode capacitance begins */
    delTemp = Temp - model.tnom;
    T0 = model.tcj * delTemp;
    if (T0 >= -1.0) {
        model_derived.SunitAreaTempJctCap = model_derived.SunitAreaJctCap * (1.0 + T0); 
        model_derived.DunitAreaTempJctCap = model_derived.DunitAreaJctCap * (1.0 + T0);
    } else {
        if (model_derived.SunitAreaJctCap > 0.0) {
            model_derived.SunitAreaTempJctCap = 0.0;
            panic!("Temperature effect has caused cjs to be negative. Cjs is clamped to zero.\n",);
        }
        if (model_derived.DunitAreaJctCap > 0.0) {
            model_derived.DunitAreaTempJctCap = 0.0;
            panic!("Temperature effect has caused cjd to be negative. Cjd is clamped to zero.\n",);
        }
    }
    T0 = model.tcjsw * delTemp;
    if (model_derived.SunitLengthSidewallJctCap < 0.0) {
        model_derived.SunitLengthSidewallJctCap = 0.0;
        panic!("CJSWS is negative. Cjsws is clamped to zero.\n");
    }
    if (model_derived.DunitLengthSidewallJctCap < 0.0) {
        model_derived.DunitLengthSidewallJctCap = 0.0;
        panic!("CJSWD is negative. Cjswd is clamped to zero.\n");
    }
    if (T0 >= -1.0) {
        model_derived.SunitLengthSidewallTempJctCap =
            model_derived.SunitLengthSidewallJctCap * (1.0 + T0);
        model_derived.DunitLengthSidewallTempJctCap =
            model_derived.DunitLengthSidewallJctCap * (1.0 + T0);
    } else {
        if (model_derived.SunitLengthSidewallJctCap > 0.0) {
            model_derived.SunitLengthSidewallTempJctCap = 0.0;
            panic!(
                "Temperature effect has caused cjsws to be negative. Cjsws is clamped to zero.\n",
            );
        }
        if (model_derived.DunitLengthSidewallJctCap > 0.0) {
            model_derived.DunitLengthSidewallTempJctCap = 0.0;
            panic!(
                "Temperature effect has caused cjswd to be negative. Cjswd is clamped to zero.\n",
            );
        }
    }
    T0 = model.tcjswg * delTemp;
    if (T0 >= -1.0) {
        model_derived.SunitLengthGateSidewallTempJctCap =
            model_derived.SunitLengthGateSidewallJctCap * (1.0 + T0);
        model_derived.DunitLengthGateSidewallTempJctCap =
            model_derived.DunitLengthGateSidewallJctCap * (1.0 + T0);
    } else {
        if (model_derived.SunitLengthGateSidewallJctCap > 0.0) {
            model_derived.SunitLengthGateSidewallTempJctCap = 0.0;
            panic!(
                "Temperature effect has caused cjswgs to be negative. Cjswgs is clamped to zero.\n",
            );
        }
        if (model_derived.DunitLengthGateSidewallJctCap > 0.0) {
            model_derived.DunitLengthGateSidewallTempJctCap = 0.0;
            panic!(
                "Temperature effect has caused cjswgd to be negative. Cjswgd is clamped to zero.\n",
            );
        }
    }

    model_derived.PhiBS = model_derived.SbulkJctPotential - model.tpb * delTemp;
    if (model_derived.PhiBS < 0.01) {
        model_derived.PhiBS = 0.01;
        panic!("Temperature effect has caused pbs to be less than 0.01. Pbs is clamped to 0.01.\n",);
    }
    model_derived.PhiBD = model_derived.DbulkJctPotential - model.tpb * delTemp;
    if (model_derived.PhiBD < 0.01) {
        model_derived.PhiBD = 0.01;
        panic!("Temperature effect has caused pbd to be less than 0.01. Pbd is clamped to 0.01.\n",);
    }

    model_derived.PhiBSWS = model_derived.SsidewallJctPotential - model.tpbsw * delTemp;
    if (model_derived.PhiBSWS <= 0.01) {
        model_derived.PhiBSWS = 0.01;
        panic!(
            "Temperature effect has caused pbsws to be less than 0.01. Pbsws is clamped to 0.01.\n",
        );
    }
    model_derived.PhiBSWD = model_derived.DsidewallJctPotential - model.tpbsw * delTemp;
    if (model_derived.PhiBSWD <= 0.01) {
        model_derived.PhiBSWD = 0.01;
        panic!(
            "Temperature effect has caused pbswd to be less than 0.01. Pbswd is clamped to 0.01.\n",
        );
    }

    model_derived.PhiBSWGS = model_derived.SGatesidewallJctPotential - model.tpbswg * delTemp;
    if (model_derived.PhiBSWGS <= 0.01) {
        model_derived.PhiBSWGS = 0.01;
        panic!("Temperature effect has caused pbswgs to be less than 0.01. Pbswgs is clamped to 0.01.\n");
    }
    model_derived.PhiBSWGD = model_derived.DGatesidewallJctPotential - model.tpbswg * delTemp;
    if (model_derived.PhiBSWGD <= 0.01) {
        model_derived.PhiBSWGD = 0.01;
        panic!("Temperature effect has caused pbswgd to be less than 0.01. Pbswgd is clamped to 0.01.\n");
    } /* End of junction capacitance */

    if (model.ijthdfwd <= 0.0) {
        // model_derived.ijthdfwd = 0.0;
        panic!("Ijthdfwd reset to %g.\n"); //model.ijthdfwd);
    }
    if (model.ijthsfwd <= 0.0) {
        // model_derived.ijthsfwd = 0.0;
        panic!("Ijthsfwd reset to %g.\n"); //model.ijthsfwd);
    }
    if (model.ijthdrev <= 0.0) {
        // model_derived.ijthdrev = 0.0;
        panic!("Ijthdrev reset to %g.\n"); //model.ijthdrev);
    }
    if (model.ijthsrev <= 0.0) {
        // model_derived.ijthsrev = 0.0;
        panic!("Ijthsrev reset to %g.\n"); //model.ijthsrev);
    }
    if ((model.xjbvd <= 0.0) && (model.diomod == 2)) {
        // model_derived.xjbvd = 0.0;
        panic!("Xjbvd reset to %g.\n"); //model.xjbvd);
    } else if ((model.xjbvd < 0.0) && (model.diomod == 0)) {
        // model_derived.xjbvd = 0.0;
        panic!("Xjbvd reset to %g.\n"); //model.xjbvd);
    }
    if (model.bvd <= 0.0) {
        // model_derived.bvd = 0.0;
        panic!("BVD reset to %g.\n"); //model.bvd);
    }
    if ((model.xjbvs <= 0.0) && (model.diomod == 2)) {
        // model_derived.xjbvs = 0.0;
        panic!("Xjbvs reset to %g.\n"); //model.xjbvs);
    } else if ((model.xjbvs < 0.0) && (model.diomod == 0)) {
        // model_derived.xjbvs = 0.0;
        panic!("Xjbvs reset to %g.\n"); //model.xjbvs);
    }
    if (model.bvs <= 0.0) {
        // model_derived.bvs = 0.0;
        panic!("BVS reset to %g.\n"); //model.bvs);
    }

    /* GEDL current reverse bias */
    T0 = (TRatio - 1.0);
    model_derived.njtsstemp = model.njts * (1.0 + model.tnjts * T0);
    model_derived.njtsswstemp = model.njtssw * (1.0 + model.tnjtssw * T0);
    model_derived.njtsswgstemp = model.njtsswg * (1.0 + model.tnjtsswg * T0);
    model_derived.njtsdtemp = model.njtsd * (1.0 + model.tnjtsd * T0);
    model_derived.njtsswdtemp = model.njtsswd * (1.0 + model.tnjtsswd * T0);
    model_derived.njtsswgdtemp = model.njtsswgd * (1.0 + model.tnjtsswgd * T0);
    /*IBM TAT*/
    if (model.jtweff < 0.0) {
        // model_derived.jtweff = 0.0;
        panic!("TAT width dependence effect is negative. Jtweff is clamped to zero.\n",);
    }

    return model_derived;
}
