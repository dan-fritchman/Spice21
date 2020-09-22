use super::bsim4defs::Bsim4ModelVals;
use super::*;
use crate::comps::consts::*;

/// BSIM4 Model
/// Derive internal parameters from specified param-values 
fn derive(model: &Bsim4ModelVals) -> Bsim4ModelDerivedParams {
    let mut Eg: f64;
    let mut Eg0: f64;
    let mut ni: f64;
    let mut T0: f64;
    let mut T1: f64;
    let mut T2: f64;
    let mut T3: f64;
    let mut delTemp: f64;
    let mut Temp: f64;
    let mut Vtm0: f64;
    let mut Tnom: f64;

    // Create our blank derived-params struct 
    let mut model_derived = Bsim4ModelDerivedParams::default();

    // This first part is not temperature-dependent,
    // but creates a few params not present in `Bims4Model` (`coxp`, `coxe`, etc.).
    model_derived.epssub = if model.mtrlmod != 0 { EPS0 * model.epsrsub } else { EPSSI };
    let epssub = model_derived.epssub;
    model_derived.coxp = if model.mtrlmod == 0 || model.mtrlcompatmod != 0 {
        model.epsrox * EPS0 / model.toxp
    } else {
        0.0
    };
    // FIXME: `coxe` is also derived elsewhere, merge when possible
    model_derived.coxe = model.epsrox * EPS0 / model.toxe;
    model_derived.vcrit = VT_REF * log(VT_REF / (SQRT2 * 1.0e-14));
    model_derived.factor1 = sqrt(epssub / (model.epsrox * EPS0) * model.toxe);

    // On to temperature dependencies
    let Temp = 300.15; // FIXME !ckt->CKTtemp;
    Tnom = model.tnom;
    model_derived.TempRatio = Temp / Tnom;

    Vtm0 = KB_OVER_Q * Tnom;
    model_derived.vtm0 = Vtm0;

    if model.mtrlmod == 0 {
        Eg0 = 1.16 - 7.02e-4 * Tnom * Tnom / (Tnom + 1108.0);
        ni = 1.45e10 * (Tnom / 300.15) * sqrt(Tnom / 300.15) * exp(21.5565981 - Eg0 / (2.0 * Vtm0));
    } else {
        Eg0 = model.bg0sub - model.tbgasub * Tnom * Tnom / (Tnom + model.tbgbsub);
        T0 = model.bg0sub - model.tbgasub * 90090.0225 / (300.15 + model.tbgbsub);
        ni = model.ni0sub * (Tnom / 300.15) * sqrt(Tnom / 300.15) * exp((T0 - Eg0) / (2.0 * Vtm0));
    }

    model_derived.Eg0 = Eg0;
    model_derived.vtm = KB_OVER_Q * Temp;
    if model.mtrlmod == 0 {
        Eg = 1.16 - 7.02e-4 * Temp * Temp / (Temp + 1108.0);
    } else {
        Eg = model.bg0sub - model.tbgasub * Temp * Temp / (Temp + model.tbgbsub);
    }
    if Temp != Tnom {
        T0 = Eg0 / Vtm0 - Eg / model_derived.vtm;
        T1 = log(Temp / Tnom);
        T2 = T0 + model.xtis * T1;
        T3 = exp(T2 / model.njs);
        model_derived.SjctTempSatCurDensity = model.jss * T3;
        model_derived.SjctSidewallTempSatCurDensity = model.jsws * T3;
        model_derived.SjctGateSidewallTempSatCurDensity = model.jswgs * T3;

        T2 = T0 + model.xtid * T1;
        T3 = exp(T2 / model.njd);
        model_derived.DjctTempSatCurDensity = model.jsd * T3;
        model_derived.DjctSidewallTempSatCurDensity = model.jswd * T3;
        model_derived.DjctGateSidewallTempSatCurDensity = model.jswgd * T3;
    } else {
        model_derived.SjctTempSatCurDensity = model.jss;
        model_derived.SjctSidewallTempSatCurDensity = model.jsws;
        model_derived.SjctGateSidewallTempSatCurDensity = model.jswgs;
        model_derived.DjctTempSatCurDensity = model.jsd;
        model_derived.DjctSidewallTempSatCurDensity = model.jswd;
        model_derived.DjctGateSidewallTempSatCurDensity = model.jswgd;
    }

    if model_derived.SjctTempSatCurDensity < 0.0 {
        model_derived.SjctTempSatCurDensity = 0.0;
    }
    if model_derived.SjctSidewallTempSatCurDensity < 0.0 {
        model_derived.SjctSidewallTempSatCurDensity = 0.0;
    }
    if model_derived.SjctGateSidewallTempSatCurDensity < 0.0 {
        model_derived.SjctGateSidewallTempSatCurDensity = 0.0;
    }
    if model_derived.DjctTempSatCurDensity < 0.0 {
        model_derived.DjctTempSatCurDensity = 0.0;
    }
    if model_derived.DjctSidewallTempSatCurDensity < 0.0 {
        model_derived.DjctSidewallTempSatCurDensity = 0.0;
    }
    if model_derived.DjctGateSidewallTempSatCurDensity < 0.0 {
        model_derived.DjctGateSidewallTempSatCurDensity = 0.0;
    }

    /* Temperature dependence of D/B and S/B diode capacitance begins */
    delTemp = Temp - model.tnom;
    T0 = model.tcj * delTemp;
    if T0 >= -1.0 {
        model_derived.SunitAreaTempJctCap = model.cjs * (1.0 + T0);
        model_derived.DunitAreaTempJctCap = model.cjd * (1.0 + T0);
    } else {
        if model.cjs > 0.0 {
            model_derived.SunitAreaTempJctCap = 0.0;
            println!("Temperature effect has caused cjs to be negative. Cjs is clamped to zero.\n",);
        }
        if model.cjd > 0.0 {
            model_derived.DunitAreaTempJctCap = 0.0;
            println!("Temperature effect has caused cjd to be negative. Cjd is clamped to zero.\n",);
        }
    }
    T0 = model.tcjsw * delTemp;
    if T0 >= -1.0 {
        model_derived.SunitLengthSidewallTempJctCap = model.cjsws * (1.0 + T0);
        model_derived.DunitLengthSidewallTempJctCap = model.cjswd * (1.0 + T0);
    } else {
        if model.cjsws > 0.0 {
            model_derived.SunitLengthSidewallTempJctCap = 0.0;
            println!("Temperature effect has caused cjsws to be negative. Cjsws is clamped to zero.\n",);
        }
        if model.cjswd > 0.0 {
            model_derived.DunitLengthSidewallTempJctCap = 0.0;
            println!("Temperature effect has caused cjswd to be negative. Cjswd is clamped to zero.\n",);
        }
    }
    T0 = model.tcjswg * delTemp;
    if T0 >= -1.0 {
        model_derived.SunitLengthGateSidewallTempJctCap = model.cjswgs * (1.0 + T0);
        model_derived.DunitLengthGateSidewallTempJctCap = model.cjswgd * (1.0 + T0);
    } else {
        if model.cjswgs > 0.0 {
            model_derived.SunitLengthGateSidewallTempJctCap = 0.0;
            println!("Temperature effect has caused cjswgs to be negative. Cjswgs is clamped to zero.\n",);
        }
        if model.cjswgd > 0.0 {
            model_derived.DunitLengthGateSidewallTempJctCap = 0.0;
            println!("Temperature effect has caused cjswgd to be negative. Cjswgd is clamped to zero.\n",);
        }
    }

    model_derived.PhiBS = model.pbs - model.tpb * delTemp;
    if model_derived.PhiBS < 0.01 {
        model_derived.PhiBS = 0.01;
        println!("Temperature effect has caused pbs to be less than 0.01. Pbs is clamped to 0.01.\n",);
    }
    model_derived.PhiBD = model.pbd - model.tpb * delTemp;
    if model_derived.PhiBD < 0.01 {
        model_derived.PhiBD = 0.01;
        println!("Temperature effect has caused pbd to be less than 0.01. Pbd is clamped to 0.01.\n",);
    }

    model_derived.PhiBSWS = model.pbsws - model.tpbsw * delTemp;
    if model_derived.PhiBSWS <= 0.01 {
        model_derived.PhiBSWS = 0.01;
        println!("Temperature effect has caused pbsws to be less than 0.01. Pbsws is clamped to 0.01.\n",);
    }
    model_derived.PhiBSWD = model.pbswd - model.tpbsw * delTemp;
    if model_derived.PhiBSWD <= 0.01 {
        model_derived.PhiBSWD = 0.01;
        println!("Temperature effect has caused pbswd to be less than 0.01. Pbswd is clamped to 0.01.\n",);
    }

    model_derived.PhiBSWGS = model.pbswgs - model.tpbswg * delTemp;
    if model_derived.PhiBSWGS <= 0.01 {
        model_derived.PhiBSWGS = 0.01;
        println!("Temperature effect has caused pbswgs to be less than 0.01. Pbswgs is clamped to 0.01.\n");
    }
    model_derived.PhiBSWGD = model.pbswgd - model.tpbswg * delTemp;
    if model_derived.PhiBSWGD <= 0.01 {
        model_derived.PhiBSWGD = 0.01;
        println!("Temperature effect has caused pbswgd to be less than 0.01. Pbswgd is clamped to 0.01.\n");
    } /* End of junction capacitance */

    /* GEDL current reverse bias */
    T0 = model_derived.TempRatio - 1.0;
    model_derived.njtsstemp = model.njts * (1.0 + model.tnjts * T0);
    model_derived.njtsswstemp = model.njtssw * (1.0 + model.tnjtssw * T0);
    model_derived.njtsswgstemp = model.njtsswg * (1.0 + model.tnjtsswg * T0);
    model_derived.njtsdtemp = model.njtsd * (1.0 + model.tnjtsd * T0);
    model_derived.njtsswdtemp = model.njtsswd * (1.0 + model.tnjtsswd * T0);
    model_derived.njtsswgdtemp = model.njtsswgd * (1.0 + model.tnjtsswgd * T0);

    return model_derived;
}
