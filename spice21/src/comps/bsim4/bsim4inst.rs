use super::bsim4defs::*;
use super::*;
use crate::comps::consts::*;
use crate::comps::mos::MosType;

const DELTA: f64 = 1e-9;

/// Derive Bsim4 Internal Instance Parameters from Model and Instance params
fn from(model: &Bsim4ModelVals, model_derived: &Bsim4ModelDerivedParams, inst: &Bsim4InstSpecs) -> Bsim4InternalParams {
    let mut tmp: f64;
    let mut tmp1: f64;
    let mut tmp2: f64;
    let mut tmp3: f64;

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

    let mut Temp: f64;
    let mut Inv_L: f64;
    let mut Inv_W: f64;
    let mut Inv_LW: f64;
    let mut Dw: f64;
    let mut Dl: f64;

    let mut dumPs: f64;
    let mut dumPd: f64;
    let mut dumAs: f64;
    let mut dumAd: f64;
    let mut PowWeffWr: f64;
    let mut DMCGeff: f64;
    let mut DMCIeff: f64;
    let mut DMDGeff: f64;
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
    let mut sceff: f64;
    let mut Wdrn: f64;
    let mut V0: f64;
    let mut lt1: f64;
    let mut ltw: f64;
    let mut Theta0: f64;
    let mut Delt_vth: f64;
    let mut Vth_NarrowW: f64;
    let mut Lpe_Vb: f64;
    let mut Vth: f64;
    let mut n: f64;
    let mut n0: f64;
    let mut Vgsteff: f64;
    let mut Vgs_eff: f64;
    let mut toxpf: f64;
    let mut toxpi: f64;
    let mut Tcen: f64;

    let mut epsrox: f64;
    let mut vddeot: f64;
    let mut vtfbphi2eot: f64;
    let mut phieot: f64;
    let mut Vtm0eot: f64;
    let mut Vtmeot: f64;
    let mut vbieot: f64;
    let mut Size_Not_Found: bool;
    let mut i: usize;

    let Temp = 300.15; // FIXME !ckt->CKTtemp;
    let delTemp = Temp - model.tnom;

    // Start with blank sets of parameters
    let mut size_params = Bsim4SizeDepParams::default();
    let mut intp = Bsim4InternalParams::default();

    // Instance default values
    intp.l = if let Some(val) = inst.l { val } else { 5.0e-6 };
    intp.w = if let Some(val) = inst.w { val } else { 5.0e-6 };
    intp.nf = if let Some(val) = inst.nf { val } else { 1.0 };
    intp.sa = if let Some(val) = inst.sa { val } else { 0.0 };
    intp.sb = if let Some(val) = inst.sb { val } else { 0.0 };
    intp.sd = if let Some(val) = inst.sd { val } else { 2.0 * model.dmcg };
    intp.sc = if let Some(val) = inst.sc { val } else { 0.0 };
    intp.ad = if let Some(val) = inst.ad { val } else { 0.0 };
    intp.r#as = if let Some(val) = inst.r#as { val } else { 0.0 }; // Note renamed from the keyword `as`
    intp.pd = if let Some(val) = inst.pd { val } else { 0.0 };
    intp.ps = if let Some(val) = inst.ps { val } else { 0.0 };
    intp.nrd = if let Some(val) = inst.nrd { val } else { 1.0 };
    intp.nrs = if let Some(val) = inst.nrs { val } else { 1.0 };
    intp.delvto = if let Some(val) = inst.delvto { val } else { 0.0 };
    // Modal instance params
    // FIXME: check ranges, or enum-ize
    intp.min = if let Some(val) = inst.min { val } else { 0 };
    intp.rgeomod = if let Some(val) = inst.rgeomod { val } else { 0 };

    // Also model parameters
    intp.rbdb = if let Some(val) = inst.rbdb { val } else { model.rbdb };
    intp.rbsb = if let Some(val) = inst.rbsb { val } else { model.rbsb };
    intp.rbpb = if let Some(val) = inst.rbpb { val } else { model.rbpb };
    intp.rbps = if let Some(val) = inst.rbps { val } else { model.rbps };
    intp.rbpd = if let Some(val) = inst.rbpd { val } else { model.rbpd };
    intp.xgw = if let Some(val) = inst.xgw { val } else { model.xgw };
    intp.ngcon = if let Some(val) = inst.ngcon { val } else { model.ngcon };

    // More modes
    // intp.trnqsmod = if let Some(val) = inst.trnqsmod { val } else { model.trnqsmod };
    // intp.acnqsmod = if let Some(val) = inst.acnqsmod { val } else { model.acnqsmod };
    // intp.rbodymod = if let Some(val) = inst.rbodymod { val } else { model.rbodymod };
    // intp.rgatemod = if let Some(val) = inst.rgatemod { val } else { model.rgatemod };
    // intp.geomod = if let Some(val) = inst.geomod { val } else { model.geomod };

    // FIXME: set up a hash table of these size params, cache it, etc.
    Size_Not_Found = true;
    // pSizeDependParamKnot = model->pSizeDependParamKnot;
    //   while ((pSizeDependParamKnot != NULL) && Size_Not_Found)
    //   {   if ((inst.l == pSizeDependParamKnot->Length)
    //       && (inst.w == pSizeDependParamKnot->Width)
    //       && (inst.nf == pSizeDependParamKnot->NFinger))
    //           {   Size_Not_Found = 0;
    //       here->pParam = pSizeDependParamKnot;
    //       pParam = here->pParam; /*bug-fix  */
    //   }
    //   else
    //   {   pLastKnot = pSizeDependParamKnot;
    //       pSizeDependParamKnot = pSizeDependParamKnot->pNext;
    //   }
    //       }

    /* stress effect */
    Ldrn = intp.l;
    Wdrn = intp.w / intp.nf;
    let Lnew = intp.l + model.xl;
    let Wnew = intp.w / intp.nf + model.xw;

    if Size_Not_Found {
        //   pParam = (struct bsim4SizeDependParam *)malloc(
        //                 sizeof(struct bsim4SizeDependParam));
        //       if (pLastKnot == NULL)
        //   {model->pSizeDependParamKnot = pParam;}
        //       else
        //   {pLastKnot->pNext = pParam;}
        //       size_params.pNext = NULL;
        //       here->pParam = pParam;

        size_params.Length = intp.l;
        size_params.Width = intp.w;
        size_params.NFinger = intp.nf;

        T0 = pow(Lnew, model.lln);
        T1 = pow(Wnew, model.lwn);
        tmp1 = model.ll / T0 + model.lw / T1 + model.lwl / (T0 * T1);
        size_params.dl = model.lint + tmp1;
        tmp2 = model.llc / T0 + model.lwc / T1 + model.lwlc / (T0 * T1);
        size_params.dlc = model.dlc + tmp2;

        T2 = pow(Lnew, model.wln);
        T3 = pow(Wnew, model.wwn);
        tmp1 = model.wl / T2 + model.ww / T3 + model.wwl / (T2 * T3);
        size_params.dw = model.wint + tmp1;
        tmp2 = model.wlc / T2 + model.wwc / T3 + model.wwlc / (T2 * T3);
        size_params.dwc = model.dwc + tmp2;
        size_params.dwj = model.dwj + tmp2;

        size_params.leff = Lnew - 2.0 * size_params.dl;
        if size_params.leff <= 0.0 {
            panic!("BSIM4: mosfet %s, model %s: Effective channel length <= 0");
        }

        size_params.weff = Wnew - 2.0 * size_params.dw;
        if size_params.weff <= 0.0 {
            panic!("BSIM4: mosfet %s, model %s: Effective channel width <= 0");
        }

        size_params.leffCV = Lnew - 2.0 * size_params.dlc;
        if size_params.leffCV <= 0.0 {
            panic!("BSIM4: mosfet %s, model %s: Effective channel length for C-V <= 0");
        }

        size_params.weffCV = Wnew - 2.0 * size_params.dwc;
        if size_params.weffCV <= 0.0 {
            panic!("BSIM4: mosfet %s, model %s: Effective channel width for C-V <= 0");
        }

        size_params.weffCJ = Wnew - 2.0 * size_params.dwj;
        if size_params.weffCJ <= 0.0 {
            panic!("BSIM4: mosfet %s, model %s: Effective channel width for S/D junctions <= 0");
        }

        if model.binunit == 1 {
            Inv_L = 1.0e-6 / size_params.leff;
            Inv_W = 1.0e-6 / size_params.weff;
            Inv_LW = 1.0e-12 / (size_params.leff * size_params.weff);
        } else {
            Inv_L = 1.0 / size_params.leff;
            Inv_W = 1.0 / size_params.weff;
            Inv_LW = 1.0 / (size_params.leff * size_params.weff);
        }
        size_params.cdsc = model.cdsc + model.lcdsc * Inv_L + model.wcdsc * Inv_W + model.pcdsc * Inv_LW;
        size_params.cdscb = model.cdscb + model.lcdscb * Inv_L + model.wcdscb * Inv_W + model.pcdscb * Inv_LW;
        size_params.cdscd = model.cdscd + model.lcdscd * Inv_L + model.wcdscd * Inv_W + model.pcdscd * Inv_LW;
        size_params.cit = model.cit + model.lcit * Inv_L + model.wcit * Inv_W + model.pcit * Inv_LW;
        size_params.nfactor = model.nfactor + model.lnfactor * Inv_L + model.wnfactor * Inv_W + model.pnfactor * Inv_LW;
        size_params.tnfactor = model.tnfactor + model.ltnfactor * Inv_L + model.wtnfactor * Inv_W + model.ptnfactor * Inv_LW;
        size_params.xj = model.xj + model.lxj * Inv_L + model.wxj * Inv_W + model.pxj * Inv_LW;
        size_params.vsat = model.vsat + model.lvsat * Inv_L + model.wvsat * Inv_W + model.pvsat * Inv_LW;
        size_params.at = model.at + model.lat * Inv_L + model.wat * Inv_W + model.pat * Inv_LW;
        size_params.a0 = model.a0 + model.la0 * Inv_L + model.wa0 * Inv_W + model.pa0 * Inv_LW;
        size_params.ags = model.ags + model.lags * Inv_L + model.wags * Inv_W + model.pags * Inv_LW;
        size_params.a1 = model.a1 + model.la1 * Inv_L + model.wa1 * Inv_W + model.pa1 * Inv_LW;
        size_params.a2 = model.a2 + model.la2 * Inv_L + model.wa2 * Inv_W + model.pa2 * Inv_LW;
        size_params.keta = model.keta + model.lketa * Inv_L + model.wketa * Inv_W + model.pketa * Inv_LW;
        size_params.nsub = model.nsub + model.lnsub * Inv_L + model.wnsub * Inv_W + model.pnsub * Inv_LW;
        size_params.ndep = model.ndep + model.lndep * Inv_L + model.wndep * Inv_W + model.pndep * Inv_LW;
        size_params.nsd = model.nsd + model.lnsd * Inv_L + model.wnsd * Inv_W + model.pnsd * Inv_LW;
        size_params.phin = model.phin + model.lphin * Inv_L + model.wphin * Inv_W + model.pphin * Inv_LW;
        size_params.ngate = model.ngate + model.lngate * Inv_L + model.wngate * Inv_W + model.pngate * Inv_LW;
        size_params.gamma1 = model.gamma1 + model.lgamma1 * Inv_L + model.wgamma1 * Inv_W + model.pgamma1 * Inv_LW;
        size_params.gamma2 = model.gamma2 + model.lgamma2 * Inv_L + model.wgamma2 * Inv_W + model.pgamma2 * Inv_LW;
        size_params.vbx = model.vbx + model.lvbx * Inv_L + model.wvbx * Inv_W + model.pvbx * Inv_LW;
        size_params.vbm = model.vbm + model.lvbm * Inv_L + model.wvbm * Inv_W + model.pvbm * Inv_LW;
        size_params.xt = model.xt + model.lxt * Inv_L + model.wxt * Inv_W + model.pxt * Inv_LW;
        size_params.vfb = model.vfb + model.lvfb * Inv_L + model.wvfb * Inv_W + model.pvfb * Inv_LW;
        size_params.k1 = model.k1 + model.lk1 * Inv_L + model.wk1 * Inv_W + model.pk1 * Inv_LW;
        size_params.kt1 = model.kt1 + model.lkt1 * Inv_L + model.wkt1 * Inv_W + model.pkt1 * Inv_LW;
        size_params.kt1l = model.kt1l + model.lkt1l * Inv_L + model.wkt1l * Inv_W + model.pkt1l * Inv_LW;
        size_params.k2 = model.k2 + model.lk2 * Inv_L + model.wk2 * Inv_W + model.pk2 * Inv_LW;
        size_params.kt2 = model.kt2 + model.lkt2 * Inv_L + model.wkt2 * Inv_W + model.pkt2 * Inv_LW;
        size_params.k3 = model.k3 + model.lk3 * Inv_L + model.wk3 * Inv_W + model.pk3 * Inv_LW;
        size_params.k3b = model.k3b + model.lk3b * Inv_L + model.wk3b * Inv_W + model.pk3b * Inv_LW;
        size_params.w0 = model.w0 + model.lw0 * Inv_L + model.ww0 * Inv_W + model.pw0 * Inv_LW;
        size_params.lpe0 = model.lpe0 + model.llpe0 * Inv_L + model.wlpe0 * Inv_W + model.plpe0 * Inv_LW;
        size_params.lpeb = model.lpeb + model.llpeb * Inv_L + model.wlpeb * Inv_W + model.plpeb * Inv_LW;
        size_params.dvtp0 = model.dvtp0 + model.ldvtp0 * Inv_L + model.wdvtp0 * Inv_W + model.pdvtp0 * Inv_LW;
        size_params.dvtp1 = model.dvtp1 + model.ldvtp1 * Inv_L + model.wdvtp1 * Inv_W + model.pdvtp1 * Inv_LW;
        size_params.dvtp2 = model.dvtp2 + model.ldvtp2 * Inv_L + model.wdvtp2 * Inv_W + model.pdvtp2 * Inv_LW;
        size_params.dvtp3 = model.dvtp3 + model.ldvtp3 * Inv_L + model.wdvtp3 * Inv_W + model.pdvtp3 * Inv_LW;
        size_params.dvtp4 = model.dvtp4 + model.ldvtp4 * Inv_L + model.wdvtp4 * Inv_W + model.pdvtp4 * Inv_LW;
        size_params.dvtp5 = model.dvtp5 + model.ldvtp5 * Inv_L + model.wdvtp5 * Inv_W + model.pdvtp5 * Inv_LW;
        size_params.dvt0 = model.dvt0 + model.ldvt0 * Inv_L + model.wdvt0 * Inv_W + model.pdvt0 * Inv_LW;
        size_params.dvt1 = model.dvt1 + model.ldvt1 * Inv_L + model.wdvt1 * Inv_W + model.pdvt1 * Inv_LW;
        size_params.dvt2 = model.dvt2 + model.ldvt2 * Inv_L + model.wdvt2 * Inv_W + model.pdvt2 * Inv_LW;
        size_params.dvt0w = model.dvt0w + model.ldvt0w * Inv_L + model.wdvt0w * Inv_W + model.pdvt0w * Inv_LW;
        size_params.dvt1w = model.dvt1w + model.ldvt1w * Inv_L + model.wdvt1w * Inv_W + model.pdvt1w * Inv_LW;
        size_params.dvt2w = model.dvt2w + model.ldvt2w * Inv_L + model.wdvt2w * Inv_W + model.pdvt2w * Inv_LW;
        size_params.drout = model.drout + model.ldrout * Inv_L + model.wdrout * Inv_W + model.pdrout * Inv_LW;
        size_params.dsub = model.dsub + model.ldsub * Inv_L + model.wdsub * Inv_W + model.pdsub * Inv_LW;
        size_params.vth0 = model.vth0 + model.lvth0 * Inv_L + model.wvth0 * Inv_W + model.pvth0 * Inv_LW;
        size_params.ua = model.ua + model.lua * Inv_L + model.wua * Inv_W + model.pua * Inv_LW;
        size_params.ua1 = model.ua1 + model.lua1 * Inv_L + model.wua1 * Inv_W + model.pua1 * Inv_LW;
        size_params.ub = model.ub + model.lub * Inv_L + model.wub * Inv_W + model.r#pub * Inv_LW;
        size_params.ub1 = model.ub1 + model.lub1 * Inv_L + model.wub1 * Inv_W + model.pub1 * Inv_LW;
        size_params.uc = model.uc + model.luc * Inv_L + model.wuc * Inv_W + model.puc * Inv_LW;
        size_params.uc1 = model.uc1 + model.luc1 * Inv_L + model.wuc1 * Inv_W + model.puc1 * Inv_LW;
        size_params.ud = model.ud + model.lud * Inv_L + model.wud * Inv_W + model.pud * Inv_LW;
        size_params.ud1 = model.ud1 + model.lud1 * Inv_L + model.wud1 * Inv_W + model.pud1 * Inv_LW;
        size_params.up = model.up + model.lup * Inv_L + model.wup * Inv_W + model.pup * Inv_LW;
        size_params.lp = model.lp + model.llp * Inv_L + model.wlp * Inv_W + model.plp * Inv_LW;
        size_params.eu = model.eu + model.leu * Inv_L + model.weu * Inv_W + model.peu * Inv_LW;
        size_params.u0 = model.u0 + model.lu0 * Inv_L + model.wu0 * Inv_W + model.pu0 * Inv_LW;
        size_params.ute = model.ute + model.lute * Inv_L + model.wute * Inv_W + model.pute * Inv_LW;
        size_params.ucs = model.ucs + model.lucs * Inv_L + model.wucs * Inv_W + model.pucs * Inv_LW;
        size_params.ucste = model.ucste + model.lucste * Inv_L + model.wucste * Inv_W + model.pucste * Inv_LW;
        size_params.voff = model.voff + model.lvoff * Inv_L + model.wvoff * Inv_W + model.pvoff * Inv_LW;
        size_params.tvoff = model.tvoff + model.ltvoff * Inv_L + model.wtvoff * Inv_W + model.ptvoff * Inv_LW;
        size_params.minv = model.minv + model.lminv * Inv_L + model.wminv * Inv_W + model.pminv * Inv_LW;
        size_params.minvcv = model.minvcv + model.lminvcv * Inv_L + model.wminvcv * Inv_W + model.pminvcv * Inv_LW;
        size_params.fprout = model.fprout + model.lfprout * Inv_L + model.wfprout * Inv_W + model.pfprout * Inv_LW;
        size_params.pdits = model.pdits + model.lpdits * Inv_L + model.wpdits * Inv_W + model.ppdits * Inv_LW;
        size_params.pditsd = model.pditsd + model.lpditsd * Inv_L + model.wpditsd * Inv_W + model.ppditsd * Inv_LW;
        size_params.delta = model.delta + model.ldelta * Inv_L + model.wdelta * Inv_W + model.pdelta * Inv_LW;
        size_params.rdsw = model.rdsw + model.lrdsw * Inv_L + model.wrdsw * Inv_W + model.prdsw * Inv_LW;
        size_params.rdw = model.rdw + model.lrdw * Inv_L + model.wrdw * Inv_W + model.prdw * Inv_LW;
        size_params.rsw = model.rsw + model.lrsw * Inv_L + model.wrsw * Inv_W + model.prsw * Inv_LW;
        size_params.prwg = model.prwg + model.lprwg * Inv_L + model.wprwg * Inv_W + model.pprwg * Inv_LW;
        size_params.prwb = model.prwb + model.lprwb * Inv_L + model.wprwb * Inv_W + model.pprwb * Inv_LW;
        size_params.prt = model.prt + model.lprt * Inv_L + model.wprt * Inv_W + model.pprt * Inv_LW;
        size_params.eta0 = model.eta0 + model.leta0 * Inv_L + model.weta0 * Inv_W + model.peta0 * Inv_LW;
        size_params.teta0 = model.teta0 + model.lteta0 * Inv_L + model.wteta0 * Inv_W + model.pteta0 * Inv_LW;
        size_params.tvoffcv = model.tvoffcv + model.ltvoffcv * Inv_L + model.wtvoffcv * Inv_W + model.ptvoffcv * Inv_LW;
        size_params.etab = model.etab + model.letab * Inv_L + model.wetab * Inv_W + model.petab * Inv_LW;
        size_params.pclm = model.pclm + model.lpclm * Inv_L + model.wpclm * Inv_W + model.ppclm * Inv_LW;
        size_params.pdibl1 = model.pdiblc1 + model.lpdiblc1 * Inv_L + model.wpdiblc1 * Inv_W + model.ppdiblc1 * Inv_LW;
        size_params.pdibl2 = model.pdiblc2 + model.lpdiblc2 * Inv_L + model.wpdiblc2 * Inv_W + model.ppdiblc2 * Inv_LW;
        size_params.pdiblb = model.pdiblcb + model.lpdiblcb * Inv_L + model.wpdiblcb * Inv_W + model.ppdiblcb * Inv_LW;
        size_params.pscbe1 = model.pscbe1 + model.lpscbe1 * Inv_L + model.wpscbe1 * Inv_W + model.ppscbe1 * Inv_LW;
        size_params.pscbe2 = model.pscbe2 + model.lpscbe2 * Inv_L + model.wpscbe2 * Inv_W + model.ppscbe2 * Inv_LW;
        size_params.pvag = model.pvag + model.lpvag * Inv_L + model.wpvag * Inv_W + model.ppvag * Inv_LW;
        size_params.wr = model.wr + model.lwr * Inv_L + model.wwr * Inv_W + model.pwr * Inv_LW;
        size_params.dwg = model.dwg + model.ldwg * Inv_L + model.wdwg * Inv_W + model.pdwg * Inv_LW;
        size_params.dwb = model.dwb + model.ldwb * Inv_L + model.wdwb * Inv_W + model.pdwb * Inv_LW;
        size_params.b0 = model.b0 + model.lb0 * Inv_L + model.wb0 * Inv_W + model.pb0 * Inv_LW;
        size_params.b1 = model.b1 + model.lb1 * Inv_L + model.wb1 * Inv_W + model.pb1 * Inv_LW;
        size_params.alpha0 = model.alpha0 + model.lalpha0 * Inv_L + model.walpha0 * Inv_W + model.palpha0 * Inv_LW;
        size_params.alpha1 = model.alpha1 + model.lalpha1 * Inv_L + model.walpha1 * Inv_W + model.palpha1 * Inv_LW;
        size_params.beta0 = model.beta0 + model.lbeta0 * Inv_L + model.wbeta0 * Inv_W + model.pbeta0 * Inv_LW;
        size_params.agidl = model.agidl + model.lagidl * Inv_L + model.wagidl * Inv_W + model.pagidl * Inv_LW;
        size_params.bgidl = model.bgidl + model.lbgidl * Inv_L + model.wbgidl * Inv_W + model.pbgidl * Inv_LW;
        size_params.cgidl = model.cgidl + model.lcgidl * Inv_L + model.wcgidl * Inv_W + model.pcgidl * Inv_LW;
        size_params.egidl = model.egidl + model.legidl * Inv_L + model.wegidl * Inv_W + model.pegidl * Inv_LW;
        size_params.rgidl = model.rgidl + model.lrgidl * Inv_L + model.wrgidl * Inv_W + model.prgidl * Inv_LW;
        size_params.kgidl = model.kgidl + model.lkgidl * Inv_L + model.wkgidl * Inv_W + model.pkgidl * Inv_LW;
        size_params.fgidl = model.fgidl + model.lfgidl * Inv_L + model.wfgidl * Inv_W + model.pfgidl * Inv_LW;
        size_params.agisl = model.agisl + model.lagisl * Inv_L + model.wagisl * Inv_W + model.pagisl * Inv_LW;
        size_params.bgisl = model.bgisl + model.lbgisl * Inv_L + model.wbgisl * Inv_W + model.pbgisl * Inv_LW;
        size_params.cgisl = model.cgisl + model.lcgisl * Inv_L + model.wcgisl * Inv_W + model.pcgisl * Inv_LW;
        size_params.egisl = model.egisl + model.legisl * Inv_L + model.wegisl * Inv_W + model.pegisl * Inv_LW;
        size_params.rgisl = model.rgisl + model.lrgisl * Inv_L + model.wrgisl * Inv_W + model.prgisl * Inv_LW;
        size_params.kgisl = model.kgisl + model.lkgisl * Inv_L + model.wkgisl * Inv_W + model.pkgisl * Inv_LW;
        size_params.fgisl = model.fgisl + model.lfgisl * Inv_L + model.wfgisl * Inv_W + model.pfgisl * Inv_LW;
        size_params.aigc = model.aigc + model.laigc * Inv_L + model.waigc * Inv_W + model.paigc * Inv_LW;
        size_params.bigc = model.bigc + model.lbigc * Inv_L + model.wbigc * Inv_W + model.pbigc * Inv_LW;
        size_params.cigc = model.cigc + model.lcigc * Inv_L + model.wcigc * Inv_W + model.pcigc * Inv_LW;
        size_params.aigsd = model.aigsd + model.laigsd * Inv_L + model.waigsd * Inv_W + model.paigsd * Inv_LW;
        size_params.bigsd = model.bigsd + model.lbigsd * Inv_L + model.wbigsd * Inv_W + model.pbigsd * Inv_LW;
        size_params.cigsd = model.cigsd + model.lcigsd * Inv_L + model.wcigsd * Inv_W + model.pcigsd * Inv_LW;
        size_params.aigs = model.aigs + model.laigs * Inv_L + model.waigs * Inv_W + model.paigs * Inv_LW;
        size_params.bigs = model.bigs + model.lbigs * Inv_L + model.wbigs * Inv_W + model.pbigs * Inv_LW;
        size_params.cigs = model.cigs + model.lcigs * Inv_L + model.wcigs * Inv_W + model.pcigs * Inv_LW;
        size_params.aigd = model.aigd + model.laigd * Inv_L + model.waigd * Inv_W + model.paigd * Inv_LW;
        size_params.bigd = model.bigd + model.lbigd * Inv_L + model.wbigd * Inv_W + model.pbigd * Inv_LW;
        size_params.cigd = model.cigd + model.lcigd * Inv_L + model.wcigd * Inv_W + model.pcigd * Inv_LW;
        size_params.aigbacc = model.aigbacc + model.laigbacc * Inv_L + model.waigbacc * Inv_W + model.paigbacc * Inv_LW;
        size_params.bigbacc = model.bigbacc + model.lbigbacc * Inv_L + model.wbigbacc * Inv_W + model.pbigbacc * Inv_LW;
        size_params.cigbacc = model.cigbacc + model.lcigbacc * Inv_L + model.wcigbacc * Inv_W + model.pcigbacc * Inv_LW;
        size_params.aigbinv = model.aigbinv + model.laigbinv * Inv_L + model.waigbinv * Inv_W + model.paigbinv * Inv_LW;
        size_params.bigbinv = model.bigbinv + model.lbigbinv * Inv_L + model.wbigbinv * Inv_W + model.pbigbinv * Inv_LW;
        size_params.cigbinv = model.cigbinv + model.lcigbinv * Inv_L + model.wcigbinv * Inv_W + model.pcigbinv * Inv_LW;
        size_params.nigc = model.nigc + model.lnigc * Inv_L + model.wnigc * Inv_W + model.pnigc * Inv_LW;
        size_params.nigbacc = model.nigbacc + model.lnigbacc * Inv_L + model.wnigbacc * Inv_W + model.pnigbacc * Inv_LW;
        size_params.nigbinv = model.nigbinv + model.lnigbinv * Inv_L + model.wnigbinv * Inv_W + model.pnigbinv * Inv_LW;
        size_params.ntox = model.ntox + model.lntox * Inv_L + model.wntox * Inv_W + model.pntox * Inv_LW;
        size_params.eigbinv = model.eigbinv + model.leigbinv * Inv_L + model.weigbinv * Inv_W + model.peigbinv * Inv_LW;
        size_params.pigcd = model.pigcd + model.lpigcd * Inv_L + model.wpigcd * Inv_W + model.ppigcd * Inv_LW;
        size_params.poxedge = model.poxedge + model.lpoxedge * Inv_L + model.wpoxedge * Inv_W + model.ppoxedge * Inv_LW;
        size_params.xrcrg1 = model.xrcrg1 + model.lxrcrg1 * Inv_L + model.wxrcrg1 * Inv_W + model.pxrcrg1 * Inv_LW;
        size_params.xrcrg2 = model.xrcrg2 + model.lxrcrg2 * Inv_L + model.wxrcrg2 * Inv_W + model.pxrcrg2 * Inv_LW;
        size_params.lambda = model.lambda + model.llambda * Inv_L + model.wlambda * Inv_W + model.plambda * Inv_LW;
        size_params.vtl = model.vtl + model.lvtl * Inv_L + model.wvtl * Inv_W + model.pvtl * Inv_LW;
        size_params.xn = model.xn + model.lxn * Inv_L + model.wxn * Inv_W + model.pxn * Inv_LW;
        size_params.vfbsdoff = model.vfbsdoff + model.lvfbsdoff * Inv_L + model.wvfbsdoff * Inv_W + model.pvfbsdoff * Inv_LW;
        size_params.tvfbsdoff = model.tvfbsdoff + model.ltvfbsdoff * Inv_L + model.wtvfbsdoff * Inv_W + model.ptvfbsdoff * Inv_LW;

        size_params.cgsl = model.cgsl + model.lcgsl * Inv_L + model.wcgsl * Inv_W + model.pcgsl * Inv_LW;
        size_params.cgdl = model.cgdl + model.lcgdl * Inv_L + model.wcgdl * Inv_W + model.pcgdl * Inv_LW;
        size_params.ckappas = model.ckappas + model.lckappas * Inv_L + model.wckappas * Inv_W + model.pckappas * Inv_LW;
        size_params.ckappad = model.ckappad + model.lckappad * Inv_L + model.wckappad * Inv_W + model.pckappad * Inv_LW;
        size_params.cf = model.cf + model.lcf * Inv_L + model.wcf * Inv_W + model.pcf * Inv_LW;
        size_params.clc = model.clc + model.lclc * Inv_L + model.wclc * Inv_W + model.pclc * Inv_LW;
        size_params.cle = model.cle + model.lcle * Inv_L + model.wcle * Inv_W + model.pcle * Inv_LW;
        size_params.vfbcv = model.vfbcv + model.lvfbcv * Inv_L + model.wvfbcv * Inv_W + model.pvfbcv * Inv_LW;
        size_params.acde = model.acde + model.lacde * Inv_L + model.wacde * Inv_W + model.pacde * Inv_LW;
        size_params.moin = model.moin + model.lmoin * Inv_L + model.wmoin * Inv_W + model.pmoin * Inv_LW;
        size_params.noff = model.noff + model.lnoff * Inv_L + model.wnoff * Inv_W + model.pnoff * Inv_LW;
        size_params.voffcv = model.voffcv + model.lvoffcv * Inv_L + model.wvoffcv * Inv_W + model.pvoffcv * Inv_LW;
        size_params.kvth0we = model.kvth0we + model.lkvth0we * Inv_L + model.wkvth0we * Inv_W + model.pkvth0we * Inv_LW;
        size_params.k2we = model.k2we + model.lk2we * Inv_L + model.wk2we * Inv_W + model.pk2we * Inv_LW;
        size_params.ku0we = model.ku0we + model.lku0we * Inv_L + model.wku0we * Inv_W + model.pku0we * Inv_LW;
        size_params.abulkCVfactor = 1.0 + pow((size_params.clc / size_params.leffCV), size_params.cle);

        T0 = model_derived.TempRatio - 1.0;

        PowWeffWr = pow(size_params.weffCJ * 1.0e6, size_params.wr) * intp.nf;

        T1 = 0.0;
        T2 = 0.0;
        T3 = 0.0;
        T4 = 0.0;
        size_params.ucs = size_params.ucs * pow(model_derived.TempRatio, size_params.ucste);
        if model.tempmod == 0 {
            size_params.ua = size_params.ua + size_params.ua1 * T0;
            size_params.ub = size_params.ub + size_params.ub1 * T0;
            size_params.uc = size_params.uc + size_params.uc1 * T0;
            size_params.ud = size_params.ud + size_params.ud1 * T0;
            size_params.vsattemp = size_params.vsat - size_params.at * T0;
            T10 = size_params.prt * T0;
            if model.rdsmod != 0 {
                /* External Rd(V) */
                T1 = size_params.rdw + T10;
                T2 = model.rdwmin + T10;
                /* External Rs(V) */
                T3 = size_params.rsw + T10;
                T4 = model.rswmin + T10;
            }
            /* Internal Rds(V) in IV */
            size_params.rds0 = (size_params.rdsw + T10) * intp.nf / PowWeffWr;
            size_params.rdswmin = (model.rdswmin + T10) * intp.nf / PowWeffWr;
        } else {
            if model.tempmod == 3 {
                size_params.ua = size_params.ua * pow(model_derived.TempRatio, size_params.ua1);
                size_params.ub = size_params.ub * pow(model_derived.TempRatio, size_params.ub1);
                size_params.uc = size_params.uc * pow(model_derived.TempRatio, size_params.uc1);
                size_params.ud = size_params.ud * pow(model_derived.TempRatio, size_params.ud1);
            } else {
                /* tempmod = 1, 2 */
                size_params.ua = size_params.ua * (1.0 + size_params.ua1 * delTemp);
                size_params.ub = size_params.ub * (1.0 + size_params.ub1 * delTemp);
                size_params.uc = size_params.uc * (1.0 + size_params.uc1 * delTemp);
                size_params.ud = size_params.ud * (1.0 + size_params.ud1 * delTemp);
            }
            size_params.vsattemp = size_params.vsat * (1.0 - size_params.at * delTemp);
            T10 = 1.0 + size_params.prt * delTemp;
            if model.rdsmod != 0 {
                /* External Rd(V) */
                T1 = size_params.rdw * T10;
                T2 = model.rdwmin * T10;

                /* External Rs(V) */
                T3 = size_params.rsw * T10;
                T4 = model.rswmin * T10;
            }
            /* Internal Rds(V) in IV */
            size_params.rds0 = size_params.rdsw * T10 * intp.nf / PowWeffWr;
            size_params.rdswmin = model.rdswmin * T10 * intp.nf / PowWeffWr;
        }

        if T1 < 0.0 {
            T1 = 0.0;
            println!("Warning: Rdw at current temperature is negative; set to 0.\n");
        }
        if T2 < 0.0 {
            T2 = 0.0;
            println!("Warning: Rdwmin at current temperature is negative; set to 0.\n");
        }
        size_params.rd0 = T1 / PowWeffWr;
        size_params.rdwmin = T2 / PowWeffWr;
        if T3 < 0.0 {
            T3 = 0.0;
            println!("Warning: Rsw at current temperature is negative; set to 0.\n");
        }
        if T4 < 0.0 {
            T4 = 0.0;
            println!("Warning: Rswmin at current temperature is negative; set to 0.\n");
        }
        size_params.rs0 = T3 / PowWeffWr;
        size_params.rswmin = T4 / PowWeffWr;

        if size_params.u0 > 1.0 {
            size_params.u0 = size_params.u0 / 1.0e4;
        }

        /* mobility channel length dependence */
        T5 = 1.0 - size_params.up * exp(-size_params.leff / size_params.lp);
        size_params.u0temp = size_params.u0 * T5 * pow(model_derived.TempRatio, size_params.ute);
        if size_params.eu < 0.0 {
            size_params.eu = 0.0;
            println!("Warning: eu has been negative; reset to 0.0.\n");
        }
        if size_params.ucs < 0.0 {
            size_params.ucs = 0.0;
            println!("Warning: ucs has been negative; reset to 0.0.\n");
        }

        size_params.vfbsdoff = size_params.vfbsdoff * (1.0 + size_params.tvfbsdoff * delTemp);
        size_params.voff = size_params.voff * (1.0 + size_params.tvoff * delTemp);

        size_params.nfactor = size_params.nfactor + size_params.tnfactor * delTemp / model.tnom;
        size_params.voffcv = size_params.voffcv * (1.0 + size_params.tvoffcv * delTemp);
        size_params.eta0 = size_params.eta0 + size_params.teta0 * delTemp / model.tnom;

        /* Source End Velocity Limit  */
        if ((model.vtlGiven) && (model.vtl > 0.0)) {
            if model.lc < 0.0 {
                size_params.lc = 0.0;
            } else {
                size_params.lc = model.lc;
            }
            T0 = size_params.leff / (size_params.xn * size_params.leff + size_params.lc);
            size_params.tfactor = (1.0 - T0) / (1.0 + T0);
        }

        size_params.cgdo = (model.cgdo + size_params.cf) * size_params.weffCV;
        size_params.cgso = (model.cgso + size_params.cf) * size_params.weffCV;
        size_params.cgbo = model.cgbo * size_params.leffCV * intp.nf;

        if !model.ndepGiven && model.gamma1Given {
            T0 = size_params.gamma1 * model_derived.coxe;
            size_params.ndep = 3.01248e22 * T0 * T0;
        }

        size_params.phi = model_derived.vtm0 * log(size_params.ndep / model_derived.ni) + size_params.phin + 0.4;

        size_params.sqrtPhi = sqrt(size_params.phi);
        size_params.phis3 = size_params.sqrtPhi * size_params.phi;

        size_params.Xdep0 = sqrt(2.0 * model_derived.epssub / (Q * size_params.ndep * 1.0e6)) * size_params.sqrtPhi;
        size_params.sqrtXdep0 = sqrt(size_params.Xdep0);

        if model.mtrlmod == 0 {
            size_params.litl = sqrt(3.0 * 3.9 / model.epsrox * size_params.xj * model.toxe);
        } else {
            size_params.litl = sqrt(model.epsrsub / model.epsrox * size_params.xj * model.toxe);
        }

        size_params.vbi = model_derived.vtm0 * log(size_params.nsd * size_params.ndep / (model_derived.ni * model_derived.ni));

        if model.mtrlmod == 0 {
            if size_params.ngate > 0.0 {
                size_params.vfbsd = model_derived.vtm0 * log(size_params.ngate / size_params.nsd);
            } else {
                size_params.vfbsd = 0.0;
            }
        } else {
            T0 = model_derived.vtm0 * log(size_params.nsd / model_derived.ni);
            T1 = 0.5 * model_derived.Eg0;
            if T0 > T1 {
                T0 = T1;
            }
            T2 = model.easub + T1 - model.p() * T0;
            size_params.vfbsd = model.phig - T2;
        }

        size_params.cdep0 = sqrt(Q * model_derived.epssub * size_params.ndep * 1.0e6 / 2.0 / size_params.phi);

        size_params.ToxRatio = exp(size_params.ntox * log(model.toxref / model.toxe)) / model.toxe / model.toxe;
        size_params.ToxRatioEdge = exp(size_params.ntox * log(model.toxref / (model.toxe * size_params.poxedge)))
            / model.toxe
            / model.toxe
            / size_params.poxedge
            / size_params.poxedge;
        size_params.Aechvb = match model.mos_type {
            MosType::NMOS => 4.97232e-7,
            MosType::PMOS => 3.42537e-7,
        };
        size_params.Bechvb = if model.p() == 1.0 {
            // FIXME: MOS enum
            7.45669e11
        } else {
            1.16645e12
        };
        size_params.AechvbEdgeS = size_params.Aechvb * size_params.weff * model.dlcig * size_params.ToxRatioEdge;
        size_params.AechvbEdgeD = size_params.Aechvb * size_params.weff * model.dlcigd * size_params.ToxRatioEdge;
        size_params.BechvbEdge = -size_params.Bechvb * model.toxe * size_params.poxedge;
        size_params.Aechvb *= size_params.weff * size_params.leff * size_params.ToxRatio;
        size_params.Bechvb *= -model.toxe;

        size_params.mstar = 0.5 + atan(size_params.minv) / PI;
        size_params.mstarcv = 0.5 + atan(size_params.minvcv) / PI;
        size_params.voffcbn = size_params.voff + model.voffl / size_params.leff;
        size_params.voffcbncv = size_params.voffcv + model.voffcvl / size_params.leff;

        size_params.ldeb = sqrt(model_derived.epssub * model_derived.vtm0 / (Q * size_params.ndep * 1.0e6)) / 3.0;
        size_params.acde *= pow((size_params.ndep / 2.0e16), -0.25);

        if model.k1Given || model.k2Given {
            if !model.k1Given {
                size_params.k1 = 0.53;
                println!("Warning: k1 should be specified with k2.\n");
            }
            if !model.k2Given {
                size_params.k2 = -0.0186;
                println!("Warning: k2 should be specified with k1.\n");
            }
            if model.nsubGiven {
                println!("Warning: nsub is ignored because k1 or k2 is given.\n",);
            }
            if model.xtGiven {
                println!("Warning: xt is ignored because k1 or k2 is given.\n",);
            }
            if model.vbxGiven {
                println!("Warning: vbx is ignored because k1 or k2 is given.\n",);
            }
            if model.gamma1Given {
                println!("Warning: gamma1 is ignored because k1 or k2 is given.\n",);
            }
            if model.gamma2Given {
                println!("Warning: gamma2 is ignored because k1 or k2 is given.\n",);
            }
        } else {
            if !model.vbxGiven {
                size_params.vbx = size_params.phi - 7.7348e-4 * size_params.ndep * size_params.xt * size_params.xt;
            }
            if size_params.vbx > 0.0 {
                size_params.vbx = -size_params.vbx;
            }
            if size_params.vbm > 0.0 {
                size_params.vbm = -size_params.vbm;
            }

            if !model.gamma1Given {
                size_params.gamma1 = 5.753e-12 * sqrt(size_params.ndep) / model_derived.coxe;
            }
            if !model.gamma2Given {
                size_params.gamma2 = 5.753e-12 * sqrt(size_params.nsub) / model_derived.coxe;
            }

            T0 = size_params.gamma1 - size_params.gamma2;
            T1 = sqrt(size_params.phi - size_params.vbx) - size_params.sqrtPhi;
            T2 = sqrt(size_params.phi * (size_params.phi - size_params.vbm)) - size_params.phi;
            size_params.k2 = T0 * T1 / (2.0 * T2 + size_params.vbm);
            size_params.k1 = size_params.gamma2 - 2.0 * size_params.k2 * sqrt(size_params.phi - size_params.vbm);
        }

        if !model.vfbGiven {
            if model.vth0Given {
                size_params.vfb = model.p() * size_params.vth0 - size_params.phi - size_params.k1 * size_params.sqrtPhi;
            } else {
                if ((model.mtrlmod != 0) && (model.phigGiven) && (model.nsubGiven)) {
                    T0 = model_derived.vtm0 * log(size_params.nsub / model_derived.ni);
                    T1 = 0.5 * model_derived.Eg0;
                    if T0 > T1 {
                        T0 = T1;
                    }
                    T2 = model.easub + T1 + model.p() * T0;
                    size_params.vfb = model.phig - T2;
                } else {
                    size_params.vfb = -1.0;
                }
            }
        }
        if !model.vth0Given {
            size_params.vth0 = model.p() * (size_params.vfb + size_params.phi + size_params.k1 * size_params.sqrtPhi);
        }

        size_params.k1ox = size_params.k1 * model.toxe / model.toxm;

        tmp = sqrt(model_derived.epssub / (model.epsrox * EPS0) * model.toxe * size_params.Xdep0);
        T0 = size_params.dsub * size_params.leff / tmp;
        if T0 < EXP_THRESHOLD {
            T1 = exp(T0);
            T2 = T1 - 1.0;
            T3 = T2 * T2;
            T4 = T3 + 2.0 * T1 * MIN_EXP;
            size_params.theta0vb0 = T1 / T4;
        } else {
            size_params.theta0vb0 = 1.0 / (MAX_EXP - 2.0);
        }

        T0 = size_params.drout * size_params.leff / tmp;
        if T0 < EXP_THRESHOLD {
            T1 = exp(T0);
            T2 = T1 - 1.0;
            T3 = T2 * T2;
            T4 = T3 + 2.0 * T1 * MIN_EXP;
            T5 = T1 / T4;
        } else {
            T5 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
        }
        size_params.thetaRout = size_params.pdibl1 * T5 + size_params.pdibl2;

        tmp = sqrt(size_params.Xdep0);
        tmp1 = size_params.vbi - size_params.phi;
        tmp2 = model_derived.factor1 * tmp;

        T0 = size_params.dvt1w * size_params.weff * size_params.leff / tmp2;
        if T0 < EXP_THRESHOLD {
            T1 = exp(T0);
            T2 = T1 - 1.0;
            T3 = T2 * T2;
            T4 = T3 + 2.0 * T1 * MIN_EXP;
            T8 = T1 / T4;
        } else {
            T8 = 1.0 / (MAX_EXP - 2.0);
        }
        T0 = size_params.dvt0w * T8;
        T8 = T0 * tmp1;

        T0 = size_params.dvt1 * size_params.leff / tmp2;
        if T0 < EXP_THRESHOLD {
            T1 = exp(T0);
            T2 = T1 - 1.0;
            T3 = T2 * T2;
            T4 = T3 + 2.0 * T1 * MIN_EXP;
            T9 = T1 / T4;
        } else {
            T9 = 1.0 / (MAX_EXP - 2.0);
        }
        T9 = size_params.dvt0 * T9 * tmp1;

        T4 = model.toxe * size_params.phi / (size_params.weff + size_params.w0);

        T0 = sqrt(1.0 + size_params.lpe0 / size_params.leff);
        if ((model.tempmod == 1) || (model.tempmod == 0)) {
            T3 = (size_params.kt1 + size_params.kt1l / size_params.leff) * (model_derived.TempRatio - 1.0);
        }
        if ((model.tempmod == 2) || (model.tempmod == 3)) {
            T3 = -size_params.kt1 * (model_derived.TempRatio - 1.0);
        }

        T5 = size_params.k1ox * (T0 - 1.0) * size_params.sqrtPhi + T3;
        size_params.vfbzbfactor = -T8 - T9 + size_params.k3 * T4 + T5 - size_params.phi - size_params.k1 * size_params.sqrtPhi;

        /* stress effect */
        T0 = pow(Lnew, model.llodku0);
        W_tmp = Wnew + model.wlod;
        T1 = pow(W_tmp, model.wlodku0);
        tmp1 = model.lku0 / T0 + model.wku0 / T1 + model.pku0 / (T0 * T1);
        size_params.ku0 = 1.0 + tmp1;

        T0 = pow(Lnew, model.llodvth);
        T1 = pow(W_tmp, model.wlodvth);
        tmp1 = model.lkvth0 / T0 + model.wkvth0 / T1 + model.pkvth0 / (T0 * T1);
        size_params.kvth0 = 1.0 + tmp1;
        size_params.kvth0 = sqrt(size_params.kvth0 * size_params.kvth0 + DELTA);

        T0 = (model_derived.TempRatio - 1.0);
        size_params.ku0temp = size_params.ku0 * (1.0 + model.tku0 * T0) + DELTA;

        Inv_saref = 1.0 / (model.saref + 0.5 * Ldrn);
        Inv_sbref = 1.0 / (model.sbref + 0.5 * Ldrn);
        size_params.inv_od_ref = Inv_saref + Inv_sbref;
        size_params.rho_ref = model.ku0 / size_params.ku0temp * size_params.inv_od_ref;

        /*high k*/
        /*Calculate VgsteffVth for mobMod=3*/
        if model.mobmod == 3 {
            /*Calculate n @ Vbs=Vds=0*/
            lt1 = model_derived.factor1 * size_params.sqrtXdep0;
            T0 = size_params.dvt1 * size_params.leff / lt1;
            if T0 < EXP_THRESHOLD {
                T1 = exp(T0);
                T2 = T1 - 1.0;
                T3 = T2 * T2;
                T4 = T3 + 2.0 * T1 * MIN_EXP;
                Theta0 = T1 / T4;
            } else {
                Theta0 = 1.0 / (MAX_EXP - 2.0);
            }

            tmp1 = model_derived.epssub / size_params.Xdep0;
            tmp2 = size_params.nfactor * tmp1;
            tmp3 = (tmp2 + size_params.cdsc * Theta0 + size_params.cit) / model_derived.coxe;
            if tmp3 >= -0.5 {
                n0 = 1.0 + tmp3;
            } else {
                T0 = 1.0 / (3.0 + 8.0 * tmp3);
                n0 = (1.0 + 3.0 * tmp3) * T0;
            }

            T0 = n0 * model_derived.vtm;
            T1 = size_params.voffcbn;
            T2 = T1 / T0;
            if T2 < -EXP_THRESHOLD {
                T3 = model_derived.coxe * MIN_EXP / size_params.cdep0;
                T4 = size_params.mstar + T3 * n0;
            } else if T2 > EXP_THRESHOLD {
                T3 = model_derived.coxe * MAX_EXP / size_params.cdep0;
                T4 = size_params.mstar + T3 * n0;
            } else {
                T3 = exp(T2) * model_derived.coxe / size_params.cdep0;
                T4 = size_params.mstar + T3 * n0;
            }
            size_params.VgsteffVth = T0 * log(2.0) / T4;
        }

        /* New DITS term added in 4.7 */
        T0 = -size_params.dvtp3 * log(size_params.leff);
        T1 = dexpb(T0);
        size_params.dvtp2factor = size_params.dvtp5 + size_params.dvtp2 * T1;
    } /* End of SizeNotFound */

    /*  stress effect */
    if ((intp.sa > 0.0) && (intp.sb > 0.0) && ((intp.nf == 1.0) || ((intp.nf > 1.0) && (intp.sd > 0.0)))) {
        Inv_sa = 0.0;
        Inv_sb = 0.0;

        kvsat = model.kvsat;
        if model.kvsat < -1.0 {
            kvsat = -1.0;
            println!("Warning: KVSAT = {} is too small; -1.0 is used.\n", model.kvsat,);
        }
        if model.kvsat > 1.0 {
            kvsat = 1.0;
            println!("Warning: KVSAT = {} is too big; 1.0 is used.\n", model.kvsat,);
        }
        let nfi = intp.nf as usize;
        for i in 0..nfi {
            let t0 = 1.0 / intp.nf / (intp.sa + 0.5 * Ldrn + i as f64 * (intp.sd + Ldrn));
            let t1 = 1.0 / intp.nf / (intp.sb + 0.5 * Ldrn + i as f64 * (intp.sd + Ldrn));
            Inv_sa += t0;
            Inv_sb += t1;
        }
        Inv_ODeff = Inv_sa + Inv_sb;
        rho = model.ku0 / size_params.ku0temp * Inv_ODeff;
        T0 = (1.0 + rho) / (1.0 + size_params.rho_ref);
        intp.u0temp = size_params.u0temp * T0;

        T1 = (1.0 + kvsat * rho) / (1.0 + kvsat * size_params.rho_ref);
        intp.vsattemp = size_params.vsattemp * T1;

        OD_offset = Inv_ODeff - size_params.inv_od_ref;
        dvth0_lod = model.kvth0 / size_params.kvth0 * OD_offset;
        dk2_lod = model.stk2 / pow(size_params.kvth0, model.lodk2) * OD_offset;
        deta0_lod = model.steta0 / pow(size_params.kvth0, model.lodeta0) * OD_offset;
        intp.vth0 = size_params.vth0 + dvth0_lod;

        intp.eta0 = size_params.eta0 + deta0_lod;
        intp.k2 = size_params.k2 + dk2_lod;
    } else {
        intp.u0temp = size_params.u0temp;
        intp.vth0 = size_params.vth0;
        intp.vsattemp = size_params.vsattemp;
        intp.eta0 = size_params.eta0;
        intp.k2 = size_params.k2;
    }

    /*  Well Proximity Effect  */
    if model.wpemod != 0 {
        // FIXME: all the WPE parameter effects are done here;
        // scX can become locals
        intp.sca = if let Some(val) = inst.sca { val } else { 0.0 };
        intp.scb = if let Some(val) = inst.scb { val } else { 0.0 };
        intp.scc = if let Some(val) = inst.scc { val } else { 0.0 };

        if ((!inst.sca.is_some()) && (!inst.scb.is_some()) && (!inst.scc.is_some())) {
            if ((inst.sc.is_some()) && (intp.sc > 0.0)) {
                T1 = intp.sc + Wdrn;
                T2 = 1.0 / model.scref;
                intp.sca = model.scref * model.scref / (intp.sc * T1);
                intp.scb = ((0.1 * intp.sc + 0.01 * model.scref) * exp(-10.0 * intp.sc * T2)
                    - (0.1 * T1 + 0.01 * model.scref) * exp(-10.0 * T1 * T2))
                    / Wdrn;
                intp.scc = ((0.05 * intp.sc + 0.0025 * model.scref) * exp(-20.0 * intp.sc * T2)
                    - (0.05 * T1 + 0.0025 * model.scref) * exp(-20.0 * T1 * T2))
                    / Wdrn;
            } else {
                println!("Warning: No WPE as none of SCA, SCB, SCC, SC is given and/or SC not positive.\n");
            }
        }
        if intp.sca < 0.0 {
            println!("Warning: SCA = {} is negative. Set to 0.0.\n", intp.sca);
            intp.sca = 0.0;
        }
        if intp.scb < 0.0 {
            println!("Warning: SCB = {} is negative. Set to 0.0.\n", intp.scb);
            intp.scb = 0.0;
        }
        if intp.scc < 0.0 {
            println!("Warning: SCC = {} is negative. Set to 0.0.\n", intp.scc);
            intp.scc = 0.0;
        }
        if intp.sc < 0.0 {
            println!("Warning: SC = {} is negative. Set to 0.0.\n", intp.sc);
            intp.sc = 0.0;
        }
        sceff = intp.sca + model.web * intp.scb + model.wec * intp.scc;
        intp.vth0 += size_params.kvth0we * sceff;
        intp.k2 += size_params.k2we * sceff;
        T3 = 1.0 + size_params.ku0we * sceff;
        if T3 <= 0.0 {
            T3 = 0.0;
            println!("Warning: ku0we = {} is negatively too high. Negative mobility! \n", size_params.ku0we,);
        }
        intp.u0temp *= T3;
    }

    /* adding delvto  */
    intp.vth0 += intp.delvto;
    intp.vfb = size_params.vfb + model.p() * intp.delvto;

    /* Instance variables calculation  */
    T3 = model.p() * intp.vth0 - intp.vfb - size_params.phi;
    T4 = T3 + T3;
    T5 = 2.5 * T3;
    intp.vtfbphi1 = if model.p() > 0.0 { T4 } else { T5 };
    if intp.vtfbphi1 < 0.0 {
        intp.vtfbphi1 = 0.0;
    }

    intp.vtfbphi2 = 4.0 * T3;
    if intp.vtfbphi2 < 0.0 {
        intp.vtfbphi2 = 0.0;
    }

    if intp.k2 < 0.0 {
        T0 = 0.5 * size_params.k1 / intp.k2;
        intp.vbsc = 0.9 * (size_params.phi - T0 * T0);
        if intp.vbsc > -3.0 {
            intp.vbsc = -3.0;
        } else if intp.vbsc < -30.0 {
            intp.vbsc = -30.0;
        }
    } else {
        intp.vbsc = -30.0;
    }
    if intp.vbsc > size_params.vbm {
        intp.vbsc = size_params.vbm;
    }
    intp.k2ox = intp.k2 * model.toxe / model.toxm;

    intp.vfbzb = size_params.vfbzbfactor + model.p() * intp.vth0;

    // FIXME! whether to include
    // intp.cgso = size_params.cgso;
    // intp.cgdo = size_params.cgdo;

    lnl = log(size_params.leff * 1.0e6);
    lnw = log(size_params.weff * 1.0e6);
    lnnf = log(intp.nf);

    if model.rbodymod == 2 {
        if model.bodymode == 5 {
            rbsbx = model.rbsbx0 * exp(model.rbsdbxl * lnl + model.rbsdbxw * lnw + model.rbsdbxnf * lnnf);
            rbsby = model.rbsby0 * exp(model.rbsdbyl * lnl + model.rbsdbyw * lnw + model.rbsdbynf * lnnf);
            intp.rbsb = rbsbx * rbsby / (rbsbx + rbsby);

            rbdbx = model.rbdbx0 * exp(model.rbsdbxl * lnl + model.rbsdbxw * lnw + model.rbsdbxnf * lnnf);
            rbdby = model.rbdby0 * exp(model.rbsdbyl * lnl + model.rbsdbyw * lnw + model.rbsdbynf * lnnf);

            intp.rbdb = rbdbx * rbdby / (rbdbx + rbdby);
        }

        if ((model.bodymode == 3) || (model.bodymode == 5)) {
            intp.rbps = model.rbps0 * exp(model.rbpsl * lnl + model.rbpsw * lnw + model.rbpsnf * lnnf);
            intp.rbpd = model.rbpd0 * exp(model.rbpdl * lnl + model.rbpdw * lnw + model.rbpdnf * lnnf);
        }
        rbpbx = model.rbpbx0 * exp(model.rbpbxl * lnl + model.rbpbxw * lnw + model.rbpbxnf * lnnf);
        rbpby = model.rbpby0 * exp(model.rbpbyl * lnl + model.rbpbyw * lnw + model.rbpbynf * lnnf);

        intp.rbpb = rbpbx * rbpby / (rbpbx + rbpby);
    }

    if ((model.rbodymod == 1) || ((model.rbodymod == 2) && (model.bodymode == 5))) {
        if intp.rbdb < 1.0e-3 {
            intp.grbdb = 1.0e3; /* in mho */
        } else {
            intp.grbdb = model.gbmin + 1.0 / intp.rbdb;
        }
        if intp.rbpb < 1.0e-3 {
            intp.grbpb = 1.0e3;
        } else {
            intp.grbpb = model.gbmin + 1.0 / intp.rbpb;
        }
        if intp.rbps < 1.0e-3 {
            intp.grbps = 1.0e3;
        } else {
            intp.grbps = model.gbmin + 1.0 / intp.rbps;
        }
        if intp.rbsb < 1.0e-3 {
            intp.grbsb = 1.0e3;
        } else {
            intp.grbsb = model.gbmin + 1.0 / intp.rbsb;
        }
        if intp.rbpd < 1.0e-3 {
            intp.grbpd = 1.0e3;
        } else {
            intp.grbpd = model.gbmin + 1.0 / intp.rbpd;
        }
    }

    if ((model.rbodymod == 2) && (model.bodymode == 3)) {
        intp.grbdb = model.gbmin;
        intp.grbsb = model.gbmin;
        if intp.rbpb < 1.0e-3 {
            intp.grbpb = 1.0e3;
        } else {
            intp.grbpb = model.gbmin + 1.0 / intp.rbpb;
        }
        if intp.rbps < 1.0e-3 {
            intp.grbps = 1.0e3;
        } else {
            intp.grbps = model.gbmin + 1.0 / intp.rbps;
        }
        if intp.rbpd < 1.0e-3 {
            intp.grbpd = 1.0e3;
        } else {
            intp.grbpd = model.gbmin + 1.0 / intp.rbpd;
        }
    }

    if ((model.rbodymod == 2) && (model.bodymode == 1)) {
        intp.grbdb = model.gbmin;
        intp.grbsb = model.gbmin;
        intp.grbps = 1.0e3;
        intp.grbpd = 1.0e3;
        if intp.rbpb < 1.0e-3 {
            intp.grbpb = 1.0e3;
        } else {
            intp.grbpb = model.gbmin + 1.0 / intp.rbpb;
        }
    }

    /*
     * Process geomertry dependent parasitics
     */

    intp.grgeltd = model.rshg * (intp.xgw + size_params.weffCJ / 3.0 / intp.ngcon) / (intp.ngcon * intp.nf * (Lnew - model.xgl));
    if intp.grgeltd > 0.0 {
        intp.grgeltd = 1.0 / intp.grgeltd;
    } else {
        intp.grgeltd = 1.0e3; /* mho */
        if model.rgatemod != 0 {
            println!("Warning: The gate conductance reset to 1.0e3 mho.\n");
        }
    }

    DMCGeff = model.dmcg - model.dmcgt;
    DMCIeff = model.dmci;
    DMDGeff = model.dmdg - model.dmcgt;
    let (ps_calc, pd_calc, as_calc, ad_calc) = BSIM4PAeffGeo(intp.nf, model.geomod, intp.min, size_params.weffCJ, DMCGeff, DMCIeff, DMDGeff);

    intp.Pseff = if let Some(val) = inst.ps {
        if val < 0.0 {
            println!("Warning: Source Perimeter is specified as negative, it is set to zero.\n");
            0.0
        } else if model.permod == 0 {
            val
        } else {
            val - size_params.weffCJ * intp.nf
        }
    } else {
        ps_calc
    };
    if intp.Pseff < 0.0 {
        intp.Pseff = 0.0;
        println!("Warning: Pseff is negative, it is set to zero.\n");
    }
    intp.Pdeff = if let Some(val) = inst.pd {
        if val < 0.0 {
            println!("Warning: Drain Perimeter is specified as negative, it is set to zero.\n");
            0.0
        } else if model.permod == 0 {
            val
        } else {
            val - size_params.weffCJ * intp.nf
        }
    } else {
        pd_calc
    };
    if intp.Pdeff < 0.0 {
        intp.Pdeff = 0.0;
        println!("Warning: Pdeff is negative, it is set to zero.\n");
    }

    intp.Aseff = if let Some(val) = inst.r#as { val } else { as_calc };
    if intp.Aseff < 0.0 {
        intp.Aseff = 0.0;
        println!("Warning: Aseff is negative, it is set to zero.\n");
    }
    intp.Adeff = if let Some(val) = inst.ad { val } else { ad_calc };
    if intp.Adeff < 0.0 {
        intp.Adeff = 0.0;
        println!("Warning: Adeff is negative, it is set to zero.\n");
    }

    // // FIXME: these nodes won't be available here, figure out when & where to do this
    // // probably just do it here per geometry, and each instance will decide whether to use it.
    // /* Processing S/D resistance and conductance below */
    // if intp.sNodePrime != intp.sNode {
    //     intp.sourceConductance = 0.0;
    //     if inst.nrsGiven {
    //         intp.sourceConductance = model.sheetResistance * inst.nrs;
    //     } else if model.rgeomod > 0 {
    //         BSIM4RdseffGeo(
    //             intp.nf,
    //             model.geomod,
    //             model.rgeomod,
    //             intp.min,
    //             size_params.weffCJ,
    //             model.sheetResistance,
    //             DMCGeff,
    //             DMCIeff,
    //             DMDGeff,
    //             1,
    //             &(intp.sourceConductance),
    //         );
    //     } else {
    //         intp.sourceConductance = 0.0;
    //     }
    //     if intp.sourceConductance > 0.0 {
    //         intp.sourceConductance = 1.0 / intp.sourceConductance;
    //     } else {
    //         intp.sourceConductance = 1.0e3; /* mho */
    //         println!("Warning: Source conductance reset to 1.0e3 mho.\n");
    //     }
    // } else {
    //     intp.sourceConductance = 0.0;
    // }
    // if intp.dNodePrime != intp.dNode {
    //     intp.drainConductance = 0.0;
    //     if inst.nrdGiven {
    //         intp.drainConductance = model.sheetResistance * inst.nrd;
    //     } else if model.rgeomod > 0 {
    //         BSIM4RdseffGeo(
    //             intp.nf,
    //             model.geomod,
    //             model.rgeomod,
    //             intp.min,
    //             size_params.weffCJ,
    //             model.sheetResistance,
    //             DMCGeff,
    //             DMCIeff,
    //             DMDGeff,
    //             0,
    //             &(intp.drainConductance),
    //         );
    //     } else {
    //         intp.drainConductance = 0.0;
    //     }
    //     if intp.drainConductance > 0.0 {
    //         intp.drainConductance = 1.0 / intp.drainConductance;
    //     } else {
    //         intp.drainConductance = 1.0e3; /* mho */
    //         println!("Warning: Drain conductance reset to 1.0e3 mho.\n");
    //     }
    // } else {
    //     intp.drainConductance = 0.0;
    // }
    // /* End of Rsd processing */
    let Nvtms = model_derived.vtm * model.njs;
    if ((intp.Aseff <= 0.0) && (intp.Pseff <= 0.0)) {
        SourceSatCurrent = 0.0;
    } else {
        SourceSatCurrent = intp.Aseff * model_derived.SjctTempSatCurDensity
            + intp.Pseff * model_derived.SjctSidewallTempSatCurDensity
            + size_params.weffCJ * intp.nf * model_derived.SjctGateSidewallTempSatCurDensity;
    }
    if SourceSatCurrent > 0.0 {
        match model.diomod {
            0 => {
                if ((model.bvs / Nvtms) > EXP_THRESHOLD) {
                    intp.XExpBVS = model.xjbvs * MIN_EXP;
                } else {
                    intp.XExpBVS = model.xjbvs * exp(-model.bvs / Nvtms);
                }
            }
            1 => {
                intp.vjsmFwd = BSIM4DioIjthVjmEval(Nvtms, model.ijthsfwd, SourceSatCurrent, 0.0);
                intp.IVjsmFwd = SourceSatCurrent * exp(intp.vjsmFwd / Nvtms);
            }
            2 => {
                if ((model.bvs / Nvtms) > EXP_THRESHOLD) {
                    intp.XExpBVS = model.xjbvs * MIN_EXP;
                    tmp = MIN_EXP;
                } else {
                    intp.XExpBVS = exp(-model.bvs / Nvtms);
                    tmp = intp.XExpBVS;
                    intp.XExpBVS *= model.xjbvs;
                }

                intp.vjsmFwd = BSIM4DioIjthVjmEval(Nvtms, model.ijthsfwd, SourceSatCurrent, intp.XExpBVS);
                T0 = exp(intp.vjsmFwd / Nvtms);
                intp.IVjsmFwd = SourceSatCurrent * (T0 - intp.XExpBVS / T0 + intp.XExpBVS - 1.0);
                intp.SslpFwd = SourceSatCurrent * (T0 + intp.XExpBVS / T0) / Nvtms;

                T2 = model.ijthsrev / SourceSatCurrent;
                if T2 < 1.0 {
                    T2 = 10.0;
                    println!("Warning: ijthsrev too small and set to 10 times IsbSat.\n",);
                }
                intp.vjsmRev = -model.bvs - Nvtms * log((T2 - 1.0) / model.xjbvs);
                T1 = model.xjbvs * exp(-(model.bvs + intp.vjsmRev) / Nvtms);
                intp.IVjsmRev = SourceSatCurrent * (1.0 + T1);
                intp.SslpRev = -SourceSatCurrent * T1 / Nvtms;
            }
            _ => println!("Specified dioMod {} not matched\n", model.diomod),
        }
    }

    let Nvtmd = model_derived.vtm * model.njd;
    if ((intp.Adeff <= 0.0) && (intp.Pdeff <= 0.0)) {
        DrainSatCurrent = 0.0;
    } else {
        DrainSatCurrent = intp.Adeff * model_derived.DjctTempSatCurDensity
            + intp.Pdeff * model_derived.DjctSidewallTempSatCurDensity
            + size_params.weffCJ * intp.nf * model_derived.DjctGateSidewallTempSatCurDensity;
    }
    if DrainSatCurrent > 0.0 {
        match model.diomod {
            0 => {
                if ((model.bvd / Nvtmd) > EXP_THRESHOLD) {
                    intp.XExpBVD = model.xjbvd * MIN_EXP;
                } else {
                    intp.XExpBVD = model.xjbvd * exp(-model.bvd / Nvtmd);
                }
            }
            1 => {
                intp.vjdmFwd = BSIM4DioIjthVjmEval(Nvtmd, model.ijthdfwd, DrainSatCurrent, 0.0);
                intp.IVjdmFwd = DrainSatCurrent * exp(intp.vjdmFwd / Nvtmd);
            }
            2 => {
                if ((model.bvd / Nvtmd) > EXP_THRESHOLD) {
                    intp.XExpBVD = model.xjbvd * MIN_EXP;
                    tmp = MIN_EXP;
                } else {
                    intp.XExpBVD = exp(-model.bvd / Nvtmd);
                    tmp = intp.XExpBVD;
                    intp.XExpBVD *= model.xjbvd;
                }

                intp.vjdmFwd = BSIM4DioIjthVjmEval(Nvtmd, model.ijthdfwd, DrainSatCurrent, intp.XExpBVD);
                T0 = exp(intp.vjdmFwd / Nvtmd);
                intp.IVjdmFwd = DrainSatCurrent * (T0 - intp.XExpBVD / T0 + intp.XExpBVD - 1.0);
                intp.DslpFwd = DrainSatCurrent * (T0 + intp.XExpBVD / T0) / Nvtmd;

                T2 = model.ijthdrev / DrainSatCurrent;
                if T2 < 1.0 {
                    T2 = 10.0;
                    println!("Warning: ijthdrev too small and set to 10 times IdbSat.\n",);
                }
                intp.vjdmRev = -model.bvd - Nvtmd * log((T2 - 1.0) / model.xjbvd); /* bugfix */
                T1 = model.xjbvd * exp(-(model.bvd + intp.vjdmRev) / Nvtmd);
                intp.IVjdmRev = DrainSatCurrent * (1.0 + T1);
                intp.DslpRev = -DrainSatCurrent * T1 / Nvtmd;
            }
            _ => println!("Specified dioMod {} not matched\n", model.diomod),
        }
    }

    T0 = (model_derived.TempRatio - 1.0);
    T7 = model_derived.Eg0 / model_derived.vtm * T0;
    T9 = model.xtss * T7;
    T1 = dexpb(T9);
    T9 = model.xtsd * T7;
    T2 = dexpb(T9);
    T9 = model.xtssws * T7;
    T3 = dexpb(T9);
    T9 = model.xtsswd * T7;
    T4 = dexpb(T9);
    T9 = model.xtsswgs * T7;
    T5 = dexpb(T9);
    T9 = model.xtsswgd * T7;
    T6 = dexpb(T9);
    T11 = sqrt(model.jtweff / size_params.weffCJ) + 1.0;

    T10 = size_params.weffCJ * intp.nf;
    intp.SjctTempRevSatCur = T1 * intp.Aseff * model.jtss;
    intp.DjctTempRevSatCur = T2 * intp.Adeff * model.jtsd;
    intp.SswTempRevSatCur = T3 * intp.Pseff * model.jtssws;
    intp.DswTempRevSatCur = T4 * intp.Pdeff * model.jtsswd;
    intp.SswgTempRevSatCur = T5 * T10 * T11 * model.jtsswgs;
    intp.DswgTempRevSatCur = T6 * T10 * T11 * model.jtsswgd;

    if model.mtrlmod != 0 && model.mtrlcompatmod == 0 {
        /* Calculate TOXP from EOT */
        /* Calculate Vgs_eff @ Vgs = VDD with Poly Depletion Effect */
        Vtm0eot = KB_OVER_Q * model.tempeot;
        Vtmeot = Vtm0eot;
        vbieot = Vtm0eot * log(size_params.nsd * size_params.ndep / (model_derived.ni * model_derived.ni));
        phieot = Vtm0eot * log(size_params.ndep / model_derived.ni) + size_params.phin + 0.4;
        tmp2 = intp.vfb + phieot;
        vddeot = model.p() * model.vddeot;
        T0 = model.epsrgate * EPS0;
        if ((size_params.ngate > 1.0e18) && (size_params.ngate < 1.0e25) && (vddeot > tmp2) && (T0 != 0.0)) {
            T1 = 1.0e6 * Q * T0 * size_params.ngate / (model_derived.coxe * model_derived.coxe);
            T8 = vddeot - tmp2;
            T4 = sqrt(1.0 + 2.0 * T8 / T1);
            T2 = 2.0 * T8 / (T4 + 1.0);
            T3 = 0.5 * T2 * T2 / T1;
            T7 = 1.12 - T3 - 0.05;
            T6 = sqrt(T7 * T7 + 0.224);
            T5 = 1.12 - 0.5 * (T7 + T6);
            Vgs_eff = vddeot - T5;
        } else {
            Vgs_eff = vddeot;
        }

        /* Calculate Vth @ Vds=Vbs=0 */

        V0 = vbieot - phieot;
        lt1 = model_derived.factor1 * size_params.sqrtXdep0;
        ltw = lt1;
        T0 = size_params.dvt1 * model.leffeot / lt1;
        if T0 < EXP_THRESHOLD {
            T1 = exp(T0);
            T2 = T1 - 1.0;
            T3 = T2 * T2;
            T4 = T3 + 2.0 * T1 * MIN_EXP;
            Theta0 = T1 / T4;
        } else {
            Theta0 = 1.0 / (MAX_EXP - 2.0);
        }
        Delt_vth = size_params.dvt0 * Theta0 * V0;
        T0 = size_params.dvt1w * model.weffeot * model.leffeot / ltw;
        if T0 < EXP_THRESHOLD {
            T1 = exp(T0);
            T2 = T1 - 1.0;
            T3 = T2 * T2;
            T4 = T3 + 2.0 * T1 * MIN_EXP;
            T5 = T1 / T4;
        } else {
            T5 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
        }
        T2 = size_params.dvt0w * T5 * V0;
        let TempRatioeot = model.tempeot / model.tnom - 1.0;
        T0 = sqrt(1.0 + size_params.lpe0 / model.leffeot);
        T1 = size_params.k1ox * (T0 - 1.0) * sqrt(phieot) + (size_params.kt1 + size_params.kt1l / model.leffeot) * TempRatioeot;
        Vth_NarrowW = model.toxe * phieot / (model.weffeot + size_params.w0);
        Lpe_Vb = sqrt(1.0 + size_params.lpeb / model.leffeot);
        Vth = model.p() * intp.vth0 + (size_params.k1ox - size_params.k1) * sqrt(phieot) * Lpe_Vb - Delt_vth - T2 + size_params.k3 * Vth_NarrowW + T1;

        /* Calculate n */
        tmp1 = model_derived.epssub / size_params.Xdep0;
        tmp2 = size_params.nfactor * tmp1;
        tmp3 = (tmp2 + size_params.cdsc * Theta0 + size_params.cit) / model_derived.coxe;
        if tmp3 >= -0.5 {
            n = 1.0 + tmp3;
        } else {
            T0 = 1.0 / (3.0 + 8.0 * tmp3);
            n = (1.0 + 3.0 * tmp3) * T0;
        }

        /* Vth correction for Pocket implant */
        if size_params.dvtp0 > 0.0 {
            let t_ = model.leffeot + size_params.dvtp0 * 2.0;
            let v_ = if model.tempmod < 2 { Vtmeot } else { Vtm0eot };
            Vth -= n * v_ * log(model.leffeot / t_);
        }
        Vgsteff = Vgs_eff - Vth;
        /* calculating Toxp */
        T3 = model.p() * intp.vth0 - intp.vfb - phieot;
        T4 = T3 + T3;
        T5 = 2.5 * T3;

        vtfbphi2eot = 4.0 * T3;
        if vtfbphi2eot < 0.0 {
            vtfbphi2eot = 0.0;
        }

        let niter = 0;
        toxpf = model.toxe;
        for i in 0..4 {
            toxpi = toxpf;
            tmp2 = 2.0e8 * toxpf;
            T0 = (Vgsteff + vtfbphi2eot) / tmp2;
            T1 = 1.0 + exp(model.bdos * 0.7 * log(T0));
            Tcen = model.ados * 1.9e-9 / T1;
            toxpf = model.toxe - model.epsrox / model.epsrsub * Tcen;
        }
        intp.toxp = toxpf;
        intp.coxp = model.epsrox * EPS0 / model.toxp;
    } else {
        intp.toxp = model.toxp;
        intp.coxp = model_derived.coxp;
    }

    //   if (BSIM4checkModel(model, here, ckt))
    //   {   IFuid namarray[2];
    //       namarray[0] = model.modName;
    //       namarray[1] = intp.name;
    //       (*(SPfrontEnd->IFerror)) (ERR_FATAL, "Fatal error(s) detected during BSIM4.6.0 parameter checking for %s in model %s", namarray);
    //       return(E_BADPARM);
    //   }

    // FIXME: linking these two, or returning both
    // intp.size_params = size_params;
    return intp;
}

fn BSIM4DioIjthVjmEval(Nvtm: f64, Ijth: f64, Isb: f64, XExpBV: f64) -> f64 {
    let Tc = XExpBV;
    let Tb = 1.0 + Ijth / Isb - Tc;
    let EVjmovNv = 0.5 * (Tb + sqrt(Tb * Tb + 4.0 * Tc));
    return Nvtm * log(EVjmovNv);
}

fn BSIM4PAeffGeo(nf: f64, geo: usize, minSD: usize, Weffcj: f64, DMCG: f64, DMCI: f64, DMDG: f64) -> (f64, f64, f64, f64) {
    let (nuIntD, nuEndD, nuIntS, nuEndS) = if geo < 9 {
        BSIM4NumFingerDiff(nf, minSD)
    } else {
        // For geo = 9 and 10, the numbers of S/D diffusions already known
        (0.0, 0.0, 0.0, 0.0)
    };
    let tmp = DMCG + DMCI;
    let p_iso = tmp + tmp + Weffcj;
    let p_sha = DMCG + DMCG;
    let p_mer = DMDG + DMDG;
    let a_iso = tmp * Weffcj;
    let a_sha = DMCG * Weffcj;
    let a_mer = DMDG * Weffcj;

    match geo {
        0 => {
            let Ps = nuEndS * p_iso + nuIntS * p_sha;
            let Pd = nuEndD * p_iso + nuIntD * p_sha;
            let As = nuEndS * a_iso + nuIntS * a_sha;
            let Ad = nuEndD * a_iso + nuIntD * a_sha;

            return (Ps, Pd, As, Ad);
        }
        1 => {
            let Ps = nuEndS * p_iso + nuIntS * p_sha;
            let Pd = (nuEndD + nuIntD) * p_sha;
            let As = nuEndS * a_iso + nuIntS * a_sha;
            let Ad = (nuEndD + nuIntD) * a_sha;
            return (Ps, Pd, As, Ad);
        }
        2 => {
            let Ps = (nuEndS + nuIntS) * p_sha;
            let Pd = nuEndD * p_iso + nuIntD * p_sha;
            let As = (nuEndS + nuIntS) * a_sha;
            let Ad = nuEndD * a_iso + nuIntD * a_sha;
            return (Ps, Pd, As, Ad);
        }
        3 => {
            let Ps = (nuEndS + nuIntS) * p_sha;
            let Pd = (nuEndD + nuIntD) * p_sha;
            let As = (nuEndS + nuIntS) * a_sha;
            let Ad = (nuEndD + nuIntD) * a_sha;
            return (Ps, Pd, As, Ad);
        }
        4 => {
            let Ps = nuEndS * p_iso + nuIntS * p_sha;
            let Pd = nuEndD * p_mer + nuIntD * p_sha;
            let As = nuEndS * a_iso + nuIntS * a_sha;
            let Ad = nuEndD * a_mer + nuIntD * a_sha;
            return (Ps, Pd, As, Ad);
        }
        5 => {
            let Ps = (nuEndS + nuIntS) * p_sha;
            let Pd = nuEndD * p_mer + nuIntD * p_sha;
            let As = (nuEndS + nuIntS) * a_sha;
            let Ad = nuEndD * a_mer + nuIntD * a_sha;
            return (Ps, Pd, As, Ad);
        }
        6 => {
            let Ps = nuEndS * p_mer + nuIntS * p_sha;
            let Pd = nuEndD * p_iso + nuIntD * p_sha;
            let As = nuEndS * a_mer + nuIntS * a_sha;
            let Ad = nuEndD * a_iso + nuIntD * a_sha;
            return (Ps, Pd, As, Ad);
        }
        7 => {
            let Ps = nuEndS * p_mer + nuIntS * p_sha;
            let Pd = (nuEndD + nuIntD) * p_sha;
            let As = nuEndS * a_mer + nuIntS * a_sha;
            let Ad = (nuEndD + nuIntD) * a_sha;
            return (Ps, Pd, As, Ad);
        }
        8 => {
            let Ps = nuEndS * p_mer + nuIntS * p_sha;
            let Pd = nuEndD * p_mer + nuIntD * p_sha;
            let As = nuEndS * a_mer + nuIntS * a_sha;
            let Ad = nuEndD * a_mer + nuIntD * a_sha;
            return (Ps, Pd, As, Ad);
        }
        9 => {
            /* geo = 9 and 10 happen only when nf = even */
            let Ps = p_iso + (nf - 1.0) * p_sha;
            let Pd = nf * p_sha;
            let As = a_iso + (nf - 1.0) * a_sha;
            let Ad = nf * a_sha;
            return (Ps, Pd, As, Ad);
        }
        10 => {
            let Ps = nf * p_sha;
            let Pd = p_iso + (nf - 1.0) * p_sha;
            let As = nf * a_sha;
            let Ad = a_iso + (nf - 1.0) * a_sha;
            return (Ps, Pd, As, Ad);
        }
        _ => {
            println!("Warning: Specified GEO = {} not matched\n", geo);
            return (0.0, 0.0, 0.0, 0.0);
        }
    }
}

fn BSIM4NumFingerDiff(nf: f64, minSD: usize) -> (f64, f64, f64, f64) {
    let NF = nf as usize;
    if (NF % 2) != 0 {
        let nint = 2.0 * MAX((nf - 1.0) / 2.0, 0.0);
        let nend = 1.0;
        return (nint, nend, nint, nend);
    } else if minSD == 1 {
        // minimize # of source
        let nuEndD = 2.0;
        let nuIntD = 2.0 * max((nf / 2.0 - 1.0), 0.0);
        let nuEndS = 0.0;
        let nuIntS = nf;
        return (nuIntD, nuEndD, nuIntS, nuEndS);
    } else {
        let nuEndD = 0.0;
        let nuIntD = nf;
        let nuEndS = 2.0;
        let nuIntS = 2.0 * max((nf / 2.0 - 1.0), 0.0);
        return (nuIntD, nuEndD, nuIntS, nuEndS);
    }
}
