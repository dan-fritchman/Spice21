use super::consts::*;
use super::bsim4::*;
use super::bsim4defs::*;

/// Derive Bsim4 Internal Parameters from Model and Instance params
fn from(
    model: &Bsim4Model,
    model_derived: &Bsim4ModelDerivedParams,
    inst: &Bsim4Inst,
) -> Bsim4InternalParams {
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
    let mut bodymode: usize;
    let mut Size_Not_Found: bool;
    let mut i: usize;

    let mut size_params = Bsim4SizeDepParams::default();
    let mut intp = Bsim4InternalParams::default();

    // FIXME: set up a hash table of these size params, cache it etc.
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
    Ldrn = inst.l;
    Wdrn = inst.w / inst.nf;

    if (Size_Not_Found) {
        //   pParam = (struct bsim4SizeDependParam *)malloc(
        //                 sizeof(struct bsim4SizeDependParam));
        //       if (pLastKnot == NULL)
        //   {model->pSizeDependParamKnot = pParam;}
        //       else
        //   {pLastKnot->pNext = pParam;}
        //       size_params.pNext = NULL;
        //       here->pParam = pParam;

        size_params.Length = inst.l;
        size_params.Width = inst.w;
        size_params.NFinger = inst.nf;

        Lnew = inst.l + model.xl;
        Wnew = inst.w / inst.nf + model.xw;

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
        if (size_params.leff <= 0.0) {
            panic!("BSIM4: mosfet %s, model %s: Effective channel length <= 0");
        }

        size_params.weff = Wnew - 2.0 * size_params.dw;
        if (size_params.weff <= 0.0) {
            panic!("BSIM4: mosfet %s, model %s: Effective channel width <= 0");
        }

        size_params.leffCV = Lnew - 2.0 * size_params.dlc;
        if (size_params.leffCV <= 0.0) {
            panic!("BSIM4: mosfet %s, model %s: Effective channel length for C-V <= 0");
        }

        size_params.weffCV = Wnew - 2.0 * size_params.dwc;
        if (size_params.weffCV <= 0.0) {
            panic!("BSIM4: mosfet %s, model %s: Effective channel width for C-V <= 0");
        }

        size_params.weffCJ = Wnew - 2.0 * size_params.dwj;
        if (size_params.weffCJ <= 0.0) {
            panic!("BSIM4: mosfet %s, model %s: Effective channel width for S/D junctions <= 0");
        }

        if (model.binunit == 1) {
            Inv_L = 1.0e-6 / size_params.leff;
            Inv_W = 1.0e-6 / size_params.weff;
            Inv_LW = 1.0e-12 / (size_params.leff * size_params.weff);
        } else {
            Inv_L = 1.0 / size_params.leff;
            Inv_W = 1.0 / size_params.weff;
            Inv_LW = 1.0 / (size_params.leff * size_params.weff);
        }
        size_params.cdsc =
            model.cdsc + model.lcdsc * Inv_L + model.wcdsc * Inv_W + model.pcdsc * Inv_LW;
        size_params.cdscb =
            model.cdscb + model.lcdscb * Inv_L + model.wcdscb * Inv_W + model.pcdscb * Inv_LW;

        size_params.cdscd =
            model.cdscd + model.lcdscd * Inv_L + model.wcdscd * Inv_W + model.pcdscd * Inv_LW;

        size_params.cit = model.cit + model.lcit * Inv_L + model.wcit * Inv_W + model.pcit * Inv_LW;
        size_params.nfactor = model.nfactor
            + model.lnfactor * Inv_L
            + model.wnfactor * Inv_W
            + model.pnfactor * Inv_LW;
        size_params.tnfactor = model.tnfactor 
				       + model.ltnfactor * Inv_L
				       + model.wtnfactor * Inv_W
				       + model.ptnfactor * Inv_LW;
        size_params.xj = model.xj + model.lxj * Inv_L + model.wxj * Inv_W + model.pxj * Inv_LW;
        size_params.vsat =
            model.vsat + model.lvsat * Inv_L + model.wvsat * Inv_W + model.pvsat * Inv_LW;
        size_params.at = model.at + model.lat * Inv_L + model.wat * Inv_W + model.pat * Inv_LW;
        size_params.a0 = model.a0 + model.la0 * Inv_L + model.wa0 * Inv_W + model.pa0 * Inv_LW;

        size_params.ags = model.ags + model.lags * Inv_L + model.wags * Inv_W + model.pags * Inv_LW;

        size_params.a1 = model.a1 + model.la1 * Inv_L + model.wa1 * Inv_W + model.pa1 * Inv_LW;
        size_params.a2 = model.a2 + model.la2 * Inv_L + model.wa2 * Inv_W + model.pa2 * Inv_LW;
        size_params.keta =
            model.keta + model.lketa * Inv_L + model.wketa * Inv_W + model.pketa * Inv_LW;
        size_params.nsub =
            model.nsub + model.lnsub * Inv_L + model.wnsub * Inv_W + model.pnsub * Inv_LW;
        size_params.ndep =
            model.ndep + model.lndep * Inv_L + model.wndep * Inv_W + model.pndep * Inv_LW;
        size_params.nsd = model.nsd + model.lnsd * Inv_L + model.wnsd * Inv_W + model.pnsd * Inv_LW;
        size_params.phin =
            model.phin + model.lphin * Inv_L + model.wphin * Inv_W + model.pphin * Inv_LW;
        size_params.ngate =
            model.ngate + model.lngate * Inv_L + model.wngate * Inv_W + model.pngate * Inv_LW;
        size_params.gamma1 =
            model.gamma1 + model.lgamma1 * Inv_L + model.wgamma1 * Inv_W + model.pgamma1 * Inv_LW;
        size_params.gamma2 =
            model.gamma2 + model.lgamma2 * Inv_L + model.wgamma2 * Inv_W + model.pgamma2 * Inv_LW;
        size_params.vbx = model.vbx + model.lvbx * Inv_L + model.wvbx * Inv_W + model.pvbx * Inv_LW;
        size_params.vbm = model.vbm + model.lvbm * Inv_L + model.wvbm * Inv_W + model.pvbm * Inv_LW;
        size_params.xt = model.xt + model.lxt * Inv_L + model.wxt * Inv_W + model.pxt * Inv_LW;
        size_params.vfb = model.vfb + model.lvfb * Inv_L + model.wvfb * Inv_W + model.pvfb * Inv_LW;
        size_params.k1 = model.k1 + model.lk1 * Inv_L + model.wk1 * Inv_W + model.pk1 * Inv_LW;
        size_params.kt1 = model.kt1 + model.lkt1 * Inv_L + model.wkt1 * Inv_W + model.pkt1 * Inv_LW;
        size_params.kt1l =
            model.kt1l + model.lkt1l * Inv_L + model.wkt1l * Inv_W + model.pkt1l * Inv_LW;
        size_params.k2 = model.k2 + model.lk2 * Inv_L + model.wk2 * Inv_W + model.pk2 * Inv_LW;
        size_params.kt2 = model.kt2 + model.lkt2 * Inv_L + model.wkt2 * Inv_W + model.pkt2 * Inv_LW;
        size_params.k3 = model.k3 + model.lk3 * Inv_L + model.wk3 * Inv_W + model.pk3 * Inv_LW;
        size_params.k3b = model.k3b + model.lk3b * Inv_L + model.wk3b * Inv_W + model.pk3b * Inv_LW;
        size_params.w0 = model.w0 + model.lw0 * Inv_L + model.ww0 * Inv_W + model.pw0 * Inv_LW;
        size_params.lpe0 =
            model.lpe0 + model.llpe0 * Inv_L + model.wlpe0 * Inv_W + model.plpe0 * Inv_LW;
        size_params.lpeb =
            model.lpeb + model.llpeb * Inv_L + model.wlpeb * Inv_W + model.plpeb * Inv_LW;
        size_params.dvtp0 =
            model.dvtp0 + model.ldvtp0 * Inv_L + model.wdvtp0 * Inv_W + model.pdvtp0 * Inv_LW;
        size_params.dvtp1 =
            model.dvtp1 + model.ldvtp1 * Inv_L + model.wdvtp1 * Inv_W + model.pdvtp1 * Inv_LW;
        size_params.dvtp2 = model.dvtp2 		/* v4.7  */
                                     + model.ldvtp2 * Inv_L
                                     + model.wdvtp2 * Inv_W
                                     + model.pdvtp2 * Inv_LW;
        size_params.dvtp3 = model.dvtp3 		/* v4.7  */
                                     + model.ldvtp3 * Inv_L
                                     + model.wdvtp3 * Inv_W
                                     + model.pdvtp3 * Inv_LW;
        size_params.dvtp4 = model.dvtp4 		/* v4.7  */
                                     + model.ldvtp4 * Inv_L
                                     + model.wdvtp4 * Inv_W
                                     + model.pdvtp4 * Inv_LW;
        size_params.dvtp5 = model.dvtp5 		/* v4.7  */
                                     + model.ldvtp5 * Inv_L
                                     + model.wdvtp5 * Inv_W
                                     + model.pdvtp5 * Inv_LW;
        size_params.dvt0 =
            model.dvt0 + model.ldvt0 * Inv_L + model.wdvt0 * Inv_W + model.pdvt0 * Inv_LW;
        size_params.dvt1 =
            model.dvt1 + model.ldvt1 * Inv_L + model.wdvt1 * Inv_W + model.pdvt1 * Inv_LW;
        size_params.dvt2 =
            model.dvt2 + model.ldvt2 * Inv_L + model.wdvt2 * Inv_W + model.pdvt2 * Inv_LW;
        size_params.dvt0w =
            model.dvt0w + model.ldvt0w * Inv_L + model.wdvt0w * Inv_W + model.pdvt0w * Inv_LW;
        size_params.dvt1w =
            model.dvt1w + model.ldvt1w * Inv_L + model.wdvt1w * Inv_W + model.pdvt1w * Inv_LW;
        size_params.dvt2w =
            model.dvt2w + model.ldvt2w * Inv_L + model.wdvt2w * Inv_W + model.pdvt2w * Inv_LW;
        size_params.drout =
            model.drout + model.ldrout * Inv_L + model.wdrout * Inv_W + model.pdrout * Inv_LW;
        size_params.dsub =
            model.dsub + model.ldsub * Inv_L + model.wdsub * Inv_W + model.pdsub * Inv_LW;
        size_params.vth0 =
            model.vth0 + model.lvth0 * Inv_L + model.wvth0 * Inv_W + model.pvth0 * Inv_LW;
        size_params.ua = model.ua + model.lua * Inv_L + model.wua * Inv_W + model.pua * Inv_LW;
        size_params.ua1 = model.ua1 + model.lua1 * Inv_L + model.wua1 * Inv_W + model.pua1 * Inv_LW;
        size_params.ub = model.ub + model.lub * Inv_L + model.wub * Inv_W + model.pub_ * Inv_LW;
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
        /*high k mobility*/
        size_params.ucs = model.ucs + model.lucs * Inv_L + model.wucs * Inv_W + model.pucs * Inv_LW;
        size_params.ucste =
            model.ucste + model.lucste * Inv_L + model.wucste * Inv_W + model.pucste * Inv_LW;

        size_params.voff =
            model.voff + model.lvoff * Inv_L + model.wvoff * Inv_W + model.pvoff * Inv_LW;
        size_params.tvoff =
            model.tvoff + model.ltvoff * Inv_L + model.wtvoff * Inv_W + model.ptvoff * Inv_LW;
        size_params.minv =
            model.minv + model.lminv * Inv_L + model.wminv * Inv_W + model.pminv * Inv_LW;
        size_params.minvcv =
            model.minvcv + model.lminvcv * Inv_L + model.wminvcv * Inv_W + model.pminvcv * Inv_LW;
        size_params.fprout =
            model.fprout + model.lfprout * Inv_L + model.wfprout * Inv_W + model.pfprout * Inv_LW;
        size_params.pdits =
            model.pdits + model.lpdits * Inv_L + model.wpdits * Inv_W + model.ppdits * Inv_LW;
        size_params.pditsd =
            model.pditsd + model.lpditsd * Inv_L + model.wpditsd * Inv_W + model.ppditsd * Inv_LW;
        size_params.delta =
            model.delta + model.ldelta * Inv_L + model.wdelta * Inv_W + model.pdelta * Inv_LW;
        size_params.rdsw =
            model.rdsw + model.lrdsw * Inv_L + model.wrdsw * Inv_W + model.prdsw * Inv_LW;
        size_params.rdw = model.rdw + model.lrdw * Inv_L + model.wrdw * Inv_W + model.prdw * Inv_LW;
        size_params.rsw = model.rsw + model.lrsw * Inv_L + model.wrsw * Inv_W + model.prsw * Inv_LW;
        size_params.prwg =
            model.prwg + model.lprwg * Inv_L + model.wprwg * Inv_W + model.pprwg * Inv_LW;
        size_params.prwb =
            model.prwb + model.lprwb * Inv_L + model.wprwb * Inv_W + model.pprwb * Inv_LW;
        size_params.prt = model.prt + model.lprt * Inv_L + model.wprt * Inv_W + model.pprt * Inv_LW;
        size_params.eta0 =
            model.eta0 + model.leta0 * Inv_L + model.weta0 * Inv_W + model.peta0 * Inv_LW;
        size_params.teta0 = model.teta0 		/* v4.7  */
				    + model.lteta0 * Inv_L
				    + model.wteta0 * Inv_W
				    + model.pteta0 * Inv_LW;
        size_params.tvoffcv = model.tvoffcv 	/* v4.8.0  */
				    + model.ltvoffcv * Inv_L
				    + model.wtvoffcv * Inv_W
				    + model.ptvoffcv * Inv_LW;
        size_params.etab =
            model.etab + model.letab * Inv_L + model.wetab * Inv_W + model.petab * Inv_LW;
        size_params.pclm =
            model.pclm + model.lpclm * Inv_L + model.wpclm * Inv_W + model.ppclm * Inv_LW;
        size_params.pdibl1 = model.pdiblc1
            + model.lpdiblc1 * Inv_L
            + model.wpdiblc1 * Inv_W
            + model.ppdiblc1 * Inv_LW;
        size_params.pdibl2 = model.pdiblc2
            + model.lpdiblc2 * Inv_L
            + model.wpdiblc2 * Inv_W
            + model.ppdiblc2 * Inv_LW;
        size_params.pdiblb = model.pdiblcb
            + model.lpdiblcb * Inv_L
            + model.wpdiblcb * Inv_W
            + model.ppdiblcb * Inv_LW;
        size_params.pscbe1 =
            model.pscbe1 + model.lpscbe1 * Inv_L + model.wpscbe1 * Inv_W + model.ppscbe1 * Inv_LW;
        size_params.pscbe2 =
            model.pscbe2 + model.lpscbe2 * Inv_L + model.wpscbe2 * Inv_W + model.ppscbe2 * Inv_LW;
        size_params.pvag =
            model.pvag + model.lpvag * Inv_L + model.wpvag * Inv_W + model.ppvag * Inv_LW;
        size_params.wr = model.wr + model.lwr * Inv_L + model.wwr * Inv_W + model.pwr * Inv_LW;
        size_params.dwg = model.dwg + model.ldwg * Inv_L + model.wdwg * Inv_W + model.pdwg * Inv_LW;
        size_params.dwb = model.dwb + model.ldwb * Inv_L + model.wdwb * Inv_W + model.pdwb * Inv_LW;
        size_params.b0 = model.b0 + model.lb0 * Inv_L + model.wb0 * Inv_W + model.pb0 * Inv_LW;
        size_params.b1 = model.b1 + model.lb1 * Inv_L + model.wb1 * Inv_W + model.pb1 * Inv_LW;
        size_params.alpha0 =
            model.alpha0 + model.lalpha0 * Inv_L + model.walpha0 * Inv_W + model.palpha0 * Inv_LW;
        size_params.alpha1 =
            model.alpha1 + model.lalpha1 * Inv_L + model.walpha1 * Inv_W + model.palpha1 * Inv_LW;
        size_params.beta0 =
            model.beta0 + model.lbeta0 * Inv_L + model.wbeta0 * Inv_W + model.pbeta0 * Inv_LW;
        size_params.agidl =
            model.agidl + model.lagidl * Inv_L + model.wagidl * Inv_W + model.pagidl * Inv_LW;
        size_params.bgidl =
            model.bgidl + model.lbgidl * Inv_L + model.wbgidl * Inv_W + model.pbgidl * Inv_LW;
        size_params.cgidl =
            model.cgidl + model.lcgidl * Inv_L + model.wcgidl * Inv_W + model.pcgidl * Inv_LW;
        size_params.egidl =
            model.egidl + model.legidl * Inv_L + model.wegidl * Inv_W + model.pegidl * Inv_LW;
        size_params.rgidl = model.rgidl		/* v4.7 New GIDL/GISL */
                                     + model.lrgidl * Inv_L
                                     + model.wrgidl * Inv_W
                                     + model.prgidl * Inv_LW;
        size_params.kgidl = model.kgidl		/* v4.7 New GIDL/GISL */
                                     + model.lkgidl * Inv_L
                                     + model.wkgidl * Inv_W
                                     + model.pkgidl * Inv_LW;
        size_params.fgidl = model.fgidl		/* v4.7 New GIDL/GISL */
                                     + model.lfgidl * Inv_L
                                     + model.wfgidl * Inv_W
                                     + model.pfgidl * Inv_LW;
        size_params.agisl =
            model.agisl + model.lagisl * Inv_L + model.wagisl * Inv_W + model.pagisl * Inv_LW;
        size_params.bgisl =
            model.bgisl + model.lbgisl * Inv_L + model.wbgisl * Inv_W + model.pbgisl * Inv_LW;
        size_params.cgisl =
            model.cgisl + model.lcgisl * Inv_L + model.wcgisl * Inv_W + model.pcgisl * Inv_LW;
        size_params.egisl =
            model.egisl + model.legisl * Inv_L + model.wegisl * Inv_W + model.pegisl * Inv_LW;
        size_params.rgisl = model.rgisl		/* v4.7 New GIDL/GISL */
                                     + model.lrgisl * Inv_L
                                     + model.wrgisl * Inv_W
                                     + model.prgisl * Inv_LW;
        size_params.kgisl = model.kgisl		/* v4.7 New GIDL/GISL */
                                     + model.lkgisl * Inv_L
                                     + model.wkgisl * Inv_W
                                     + model.pkgisl * Inv_LW;
        size_params.fgisl = model.fgisl		/* v4.7 New GIDL/GISL */
                                     + model.lfgisl * Inv_L
                                     + model.wfgisl * Inv_W
                                     + model.pfgisl * Inv_LW;
        size_params.aigc =
            model.aigc + model.laigc * Inv_L + model.waigc * Inv_W + model.paigc * Inv_LW;
        size_params.bigc =
            model.bigc + model.lbigc * Inv_L + model.wbigc * Inv_W + model.pbigc * Inv_LW;
        size_params.cigc =
            model.cigc + model.lcigc * Inv_L + model.wcigc * Inv_W + model.pcigc * Inv_LW;
        size_params.aigsd =
            model.aigsd + model.laigsd * Inv_L + model.waigsd * Inv_W + model.paigsd * Inv_LW;
        size_params.bigsd =
            model.bigsd + model.lbigsd * Inv_L + model.wbigsd * Inv_W + model.pbigsd * Inv_LW;
        size_params.cigsd =
            model.cigsd + model.lcigsd * Inv_L + model.wcigsd * Inv_W + model.pcigsd * Inv_LW;
        size_params.aigs =
            model.aigs + model.laigs * Inv_L + model.waigs * Inv_W + model.paigs * Inv_LW;
        size_params.bigs =
            model.bigs + model.lbigs * Inv_L + model.wbigs * Inv_W + model.pbigs * Inv_LW;
        size_params.cigs =
            model.cigs + model.lcigs * Inv_L + model.wcigs * Inv_W + model.pcigs * Inv_LW;
        size_params.aigd =
            model.aigd + model.laigd * Inv_L + model.waigd * Inv_W + model.paigd * Inv_LW;
        size_params.bigd =
            model.bigd + model.lbigd * Inv_L + model.wbigd * Inv_W + model.pbigd * Inv_LW;
        size_params.cigd =
            model.cigd + model.lcigd * Inv_L + model.wcigd * Inv_W + model.pcigd * Inv_LW;
        size_params.aigbacc = model.aigbacc
            + model.laigbacc * Inv_L
            + model.waigbacc * Inv_W
            + model.paigbacc * Inv_LW;
        size_params.bigbacc = model.bigbacc
            + model.lbigbacc * Inv_L
            + model.wbigbacc * Inv_W
            + model.pbigbacc * Inv_LW;
        size_params.cigbacc = model.cigbacc
            + model.lcigbacc * Inv_L
            + model.wcigbacc * Inv_W
            + model.pcigbacc * Inv_LW;
        size_params.aigbinv = model.aigbinv
            + model.laigbinv * Inv_L
            + model.waigbinv * Inv_W
            + model.paigbinv * Inv_LW;
        size_params.bigbinv = model.bigbinv
            + model.lbigbinv * Inv_L
            + model.wbigbinv * Inv_W
            + model.pbigbinv * Inv_LW;
        size_params.cigbinv = model.cigbinv
            + model.lcigbinv * Inv_L
            + model.wcigbinv * Inv_W
            + model.pcigbinv * Inv_LW;
        size_params.nigc =
            model.nigc + model.lnigc * Inv_L + model.wnigc * Inv_W + model.pnigc * Inv_LW;
        size_params.nigbacc = model.nigbacc
            + model.lnigbacc * Inv_L
            + model.wnigbacc * Inv_W
            + model.pnigbacc * Inv_LW;
        size_params.nigbinv = model.nigbinv
            + model.lnigbinv * Inv_L
            + model.wnigbinv * Inv_W
            + model.pnigbinv * Inv_LW;
        size_params.ntox =
            model.ntox + model.lntox * Inv_L + model.wntox * Inv_W + model.pntox * Inv_LW;
        size_params.eigbinv = model.eigbinv
            + model.leigbinv * Inv_L
            + model.weigbinv * Inv_W
            + model.peigbinv * Inv_LW;
        size_params.pigcd =
            model.pigcd + model.lpigcd * Inv_L + model.wpigcd * Inv_W + model.ppigcd * Inv_LW;
        size_params.poxedge = model.poxedge
            + model.lpoxedge * Inv_L
            + model.wpoxedge * Inv_W
            + model.ppoxedge * Inv_LW;
        size_params.xrcrg1 =
            model.xrcrg1 + model.lxrcrg1 * Inv_L + model.wxrcrg1 * Inv_W + model.pxrcrg1 * Inv_LW;
        size_params.xrcrg2 =
            model.xrcrg2 + model.lxrcrg2 * Inv_L + model.wxrcrg2 * Inv_W + model.pxrcrg2 * Inv_LW;
        size_params.lambda =
            model.lambda + model.llambda * Inv_L + model.wlambda * Inv_W + model.plambda * Inv_LW;
        size_params.vtl = model.vtl + model.lvtl * Inv_L + model.wvtl * Inv_W + model.pvtl * Inv_LW;
        size_params.xn = model.xn + model.lxn * Inv_L + model.wxn * Inv_W + model.pxn * Inv_LW;
        size_params.vfbsdoff = model.vfbsdoff
            + model.lvfbsdoff * Inv_L
            + model.wvfbsdoff * Inv_W
            + model.pvfbsdoff * Inv_LW;
        size_params.tvfbsdoff = model.tvfbsdoff
            + model.ltvfbsdoff * Inv_L
            + model.wtvfbsdoff * Inv_W
            + model.ptvfbsdoff * Inv_LW;

        size_params.cgsl =
            model.cgsl + model.lcgsl * Inv_L + model.wcgsl * Inv_W + model.pcgsl * Inv_LW;
        size_params.cgdl =
            model.cgdl + model.lcgdl * Inv_L + model.wcgdl * Inv_W + model.pcgdl * Inv_LW;
        size_params.ckappas = model.ckappas
            + model.lckappas * Inv_L
            + model.wckappas * Inv_W
            + model.pckappas * Inv_LW;
        size_params.ckappad = model.ckappad
            + model.lckappad * Inv_L
            + model.wckappad * Inv_W
            + model.pckappad * Inv_LW;
        size_params.cf = model.cf + model.lcf * Inv_L + model.wcf * Inv_W + model.pcf * Inv_LW;
        size_params.clc = model.clc + model.lclc * Inv_L + model.wclc * Inv_W + model.pclc * Inv_LW;
        size_params.cle = model.cle + model.lcle * Inv_L + model.wcle * Inv_W + model.pcle * Inv_LW;
        size_params.vfbcv =
            model.vfbcv + model.lvfbcv * Inv_L + model.wvfbcv * Inv_W + model.pvfbcv * Inv_LW;
        size_params.acde =
            model.acde + model.lacde * Inv_L + model.wacde * Inv_W + model.pacde * Inv_LW;
        size_params.moin =
            model.moin + model.lmoin * Inv_L + model.wmoin * Inv_W + model.pmoin * Inv_LW;
        size_params.noff =
            model.noff + model.lnoff * Inv_L + model.wnoff * Inv_W + model.pnoff * Inv_LW;
        size_params.voffcv =
            model.voffcv + model.lvoffcv * Inv_L + model.wvoffcv * Inv_W + model.pvoffcv * Inv_LW;
        size_params.kvth0we = model.kvth0we
            + model.lkvth0we * Inv_L
            + model.wkvth0we * Inv_W
            + model.pkvth0we * Inv_LW;
        size_params.k2we =
            model.k2we + model.lk2we * Inv_L + model.wk2we * Inv_W + model.pk2we * Inv_LW;
        size_params.ku0we =
            model.ku0we + model.lku0we * Inv_L + model.wku0we * Inv_W + model.pku0we * Inv_LW;

        size_params.abulkCVfactor =
            1.0 + pow((size_params.clc / size_params.leffCV), size_params.cle);

        T0 = (TRatio - 1.0);

        PowWeffWr = pow(size_params.weffCJ * 1.0e6, size_params.wr) * intp.nf;

        T1 = 0.0;
        T2 = 0.0;
        T3 = 0.0;
        T4 = 0.0;
        size_params.ucs = size_params.ucs * pow(TRatio, size_params.ucste);
        if (model.tempmod == 0) {
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
            if (model.tempmod == 3) {
                size_params.ua = size_params.ua * pow(TRatio, size_params.ua1);
                size_params.ub = size_params.ub * pow(TRatio, size_params.ub1);
                size_params.uc = size_params.uc * pow(TRatio, size_params.uc1);
                size_params.ud = size_params.ud * pow(TRatio, size_params.ud1);
            } else {
                /* tempMod = 1, 2 */
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

        if (T1 < 0.0) {
            T1 = 0.0;
            panic!("Warning: Rdw at current temperature is negative; set to 0.\n");
        }
        if (T2 < 0.0) {
            T2 = 0.0;
            panic!("Warning: Rdwmin at current temperature is negative; set to 0.\n");
        }
        size_params.rd0 = T1 / PowWeffWr;
        size_params.rdwmin = T2 / PowWeffWr;
        if (T3 < 0.0) {
            T3 = 0.0;
            panic!("Warning: Rsw at current temperature is negative; set to 0.\n");
        }
        if (T4 < 0.0) {
            T4 = 0.0;
            panic!("Warning: Rswmin at current temperature is negative; set to 0.\n");
        }
        size_params.rs0 = T3 / PowWeffWr;
        size_params.rswmin = T4 / PowWeffWr;

        if (size_params.u0 > 1.0) {
            size_params.u0 = size_params.u0 / 1.0e4;
        }

        /* mobility channel length dependence */
        T5 = 1.0 - size_params.up * exp(-size_params.leff / size_params.lp);
        size_params.u0temp = size_params.u0 * T5 * pow(TRatio, size_params.ute);
        if (size_params.eu < 0.0) {
            size_params.eu = 0.0;
            panic!("Warning: eu has been negative; reset to 0.0.\n");
        }
        if (size_params.ucs < 0.0) {
            size_params.ucs = 0.0;
            panic!("Warning: ucs has been negative; reset to 0.0.\n");
        }

        size_params.vfbsdoff = size_params.vfbsdoff * (1.0 + size_params.tvfbsdoff * delTemp);
        size_params.voff = size_params.voff * (1.0 + size_params.tvoff * delTemp);

        size_params.nfactor = size_params.nfactor + size_params.tnfactor * delTemp / Tnom; /* v4.7 temp dep of leakage currents */
        size_params.voffcv = size_params.voffcv * (1.0 + size_params.tvoffcv * delTemp); /*	 v4.7 temp dep of leakage currents */
        size_params.eta0 = size_params.eta0 + size_params.teta0 * delTemp / Tnom; /*	 v4.7 temp dep of leakage currents */

        /* Source End Velocity Limit  */
        if ((model.vtlGiven) && (model.vtl > 0.0)) {
            if (model.lc < 0.0) {
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

        if (!model.ndepGiven && model.gamma1Given) {
            T0 = size_params.gamma1 * model_derived.coxe;
            size_params.ndep = 3.01248e22 * T0 * T0;
        }

        size_params.phi = Vtm0 * log(size_params.ndep / ni) + size_params.phin + 0.4;

        size_params.sqrtPhi = sqrt(size_params.phi);
        size_params.phis3 = size_params.sqrtPhi * size_params.phi;

        size_params.Xdep0 =
            sqrt(2.0 * epssub / (Q * size_params.ndep * 1.0e6)) * size_params.sqrtPhi;
        size_params.sqrtXdep0 = sqrt(size_params.Xdep0);

        if (model.mtrlmod == 0) {
            size_params.litl = sqrt(3.0 * 3.9 / epsrox * size_params.xj * toxe);
        } else {
            size_params.litl = sqrt(model.epsrsub / epsrox * size_params.xj * toxe);
        }

        size_params.vbi = Vtm0 * log(size_params.nsd * size_params.ndep / (ni * ni));

        if (model.mtrlmod == 0) {
            if (size_params.ngate > 0.0) {
                size_params.vfbsd = Vtm0 * log(size_params.ngate / size_params.nsd);
            } else {
                size_params.vfbsd = 0.0;
            }
        } else {
            T0 = Vtm0 * log(size_params.nsd / ni);
            T1 = 0.5 * Eg0;
            if (T0 > T1) {
                T0 = T1;
            }
            T2 = model.easub + T1 - model.p() * T0;
            size_params.vfbsd = model.phig - T2;
        }

        size_params.cdep0 = sqrt(Q * epssub * size_params.ndep * 1.0e6 / 2.0 / size_params.phi);

        size_params.ToxRatio = exp(size_params.ntox * log(model.toxref / toxe)) / toxe / toxe;
        size_params.ToxRatioEdge =
            exp(size_params.ntox * log(model.toxref / (toxe * size_params.poxedge)))
                / toxe
                / toxe
                / size_params.poxedge
                / size_params.poxedge;
        size_params.Aechvb = if model.p() == 1.0 {
            // FIXME: MOS enum
            4.97232e-7
        } else {
            3.42537e-7
        };
        size_params.Bechvb = if model.p() == 1.0 {
            // FIXME: MOS enum
            7.45669e11
        } else {
            1.16645e12
        };
        size_params.AechvbEdgeS =
            size_params.Aechvb * size_params.weff * model.dlcig * size_params.ToxRatioEdge;
        size_params.AechvbEdgeD =
            size_params.Aechvb * size_params.weff * model.dlcigd * size_params.ToxRatioEdge;
        size_params.BechvbEdge = -size_params.Bechvb * toxe * size_params.poxedge;
        size_params.Aechvb *= size_params.weff * size_params.leff * size_params.ToxRatio;
        size_params.Bechvb *= -toxe;

        size_params.mstar = 0.5 + atan(size_params.minv) / PI;
        size_params.mstarcv = 0.5 + atan(size_params.minvcv) / PI;
        size_params.voffcbn = size_params.voff + model.voffl / size_params.leff;
        size_params.voffcbncv = size_params.voffcv + model.voffcvl / size_params.leff;

        size_params.ldeb = sqrt(epssub * Vtm0 / (Q * size_params.ndep * 1.0e6)) / 3.0;
        size_params.acde *= pow((size_params.ndep / 2.0e16), -0.25);

        if (model.k1Given || model.k2Given) {
            if (!model.k1Given) {
                size_params.k1 = 0.53;
                panic!("Warning: k1 should be specified with k2.\n");
            }
            if (!model.k2Given) {
                size_params.k2 = -0.0186;
                panic!("Warning: k2 should be specified with k1.\n");
            }
            if (model.nsubGiven) {
                panic!("Warning: nsub is ignored because k1 or k2 is given.\n",);
            }
            if (model.xtGiven) {
                panic!("Warning: xt is ignored because k1 or k2 is given.\n",);
            }
            if (model.vbxGiven) {
                panic!("Warning: vbx is ignored because k1 or k2 is given.\n",);
            }
            if (model.gamma1Given) {
                panic!("Warning: gamma1 is ignored because k1 or k2 is given.\n",);
            }
            if (model.gamma2Given) {
                panic!("Warning: gamma2 is ignored because k1 or k2 is given.\n",);
            }
        } else {
            if (!model.vbxGiven) {
                size_params.vbx = size_params.phi
                    - 7.7348e-4 * size_params.ndep * size_params.xt * size_params.xt;
            }
            if (size_params.vbx > 0.0) {
                size_params.vbx = -size_params.vbx;
            }
            if (size_params.vbm > 0.0) {
                size_params.vbm = -size_params.vbm;
            }

            if (!model.gamma1Given) {
                size_params.gamma1 = 5.753e-12 * sqrt(size_params.ndep) / model_derived.coxe;
            }
            if (!model.gamma2Given) {
                size_params.gamma2 = 5.753e-12 * sqrt(size_params.nsub) / model_derived.coxe;
            }

            T0 = size_params.gamma1 - size_params.gamma2;
            T1 = sqrt(size_params.phi - size_params.vbx) - size_params.sqrtPhi;
            T2 = sqrt(size_params.phi * (size_params.phi - size_params.vbm)) - size_params.phi;
            size_params.k2 = T0 * T1 / (2.0 * T2 + size_params.vbm);
            size_params.k1 =
                size_params.gamma2 - 2.0 * size_params.k2 * sqrt(size_params.phi - size_params.vbm);
        }

        if (!model.vfbGiven) {
            if (model.vth0Given) {
                size_params.vfb = model.p() * size_params.vth0
                    - size_params.phi
                    - size_params.k1 * size_params.sqrtPhi;
            } else {
                if ((model.mtrlmod) && (model.phigGiven) && (model.nsubGiven)) {
                    T0 = Vtm0 * log(size_params.nsub / ni);
                    T1 = 0.5 * Eg0;
                    if (T0 > T1) {
                        T0 = T1;
                    }
                    T2 = model.easub + T1 + model.p() * T0;
                    size_params.vfb = model.phig - T2;
                } else {
                    size_params.vfb = -1.0;
                }
            }
        }
        if (!model.vth0Given) {
            size_params.vth0 = model.p()
                * (size_params.vfb + size_params.phi + size_params.k1 * size_params.sqrtPhi);
        }

        size_params.k1ox = size_params.k1 * toxe / model.toxm;

        tmp = sqrt(epssub / (epsrox * EPS0) * toxe * size_params.Xdep0);
        T0 = size_params.dsub * size_params.leff / tmp;
        if (T0 < EXP_THRESHOLD) {
            T1 = exp(T0);
            T2 = T1 - 1.0;
            T3 = T2 * T2;
            T4 = T3 + 2.0 * T1 * MIN_EXP;
            size_params.theta0vb0 = T1 / T4;
        } else {
            size_params.theta0vb0 = 1.0 / (MAX_EXP - 2.0);
        }

        T0 = size_params.drout * size_params.leff / tmp;
        if (T0 < EXP_THRESHOLD) {
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
        if (T0 < EXP_THRESHOLD) {
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
        if (T0 < EXP_THRESHOLD) {
            T1 = exp(T0);
            T2 = T1 - 1.0;
            T3 = T2 * T2;
            T4 = T3 + 2.0 * T1 * MIN_EXP;
            T9 = T1 / T4;
        } else {
            T9 = 1.0 / (MAX_EXP - 2.0);
        }
        T9 = size_params.dvt0 * T9 * tmp1;

        T4 = toxe * size_params.phi / (size_params.weff + size_params.w0);

        T0 = sqrt(1.0 + size_params.lpe0 / size_params.leff);
        if ((model.tempMod == 1) || (model.tempmod == 0)) {
            T3 = (size_params.kt1 + size_params.kt1l / size_params.leff) * (TRatio - 1.0);
        }
        if ((model.tempMod == 2) || (model.tempmod == 3)) {
            T3 = -size_params.kt1 * (TRatio - 1.0);
        }

        T5 = size_params.k1ox * (T0 - 1.0) * size_params.sqrtPhi + T3;
        size_params.vfbzbfactor = -T8 - T9 + size_params.k3 * T4 + T5
            - size_params.phi
            - size_params.k1 * size_params.sqrtPhi;

        /* stress effect */

        wlod = model.wlod;
        if (model.wlod < 0.0) {
            panic!(
                "Warning: WLOD = %g is less than 0. 0.0 is used\n",
                model.wlod,
            );
            wlod = 0.0;
        }
        T0 = pow(Lnew, model.llodku0);
        W_tmp = Wnew + wlod;
        T1 = pow(W_tmp, model.wlodku0);
        tmp1 = model.lku0 / T0 + model.wku0 / T1 + model.pku0 / (T0 * T1);
        size_params.ku0 = 1.0 + tmp1;

        T0 = pow(Lnew, model.llodvth);
        T1 = pow(W_tmp, model.wlodvth);
        tmp1 = model.lkvth0 / T0 + model.wkvth0 / T1 + model.pkvth0 / (T0 * T1);
        size_params.kvth0 = 1.0 + tmp1;
        size_params.kvth0 = sqrt(size_params.kvth0 * size_params.kvth0 + DELTA);

        T0 = (TRatio - 1.0);
        size_params.ku0temp = size_params.ku0 * (1.0 + model.tku0 * T0) + DELTA;

        Inv_saref = 1.0 / (model.saref + 0.5 * Ldrn);
        Inv_sbref = 1.0 / (model.sbref + 0.5 * Ldrn);
        size_params.inv_od_ref = Inv_saref + Inv_sbref;
        size_params.rho_ref = model.ku0 / size_params.ku0temp * size_params.inv_od_ref;

        /*high k*/
        /*Calculate VgsteffVth for mobMod=3*/
        if (model.mobmod == 3) {
            /*Calculate n @ Vbs=Vds=0*/
            lt1 = model_derived.factor1 * size_params.sqrtXdep0;
            T0 = size_params.dvt1 * size_params.leff / lt1;
            if (T0 < EXP_THRESHOLD) {
                T1 = exp(T0);
                T2 = T1 - 1.0;
                T3 = T2 * T2;
                T4 = T3 + 2.0 * T1 * MIN_EXP;
                Theta0 = T1 / T4;
            } else {
                Theta0 = 1.0 / (MAX_EXP - 2.0);
            }

            tmp1 = epssub / size_params.Xdep0;
            tmp2 = size_params.nfactor * tmp1;
            tmp3 = (tmp2 + size_params.cdsc * Theta0 + size_params.cit) / model_derived.coxe;
            if (tmp3 >= -0.5) {
                n0 = 1.0 + tmp3;
            } else {
                T0 = 1.0 / (3.0 + 8.0 * tmp3);
                n0 = (1.0 + 3.0 * tmp3) * T0;
            }

            T0 = n0 * model_derived.vtm;
            T1 = size_params.voffcbn;
            T2 = T1 / T0;
            if (T2 < -EXP_THRESHOLD) {
                T3 = model_derived.coxe * MIN_EXP / size_params.cdep0;
                T4 = size_params.mstar + T3 * n0;
            } else if (T2 > EXP_THRESHOLD) {
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
    if ((intp.sa > 0.0)
        && (intp.sb > 0.0)
        && ((intp.nf == 1.0) || ((intp.nf > 1.0) && (intp.sd > 0.0))))
    {
        Inv_sa = 0.0;
        Inv_sb = 0.0;

        kvsat = model.kvsat;
        if (model.kvsat < -1.0) {
            panic!(
                "Warning: KVSAT = %g is too small; -1.0 is used.\n",
                model.kvsat,
            );
            kvsat = -1.0;
        }
        if (model.kvsat > 1.0) {
            panic!(
                "Warning: KVSAT = %g is too big; 1.0 is used.\n",
                model.kvsat,
            );
            kvsat = 1.0;
        }
        let nfi = intp.nf as usize;
        for i in 0..nfi {
            T0 = 1.0 / intp.nf / (intp.sa + 0.5 * Ldrn + i as f64 * (intp.sd + Ldrn));
            T1 = 1.0 / intp.nf / (intp.sb + 0.5 * Ldrn + i as f64 * (intp.sd + Ldrn));
            Inv_sa += T0;
            Inv_sb += T1;
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
        if ((!intp.scaGiven) && (!intp.scbGiven) && (!intp.sccGiven)) {
            if ((intp.scGiven) && (intp.sc > 0.0)) {
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
                panic!("Warning: No WPE as none of SCA, SCB, SCC, SC is given and/or SC not positive.\n");
            }
        }

        if (intp.sca < 0.0) {
            intp.sca = 0.0;
            panic!("Warning: SCA = %g is negative. Set to 0.0.\n", intp.sca);
        }
        if (intp.scb < 0.0) {
            intp.scb = 0.0;
            panic!("Warning: SCB = %g is negative. Set to 0.0.\n", intp.scb);
        }
        if (intp.scc < 0.0) {
            intp.scc = 0.0;
            panic!("Warning: SCC = %g is negative. Set to 0.0.\n", intp.scc);
        }
        if (intp.sc < 0.0) {
            intp.sc = 0.0;
            panic!("Warning: SC = %g is negative. Set to 0.0.\n", intp.sc);
        }
        /*4.6.2*/
        sceff = intp.sca + model.web * intp.scb + model.wec * intp.scc;
        intp.vth0 += size_params.kvth0we * sceff;
        intp.k2 += size_params.k2we * sceff;
        T3 = 1.0 + size_params.ku0we * sceff;
        if (T3 <= 0.0) {
            T3 = 0.0;
            panic!(
                "Warning: ku0we = %g is negatively too high. Negative mobility! \n",
                size_params.ku0we,
            );
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
    if (intp.vtfbphi1 < 0.0) {
        intp.vtfbphi1 = 0.0;
    }

    intp.vtfbphi2 = 4.0 * T3;
    if (intp.vtfbphi2 < 0.0) {
        intp.vtfbphi2 = 0.0;
    }

    if (intp.k2 < 0.0) {
        T0 = 0.5 * size_params.k1 / intp.k2;
        intp.vbsc = 0.9 * (size_params.phi - T0 * T0);
        if (intp.vbsc > -3.0) {
            intp.vbsc = -3.0;
        } else if (intp.vbsc < -30.0) {
            intp.vbsc = -30.0;
        }
    } else {
        intp.vbsc = -30.0;
    }
    if (intp.vbsc > size_params.vbm) {
        intp.vbsc = size_params.vbm;
    }
    intp.k2ox = intp.k2 * toxe / model.toxm;

    intp.vfbzb = size_params.vfbzbfactor + model.p() * intp.vth0;

    // intp.cgso = size_params.cgso;
    // intp.cgdo = size_params.cgdo;

    lnl = log(size_params.leff * 1.0e6);
    lnw = log(size_params.weff * 1.0e6);
    lnnf = log(intp.nf);

    bodymode = 5;
    if ((!model.rbps0Given) || (!model.rbpd0Given)) {
        bodymode = 1;
    } else if ((!model.rbsbx0Given && !model.rbsby0Given)
        || (!model.rbdbx0Given && !model.rbdby0Given))
    {
        bodymode = 3;
    }

    if (model.rbodymod == 2) {
        if (bodymode == 5) {
            rbsbx = model.rbsbx0
                * exp(model.rbsdbxl * lnl + model.rbsdbxw * lnw + model.rbsdbxnf * lnnf);
            rbsby = model.rbsby0
                * exp(model.rbsdbyl * lnl + model.rbsdbyw * lnw + model.rbsdbynf * lnnf);
            intp.rbsb = rbsbx * rbsby / (rbsbx + rbsby);

            rbdbx = model.rbdbx0
                * exp(model.rbsdbxl * lnl + model.rbsdbxw * lnw + model.rbsdbxnf * lnnf);
            rbdby = model.rbdby0
                * exp(model.rbsdbyl * lnl + model.rbsdbyw * lnw + model.rbsdbynf * lnnf);

            intp.rbdb = rbdbx * rbdby / (rbdbx + rbdby);
        }

        if ((bodymode == 3) || (bodymode == 5)) {
            intp.rbps =
                model.rbps0 * exp(model.rbpsl * lnl + model.rbpsw * lnw + model.rbpsnf * lnnf);
            intp.rbpd =
                model.rbpd0 * exp(model.rbpdl * lnl + model.rbpdw * lnw + model.rbpdnf * lnnf);
        }
        rbpbx = model.rbpbx0 * exp(model.rbpbxl * lnl + model.rbpbxw * lnw + model.rbpbxnf * lnnf);
        rbpby = model.rbpby0 * exp(model.rbpbyl * lnl + model.rbpbyw * lnw + model.rbpbynf * lnnf);

        intp.rbpb = rbpbx * rbpby / (rbpbx + rbpby);
    }

    if ((model.rbodymod == 1) || ((model.rbodymod == 2) && (bodymode == 5))) {
        if (intp.rbdb < 1.0e-3) {
            intp.grbdb = 1.0e3; /* in mho */
        } else {
            intp.grbdb = model.gbmin + 1.0 / intp.rbdb;
        }
        if (intp.rbpb < 1.0e-3) {
            intp.grbpb = 1.0e3;
        } else {
            intp.grbpb = model.gbmin + 1.0 / intp.rbpb;
        }
        if (intp.rbps < 1.0e-3) {
            intp.grbps = 1.0e3;
        } else {
            intp.grbps = model.gbmin + 1.0 / intp.rbps;
        }
        if (intp.rbsb < 1.0e-3) {
            intp.grbsb = 1.0e3;
        } else {
            intp.grbsb = model.gbmin + 1.0 / intp.rbsb;
        }
        if (intp.rbpd < 1.0e-3) {
            intp.grbpd = 1.0e3;
        } else {
            intp.grbpd = model.gbmin + 1.0 / intp.rbpd;
        }
    }

    if ((model.rbodymod == 2) && (bodymode == 3)) {
        intp.grbdb = model.gbmin;
        intp.grbsb = model.gbmin;
        if (intp.rbpb < 1.0e-3) {
            intp.grbpb = 1.0e3;
        } else {
            intp.grbpb = model.gbmin + 1.0 / intp.rbpb;
        }
        if (intp.rbps < 1.0e-3) {
            intp.grbps = 1.0e3;
        } else {
            intp.grbps = model.gbmin + 1.0 / intp.rbps;
        }
        if (intp.rbpd < 1.0e-3) {
            intp.grbpd = 1.0e3;
        } else {
            intp.grbpd = model.gbmin + 1.0 / intp.rbpd;
        }
    }

    if ((model.rbodymod == 2) && (bodymode == 1)) {
        intp.grbdb = model.gbmin;
        intp.grbsb = model.gbmin;
        intp.grbps = 1.0e3;
        intp.grbpd = 1.0e3;
        if (intp.rbpb < 1.0e-3) {
            intp.grbpb = 1.0e3;
        } else {
            intp.grbpb = model.gbmin + 1.0 / intp.rbpb;
        }
    }

    /*
     * Process geomertry dependent parasitics
     */

    intp.grgeltd = model.rshg * (intp.xgw + size_params.weffCJ / 3.0 / intp.ngcon)
        / (intp.ngcon * intp.nf * (Lnew - model.xgl));
    if (intp.grgeltd > 0.0) {
        intp.grgeltd = 1.0 / intp.grgeltd;
    } else {
        intp.grgeltd = 1.0e3; /* mho */
        if (model.rgatemod != 0) {
            panic!("Warning: The gate conductance reset to 1.0e3 mho.\n");
        }
    }

    DMCGeff = model.dmcg - model.dmcgt;
    DMCIeff = model.dmci;
    DMDGeff = model.dmdg - model.dmcgt;

    /* New Diode Model v4.7*/
    if (intp.sourcePerimeterGiven) {
        /* given */
        if (intp.sourcePerimeter == 0.0) {
            intp.Pseff = 0.0;
        } else if (intp.sourcePerimeter < 0.0) {
            panic!("Warning: Source Perimeter is specified as negative, it is set to zero.\n");
            intp.Pseff = 0.0;
        } else {
            if (model.permod == 0) {
                intp.Pseff = intp.sourcePerimeter;
            } else {
                intp.Pseff = intp.sourcePerimeter - size_params.weffCJ * intp.nf;
            }
        }
    } else
    /* not given */
    {
        BSIM4PAeffGeo(
            intp.nf,
            model.geomod,
            intp.min,
            size_params.weffCJ,
            DMCGeff,
            DMCIeff,
            DMDGeff,
            &(intp.Pseff),
            &dumPd,
            &dumAs,
            &dumAd,
        );
    }

    if (intp.Pseff < 0.0) {
        intp.Pseff = 0.0;
        panic!("Warning: Pseff is negative, it is set to zero.\n");
    }

    if (intp.drainPerimeterGiven) {
        /* given */
        if (intp.drainPerimeter == 0.0) {
            intp.Pdeff = 0.0;
        } else if (intp.drainPerimeter < 0.0) {
            panic!("Warning: Drain Perimeter is specified as negative, it is set to zero.\n");
            intp.Pdeff = 0.0;
        } else {
            if (model.permod == 0) {
                intp.Pdeff = intp.drainPerimeter;
            } else {
                intp.Pdeff = intp.drainPerimeter - size_params.weffCJ * intp.nf;
            }
        }
    } else
    /* not given */
    {
        BSIM4PAeffGeo(
            intp.nf,
            model.geomod,
            intp.min,
            size_params.weffCJ,
            DMCGeff,
            DMCIeff,
            DMDGeff,
            &dumPs,
            &(intp.Pdeff),
            &dumAs,
            &dumAd,
        );
    }

    if (intp.Pdeff < 0.0) {
        intp.Pdeff = 0.0; /*New Diode v4.7*/
        panic!("Warning: Pdeff is negative, it is set to zero.\n");
    }
    if (intp.sourceAreaGiven) {
        intp.Aseff = intp.sourceArea;
    } else {
        BSIM4PAeffGeo(
            intp.nf,
            model.geomod,
            intp.min,
            size_params.weffCJ,
            DMCGeff,
            DMCIeff,
            DMDGeff,
            &dumPs,
            &dumPd,
            &(intp.Aseff),
            &dumAd,
        );
    }
    if (intp.Aseff < 0.0) {
        intp.Aseff = 0.0; /* v4.7 */
        panic!("Warning: Aseff is negative, it is set to zero.\n");
    }
    if (intp.drainAreaGiven) {
        intp.Adeff = intp.drainArea;
    } else {
        BSIM4PAeffGeo(
            intp.nf,
            model.geomod,
            intp.min,
            size_params.weffCJ,
            DMCGeff,
            DMCIeff,
            DMDGeff,
            &dumPs,
            &dumPd,
            &dumAs,
            &(intp.Adeff),
        );
    }
    if (intp.Adeff < 0.0) {
        intp.Adeff = 0.0; /* v4.7 */
        panic!("Warning: Adeff is negative, it is set to zero.\n");
    }

    // // FIXME: these nodes won't be available here, figure out when & where to do this
    // // probably just do it here per geometry, and each instance will decide whether to use it.
    // /* Processing S/D resistance and conductance below */
    // if (intp.sNodePrime != intp.sNode) {
    //     intp.sourceConductance = 0.0;
    //     if (intp.sourceSquaresGiven) {
    //         intp.sourceConductance = model.sheetResistance * intp.sourceSquares;
    //     } else if (model.rgeomod > 0) {
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
    //     if (intp.sourceConductance > 0.0) {
    //         intp.sourceConductance = 1.0 / intp.sourceConductance;
    //     } else {
    //         intp.sourceConductance = 1.0e3; /* mho */
    //         panic!("Warning: Source conductance reset to 1.0e3 mho.\n");
    //     }
    // } else {
    //     intp.sourceConductance = 0.0;
    // }
    // if (intp.dNodePrime != intp.dNode) {
    //     intp.drainConductance = 0.0;
    //     if (intp.drainSquaresGiven) {
    //         intp.drainConductance = model.sheetResistance * intp.drainSquares;
    //     } else if (model.rgeomod > 0) {
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
    //     if (intp.drainConductance > 0.0) {
    //         intp.drainConductance = 1.0 / intp.drainConductance;
    //     } else {
    //         intp.drainConductance = 1.0e3; /* mho */
    //         panic!("Warning: Drain conductance reset to 1.0e3 mho.\n");
    //     }
    // } else {
    //     intp.drainConductance = 0.0;
    // }
    // /* End of Rsd processing */
    Nvtms = model_derived.vtm * model_derived.SjctEmissionCoeff;
    if ((intp.Aseff <= 0.0) && (intp.Pseff <= 0.0)) {
        SourceSatCurrent = 0.0;
    } else {
        SourceSatCurrent = intp.Aseff * model_derived.SjctTempSatCurDensity
            + intp.Pseff * model_derived.SjctSidewallTempSatCurDensity
            + size_params.weffCJ * intp.nf * model_derived.SjctGateSidewallTempSatCurDensity;
    }
    if (SourceSatCurrent > 0.0) {
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

                intp.vjsmFwd =
                    BSIM4DioIjthVjmEval(Nvtms, model.ijthsfwd, SourceSatCurrent, intp.XExpBVS);
                T0 = exp(intp.vjsmFwd / Nvtms);
                intp.IVjsmFwd = SourceSatCurrent * (T0 - intp.XExpBVS / T0 + intp.XExpBVS - 1.0);
                intp.SslpFwd = SourceSatCurrent * (T0 + intp.XExpBVS / T0) / Nvtms;

                T2 = model.ijthsrev / SourceSatCurrent;
                if (T2 < 1.0) {
                    T2 = 10.0;
                    panic!("Warning: ijthsrev too small and set to 10 times IsbSat.\n",);
                }
                intp.vjsmRev = -model.bvs - Nvtms * log((T2 - 1.0) / model.xjbvs);
                T1 = model.xjbvs * exp(-(model.bvs + intp.vjsmRev) / Nvtms);
                intp.IVjsmRev = SourceSatCurrent * (1.0 + T1);
                intp.SslpRev = -SourceSatCurrent * T1 / Nvtms;
            }
            _ => panic!("Specified dioMod %d not matched\n", model.diomod),
        }
    }

    Nvtmd = model_derived.vtm * model_derived.DjctEmissionCoeff;
    if ((intp.Adeff <= 0.0) && (intp.Pdeff <= 0.0)) {
        /* DrainSatCurrent = 1.0e-14; 	v4.7 */
        DrainSatCurrent = 0.0;
    } else {
        DrainSatCurrent = intp.Adeff * model_derived.DjctTempSatCurDensity
            + intp.Pdeff * model_derived.DjctSidewallTempSatCurDensity
            + size_params.weffCJ * intp.nf * model_derived.DjctGateSidewallTempSatCurDensity;
    }
    if (DrainSatCurrent > 0.0) {
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

                intp.vjdmFwd =
                    BSIM4DioIjthVjmEval(Nvtmd, model.ijthdfwd, DrainSatCurrent, intp.XExpBVD);
                T0 = exp(intp.vjdmFwd / Nvtmd);
                intp.IVjdmFwd = DrainSatCurrent * (T0 - intp.XExpBVD / T0 + intp.XExpBVD - 1.0);
                intp.DslpFwd = DrainSatCurrent * (T0 + intp.XExpBVD / T0) / Nvtmd;

                T2 = model.ijthdrev / DrainSatCurrent;
                if (T2 < 1.0) {
                    T2 = 10.0;
                    panic!("Warning: ijthdrev too small and set to 10 times IdbSat.\n",);
                }
                intp.vjdmRev = -model.bvd - Nvtmd * log((T2 - 1.0) / model.xjbvd); /* bugfix */
                T1 = model.xjbvd * exp(-(model.bvd + intp.vjdmRev) / Nvtmd);
                intp.IVjdmRev = DrainSatCurrent * (1.0 + T1);
                intp.DslpRev = -DrainSatCurrent * T1 / Nvtmd;
            }
            _ => panic!("Specified dioMod %d not matched\n", model.diomod),
        }
    }

    T7 = Eg0 / model_derived.vtm * T0;
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

    if (model.mtrlmod != 0 && model.mtrlcompatmod == 0) {
        /* Calculate TOXP from EOT */
        /* Calculate Vgs_eff @ Vgs = VDD with Poly Depletion Effect */
        Vtm0eot = KB_OVER_Q * model.tempeot;
        Vtmeot = Vtm0eot;
        vbieot = Vtm0eot * log(size_params.nsd * size_params.ndep / (ni * ni));
        phieot = Vtm0eot * log(size_params.ndep / ni) + size_params.phin + 0.4;
        tmp2 = intp.vfb + phieot;
        vddeot = model.p() * model.vddeot;
        T0 = model.epsrgate * EPS0;
        if ((size_params.ngate > 1.0e18)
            && (size_params.ngate < 1.0e25)
            && (vddeot > tmp2)
            && (T0 != 0.0))
        {
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
        if (T0 < EXP_THRESHOLD) {
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
        if (T0 < EXP_THRESHOLD) {
            T1 = exp(T0);
            T2 = T1 - 1.0;
            T3 = T2 * T2;
            T4 = T3 + 2.0 * T1 * MIN_EXP;
            T5 = T1 / T4;
        } else {
            T5 = 1.0 / (MAX_EXP - 2.0); /* 3.0 * MIN_EXP omitted */
        }
        T2 = size_params.dvt0w * T5 * V0;
        TempRatioeot = model.tempeot / model.tnom - 1.0;
        T0 = sqrt(1.0 + size_params.lpe0 / model.leffeot);
        T1 = size_params.k1ox * (T0 - 1.0) * sqrt(phieot)
            + (size_params.kt1 + size_params.kt1l / model.leffeot) * TempRatioeot;
        Vth_NarrowW = toxe * phieot / (model.weffeot + size_params.w0);
        Lpe_Vb = sqrt(1.0 + size_params.lpeb / model.leffeot);
        Vth = model.p() * intp.vth0 + (size_params.k1ox - size_params.k1) * sqrt(phieot) * Lpe_Vb
            - Delt_vth
            - T2
            + size_params.k3 * Vth_NarrowW
            + T1;

        /* Calculate n */
        tmp1 = epssub / size_params.Xdep0;
        tmp2 = size_params.nfactor * tmp1;
        tmp3 = (tmp2 + size_params.cdsc * Theta0 + size_params.cit) / model_derived.coxe;
        if (tmp3 >= -0.5) {
            n = 1.0 + tmp3;
        } else {
            T0 = 1.0 / (3.0 + 8.0 * tmp3);
            n = (1.0 + 3.0 * tmp3) * T0;
        }

        /* Vth correction for Pocket implant */
        if (size_params.dvtp0 > 0.0) {
            T3 = model.leffeot + size_params.dvtp0 * 2.0;
            if (model.tempmod < 2) {
                T4 = Vtmeot * log(model.leffeot / T3);
            } else {
                T4 = Vtm0eot * log(model.leffeot / T3);
            }
            Vth -= n * T4;
        }
        Vgsteff = Vgs_eff - Vth;
        /* calculating Toxp */
        T3 = model.p() * intp.vth0 - intp.vfb - phieot;
        T4 = T3 + T3;
        T5 = 2.5 * T3;

        vtfbphi2eot = 4.0 * T3;
        if (vtfbphi2eot < 0.0) {
            vtfbphi2eot = 0.0;
        }

        let niter = 0;
        toxpf = toxe;
        // loop  {
        // toxpi = toxpf;
        // tmp2 = 2.0e8 * toxpf;
        // T0 = (Vgsteff + vtfbphi2eot) / tmp2;
        // T1 = 1.0 + exp(model.bdos * 0.7 * log(T0));
        // Tcen = model.ados * 1.9e-9 / T1;
        // toxpf = toxe - epsrox/model.epsrsub * Tcen;
        // niter++;
        //   } while ((niter<=4)&&(abs(toxpf-toxpi)>1e-12));
        intp.toxp = toxpf;
        intp.coxp = epsrox * EPS0 / model.toxp;
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
    intp.size_params = size_params;
    return intp;
}

fn BSIM4DioIjthVjmEval(Nvtm: f64, Ijth: f64, Isb: f64, XExpBV: f64) -> f64 {
    let Tc = XExpBV;
    let Tb = 1.0 + Ijth / Isb - Tc;
    let EVjmovNv = 0.5 * (Tb + sqrt(Tb * Tb + 4.0 * Tc));
    return Nvtm * log(EVjmovNv);
}
