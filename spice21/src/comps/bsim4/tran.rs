use super::{Bsim4, Bsim4OpPoint};
use crate::analysis::TranState;

impl Bsim4 {
    ///
    /// Given the DC portion of an operating point,
    /// add transient impedances and currents for all time-varying components
    /// Argument `newop` is the partially-completed operating point,
    /// with its (admittedly ill-defined) DC and charge parameters set.
    ///
    pub(crate) fn tran_op(&self, newop: &mut Bsim4OpPoint, tran_state: &TranState) {
        let nqs_scaling_factor = 1.0e-9; // FIXME: lump this in an NQS thing

        // Initially zero all capacitances and their impedances
        // Many complicated paths through the code below do not assure they are otherwise initialized.
        let mut qgmid = 0.0;

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
        let mut dxpart = if newop.mode > 0 { 0.4 } else { 0.6 };
        let mut sxpart = 1.0 - dxpart;
        let mut ddxpart_dVd = 0.0;
        let mut ddxpart_dVg = 0.0;
        let mut ddxpart_dVb = 0.0;
        let mut ddxpart_dVs = 0.0;
        let mut dsxpart_dVd = 0.0;
        let mut dsxpart_dVg = 0.0;
        let mut dsxpart_dVb = 0.0;
        let mut dsxpart_dVs = 0.0;

        let Bsim4OpPoint {
            mut qbulk,
            mut qgate,
            mut qsrc,
            mut qdrn,
            ..
        } = *newop; // FIXME: the story with mutating these
        let Bsim4OpPoint { qgdo, qgso, qdef, .. } = *newop;
        let Bsim4OpPoint { cgdo, cgso, .. } = *newop;
        let Bsim4OpPoint {
            vbs,
            vbd,
            vgb,
            vgmb,
            vbs_jct,
            vbd_jct,
            ..
        } = *newop;

        // TODO: the BSIM4 reference implementation essentially bakes numerical integration in here,
        // ignoring the circuit/ analysis integration method.
        // All of these impedances are calculated as g = C/dt, e.g. using Backward Euler.
        // Figure out whether this is the implementation intent, or just for reference.

        let ag0 = 1.0 / tran_state.dt;
        if newop.mode > 0 {
            if self.model.trnqsmod == 0 {
                qdrn -= qgdo;
                if self.model.rgatemod == 3 {
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

                    let qgmb = self.size_params.cgbo * vgmb;
                    qgmid = qgdo + qgso + qgmb;
                    qbulk -= qgmb;
                    qsrc = -(qgate + qgmid + qbulk + qdrn);
                } else { // Default rgatemod==0
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

                    let qgb = self.size_params.cgbo * vgb;
                    qgate += qgdo + qgso + qgb;
                    qbulk -= qgb;
                    qsrc = -(qgate + qbulk + qdrn);
                }
                gcddb = (newop.cddb + newop.capbd + cgdo) * ag0;
                gcdsb = newop.cdsb * ag0;

                gcsdb = -(newop.cgdb + newop.cbdb + newop.cddb) * ag0;
                gcssb = (newop.capbs + cgso - (newop.cgsb + newop.cbsb + newop.cdsb)) * ag0;

                if self.model.rbodymod == 0 { // Default model 
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
                let qcheq = newop.qchqs;
                let CoxWL = self.model_derived.coxe * self.size_params.weffCV * self.intp.nf * self.size_params.leffCV;
                let T0 = qdef * nqs_scaling_factor / CoxWL;

                ggtg = T0 * newop.gcrgg;
                newop.gtg = ggtg;
                ggtd = T0 * newop.gcrgd;
                newop.gtd = ggtd;
                ggts = T0 * newop.gcrgs;
                newop.gts = ggts;
                ggtb = T0 * newop.gcrgb;
                newop.gtb = ggtb;
                gqdef = nqs_scaling_factor * ag0;

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
                    let Cdd = newop.cddb;
                    let Csd = -(newop.cgdb + newop.cddb + newop.cbdb);
                    ddxpart_dVd = (Cdd - dxpart * (Cdd + Csd)) / qcheq;
                    let Cdg = newop.cdgb;
                    let Csg = -(newop.cggb + newop.cdgb + newop.cbgb);
                    ddxpart_dVg = (Cdg - dxpart * (Cdg + Csg)) / qcheq;

                    let Cds = newop.cdsb;
                    let Css = -(newop.cgsb + newop.cdsb + newop.cbsb);
                    ddxpart_dVs = (Cds - dxpart * (Cds + Css)) / qcheq;

                    ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);
                }
                sxpart = 1.0 - dxpart;
                dsxpart_dVd = -ddxpart_dVd;
                dsxpart_dVg = -ddxpart_dVg;
                dsxpart_dVs = -ddxpart_dVs;
                dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);

                if self.model.rgatemod == 3 {
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

                    let qgmb = self.size_params.cgbo * vgmb;
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

                    let qgb = self.size_params.cgbo * vgb;
                    qgate = qgdo + qgso + qgb;
                    qbulk = -qgb;
                    qdrn = -qgdo;
                    qsrc = -(qgate + qbulk + qdrn);
                }

                gcddb = (newop.capbd + cgdo) * ag0;
                gcdsb = 0.0;
                gcsdb = 0.0;
                gcssb = (newop.capbs + cgso) * ag0;

                if self.model.rbodymod == 0 {
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
                if self.model.rgatemod == 3 {
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

                    let qgmb = self.size_params.cgbo * vgmb;
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

                    let qgb = self.size_params.cgbo * vgb;
                    qgate += qgdo + qgso + qgb;
                    qbulk -= qgb;
                    qdrn = -(qgate + qbulk + qsrc);
                }
                gcddb = (newop.capbd + cgdo - (newop.cgsb + newop.cbsb + newop.cdsb)) * ag0;
                gcdsb = -(newop.cgdb + newop.cbdb + newop.cddb) * ag0;

                gcsdb = newop.cdsb * ag0;
                gcssb = (newop.cddb + newop.capbs + cgso) * ag0;

                if self.model.rbodymod == 0 {
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
                let qcheq = newop.qchqs;
                let CoxWL = self.model_derived.coxe * self.size_params.weffCV * self.intp.nf * self.size_params.leffCV;
                let T0 = qdef * nqs_scaling_factor / CoxWL;
                ggtg = T0 * newop.gcrgg;
                newop.gtg = ggtg;
                ggts = T0 * newop.gcrgd;
                newop.gts = ggts;
                ggtd = T0 * newop.gcrgs;
                newop.gtd = ggtd;
                ggtb = T0 * newop.gcrgb;
                newop.gtb = ggtb;
                gqdef = nqs_scaling_factor * ag0;

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
                    let Css = newop.cddb;
                    let Cds = -(newop.cgdb + newop.cddb + newop.cbdb);
                    dsxpart_dVs = (Css - sxpart * (Css + Cds)) / qcheq;
                    let Csg = newop.cdgb;
                    let Cdg = -(newop.cggb + newop.cdgb + newop.cbgb);
                    dsxpart_dVg = (Csg - sxpart * (Csg + Cdg)) / qcheq;

                    let Csd = newop.cdsb;
                    let Cdd = -(newop.cgsb + newop.cdsb + newop.cbs);
                    dsxpart_dVd = (Csd - sxpart * (Csd + Cdd)) / qcheq;

                    dsxpart_dVb = -(dsxpart_dVd + dsxpart_dVg + dsxpart_dVs);
                }
                dxpart = 1.0 - sxpart;
                ddxpart_dVd = -dsxpart_dVd;
                ddxpart_dVg = -dsxpart_dVg;
                ddxpart_dVs = -dsxpart_dVs;
                ddxpart_dVb = -(ddxpart_dVd + ddxpart_dVg + ddxpart_dVs);

                if self.model.rgatemod == 3 {
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

                    let qgmb = self.size_params.cgbo * vgmb;
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

                    let qgb = self.size_params.cgbo * vgb;
                    qgate = qgdo + qgso + qgb;
                    qbulk = -qgb;
                    qdrn = -qgdo;
                    qsrc = -qgso;
                }

                gcddb = (newop.capbd + cgdo) * ag0;
                gcdsb = 0.0;
                gcsdb = 0.0;
                gcssb = (newop.capbs + cgso) * ag0;
                if self.model.rbodymod == 0 {
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

        newop.qg = qgate;
        newop.qd = qdrn - newop.qbd;
        newop.qs = qsrc - newop.qbs;
        if self.model.rgatemod == 3 {
            newop.qgmid = qgmid;
        }

        if self.model.rbodymod != 0 {
            newop.qb = qbulk + newop.qbd + newop.qbs;
        } else {
            newop.qb = qbulk;
        }
        // FIXME! Something got very lost in translation here
        // if let AnalysisInfo::TRAN(_, state) = an {
        //     // Transient, Do Numerical Integration
        //     if self.model.trnqsmod != 0 {
        //         newop.qcdump = qdef * nqs_scaling_factor;
        //         let (_g, i, _r) = state.integrate(newop.qcdump - self.op.qcdump, 0.0, 0.0, self.op.cqcdump);
        //         newop.cqcdump = i;
        //     }

        //     let (_g, i, _r) = state.integrate(newop.qb - self.op.qb, 0.0, 0.0, self.op.cqb);
        //     newop.cqb = i;
        //     let (_g, i, _r) = state.integrate(newop.qg - self.op.qg, 0.0, 0.0, self.op.cqg);
        //     newop.cqg = i;
        //     let (_g, i, _r) = state.integrate(newop.qd - self.op.qd, 0.0, 0.0, self.op.cqd);
        //     newop.cqd = i;

        //     if self.model.rgatemod == 3 {
        //         let (_g, i, _r) = state.integrate(newop.qgmid - self.op.qgmid, 0.0, 0.0, self.op.cqgmid);
        //         newop.cqgmid = i;
        //     }

        //     if self.model.rbodymod != 0 {
        //         let (_g, i, _r) = state.integrate(newop.qbs - self.op.qbs, 0.0, 0.0, self.op.cqbs);
        //         newop.cqbs = i;
        //         let (_g, i, _r) = state.integrate(newop.qbd - self.op.qbd, 0.0, 0.0, self.op.cqbd);
        //         newop.cqbd = i;
        //     }
        // }

        /* Calculate equivalent charge current */

        let cqgate = newop.cqg;
        let cqbody = newop.cqb;
        let cqdrn = newop.cqd;

        ceqqg = cqgate - gcggb * vgb + gcgdb * vbd + gcgsb * vbs;
        ceqqd = cqdrn - gcdgb * vgb - gcdgmb * vgmb + (gcddb + gcdbdb) * vbd - gcdbdb * vbd_jct + gcdsb * vbs;
        ceqqb = cqbody - gcbgb * vgb - gcbgmb * vgmb + gcbdb * vbd + gcbsb * vbs;

        newop.ceqqgmid = if self.model.rgatemod == 3 {
            newop.cqgmid + gcgmdb * vbd + gcgmsb * vbs - gcgmgmb * vgmb
        } else { // Default rgatemod==0
            0.0
        };

        if self.model.rbodymod != 0 {
            ceqqjs = newop.cqbs + gcsbsb * vbs_jct;
            ceqqjd = newop.cqbd + gcdbdb * vbd_jct;
        }

        if self.model.trnqsmod != 0 {
            let (_g, i, _r) = tran_state.integrate(newop.qcheq - self.op.qcheq, 0.0, 0.0, self.op.cqcheq);
            newop.cqcheq = i;

            let T0 = ggtg * vgb - ggtd * vbd - ggts * vbs;
            ceqqg += T0;
            let T1 = qdef * newop.gtau;
            ceqqd -= dxpart * T0 + T1 * (ddxpart_dVg * vgb - ddxpart_dVd * vbd - ddxpart_dVs * vbs);
            cqdef = newop.cqcdump - gqdef * qdef;
            cqcheq = newop.cqcheq - (gcqgb * vgb - gcqdb * vbd - gcqsb * vbs) + T0;
        }

        newop.gqdef = gqdef;
        newop.ggtg = ggtg;
        newop.ggtd = ggtd;
        newop.ggts = ggts;
        newop.ggtb = ggtb;
        newop.gcsbsb = gcsbsb;
        newop.gcsbsb = gcsbsb;
        newop.gcqsb = gcqsb;
        newop.gcqdb = gcqdb;
        newop.gcqgb = gcqgb;
        newop.cqcheq = cqcheq;
        newop.cqdef = cqdef;
        newop.gcqbb = gcqbb;
        newop.ceqqg = ceqqg;
        newop.ceqqd = ceqqd;
        newop.ceqqb = ceqqb;
        newop.ceqqjs = ceqqjs;
        newop.ceqqjd = ceqqjd;
        newop.gcgbb = gcgbb;
        newop.gcsbb = gcsbb;
        newop.gcssb = gcssb;
        newop.gcsgb = gcsgb;
        newop.gcsdb = gcsdb;
        newop.gcdbb = gcdbb;
        newop.gcggb = gcggb;
        newop.gcgdb = gcgdb;
        newop.gcgsb = gcgsb;
        newop.gcgsb = gcgsb;
        newop.gcgmgmb = gcgmgmb;
        newop.gcgmsb = gcgmsb;
        newop.gcgmdb = gcgmdb;
        newop.gcgmbb = gcgmbb;
        newop.gcdgmb = gcdgmb;
        newop.gcsgmb = gcsgmb;
        newop.ddxpart_dVd = ddxpart_dVd;
        newop.ddxpart_dVs = ddxpart_dVs;
        newop.ddxpart_dVb = ddxpart_dVb;
        newop.dsxpart_dVb = dsxpart_dVb;
        newop.dsxpart_dVs = dsxpart_dVs;
        newop.dsxpart_dVg = dsxpart_dVg;
        newop.dsxpart_dVd = dsxpart_dVd;
        newop.ddxpart_dVg = ddxpart_dVg;
        newop.gcgmdb = gcgmdb;
        newop.gcbgmb = gcbgmb;
        newop.gcddb = gcddb;
        newop.gcdgb = gcdgb;
        newop.gcdsb = gcdsb;
        newop.gcbdb = gcbdb;
        newop.gcbgb = gcbgb;
        newop.gcbsb = gcbsb;
        newop.gcbbb = gcbbb;
        newop.dxpart = dxpart;
        newop.sxpart = sxpart;
    }
}
