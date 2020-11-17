use super::{Bsim4, Bsim4OpPoint};
use crate::analysis::{Stamps, VarIndex};
use crate::comps::mos::MosType;
use crate::sparse21::{Eindex, Matrix};

impl Bsim4 {
    /// Convert operating-point into matrix stamps
    pub(crate) fn stamp(&self) -> Stamps<f64> {
        // Extract the operating-point object
        let newop = &self.guess;

        let mut ceqdrn: f64;
        let mut ceqbd: f64;
        let mut ceqbs: f64;
        let mut Gm: f64;
        let mut Gmbs: f64;
        let mut FwdSum: f64;
        let mut RevSum: f64;
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

        let mut ceqgcrg = 0.0;
        let mut gcrg = 0.0;
        let mut gcrgd = 0.0;
        let mut gcrgg = 0.0;
        let mut gcrgs = 0.0;
        let mut gcrgb = 0.0;

        if newop.mode >= 0 {
            Gm = newop.gm;
            Gmbs = newop.gmbs;
            FwdSum = Gm + Gmbs;
            RevSum = 0.0;

            ceqdrn = self.model.p() * (newop.cd - newop.gds * newop.vds - Gm * newop.vgs - Gmbs * newop.vbs);
            ceqbd = self.model.p()
                * (newop.csub + newop.Igidl
                    - (newop.gbds + newop.ggidld) * newop.vds
                    - (newop.gbgs + newop.ggidlg) * newop.vgs
                    - (newop.gbbs + newop.ggidlb) * newop.vbs);
            ceqbs = self.model.p() * (newop.Igisl + newop.ggisls * newop.vds - newop.ggislg * newop.vgd - newop.ggislb * newop.vbd);

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
                Istoteq = self.model.p() * (newop.Igs + newop.Igcs - gIstotg * newop.vgs - newop.gIgcsd * newop.vds - newop.gIgcsb * newop.vbs);

                gIdtotg = newop.gIgdg + newop.gIgcdg;
                gIdtotd = newop.gIgdd + newop.gIgcdd;
                gIdtots = newop.gIgcds;
                gIdtotb = newop.gIgcdb;
                Idtoteq = self.model.p()
                    * (newop.Igd + newop.Igcd
                        - newop.gIgdg * newop.vgd
                        - newop.gIgcdg * newop.vgs
                        - newop.gIgcdd * newop.vds
                        - newop.gIgcdb * newop.vbs);
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
                Ibtoteq = self.model.p() * (newop.Igb - newop.gIgbg * newop.vgs - newop.gIgbd * newop.vds - newop.gIgbb * newop.vbs);
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

            if self.model.rgatemod > 1 {
                let tmp = if self.model.rgatemod == 2 {
                    newop.vges - newop.vgs
                } else {
                    // rgatemod == 3
                    newop.vgms - newop.vgs
                };
                gcrgd = newop.gcrgd * tmp;
                gcrgg = newop.gcrgg * tmp;
                gcrgs = newop.gcrgs * tmp;
                gcrgb = newop.gcrgb * tmp;
                ceqgcrg = -(gcrgd * newop.vds + gcrgg * newop.vgs + gcrgb * newop.vbs);
                gcrgg -= newop.gcrg;
                gcrg = newop.gcrg;
            }
        } else {
            Gm = -newop.gm;
            Gmbs = -newop.gmbs;
            FwdSum = 0.0;
            RevSum = -(Gm + Gmbs);

            ceqdrn = -self.model.p() * (newop.cd + newop.gds * newop.vds + Gm * newop.vgd + Gmbs * newop.vbd);

            ceqbs = self.model.p()
                * (newop.csub + newop.Igisl + (newop.gbds + newop.ggisls) * newop.vds
                    - (newop.gbgs + newop.ggislg) * newop.vgd
                    - (newop.gbbs + newop.ggislb) * newop.vbd);
            ceqbd = self.model.p() * (newop.Igidl - newop.ggidld * newop.vds - newop.ggidlg * newop.vgs - newop.ggidlb * newop.vbs);

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
                    * (newop.Igs + newop.Igcd - newop.gIgsg * newop.vgs - newop.gIgcdg * newop.vgd + newop.gIgcdd * newop.vds
                        - newop.gIgcdb * newop.vbd);

                gIdtotg = newop.gIgdg + newop.gIgcsg;
                gIdtotd = newop.gIgdd + newop.gIgcss;
                gIdtots = newop.gIgcsd;
                gIdtotb = newop.gIgcsb;
                Idtoteq = self.model.p()
                    * (newop.Igd + newop.Igcs - (newop.gIgdg + newop.gIgcsg) * newop.vgd + newop.gIgcsd * newop.vds - newop.gIgcsb * newop.vbd);
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
                Ibtoteq = self.model.p() * (newop.Igb - newop.gIgbg * newop.vgd + newop.gIgbd * newop.vds - newop.gIgbb * newop.vbd);
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

            if self.model.rgatemod > 1 {
                let tmp = if self.model.rgatemod == 2 {
                    newop.vges - newop.vgs
                } else {
                    // rgatemod == 3
                    newop.vgms - newop.vgs
                };
                gcrgd = newop.gcrgs * tmp;
                gcrgg = newop.gcrgg * tmp;
                gcrgs = newop.gcrgd * tmp;
                gcrgb = newop.gcrgb * tmp;
                ceqgcrg = -(gcrgg * newop.vgd - gcrgs * newop.vds + gcrgb * newop.vbd);
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

        let mut gstot = 0.0;
        let mut gstotd = 0.0;
        let mut gstotg = 0.0;
        let mut gstots = 0.0;
        let mut gstotb = 0.0;
        let mut ceqgstot = 0.0;
        let mut gdtot = 0.0;
        let mut gdtotd = 0.0;
        let mut gdtotg = 0.0;
        let mut gdtots = 0.0;
        let mut gdtotb = 0.0;
        let mut ceqgdtot = 0.0;

        if self.model.rdsmod == 1 {
            ceqgstot = self.model.p() * (newop.gstotd * newop.vds + newop.gstotg * newop.vgs + newop.gstotb * newop.vbs);
            gstot = newop.gstot;
            gstotd = newop.gstotd;
            gstotg = newop.gstotg;
            gstots = newop.gstots - gstot;
            gstotb = newop.gstotb;

            ceqgdtot = -self.model.p() * (newop.gdtotd * newop.vds + newop.gdtotg * newop.vgs + newop.gdtotb * newop.vbs);
            gdtot = newop.gdtot;
            gdtotd = newop.gdtotd - gdtot;
            gdtotg = newop.gdtotg;
            gdtots = newop.gdtots;
            gdtotb = newop.gdtotb;
        }

        let (ceqjs, ceqjd) = match self.model.mos_type {
            MosType::NMOS => (newop.cbs - newop.gbs * newop.vbs_jct, newop.cbd - newop.gbd * newop.vbd_jct),
            MosType::PMOS => (-(newop.cbs - newop.gbs * newop.vbs_jct), -(newop.cbd - newop.gbd * newop.vbd_jct)),
        };

        let Bsim4OpPoint {
            mut ceqqg,
            mut ceqqd,
            mut ceqqb,
            mut cqdef,
            mut cqcheq,
            mut ceqqjs,
            mut ceqqjd,
            mut ceqqgmid,
            ..
        } = newop;

        if self.model.p() < 0.0 {
            ceqqg = -ceqqg;
            ceqqd = -ceqqd;
            ceqqb = -ceqqb;
            ceqgcrg = -ceqgcrg;

            if self.model.trnqsmod != 0 {
                cqdef = -cqdef;
                cqcheq = -cqcheq;
            }
            if self.model.rbodymod != 0 {
                ceqqjs = -ceqqjs;
                ceqqjd = -ceqqjd;
            }
            if self.model.rgatemod == 3 {
                ceqqgmid = -ceqqgmid;
            }
        }

        // Gather up RHS current-vector terms
        let mut b: Vec<(Option<VarIndex>, f64)> = vec![];

        b.push((self.ports.dNodePrime, (ceqjd - ceqbd + ceqgdtot - ceqdrn - ceqqd + Idtoteq)));
        b.push((self.ports.gNodePrime, -(ceqqg - ceqgcrg + Igtoteq)));

        if self.model.rgatemod == 2 {
            b.push((self.ports.gNodeExt, -ceqgcrg));
        } else if self.model.rgatemod == 3 {
            b.push((self.ports.gNodeMid, -(ceqqgmid + ceqgcrg)));
        }

        if self.model.rbodymod == 0 {
            b.push((self.ports.bNodePrime, (ceqbd + ceqbs - ceqjd - ceqjs - ceqqb + Ibtoteq)));
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
                (ceqdrn - ceqbs + ceqjs + ceqqd + ceqqg + ceqqb + ceqqjd + ceqqjs + ceqqgmid - ceqgstot + Istoteq),
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

        let (gjbd, gjbs) = if self.model.rbodymod != 0 { (newop.gbd, newop.gbs) } else { (0.0, 0.0) };
        let (gdpr, gspr) = if self.model.rdsmod != 0 {
            (self.intp.drainConductance, self.intp.sourceConductance)
        } else {
            (0.0, 0.0)
        };

        // FIXME: it appears all these gate currents go to "ground", both here and in the BSIM4 reference. Is that the idea?
        if self.model.rgatemod == 1 {
            let geltd = self.intp.grgeltd;
            j.push((self.matps.GEgePtr, geltd));
            j.push((self.matps.GPgePtr, -(geltd)));
            j.push((self.matps.GEgpPtr, -(geltd)));
            j.push((self.matps.GPgpPtr, newop.gcggb + geltd - newop.ggtg + gIgtotg));
            j.push((self.matps.GPdpPtr, newop.gcgdb - newop.ggtd + gIgtotd));
            j.push((self.matps.GPspPtr, newop.gcgsb - newop.ggts + gIgtots));
            j.push((self.matps.GPbpPtr, newop.gcgbb - newop.ggtb + gIgtotb));
        } else if self.model.rgatemod == 2 {
            j.push((self.matps.GEgePtr, gcrg));
            j.push((self.matps.GEgpPtr, gcrgg));
            j.push((self.matps.GEdpPtr, gcrgd));
            j.push((self.matps.GEspPtr, gcrgs));
            j.push((self.matps.GEbpPtr, gcrgb));

            j.push((self.matps.GPgePtr, -(gcrg)));
            j.push((self.matps.GPgpPtr, newop.gcggb - gcrgg - newop.ggtg + gIgtotg));
            j.push((self.matps.GPdpPtr, newop.gcgdb - gcrgd - newop.ggtd + gIgtotd));
            j.push((self.matps.GPspPtr, newop.gcgsb - gcrgs - newop.ggts + gIgtots));
            j.push((self.matps.GPbpPtr, newop.gcgbb - gcrgb - newop.ggtb + gIgtotb));
        } else if self.model.rgatemod == 3 {
            let geltd = self.intp.grgeltd;
            j.push((self.matps.GEgePtr, geltd));
            j.push((self.matps.GEgmPtr, -(geltd)));
            j.push((self.matps.GMgePtr, -(geltd)));
            j.push((self.matps.GMgmPtr, geltd + gcrg + newop.gcgmgmb));

            j.push((self.matps.GMdpPtr, gcrgd + newop.gcgmdb));
            j.push((self.matps.GMgpPtr, gcrgg));
            j.push((self.matps.GMspPtr, gcrgs + newop.gcgmsb));
            j.push((self.matps.GMbpPtr, gcrgb + newop.gcgmbb));

            j.push((self.matps.DPgmPtr, newop.gcdgmb));
            j.push((self.matps.GPgmPtr, -(gcrg)));
            j.push((self.matps.SPgmPtr, newop.gcsgmb));
            j.push((self.matps.BPgmPtr, newop.gcbgmb));

            j.push((self.matps.GPgpPtr, newop.gcggb - gcrgg - newop.ggtg + gIgtotg));
            j.push((self.matps.GPdpPtr, newop.gcgdb - gcrgd - newop.ggtd + gIgtotd));
            j.push((self.matps.GPspPtr, newop.gcgsb - gcrgs - newop.ggts + gIgtots));
            j.push((self.matps.GPbpPtr, newop.gcgbb - gcrgb - newop.ggtb + gIgtotb));
        } else {
            // Default case: rgatemode == 0
            j.push((self.matps.GPgpPtr, newop.gcggb - newop.ggtg + gIgtotg));
            j.push((self.matps.GPdpPtr, newop.gcgdb - newop.ggtd + gIgtotd));
            j.push((self.matps.GPspPtr, newop.gcgsb - newop.ggts + gIgtots));
            j.push((self.matps.GPbpPtr, newop.gcgbb - newop.ggtb + gIgtotb));
        }

        if self.model.rdsmod != 0 {
            j.push((self.matps.DgpPtr, gdtotg));
            j.push((self.matps.DspPtr, gdtots));
            j.push((self.matps.DbpPtr, gdtotb));
            j.push((self.matps.SdpPtr, gstotd));
            j.push((self.matps.SgpPtr, gstotg));
            j.push((self.matps.SbpPtr, gstotb));
        }
        {
            let tmp1 = newop.qdef * newop.gtau;
            j.push((
                self.matps.DPdpPtr,
                gdpr + newop.gds + newop.gbd + tmp1 * newop.ddxpart_dVd - gdtotd + RevSum + newop.gcddb + gbdpdp + newop.dxpart * newop.ggtd
                    - gIdtotd,
            ));
            j.push((self.matps.DPdPtr, -(gdpr + gdtot)));
            j.push((
                self.matps.DPgpPtr,
                Gm + newop.gcdgb - gdtotg + gbdpg - gIdtotg + newop.dxpart * newop.ggtg + tmp1 * newop.ddxpart_dVg,
            ));
            j.push((
                self.matps.DPspPtr,
                -(newop.gds + gdtots - newop.dxpart * newop.ggts + gIdtots - tmp1 * newop.ddxpart_dVs + FwdSum - newop.gcdsb - gbdpsp),
            ));
            j.push((
                self.matps.DPbpPtr,
                -(gjbd + gdtotb - Gmbs - newop.gcdbb - gbdpb + gIdtotb - tmp1 * newop.ddxpart_dVb - newop.dxpart * newop.ggtb),
            ));

            j.push((self.matps.DdpPtr, -(gdpr - gdtotd)));
            j.push((self.matps.DdPtr, gdpr + gdtot));

            j.push((
                self.matps.SPdpPtr,
                -(newop.gds + gstotd + RevSum - newop.gcsdb - gbspdp - tmp1 * newop.dsxpart_dVd - newop.sxpart * newop.ggtd + gIstotd),
            ));
            j.push((
                self.matps.SPgpPtr,
                newop.gcsgb - Gm - gstotg + gbspg + newop.sxpart * newop.ggtg + tmp1 * newop.dsxpart_dVg - gIstotg,
            ));
            j.push((
                self.matps.SPspPtr,
                gspr + newop.gds + newop.gbs + tmp1 * newop.dsxpart_dVs - gstots + FwdSum + newop.gcssb + gbspsp + newop.sxpart * newop.ggts
                    - gIstots,
            ));
            j.push((self.matps.SPsPtr, -(gspr + gstot)));
            j.push((
                self.matps.SPbpPtr,
                -(gjbs + gstotb + Gmbs - newop.gcsbb - gbspb - newop.sxpart * newop.ggtb - tmp1 * newop.dsxpart_dVb + gIstotb),
            ));
        }

        j.push((self.matps.SspPtr, -(gspr - gstots)));
        j.push((self.matps.SsPtr, gspr + gstot));

        j.push((self.matps.BPdpPtr, newop.gcbdb - gjbd + gbbdp - gIbtotd));
        j.push((self.matps.BPgpPtr, newop.gcbgb - newop.gbgs - gIbtotg));
        j.push((self.matps.BPspPtr, newop.gcbsb - gjbs + gbbsp - gIbtots));
        j.push((self.matps.BPbpPtr, gjbd + gjbs + newop.gcbbb - newop.gbbs - gIbtotb));
        {
            // GIDL & GISL Section
            let ggidld = newop.ggidld;
            let ggidlg = newop.ggidlg;
            let ggidlb = newop.ggidlb;
            let ggislg = newop.ggislg;
            let ggisls = newop.ggisls;
            let ggislb = newop.ggislb;

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
        }
        // FIXME: are gbs, gd only balanced out for rbodymod != 0 ?
        if self.model.rbodymod != 0 {
            j.push((self.matps.DPdbPtr, newop.gcdbdb - newop.gbd));
            j.push((self.matps.SPsbPtr, -(newop.gbs - newop.gcsbsb)));

            j.push((self.matps.DBdpPtr, newop.gcdbdb - newop.gbd));
            j.push((self.matps.DBdbPtr, newop.gbd - newop.gcdbdb + self.intp.grbpd + self.intp.grbdb));
            j.push((self.matps.DBbpPtr, -(self.intp.grbpd)));
            j.push((self.matps.DBbPtr, -(self.intp.grbdb)));

            j.push((self.matps.BPdbPtr, -(self.intp.grbpd)));
            j.push((self.matps.BPbPtr, -(self.intp.grbpb)));
            j.push((self.matps.BPsbPtr, -(self.intp.grbps)));
            j.push((self.matps.BPbpPtr, self.intp.grbpd + self.intp.grbps + self.intp.grbpb));

            j.push((self.matps.SBspPtr, newop.gcsbsb - newop.gbs));
            j.push((self.matps.SBbpPtr, -(self.intp.grbps)));
            j.push((self.matps.SBbPtr, -(self.intp.grbsb)));
            j.push((self.matps.SBsbPtr, newop.gbs - newop.gcsbsb + self.intp.grbps + self.intp.grbsb));

            j.push((self.matps.BdbPtr, -(self.intp.grbdb)));
            j.push((self.matps.BbpPtr, -(self.intp.grbpb)));
            j.push((self.matps.BsbPtr, -(self.intp.grbsb)));
            j.push((self.matps.BbPtr, self.intp.grbsb + self.intp.grbdb + self.intp.grbpb));
        }

        if self.model.trnqsmod != 0 {
            j.push((self.matps.QqPtr, newop.gqdef + newop.gtau));
            j.push((self.matps.QgpPtr, newop.ggtg - newop.gcqgb));
            j.push((self.matps.QdpPtr, newop.ggtd - newop.gcqdb));
            j.push((self.matps.QspPtr, newop.ggts - newop.gcqsb));
            j.push((self.matps.QbpPtr, newop.ggtb - newop.gcqbb));

            j.push((self.matps.DPqPtr, newop.dxpart * newop.gtau));
            j.push((self.matps.SPqPtr, newop.sxpart * newop.gtau));
            j.push((self.matps.GPqPtr, -(newop.gtau)));
        }
        // And return our matrix stamps
        return Stamps { g: j, b };
    }
}
