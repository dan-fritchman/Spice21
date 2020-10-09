use super::bsim4defs::*;
use super::Bsim4InternalParams;
use crate::analysis::{Solver, VarIndex, VarKind};
use crate::comps::mos::{MosTerm, MosPorts};
use crate::SpNum;
use serde::{Deserialize, Serialize};

#[derive(Default, Debug, Deserialize, Serialize)]
pub(crate) struct Bsim4Ports<T> {
    pub(crate) dNode: T,
    pub(crate) dNodePrime: T,
    pub(crate) sNode: T,
    pub(crate) sNodePrime: T,
    pub(crate) gNodeExt: T,
    pub(crate) gNodePrime: T,
    pub(crate) gNodeMid: T,
    pub(crate) bNode: T,
    pub(crate) bNodePrime: T,
    pub(crate) dbNode: T,
    pub(crate) sbNode: T,
    pub(crate) qNode: T,
}

impl Bsim4Ports<Option<VarIndex>> {
    pub(crate) fn from<T: SpNum>(
        iname: String,
        terms: &MosPorts<Option<VarIndex>>,
        model: &Bsim4ModelVals,
        intp: &Bsim4InternalParams,
        solver: &mut Solver<T>,
    ) -> Self {
        use MosTerm::{B, D, G, S};

        let dNode = terms[D];
        let sNode = terms[S];
        let gNodeExt = terms[G];
        let bNode = terms[B];

        // Terminal resistances
        // FIXME: can potentially filter our the noise-mode depenedence only when doing noise analysis
        let dNodePrime = if model.rdsmod != 0 || (model.tnoimod == 1) {
            let mut name = iname.clone();
            name.push_str("_drain");
            Some(solver.vars.add(name, VarKind::V))
        } else {
            dNode
        };
        let sNodePrime = if model.rdsmod != 0 || (model.tnoimod == 1) {
            let mut name = iname.clone();
            name.push_str("_source");
            Some(solver.vars.add(name, VarKind::V))
        } else {
            sNode
        };

        // Gate Resistance
        let gNodePrime = if intp.rgatemod > 0 {
            let mut name = iname.clone();
            name.push_str("_gate");
            Some(solver.vars.add(name, VarKind::V))
        } else {
            gNodeExt
        };

        let gNodeMid = if intp.rgatemod == 3 {
            let mut name = iname.clone();
            name.push_str("_midgate");
            Some(solver.vars.add(name, VarKind::V))
        } else {
            gNodeExt
        };

        /* internal body nodes for body resistance model */
        let (dbNode, bNodePrime, sbNode) = if intp.rbodymod == 1 || intp.rbodymod == 2 {
            let mut name = iname.clone();
            name.push_str("_dbody");
            let dbNode = Some(solver.vars.add(name, VarKind::V));

            let mut name = iname.clone();
            name.push_str("_body");
            let bNodePrime = Some(solver.vars.add(name, VarKind::V));

            let mut name = iname.clone();
            name.push_str("_sbody");
            let sbNode = Some(solver.vars.add(name, VarKind::V));

            (dbNode, bNodePrime, sbNode)
        } else {
            (bNode, bNode, bNode)
        };

        // NQS charge node
        let qNode = if intp.trnqsmod != 0 {
            let mut name = iname.clone();
            name.push_str("_charge");
            Some(solver.vars.add(name, VarKind::Q))
        } else {
            None
        };

        Bsim4Ports {
            dNode,
            dNodePrime,
            sNode,
            sNodePrime,
            gNodeExt,
            gNodePrime,
            gNodeMid,
            bNode,
            bNodePrime,
            dbNode,
            sbNode,
            qNode,
        }
    }
}
