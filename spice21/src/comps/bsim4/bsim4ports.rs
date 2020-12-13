use super::{Bsim4InstVals, Bsim4ModelVals};
use super::Bsim4InternalParams;
use crate::analysis::{Variables, VarIndex, VarKind};
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
    /// Create Bsim4Ports from all the required device info
    /// * Terminal connections `terms` include the primary/ external G,D,S,B nodes
    /// * Model dictates whether several internal nodes are created
    /// * Internal params `intp` dictates several more
    /// * If internal variables are to be created, they are added to mutable Variables `vars`. 
    /// Instance path is required for generation of internal node-names. 
    pub(crate) fn from<T: SpNum>(
        path: String,
        terms: &MosPorts<Option<VarIndex>>,
        model: &Bsim4ModelVals,
        intp: &Bsim4InternalParams,
        vars: &mut Variables<T>,
    ) -> Self {
        use MosTerm::{B, D, G, S};

        let dNode = terms[D];
        let sNode = terms[S];
        let gNodeExt = terms[G];
        let bNode = terms[B];

        // Terminal resistances
        // FIXME: can potentially filter our the noise-mode dependence only when doing noise analysis
        let dNodePrime = if model.rdsmod != 0 || (model.tnoimod == 1) {
            let name = format!("{}.{}", path, "drain");
            Some(vars.add(name, VarKind::V))
        } else {
            dNode
        };
        let sNodePrime = if model.rdsmod != 0 || (model.tnoimod == 1) {
            let name = format!("{}.{}", path, "source");
            Some(vars.add(name, VarKind::V))
        } else {
            sNode
        };
        // Gate Resistance
        let gNodePrime = if intp.rgatemod > 0 {
            let name = format!("{}.{}", path, "gate");
            Some(vars.add(name, VarKind::V))
        } else {
            gNodeExt
        };

        let gNodeMid = if intp.rgatemod == 3 {
            let name = format!("{}.{}", path, "midgate");
            Some(vars.add(name, VarKind::V))
        } else {
            gNodeExt
        };

        /* internal body nodes for body resistance model */
        let (dbNode, bNodePrime, sbNode) = if intp.rbodymod == 1 || intp.rbodymod == 2 {
            let name = format!("{}.{}", path, "dbody");
            let dbNode = Some(vars.add(name, VarKind::V));

            let name = format!("{}.{}", path, "body");
            let bNodePrime = Some(vars.add(name, VarKind::V));

            let name = format!("{}.{}", path, "sbody");
            let sbNode = Some(vars.add(name, VarKind::V));

            (dbNode, bNodePrime, sbNode)
        } else {
            (bNode, bNode, bNode)
        };

        // NQS charge node
        let qNode = if intp.trnqsmod != 0 {
            let name = format!("{}.{}", path, "charge");
            Some(vars.add(name, VarKind::Q))
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
