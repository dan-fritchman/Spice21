pub(crate) mod vals;
pub(crate) use super::bsim4defs::Bsim4ModelSpecs;
pub(crate) use super::bsim4defs::Bsim4ModelVals;
pub(crate) use crate::proto::Bsim4Model as Bsim4ModelProto;
pub(crate) use vals::*;

use crate::comps::mos::MosType;

impl Bsim4ModelSpecs {
    /// Create a new Bsim4Model of MosType `t`.
    /// All other parameters take on their default values
    pub(crate) fn new (t: MosType) -> Self {
        Self {
            mos_type: Some(t), 
            ..Default::default()
        }
    }
    pub(crate) fn from(proto: &Bsim4ModelProto) -> Self {
        let mos_type = match proto.mos_type {
            1 => MosType::PMOS,
            _ => MosType::NMOS,
        };
        Self {
            mos_type:Some(mos_type),
            ..Default::default()
        }
    }
}

impl Bsim4ModelVals {
    /// Polarity function
    pub(crate) fn p(&self) -> f64 {
        self.mos_type.p()
    }
}
