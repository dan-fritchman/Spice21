pub(crate) mod vals;
pub(crate) use super::bsim4defs::Bsim4ModelVals;
pub(crate) use crate::proto::Bsim4Model as Bsim4ModelSpecs;
pub(crate) use vals::*;

use crate::comps::mos::MosType;

impl Bsim4ModelSpecs {
    /// Create a new Bsim4Model of MosType `t`.
    /// All other parameters take on their default values
    pub(crate) fn new(t: MosType) -> Self {
        Self {
            mos_type: t as i32, // Convert into protobuf's enum repr, which is just i32
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
