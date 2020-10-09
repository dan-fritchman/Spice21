use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use super::inst::{Bsim4InstSpecs, Bsim4InstVals};
use super::model::{Bsim4ModelSpecs, Bsim4ModelVals};
use super::{Bsim4InternalParams, Bsim4ModelDerivedParams, Bsim4SizeDepParams};

#[derive(Clone)]
pub(crate) struct Bsim4ModelEntry {
    pub(crate) specs: Bsim4ModelSpecs,
    pub(crate) vals: Bsim4ModelVals,
    pub(crate) derived: Bsim4ModelDerivedParams,
    pub(crate) insts: Vec<Bsim4InstEntry>,
}

impl Bsim4ModelEntry {
    fn new(specs: Bsim4ModelSpecs) -> Self {
        use super::bsim4derive::derive;
        use super::model::vals::resolve;

        let vals = resolve(&specs);
        let derived = derive(&vals);
        Self {
            specs,
            vals,
            derived,
            insts: vec![],
        }
    }
}

#[derive(Clone)]
pub(crate) struct Bsim4InstEntry {
    pub(crate) specs: Bsim4InstSpecs,
    pub(crate) intp: Bsim4InternalParams,
    pub(crate) size_params: Bsim4SizeDepParams,
}

impl Bsim4InstEntry {
    fn new(specs: Bsim4InstSpecs, model: &Bsim4ModelEntry) -> Self {
        use super::bsim4inst::from;
        let (intp, size_params) = from(&model.vals, &model.derived, &specs);
        Self { specs, intp, size_params }
    }
}

pub(crate) struct Bsim4ModelCache(HashMap<String, Bsim4ModelEntry>);

impl Bsim4ModelCache {
    pub(crate) fn new() -> Self {
        Self(HashMap::new())
    }
    pub(crate) fn add<S: Into<String>>(&mut self, name: S, specs: Bsim4ModelSpecs) {
        let entry = Bsim4ModelEntry::new(specs);
        self.0.insert(name.into(), entry);
    }
    pub(crate) fn model(&mut self, model_name: &String) -> Option<&mut Bsim4ModelEntry> {
        self.0.get_mut(model_name)
    }
    pub(crate) fn inst(&mut self, model_name: &String, specs: Bsim4InstSpecs) -> Option<(Bsim4ModelEntry, Bsim4InstEntry)> {
        // FIXME: actually check whether these things are already in the cache!
        let model: &mut Bsim4ModelEntry = self.0.get_mut(model_name)?;
        let inst = Bsim4InstEntry::new(specs, &model);
        model.insts.push(inst.clone());
        // FIXME: stop cloning, return references or pointers
        // Some((&*model, &model.insts[model.insts.len() - 1]))
        Some((model.clone(), inst))
    }
}
