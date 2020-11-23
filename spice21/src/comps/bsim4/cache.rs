///
/// Bsim4 Model, Instance, and Internal Parameter Registry and Cache
///
/// Combining and merging model and instance parameters is also performed here,
/// on-demand as combinations are requested.
/// The `ModelEntry` and `InstEntry` structs hold these combinations.
///
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use super::inst::{Bsim4InstSpecs, Bsim4InstVals};
use super::model::{Bsim4ModelSpecs, Bsim4ModelVals};
use super::{Bsim4InternalParams, Bsim4ModelDerivedParams, Bsim4SizeDepParams};
use crate::SpResult;

/// Entries of Derived Model Parameters
#[derive(Clone)]
pub(crate) struct Bsim4ModelEntry {
    pub(crate) vals: Bsim4ModelVals,
    pub(crate) derived: Bsim4ModelDerivedParams,
}
impl Bsim4ModelEntry {
    fn new(specs: &Bsim4ModelSpecs) -> Self {
        use super::bsim4derive::derive;
        use super::model::vals::resolve;

        let vals = resolve(specs);
        let derived = derive(&vals);
        Self { vals, derived }
    }
}

/// Entries of Derived Instance Parameters
#[derive(Clone)]
pub(crate) struct Bsim4InstEntry {
    pub(crate) intp: Bsim4InternalParams,
    pub(crate) size_params: Bsim4SizeDepParams,
}
impl Bsim4InstEntry {
    fn new(specs: &Bsim4InstSpecs, model: &Bsim4ModelEntry) -> Self {
        use super::bsim4inst::from;
        let (intp, size_params) = from(&model.vals, &model.derived, specs);
        Self { intp, size_params }
    }
}

/// Model, Instance, and Combination Registries
pub(crate) struct Bsim4Cache {
    pub(crate) models: HashMap<String, Bsim4ModelSpecs>,
    pub(crate) insts: HashMap<String, Bsim4InstSpecs>,
    cache: HashMap<(String, String), (Bsim4ModelEntry, Bsim4InstEntry)>,
}
impl Bsim4Cache {
    pub(crate) fn new() -> Self {
        Self {
            models: HashMap::new(),
            insts: HashMap::new(),
            cache: HashMap::new(),
        }
    }
    pub(crate) fn add_model(&mut self, name:&str, specs: Bsim4ModelSpecs) {
        self.models.insert(name.to_string(), specs);
    }
    pub(crate) fn add_inst(&mut self, inst: Bsim4InstSpecs) {
        self.insts.insert(inst.name.clone(), inst);
    }
    pub(crate) fn get(&mut self, model_name: &String, inst_name: &String) -> Option<(Bsim4ModelEntry, Bsim4InstEntry)> {
        if let Some(e) = self.cache.get(&(model_name.clone(), inst_name.clone())) {
            return Some(e.clone()); // FIXME: pointers
        }
        // Not in cache, create anew and insert 
        let model = self.models.get(model_name)?;
        let me = Bsim4ModelEntry::new(model);
        let inst = self.insts.get(inst_name)?;
        let ie = Bsim4InstEntry::new(inst, &me);
        self.cache.insert((model_name.clone(), inst_name.clone()), (me.clone(), ie.clone()));

        // FIXME: stop cloning, return references or pointers
        Some((me.clone(), ie.clone()))
    }
}
