///
/// # Spice21 Circuit-Definitions Depots
///
use std::collections::HashMap;
use std::sync::{Arc, RwLock, RwLockReadGuard};

use crate::analysis;

///
/// # Definition Pointer
///
/// The typical type of thread-shared pointers to Definitions
///
#[derive(Clone, Default)]
pub struct DefPtr<T>(Arc<RwLock<T>>);

impl<T> DefPtr<T> {
    /// Definition-Pointer Constructor
    pub fn new(i: T) -> Self {
        Self(Arc::new(RwLock::new(i)))
    }
    /// Read our definition.
    /// Panics if read fails.
    /// Typical usage requires a deref-and-re-ref on the result, e.g.
    /// `&*ptr.read()`
    /// That `&*` is ugly, and hopefully can eventually be pulled in here.
    pub fn read(&self) -> RwLockReadGuard<T> {
        self.0.read().unwrap()
    }
    pub fn clone(i: &Self) -> Self {
        Self(Arc::clone(&i.0))
    }
}

// The Module-Definition depot is defined here.
// All others (models, instance parameters, etc.)
// are imported from their comp-specific modules.

use crate::proto::Module as ModuleDef;
///
/// # Module Definitions Depot
///
#[derive(Default)]
pub struct ModuleDefs {
    pub store: HashMap<String, DefPtr<ModuleDef>>,
}
impl ModuleDefs {
    pub(crate) fn add(&mut self, m: ModuleDef) {
        self.store.insert(m.name.clone(), DefPtr::new(m));
    }
    pub(crate) fn get(&mut self, name: &str) -> Option<DefPtr<ModuleDef>> {
        match self.store.get(name) {
            Some(ptr) => Some(DefPtr::clone(ptr)),
            None => None,
        }
    }
}

//
// Shared Model & Instance Caching Infrastructure
//
// Each of the components with `model`/`instance` parameter separation uses
// a version of `ModelInstanceCache` to store its model-parameters, instance-parameters,
// and combinations thereof, typically including some derived information
// and simulator option-specifics.
//
// The `ModelInstanceCache` struct and `CacheEntry` trait comprise most
// code for searching these caches, and where necessary creating
// and adding the derived data.
//

pub trait CacheEntry: Clone {
    type Model;
    type Instance;
    fn new(model: &DefPtr<Self::Model>, inst: &DefPtr<Self::Instance>, opts: &analysis::Options) -> Self;
}

#[derive(Default)]
pub struct ModelInstanceCache<Model, Instance, Entry>
where
    Entry: CacheEntry<Model = Model, Instance = Instance>,
{
    pub(crate) models: HashMap<String, DefPtr<Model>>,
    pub(crate) insts: HashMap<String, DefPtr<Instance>>,
    pub(crate) cache: HashMap<(String, String), Entry>,
}
impl<Model, Instance, Entry> ModelInstanceCache<Model, Instance, Entry>
where
    Entry: CacheEntry<Model = Model, Instance = Instance>,
{
    pub(crate) fn add_model(&mut self, name: &str, model: Model) {
        self.models.insert(name.to_string(), DefPtr::new(model));
    }
    pub(crate) fn add_inst(&mut self, name: &str, inst: Instance) {
        self.insts.insert(name.to_string(), DefPtr::new(inst));
    }
    pub(crate) fn get(&mut self, inst: &str, model: &str, opts: &analysis::Options) -> Option<Entry> {
        // If we've already derived these parameters, clone a new pointer to them
        if let Some(e) = self.cache.get(&(inst.to_string(), model.to_string())) {
            return Some(e.clone());
        }

        // Not in cache, check whether we have definitions.
        let instptr = self.insts.get(inst)?;
        let modelptr = self.models.get(model)?;

        // If we get here, we found definitions of both instance and model params.
        // Now derive the internal ones, including any circuit options.
        let e = Entry::new(modelptr, instptr, opts);

        // Insert a copy in our cache, and return the original
        self.cache.insert((inst.to_string(), model.to_string()), e.clone());
        Some(e)
    }
}

// Collect up device-type-specific depots/ caches
use crate::comps::{bsim4, diode, mos};

///
/// # Definitions Struct
///
/// Central repository for top-level circuit definitions, including:
/// * Module definitions
/// * Models
/// * Instance parameter-sets
///
#[derive(Default)]
pub struct Defs {
    pub(crate) modules: ModuleDefs,
    pub(crate) mos0: HashMap<String, mos::MosType>,
    pub(crate) mos1: mos::Mos1Defs,
    pub(crate) bsim4: bsim4::Bsim4Cache,
    pub(crate) diodes: diode::DiodeDefs,
}
