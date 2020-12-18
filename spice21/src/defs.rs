///
/// # Spice21 Circuit-Definitions Depots
///
use std::collections::HashMap;
use std::sync::{Arc, RwLock, RwLockReadGuard};

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
