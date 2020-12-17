///
/// # Spice21 Circuit-Definitions Depots
///
use std::collections::HashMap;
use std::sync::{Arc, RwLock};

/// The typical type of thread-shared pointers to Definitions
pub(crate) type DefPtr<T> = Arc<RwLock<T>>;
/// DefPtr "Constructor"
/// Note we cannot add `impl` methods, as `DefPtr` is externally defined.
pub(crate) fn defptr<T>(i: T) -> DefPtr<T> {
    Arc::new(RwLock::new(i))
}

// Collect up device-type-specific depots/ caches 
use crate::comps::{bsim4, diode, mos};

// The Module-Definition depot is defined here
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
        self.store.insert(m.name.clone(), defptr(m));
    } 
    pub(crate) fn get(&mut self, name: &str) -> Option<DefPtr<ModuleDef>> {
        match self.store.get(name) {
            Some(ptr) => Some(DefPtr::clone(ptr)),
            None => None
        } 
    }
}

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
// impl Defs {
//     pub fn new() -> Self {
//         Self {
//             modules: HashMap::new(),
//             mos0: HashMap::new(),
//             mos1: mos::Mos1Defs::default(),
//             bsim4: bsim4::Bsim4Cache::default(),
//             diodes: diode::DiodeDefs::default(),
//         }
//     }
// }
