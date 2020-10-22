//!
//! # Spice21 Circuit
//!
//! Roughly the "middle-end" circuit description.
//!
//! Typically created from either:
//! * (a) Protobuf via `spice21::proto`, or
//! * (b) Directly for testing or in Rust use-cases
//!

use std::collections::HashMap;

use super::comps::{DiodeInstParams, DiodeModel};
use super::comps::{Mos1InstanceParams, Mos1Model, MosType};
use super::proto::def::Defines as DefProto;
use super::proto::instance::Comp as CompProto;
use super::proto::Circuit as CircuitProto;
use crate::comps::bsim4::Bsim4Cache;
use crate::comps::mos::MosPorts;
use crate::{SpError, SpResult};

/// Node Reference
#[derive(Debug, Clone)]
pub enum NodeRef {
    Gnd,
    Num(usize),
    Name(String),
}
use NodeRef::Gnd;

/// Conversion to create Nodes from string-refs
impl From<&str> for NodeRef {
    fn from(f: &str) -> Self {
        n(f)
    }
}
/// Create a Node from anything convertible into String
/// Empty string is a cardinal value for creating Gnd
pub(crate) fn n<S: Into<String>>(name: S) -> NodeRef {
    let s: String = name.into();
    if s.len() == 0 {
        NodeRef::Gnd
    } else {
        NodeRef::Name(s)
    }
}

/// Convert anything convertible into String
pub(crate) fn s<S: Into<String>>(from: S) -> String {
    from.into()
}

/// Voltage Source
pub struct Vs {
    pub name: String,
    pub vdc: f64,
    pub acm: f64,
    pub p: NodeRef,
    pub n: NodeRef,
}

impl From<Bsim4i> for Comp {
    fn from(x: Bsim4i) -> Self {
        Comp::Bsim4(x)
    }
}

impl From<Vs> for Comp {
    fn from(x: Vs) -> Self {
        Comp::V(x)
    }
}

/// Diode
pub struct Ds {
    pub name: String,
    pub model: DiodeModel,
    pub inst: DiodeInstParams,
    pub p: NodeRef,
    pub n: NodeRef,
}

impl Ds {
    pub fn new<S: Into<String>>(name: S, p: NodeRef, n: NodeRef) -> Self {
        Self {
            p,
            n,
            name: name.into(),
            inst: DiodeInstParams::default(),
            model: DiodeModel::default(),
        }
    }
}

pub struct Mos0i {
    pub(crate) name: String,
    pub(crate) mos_type: MosType,
    pub(crate) ports: MosPorts<NodeRef>,
}
pub struct Mos1i {
    pub(crate) name: String,
    pub(crate) ports: MosPorts<NodeRef>,
    pub(crate) model: String,
    pub(crate) params: String,
}
pub struct Bsim4i {
    pub(crate) name: String,
    pub(crate) ports: MosPorts<NodeRef>,
    pub(crate) model: String,
    pub(crate) params: String,
}

/// Component Enum.
/// Circuits are mostly a list of these variants.
pub enum Comp {
    V(Vs),
    I(f64, NodeRef, NodeRef),
    R(f64, NodeRef, NodeRef),
    C(f64, NodeRef, NodeRef),
    D(Ds),
    Mos0(Mos0i),
    Mos1(Mos1i),
    Bsim4(Bsim4i),
}

impl Comp {
    /// Replacement for deprecated original `V` enum variant
    pub fn vdc<S: Into<String>>(name: S, vdc: f64, p: NodeRef, n: NodeRef) -> Comp {
        Comp::V(Vs {
            name: name.into(),
            vdc,
            acm: 0.0,
            p,
            n,
        })
    }
    // Convert from protobuf-generated classes
    pub fn from(c: CompProto, models: &ModelCache) -> Self {
        match c {
            CompProto::I(i) => Comp::I(i.dc, n(i.p), n(i.n)),
            CompProto::R(r) => Comp::R(r.g, n(r.p), n(r.n)),
            CompProto::D(d) => {
                let d1 = Ds::new(d.name, n(d.p), n(d.n));
                Comp::D(d1)
            }
            CompProto::V(v) => {
                let vs = Vs {
                    name: v.name,
                    p: n(v.p),
                    n: n(v.n),
                    vdc: v.dc,
                    acm: v.acm,
                };
                Comp::V(vs)
            }
            CompProto::C(c) => Comp::C(c.c, n(c.p), n(c.n)),
            CompProto::M(m) => {
                let ports: MosPorts<NodeRef> = match m.ports {
                    Some(p) => MosPorts {
                        g: n(p.g),
                        d: n(p.d),
                        s: n(p.s),
                        b: n(p.b),
                    },
                    None => MosPorts {
                        g: Gnd,
                        d: Gnd,
                        s: Gnd,
                        b: Gnd,
                    }, // FIXME: whether to default or not
                };
                // Mos instances break out to their solver-types here
                if let Some(_model) = models.bsim4.models.get(&m.model) {
                    Comp::Bsim4(Bsim4i {
                        name: m.name,
                        model: m.model.clone(),
                        params: m.params.clone(),
                        ports,
                    })
                } else {
                    Comp::Mos1(Mos1i {
                        name: m.name,
                        model: m.model.clone(),
                        params: m.params.clone(),
                        ports,
                    })
                }
            }
        }
    }
}

impl From<Ds> for Comp {
    fn from(d: Ds) -> Self {
        Comp::D(d)
    }
}

use crate::proto::Mos1InstParams as Mos1InstSpecs;
use crate::proto::Mos1Model as Mos1ModelSpecs;
#[derive(Default)]
pub struct Mos1Defs {
    pub(crate) models: HashMap<String, Mos1ModelSpecs>,
    pub(crate) insts: HashMap<String, Mos1InstSpecs>,
}
pub struct ModelCache {
    pub(crate) mos1: Mos1Defs,
    pub(crate) bsim4: Bsim4Cache,
}
impl ModelCache {
    pub(crate) fn new() -> Self {
        Self {
            mos1: Mos1Defs::default(),
            bsim4: Bsim4Cache::new(),
        }
    }
}

/// Primary Circuit Structure
pub struct Ckt {
    pub comps: Vec<Comp>,
    pub models: ModelCache,
}

impl Ckt {
    /// Create a new, empty Circuit
    pub fn new() -> Self {
        Self {
            comps: vec![],
            models: ModelCache::new(),
        }
    }
    pub fn from_comps(comps: Vec<Comp>) -> Self {
        Self {
            comps: comps,
            models: ModelCache::new(),
        }
    }
    /// Create from a protobuf-generated circuit
    pub fn from(c: CircuitProto) -> SpResult<Self> {
        let CircuitProto { comps, defs, .. } = c;
        let mut models = ModelCache::new();

        // Step through all definitions
        for def in defs.into_iter() {
            match def.defines.unwrap() {
                DefProto::Bsim4inst(x) => {
                    models.bsim4.add_inst(x);
                }
                DefProto::Bsim4model(x) => {
                    models.bsim4.add_model(x);
                }
                DefProto::Mos1model(x) => {
                    models.mos1.models.insert(x.name.clone(), x);
                }
                DefProto::Mos1inst(x) => {
                    models.mos1.insts.insert(x.name.clone(), x);
                }
                // DefProto::Subckt(_x),
                // DefProto::Lib(_x),
                // DefProto::Diodemodel(_x),
                _ => {
                    return Err(SpError::new("Unsupported Definition"));
                }
            }
        }
        // And step through all instances
        let mut cs: Vec<Comp> = vec![];
        for opt in comps.into_iter() {
            if let Some(c) = opt.comp {
                cs.push(Comp::from(c, &models));
            }
        }
        Ok(Self { comps: cs, models })
    }
    /// Decode from bytes, via proto definitions
    pub fn decode(bytes_: &[u8]) -> SpResult<Self> {
        use prost::Message;
        use std::io::Cursor;

        // Decode the protobuf version
        let ckt_proto = CircuitProto::decode(&mut Cursor::new(bytes_))?;
        // And convert into a Circuit
        Self::from(ckt_proto)
    }
    /// Add anything convertible into `Comp`,
    /// typically the enum-associated structs `Vsrc` et al.
    pub fn add<C: Into<Comp>>(&mut self, comp: C) {
        self.comps.push(comp.into());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spresult::TestResult;

    #[test]
    fn test_ckt_parse() -> TestResult {
        Ckt {
            comps: vec![Comp::I(1e-3, NodeRef::Num(0), NodeRef::Gnd), Comp::R(1e-3, NodeRef::Num(0), NodeRef::Gnd)],
            models: ModelCache::new(),
        };
        Ok(())
    }
}
