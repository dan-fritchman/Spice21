//!
//! # Spice21 Circuit
//!
//! Roughly the "middle-end" circuit description.
//!
//! Typically created from either:
//! * (a) Protobuf via `spice21::proto`, or
//! * (b) Directly for testing or in Rust use-cases
//!

use super::comps::{DiodeInstParams, DiodeModel};
use super::comps::{Mos1InstanceParams, Mos1Model, MosType};
use super::proto::instance::Comp as CompProto;
use super::proto::Circuit as CircuitProto;
use crate::SpResult;

use crate::comps::bsim4::Bsim4InstSpecs;
use crate::comps::mos::MosPorts;

use crate::comps::bsim4::Bsim4ModelCache;

/// Node Reference
#[derive(Debug, Clone)]
pub enum NodeRef {
    Gnd,
    Num(usize),
    Name(String),
}
/// Conversion to create Nodes from string-refs
impl From<&str> for NodeRef {
    fn from(f: &str) -> Self { n(f) }
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

pub struct Bsim4i {
    pub(crate) name: String,
    pub(crate) ports: MosPorts<NodeRef>,
    pub(crate) model: String,
    pub(crate) params: Bsim4InstSpecs,
}

pub struct Mos0i {
    pub(crate) name: String,
    pub(crate) mos_type: MosType,
    pub(crate) ports: MosPorts<NodeRef>,
}
pub struct Mos1i {
    pub(crate) name: String,
    pub(crate) model: Mos1Model,
    pub(crate) params: Mos1InstanceParams,
    pub(crate) ports: MosPorts<NodeRef>,
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
    pub fn from(c: CompProto) -> Self {
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
            CompProto::M(m) => Comp::Mos1(Mos1i { name: m.name, model:Mos1Model::default(), params:Mos1InstanceParams::default(), ports: MosPorts { g: n(m.g), d: n(m.d), s: n(m.s), b: n(m.b) } }),
        }
    }
}

impl From<Ds> for Comp {
    fn from(d: Ds) -> Self {
        Comp::D(d)
    }
}

pub struct ModelCache {
    pub(crate) bsim4: Bsim4ModelCache,
}
impl ModelCache {
    pub(crate) fn new() -> Self {
        Self {
            bsim4: Bsim4ModelCache::new(),
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
    pub(crate) fn from_comps(comps: Vec<Comp>) -> Self {
        Self {
            comps: comps,
            models: ModelCache::new(),
        }
    }
    /// Create from a protobuf-generated circuit
    pub fn from(c: CircuitProto) -> Self {
        let CircuitProto { name, comps, .. } = c;

        let mut cs: Vec<Comp> = vec![];
        for opt in comps.into_iter() {
            if let Some(c) = opt.comp {
                cs.push(Comp::from(c));
            }
        }
        Self {
            comps: cs,
            models: ModelCache::new(),
        }
    }
    /// Decode from bytes
    pub fn decode(bytes_: &[u8]) -> SpResult<Self> {
        use prost::Message;
        use std::io::Cursor;

        // Decode the protobuf version
        let c = CircuitProto::decode(&mut Cursor::new(bytes_));
        // Unfortunately these conversion errors don't convert to our Result
        let ckt_proto = match c {
            Ok(ckt) => ckt,
            Err(e) => panic!("Circuit Decode Failed"),
        };
        // And convert into a Circuit
        let c = Self::from(ckt_proto);
        Ok(c)
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
