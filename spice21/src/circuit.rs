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

// use crate::comps::mos::MosTerminals;
// use crate::comps::bsim4::{Bsim4ModelSpecs, Bsim4InstSpecs};

/// Node Reference
#[derive(Debug, Clone)]
pub enum NodeRef {
    Gnd,
    Num(usize),
    Name(String),
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

// impl From<Bsim4i> for Comp {
//     fn from(x: Bsim4i) -> Self {
//         Comp::Bsim4(x)
//     }
// }

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


// pub struct Bsim4i {
//     pub(crate) name: String,
//     pub(crate) ports: MosTerminals<NodeRef>,
//     pub(crate) model: Bsim4ModelSpecs,
//     pub(crate) params: Bsim4InstSpecs,
// }
/// Component Enum.
/// Circuits are mostly a list of these variants.
pub enum Comp {
    V(Vs),
    I(f64, NodeRef, NodeRef),
    R(f64, NodeRef, NodeRef),
    C(f64, NodeRef, NodeRef),
    D(Ds),

    Mos0(MosType, NodeRef, NodeRef, NodeRef, NodeRef),
    Mos1(
        Mos1Model,
        Mos1InstanceParams,
        NodeRef,
        NodeRef,
        NodeRef,
        NodeRef,
    ),
    // Bsim4(Bsim4i),
}

impl Comp {
    /// Replacement for deprecated original `V` enum variant
    pub fn vdc(vdc: f64, p: NodeRef, n: NodeRef) -> Comp {
        Comp::V(Vs {
            name: s("tbd"),
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
            CompProto::Mos(m) => Comp::Mos1(
                Mos1Model::default(),
                Mos1InstanceParams::default(),
                n(m.g),
                n(m.d),
                n(m.s),
                n(m.b),
            ),
        }
    }
}

impl From<Ds> for Comp {
    fn from(d: Ds) -> Self {
        Comp::D(d)
    }
}

/// Primary Circuit Structure
pub struct Ckt {
    pub comps: Vec<Comp>,
}

impl Ckt {
    /// Create a new, empty Circuit
    pub fn new() -> Self {
        Self { comps: vec![] }
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
        Self { comps: cs }
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
            comps: vec![
                Comp::I(1e-3, NodeRef::Num(0), NodeRef::Gnd),
                Comp::R(1e-3, NodeRef::Num(0), NodeRef::Gnd),
            ],
        };
        Ok(())
    }
}
