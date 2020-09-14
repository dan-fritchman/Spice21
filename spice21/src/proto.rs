use super::comps::{Mos1InstanceParams, Mos1Model, MosType};
use crate::SpResult;
use prost::Message;

pub mod proto {
    // Include the prost-expanded proto-file content
    include!(concat!(env!("OUT_DIR"), "/spice21.proto.rs"));

    // fn some() -> Circuit {
    //     use super::s;
    //     use instance::Comp::{Mos, C, D, I, R, V};
    //     return Circuit {
    //         name: String::from("tbd"),
    //         comps: vec![Instance {
    //             comp: Some(D(Diode {
    //                 name: Some(s("???")),
    //                 p: Some(s("p")),
    //                 n: Some(s("N")),
    //                 model: Some(DiodeModel::default()),
    //                 params: Some(DiodeInstParams::default()),
    //             })),
    //         }],
    //         defs: vec![],
    //     };
    // }
}

#[derive(Debug, Clone)]
pub enum NodeRef {
    Gnd,
    Num(usize),
    Name(String),
}

/// Create a Node from anything convertible into String
/// Empty string is a cardinal value for creating Gnd
pub fn n<S: Into<String>>(name: S) -> NodeRef {
    let s: String = name.into();
    if s.len() == 0 {
        NodeRef::Gnd
    } else {
        NodeRef::Name(s)
    }
}

/// Convert anything convertible into String
pub fn s<S: Into<String>>(from: S) -> String {
    from.into()
}

pub struct Vs {
    pub name: String,
    pub vdc: f64,
    pub acm: f64,
    pub p: NodeRef,
    pub n: NodeRef,
}

impl From<Vs> for CompParse {
    fn from(x: Vs) -> Self {
        CompParse::Vb(x)
    }
}

use super::comps::{DiodeInstParams, DiodeModel};

pub struct D1 {
    pub name: String,
    pub model: DiodeModel,
    pub inst: DiodeInstParams,
    pub p: NodeRef,
    pub n: NodeRef,
}

impl D1 {
    pub fn new<S: Into<String>>(name: S, p: NodeRef, n: NodeRef) -> Self {
        D1 {
            p,
            n,
            name: name.into(),
            inst: DiodeInstParams::default(),
            model: DiodeModel::default(),
        }
    }
}

pub enum CompParse {
    Vb(Vs),
    I(f64, NodeRef, NodeRef),
    R(f64, NodeRef, NodeRef),
    C(f64, NodeRef, NodeRef),
    D1(D1),
    Mos0(MosType, NodeRef, NodeRef, NodeRef, NodeRef),
    Mos1(
        Mos1Model,
        Mos1InstanceParams,
        NodeRef,
        NodeRef,
        NodeRef,
        NodeRef,
    ),
}

use proto::instance::Comp;

impl CompParse {
    /// Replacement for deprecated `V` enum variant
    pub fn V(vdc: f64, p: NodeRef, n: NodeRef) -> CompParse {
        CompParse::Vb(Vs {
            name: s("tbd"),
            vdc,
            acm: 0.0,
            p,
            n,
        })
    }
    // Convert from protobuf-generated classes
    pub fn from(c: Comp) -> Self {
        match c {
            Comp::I(i) => CompParse::I(i.dc, n(i.p), n(i.n)),
            Comp::R(r) => CompParse::R(r.g, n(r.p), n(r.n)),
            Comp::D(d) => {
                let d1 = D1::new(d.name, n(d.p), n(d.n));
                CompParse::D1(d1)
            }
            Comp::V(v) => {
                let vs = Vs {
                    name: v.name,
                    p: n(v.p),
                    n: n(v.n),
                    vdc: v.dc,
                    acm: v.acm,
                };
                CompParse::Vb(vs)
            }
            Comp::C(c) => CompParse::C(c.c, n(c.p), n(c.n)),
            Comp::Mos(m) => CompParse::Mos1(
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

impl From<D1> for CompParse {
    fn from(d: D1) -> Self {
        CompParse::D1(d)
    }
}

pub struct CktParse {
    pub nodes: usize,
    pub comps: Vec<CompParse>,
}

impl CktParse {
    pub fn new() -> Self {
        Self {
            nodes: 0,
            comps: vec![],
        }
    }
    // Create from a protobuf-generated circuit
    pub fn from(c: proto::Circuit) -> Self {
        let proto::Circuit { name, comps, .. } = c;

        let mut cs: Vec<CompParse> = vec![];
        for opt in comps.into_iter() {
            if let Some(c) = opt.comp {
                cs.push(CompParse::from(c));
            }
        }
        Self {
            nodes: 0,
            comps: cs,
        }
    }
    /// Decode from bytes
    pub fn decode(bytes_: &[u8]) -> SpResult<Self> {
        use proto::Circuit;
        use std::io::Cursor;

        let c = Circuit::decode(&mut Cursor::new(bytes_));
        // Unfortunately these conversion errors don't convert to our Result
        let ckt_proto = match c {
            Ok(ckt) => ckt,
            Err(e) => panic!("Circuit Decode Failed"),
        };
        // Create the circuit-parse analysis needs
        let c = Self::from(ckt_proto);
        Ok(c)
    }
    /// Add anything convertible into `CompParse`,
    /// typically the enum-associated structs `Vsrc` et al.
    pub fn add<C: Into<CompParse>>(&mut self, comp: C) {
        self.comps.push(comp.into());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spresult::TestResult;

    #[test]
    fn test_ckt_parse() -> TestResult {
        CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(1e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Gnd),
            ],
        };
        Ok(())
    }
}
