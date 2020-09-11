use super::comps::{Mos1InstanceParams, Mos1Model, MosType};
// use prost::Message;

pub mod proto {
    // Include the prost-expanded proto-file content
    include!(concat!(env!("OUT_DIR"), "/spice21.proto.rs"));

    // fn some() -> Circuit {
    //     return Circuit {
    //         name: String::from("tbd"),
    //         statements: vec![Statement {
    //             // aint look great
    //             statement: Some(statement::Statement::Instance(Instance {
    //                 name: String::from("???"),
    //             })),
    //         }],
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
pub fn n<S: Into<String>>(name: S) -> NodeRef {
    NodeRef::Name(name.into())
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
    // D(f64, f64, NodeRef, NodeRef),
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
