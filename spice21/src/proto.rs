use super::comps::{Mos1InstanceParams, Mos1Model, MosType};
use prost::Message;

// Include the `items` module, which is generated from items.proto.
pub mod proto {
    include!(concat!(env!("OUT_DIR"), "/spice21.proto.rs"));
    fn some() -> Circuit {
        return Circuit {
            name: String::from("tbd"),
            statements: vec![Statement {
                // aint look great
                statement: Some(statement::Statement::Instance(Instance {
                    name: String::from("???"),
                })),
            }],
        };
    }
}

#[derive(Debug, Copy, Clone)]
pub enum NodeRef {
    Gnd,
    Num(usize),
}

pub struct Vs {
    pub vdc: f64,
    pub acm: f64,
    pub p: NodeRef,
    pub n: NodeRef,
}

pub enum CompParse {
    Vb(Vs),
    V(f64, NodeRef, NodeRef),
    I(f64, NodeRef, NodeRef),
    R(f64, NodeRef, NodeRef),
    C(f64, NodeRef, NodeRef),
    D(f64, f64, NodeRef, NodeRef),
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

pub struct CktParse {
    pub nodes: usize,
    pub comps: Vec<CompParse>,
}

#[cfg(test)]
mod tests {
    use super::NodeRef::{Gnd, Num};
    use super::*;
    use crate::spresult::TestResult;

    /// Create a very basic Circuit
    fn parse_ckt() -> CktParse {
        CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(1e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Gnd),
            ],
        }
    }

    #[test]
    fn test_ckt_parse() -> TestResult {
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(1e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Gnd),
            ],
        };
        Ok(())
    }
}
