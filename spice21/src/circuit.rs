//!
//! # Spice21 Circuit
//!
//! Roughly the "middle-end" circuit description.
//!
//! Typically created from either:
//! * (a) Protobuf via `spice21::proto`, or
//! * (b) Directly for testing or in Rust use-cases
//!

use enum_dispatch::enum_dispatch;

use super::comps::mos::MosPorts;
use super::defs::Defs;
use crate::{SpError, SpResult};

use super::proto;
use super::proto::instance::Comp as CompProto;
use super::proto::Circuit as CircuitProto;

// Re-exports
pub use super::proto::Diode as DiodeI;
pub use super::proto::Module as ModuleDef;
pub use super::proto::ModuleInstance as ModuleI;

/// Node Reference
#[derive(Debug, Clone)]
pub enum NodeRef {
    Gnd,
    Num(usize),
    Name(String),
}
use NodeRef::Gnd;
impl NodeRef {
    pub(crate) fn to_string(&self) -> String {
        match self {
            NodeRef::Name(s) => s.clone(),
            NodeRef::Num(s) => s.to_string(),
            NodeRef::Gnd => "".into(),
        }
    }
}
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

/// Utility method: convert anything `Into<String>`-convertible into `String`.
pub(crate) fn s<S: Into<String>>(from: S) -> String {
    from.into()
}
/// Voltage Source Instance
pub struct Vi {
    pub name: String,
    pub vdc: f64,
    pub acm: f64,
    pub p: NodeRef,
    pub n: NodeRef,
}
/// Current Source Instance
pub struct Ii {
    pub name: String,
    pub dc: f64,
    pub acm: f64,
    pub p: NodeRef,
    pub n: NodeRef,
}
/// Resistance (really conductance) Instance
pub struct Ri {
    pub name: String,
    pub g: f64,
    pub p: NodeRef,
    pub n: NodeRef,
}
/// Capacitor Instance
pub struct Ci {
    pub name: String,
    pub c: f64,
    pub p: NodeRef,
    pub n: NodeRef,
}

/// Mos Instance
pub struct Mosi {
    pub(crate) name: String,             // Instance Name
    pub(crate) model: String,            // Model Name
    pub(crate) params: String,           // Instance Param-Set Name
    pub(crate) ports: MosPorts<NodeRef>, // Port Connections
}
///
/// # Component Enum
///
/// Circuits are primarily composed of a vector of these.
/// From and Into methods for each variant are generated by `enum_dispatch` macros.
///
#[enum_dispatch]
pub enum Comp {
    V(Vi),
    I(Ii),
    R(Ri),
    C(Ci),
    D(DiodeI),
    Mos(Mosi),
    Module(ModuleI),
}
// The empty `CompTrait` allows the `enum_dispatch` macros to generate `From` and `Into`
// between the Comp enum and each of its variants.
#[enum_dispatch(Comp)]
trait CompTrait {}

impl Comp {
    /// Replacement for deprecated original `V` enum variant
    pub fn vdc<S: Into<String>>(name: S, vdc: f64, p: NodeRef, n: NodeRef) -> Comp {
        Comp::V(Vi {
            name: name.into(),
            vdc,
            acm: 0.0,
            p,
            n,
        })
    }
    pub fn idc<S: Into<String>>(name: S, dc: f64, p: NodeRef, n: NodeRef) -> Comp {
        Comp::I(Ii {
            name: name.into(),
            dc,
            acm: 0.0,
            p,
            n,
        })
    }
    pub fn r<S: Into<String>>(name: S, g: f64, p: NodeRef, n: NodeRef) -> Comp {
        Comp::R(Ri { name: name.into(), g, p, n })
    }
    pub fn c<S: Into<String>>(name: S, c: f64, p: NodeRef, n: NodeRef) -> Comp {
        Comp::C(Ci { name: name.into(), c, p, n })
    }
    /// Convert from protobuf-generated classes
    pub fn from(c: CompProto) -> Self {
        match c {
            CompProto::I(i) => {
                let x = Ii {
                    name: i.name.into(),
                    p: n(i.p),
                    n: n(i.n),
                    dc: i.dc,
                    acm: 0.0, // FIXME: no value on proto yet
                };
                Comp::I(x)
            }
            CompProto::R(r) => {
                let x = Ri {
                    name: r.name.into(),
                    p: n(r.p),
                    n: n(r.n),
                    g: r.g,
                };
                Comp::R(x)
            }
            CompProto::C(c) => {
                let x = Ci {
                    name: c.name.into(),
                    p: n(c.p),
                    n: n(c.n),
                    c: c.c,
                };
                Comp::C(x)
            }
            CompProto::V(v) => {
                let vs = Vi {
                    name: v.name,
                    p: n(v.p),
                    n: n(v.n),
                    vdc: v.dc,
                    acm: v.acm,
                };
                Comp::V(vs)
            }
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
                Comp::Mos(Mosi {
                    name: m.name,
                    model: m.model.clone(),
                    params: m.params.clone(),
                    ports,
                })
            }
            CompProto::D(x) => Comp::D(x),
            CompProto::X(x) => Comp::Module(x),
        }
    }
}

///
/// # Primary Circuit Structure
///
#[derive(Default)]
pub struct Ckt {
    pub name: String,
    pub signals: Vec<String>,
    pub comps: Vec<Comp>,
    pub defs: Defs,
}
impl Ckt {
    /// Create a new, empty Circuit
    pub fn new() -> Self {
        Self {
            name: String::from(""),
            signals: Vec::new(),
            comps: Vec::new(),
            defs: Defs::default(),
        }
    }
    /// Create a Circuit from a vector of Components
    pub fn from_comps(comps: Vec<Comp>) -> Self {
        Self {
            name: String::from(""),
            signals: Vec::new(),
            comps: comps,
            defs: Defs::default(),
        }
    }
    /// Decode from bytes, via proto definitions
    pub fn decode(bytes_: &[u8]) -> SpResult<Self> { 
        use prost::Message;
        use std::io::Cursor;

        // Decode the protobuf version
        let proto = CircuitProto::decode(&mut Cursor::new(bytes_))?;
        // And convert into a Circuit
        Self::from_proto(proto)
    }
    /// Add anything convertible into `Comp`,
    /// typically the enum-associated structs `Vi` et al.
    pub fn add<C: Into<Comp>>(&mut self, comp: C) {
        self.comps.push(comp.into());
    }
    /// Convert from YAML string  
    pub fn from_yaml(y: &str) -> SpResult<Self> {
        use textwrap::dedent;
        let proto: CircuitProto = serde_yaml::from_str(&dedent(y)).unwrap();
        Self::from_proto(proto)
    }
    /// Convert from TOML string  
    pub fn from_toml(y: &str) -> SpResult<Self> {
        use textwrap::dedent;
        let proto: CircuitProto = toml::from_str(&dedent(y)).unwrap();
        Self::from_proto(proto)
    }
    /// Convert from JSON string  
    pub fn from_json(y: &str) -> SpResult<Self> {
        use textwrap::dedent;
        let proto: CircuitProto = serde_json::from_str(&dedent(y)).unwrap();
        Self::from_proto(proto)
    }
    /// Create from a protobuf-generated circuit
    pub fn from_proto(c: proto::Circuit) -> SpResult<Ckt> {
        let CircuitProto {
            name,
            comps: cs_,
            defs: ds_,
            signals,
        } = c;
        let mut defs = Defs::default();

        use super::proto::def::Defines as DefProto;
        // Step through all definitions
        for def in ds_.into_iter() {
            match def.defines.unwrap() {
                DefProto::Bsim4inst(x) => {
                    defs.bsim4.add_inst(x);
                }
                DefProto::Bsim4model(x) => {
                    use crate::comps::bsim4::Bsim4ModelSpecs;
                    let specs = Bsim4ModelSpecs::from(&x);
                    defs.bsim4.add_model(&x.name, specs);
                }
                DefProto::Mos1model(x) => {
                    use crate::comps::mos::Mos1Model;
                    defs.mos1.add_model(&x.name.clone(), Mos1Model::resolve(&x));
                }
                DefProto::Mos1inst(x) => {
                    use crate::comps::mos::Mos1InstanceParams;
                    defs.mos1.add_inst(&x.name.clone(), Mos1InstanceParams::resolve(&x));
                }
                DefProto::Diodemodel(x) => {
                    use crate::comps::diode::DiodeModel;
                    defs.diodes.add_model(&x.name.clone(), DiodeModel::from(x))
                }
                DefProto::Diodeinst(x) => defs.diodes.add_inst(&x.name.clone(), x),
                DefProto::Module(x) => {
                    defs.modules.add(x);
                }
            }
        }
        // And step through all instances
        let mut comps: Vec<Comp> = vec![];
        for opt in cs_.into_iter() {
            if let Some(c) = opt.comp {
                comps.push(Comp::from(c));
            } else {
                return Err(SpError::new("Invalid Component"));
            }
        }
        Ok(Ckt { comps, defs, name, signals })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert::assert;
    use crate::spresult::TestResult;

    #[test]
    fn test_ckt_parse() -> TestResult {
        let ckt = Ckt::from_comps(vec![
            Comp::idc("i1", 1e-3, NodeRef::Num(0), NodeRef::Gnd),
            Comp::r("r1", 1e-3, NodeRef::Num(0), NodeRef::Gnd),
        ]);
        assert(ckt.comps.len()).eq(2)?;
        Ok(())
    }
    #[test]
    fn test_from_yaml() -> TestResult {
        let ckt = Ckt::from_yaml(
            r#"
            name: tbd
            defs: []
            comps:
              - {type: R, name: r1, p: a, n: "", g: 0.001 }
              - {type: C, name: r2, p: a, n: "", c: 0.001 }
              - {type: M, name: mq, params: noparams, model: nomodel, ports: {d: b, g: a } }
                "#,
        )?;
        assert(ckt.comps.len()).eq(3)?;
        Ok(())
    }
    #[test]
    fn test_from_toml() -> TestResult {
        let ckt = Ckt::from_toml(
            r#"
            name = "tbd"
            defs = []
            comps = [
              {type="R", name= "r1", p= "a", n= "", g= 0.001 },
              {type="C", name= "r2", p= "a", n= "", c= 0.001 }
            ]"#,
        )?;
        assert(ckt.comps.len()).eq(2)?;
        Ok(())
    }
    #[test]
    fn test_from_json() -> TestResult {
        let ckt = Ckt::from_json(
            r#"
            {
                "name": "tbd",
                "comps": [
                    {
                        "type": "I",
                        "name": "ii",
                        "p": "ip",
                        "n": "in",
                        "dc": 1e-12
                    },
                    {
                        "type": "V",
                        "name": "vv",
                        "p": "vp",
                        "n": "vn",
                        "dc": 1e-12,
                        "acm": 0.0
                    },
                    {
                        "type": "C",
                        "name": "cccc",
                        "p": "ac",
                        "n": "bc",
                        "c": 1e-12
                    },
                    {
                        "type": "R",
                        "name": "dtbd",
                        "p": "a",
                        "n": "b",
                        "g": 0.001
                    },
                    {
                        "type": "D",
                        "name": "dtbd",
                        "p": "a",
                        "n": "b",
                        "model": "default",
                        "params": "default"
                    },
                    {
                        "type": "M",
                        "name": "mq",
                        "model": "nomodel",
                        "params": "noparams",
                        "ports": {
                            "d": "b",
                            "g": "a" 
                        }
                    }
                ],
                "defs": []
            }"#,
        )?;
        assert(ckt.comps.len()).eq(6)?;
        Ok(())
    }
}
