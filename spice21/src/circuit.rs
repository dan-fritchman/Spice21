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
use std::collections::HashMap;

use super::comps::bsim4::Bsim4Cache;
use super::comps::mos::{MosPorts, MosType};
use super::comps::{DiodeInstParams, DiodeModel};
use super::proto::def::Defines as DefProto;
use super::proto::instance::Comp as CompProto;
use super::proto::Circuit as CircuitProto;

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
impl From<Vs> for Comp {
    fn from(x: Vs) -> Self {
        Comp::V(x)
    }
}
/// Current Source
pub struct Ii {
    pub name: String,
    pub dc: f64,
    pub acm: f64,
    pub p: NodeRef,
    pub n: NodeRef,
}
impl From<Ii> for Comp {
    fn from(x: Ii) -> Self {
        Comp::I(x)
    }
}
/// Resistance (really conductance)
pub struct Ri {
    pub name: String,
    pub g: f64,
    pub p: NodeRef,
    pub n: NodeRef,
}
impl From<Ri> for Comp {
    fn from(x: Ri) -> Self {
        Comp::R(x)
    }
}
/// Capacitor
pub struct Ci {
    pub name: String,
    pub c: f64,
    pub p: NodeRef,
    pub n: NodeRef,
}
impl From<Ci> for Comp {
    fn from(x: Ci) -> Self {
        Comp::C(x)
    }
}


/// Diode Instance
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

/// Mos Instance
pub struct Mosi {
    pub(crate) name: String,             // Instance Name
    pub(crate) model: String,            // Model Name
    pub(crate) params: String,           // Instance Param-Set Name
    pub(crate) ports: MosPorts<NodeRef>, // Port Connections
}
/// "Level Zero" MOS Instance
pub struct Mos0i {
    pub(crate) name: String,
    pub(crate) mos_type: MosType,
    pub(crate) ports: MosPorts<NodeRef>,
}
/// Component Enum.
/// Circuits are mostly a list of these variants.
// FIXME: #[enum_dispatch]
pub enum Comp {
    V(Vs),
    I(Ii),
    R(Ri),
    C(Ci),
    D(Ds),
    Mos(Mosi),
    Mos0(Mos0i),
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
        Comp::R(Ri {
            name: name.into(),
            g,
            p,
            n,
        })
    }
    pub fn c<S: Into<String>>(name: S, c: f64, p: NodeRef, n: NodeRef) -> Comp {
        Comp::C(Ci {
            name: name.into(),
            c,
            p,
            n,
        })
    }
    /// Convert from protobuf-generated classes
    // FIXME: remove the unused
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
                    g: r.g
                };
                Comp::R(x)
            }
            CompProto::C(c) => {
                let x = Ci {
                    name: c.name.into(),
                    p: n(c.p),
                    n: n(c.n),
                    c: c.c 
                };
                Comp::C(x)
            }
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
        }
    }
}

impl From<Mosi> for Comp {
    fn from(m: Mosi) -> Self {
        Comp::Mos(m)
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
    pub(crate) mos0: HashMap<String, MosType>,
    pub(crate) mos1: Mos1Defs,
    pub(crate) bsim4: Bsim4Cache,
}
impl ModelCache {
    pub(crate) fn new() -> Self {
        Self {
            mos0: HashMap::new(),
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
                    use crate::comps::bsim4::Bsim4ModelSpecs;
                    let specs = Bsim4ModelSpecs::from(&x);
                    models.bsim4.add_model(&x.name, specs);
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
                cs.push(Comp::from(c));
            } else {
                return Err(SpError::new("Invalid Component"));
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
    /// Convert from YAML string  
    pub fn from_yaml(y: &str) -> SpResult<Self> {
        use textwrap::dedent;
        let proto: CircuitProto = serde_yaml::from_str(&dedent(y)).unwrap();
        Self::from(proto)
    }
    /// Convert from TOML string  
    pub fn from_toml(y: &str) -> SpResult<Self> {
        // use serde_yaml;
        use textwrap::dedent;

        let proto: CircuitProto = toml::from_str(&dedent(y)).unwrap();
        Self::from(proto)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert::assert;
    use crate::spresult::TestResult;

    #[test]
    fn test_ckt_parse() -> TestResult {
        let ckt = Ckt {
            comps: vec![
                Comp::idc("i1", 1e-3, NodeRef::Num(0), NodeRef::Gnd),
                Comp::r("r1", 1e-3, NodeRef::Num(0), NodeRef::Gnd),
            ],
            models: ModelCache::new(),
        };
        assert(ckt.comps.len()).eq(2)?;
        Ok(())
    }
    #[test]
    fn test_from_yaml() -> TestResult {
        use serde::{Deserialize, Serialize};
        #[derive(Serialize, Deserialize)]
        struct E {
            f: String,
            h: String,
        }
        #[derive(Serialize, Deserialize)]
        struct A {
            a: String,
            c: String,
            e: E,
        }
        let y: E = serde_yaml::from_str("{f: good, h: great}").unwrap();
        let y: A = serde_yaml::from_str("{ a: b, c: d, e: {f: good, h: great } }").unwrap();
        let ckt = Ckt::from_yaml(
            r#"
            name: tbd
            defs: []
            comps:
              - {type: R, name: r1, p: a, n: "", g: 0.001 }
              - {type: C, name: r2, p: a, n: "", c: 0.001 }
              - {type: M, name: mq, params: noparams, model: nomodel, ports: {d: b, g: a, s: c, b: d} }
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
}
