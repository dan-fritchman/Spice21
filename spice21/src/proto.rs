//!
//! # Spice21 Protobuf Definitions
//!

#[allow(unused_imports)]
use prost::Message;
#[allow(unused_imports)] // These are used by the macro-expanded code
use serde::{Deserialize, Serialize};

// Error Conversion
use crate::SpError;
impl From<prost::DecodeError> for SpError {
    fn from(_e: prost::DecodeError) -> Self {
        SpError::new("Error Decoding Circuit Binary")
    }
}

// Include the prost-expanded proto-file content
include!(concat!(env!("OUT_DIR"), "/spice21.rs"));

#[cfg(test)]
mod tests {
    use super::def::Defines;
    use super::instance::Comp;
    use super::*;
    use crate::assert::assert;
    use crate::circuit::s;
    use crate::*;
    use std::collections::HashMap;

    /// Proto-defined struct creation, and round-tripping via alternate serializations
    #[test]
    fn test_ckt_proto() -> TestResult {
        let mut p: HashMap<String, f64> = HashMap::new();
        p.insert("p1".into(), 1.0);
        let mut conns: HashMap<String, String> = HashMap::new();
        conns.insert("a".into(), "a1".into());
        conns.insert("b".into(), "b1".into());
        conns.insert("c".into(), "c1".into());
        conns.insert("d".into(), "d1".into());

        let c = Circuit {
            name: s("tbd"),
            signals: vec!["a".into(), "b".into()],
            comps: vec![
                Instance {
                    comp: Some(Comp::I(Isrc {
                        name: s("ii"),
                        p: s("ip"),
                        n: s("in"),
                        dc: 1e-12,
                    })),
                },
                Instance {
                    comp: Some(Comp::V(Vsrc {
                        name: s("vv"),
                        p: s("vp"),
                        n: s("vn"),
                        dc: 1e-12,
                        acm: 0.0,
                    })),
                },
                Instance {
                    comp: Some(Comp::C(Capacitor {
                        name: s("cccc"),
                        p: s("ac"),
                        n: s("bc"),
                        c: 1e-12,
                    })),
                },
                Instance {
                    comp: Some(Comp::R(Resistor {
                        name: s("dtbd"),
                        p: s("a"),
                        n: s("b"),
                        g: 1e-3,
                    })),
                },
                Instance {
                    comp: Some(Comp::D(Diode {
                        name: s("dtbd"),
                        p: s("a"),
                        n: s("b"),
                        area: 1.0,
                        temp: 300.0,
                    })),
                },
                Instance {
                    comp: Some(Comp::M(Mos {
                        name: s("mq"),
                        model: s("nomodel"),
                        params: s("noparams"),
                        ports: Some(MosPorts {
                            g: s("a"),
                            d: s("b"),
                            s: s("c"),
                            b: s("d"),
                        }),
                    })),
                },
                Instance {
                    comp: Some(Comp::X(ModuleInstance {
                        name: s("xxx"),
                        module: s("good_luck"),
                        params: p.clone(),
                        ports: conns.clone(),
                    })),
                },
            ],
            defs: vec![
                Def {
                    defines: Some(Defines::Mos1model(Mos1Model::default())),
                },
                Def {
                    defines: Some(Defines::Module(Module {
                        name: "good_luck".into(),
                        ports: vec!["a".into(), "b".into(), "c".into(), "d".into()],
                        params: p.clone(),
                        signals: vec!["i1".into(), "i2".into()],
                        comps: vec![
                            Instance {
                                comp: Some(Comp::I(Isrc {
                                    name: s("ii"),
                                    p: s("ip"),
                                    n: s("in"),
                                    dc: 1e-12,
                                })),
                            },
                            Instance {
                                comp: Some(Comp::V(Vsrc {
                                    name: s("vv"),
                                    p: s("vp"),
                                    n: s("vn"),
                                    dc: 1e-12,
                                    acm: 0.0,
                                })),
                            },
                            Instance {
                                comp: Some(Comp::C(Capacitor {
                                    name: s("cccc"),
                                    p: s("ac"),
                                    n: s("bc"),
                                    c: 1e-12,
                                })),
                            },
                            Instance {
                                comp: Some(Comp::R(Resistor {
                                    name: s("dtbd"),
                                    p: s("a"),
                                    n: s("b"),
                                    g: 1e-3,
                                })),
                            },
                            Instance {
                                comp: Some(Comp::D(Diode {
                                    name: s("dtbd"),
                                    p: s("a"),
                                    n: s("b"),
                                    area: 1.0,
                                    temp: 300.0,
                                })),
                            },
                            Instance {
                                comp: Some(Comp::M(Mos {
                                    name: s("mq"),
                                    model: s("nomodel"),
                                    params: s("noparams"),
                                    ports: Some(MosPorts {
                                        g: s("a"),
                                        d: s("b"),
                                        s: s("c"),
                                        b: s("d"),
                                    }),
                                })),
                            },
                        ],
                    })),
                },
            ],
        };

        use std::fs::{create_dir_all, File};
        use std::io::prelude::*;
        use std::io::Cursor;

        // Make sure we have a `scratch/` directory
        create_dir_all("scratch").unwrap();

        // Prost/ Protobuf Serialization Round-Trip
        let mut buf = Vec::<u8>::new();
        buf.reserve(c.encoded_len());
        c.encode(&mut buf).unwrap();
        let mut file = File::create("scratch/ckt.bin").unwrap();
        file.write_all(&buf).unwrap();
        let d = Circuit::decode(&mut Cursor::new(buf)).unwrap();
        assert(&c.name).eq(&d.name)?;

        // Serde-JSON Round-Trip
        let s = serde_json::to_string(&c).unwrap();
        let mut rfj = File::create("scratch/ckt.json").unwrap();
        rfj.write_all(s.as_bytes()).unwrap();
        let d: Circuit = serde_json::from_str(&s).unwrap();
        assert(&c.name).eq(&d.name)?;

        // Serde-YAML Round-Trip
        let s = serde_yaml::to_string(&c).unwrap();
        let mut rfj = File::create("scratch/ckt.yaml").unwrap();
        rfj.write_all(s.as_bytes()).unwrap();
        let d: Circuit = serde_yaml::from_str(&s).unwrap();
        assert(&c.name).eq(&d.name)?;

        // Serde-TOML Round-Trip
        let s = toml::to_string(&c).unwrap();
        let mut rfj = File::create("scratch/ckt.toml").unwrap();
        rfj.write_all(s.as_bytes()).unwrap();
        let d: Circuit = toml::from_str(&s).unwrap();

        Ok(())
    }
}
