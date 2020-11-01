//!
//! # Spice21 Protobuf Definitions
//!

#[allow(unused_imports)] // These are used by the macro-expanded code
use serde::{Deserialize, Serialize};
#[allow(unused_imports)]
use prost::Message;

// Include the prost-expanded proto-file content
include!(concat!(env!("OUT_DIR"), "/spice21.rs"));


#[cfg(test)]
mod tests {
    use super::*;
    use crate::*;

    #[test]
    fn test_ckt_proto() -> TestResult {
        use crate::assert::assert;
        use circuit::s;
        use instance::Comp::{C, D, I, M, R, V};
        use def::Defines;

        let c = Circuit {
            name: String::from("tbd"),
            comps: vec![
                Instance {
                    comp: Some(C(Capacitor {
                        name: s("cccc"),
                        p: s("ac"),
                        n: s("bc"),
                        c: 1e-12,
                    })),
                },
                Instance {
                    comp: Some(R(Resistor {
                        name: s("dtbd"),
                        p: s("a"),
                        n: s("b"),
                        g: 1e-3,
                    })),
                },
                Instance {
                    comp: Some(D(Diode {
                        name: s("dtbd"),
                        p: s("a"),
                        n: s("b"),
                        area: 1.0,
                        temp: 300.0,
                    })),
                },
                Instance {
                    comp: Some(M(Mos {
                        name: s("mq"),
                        ports: Some(MosPorts {
                            g: s("a"),
                            d: s("b"),
                            s: s("c"),
                            b: s("d"),
                        }),
                        model: s("nomodel"),
                        params: s("noparams"),
                    })),
                },
            ],
            defs: vec![
                Def {
                    defines: Some(Defines::Mos1model(
                        Mos1Model::default()
                    ))
                }
            ],
        };

        use std::io::Cursor;
        use std::fs::File;
        use std::io::prelude::*;

        // Prost/ Protobuf Serialization Round-Trip
        let mut buf = Vec::<u8>::new();
        buf.reserve(c.encoded_len());
        c.encode(&mut buf).unwrap();
        // let mut file = File::create("ckt.bin").unwrap();
        // file.write_all(&buf).unwrap();
        let d = Circuit::decode(&mut Cursor::new(buf)).unwrap();
        assert(&c.name).eq(&d.name)?;

        // Serde-JSON Serialization Round-Trip
        let s = serde_json::to_string(&c).unwrap();
        let mut rfj = File::create("ckt.json").unwrap();
        rfj.write_all(s.as_bytes()).unwrap();
        let d: Circuit = serde_json::from_str(&s).unwrap();
        assert(&c.name).eq(&d.name)?;
        // Serde-YAML Round-Trip
        // let s = serde_yaml::to_string(&c).unwrap();
        // let mut rfj = File::create("ckt.yaml").unwrap();
        // rfj.write_all(s.as_bytes()).unwrap();
        
        Ok(())
    }
}