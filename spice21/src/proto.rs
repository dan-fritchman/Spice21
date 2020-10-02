//!
//! # Spice21 Protobuf Definitions
//!
use prost::Message;
use serde::{Deserialize, Serialize};

// Include the prost-expanded proto-file content
include!(concat!(env!("OUT_DIR"), "/spice21.proto.rs"));

#[cfg(test)]
mod tests {
    use super::*;
    use crate::*;

    #[test]
    fn test_ckt_proto() -> TestResult {
        use crate::assert::assert;
        use circuit::s;
        use instance::Comp::{C, D, I, M, R, V};

        let c = Circuit {
            name: String::from("tbd"),
            comps: vec![
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
                        g: s("a"),
                        d: s("b"),
                        s: s("b"),
                        b: s("b"),
                    })),
                },
            ],
            defs: vec![],
        };

        use std::io::Cursor;
        // use std::fs::File;
        // use std::io::prelude::*;

        // Prost/ Protobuf Serialization Round-Trip
        let mut buf = Vec::<u8>::new();
        buf.reserve(c.encoded_len());
        c.encode(&mut buf).unwrap();
        // let mut file = File::create("ckt.bin").unwrap();
        // file.write_all(&buf).unwrap();
        let d = Circuit::decode(&mut Cursor::new(buf)).unwrap();
        assert(&c.name).eq(&d.name)?;

        // Serde-JSON Serialization Round-Trip
        use serde::{Deserialize, Serialize};
        let s = serde_json::to_string(&c).unwrap();
        // let mut rfj = File::create("ckt.json").unwrap();
        // rfj.write_all(s.as_bytes()).unwrap();
        let d: Circuit = serde_json::from_str(&s).unwrap();
        assert(&c.name).eq(&d.name)?;

        Ok(())
    }
}
