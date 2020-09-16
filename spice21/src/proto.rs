//!
//! # Spice21 Protobuf Definitions
//!
use prost::Message;

// Include the prost-expanded proto-file content
include!(concat!(env!("OUT_DIR"), "/spice21.proto.rs"));

#[cfg(test)]
mod tests {
    use super::*;
    use crate::*;

    #[test]
    fn test1() -> TestResult {
        use circuit::s;
        use instance::Comp::{Mos, C, D, I, R, V};
        let c = Circuit {
            name: String::from("tbd"),
            comps: vec![Instance {
                comp: Some(D(Diode {
                    name: s("dtbd"),
                    p: s("a"),
                    n: s("b"),
                    area: 1.0,
                    temp: 300.0,
                })),
            }],
            defs: vec![],
        };
        Ok(())
    }
}
