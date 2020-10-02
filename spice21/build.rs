//! 
//! Spice21 Build Module
//! 
//! Primarily expand protobuf definitions
//! 
extern crate prost_build;

fn main() {
    println!("Spice21 build.rs");
    let mut config = prost_build::Config::new();
    // Add serde traits 
    config.type_attribute(".", "#[derive(serde_derive::Serialize, serde_derive::Deserialize)]");
    config.compile_protos(&["src/spice21.proto"], &["src/"]).unwrap();
}
