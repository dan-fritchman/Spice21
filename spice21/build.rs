extern crate prost_build;

fn main() {
    println!("Spice21 build.rs");
    prost_build::compile_protos(&["src/spice21.proto"], &["src/"]).unwrap();
}
