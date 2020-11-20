//!
//! Spice21 Build Module
//!
//! Primarily expands protobuf definitions
//!
extern crate prost_build;

fn main() {
    let mut config = prost_build::Config::new();
    // Add serde traits
    config.type_attribute(".", "#[derive(serde_derive::Serialize, serde_derive::Deserialize)]");
    config.field_attribute("spice21.MosPorts.d", "#[serde(default)]");
    config.field_attribute("spice21.MosPorts.g", "#[serde(default)]");
    config.field_attribute("spice21.MosPorts.s", "#[serde(default)]");
    config.field_attribute("spice21.MosPorts.b", "#[serde(default)]");
    // Nicen up our repeated and enum fields
    config.type_attribute("spice21.Instance.comp", "#[serde(tag = \"type\")]");
    config.field_attribute("spice21.Instance.comp", "#[serde(flatten)]");
    config.type_attribute("spice21.Def.defines", "#[serde(tag = \"type\")]");
    config.field_attribute("spice21.Def.defines", "#[serde(flatten)]");
    // And build!
    config
        .compile_protos(&["src/spice21.proto", "src/mos.proto", "src/bsim4.proto"], &["src/"])
        .unwrap();
}
