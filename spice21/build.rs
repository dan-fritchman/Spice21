//!
//! # Spice21 Build Module
//!
//! Primarily expands protobuf definitions,
//! adding a number of annotations.
//!
use prost_build;

fn main() {
    let mut config = prost_build::Config::new();
    // Add our custom trait
    config.type_attribute("spice21.Circuit", "#[derive(SpProto)]");
    config.type_attribute("spice21.Module", "#[derive(SpProto)]");
    config.type_attribute("spice21.Def", "#[derive(SpProto)]");
    config.type_attribute("spice21.Defs", "#[derive(SpProto)]");
    config.type_attribute("spice21.Sim", "#[derive(SpProto)]");
    config.type_attribute("spice21.SimOptions", "#[derive(SpProto)]");
    config.type_attribute("spice21.SimResult", "#[derive(SpProto)]");
    config.type_attribute("spice21.Op", "#[derive(SpProto)]");
    config.type_attribute("spice21.OpResult", "#[derive(SpProto)]");
    config.type_attribute("spice21.Tran", "#[derive(SpProto)]");
    config.type_attribute("spice21.TranResult", "#[derive(SpProto)]");
    config.type_attribute("spice21.Ac", "#[derive(SpProto)]");
    config.type_attribute("spice21.AcResult", "#[derive(SpProto)]");

    // Add serde traits
    config.type_attribute(".", "#[derive(serde_derive::Serialize, serde_derive::Deserialize)]");

    // Add field attributes
    // Ideally, we'll figure out how to add more than one at a time
    config.field_attribute("spice21.Circuit.name", "#[serde(default)]");
    config.field_attribute("spice21.Circuit.signals", "#[serde(default)]");
    config.field_attribute("spice21.Circuit.params", "#[serde(default)]");
    config.field_attribute("spice21.Circuit.defs", "#[serde(default)]");
    config.field_attribute("spice21.Circuit.comps", "#[serde(default)]");

    config.field_attribute("spice21.MosPorts.g", "#[serde(default)]");
    config.field_attribute("spice21.MosPorts.s", "#[serde(default)]");
    config.field_attribute("spice21.MosPorts.b", "#[serde(default)]");
    config.field_attribute("spice21.MosPorts.d", "#[serde(default)]");

    // Nicen up our repeated and enum fields
    config.type_attribute("spice21.Instance.comp", "#[serde(tag = \"type\")]");
    config.field_attribute("spice21.Instance.comp", "#[serde(flatten)]");
    config.type_attribute("spice21.Def.defines", "#[serde(tag = \"type\")]");
    config.field_attribute("spice21.Def.defines", "#[serde(flatten)]");
    // And build!
    config.compile_protos(&["protos/spice21.proto"], &["protos/"]).unwrap();
}
