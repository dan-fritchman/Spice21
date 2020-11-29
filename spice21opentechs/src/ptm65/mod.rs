use serde_yaml;
use spice21::circuit::Defs;
use spice21::comps::bsim4::Bsim4ModelSpecs;

pub fn defs() -> Defs {
    let mut defs = Defs::new();
    let pmos: Bsim4ModelSpecs = serde_yaml::from_str(include_str!("pmos.yaml")).unwrap();
    let nmos: Bsim4ModelSpecs = serde_yaml::from_str(include_str!("nmos.yaml")).unwrap();

    // FIXME: actual public addition method
    // defs.bsim4.models.insert("nmos".into(), nmos);
    // defs.bsim4.models.insert("pmos".into(), pmos);
    
    defs
}

#[cfg(test)]
mod tests;
