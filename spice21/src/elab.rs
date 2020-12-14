use std::collections::HashMap;

use crate::analysis::{Options, VarIndex, Variables};
use crate::circuit;
use crate::circuit::{Comp, NodeRef};
use crate::comps::ComponentSolver;
use crate::SpNum;

///
/// # Hierarchy Elaborator
///
/// Transforms "folded" circuits with sub-module definitions into "flat"
/// sets of `ComponentSolvers`.
/// The primary constructor argument for each `Elaborator` is a `Ckt`,
/// which becomes owned by the `Elaborator`, and is broken into its
/// definitions and instances.
/// Instances are visited depth-first, and primitive instances
/// transformed into `ComponentSolvers`.
/// All circuit definitions (Modules, Models, and the like)
/// are retained for simulation use.
///
pub(crate) struct Elaborator<'a, NumT: SpNum> {
    pub(crate) comps: Vec<ComponentSolver<'a>>,
    pub(crate) vars: Variables<NumT>,
    pub(crate) defs: circuit::Defs,
    pub(crate) path: Vec<String>,
    pub(crate) opts: Options,
}
impl<'a, NumT: SpNum> Elaborator<'a, NumT> {
    /// Get or create a Variable for Node `node`.
    /// Behavior *heavily* depends on the boolean parameter `autonode`.
    /// For `autonode=0`, no variables are created, only this in namespace `ns` are returned.
    /// For `autonode=1`, new nodes are created for any identifier not previously encountered (ala SPICE netlists)
    /// This is (hopefully) a temporary measure.
    fn node_var(&mut self, node: NodeRef, autonode: bool, ns: &mut HashMap<String, Option<VarIndex>>) -> Option<VarIndex> {
        if autonode {
            if let NodeRef::Gnd = node {
                return None;
            }
            self.path.push(node.to_string());
            let pathname = self.path.join(".");
            let var = self.vars.find_or_create(NodeRef::Name(pathname)).clone();
            ns.insert(node.to_string(), var.clone());
            self.path.pop();
            var
        } else {
            match ns.get(&node.to_string()) {
                Some(n) => n.clone(),
                None => panic!("!!!"),
            }
        }
    }
    /// Elaborate a Module or Component Instance
    /// Dispatches based on circuit::Comp variants.
    pub(crate) fn elaborate_instance(&mut self, inst: Comp, ns: &mut HashMap<String, Option<VarIndex>>, autonode: bool) {
        // FIXME: port/signal-name paths
        match inst {
            Comp::R(r) => {
                let circuit::Ri { g, p, n, .. } = r;
                use crate::comps::Resistor;
                let pvar = self.node_var(p, autonode, ns);
                let nvar = self.node_var(n, autonode, ns);
                self.comps.push(Resistor::new(g, pvar.clone(), nvar.clone()).into());
            }
            Comp::C(c) => {
                let circuit::Ci { c, p, n, .. } = c;
                use crate::comps::Capacitor;
                let pvar = self.node_var(p, autonode, ns);
                let nvar = self.node_var(n, autonode, ns);
                self.comps.push(Capacitor::new(c, pvar.clone(), nvar.clone()).into());
            }
            Comp::I(i) => {
                let circuit::Ii { dc, p, n, .. } = i;
                use crate::comps::Isrc;
                let pvar = self.node_var(p, autonode, ns);
                let nvar = self.node_var(n, autonode, ns);
                self.comps.push(Isrc::new(dc, pvar.clone(), nvar.clone()).into());
            }
            Comp::V(x) => self.elaborate_vsrc(x, ns),
            Comp::D(x) => self.elaborate_diode(x, ns),
            Comp::Mos(x) => self.elaborate_mos(x, ns),
            Comp::Module(x) => self.elaborate_module_inst(x, ns),
        }
    }
    pub(crate) fn elaborate_diode(&mut self, d: circuit::DiodeI, ns: &mut HashMap<String, Option<VarIndex>>) {
        use crate::comps::diode;
        // Destruct the key parser-diode attributes
        let circuit::DiodeI { name, params, p, n, .. } = d;
        // Create or retrive the solver node-variables
        // FIXME: note no auto-noding here
        let pvar = ns.get(&p).unwrap().clone();
        let nvar = ns.get(&n).unwrap().clone();
        self.path.push(name.clone());
        // Get our model and params from definitions
        let ddef = match self.defs.diodes.get(&params.clone(), &self.opts) {
            Some(e) => e,
            None => panic!(format!("Parameters not defined: {}", params)),
        };
        let diode::DiodeCacheEntry { model, intp } = ddef;
        // Derive internal params
        let ports = {
            let m = &*(model.read().unwrap());
            diode::DiodePorts::from(self.pathstr(), m, pvar, nvar, &mut self.vars)
        };
        // And create our solver
        let d = diode::Diode {
            ports,
            model,
            intp,
            ..Default::default()
        };
        self.path.pop();
        self.comps.push(d.into());
    }
    pub(crate) fn elaborate_vsrc(&mut self, vi: circuit::Vi, ns: &mut HashMap<String, Option<VarIndex>>) {
        use crate::comps::Vsrc;
        let circuit::Vi { name, p, n, vdc, acm } = vi;
        // Note order of ops here is, as in many cases,
        // effected by the `autonode`-ing
        // Create or retrieve our node-variables
        let pvar = self.node_var(p, self.on_top(), ns);
        let nvar = self.node_var(n, self.on_top(), ns);
        // Create the current variable, named `self.path`
        self.path.push(name);
        let ivar = self.vars.addi(self.pathstr());
        self.path.pop();
        // And create our solver
        self.comps.push(Vsrc::new(vdc, acm, pvar, nvar, ivar).into());
    }
    pub(crate) fn elaborate_mos(&mut self, m: circuit::Mosi, ns: &mut HashMap<String, Option<VarIndex>>) {
        use crate::comps::mos::MosPorts;
        let circuit::Mosi { name, ports, model, params } = m;
        let MosPorts { d, g, s: s_, b } = ports;
        let ports: MosPorts<Option<VarIndex>> = [
            self.node_var(d, self.on_top(), ns),
            self.node_var(g, self.on_top(), ns),
            self.node_var(s_, self.on_top(), ns),
            self.node_var(b, self.on_top(), ns),
        ]
        .into();
        // Add the instance-name to our path
        self.path.push(name);
        // Determine solver-type from defined models
        let c: ComponentSolver = if let Some(m_) = self.defs.bsim4.models.get(&model) {
            use crate::comps::bsim4::bsim4ports::Bsim4Ports;
            use crate::comps::bsim4::Bsim4;
            let (model, inst) = self.defs.bsim4.get(&model, &params).unwrap();
            let ports = Bsim4Ports::from(self.pathstr(), &ports, &model.vals, &inst.intp, &mut self.vars);
            Bsim4::new(ports, model, inst).into()
        } else if let Some(model) = self.defs.mos1.models.get(&model) {
            use crate::comps::{Mos1, Mos1InstanceParams, Mos1Model};
            let params = match self.defs.mos1.insts.get(&params) {
                Some(p) => p.clone(),
                None => panic!(format!("Parameters not defined: {}", params)),
            };
            Mos1::new(Mos1Model::resolve(model.clone()), Mos1InstanceParams::resolve(params), ports.into()).into()
        } else if let Some(mos_type) = self.defs.mos0.get(&model) {
            // Mos0 has no instance params, and only the PMOS/NMOS type as a "model"
            use crate::comps::Mos0;
            Mos0::new(ports.into(), mos_type.clone()).into()
        } else {
            panic!(format!("Model not defined: {}", model));
        };
        self.comps.push(c);
        self.path.pop();
    }
    /// Concatenate our path into a period-separated string
    fn pathstr(&self) -> String {
        self.path.join(".")
    }
    /// Boolean helper function, indicating whether we are currently at top-level
    fn on_top(&self) -> bool {
        self.path.len() == 0
    }
    pub(crate) fn elaborate_module_inst(&mut self, m: circuit::ModuleI, ns: &mut HashMap<String, Option<VarIndex>>) {
        let circuit::ModuleI { name, module, ports, params } = m;
        // FIXME: parameter handling

        let mdef = match self.defs.modules.get(&module) {
            Some(md) => md,
            None => panic!("ModuleDef not found: {}", module),
        };
        // Each Module instance generates a new namespace.
        // Initialize it by grabbing the variables corresponding to each port.
        // By the time we get here, each value in the `m.ports` map
        // must correspond to an existing variable, or elaboration fails.
        // (This is essentially where connections are made.)
        // This variable-map `inst_ns` seeds the module-innards namespace.
        let mut inst_ns: HashMap<String, Option<VarIndex>> = HashMap::new();
        for (k, v) in &ports {
            let var = ns.get(v).unwrap().clone();
            inst_ns.insert(k.clone(), var);
        }
        self.path.push(name);
        if self.path.len() > 1024 {
            panic!("Elaboration Error: Too deep a hierarchy (for now)!");
        }
        self.elaborate_module(mdef.clone(), &mut inst_ns); // FIXME: stop cloning here please!
        self.path.pop();
    }
    /// Elaborate the content of `ModuleDef` `m`.
    pub(crate) fn elaborate_module(&mut self, m: circuit::ModuleDef, ns: &mut HashMap<String, Option<VarIndex>>) {
        let circuit::ModuleDef { name, signals, comps, .. } = m;
        // FIXME: parameter handling

        // Create new Variables for each internal Signal, and add them to the Variable namespace
        for signame in signals.into_iter() {
            self.elaborate_signal(signame, ns);
        }
        // FIXME: check port/ param compatibility
        for inst in comps.into_iter() {
            let comp = if let Some(i) = inst.comp {
                circuit::Comp::from(i)
            } else {
                panic!("Invalid Comp!!!")
            };
            self.elaborate_instance(comp, ns, false);
        }
    }
    /// Create a new Signal at `self.path.signame`, and append it to `ns`.
    pub(crate) fn elaborate_signal(&mut self, signame: String, ns: &mut HashMap<String, Option<VarIndex>>) {
        // FIXME: add checks for name collisions
        self.path.push(signame.clone());
        let pathname = self.path.join(".");
        let var = self.vars.addv(pathname);
        ns.insert(signame, Some(var));
        self.path.pop();
    }
}
/// Elaborate a top-level circuit
/// Returns the generated `Elaborator`, including its flattened `ComponentSolvers`
/// and all definitions carried over from `ckt`.
pub(crate) fn elaborate<'a, T: SpNum>(ckt: circuit::Ckt, opts: Options) -> Elaborator<'a, T> {
    let circuit::Ckt { name, comps, defs, signals } = ckt;
    let mut e = Elaborator {
        comps: Vec::new(),
        vars: Variables::new(),
        defs,
        path: Vec::new(),
        opts,
    };
    // Initialize the top-level namespace with Gnd
    let mut ns: HashMap<String, Option<VarIndex>> = HashMap::new();
    ns.insert("".into(), None);
    // Add Variables for each top-level Signal
    for signame in signals.into_iter() {
        e.elaborate_signal(signame, &mut ns);
    }
    // Visit all of our components
    for inst in comps.into_iter() {
        e.elaborate_instance(inst, &mut ns, true); // FIXME: autonode'ing top-level instances
    }
    e
}
