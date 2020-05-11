// struct Circuit {
//     comps: Vec<Box<dyn Component>>,
//     nodes: Vec<Node>,
//     vars: usize, // FIXME: more elaborate non-node variable handling
//     node0: Node,
// }

// impl Circuit {
//     fn new() -> Circuit {
//         Circuit {
//             comps: vec![],
//             nodes: vec![],
//             vars: 0,
//             node0: Node {
//                 rf: NodeRef::Gnd,
//                 solve: false,
//             },
//         }
//     }
//     fn add_node(&mut self) -> NodeRef {
//         let rf = NodeRef::Num(self.nodes.len());
//         let node = Node {
//             rf: rf.clone(),
//             solve: true,
//         };
//         self.nodes.push(node);
//         return rf;
//     }
//     fn add_var(&mut self) -> usize {
//         let nv = self.vars;
//         self.vars += 1;
//         return nv;
//     }
//     fn add_comp<C: Component + 'static>(&mut self, comp: C) {
//         comp.on_add(&self);
//         self.comps.push(Box::new(comp));
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn test_add_node() -> TestResult {
    //     let mut c = Circuit::new();
    //     let rf = c.add_node();
    //     match rf {
    //         NodeRef::Gnd => return Err("gnd?"),
    //         NodeRef::Num(n) => assert_eq!(n, 0),
    //     };
    //     assert_eq!(c.nodes.len(), 1);
    //     assert_eq!(c.comps.len(), 0);
    //     Ok(())
    // }

    // #[test]
    // fn test_add_res() -> TestResult {
    //     let mut c = Circuit::new();
    //     let rf = c.add_node();
    //     let r = Resistor {
    //         g: 1e-3,
    //         p: rf,
    //         n: NodeRef::Gnd,
    //         pp: None,
    //         nn: None,
    //         np: None,
    //         pn: None,
    //     };
    //     c.add_comp(r);
    //     Ok(())
    // }
}
