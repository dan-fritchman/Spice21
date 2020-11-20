/// "Integration" Tests
#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use crate::analysis::*;
    use crate::assert::*;
    use crate::circuit::NodeRef::{Gnd, Num};
    use crate::circuit::*;
    use crate::comps::*;
    use crate::spresult::*;

    /// Create a very basic Circuit
    #[test]
    fn test_ckt_parse() -> TestResult {
        Ckt::from_comps(vec![
            Comp::idc("i1", 1e-3, NodeRef::Name(s("0")), NodeRef::Gnd),
            Comp::r("r1", 1e-3, NodeRef::Name(s("0")), NodeRef::Gnd),
        ]);
        Ok(())
    }
    /// R-Only DCOP
    #[test]
    fn test_dcop1() -> TestResult {
        let ckt = Ckt::from_yaml(
            r#"
            name: tbd
            defs: []
            comps:
              - {type: R, name: r1, p: a, n: "", g: 0.001 }"#,
        )?;
        let soln = dcop(ckt)?;
        assert_eq!(soln.values, vec![0.0]);
        Ok(())
    }
    /// I-R DCOP
    #[test]
    fn test_dcop2() -> TestResult {
        let ckt = Ckt::from_comps(vec![
            Comp::idc("i1", 1e-3, NodeRef::Num(0), NodeRef::Gnd),
            Comp::r("r1", 1e-3, NodeRef::Num(0), NodeRef::Gnd),
        ]);

        let soln = dcop(ckt)?;
        assert_eq!(soln.values, vec![1.0,]);
        Ok(())
    }
    /// I - R - R divider
    #[test]
    fn test_dcop3() -> TestResult {
        let ckt = Ckt::from_comps(vec![
            Comp::r("r1", 1e-3, NodeRef::Num(0), NodeRef::Gnd),
            Comp::r("r1", 1e-3, NodeRef::Num(1), NodeRef::Num(0)),
            Comp::idc("i1", 1e-3, NodeRef::Num(1), NodeRef::Gnd),
        ]);

        let soln = dcop(ckt)?;
        assert!((soln[0] - 1.0).abs() < 1e-4);
        assert!((soln[1] - 2.0).abs() < 1e-4);
        Ok(())
    }

    /// V - R - R divider
    #[test]
    fn test_dcop4() -> TestResult {
        use Comp::R;
        use NodeRef::Gnd;
        let ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::r("r1", 2e-3, n("vdd"), n("div")),
            Comp::r("r2", 2e-3, n("div"), Gnd),
        ]);
        let soln = dcop(ckt)?;
        assert(soln.values).eq(vec![-1e-3, 1.0, 0.5])?;
        Ok(())
    }

    /// Diode DcOp Tests
    /// Voltage & Current-Biased
    #[test]
    fn test_dcop5() -> TestResult {
        // I - R - Diode
        use crate::circuit::{Di, Vi};
        use NodeRef::Gnd;

        // Voltage-biased Diode
        let v = 0.70;
        let mut ckt = Ckt::new();
        ckt.add(Di::new("dd", n("p"), Gnd));
        ckt.add(Vi {
            name: s("vin"),
            p: n("p"),
            n: Gnd,
            vdc: v,
            acm: 0.0,
        });

        let soln = dcop(ckt)?;
        let i = soln.map.get(&s("vin")).ok_or("Failed to measure current")?.abs();
        // Some broad bounds checks
        assert(i).gt(1e-3)?;
        assert(i).lt(100e-3)?;

        // Current-biased Diode, with the measured current
        let mut ckt = Ckt::new();
        ckt.add(Di::new("dd", n("p"), Gnd));
        ckt.add(Comp::idc("i1", i, n("p"), Gnd));

        // Check the voltage matches our initial v-bias
        let soln = dcop(ckt)?;
        assert(soln[0]).isclose(v, 1e-3)?;
        assert(soln[0] - v).abs().lt(1e-3)?; // (same thing really)
        Ok(())
    }

    #[test]
    fn test_dcop6() -> TestResult {
        // NMOS Char
        use NodeRef::{Gnd, Num};
        let mut ckt = Ckt::from_comps(vec![
            Comp::Mos(Mosi {
                name: s("m"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Num(1),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::vdc("v1", 1.0, NodeRef::Num(0), NodeRef::Gnd),
            Comp::vdc("v1", 1.0, NodeRef::Num(1), NodeRef::Gnd),
        ]);
        add_mos0_defaults(&mut ckt);

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], 1.0);
        assert_eq!(soln[1], 1.0);
        assert_eq!(soln[2], 0.0);
        assert!((soln[3] + 14.1e-3).abs() < 1e-4);
        Ok(())
    }

    #[test]
    fn test_dcop7() -> TestResult {
        // PMOS Char
        use NodeRef::{Gnd, Num};
        let mut ckt = Ckt::from_comps(vec![
            Comp::Mos(Mosi {
                name: s("m"),
                model: "pmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Num(1),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::vdc("v1", -1.0, NodeRef::Num(0), NodeRef::Gnd),
            Comp::vdc("v1", -1.0, NodeRef::Num(1), NodeRef::Gnd),
        ]);
        add_mos0_defaults(&mut ckt);

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], -1.0);
        assert_eq!(soln[1], -1.0);
        assert_eq!(soln[2], 0.0);
        assert!((soln[3] - 14.1e-3).abs() < 1e-4);
        Ok(())
    }

    #[test]
    fn test_dcop8() -> TestResult {
        // Diode NMOS
        use NodeRef::{Gnd, Num};
        let mut ckt = Ckt::from_comps(vec![
            Comp::idc("i1", 5e-3, Num(0), Gnd),
            Comp::Mos(Mosi {
                name: s("m"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::r("r1", 1e-12, Num(0), Gnd), // "gmin"
        ]);
        add_mos0_defaults(&mut ckt);

        let soln = dcop(ckt)?;
        assert!((soln[0] - 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_diode_nmos_tran() -> TestResult {
        // Diode NMOS Tran
        use NodeRef::{Gnd, Num};
        let mut ckt = Ckt::from_comps(vec![
            Comp::idc("i1", 5e-3, Num(0), Gnd),
            Comp::Mos(Mosi {
                name: s("m"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::r("r1", 1e-12, Num(0), Gnd), // "gmin"
        ]);
        add_mos0_defaults(&mut ckt);
        let opts = TranOptions {
            tstep: 1e-12,
            tstop: 100e-12,
            ..Default::default()
        };
        let soln = tran(ckt, opts)?;
        for point in soln.data.iter() {
            assert!((point[0] - 0.697).abs() < 1e-3);
        }
        Ok(())
    }

    #[test]
    fn test_dcop8b() -> TestResult {
        // Diode NMOS, S/D Swapped
        use NodeRef::{Gnd, Num};
        let mut ckt = Ckt::from_comps(vec![
            Comp::idc("i1", 5e-3, Num(0), Gnd),
            Comp::Mos(Mosi {
                name: s("m"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Gnd,
                    s: Num(0),
                    b: Gnd,
                },
            }),
            Comp::r("r1", 1e-12, Num(0), Gnd), // "gmin"
        ]);
        add_mos0_defaults(&mut ckt);

        let soln = dcop(ckt)?;
        assert!((soln[0] - 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_diode_pmos_dcop() -> TestResult {
        // Diode PMOS
        use NodeRef::{Gnd, Num};
        let mut ckt = Ckt::from_comps(vec![
            Comp::idc("i1", -5e-3, Num(0), Gnd),
            Comp::Mos(Mosi {
                name: s("m"),
                model: "pmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::r("r1", 1e-12, Num(0), Gnd), // "gmin"
        ]);
        add_mos0_defaults(&mut ckt);

        let soln = dcop(ckt)?;
        assert!((soln[0] + 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_diode_pmos_tran() -> TestResult {
        // Diode PMOS Tran
        use NodeRef::{Gnd, Num};
        let mut ckt = Ckt::from_comps(vec![
            Comp::idc("i1", -5e-3, Num(0), Gnd),
            Comp::Mos(Mosi {
                name: s("m"),
                model: "pmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::r("r1", 1e-12, Num(0), Gnd), // "gmin"
        ]);
        add_mos0_defaults(&mut ckt);

        let opts = TranOptions {
            tstep: 1e-12,
            tstop: 100e-12,
            ..Default::default()
        };
        let soln = tran(ckt, opts)?;
        for point in soln.data.iter() {
            assert!((point[0] + 0.697).abs() < 1e-3);
        }
        Ok(())
    }

    #[test]
    fn test_dcop8d() -> TestResult {
        // Diode PMOS, S/D Swapped
        let mut ckt = Ckt::from_comps(vec![
            Comp::idc("i1", -5e-3, Num(0), Gnd),
            Comp::Mos(Mosi {
                name: s("m"),
                model: "pmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Gnd,
                    s: Num(0),
                    b: Gnd,
                },
            }),
            Comp::r("r1", 1e-12, Num(0), Gnd), // "gmin"
        ]);
        add_mos0_defaults(&mut ckt);

        let soln = dcop(ckt)?;
        assert!((soln[0] + 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_dcop9() -> TestResult {
        // NMOS-R, "Grounded"
        use NodeRef::{Gnd, Num};
        let mut ckt = Ckt::from_comps(vec![
            Comp::r("r1", 1e-3, Num(0), Gnd),
            Comp::Mos(Mosi {
                name: s("m"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
        ]);
        add_mos0_defaults(&mut ckt);

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], 0.0);
        Ok(())
    }

    #[test]
    fn test_dcop9b() -> TestResult {
        // NMOS-R, "Grounded", S/D Swapped
        use NodeRef::{Gnd, Num};
        let mut ckt = Ckt::from_comps(vec![
            Comp::r("r1", 1e-3, Num(0), Gnd),
            Comp::Mos(Mosi {
                name: s("m"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Gnd,
                    s: Num(0),
                    b: Gnd,
                },
            }),
        ]);
        add_mos0_defaults(&mut ckt);

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], 0.0);
        Ok(())
    }

    /// PMOS-R, "Grounded"
    #[test]
    fn test_dcop9c() -> TestResult {
        use NodeRef::{Gnd, Num};
        let mut ckt = Ckt::from_comps(vec![
            Comp::r("r1", 1e-3, Num(0), Gnd),
            Comp::Mos(Mosi {
                name: s("m"),
                model: "pmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
        ]);
        add_mos0_defaults(&mut ckt);

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], 0.0);
        Ok(())
    }

    /// PMOS-R, "Grounded", S/D Swapped
    #[test]
    fn test_dcop9d() -> TestResult {
        use NodeRef::{Gnd, Num};
        let mut ckt = Ckt::from_comps(vec![
            Comp::r("r1", 1e-3, Num(0), Gnd),
            Comp::Mos(Mosi {
                name: s("m"),
                model: "pmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Gnd,
                    s: Num(0),
                    b: Gnd,
                },
            }),
        ]);
        add_mos0_defaults(&mut ckt);

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], 0.0);
        Ok(())
    }

    /// NMOS-R Inverter
    #[test]
    fn test_dcop10() -> TestResult {
        let mut ckt = Ckt::from_comps(vec![
            Comp::r("r1", 1e-3, n("vdd"), n("d")),
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::Mos(Mosi {
                name: s("m"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: n("vdd"),
                    d: n("d"),
                    s: Gnd,
                    b: Gnd,
                },
            }),
        ]);
        add_mos0_defaults(&mut ckt);

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], 1.0);
        assert!(soln[1] < 50e-3);
        assert!((soln[2] + 1e-3).abs() < 0.1e-3);
        Ok(())
    }

    /// PMOS-R Inverter
    #[test]
    fn test_dcop10b() -> TestResult {
        let mut ckt = Ckt::from_comps(vec![
            Comp::r("r1", 1e-3, n("g"), n("d")),
            Comp::vdc("v1", -1.0, n("g"), Gnd),
            Comp::Mos(Mosi {
                name: s("m"),
                model: "pmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: n("g"),
                    d: n("d"),
                    s: Gnd,
                    b: Gnd,
                },
            }),
        ]);
        add_mos0_defaults(&mut ckt);

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], -1.0);
        assert!(soln[1].abs() < 50e-3);
        assert!((soln[2] - 1e-3).abs() < 0.1e-3);
        Ok(())
    }
    /// Mos0 CMOS Inverter DC-Op, Vin=Vdd
    #[test]
    fn test_dcop11() -> TestResult {
        use MosType::{NMOS, PMOS};
        use NodeRef::Gnd;

        let mut ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::r("r1", 1e-9, n("d"), Gnd), // "gmin"
            Comp::Mos(Mosi {
                name: s("p"),
                model: "pmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: n("vdd"),
                    d: n("d"),
                    s: n("vdd"),
                    b: n("vdd"),
                },
            }),
            Comp::Mos(Mosi {
                name: s("n"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: n("vdd"),
                    d: n("d"),
                    s: Gnd,
                    b: Gnd,
                },
            }),
        ]);
        add_mos0_defaults(&mut ckt);

        let soln = dcop(ckt)?;
        assert(soln.values).eq(vec![0.0, 1.0, 0.0])?;
        Ok(())
    }
    /// Mos0 CMOS Inverter DC-Op, Vin=Vss
    #[test]
    fn test_dcop11b() -> TestResult {
        use MosType::{NMOS, PMOS};
        use NodeRef::Gnd;

        let mut ckt = Ckt::from_comps(vec![
            Comp::Mos(Mosi {
                name: s("p"),
                model: "pmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Gnd,
                    d: n("d"),
                    s: n("vdd"),
                    b: n("vdd"),
                },
            }),
            Comp::Mos(Mosi {
                name: s("n"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Gnd,
                    d: n("d"),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::r("r1", 1e-9, n("d"), n("vdd")), // "gmin"
        ]);
        add_mos0_defaults(&mut ckt);

        let soln = dcop(ckt)?;
        assert(soln.values).eq(vec![1.0, 1.0, 0.0])?;
        Ok(())
    }
    /// DCOP, Several Series CMOS Inverters
    #[test]
    fn test_dcop12() -> TestResult {
        use NodeRef::{Gnd, Num};
        let mut ckt = Ckt::from_comps(vec![
            Comp::r("r1", 1e-9, Num(0), Gnd),
            Comp::r("r1", 1e-9, Num(1), Gnd),
            Comp::r("r1", 1e-9, Num(2), Gnd),
            Comp::r("r1", 1e-9, Num(3), Gnd),
            Comp::r("r1", 1e-9, Num(4), Gnd),
            Comp::Mos(Mosi {
                name: s("p1"),
                model: "pmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Num(1),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos(Mosi {
                name: s("n1"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Num(1),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::Mos(Mosi {
                name: s("p2"),
                model: "pmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos(Mosi {
                name: s("n2"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::Mos(Mosi {
                name: s("p3"),
                model: "pmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos(Mosi {
                name: s("n3"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::Mos(Mosi {
                name: s("p4"),
                model: "pmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(3),
                    d: Num(4),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos(Mosi {
                name: s("n4"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(3),
                    d: Num(4),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::vdc("v1", 1.0, Num(0), Gnd),
        ]);
        add_mos0_defaults(&mut ckt);

        let soln = dcop(ckt)?;
        assert(soln[0]).eq(1.0)?;
        assert!(soln[1].abs() < 1e-3);
        assert!((soln[2] - 1.0).abs() < 1e-3);
        assert!(soln[3].abs() < 1e-3);
        assert!((soln[4] - 1.0).abs() < 1e-3);
        assert!(soln[5].abs() < 1e-6);
        Ok(())
    }

    /// RC Low-Pass Filter DcOp
    #[test]
    fn test_dcop13() -> TestResult {
        use NodeRef::{Gnd, Num};

        let ckt = Ckt::from_comps(vec![
            Comp::r("r1", 1e-3, Num(1), Num(0)),
            Comp::c("c1", 1e-9, Num(1), Gnd),
            Comp::vdc("v1", 1.0, Num(0), Gnd),
        ]);

        let soln = dcop(ckt)?;
        assert_eq!(soln.values, vec![1.0, 1.0, 0.0]);
        Ok(())
    }

    // RC High-Pass Filter DcOp
    #[test]
    fn test_dcop13b() -> TestResult {
        let ckt = Ckt::from_comps(vec![
            Comp::c("c1", 1e-9, n("i"), n("o")),
            Comp::r("r1", 1e-3, n("o"), Gnd),
            Comp::vdc("v1", 1.0, n("i"), Gnd),
        ]);

        let soln = dcop(ckt)?;
        assert_eq!(soln.values, vec![1.0, 0.0, 0.0]);
        Ok(())
    }

    /// RC Low-Pass Filter Tran
    #[test]
    fn test_tran1() -> TestResult {
        // Circuit
        let ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, n("inp"), Gnd),
            Comp::r("r1", 1e-3, n("inp"), n("out")),
            Comp::c("c1", 1e-9, n("out"), Gnd),
        ]);
        // Simulate
        let opts = TranOptions {
            tstep: 10e-9,
            tstop: 10e-6,
            ic: vec![(n("out"), 0.0)],
        };
        let soln = tran(ckt, opts)?;
        // Checks
        let inp = soln.get("inp")?;
        assert(inp).is().constant(1.0)?;
        let out = soln.get("out")?;
        assert(out[0]).abs().lt(1e-3)?;
        assert(out[out.len() - 1]).isclose(1.0, 1e-3)?;
        assert(out).is().increasing()?;
        Ok(())
    }

    /// I-C Integrator with Initial Condition
    #[test]
    #[ignore]
    fn test_tran2() -> TestResult {
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![Comp::idc("i1", 1e-3, Num(0), Gnd), Comp::c("c1", 4e-12, Num(0), Gnd)]);

        let opts = TranOptions {
            tstep: 1e-18,
            tstop: 1e-15,
            ..Default::default()
        };
        let mut tran = Tran::new(ckt, opts);
        tran.ic(Num(0), 0.0);
        let soln = tran.solve()?;

        assert(soln[0][0]).eq(5e-3)?;
        assert(soln[0][1]).eq(0.0)?;
        assert(soln[0][2]).eq(5e-3)?;
        for k in 1..soln.len() {
            assert((soln[k][0] - soln[k - 1][0] - 5e-3).abs()).lt(1e-6)?;
            assert(soln[k][1]).eq(0.0)?;
            assert(soln[k][2]).lt(1e-6)?;
        }
        Ok(())
    }
    /// I-C Integrator with Initial Condition
    #[test]
    #[ignore] // FIXME: failing values to be debugged
    fn test_tran2b() -> TestResult {
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![Comp::idc("i1", 1e-6, Num(0), Gnd), Comp::c("c1", 100e-9, Num(0), Gnd)]);

        let opts = TranOptions {
            tstep: 1e-21,
            tstop: 1e-18,
            ..Default::default()
        };
        let mut tran = Tran::new(ckt, opts);
        tran.ic(Num(0), 0.0);
        let soln = tran.solve()?;

        assert(soln[0][0]).eq(5e-3)?;
        assert(soln[0][1]).eq(0.0)?;
        assert(soln[0][2]).eq(5e-3)?;
        for k in 1..soln.len() {
            assert((soln[k][0] - soln[k - 1][0] - 5e-3).abs()).lt(1e-6)?;
            assert(soln[k][1]).eq(0.0)?;
            assert(soln[k][2]).lt(1e-6)?;
        }
        Ok(())
    }

    /// Mos0 Ring Oscillator (a very fast one)
    #[test]
    fn test_mos0_cmos_ro_tran() -> TestResult {
        let opts = TranOptions {
            tstep: 1e-50, // Told you its fast
            tstop: 1e-47, // See?
            ic: vec![(Num(1), 0.0)],
        };
        let c = 2.5e-50; // Yeah especially this part
        let mut ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, Num(0), Gnd),
            Comp::r("r1", 1e-3, Num(0), Gnd),
            Comp::Mos(Mosi {
                name: s("p1"),
                model: "pmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(3),
                    d: Num(1),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos(Mosi {
                name: s("n1"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(3),
                    d: Num(1),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::r("r1", 1e-9, Num(1), Gnd),
            Comp::c("c1", c, Num(1), Gnd),
            Comp::Mos(Mosi {
                name: s("p2"),
                model: "pmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos(Mosi {
                name: s("n2"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::r("r1", 1e-9, Num(2), Gnd),
            Comp::c("c1", c, Num(2), Gnd),
            Comp::Mos(Mosi {
                name: s("p3"),
                model: "pmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos(Mosi {
                name: s("n3"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::r("r1", 1e-9, Num(3), Gnd),
            Comp::c("c1", c, Num(3), Gnd),
        ]);
        add_mos0_defaults(&mut ckt);
        // Simulate
        let soln = tran(ckt, opts)?;
        // Checks
        //to_file(&soln, "test_mos0_cmos_ro_tran.json"); // Writes new golden data
        let golden = load_golden("test_mos0_cmos_ro_tran.json");
        assert(&soln.map).isclose(golden, 1e-6)?;
        Ok(())
    }
    #[test]
    fn test_mos1_op() -> TestResult {
        let mut ckt = Ckt::from_comps(vec![
            Comp::Mos(Mosi {
                name: s("m"),
                model: "default".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::vdc("v1", 1.0, Num(0), Gnd),
        ]);
        add_mos1_defaults(&mut ckt);
        // Simulate
        let soln = dcop(ckt)?;
        // Checks
        assert(soln[0]).eq(1.0)?;
        assert(soln[1]).lt(0.0)?;
        assert(soln[1]).gt(-1e-3)?;
        Ok(())
    }

    #[test]
    fn test_mos1_tran() -> TestResult {
        let mut ckt = Ckt::from_comps(vec![
            Comp::Mos(Mosi {
                name: s("m"),
                model: "default".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::vdc("v1", 1.0, Num(0), Gnd),
        ]);
        add_mos1_defaults(&mut ckt);
        // Simulate
        let opts = TranOptions {
            tstep: 1e-9,
            tstop: 100e-9,
            ..Default::default()
        };
        let soln = tran(ckt, opts)?;
        // Checks
        for k in 1..soln.len() {
            assert(soln[k][0]).eq(1.0)?;
            assert(soln[k][1]).lt(0.0)?;
            assert(soln[k][1]).gt(-1e-3)?;
        }
        Ok(())
    }
    /// Mos1 Inverter DCOP
    #[test]
    fn test_mos1_inv_dcop() -> TestResult {
        let mut ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, Num(0), Gnd),
            Comp::Mos(Mosi {
                name: s("p"),
                model: "pmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Gnd,
                    d: Num(1),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos(Mosi {
                name: s("n"),
                model: "nmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Gnd,
                    d: Num(1),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::r("r1", 1e-4, Num(1), Gnd),
        ]);
        // Define our models & params
        add_mos1_defaults(&mut ckt);

        dcop(ckt)?;
        // FIXME: checks on solution
        Ok(())
    }
    /// Create a three-stage CMOS RO,
    /// with device-models named `pmos` and `nmos`,
    /// and instance-parameter-sets named `default.
    fn cmos_ro3() -> Ckt {
        Ckt::from_yaml(
            r#"
            name: ro
            defs: []
            comps:
              - {type: V, name: v1, p: vdd, n: "", dc: 1.0, acm: 0.0 }
              - {type: M, name: p1, params: default, model: pmos, ports: {d: "2", g: "1", s: vdd, b: vdd} }
              - {type: M, name: n1, params: default, model: nmos, ports: {d: "2", g: "1", s: "", b: ""} }
              - {type: M, name: p2, params: default, model: pmos, ports: {d: "3", g: "2", s: vdd, b: vdd} }
              - {type: M, name: n2, params: default, model: nmos, ports: {d: "3", g: "2", s: "", b: ""} }
              - {type: M, name: p3, params: default, model: pmos, ports: {d: "1", g: "3", s: vdd, b: vdd} }
              - {type: M, name: n3, params: default, model: nmos, ports: {d: "1", g: "3", s: "", b: ""} }
                "#,
        )
        .unwrap()
    }
    /// Mos1 CMOS Ring Oscillator Dc Op
    #[test]
    fn test_mos1_cmos_ro_dcop() -> TestResult {
        let mut ckt = cmos_ro3(); // Shared Circuit
        add_mos1_defaults(&mut ckt); // Add Mos1 Models & Params
                                     // Simulate
        let soln = dcop(ckt)?;
        // Checks
        assert(soln.get("vdd")?).eq(1.0)?;
        for k in ["1", "2", "3"].into_iter() {
            assert(soln.get(*k)?).gt(0.45)?;
            assert(soln.get(*k)?).lt(0.55)?;
        }
        Ok(())
    }
    /// Mos1 CMOS Ring Oscillator Tran
    #[test]
    fn test_mos1_cmos_ro_tran() -> TestResult {
        let mut ckt = cmos_ro3(); // Shared Circuit
        add_mos1_defaults(&mut ckt); // Add Mos1 Models & Params
                                     // Simulate
        let opts = TranOptions {
            tstep: 1e-11,
            tstop: 1e-8,
            ic: vec![(Num(1), 0.0)],
        };
        let soln = tran(ckt, opts)?;
        // to_file(&soln, "test_mos1_cmos_ro_tran.json"); // Writes new golden data
        // Checks
        let golden = load_golden("test_mos1_cmos_ro_tran.json");
        assert(&soln.map).isclose(golden, 1e-6)?;
        Ok(())
    }
    /// Bsim4 CMOS Ring Oscillator Tran
    #[test]
    fn test_bsim4_cmos_ro_tran() -> TestResult {
        let mut ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::Mos(Mosi {
                name: s("p1"),
                model: "pmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(3),
                    d: Num(1),
                    s: n("vdd"),
                    b: n("vdd"),
                },
            }),
            Comp::Mos(Mosi {
                name: s("n1"),
                model: "nmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(3),
                    d: Num(1),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::Mos(Mosi {
                name: s("p2"),
                model: "pmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: n("vdd"),
                    b: n("vdd"),
                },
            }),
            Comp::Mos(Mosi {
                name: s("n2"),
                model: "nmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::Mos(Mosi {
                name: s("p3"),
                model: "pmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: n("vdd"),
                    b: n("vdd"),
                },
            }),
            Comp::Mos(Mosi {
                name: s("n3"),
                model: "nmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Gnd,
                    b: Gnd,
                },
            }),
        ]);
        add_bsim4_defaults(&mut ckt); // Add Bsim4 Models & Params
                                      // Simulate
        let opts = TranOptions {
            tstep: 1e-10,
            tstop: 3e-7,
            ic: vec![(Num(1), 0.0)],
        };
        let soln = tran(ckt, opts)?;
        // to_file(&soln, "test_bsim4_cmos_ro_tran.json"); // Writes new golden data
        // Checks
        let golden = load_golden("test_bsim4_cmos_ro_tran.json");
        assert(&soln.map).isclose(golden, 1e-6)?;
        Ok(())
    }

    // Mos1 NMOS-R Oscillator Tran
    #[test]
    fn test_mos1_nmos_ro_tran() -> TestResult {
        let gl = 1e-6;
        let mut ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::Mos(Mosi {
                name: s("n1"),
                model: "nmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(3),
                    d: Num(1),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::Mos(Mosi {
                name: s("n2"),
                model: "nmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::Mos(Mosi {
                name: s("n3"),
                model: "nmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::r("r1", gl, Num(1), n("vdd")),
            Comp::r("r1", gl, Num(2), n("vdd")),
            Comp::r("r1", gl, Num(3), n("vdd")),
        ]);
        // Simulate
        let opts = TranOptions {
            tstep: 1e-11,
            tstop: 1e-8,
            ic: vec![(Num(1), 0.0)],
        };
        add_mos1_defaults(&mut ckt);
        // Simulate
        let opts = TranOptions {
            tstep: 1e-9,
            tstop: 10e-6,
            ..Default::default()
        };
        let soln = tran(ckt, opts)?;
        //to_file(&soln, "test_mos1_nmos_ro_tran.json"); // Writes new golden data
        // Checks
        let golden = load_golden("test_mos1_nmos_ro_tran.json");
        assert(&soln.map).isclose(golden, 1e-6)?;
        Ok(())
    }

    // Mos1 PMOS-R Oscillator Tran
    #[test]
    fn test_mos1_pmos_ro_tran() -> TestResult {
        let cl = 1e-16;
        let gl = 1e-6;
        let mut ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::Mos(Mosi {
                name: s("p1"),
                model: "pmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(3),
                    d: Num(1),
                    s: n("vdd"),
                    b: n("vdd"),
                },
            }),
            Comp::Mos(Mosi {
                name: s("p2"),
                model: "pmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: n("vdd"),
                    b: n("vdd"),
                },
            }),
            Comp::Mos(Mosi {
                name: s("p3"),
                model: "pmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: n("vdd"),
                    b: n("vdd"),
                },
            }),
            Comp::r("r1", gl, Num(1), Gnd),
            Comp::r("r1", gl, Num(2), Gnd),
            Comp::r("r1", gl, Num(3), Gnd),
            Comp::c("c1", cl, Num(1), Gnd),
            Comp::c("c1", cl, Num(2), Gnd),
            Comp::c("c1", cl, Num(3), Gnd),
        ]);
        add_mos1_defaults(&mut ckt);
        // Simulate
        let opts = TranOptions {
            tstep: 1e-11,
            tstop: 1e-8,
            ic: vec![(Num(1), 0.0)],
        };
        let soln = tran(ckt, opts)?;
        //to_file(&soln, "test_mos1_pmos_ro_tran.json"); // Writes new golden data
        // Checks
        let golden = load_golden("test_mos1_pmos_ro_tran.json");
        assert(&soln.map).isclose(golden, 1e-6)?;
        Ok(())
    }

    // Bsim4 PMOS-R Oscillator Tran
    #[test]
    fn test_bsim4_pmos_ro_tran() -> TestResult {
        let cl = 1e-16;
        let gl = 1e-6;
        let mut ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::Mos(Mosi {
                name: s("p1"),
                model: "pmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(3),
                    d: Num(1),
                    s: n("vdd"),
                    b: n("vdd"),
                },
            }),
            Comp::Mos(Mosi {
                name: s("p2"),
                model: "pmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: n("vdd"),
                    b: n("vdd"),
                },
            }),
            Comp::Mos(Mosi {
                name: s("p3"),
                model: "pmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: n("vdd"),
                    b: n("vdd"),
                },
            }),
            Comp::r("r1", gl, Num(1), Gnd),
            Comp::r("r1", gl, Num(2), Gnd),
            Comp::r("r1", gl, Num(3), Gnd),
            Comp::c("c1", cl, Num(1), Gnd),
            Comp::c("c1", cl, Num(2), Gnd),
            Comp::c("c1", cl, Num(3), Gnd),
        ]);
        add_bsim4_defaults(&mut ckt);
        // Simulate
        let opts = TranOptions {
            tstep: 1e-11,
            tstop: 1e-8,
            ic: vec![(Num(1), 0.0)],
        };
        let soln = tran(ckt, opts)?;
        // to_file(&soln, "test_bsim4_pmos_ro_tran.json"); // Writes new golden data
        // Checks
        let golden = load_golden("test_bsim4_pmos_ro_tran.json");
        assert(&soln.map).isclose(golden, 1e-6)?;
        Ok(())
    }

    /// Mos1 PMOS-R Amp, Tran Initial Condition Decay
    #[test]
    fn test_mos1_pmos_rload_tran() -> TestResult {
        let gl = 1e-6;
        let mut ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::Mos(Mosi {
                name: s("p1"),
                model: "pmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: n("inp"),
                    d: n("out"),
                    s: n("vdd"),
                    b: n("vdd"),
                },
            }),
            Comp::r("r1", gl, n("inp"), n("vdd")),
            Comp::r("r1", gl, n("out"), Gnd),
        ]);
        add_mos1_defaults(&mut ckt);
        // Simulate
        let opts = TranOptions {
            tstep: 1e-11,
            tstop: 1e-8,
            ic: vec![(n("inp"), 0.0)],
        };
        let soln = tran(ckt, opts)?;
        // Checks
        let inp = soln.get("inp")?;
        assert(inp[0]).isclose(0.0, 1e-6)?;
        assert(inp[inp.len() - 1]).isclose(1.0, 5e-3)?;
        assert(inp).is().nondecreasing()?;
        let out = soln.get("out")?;
        assert(out[0] - 1.0).abs().le(0.1)?;
        assert(out[out.len() - 1]).abs().le(1e-3)?;
        // assert(out).is().decreasing()?; // FIXME: time step 0-1 increases slightly
        Ok(())
    }
    // Mos1 PMOS-R Tran
    #[test]
    fn test_mos1_pmos_rg_tran() -> TestResult {
        let gl = 1e-6;
        let mut ckt = Ckt::from_comps(vec![
            Comp::Mos(Mosi {
                name: s("p1"),
                model: "pmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: n("g"),
                    d: Gnd,
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::r("r1", gl, n("g"), Gnd),
        ]);
        add_mos1_defaults(&mut ckt);
        // Simulate
        let opts = TranOptions {
            tstep: 1e-11,
            tstop: 1e-8,
            ic: vec![(n("g"), -1.0)],
        };
        let soln = tran(ckt, opts)?;
        // Checks
        let g = soln.get("g")?;
        assert(g[0]).isclose(-1.0, 1e-3)?;
        assert(g).last().isclose(0.0, 1e-3)?;
        Ok(())
    }
    // Bsim4 PMOS-R Tran
    #[test]
    fn test_bsim4_pmos_rg_tran() -> TestResult {
        let gl = 1e-5;
        let mut ckt = Ckt::from_comps(vec![
            Comp::Mos(Mosi {
                name: s("p1"),
                model: "pmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: n("g"),
                    d: Gnd,
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::r("r1", gl, n("g"), Gnd),
            Comp::c("c1", 1e-18, n("g"), Gnd),
        ]);
        add_bsim4_defaults(&mut ckt);
        // Simulate
        let opts = TranOptions {
            tstep: 1e-10,
            tstop: 1e-7,
            ic: vec![(n("g"), -1.0)],
        };
        let soln = tran(ckt, opts)?;
        // to_file(&soln, "rg.json")?;
        // Checks
        let g = soln.get("g")?;
        assert(g[0]).isclose(-1.0, 1e-3)?;
        assert(g).last().isclose(0.0, 1e-3)?;
        Ok(())
    }
    /// Mos1 NMOS-R Tran
    #[test]
    fn test_mos1_nmos_rg_tran() -> TestResult {
        let gl = 1e-6;
        let mut ckt = Ckt::from_comps(vec![
            Comp::Mos(Mosi {
                name: s("p1"),
                model: "nmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: n("g"),
                    d: Gnd,
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::r("r1", gl, n("g"), Gnd),
        ]);
        add_mos1_defaults(&mut ckt);
        // Simulate
        let opts = TranOptions {
            tstep: 1e-11,
            tstop: 1e-8,
            ic: vec![(n("g"), 1.0)],
        };
        let soln = tran(ckt, opts)?;
        // Checks
        let g = soln.get("g")?;
        assert(g[0]).isclose(1.0, 1e-3)?;
        assert(g).last().isclose(0.0, 1e-3)?;
        Ok(())
    }
    /// Bsim4 NMOS-R Tran
    #[test]
    fn test_bsim4_nmos_rg_tran() -> TestResult {
        let gl = 1e-5;
        let mut ckt = Ckt::from_comps(vec![
            Comp::Mos(Mosi {
                name: s("p1"),
                model: "nmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: n("g"),
                    d: Gnd,
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::r("r1", gl, n("g"), Gnd),
            // Comp::c("c1", 1e-15, n("g"), Gnd),
        ]);
        add_bsim4_defaults(&mut ckt);
        // Simulate
        let opts = TranOptions {
            tstep: 1e-10,
            tstop: 1e-7,
            ic: vec![(n("g"), 1.0)],
        };
        let soln = tran(ckt, opts)?;
        // to_file(&soln, "rg.json")?;
        // Checks
        let g = soln.get("g")?;
        assert(g[0]).isclose(1.0, 1e-3)?;
        assert(g).last().isclose(0.0, 1e-3)?;
        Ok(())
    }
    #[test]
    fn test_ac1() -> TestResult {
        let ckt = Ckt::from_comps(vec![Comp::r("r1", 1.0, Num(0), Gnd)]);
        ac(ckt, AcOptions::default())?;
        // FIXME: checks on solution
        Ok(())
    }

    #[test]
    fn test_ac2() -> TestResult {
        use crate::circuit::Vi;
        let ckt = Ckt::from_comps(vec![
            Comp::r("r1", 1e-3, Num(0), Num(1)),
            Comp::c("c1", 1e-9, Num(1), Gnd),
            Comp::V(Vi {
                name: s("vi"),
                vdc: 1.0,
                acm: 1.0,
                p: Num(0),
                n: Gnd,
            }),
        ]);
        ac(ckt, AcOptions::default())?;
        // FIXME: checks on solution
        Ok(())
    }

    #[test]
    #[ignore] // FIXME: aint no Mos0 AC!
    fn test_ac3() -> TestResult {
        let mut ckt = Ckt::from_comps(vec![
            Comp::r("r1", 1e-3, Num(0), Num(1)),
            Comp::c("c1", 1e-9, Num(1), Gnd),
            Comp::vdc("v1", 1.0, Num(0), Gnd),
            Comp::Mos(Mosi {
                name: s("m"),
                model: "nmos".into(),
                params: "".into(),
                ports: MosPorts {
                    g: Num(1),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
        ]);
        add_mos0_defaults(&mut ckt);
        ac(ckt, AcOptions::default())?;
        // FIXME: checks on solution
        Ok(())
    }

    /// NMOS Common-Source Amp
    #[test]
    fn test_ac4() -> TestResult {
        let mut ckt = Ckt::from_comps(vec![
            Comp::c("c1", 1e-9, n("d"), Gnd),
            Comp::Mos(Mosi {
                name: s("m"),
                model: "default".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: n("g"),
                    d: n("d"),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::V(Vi {
                name: s("vg"),
                vdc: 0.7,
                acm: 1.0,
                p: n("g"),
                n: Gnd,
            }),
        ]);
        // Define our models & params
        use crate::proto::{Mos1InstParams, Mos1Model};
        let nmos = Mos1Model {
            mos_type: MosType::NMOS as i32,
            ..Mos1Model::default()
        };
        ckt.models.mos1.models.insert("default".into(), nmos);
        let params = Mos1InstParams::default();
        ckt.models.mos1.insts.insert("default".into(), params);
        ac(ckt, AcOptions::default())?;
        // FIXME: checks on solution
        Ok(())
    }

    /// Diode-Connected NMOS AC
    #[test]
    fn test_ac5() -> TestResult {
        use crate::circuit::Vi;

        let mut ckt = Ckt::from_comps(vec![
            Comp::V(Vi {
                name: s("vd"),
                vdc: 0.5,
                acm: 1.0,
                p: Num(0),
                n: Gnd,
            }),
            Comp::Mos(Mosi {
                name: s("m"),
                model: "default".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
        ]);

        // Define our models & params
        use crate::proto::{Mos1InstParams, Mos1Model};
        let nmos = Mos1Model {
            mos_type: MosType::NMOS as i32,
            ..Mos1Model::default()
        };
        ckt.models.mos1.models.insert("default".into(), nmos);
        let params = Mos1InstParams::default();
        ckt.models.mos1.insts.insert("default".into(), params);
        ac(ckt, AcOptions::default())?;
        // FIXME: checks on solution
        Ok(())
    }

    // Bsim4 NMOS-R Oscillator Tran
    #[test]
    fn test_bsim4_nmos_ro_tran() -> TestResult {
        let gl = 1e-6;
        let cl = 1e-18;
        let mut ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::Mos(Mosi {
                name: s("n1"),
                model: "nmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(3),
                    d: Num(1),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::Mos(Mosi {
                name: s("n2"),
                model: "nmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::Mos(Mosi {
                name: s("n3"),
                model: "nmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::r("r1", gl, Num(1), n("vdd")),
            Comp::r("r1", gl, Num(2), n("vdd")),
            Comp::r("r1", gl, Num(3), n("vdd")),
            Comp::c("c1", cl, Num(1), n("vdd")),
            Comp::c("c1", cl, Num(2), n("vdd")),
            Comp::c("c1", cl, Num(3), n("vdd")),
        ]);
        // Simulate
        let opts = TranOptions {
            tstep: 1e-9,
            tstop: 1e-6,
            ic: vec![(Num(1), 0.0)],
        };
        add_bsim4_defaults(&mut ckt);
        let soln = tran(ckt, opts)?;
        // to_file(&soln, "test_bsim4_nmos_ro_tran.json"); // Writes new golden data
        // Checks
        let golden = load_golden("test_bsim4_nmos_ro_tran.json");
        assert(&soln.map).isclose(golden, 1e-6)?;
        Ok(())
    }

    /// Test-helper to write results to JSON file
    #[allow(dead_code)]
    fn to_file(soln: &TranResult, fname: &str) -> SpResult<()> {
        use std::fs::File;
        use std::io::prelude::*;

        #[allow(unused_imports)] // Need these traits in scope
        use serde::ser::{SerializeSeq, Serializer};

        let mut rfj = File::create(fname).unwrap();
        let s = serde_json::to_string(&soln.map).unwrap();
        rfj.write_all(s.as_bytes()).unwrap();
        Ok(())
    }
    /// Read golden results
    fn load_golden(fname: &str) -> HashMap<String, Vec<f64>> {
        use std::fs::File;
        use std::io::BufReader;
        use std::path::PathBuf;

        let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        d.push("resources/");
        d.push(fname);
        let file = File::open(d).unwrap();
        let reader = BufReader::new(file);
        let golden: HashMap<String, Vec<f64>> = serde_json::from_reader(reader).unwrap();
        golden
    }
    /// Helper. Modifies `ckt` adding Mos0 defaults
    fn add_mos0_defaults(ckt: &mut Ckt) {
        ckt.models.mos0.insert("default".into(), MosType::NMOS);
        ckt.models.mos0.insert("nmos".into(), MosType::NMOS);
        ckt.models.mos0.insert("pmos".into(), MosType::PMOS);
    }
    /// Helper. Modifies `ckt` adding Mos1 default instance-params, plus default NMOS and PMOS
    fn add_mos1_defaults(ckt: &mut Ckt) {
        use crate::proto::{Mos1InstParams, Mos1Model};
        let nmos = Mos1Model {
            mos_type: MosType::NMOS as i32,
            ..Mos1Model::default()
        };
        let default = nmos.clone();
        ckt.models.mos1.models.insert("default".into(), default);
        ckt.models.mos1.models.insert("nmos".into(), nmos);
        let pmos = Mos1Model {
            mos_type: MosType::PMOS as i32,
            ..Mos1Model::default()
        };
        ckt.models.mos1.models.insert("pmos".into(), pmos);
        let params = Mos1InstParams::default();
        ckt.models.mos1.insts.insert("default".into(), params);
    }
    /// Helper. Modifies `ckt` adding Bsim4 default instance-params, plus default NMOS and PMOS
    fn add_bsim4_defaults(ckt: &mut Ckt) {
        use crate::comps::bsim4::{Bsim4InstSpecs, Bsim4ModelSpecs};
        let nmos = Bsim4ModelSpecs::new(MosType::NMOS);
        let default = nmos.clone();
        ckt.models.bsim4.models.insert("default".into(), default);
        ckt.models.bsim4.models.insert("nmos".into(), nmos);
        let pmos = Bsim4ModelSpecs::new(MosType::PMOS);
        ckt.models.bsim4.models.insert("pmos".into(), pmos);
        let params = Bsim4InstSpecs::default();
        ckt.models.bsim4.insts.insert("default".into(), params);
    }
}
