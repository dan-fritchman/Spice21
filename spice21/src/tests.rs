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
            Comp::I(1e-3, NodeRef::Name(s("0")), NodeRef::Gnd),
            Comp::R(1e-3, NodeRef::Name(s("0")), NodeRef::Gnd),
        ]);
        Ok(())
    }
    /// R-Only DCOP
    #[test]
    fn test_dcop1() -> TestResult {
        let ckt = Ckt::from_comps(vec![Comp::R(1e-3, NodeRef::Num(0), NodeRef::Gnd)]);

        let soln = dcop(ckt)?;
        assert_eq!(soln.values, vec![0.0]);
        Ok(())
    }
    /// I-R DCOP
    #[test]
    fn test_dcop2() -> TestResult {
        let ckt = Ckt::from_comps(vec![
            Comp::I(1e-3, NodeRef::Num(0), NodeRef::Gnd),
            Comp::R(1e-3, NodeRef::Num(0), NodeRef::Gnd),
        ]);

        let soln = dcop(ckt)?;
        assert_eq!(soln.values, vec![1.0,]);
        Ok(())
    }
    /// I - R - R divider
    #[test]
    fn test_dcop3() -> TestResult {
        let ckt = Ckt::from_comps(vec![
            Comp::R(1e-3, NodeRef::Num(0), NodeRef::Gnd),
            Comp::R(1e-3, NodeRef::Num(1), NodeRef::Num(0)),
            Comp::I(1e-3, NodeRef::Num(1), NodeRef::Gnd),
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
            R(2e-3, n("vdd"), n("div")),
            R(2e-3, n("div"), Gnd),
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
        use crate::circuit::{Ds, Vs};
        use Comp::I;
        use NodeRef::Gnd;

        // Voltage-biased Diode
        let v = 0.70;
        let mut ckt = Ckt::new();
        ckt.add(Ds::new("dd", n("p"), Gnd));
        ckt.add(Vs {
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
        ckt.add(Ds::new("dd", n("p"), Gnd));
        ckt.add(I(i, n("p"), Gnd));

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
        let ckt = Ckt::from_comps(vec![
            Comp::Mos0(Mos0i {
                name: s("m"),
                mos_type: MosType::NMOS,
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
        let ckt = Ckt::from_comps(vec![
            Comp::Mos0(Mos0i {
                name: s("m"),
                mos_type: MosType::PMOS,
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
        let ckt = Ckt::from_comps(vec![
            Comp::I(5e-3, Num(0), Gnd),
            Comp::Mos0(Mos0i {
                name: s("m"),
                mos_type: MosType::NMOS,
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::R(1e-12, Num(0), Gnd), // "gmin"
        ]);

        let soln = dcop(ckt)?;
        assert!((soln[0] - 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_diode_nmos_tran() -> TestResult {
        // Diode NMOS Tran
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![
            Comp::I(5e-3, Num(0), Gnd),
            Comp::Mos0(Mos0i {
                name: s("m"),
                mos_type: MosType::NMOS,
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::R(1e-12, Num(0), Gnd), // "gmin"
        ]);
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
        let ckt = Ckt::from_comps(vec![
            Comp::I(5e-3, Num(0), Gnd),
            Comp::Mos0(Mos0i {
                name: s("m"),
                mos_type: MosType::NMOS,
                ports: MosPorts {
                    g: Num(0),
                    d: Gnd,
                    s: Num(0),
                    b: Gnd,
                },
            }),
            Comp::R(1e-12, Num(0), Gnd), // "gmin"
        ]);

        let soln = dcop(ckt)?;
        assert!((soln[0] - 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_diode_pmos_dcop() -> TestResult {
        // Diode PMOS
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![
            Comp::I(-5e-3, Num(0), Gnd),
            Comp::Mos0(Mos0i {
                name: s("m"),
                mos_type: MosType::PMOS,
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::R(1e-12, Num(0), Gnd), // "gmin"
        ]);

        let soln = dcop(ckt)?;
        assert!((soln[0] + 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_diode_pmos_tran() -> TestResult {
        // Diode PMOS Tran
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![
            Comp::I(-5e-3, Num(0), Gnd),
            Comp::Mos0(Mos0i {
                name: s("m"),
                mos_type: MosType::PMOS,
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::R(1e-12, Num(0), Gnd), // "gmin"
        ]);

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
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![
            Comp::I(-5e-3, Num(0), Gnd),
            Comp::Mos0(Mos0i {
                name: s("m"),
                mos_type: MosType::PMOS,
                ports: MosPorts {
                    g: Num(0),
                    d: Gnd,
                    s: Num(0),
                    b: Gnd,
                },
            }),
            Comp::R(1e-12, Num(0), Gnd), // "gmin"
        ]);

        let soln = dcop(ckt)?;
        assert!((soln[0] + 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_dcop9() -> TestResult {
        // NMOS-R, "Grounded"
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![
            Comp::R(1e-3, Num(0), Gnd),
            Comp::Mos0(Mos0i {
                name: s("m"),
                mos_type: MosType::NMOS,
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
        ]);

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], 0.0);
        Ok(())
    }

    #[test]
    fn test_dcop9b() -> TestResult {
        // NMOS-R, "Grounded", S/D Swapped
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![
            Comp::R(1e-3, Num(0), Gnd),
            Comp::Mos0(Mos0i {
                name: s("m"),
                mos_type: MosType::NMOS,
                ports: MosPorts {
                    g: Num(0),
                    d: Gnd,
                    s: Num(0),
                    b: Gnd,
                },
            }),
        ]);

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], 0.0);
        Ok(())
    }

    /// PMOS-R, "Grounded"
    #[test]
    fn test_dcop9c() -> TestResult {
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![
            Comp::R(1e-3, Num(0), Gnd),
            Comp::Mos0(Mos0i {
                name: s("m"),
                mos_type: MosType::PMOS,
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
        ]);

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], 0.0);
        Ok(())
    }

    /// PMOS-R, "Grounded", S/D Swapped
    #[test]
    fn test_dcop9d() -> TestResult {
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![
            Comp::R(1e-3, Num(0), Gnd),
            Comp::Mos0(Mos0i {
                name: s("m"),
                mos_type: MosType::PMOS,
                ports: MosPorts {
                    g: Num(0),
                    d: Gnd,
                    s: Num(0),
                    b: Gnd,
                },
            }),
        ]);

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], 0.0);
        Ok(())
    }

    /// NMOS-R Inverter
    #[test]
    fn test_dcop10() -> TestResult {
        let ckt = Ckt::from_comps(vec![
            Comp::R(1e-3, n("vdd"), n("d")),
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::Mos0(Mos0i {
                name: s("m"),
                mos_type: MosType::NMOS,
                ports: MosPorts {
                    g: n("vdd"),
                    d: n("d"),
                    s: Gnd,
                    b: Gnd,
                },
            }),
        ]);

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], 1.0);
        assert!(soln[1] < 50e-3);
        assert!((soln[2] + 1e-3).abs() < 0.1e-3);
        Ok(())
    }

    /// PMOS-R Inverter
    #[test]
    fn test_dcop10b() -> TestResult {
        let ckt = Ckt::from_comps(vec![
            Comp::R(1e-3, n("g"), n("d")),
            Comp::vdc("v1", -1.0, n("g"), Gnd),
            Comp::Mos0(Mos0i {
                name: s("m"),
                mos_type: MosType::PMOS,
                ports: MosPorts {
                    g: n("g"),
                    d: n("d"),
                    s: Gnd,
                    b: Gnd,
                },
            }),
        ]);

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

        let ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::R(1e-9, n("d"), Gnd), // "gmin"
            Comp::Mos0(Mos0i {
                name: s("p"),
                mos_type: PMOS,
                ports: MosPorts {
                    g: n("vdd"),
                    d: n("d"),
                    s: n("vdd"),
                    b: n("vdd"),
                },
            }),
            Comp::Mos0(Mos0i {
                name: s("n"),
                mos_type: NMOS,
                ports: MosPorts {
                    g: n("vdd"),
                    d: n("d"),
                    s: Gnd,
                    b: Gnd,
                },
            }),
        ]);

        let soln = dcop(ckt)?;
        assert(soln.values).eq(vec![0.0, 1.0, 0.0])?;
        Ok(())
    }
    /// Mos0 CMOS Inverter DC-Op, Vin=Vss
    #[test]
    fn test_dcop11b() -> TestResult {
        use MosType::{NMOS, PMOS};
        use NodeRef::Gnd;

        let ckt = Ckt::from_comps(vec![
            Comp::Mos0(Mos0i {
                name: s("p"),
                mos_type: PMOS,
                ports: MosPorts {
                    g: Gnd,
                    d: n("d"),
                    s: n("vdd"),
                    b: n("vdd"),
                },
            }),
            Comp::Mos0(Mos0i {
                name: s("n"),
                mos_type: NMOS,
                ports: MosPorts {
                    g: Gnd,
                    d: n("d"),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::R(1e-9, n("d"), n("vdd")), // "gmin"
        ]);

        let soln = dcop(ckt)?;
        assert(soln.values).eq(vec![1.0, 1.0, 0.0])?;
        Ok(())
    }
    /// DCOP, Several Series CMOS Inverters
    #[test]
    fn test_dcop12() -> TestResult {
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![
            Comp::R(1e-9, Num(0), Gnd),
            Comp::R(1e-9, Num(1), Gnd),
            Comp::R(1e-9, Num(2), Gnd),
            Comp::R(1e-9, Num(3), Gnd),
            Comp::R(1e-9, Num(4), Gnd),
            Comp::Mos0(Mos0i {
                name: s("p1"),
                mos_type: MosType::PMOS,
                ports: MosPorts {
                    g: Num(0),
                    d: Num(1),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos0(Mos0i {
                name: s("n1"),
                mos_type: MosType::NMOS,
                ports: MosPorts {
                    g: Num(0),
                    d: Num(1),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::Mos0(Mos0i {
                name: s("p2"),
                mos_type: MosType::PMOS,
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos0(Mos0i {
                name: s("n2"),
                mos_type: MosType::NMOS,
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::Mos0(Mos0i {
                name: s("p3"),
                mos_type: MosType::PMOS,
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos0(Mos0i {
                name: s("n3"),
                mos_type: MosType::NMOS,
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::Mos0(Mos0i {
                name: s("p4"),
                mos_type: MosType::PMOS,
                ports: MosPorts {
                    g: Num(3),
                    d: Num(4),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos0(Mos0i {
                name: s("n4"),
                mos_type: MosType::NMOS,
                ports: MosPorts {
                    g: Num(3),
                    d: Num(4),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::vdc("v1", 1.0, Num(0), Gnd),
        ]);

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
            Comp::R(1e-3, Num(1), Num(0)),
            Comp::C(1e-9, Num(1), Gnd),
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
            Comp::C(1e-9, n("i"), n("o")),
            Comp::R(1e-3, n("o"), Gnd),
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
            Comp::R(1e-3, n("inp"), n("out")),
            Comp::C(1e-9, n("out"), Gnd),
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
        let ckt = Ckt::from_comps(vec![Comp::I(1e-3, Num(0), Gnd), Comp::C(4e-12, Num(0), Gnd)]);

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
        let ckt = Ckt::from_comps(vec![Comp::I(1e-6, Num(0), Gnd), Comp::C(100e-9, Num(0), Gnd)]);

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
        let ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, Num(0), Gnd),
            Comp::R(1e-3, Num(0), Gnd),
            Comp::Mos0(Mos0i {
                name: s("p1"),
                mos_type: MosType::PMOS,
                ports: MosPorts {
                    g: Num(3),
                    d: Num(1),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos0(Mos0i {
                name: s("n1"),
                mos_type: MosType::NMOS,
                ports: MosPorts {
                    g: Num(3),
                    d: Num(1),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::R(1e-9, Num(1), Gnd),
            Comp::C(c, Num(1), Gnd),
            Comp::Mos0(Mos0i {
                name: s("p2"),
                mos_type: MosType::PMOS,
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos0(Mos0i {
                name: s("n2"),
                mos_type: MosType::NMOS,
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::R(1e-9, Num(2), Gnd),
            Comp::C(c, Num(2), Gnd),
            Comp::Mos0(Mos0i {
                name: s("p3"),
                mos_type: MosType::PMOS,
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos0(Mos0i {
                name: s("n3"),
                mos_type: MosType::NMOS,
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::R(1e-9, Num(3), Gnd),
            Comp::C(c, Num(3), Gnd),
        ]);
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
            Comp::Mos1(Mos1i {
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
            Comp::Mos1(Mos1i {
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
            Comp::Mos1(Mos1i {
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
            Comp::Mos1(Mos1i {
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
            Comp::R(1e-4, Num(1), Gnd),
        ]);
        // Define our models & params
        add_mos1_defaults(&mut ckt);

        dcop(ckt)?;
        // FIXME: checks on solution
        Ok(())
    }
    /// Mos1 CMOS Ring Oscillator Dc Op
    #[test]
    fn test_mos1_cmos_ro_dcop() -> TestResult {
        let mut ckt = Ckt::from_comps(vec![
            Comp::R(1e-9, Num(0), Gnd), // "gmin"-ish
            Comp::Mos1(Mos1i {
                name: s("p1"),
                model: "pmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos1(Mos1i {
                name: s("n1"),
                model: "nmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::Mos1(Mos1i {
                name: s("p2"),
                model: "pmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos1(Mos1i {
                name: s("n2"),
                model: "nmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::Mos1(Mos1i {
                name: s("p3"),
                model: "pmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(3),
                    d: Num(1),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos1(Mos1i {
                name: s("n3"),
                model: "nmos".into(),
                params: "default".into(),
                ports: MosPorts {
                    g: Num(3),
                    d: Num(1),
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
        assert(soln.get("0")?).eq(1.0)?;
        for k in ["1", "2", "3"].into_iter() {
            assert(soln.get(*k)?).gt(0.45)?;
            assert(soln.get(*k)?).lt(0.55)?;
        }
        Ok(())
    }
    /// Mos1 CMOS Ring Oscillator Tran
    #[test]
    fn test_mos1_cmos_ro_tran() -> TestResult {
        let mut ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::Mos1(Mos1i {
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
            Comp::Mos1(Mos1i {
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
            Comp::Mos1(Mos1i {
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
            Comp::Mos1(Mos1i {
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
            Comp::Mos1(Mos1i {
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
            Comp::Mos1(Mos1i {
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
        add_mos1_defaults(&mut ckt);
        // Simulate
        let opts = TranOptions {
            tstep: 1e-11,
            tstop: 1e-8,
            ic: vec![(Num(1), 0.0)],
        };
        let soln = tran(ckt, opts)?;
        //to_file(&soln, "test_mos1_cmos_ro_tran.json"); // Writes new golden data
        // Checks
        let golden = load_golden("test_mos1_cmos_ro_tran.json");
        assert(&soln.map).isclose(golden, 1e-6)?;
        Ok(())
    }

    // Mos1 NMOS-R Oscillator Tran
    #[test]
    fn test_mos1_nmos_ro_tran() -> TestResult {
        let gl = 1e-6;
        let mut ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::Mos1(Mos1i {
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
            Comp::Mos1(Mos1i {
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
            Comp::Mos1(Mos1i {
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
            Comp::R(gl, Num(1), n("vdd")),
            Comp::R(gl, Num(2), n("vdd")),
            Comp::R(gl, Num(3), n("vdd")),
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
            Comp::Mos1(Mos1i {
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
            Comp::Mos1(Mos1i {
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
            Comp::Mos1(Mos1i {
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
            Comp::R(gl, Num(1), Gnd),
            Comp::R(gl, Num(2), Gnd),
            Comp::R(gl, Num(3), Gnd),
            Comp::C(cl, Num(1), Gnd),
            Comp::C(cl, Num(2), Gnd),
            Comp::C(cl, Num(3), Gnd),
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

    /// Mos1 PMOS-R Amp, Tran Initial Condition Decay
    #[test]
    fn test_mos1_pmos_rload_tran() -> TestResult {
        let gl = 1e-6;
        let mut ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, n("vdd"), Gnd),
            Comp::Mos1(Mos1i {
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
            Comp::R(gl, n("inp"), n("vdd")),
            Comp::R(gl, n("out"), Gnd),
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
            Comp::Mos1(Mos1i {
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
            Comp::R(gl, n("g"), Gnd),
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
    /// Mos1 NMOS-R Tran
    #[test]
    fn test_mos1_nmos_rg_tran() -> TestResult {
        let gl = 1e-6;
        let mut ckt = Ckt::from_comps(vec![
            Comp::Mos1(Mos1i {
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
            Comp::R(gl, n("g"), Gnd),
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

    #[test]
    fn test_ac1() -> TestResult {
        let ckt = Ckt::from_comps(vec![Comp::R(1.0, Num(0), Gnd)]);
        ac(ckt, AcOptions::default())?;
        // FIXME: checks on solution
        Ok(())
    }

    #[test]
    fn test_ac2() -> TestResult {
        use crate::circuit::Vs;
        use Comp::{C, R, V};
        let ckt = Ckt::from_comps(vec![
            R(1e-3, Num(0), Num(1)),
            C(1e-9, Num(1), Gnd),
            V(Vs {
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
        let ckt = Ckt::from_comps(vec![
            Comp::R(1e-3, Num(0), Num(1)),
            Comp::C(1e-9, Num(1), Gnd),
            Comp::vdc("v1", 1.0, Num(0), Gnd),
            Comp::Mos0(Mos0i {
                name: s("m"),
                mos_type: MosType::NMOS,
                ports: MosPorts {
                    g: Num(1),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
        ]);
        ac(ckt, AcOptions::default())?;
        // FIXME: checks on solution
        Ok(())
    }

    /// NMOS Common-Source Amp
    #[test]
    fn test_ac4() -> TestResult {
        let mut ckt = Ckt::from_comps(vec![
            Comp::C(1e-9, n("d"), Gnd),
            Comp::Mos1(Mos1i {
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
            Comp::V(Vs {
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
        use crate::circuit::Vs;

        let mut ckt = Ckt::from_comps(vec![
            Comp::V(Vs {
                name: s("vd"),
                vdc: 0.5,
                acm: 1.0,
                p: Num(0),
                n: Gnd,
            }),
            Comp::Mos1(Mos1i {
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
}
