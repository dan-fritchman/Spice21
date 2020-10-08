/// "Integration" Tests
#[cfg(test)]
mod tests {
    use crate::analysis::*;
    use crate::assert::*;
    use crate::circuit::NodeRef::{Gnd, Name, Num};
    use crate::circuit::*;
    use crate::comps::*;
    use crate::spresult::*;

    /// Create a very basic Circuit
    #[test]
    fn test_ckt_parse() -> TestResult {
        let ckt = Ckt::from_comps(vec![
            Comp::I(1e-3, NodeRef::Name(s("0")), NodeRef::Gnd),
            Comp::R(1e-3, NodeRef::Name(s("0")), NodeRef::Gnd),
        ]);
        Ok(())
    }

    #[test]
    fn test_dcop1() -> TestResult {
        let ckt = Ckt::from_comps(vec![Comp::R(1e-3, NodeRef::Num(0), NodeRef::Gnd)]);

        let soln = dcop(ckt)?;
        assert_eq!(soln.values, vec![0.0]);
        Ok(())
    }

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

    #[test]
    fn test_dcop3() -> TestResult {
        // I - R - R divider
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
        let ckt = Ckt::from_comps(vec![Comp::vdc("v1", 1.0, n("vdd"), Gnd), R(2e-3, n("vdd"), n("div")), R(2e-3, n("div"), Gnd)]);
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

    #[test]
    fn test_dcop9c() -> TestResult {
        // PMOS-R, "Grounded"
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

    #[test]
    fn test_dcop9d() -> TestResult {
        // PMOS-R, "Grounded", S/D Swapped
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

    #[test]
    fn test_dcop12() -> TestResult {
        // Several CMOS Inverters
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
        use NodeRef::{Gnd, Num};
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
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![
            Comp::R(1e-3, Num(0), Num(1)),
            Comp::C(1e-9, Num(1), Gnd),
            Comp::vdc("v1", 1.0, Num(0), Gnd),
        ]);

        let opts = TranOptions {
            tstep: 1e-11,
            tstop: 100e-11,
            ..Default::default()
        };
        let soln = tran(ckt, opts)?;
        for point in soln.data.into_iter() {
            assert(point).eq(vec![1.0, 1.0, 0.0])?;
        }
        Ok(())
    }

    #[test]
    fn test_tran2() -> TestResult {
        // I-C Integrator, with Initial Condition
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![Comp::I(5e-3, Num(0), Gnd), Comp::C(1e-9, Num(0), Gnd)]);

        let opts = TranOptions {
            tstep: 1e-9,
            tstop: 100e-9,
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

    #[test]
    fn test_tran3() -> TestResult {
        // Ring Oscillator
        use NodeRef::{Gnd, Num};
        let c = 1e-10;
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
            Comp::R(1e-5, Num(1), Gnd),
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
            Comp::R(1e-5, Num(2), Gnd),
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
            Comp::R(1e-5, Num(3), Gnd),
            Comp::C(c, Num(3), Gnd),
        ]);

        let opts = TranOptions {
            tstep: 1e-12,
            tstop: 100e-12,
            ..Default::default()
        };
        let mut models = ModelCache::new();
        let mut tran = Tran::new(ckt, opts);
        tran.ic(Num(1), 0.0);
        let soln = tran.solve()?;
        // FIXME: dream up some checks
        Ok(())
    }

    #[test]
    fn test_mos1_op() -> TestResult {
        let model = Mos1Model::default();
        let params = Mos1InstanceParams::default();
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![
            Comp::Mos1(Mos1i {
                name: s("m"),
                model: model,
                params: params,
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::vdc("v1", 1.0, Num(0), Gnd),
        ]);
        let soln = dcop(ckt)?;
        assert(soln[0]).eq(1.0)?;
        assert(soln[1]).lt(0.0)?;
        assert(soln[1]).gt(-1e-3)?;
        Ok(())
    }

    #[test]
    fn test_mos1_tran() -> TestResult {
        let model = Mos1Model::default();
        let params = Mos1InstanceParams::default();
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![
            Comp::Mos1(Mos1i {
                name: s("m"),
                model: model,
                params: params,
                ports: MosPorts {
                    g: Num(0),
                    d: Num(0),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::vdc("v1", 1.0, Num(0), Gnd),
        ]);

        let opts = TranOptions {
            tstep: 1e-9,
            tstop: 100e-9,
            ..Default::default()
        };
        let soln = tran(ckt, opts)?;
        for k in 1..soln.len() {
            assert(soln[k][0]).eq(1.0)?;
            assert(soln[k][1]).lt(0.0)?;
            assert(soln[k][1]).gt(-1e-3)?;
        }
        Ok(())
    }

    #[test]
    fn test_mos1_inv_dcop() -> TestResult {
        // Mos1 Inverter DCOP

        let nmos = Mos1Model::default();
        let pmos = Mos1Model {
            mos_type: MosType::PMOS,
            ..Mos1Model::default()
        };
        let params = Mos1InstanceParams::default();
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, Num(0), Gnd),
            Comp::Mos1(Mos1i {
                name: s("p"),
                model: pmos.clone(),
                params: params,
                ports: MosPorts {
                    g: Gnd,
                    d: Num(1),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos1(Mos1i {
                name: s("n"),
                model: nmos.clone(),
                params: params,
                ports: MosPorts {
                    g: Gnd,
                    d: Num(1),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::R(1e-4, Num(1), Gnd),
        ]);
        let soln = dcop(ckt)?;
        Ok(())
    }

    /// Mos1 Ring Oscillator Dc Op
    #[test]
    fn test_mos1_ro_dcop() -> TestResult {
        let nmos = Mos1Model::default();
        let pmos = Mos1Model {
            mos_type: MosType::PMOS,
            ..Mos1Model::default()
        };
        let params = Mos1InstanceParams::default();
        use NodeRef::{Gnd, Num};
        let ckt = Ckt::from_comps(vec![
            Comp::R(1e-9, Num(0), Gnd), // "gmin"-ish
            Comp::Mos1(Mos1i {
                name: s("p1"),
                model: pmos.clone(),
                params: params,
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos1(Mos1i {
                name: s("n1"),
                model: nmos.clone(),
                params: params,
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::Mos1(Mos1i {
                name: s("p2"),
                model: pmos.clone(),
                params: params,
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos1(Mos1i {
                name: s("n2"),
                model: nmos.clone(),
                params: params,
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::Mos1(Mos1i {
                name: s("p3"),
                model: pmos.clone(),
                params: params,
                ports: MosPorts {
                    g: Num(3),
                    d: Num(1),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos1(Mos1i {
                name: s("n3"),
                model: nmos.clone(),
                params: params,
                ports: MosPorts {
                    g: Num(3),
                    d: Num(1),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::vdc("v1", 1.0, Num(0), Gnd),
        ]);

        let soln = dcop(ckt)?;
        assert(soln[0]).eq(1.0)?;
        for k in 1..3 {
            assert(soln[k]).gt(0.45)?;
            assert(soln[k]).lt(0.55)?;
        }
        Ok(())
    }

    #[test]
    fn test_mos1_ro_tran() -> TestResult {
        // Mos1 Ring Oscillator Tran

        let nmos = Mos1Model::default();
        let pmos = Mos1Model {
            mos_type: MosType::PMOS,
            ..Mos1Model::default()
        };
        let params = Mos1InstanceParams::default();
        use NodeRef::{Gnd, Num};
        let c = 1e-12;
        let ckt = Ckt::from_comps(vec![
            Comp::vdc("v1", 1.0, Num(0), Gnd),
            Comp::Mos1(Mos1i {
                name: s("p1"),
                model: pmos.clone(),
                params: params,
                ports: MosPorts {
                    g: Num(3),
                    d: Num(1),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos1(Mos1i {
                name: s("n1"),
                model: nmos.clone(),
                params: params,
                ports: MosPorts {
                    g: Num(3),
                    d: Num(1),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::C(c, Num(1), Gnd),
            Comp::R(1e-9, Num(1), Gnd),
            Comp::Mos1(Mos1i {
                name: s("p2"),
                model: pmos.clone(),
                params: params,
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos1(Mos1i {
                name: s("n2"),
                model: nmos.clone(),
                params: params,
                ports: MosPorts {
                    g: Num(1),
                    d: Num(2),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::C(c, Num(2), Gnd),
            Comp::R(1e-9, Num(2), Gnd),
            Comp::Mos1(Mos1i {
                name: s("p3"),
                model: pmos.clone(),
                params: params,
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Num(0),
                    b: Num(0),
                },
            }),
            Comp::Mos1(Mos1i {
                name: s("n3"),
                model: nmos.clone(),
                params: params,
                ports: MosPorts {
                    g: Num(2),
                    d: Num(3),
                    s: Gnd,
                    b: Gnd,
                },
            }),
            Comp::C(c, Num(3), Gnd),
            Comp::R(1e-9, Num(3), Gnd),
        ]);

        let opts = TranOptions {
            tstep: 5e-9,
            tstop: 1000e-9,
            ..Default::default()
        };
        let mut tran = Tran::new(ckt, opts);
        tran.ic(Num(1), 0.0);
        let soln = tran.solve()?;
        // FIXME: dream up some checks
        Ok(())
    }
}
