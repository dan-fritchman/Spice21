pub mod analysis;
pub mod assert;
pub mod comps;
pub mod proto;
pub mod sparse21;
pub mod spresult;

use num::traits::NumAssignOps;
use num::{Complex, Float, Num, One, Zero};
use std::cmp::PartialOrd;
use std::fmt;

// This long list of traits describes our required behavior for numeric types.
pub trait SpNum:
    Clone + Copy + NumAssignOps + Zero + Num + Abs + fmt::Display + fmt::Debug
{
}

impl<T> SpNum for T where
    T: Clone + Copy + NumAssignOps + Zero + Num + Abs + fmt::Display + fmt::Debug
{
}

/// Absolute Value Trait for numeric types
pub trait Abs {
    fn absv(&self) -> f64;
}

impl Abs for f64 {
    fn absv(&self) -> f64 {
        self.abs()
    }
}

impl Abs for Complex<f64> {
    fn absv(&self) -> f64 {
        self.norm()
    }
}

/// "Integration" Tests
#[cfg(test)]
mod tests {
    use super::analysis::*;
    use super::assert::*;
    use super::comps::*;
    use super::proto::*;
    use super::sparse21::*;
    use super::spresult::*;
    use super::*;

    /// Create a very basic Circuit
    fn parse_ckt() -> CktParse {
        CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(1e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Gnd),
            ],
        }
    }

    #[test]
    fn test_ckt_parse() -> TestResult {
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(1e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Gnd),
            ],
        };
        Ok(())
    }

    #[test]
    fn test_dcop1() -> TestResult {
        let ckt = CktParse {
            nodes: 1,
            comps: vec![CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Gnd)],
        };

        let soln = dcop(ckt)?;
        assert_eq!(soln, vec![0.0]);
        Ok(())
    }

    #[test]
    fn test_dcop2() -> TestResult {
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(1e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Gnd),
            ],
        };

        let soln = dcop(ckt)?;
        assert_eq!(soln, vec![1.0,]);
        Ok(())
    }

    #[test]
    fn test_dcop3() -> TestResult {
        // I - R - R divider
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::R(1e-3, NodeRef::Num(1), NodeRef::Num(0)),
                CompParse::I(1e-3, NodeRef::Num(1), NodeRef::Gnd),
            ],
        };

        let soln = dcop(ckt)?;
        assert!((soln[0] - 1.0).abs() < 1e-4);
        assert!((soln[1] - 2.0).abs() < 1e-4);
        Ok(())
    }

    #[test]
    fn test_dcop4() -> TestResult {
        // V - R - R divider
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, NodeRef::Num(1), NodeRef::Gnd),
                CompParse::R(2e-3, NodeRef::Num(1), NodeRef::Num(0)),
                CompParse::R(2e-3, NodeRef::Num(0), NodeRef::Gnd),
            ],
        };

        let soln = dcop(ckt)?;
        assert_eq!(soln, vec![0.5, 1.0, -1e-3]);
        Ok(())
    }

    #[test]
    fn test_dcop5() -> TestResult {
        // I - R - Diode
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::D(2e-16, 25e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::I(1e-3, NodeRef::Num(0), NodeRef::Gnd),
            ],
        };

        let soln = dcop(ckt)?;
        assert!((soln[0] - 0.7).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_dcop5b() -> TestResult {
        // V - R - Diode
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::D(2e-16, 25e-3, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::R(1e-3, NodeRef::Num(0), NodeRef::Num(1)),
                CompParse::V(1.0, NodeRef::Num(1), NodeRef::Gnd),
            ],
        };

        let soln = dcop(ckt)?;
        assert!((soln[0] - 0.7).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_dcop6() -> TestResult {
        // NMOS Char
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::V(1.0, NodeRef::Num(1), NodeRef::Gnd),
                CompParse::Mos0(MosType::NMOS, Num(0), Num(1), Gnd, Gnd),
            ],
        };

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
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(-1.0, NodeRef::Num(0), NodeRef::Gnd),
                CompParse::V(-1.0, NodeRef::Num(1), NodeRef::Gnd),
                CompParse::Mos0(MosType::PMOS, Num(0), Num(1), Gnd, Gnd),
            ],
        };

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
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(5e-3, Num(0), Gnd),
                CompParse::Mos0(MosType::NMOS, Num(0), Num(0), Gnd, Gnd),
                CompParse::R(1e-12, Num(0), Gnd), // "gmin"
            ],
        };

        let soln = dcop(ckt)?;
        assert!((soln[0] - 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_diode_nmos_tran() -> TestResult {
        // Diode NMOS Tran
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(5e-3, Num(0), Gnd),
                CompParse::Mos0(MosType::NMOS, Num(0), Num(0), Gnd, Gnd),
                CompParse::R(1e-12, Num(0), Gnd), // "gmin"
            ],
        };

        let opts = TranOptions {
            tstep: 1e-12,
            tstop: 100e-12,
        };
        let soln = tran(ckt, opts)?;
        for point in soln.iter() {
            assert!((point[0] - 0.697).abs() < 1e-3);
        }
        Ok(())
    }

    #[test]
    fn test_dcop8b() -> TestResult {
        // Diode NMOS, S/D Swapped
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(5e-3, Num(0), Gnd),
                CompParse::Mos0(MosType::NMOS, Num(0), Gnd, Num(0), Gnd),
                CompParse::R(1e-12, Num(0), Gnd), // "gmin"
            ],
        };

        let soln = dcop(ckt)?;
        assert!((soln[0] - 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_diode_pmos_dcop() -> TestResult {
        // Diode PMOS
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(-5e-3, Num(0), Gnd),
                CompParse::Mos0(MosType::PMOS, Num(0), Num(0), Gnd, Gnd),
                CompParse::R(1e-12, Num(0), Gnd), // "gmin"
            ],
        };

        let soln = dcop(ckt)?;
        assert!((soln[0] + 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_diode_pmos_tran() -> TestResult {
        // Diode PMOS Tran
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(-5e-3, Num(0), Gnd),
                CompParse::Mos0(MosType::PMOS, Num(0), Num(0), Gnd, Gnd),
                CompParse::R(1e-12, Num(0), Gnd), // "gmin"
            ],
        };

        let opts = TranOptions {
            tstep: 1e-12,
            tstop: 100e-12,
        };
        let soln = tran(ckt, opts)?;
        for point in soln.iter() {
            assert!((point[0] + 0.697).abs() < 1e-3);
        }
        Ok(())
    }

    #[test]
    fn test_dcop8d() -> TestResult {
        // Diode PMOS, S/D Swapped
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(-5e-3, Num(0), Gnd),
                CompParse::Mos0(MosType::PMOS, Num(0), Gnd, Num(0), Gnd),
                CompParse::R(1e-12, Num(0), Gnd), // "gmin"
            ],
        };

        let soln = dcop(ckt)?;
        assert!((soln[0] + 0.697).abs() < 1e-3);
        Ok(())
    }

    #[test]
    fn test_dcop9() -> TestResult {
        // NMOS-R, "Grounded"
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::R(1e-3, Num(0), Gnd),
                CompParse::Mos0(MosType::NMOS, Num(0), Num(0), Gnd, Gnd),
            ],
        };

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], 0.0);
        Ok(())
    }

    #[test]
    fn test_dcop9b() -> TestResult {
        // NMOS-R, "Grounded", S/D Swapped
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::R(1e-3, Num(0), Gnd),
                CompParse::Mos0(MosType::NMOS, Num(0), Gnd, Num(0), Gnd),
            ],
        };

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], 0.0);
        Ok(())
    }

    #[test]
    fn test_dcop9c() -> TestResult {
        // PMOS-R, "Grounded"
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::R(1e-3, Num(0), Gnd),
                CompParse::Mos0(MosType::PMOS, Num(0), Num(0), Gnd, Gnd),
            ],
        };

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], 0.0);
        Ok(())
    }

    #[test]
    fn test_dcop9d() -> TestResult {
        // PMOS-R, "Grounded", S/D Swapped
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::R(1e-3, Num(0), Gnd),
                CompParse::Mos0(MosType::PMOS, Num(0), Gnd, Num(0), Gnd),
            ],
        };

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], 0.0);
        Ok(())
    }

    #[test]
    fn test_dcop10() -> TestResult {
        // NMOS-R Inverter
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::R(1e-3, Num(1), Num(0)),
                CompParse::Mos0(MosType::NMOS, Num(0), Num(1), Gnd, Gnd),
            ],
        };

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], 1.0);
        assert!(soln[1] < 50e-3);
        assert!((soln[2] + 1e-3).abs() < 0.1e-3);
        Ok(())
    }

    #[test]
    fn test_dcop10b() -> TestResult {
        // PMOS-R Inverter
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(-1.0, Num(0), Gnd),
                CompParse::R(1e-3, Num(1), Num(0)),
                CompParse::Mos0(MosType::PMOS, Num(0), Num(1), Gnd, Gnd),
            ],
        };

        let soln = dcop(ckt)?;
        assert_eq!(soln[0], -1.0);
        assert!(soln[1].abs() < 50e-3);
        assert!((soln[2] - 1e-3).abs() < 0.1e-3);
        Ok(())
    }

    #[test]
    fn test_dcop11() -> TestResult {
        // CMOS Inverter
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::Mos0(MosType::PMOS, Num(0), Num(1), Num(0), Num(0)),
                CompParse::Mos0(MosType::NMOS, Num(0), Num(1), Gnd, Gnd),
                CompParse::R(1e-9, Num(1), Gnd),
            ],
        };

        let soln = dcop(ckt)?;
        assert_eq!(soln, vec![1.0, 0.0, 0.0]);
        Ok(())
    }

    #[test]
    fn test_dcop11b() -> TestResult {
        // CMOS Inverter
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::Mos0(MosType::PMOS, Gnd, Num(1), Num(0), Num(0)),
                CompParse::Mos0(MosType::NMOS, Gnd, Num(1), Gnd, Gnd),
                CompParse::R(1e-9, Num(1), Num(0)),
            ],
        };

        let soln = dcop(ckt)?;
        assert_eq!(soln, vec![1.0, 1.0, 0.0]);
        Ok(())
    }

    #[test]
    fn test_dcop12() -> TestResult {
        // Several CMOS Inverters
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 5,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::Mos0(MosType::PMOS, Num(0), Num(1), Num(0), Num(0)),
                CompParse::Mos0(MosType::NMOS, Num(0), Num(1), Gnd, Gnd),
                CompParse::R(1e-9, Num(1), Gnd),
                CompParse::Mos0(MosType::PMOS, Num(1), Num(2), Num(0), Num(0)),
                CompParse::Mos0(MosType::NMOS, Num(1), Num(2), Gnd, Gnd),
                CompParse::R(1e-9, Num(2), Gnd),
                CompParse::Mos0(MosType::PMOS, Num(2), Num(3), Num(0), Num(0)),
                CompParse::Mos0(MosType::NMOS, Num(2), Num(3), Gnd, Gnd),
                CompParse::R(1e-9, Num(3), Gnd),
                CompParse::Mos0(MosType::PMOS, Num(3), Num(4), Num(0), Num(0)),
                CompParse::Mos0(MosType::NMOS, Num(3), Num(4), Gnd, Gnd),
                CompParse::R(1e-9, Num(4), Gnd),
            ],
        };

        let soln = dcop(ckt)?;
        assert(soln[0]).eq(1.0)?;
        assert!(soln[1].abs() < 1e-3);
        assert!((soln[2] - 1.0).abs() < 1e-3);
        assert!(soln[3].abs() < 1e-3);
        assert!((soln[4] - 1.0).abs() < 1e-3);
        assert!(soln[5].abs() < 1e-6);
        Ok(())
    }

    #[test]
    fn test_dcop13() -> TestResult {
        // RC Low-Pass Filter
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::R(1e-3, Num(1), Num(0)),
                CompParse::C(1e-9, Num(1), Gnd),
            ],
        };

        let soln = dcop(ckt)?;
        assert_eq!(soln, vec![1.0, 1.0, 0.0]);
        Ok(())
    }

    #[test]
    fn test_dcop13b() -> TestResult {
        // RC High-Pass Filter
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::R(1e-3, Num(1), Gnd),
                CompParse::C(1e-9, Num(1), Num(0)),
            ],
        };

        let soln = dcop(ckt)?;
        assert_eq!(soln, vec![1.0, 0.0, 0.0]);
        Ok(())
    }

    #[test]
    fn test_tran1() -> TestResult {
        // RC Low-Pass Filter
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::R(1e-3, Num(1), Num(0)),
                CompParse::C(1e-9, Num(1), Gnd),
            ],
        };

        let opts = TranOptions {
            tstep: 1e-11,
            tstop: 100e-11,
        };
        let soln = tran(ckt, opts)?;
        for point in soln.into_iter() {
            assert(point).eq(vec![1.0, 1.0, 0.0])?;
        }
        Ok(())
    }

    #[test]
    fn test_tran2() -> TestResult {
        // I-C Integrator, with Initial Condition
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::I(5e-3, Num(0), Gnd),
                CompParse::C(1e-9, Num(0), Gnd),
            ],
        };

        let opts = TranOptions {
            tstep: 1e-9,
            tstop: 100e-9,
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
        let ckt = CktParse {
            nodes: 4,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::R(1e-3, Num(0), Gnd),
                CompParse::Mos0(MosType::PMOS, Num(3), Num(1), Num(0), Num(0)),
                CompParse::Mos0(MosType::NMOS, Num(3), Num(1), Gnd, Gnd),
                CompParse::R(1e-5, Num(1), Gnd),
                CompParse::C(c, Num(1), Gnd),
                CompParse::Mos0(MosType::PMOS, Num(1), Num(2), Num(0), Num(0)),
                CompParse::Mos0(MosType::NMOS, Num(1), Num(2), Gnd, Gnd),
                CompParse::R(1e-5, Num(2), Gnd),
                CompParse::C(c, Num(2), Gnd),
                CompParse::Mos0(MosType::PMOS, Num(2), Num(3), Num(0), Num(0)),
                CompParse::Mos0(MosType::NMOS, Num(2), Num(3), Gnd, Gnd),
                CompParse::R(1e-5, Num(3), Gnd),
                CompParse::C(c, Num(3), Gnd),
            ],
        };

        let opts = TranOptions {
            tstep: 1e-12,
            tstop: 100e-12,
        };
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
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::Mos1(model, params, Num(0), Num(0), Gnd, Gnd),
                CompParse::V(1.0, Num(0), Gnd),
            ],
        };
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
        let ckt = CktParse {
            nodes: 1,
            comps: vec![
                CompParse::Mos1(model, params, Num(0), Num(0), Gnd, Gnd),
                CompParse::V(1.0, Num(0), Gnd),
            ],
        };

        let opts = TranOptions {
            tstep: 1e-9,
            tstop: 100e-9,
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
        let c = 1e-10;
        let ckt = CktParse {
            nodes: 2,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::Mos1(pmos.clone(), params, Gnd, Num(1), Num(0), Num(0)),
                CompParse::Mos1(nmos.clone(), params, Gnd, Num(1), Gnd, Gnd),
                CompParse::R(1e-4, Num(1), Gnd),
            ],
        };
        let soln = dcop(ckt)?;
        Ok(())
    }

    #[test]
    fn test_mos1_ro_dcop() -> TestResult {
        // Mos1 Ring Oscillator Dc Op

        let nmos = Mos1Model::default();
        let pmos = Mos1Model {
            mos_type: MosType::PMOS,
            ..Mos1Model::default()
        };
        let params = Mos1InstanceParams::default();
        use NodeRef::{Gnd, Num};
        let ckt = CktParse {
            nodes: 4,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::Mos1(pmos.clone(), params, Num(3), Num(1), Num(0), Num(0)),
                CompParse::Mos1(nmos.clone(), params, Num(3), Num(1), Gnd, Gnd),
                CompParse::Mos1(pmos.clone(), params, Num(1), Num(2), Num(0), Num(0)),
                CompParse::Mos1(nmos.clone(), params, Num(1), Num(2), Gnd, Gnd),
                CompParse::Mos1(pmos.clone(), params, Num(2), Num(3), Num(0), Num(0)),
                CompParse::Mos1(nmos.clone(), params, Num(2), Num(3), Gnd, Gnd),
            ],
        };

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
        let ckt = CktParse {
            nodes: 4,
            comps: vec![
                CompParse::V(1.0, Num(0), Gnd),
                CompParse::Mos1(pmos.clone(), params, Num(3), Num(1), Num(0), Num(0)),
                CompParse::Mos1(nmos.clone(), params, Num(3), Num(1), Gnd, Gnd),
                CompParse::C(c, Num(1), Gnd),
                CompParse::R(1e-9, Num(1), Gnd),
                CompParse::Mos1(pmos.clone(), params, Num(1), Num(2), Num(0), Num(0)),
                CompParse::Mos1(nmos.clone(), params, Num(1), Num(2), Gnd, Gnd),
                CompParse::C(c, Num(2), Gnd),
                CompParse::R(1e-9, Num(2), Gnd),
                CompParse::Mos1(pmos.clone(), params, Num(2), Num(3), Num(0), Num(0)),
                CompParse::Mos1(nmos.clone(), params, Num(2), Num(3), Gnd, Gnd),
                CompParse::C(c, Num(3), Gnd),
                CompParse::R(1e-9, Num(3), Gnd),
            ],
        };

        let opts = TranOptions {
            tstep: 5e-9,
            tstop: 1000e-9,
        };
        let mut tran = Tran::new(ckt, opts);
        tran.ic(Num(1), 0.0);
        let soln = tran.solve()?;
        // FIXME: dream up some checks
        Ok(())
    }
}
