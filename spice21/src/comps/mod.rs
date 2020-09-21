//!
//! # Spice21 Component Solvers
//!
//! Primary `Component` trait and basic implementations
//!
use enum_dispatch::enum_dispatch;
use num::Complex;

use std::ops::{Index, IndexMut};

use super::analysis::{AnalysisInfo, Stamps, VarIndex, Variables};
use super::sparse21::{Eindex, Matrix};
use crate::{SpNum, SpResult};

use diode::Diode0;
pub mod diode;
pub(crate) use diode::{Diode, DiodeInstParams, DiodeModel};
pub mod mos;
pub(crate) use mos::*;
pub mod bsim4;
pub(crate) use bsim4::*;

/// Constants
pub mod consts {
    use std::f64::{MAX_EXP as MAX_EXPI, MIN_EXP as MIN_EXPI};
    pub(crate) use std::f64::consts::PI;

    pub(crate) const KB: f64 = 1.3806226e-23;
    pub(crate) const Q: f64 = 1.6021918e-19;
    pub(crate) const KB_OVER_Q: f64 = KB / Q;
    pub(crate) const KELVIN_TO_C: f64 = 273.15;
    pub(crate) const TEMP_REF: f64 = KELVIN_TO_C + 27.0;
    pub(crate) const VT_REF: f64 = KB * TEMP_REF / Q;
    pub(crate) const SIO2_PERMITTIVITY: f64 = 3.9 * 8.854214871e-12;
    pub(crate) const SQRT2: f64 = 1.4142135624;
    pub(crate) const MAX_EXP: f64 = MAX_EXPI as f64;
    pub(crate) const MIN_EXP: f64 = MIN_EXPI as f64;
    pub(crate) const MAX_EXPL: f64 = 2.688117142e+43;
    pub(crate) const MIN_EXPL: f64 = 3.720075976e-44;
    pub(crate) const EXPL_THRESHOLD: f64 = 100.0;
    pub(crate) const EXP_THRESHOLD: f64 = 34.0;
    pub(crate) const EPS0: f64 = 8.85418e-12;
    pub(crate) const EPSSI: f64 = 1.03594e-10;
}

///
/// Spice21 ComponentSolver
/// The primary enumeration of component-types supported in simulation.
///
/// Shout-out `enum_dispatch` for "dynamic" dispatching the methods of
/// the `Component` trait to each of these types.
///
#[enum_dispatch]
pub(crate) enum ComponentSolver {
    Vsrc,
    Isrc,
    Capacitor,
    Resistor,
    Diode0,
    Diode,
    Mos0,
    Mos1,
}

#[enum_dispatch(ComponentSolver)]
pub(crate) trait Component {
    /// Commit operating-point guesses to internal state
    fn commit(&mut self) {}

    /// Update values of single-valued components
    /// FIXME: prob not for every Component
    fn update(&mut self, _val: f64) {}

    /// Validation of parameter Values etc.
    fn validate(&self) -> SpResult<()> {
        Ok(())
    }

    fn load_ac(
        &mut self,
        _guess: &Variables<Complex<f64>>,
        _an: &AnalysisInfo,
    ) -> Stamps<Complex<f64>> {
        Stamps::<Complex<f64>>::new()
    }
    fn load(&mut self, guess: &Variables<f64>, an: &AnalysisInfo) -> Stamps<f64>;
    fn create_matrix_elems<T: SpNum>(&mut self, mat: &mut Matrix<T>);
}

pub struct Vsrc {
    v: f64,
    acm: f64,
    p: Option<VarIndex>,
    n: Option<VarIndex>,
    ivar: VarIndex,
    pi: Option<Eindex>,
    ip: Option<Eindex>,
    ni: Option<Eindex>,
    in_: Option<Eindex>,
}

impl Vsrc {
    pub fn new(
        vdc: f64,
        acm: f64,
        p: Option<VarIndex>,
        n: Option<VarIndex>,
        ivar: VarIndex,
    ) -> Vsrc {
        Vsrc {
            v: vdc,
            acm,
            p,
            n,
            ivar,
            pi: None,
            ip: None,
            ni: None,
            in_: None,
        }
    }
}

impl Component for Vsrc {
    fn update(&mut self, val: f64) {
        self.v = val;
    }
    fn create_matrix_elems<T: SpNum>(&mut self, mat: &mut Matrix<T>) {
        self.pi = make_matrix_elem(mat, self.p, Some(self.ivar));
        self.ip = make_matrix_elem(mat, Some(self.ivar), self.p);
        self.ni = make_matrix_elem(mat, self.n, Some(self.ivar));
        self.in_ = make_matrix_elem(mat, Some(self.ivar), self.n);
    }
    fn load(&mut self, _guess: &Variables<f64>, _an: &AnalysisInfo) -> Stamps<f64> {
        return Stamps {
            g: vec![
                (self.pi, 1.0),
                (self.ip, 1.0),
                (self.ni, -1.0),
                (self.in_, -1.0),
            ],
            b: vec![(Some(self.ivar), self.v)],
        };
    }
    fn load_ac(
        &mut self,
        _guess: &Variables<Complex<f64>>,
        _an: &AnalysisInfo,
    ) -> Stamps<Complex<f64>> {
        return Stamps {
            g: vec![
                (self.pi, Complex::new(1.0, 0.0)),
                (self.ip, Complex::new(1.0, 0.0)),
                (self.ni, Complex::new(-1.0, 0.0)),
                (self.in_, Complex::new(-1.0, 0.0)),
            ],
            b: vec![(Some(self.ivar), Complex::new(self.acm, 0.0))],
        };
    }
}

#[derive(Default)]
pub struct Capacitor {
    c: f64,
    p: Option<VarIndex>,
    n: Option<VarIndex>,
    pp: Option<Eindex>,
    nn: Option<Eindex>,
    pn: Option<Eindex>,
    np: Option<Eindex>,
    op: CapOpPoint,
    guess: CapOpPoint,
}

#[derive(Clone, Copy, Default)]
struct CapOpPoint {
    v: f64,
    q: f64,
    i: f64,
}

impl Capacitor {
    pub fn new(c: f64, p: Option<VarIndex>, n: Option<VarIndex>) -> Capacitor {
        Capacitor {
            c,
            p,
            n,
            ..Default::default()
        }
    }
    fn q(&self, v: f64) -> f64 {
        return self.c * v;
    }
    fn dq_dv(&self, _v: f64) -> f64 {
        return self.c;
    }
}

impl Component for Capacitor {
    fn create_matrix_elems<T: SpNum>(&mut self, mat: &mut Matrix<T>) {
        self.pp = make_matrix_elem(mat, self.p, self.p);
        self.pn = make_matrix_elem(mat, self.p, self.n);
        self.np = make_matrix_elem(mat, self.n, self.p);
        self.nn = make_matrix_elem(mat, self.n, self.n);
    }
    /// Load our last guess as the new operating point
    fn commit(&mut self) {
        self.op = self.guess;
    }
    fn load(&mut self, guess: &Variables<f64>, an: &AnalysisInfo) -> Stamps<f64> {
        let vd = guess.get(self.p) - guess.get(self.n);
        let q = self.q(vd);

        match *an {
            AnalysisInfo::OP => {
                // FIXME: calculating this during DCOP, so we copy cleanly afterward
                // Should probably just find a way to calculate it then
                self.guess = CapOpPoint {
                    v: vd,
                    q: q,
                    i: 0.0,
                };
                return Stamps::new();
            }
            AnalysisInfo::TRAN(_, state) => {
                let (g, i, rhs) = state.integrate(q - self.op.q, self.dq_dv(vd), vd, self.op.i);
                self.guess = CapOpPoint { v: vd, q: q, i: i };

                return Stamps {
                    g: vec![(self.pp, g), (self.nn, g), (self.pn, -g), (self.np, -g)],
                    b: vec![(self.p, -rhs), (self.n, rhs)],
                };
            }
            AnalysisInfo::AC(_o, _s) => panic!("HOW WE GET HERE?!?"),
        }
    }
    fn load_ac(
        &mut self,
        _guess: &Variables<Complex<f64>>,
        an: &AnalysisInfo,
    ) -> Stamps<Complex<f64>> {
        let an_st = match an {
            AnalysisInfo::AC(_, state) => state,
            _ => panic!("Invalid AC AnalysisInfo"),
        };
        let c = self.dq_dv(0.0);
        return Stamps {
            g: vec![
                (self.pp, Complex::new(0.0, an_st.omega * c)),
                (self.nn, Complex::new(0.0, an_st.omega * c)),
                (self.pn, Complex::new(0.0, -an_st.omega * c)),
                (self.np, Complex::new(0.0, -an_st.omega * c)),
            ],
            b: vec![],
        };
    }
}

#[derive(Clone, Copy)]
enum TwoTerm {
    P = 0,
    N = 1,
}

struct TwoTerminals([Option<VarIndex>; 2]);

impl Index<TwoTerm> for TwoTerminals {
    type Output = Option<VarIndex>;
    fn index(&self, t: TwoTerm) -> &Option<VarIndex> {
        &self.0[t as usize]
    }
}

struct TwoTermMatrixPointers([[Option<Eindex>; 2]; 2]);

impl Index<(TwoTerm, TwoTerm)> for TwoTermMatrixPointers {
    type Output = Option<Eindex>;
    fn index(&self, ts: (TwoTerm, TwoTerm)) -> &Option<Eindex> {
        &self.0[ts.0 as usize][ts.1 as usize]
    }
}

impl IndexMut<(TwoTerm, TwoTerm)> for TwoTermMatrixPointers {
    fn index_mut(&mut self, ts: (TwoTerm, TwoTerm)) -> &mut Self::Output {
        &mut self.0[ts.0 as usize][ts.1 as usize]
    }
}

pub struct Resistor {
    g: f64,
    terms: TwoTerminals,
    matps: TwoTermMatrixPointers,
}

impl Resistor {
    pub fn new(g: f64, p: Option<VarIndex>, n: Option<VarIndex>) -> Resistor {
        Resistor {
            g,
            terms: TwoTerminals([p, n]),
            matps: TwoTermMatrixPointers([[None; 2]; 2]),
        }
    }
}

impl Component for Resistor {
    fn update(&mut self, val: f64) {
        self.g = val;
    }
    fn create_matrix_elems<T: SpNum>(&mut self, mat: &mut Matrix<T>) {
        use TwoTerm::{N, P};
        for l in [P, N].iter() {
            for r in [P, N].iter() {
                self.matps[(*l, *r)] = make_matrix_elem(mat, self.terms[*l], self.terms[*r]);
            }
        }
    }
    fn load(&mut self, _guess: &Variables<f64>, _an: &AnalysisInfo) -> Stamps<f64> {
        use TwoTerm::{N, P};
        return Stamps {
            g: vec![
                (self.matps[(P, P)], self.g),
                (self.matps[(N, N)], self.g),
                (self.matps[(P, N)], -self.g),
                (self.matps[(N, P)], -self.g),
            ],
            b: vec![],
        };
    }
    fn load_ac(
        &mut self,
        _guess: &Variables<Complex<f64>>,
        _an: &AnalysisInfo,
    ) -> Stamps<Complex<f64>> {
        use TwoTerm::{N, P};
        return Stamps {
            g: vec![
                (self.matps[(P, P)], Complex::new(self.g, 0.0)),
                (self.matps[(N, N)], Complex::new(self.g, 0.0)),
                (self.matps[(P, N)], Complex::new(-self.g, 0.0)),
                (self.matps[(N, P)], Complex::new(-self.g, 0.0)),
            ],
            b: vec![],
        };
    }
}

#[derive(Default)]
pub struct Isrc {
    i: f64,
    p: Option<VarIndex>,
    n: Option<VarIndex>,
}

impl Isrc {
    pub fn new(i: f64, p: Option<VarIndex>, n: Option<VarIndex>) -> Isrc {
        Isrc { i, p, n }
    }
}

impl Component for Isrc {
    fn create_matrix_elems<T: SpNum>(&mut self, _mat: &mut Matrix<T>) {}
    fn load(&mut self, _guess: &Variables<f64>, _an: &AnalysisInfo) -> Stamps<f64> {
        return Stamps {
            g: vec![],
            b: vec![(self.p, self.i), (self.n, -self.i)],
        };
    }
}

/// Helper function to create matrix element at (row,col) if both are non-ground
fn make_matrix_elem<T: SpNum>(
    mat: &mut Matrix<T>,
    row: Option<VarIndex>,
    col: Option<VarIndex>,
) -> Option<Eindex> {
    if let (Some(r), Some(c)) = (row, col) {
        return Some(mat.make(r.0, c.0));
    }
    return None;
}
