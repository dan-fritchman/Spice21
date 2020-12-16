use std::cmp::{max, min};
use std::fmt;
use std::ops::{Index, IndexMut};
use std::usize::MAX;

#[allow(unused_imports)] // These traits must be in scope, no matter what `cargo check` thinks. 
use num::{Num, One, Zero};

use crate::assert::assert;
use crate::{sperror, SpNum, SpResult};

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct Eindex(usize);

// `Entry`s are a type alias for tuples of (row, col, val).
type Entry<T> = (usize, usize, T);

#[derive(Debug, Copy, Clone)]
enum Axis {
    ROWS = 0,
    COLS,
}

use Axis::*;

impl Axis {
    fn other(&self) -> Axis {
        match self {
            Axis::ROWS => Axis::COLS,
            Axis::COLS => Axis::ROWS,
        }
    }
}

struct AxisPair<T> {
    rows: T,
    cols: T,
}

impl<T> Index<Axis> for AxisPair<T> {
    type Output = T;

    fn index(&self, ax: Axis) -> &Self::Output {
        match ax {
            Axis::ROWS => &self.rows,
            Axis::COLS => &self.cols,
        }
    }
}

impl<T> IndexMut<Axis> for AxisPair<T> {
    fn index_mut(&mut self, ax: Axis) -> &mut Self::Output {
        match ax {
            Axis::ROWS => &mut self.rows,
            Axis::COLS => &mut self.cols,
        }
    }
}

#[derive(PartialEq, Debug, Copy, Clone)]
enum MatrixState {
    CREATED = 0,
    FACTORING,
    FACTORED,
    RESET,
}

#[derive(Debug)]
pub struct Element<T: Num> {
    index: Eindex,
    row: usize,
    col: usize,
    val: T,
    fillin: bool,
    orig: (usize, usize, T),
    next_in_row: Option<Eindex>,
    next_in_col: Option<Eindex>,
}

impl<T: SpNum> PartialEq for Element<T> {
    fn eq(&self, other: &Self) -> bool {
        return self.row == other.row && self.col == other.col && self.val == other.val;
    }
}

impl<T: SpNum> Element<T> {
    fn new(index: Eindex, row: usize, col: usize, val: T, fillin: bool) -> Element<T> {
        Element {
            index,
            row,
            col,
            val,
            fillin,
            orig: (row, col, val),
            next_in_row: None,
            next_in_col: None,
        }
    }
    fn loc(&self, ax: Axis) -> usize {
        match ax {
            Axis::ROWS => self.row,
            Axis::COLS => self.col,
        }
    }
    fn set_loc(&mut self, ax: Axis, to: usize) {
        match ax {
            Axis::ROWS => self.row = to,
            Axis::COLS => self.col = to,
        }
    }
    fn next(&self, ax: Axis) -> Option<Eindex> {
        match ax {
            Axis::ROWS => self.next_in_row,
            Axis::COLS => self.next_in_col,
        }
    }
    fn set_next(&mut self, ax: Axis, e: Option<Eindex>) {
        match ax {
            Axis::ROWS => self.next_in_row = e,
            Axis::COLS => self.next_in_col = e,
        }
    }
}

#[derive(Debug)]
struct AxisMapping {
    e2i: Vec<usize>,
    i2e: Vec<usize>,
    history: Vec<(usize, usize)>,
}

///
/// AxisMapping from "internal", row/column-swapped, to "external" indices
///
/// In converting between externally-indexed vectors `x` and internal-indexed vectors `xi`,
/// each of the following holds:
///
/// xi[k] == x[i2e[k]]
/// xi[e2i[k]] == x[k]
///
/// For matrix-level conversions, operations on `x` solution vectors typically
/// access the *row* mappings.
/// Conversions of RHS vectors typically access *column* mappings.
///
impl AxisMapping {
    fn new(size: usize) -> AxisMapping {
        AxisMapping {
            e2i: (0..size).collect(),
            i2e: (0..size).collect(),
            history: vec![],
        }
    }
    fn swap_int(&mut self, x: usize, y: usize) {
        // Swap internal indices x and y
        let tmp = self.i2e[x];
        self.i2e[x] = self.i2e[y];
        self.i2e[y] = tmp;
        self.e2i[self.i2e[x]] = x;
        self.e2i[self.i2e[y]] = y;
        self.history.push((x, y));
    }
}

struct AxisData {
    hdrs: Vec<Option<Eindex>>,
    qtys: Vec<usize>,
    markowitz: Vec<usize>,
    mapping: Option<AxisMapping>,
}

impl AxisData {
    fn new() -> AxisData {
        AxisData {
            hdrs: vec![],
            qtys: vec![],
            markowitz: vec![],
            mapping: None,
        }
    }
    fn grow(&mut self, to: usize) {
        if to <= self.hdrs.len() {
            return;
        }
        let by = to - self.hdrs.len();
        for _ in 0..by {
            self.hdrs.push(None);
            self.qtys.push(0);
            self.markowitz.push(0);
        }
    }
    fn setup_factoring(&mut self) {
        self.markowitz.copy_from_slice(&self.qtys);
        if self.mapping.is_none() {
            self.mapping = Some(AxisMapping::new(self.hdrs.len()));
        }
    }
    fn swap(&mut self, x: usize, y: usize) {
        self.hdrs.swap(x, y);
        self.qtys.swap(x, y);
        self.markowitz.swap(x, y);
        if let Some(m) = &mut self.mapping {
            m.swap_int(x, y);
        }
    }
}

struct MarkowitzConfig {
    rel_threshold: f64,
    abs_threshold: f64,
    ties_mult: usize,
}

const MARKOWITZ_CONFIG: MarkowitzConfig = MarkowitzConfig {
    rel_threshold: 1e-3,
    abs_threshold: 0.0,
    ties_mult: 5,
};

/// Sparse Matrix
pub struct Matrix<T: Num> {
    // Matrix.elements is the owner of all `Element`s.
    // Everything else gets referenced via `Eindex`es.
    state: MatrixState,
    elements: Vec<Element<T>>,
    axes: AxisPair<AxisData>,
    diag: Vec<Option<Eindex>>,
    fillins: Vec<Eindex>,
}

impl<T: SpNum> Matrix<T> {
    /// Create a new, initially empty `Matrix`
    pub fn new() -> Matrix<T> {
        Matrix {
            state: MatrixState::CREATED,
            axes: AxisPair {
                rows: AxisData::new(),
                cols: AxisData::new(),
            },
            diag: vec![],
            elements: vec![],
            fillins: vec![],
        }
    }
    /// Create a new `Matrix` from a vector of (row, col, val) `entries`.
    pub fn from_entries(entries: Vec<Entry<T>>) -> Matrix<T> {
        let mut m = Matrix::<T>::new();
        for e in entries.iter() {
            m.add_element(e.0, e.1, e.2);
        }
        return m;
    }
    /// Add an element at location `(row, col)` with value `val`.
    pub fn add_element(&mut self, row: usize, col: usize, val: T) {
        self._add_element(row, col, val, false);
    }
    /// Add elements correspoding to each triplet `(row, col, val)`
    /// Rows and columns are `usize`, and `vals` are `T`.
    pub fn add_elements(&mut self, elements: Vec<Entry<T>>) {
        for e in elements.iter() {
            self.add_element(e.0, e.1, e.2);
        }
    }
    /// Create a zero-valued element at `(row, col)`,
    /// or return existing Element index if present
    pub fn make(&mut self, row: usize, col: usize) -> Eindex {
        return match self.get_elem(row, col) {
            Some(ei) => ei,
            None => self._add_element(row, col, T::zero(), false),
        };
    }
    /// Reset all Elements to zero value.
    pub fn reset(&mut self) {
        for e in self.elements.iter_mut() {
            e.val = T::zero();
        }
        self.state = MatrixState::RESET;
    }
    /// Update `Element` `ei` by `val`
    pub fn update(&mut self, ei: Eindex, val: T) {
        let tmp = self[ei].val + val;
        self[ei].val = tmp;
    }
    /// Multiply by Vec
    pub fn vecmul(&self, x: &Vec<T>) -> SpResult<Vec<T>> {
        if x.len() != self.num_cols() {
            return Err(sperror("Invalid Dimensions"));
        }
        let mut y: Vec<T> = vec![T::zero(); self.num_rows()];
        for row in 0..self.num_rows() {
            let mut ep = self.hdr(ROWS, row);
            while let Some(ei) = ep {
                y[row] += self[ei].val * x[self[ei].col];
                ep = self[ei].next_in_row;
            }
        }
        return Ok(y);
    }
    pub fn res(&self, x: &Vec<T>, rhs: &Vec<T>) -> SpResult<Vec<T>> {
        let mut xi: Vec<T> = vec![T::zero(); self.num_cols()];
        if let Some(col_mapping) = self.axes[COLS].mapping.as_ref() {
            // If we have factored, unwind any column-swaps
            for k in 0..xi.len() {
                xi[k] = x[col_mapping.i2e[k]];
            }
        } else {
            for k in 0..xi.len() {
                xi[k] = x[k];
            }
        }

        // Take the matrix-vector product with `xi`
        let ri: Vec<T> = self.vecmul(&xi)?;

        // Add in the RHS, unwinding row-swaps along the way
        let mut res = vec![T::zero(); ri.len()];
        if let Some(row_mapping) = self.axes[ROWS].mapping.as_ref() {
            for k in 0..xi.len() {
                res[k] = rhs[k] - ri[row_mapping.e2i[k]];
            }
        } else {
            for k in 0..xi.len() {
                res[k] = rhs[k] - ri[k];
            }
        }
        // println!("RES: {:?}", res);
        return Ok(res);
    }
    fn insert(&mut self, e: &mut Element<T>) {
        let mut expanded = false;
        if e.row + 1 > self.num_rows() {
            self.axes[Axis::ROWS].grow(e.row + 1);
            expanded = true;
        }
        if e.col + 1 > self.num_cols() {
            self.axes[Axis::COLS].grow(e.col + 1);
            expanded = true;
        }
        if expanded {
            let new_diag_len = std::cmp::min(self.num_rows(), self.num_cols());
            for _ in 0..new_diag_len - self.diag.len() {
                self.diag.push(None);
            }
        }
        // Insert along each Axis
        self.insert_axis(Axis::COLS, e);
        self.insert_axis(Axis::ROWS, e);
        // Update row & col qtys
        self.axes[ROWS].qtys[e.row] += 1;
        self.axes[COLS].qtys[e.col] += 1;
        if self.state == MatrixState::FACTORING {
            self.axes[ROWS].markowitz[e.row] += 1;
            self.axes[COLS].markowitz[e.col] += 1;
        }
        // Update our special arrays
        if e.row == e.col {
            self.diag[e.row] = Some(e.index);
        }
        if e.fillin {
            self.fillins.push(e.index);
        }
    }
    fn insert_axis(&mut self, ax: Axis, e: &mut Element<T>) {
        // Insert Element `e` along Axis `ax`

        let head_ptr = self.axes[ax].hdrs[e.loc(ax)];
        let head_idx = match head_ptr {
            Some(h) => h,
            None => {
                // Adding first element in this row/col
                return self.set_hdr(ax, e.loc(ax), Some(e.index));
            }
        };
        let off_ax = ax.other();
        if self[head_idx].loc(off_ax) > e.loc(off_ax) {
            // `e` is the new first element
            e.set_next(ax, head_ptr);
            return self.set_hdr(ax, e.loc(ax), Some(e.index));
        }

        // `e` comes after at least one Element.  Search for its position.
        let mut prev = head_idx;
        while let Some(next) = self[prev].next(ax) {
            if self[next].loc(off_ax) >= e.loc(off_ax) {
                break;
            }
            prev = next;
        }
        // And splice it in-between `prev` and `nxt`
        e.set_next(ax, self[prev].next(ax));
        self[prev].set_next(ax, Some(e.index));
    }
    fn add_fillin(&mut self, row: usize, col: usize) -> Eindex {
        return self._add_element(row, col, T::zero(), true);
    }
    fn _add_element(&mut self, row: usize, col: usize, val: T, fillin: bool) -> Eindex {
        // Element creation & insertion, used by `add_fillin` and the public `add_element`.
        let index = Eindex(self.elements.len());
        let mut e = Element::new(index.clone(), row, col, val, fillin);
        self.insert(&mut e);
        self.elements.push(e);
        return index;
    }
    /// Returns the Element-index at `(row, col)` if present, or None if not.
    pub fn get_elem(&self, row: usize, col: usize) -> Option<Eindex> {
        if row >= self.num_rows() {
            return None;
        }
        if col >= self.num_cols() {
            return None;
        }

        if row == col {
            // On diagonal; easy access
            return self.diag[row];
        }
        // Off-diagonal. Search across `row`.
        let mut ep = self.hdr(ROWS, row);
        while let Some(ei) = ep {
            let e = &self[ei];
            if e.col == col {
                return Some(ei);
            } else if e.col > col {
                return None;
            }
            ep = e.next_in_row;
        }
        return None;
    }
    /// Returns the Element-value at `(row, col)` if present, or None if not.
    pub fn get(&self, row: usize, col: usize) -> Option<T> {
        return match self.get_elem(row, col) {
            None => None,
            Some(ei) => Some(self[ei].val),
        };
    }
    fn move_element(&mut self, ax: Axis, idx: Eindex, to: usize) {
        let loc = self[idx].loc(ax);
        if loc == to {
            return;
        }
        let off_ax = ax.other();
        let y = self[idx].loc(off_ax);

        if loc < to {
            let br = match self.before_loc(off_ax, y, to, Some(idx)) {
                Some(ei) => ei,
                None => panic!("ERROR"),
            };
            if br != idx {
                let be = self.prev(off_ax, idx, None);
                let nxt = self[idx].next(off_ax);
                match be {
                    None => self.set_hdr(off_ax, y, nxt),
                    Some(be) => self[be].set_next(off_ax, nxt),
                };
                let brn = self[br].next(off_ax);
                self[idx].set_next(off_ax, brn);
                self[br].set_next(off_ax, Some(idx));
            }
        } else {
            let br = self.before_loc(off_ax, y, to, None);
            let be = self.prev(off_ax, idx, None);

            if br != be {
                // We (may) need some pointer updates
                if let Some(ei) = be {
                    let nxt = self[idx].next(off_ax);
                    self[ei].set_next(off_ax, nxt);
                }
                match br {
                    None => {
                        // New first in row/col
                        let first = self.hdr(off_ax, y);
                        self[idx].set_next(off_ax, first);
                        self.axes[off_ax].hdrs[y] = Some(idx);
                    }
                    Some(br) => {
                        if br != idx {
                            // Splice `idx` in after `br`
                            let nxt = self[br].next(off_ax);
                            self[idx].set_next(off_ax, nxt);
                            self[br].set_next(off_ax, Some(idx));
                        }
                    }
                };
            }
        }

        // Update the moved-Element's location
        self[idx].set_loc(ax, to);

        if loc == y {
            // If idx was on our diagonal, remove it
            self.diag[loc] = None;
        } else if to == y {
            // Or if it's now on the diagonal, add it
            self.diag[to] = Some(idx);
        }
    }
    fn exchange_elements(&mut self, ax: Axis, ix: Eindex, iy: Eindex) {
        // Swap two elements `ax` indices.
        // Elements must be in the same off-axis vector,
        // and the first argument `ex` must be the lower-indexed off-axis.
        // E.g. exchange_elements(Axis.rows, ex, ey) exchanges the rows of ex and ey.

        let off_ax = ax.other();
        let off_loc = self[ix].loc(off_ax);

        let bx = self.prev(off_ax, ix, None);
        let by = match self.prev(off_ax, iy, Some(ix)) {
            Some(e) => e,
            None => panic!("ERROR!"),
        };

        let locx = self[ix].loc(ax);
        let locy = self[iy].loc(ax);
        self[iy].set_loc(ax, locx);
        self[ix].set_loc(ax, locy);

        match bx {
            None => {
                // If `ex` is the *first* entry in the column, replace it to our header-list
                self.set_hdr(off_ax, off_loc, Some(iy));
            }
            Some(bxe) => {
                // Otherwise patch ey into bx
                self[bxe].set_next(off_ax, Some(iy));
            }
        }

        if by == ix {
            // `ex` and `ey` are adjacent
            let tmp = self[iy].next(off_ax);
            self[iy].set_next(off_ax, Some(ix));
            self[ix].set_next(off_ax, tmp);
        } else {
            // Elements in-between `ex` and `ey`.  Update the last one.
            let xnxt = self[ix].next(off_ax);
            let ynxt = self[iy].next(off_ax);
            self[iy].set_next(off_ax, xnxt);
            self[ix].set_next(off_ax, ynxt);
            self[by].set_next(off_ax, Some(ix));
        }

        // Update our diagonal array, if necessary
        if locx == off_loc {
            self.diag[off_loc] = Some(iy);
        } else if locy == off_loc {
            self.diag[off_loc] = Some(ix);
        }
    }
    fn prev(&self, ax: Axis, idx: Eindex, hint: Option<Eindex>) -> Option<Eindex> {
        // Find the element previous to `idx` along axis `ax`.
        // If provided, `hint` *must* be before `idx`, or search will fail.
        let prev: Option<Eindex> = match hint {
            Some(_) => hint,
            None => self.hdr(ax, self[idx].loc(ax)),
        };
        let mut pi: Eindex = match prev {
            None => {
                return None;
            }
            Some(pi) if pi == idx => {
                return None;
            }
            Some(pi) => pi,
        };
        while let Some(nxt) = self[pi].next(ax) {
            if nxt == idx {
                break;
            }
            pi = nxt;
        }
        return Some(pi);
    }
    fn before_loc(&self, ax: Axis, loc: usize, before: usize, hint: Option<Eindex>) -> Option<Eindex> {
        let prev: Option<Eindex> = match hint {
            Some(_) => hint,
            None => self.hdr(ax, loc),
        };
        let off_ax = ax.other();
        let mut pi: Eindex = match prev {
            None => {
                return None;
            }
            Some(pi) if self[pi].loc(off_ax) >= before => {
                return None;
            }
            Some(pi) => pi,
        };
        while let Some(nxt) = self[pi].next(ax) {
            if self[nxt].loc(off_ax) >= before {
                break;
            }
            pi = nxt;
        }
        return Some(pi);
    }
    fn swap(&mut self, ax: Axis, a: usize, b: usize) {
        if a == b {
            return;
        }
        let x = min(a, b);
        let y = max(a, b);

        let hdrs = &self.axes[ax].hdrs;
        let mut ix = hdrs[x];
        let mut iy = hdrs[y];
        let off_ax = ax.other();

        loop {
            match (ix, iy) {
                (Some(ex), Some(ey)) => {
                    let ox = self[ex].loc(off_ax);
                    let oy = self[ey].loc(off_ax);
                    if ox < oy {
                        self.move_element(ax, ex, y);
                        ix = self[ex].next(ax);
                    } else if oy < ox {
                        self.move_element(ax, ey, x);
                        iy = self[ey].next(ax);
                    } else {
                        self.exchange_elements(ax, ex, ey);
                        ix = self[ex].next(ax);
                        iy = self[ey].next(ax);
                    }
                }
                (None, Some(ey)) => {
                    self.move_element(ax, ey, x);
                    iy = self[ey].next(ax);
                }
                (Some(ex), None) => {
                    self.move_element(ax, ex, y);
                    ix = self[ex].next(ax);
                }
                (None, None) => {
                    break;
                }
            }
        }
        // Swap all the relevant pointers & counters
        self.axes[ax].swap(x, y);
    }
    /// Updates self to S = L + U - I.
    /// Diagonal entries are those of U;
    /// L has diagonal entries equal to one.
    fn lu_factorize(&mut self) -> SpResult<()> {
        assert(self.diag.len()).gt(0)?;
        for k in 0..self.axes[ROWS].hdrs.len() {
            if self.hdr(ROWS, k).is_none() {
                return Err(sperror("Singular Matrix"));
            }
        }
        for k in 0..self.axes[COLS].hdrs.len() {
            if self.hdr(COLS, k).is_none() {
                return Err(sperror("Singular Matrix"));
            }
        }
        self.state = MatrixState::FACTORING;
        self.axes[ROWS].setup_factoring();
        self.axes[COLS].setup_factoring();

        for n in 0..self.diag.len() - 1 {
            let pivot = match self.search_for_pivot(n) {
                None => return Err(sperror("Pivot Search Fail")),
                Some(p) => p,
            };
            self.swap(ROWS, self[pivot].row, n);
            self.swap(COLS, self[pivot].col, n);
            self.row_col_elim(pivot, n)?;
        }
        self.state = MatrixState::FACTORED;
        return Ok(());
    }

    fn search_for_pivot(&self, n: usize) -> Option<Eindex> {
        let mut ei = self.markowitz_search_diagonal(n);
        if let Some(_) = ei {
            return ei;
        }
        ei = self.markowitz_search_submatrix(n);
        if let Some(_) = ei {
            return ei;
        }
        return self.find_max(n);
    }

    fn max_after(&self, ax: Axis, after: Eindex) -> Eindex {
        let mut best = after;
        let mut best_val = self[after].val.absv();
        let mut e = self[after].next(ax);

        while let Some(ei) = e {
            let val = self[ei].val.absv();
            if val > best_val {
                best = ei;
                best_val = val;
            }
            e = self[ei].next(ax);
        }
        return best;
    }
    /// Find the max (abs) value element at/after location `after_loc` in row/col `in_loc`
    fn max_after_loc(&self, ax: Axis, in_loc: usize, after_loc: usize) -> Option<Eindex> {
        let mut e = self.axes[ax].hdrs[in_loc];
        let off_ax = ax.other();
        while let Some(ei) = e {
            if self[ei].loc(off_ax) >= after_loc {
                break;
            }
            e = self[ei].next(ax);
        }
        let mut best = e?;
        let mut best_val = self[best].val.absv();
        while let Some(ei) = e {
            let val = self[ei].val.absv();
            if val > best_val {
                best = ei;
                best_val = val;
            }
            e = self[ei].next(ax);
        }
        return Some(best);
    }

    fn markowitz_product(&self, ei: Eindex) -> usize {
        let e = &self[ei];
        let mr = self[Axis::ROWS].markowitz[e.row];
        let mc = self[Axis::COLS].markowitz[e.col];
        assert!(mr > 0); // FIXME: figure out how to pipe these around
        assert!(mc > 0);
        return (mr - 1) * (mc - 1);
    }

    fn markowitz_search_diagonal(&self, n: usize) -> Option<Eindex> {
        let mut best_elem = None;
        let mut best_mark = MAX; // Actually use usize::MAX!
        let mut best_ratio = 0.0;
        let mut num_ties = 0;

        for k in n..self.diag.len() {
            let d = match self.diag[k] {
                None => {
                    continue;
                }
                Some(d) => d,
            };

            // Check whether this element meets our threshold criteria
            let max_in_col = match self.max_after_loc(COLS, k, n) {
                None => {
                    continue;
                }
                Some(d) => d,
            };
            // FIXME: abs-value criterion for max_in_col
            let threshold = MARKOWITZ_CONFIG.rel_threshold * self[max_in_col].val.absv() + MARKOWITZ_CONFIG.abs_threshold;
            if self[d].val.absv() < threshold {
                continue;
            }

            // If so, compute and compare its Markowitz product to our best
            let mark = self.markowitz_product(d);
            if mark < best_mark {
                num_ties = 0;
                best_elem = self.diag[k];
                best_mark = mark;
                best_ratio = (self[d].val / self[max_in_col].val).absv();
            } else if mark == best_mark {
                num_ties += 1;
                let ratio = (self[d].val / self[max_in_col].val).absv();
                if ratio > best_ratio {
                    best_elem = self.diag[k];
                    best_mark = mark;
                    best_ratio = ratio;
                }
                if num_ties >= best_mark * MARKOWITZ_CONFIG.ties_mult {
                    return best_elem;
                }
            }
        }
        return best_elem;
    }

    fn markowitz_search_submatrix(&self, n: usize) -> Option<Eindex> {
        let mut best_elem = None;
        let mut best_mark = MAX; // Actually use usize::MAX!
        let mut best_ratio = 0.0;
        //        let mut num_ties = 0;

        for _k in n..self.axes[COLS].hdrs.len() {
            let mut e = self.hdr(COLS, n);
            // Advance to a row ≥ n
            while let Some(ei) = e {
                if self[ei].row >= n {
                    break;
                }
                e = self[ei].next_in_col;
            }
            let ei = match e {
                None => {
                    continue;
                }
                Some(d) => d,
            };

            // Check whether this element meets our threshold criteria
            let max_in_col = self.max_after(COLS, ei);
            //            let threshold = MARKOWITZ_CONFIG.rel_threshol * self[max_in_col].val.absv() + MARKOWITZ_CONFIG.abs_threshol;

            while let Some(ei) = e {
                // If so, compute and compare its Markowitz product to our best
                let mark = self.markowitz_product(ei);
                if mark < best_mark {
                    //                    num_ties = 0;
                    best_elem = e;
                    best_mark = mark;
                    best_ratio = (self[ei].val / self[max_in_col].val).absv();
                } else if mark == best_mark {
                    //                    num_ties += 1;
                    let ratio = (self[ei].val / self[max_in_col].val).absv();
                    if ratio > best_ratio {
                        best_elem = e;
                        best_mark = mark;
                        best_ratio = ratio;
                    }
                    //                    // FIXME: do we want tie-counting in here?
                    //                    if num_ties >= best_mark * MARKOWITZ_CONFIG.ties_mult { return best_elem; }
                }
                e = self[ei].next_in_col;
            }
        }
        return best_elem;
    }
    /// Find the max (abs value) element in sub-matrix of indices ≥ `n`.
    /// Returns `None` if no elements present.
    fn find_max(&self, n: usize) -> Option<Eindex> {
        let mut max_elem = None;
        let mut max_val = 0.0;

        // Search each column ≥ n
        for k in n..self.axes[COLS].hdrs.len() {
            let mut ep = self.hdr(COLS, k);

            // Advance to a row ≥ n
            while let Some(ei) = ep {
                if self[ei].row >= n {
                    break;
                }
                ep = self[ei].next_in_col;
            }
            // And search over remaining elements
            while let Some(ei) = ep {
                let val = self[ei].val.absv();
                if val > max_val {
                    max_elem = ep;
                    max_val = val;
                }
                ep = self[ei].next_in_col;
            }
        }
        return max_elem;
    }

    fn row_col_elim(&mut self, pivot: Eindex, n: usize) -> SpResult<()> {
        let de = match self.diag[n] {
            Some(de) => de,
            None => return Err(sperror("Singular Matrix")),
        };
        assert(de).eq(pivot)?;
        let pivot_val = self[pivot].val;
        assert(pivot_val).ne(T::zero())?;

        // Divide elements in the pivot column by the pivot-value
        let mut plower = self[pivot].next_in_col;
        while let Some(ple) = plower {
            self[ple].val /= pivot_val;
            plower = self[ple].next_in_col;
        }

        let mut pupper = self[pivot].next_in_row;
        while let Some(pue) = pupper {
            let pupper_col = self[pue].col;
            plower = self[pivot].next_in_col;
            let mut psub = self[pue].next_in_col;
            while let Some(ple) = plower {
                // Walk `psub` down to the lower pointer
                while let Some(pse) = psub {
                    if self[pse].row >= self[ple].row {
                        break;
                    }
                    psub = self[pse].next_in_col;
                }
                let pse = match psub {
                    None => self.add_fillin(self[ple].row, pupper_col),
                    Some(pse) if self[pse].row > self[ple].row => self.add_fillin(self[ple].row, pupper_col),
                    Some(pse) => pse,
                };

                // Update the `psub` element value
                let v: T = self[pue].val.clone() * self[ple].val.clone();
                self[pse].val -= v;
                psub = self[pse].next_in_col;
                plower = self[ple].next_in_col;
            }
            self.axes[COLS].markowitz[pupper_col] -= 1;
            pupper = self[pue].next_in_row;
        }
        // Update remaining Markowitz counts
        self.axes[ROWS].markowitz[n] -= 1;
        self.axes[COLS].markowitz[n] -= 1;
        plower = self[pivot].next_in_col;
        while let Some(ple) = plower {
            let plower_row = self[ple].row;
            self.axes[ROWS].markowitz[plower_row] -= 1;
            plower = self[ple].next_in_col;
        }
        return Ok(());
    }
    /// Solve the system `Ax=b`, where:
    /// * `A` is `self`
    /// * `b` is argument `rhs`
    /// * `x` is the return value.
    ///
    /// Returns a `Result` containing the `Vec<T>` representing `x` if successful.
    /// Returns an `Err` if unsuccessful.
    ///
    /// Performs LU factorization, forward and backward substitution.
    pub fn solve(&mut self, rhs: Vec<T>) -> SpResult<Vec<T>> {
        if self.state != MatrixState::FACTORED {
            self.lu_factorize()?;
        }
        assert(self.state).eq(MatrixState::FACTORED)?;

        // Unwind any row-swaps
        let mut c: Vec<T> = vec![T::zero(); rhs.len()];

        if let Some(row_mapping) = self.axes[ROWS].mapping.as_ref() {
            for k in 0..c.len() {
                c[k] = rhs[row_mapping.i2e[k]];
            }
        } else {
            return Err(sperror("Missing Row Mapping"));
        }

        // Forward substitution: Lc=b
        for k in 0..self.diag.len() {
            // Walk down each column, update c
            if c[k] == T::zero() {
                continue;
            } // No updates to make on this iteration

            // c[d.row] /= d.val

            let di = match self.diag[k] {
                Some(di) => di,
                None => return Err(sperror("Singular Matrix")),
            };
            let mut e = self[di].next_in_col;
            while let Some(ei) = e {
                c[self[ei].row] = c[self[ei].row] - c[k] * self[ei].val; // FIXME: SubAssign
                e = self[ei].next_in_col;
            }
        }

        // Backward substitution: Ux=c
        for k in (0..self.diag.len()).rev() {
            // Walk each row, update c
            let di = match self.diag[k] {
                Some(di) => di,
                None => return Err(sperror("Singular Matrix")),
            };
            let mut ep = self[di].next_in_row;
            while let Some(ei) = ep {
                c[k] = c[k] - c[self[ei].col] * self[ei].val; // FIXME: SubAssign
                ep = self[ei].next_in_row;
            }
            c[k] /= self[di].val;
        }

        // Unwind any column-swaps
        let mut soln: Vec<T> = vec![T::zero(); c.len()];
        if let Some(col_mapping) = self.axes[COLS].mapping.as_ref() {
            for k in 0..c.len() {
                soln[k] = c[col_mapping.e2i[k]];
            }
        } else {
            return Err(sperror("Missing Column Mapping"));
        }
        return Ok(soln);
    }
    /// Create a row-majory dense matrix representation
    pub fn to_dense(&self) -> Vec<Vec<T>> {
        let mut res = vec![vec![T::zero(); self.num_cols()]; self.num_rows()];
        for ei in self.elements.iter() {
            res[self[ei.index].row][self[ei.index].col] = self[ei.index].val;
        }
        return res;
    }
    fn hdr(&self, ax: Axis, loc: usize) -> Option<Eindex> {
        self.axes[ax].hdrs[loc]
    }
    fn set_hdr(&mut self, ax: Axis, loc: usize, ei: Option<Eindex>) {
        self.axes[ax].hdrs[loc] = ei;
    }
    fn num_rows(&self) -> usize {
        self.axes[ROWS].hdrs.len()
    }
    fn num_cols(&self) -> usize {
        self.axes[COLS].hdrs.len()
    }
}

impl<T: SpNum + One> Matrix<T> {
    /// Create an n*n identity `Matrix`
    pub fn identity(n: usize) -> Matrix<T> {
        let mut m = Matrix::<T>::new();
        for k in 0..n {
            m.add_element(k, k, T::one());
        }
        return m;
    }
}

impl<T: SpNum> Index<Eindex> for Matrix<T> {
    type Output = Element<T>;
    fn index(&self, index: Eindex) -> &Self::Output {
        &self.elements[index.0]
    }
}

impl<T: SpNum> IndexMut<Eindex> for Matrix<T> {
    fn index_mut(&mut self, index: Eindex) -> &mut Self::Output {
        &mut self.elements[index.0]
    }
}

impl<T: SpNum> Index<Axis> for Matrix<T> {
    type Output = AxisData;
    fn index(&self, ax: Axis) -> &Self::Output {
        &self.axes[ax]
    }
}

impl<T: SpNum> IndexMut<Axis> for Matrix<T> {
    fn index_mut(&mut self, ax: Axis) -> &mut Self::Output {
        &mut self.axes[ax]
    }
}

impl<T: SpNum> fmt::Debug for Matrix<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "sparse21::Matrix (rows={}, cols={}, elems={}\n",
            self.num_rows(),
            self.num_cols(),
            self.elements.len()
        )?;
        let d = self.to_dense();
        for row in d.iter() {
            write!(f, "{:?}\n", *row)?;
        }
        for e in self.elements.iter() {
            write!(f, "({}, {}, {}) \n", e.row, e.col, e.val)?;
        }
        write!(f, "\n")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spresult::TestResult;
    use num::Complex;

    impl<T: SpNum> Matrix<T> {
        pub fn swap_rows(&mut self, x: usize, y: usize) {
            self.swap(ROWS, x, y)
        }
        pub fn swap_cols(&mut self, x: usize, y: usize) {
            self.swap(COLS, x, y)
        }
        pub fn size(&self) -> (usize, usize) {
            (self.num_rows(), self.num_cols())
        }

        pub fn checkups(&self) -> TestResult {
            // Internal consistency tests.  Probably pretty slow.

            self.check_diagonal()?;

            let mut next_in_rows: Vec<Eindex> = vec![];
            let mut next_in_cols: Vec<Eindex> = vec![];

            for n in 0..self.axes[COLS].hdrs.len() {
                let mut ep = self.hdr(COLS, n);
                while let Some(ei) = ep {
                    assert(self[ei].col).eq(n)?;

                    if let Some(nxt) = self[ei].next_in_col {
                        assert(self[nxt].row).gt(self[ei].row)?;
                        assert!(!next_in_cols.contains(&nxt));
                        next_in_cols.push(nxt);
                    }
                    if let Some(nxt) = self[ei].next_in_row {
                        assert(self[nxt].col).gt(self[ei].col)?;
                        assert!(!next_in_rows.contains(&nxt));
                        next_in_rows.push(nxt);
                    }
                    ep = self[ei].next_in_col;
                }
            }
            // Add the row/column headers to the "next" vectors
            for ep in self.axes[Axis::COLS].hdrs.iter() {
                if let Some(ei) = ep {
                    assert!(!next_in_cols.contains(ei));
                    next_in_cols.push(*ei);
                }
            }
            for ep in self.axes[Axis::ROWS].hdrs.iter() {
                if let Some(ei) = ep {
                    assert!(!next_in_rows.contains(ei));
                    next_in_rows.push(*ei);
                }
            }
            // Check that all elements are included
            assert(next_in_cols.len()).eq(self.elements.len())?;
            assert(next_in_rows.len()).eq(self.elements.len())?;
            for n in 0..self.elements.len() {
                assert!(next_in_cols.contains(&Eindex(n)));
                assert!(next_in_rows.contains(&Eindex(n)));
            }
            Ok(())
        }

        fn check_diagonal(&self) -> TestResult {
            for r in 0..self.diag.len() {
                let eo = self.get(r, r);
                if let Some(e) = eo {
                    if let Some(d) = self.diag[r] {
                        assert(self[d].val).eq(e)?;
                    } else {
                        return Err(sperror("FAIL!"));
                    }
                // FIXME: would prefer something like the previous "same element ID" testing
                // assert_eq!(e, m[self.diag[r]].val);
                //                    assert_eq!(e.index, self.diag[r]);
                //                    assert_eq!(e.row, r);
                //                    assert_eq!(e.col, r);
                } else {
                    assert(self.diag[r]).eq(None)?;
                }
            }
            Ok(())
        }
        fn assert_entries(&self, entries: Vec<Entry<T>>) -> TestResult {
            for e in entries.into_iter() {
                assert(self.get(e.0, e.1).ok_or("ElementMissing")?).eq(e.2)?;
            }
            Ok(())
        }
    }

    #[test]
    fn test_create_element() -> TestResult {
        let e = Element::new(Eindex(0), 0, 0, 1.0, false);
        assert_eq!(e.index.0, 0);
        assert_eq!(e.row, 0);
        assert_eq!(e.col, 0);
        assert_eq!(e.val, 1.0);
        assert_eq!(e.fillin, false);
        assert_eq!(e.next_in_row, None);
        assert_eq!(e.next_in_col, None);
        Ok(())
    }

    #[test]
    fn test_create_matrix() -> TestResult {
        let m = Matrix::<f64>::new();
        assert_eq!(m.state, MatrixState::CREATED);
        assert_eq!(m.diag, vec![]);
        assert_eq!(m.axes[Axis::ROWS].hdrs, vec![]);
        Ok(())
    }

    #[test]
    fn test_add_element() -> TestResult {
        let mut m = Matrix::new();

        m.add_element(0, 0, 1.0);
        assert_eq!(m.num_rows(), 1);
        assert_eq!(m.num_cols(), 1);
        assert_eq!(m.size(), (1, 1));
        assert_eq!(m.diag.len(), 1);

        m.add_element(100, 100, 1.0);
        assert_eq!(m.num_rows(), 101);
        assert_eq!(m.num_cols(), 101);
        assert_eq!(m.size(), (101, 101));
        assert_eq!(m.diag.len(), 101);
        Ok(())
    }

    #[test]
    fn test_get() -> TestResult {
        let mut m = Matrix::new();
        m.add_element(0, 0, 1.0);
        assert(m.get(0, 0).ok_or("ElementMissing")?).eq(1.0)?;
        Ok(())
    }

    #[test]
    fn test_identity() -> TestResult {
        // Check identity matrices of each (small) size
        for k in 1..10 {
            let ik = Matrix::<f64>::identity(k);

            // Basic size checks
            assert_eq!(ik.num_rows(), k);
            assert_eq!(ik.num_cols(), k);
            assert_eq!(ik.size(), (k, k));
            assert_eq!(ik.elements.len(), k);
            ik.checkups()?;

            for v in 0..k {
                // Check each row/ col head is the same element, and this element is on the diagonal
                let ro = ik.hdr(Axis::ROWS, v).unwrap();
                let co = ik.hdr(Axis::COLS, v).unwrap();
                let d0 = ik.get(v, v).ok_or("ElementMissing")?;
                assert_eq!(ro, co);
                assert_eq!(ik[ro].val, d0);
                assert_eq!(ik[co].val, d0);
            }
        }
        Ok(())
    }

    #[test]
    fn test_swap_rows0() -> TestResult {
        let mut m = Matrix::new();

        m.add_element(0, 0, 11.0);
        m.add_element(7, 0, 22.0);
        m.add_element(0, 7, 33.0);
        m.add_element(7, 7, 44.0);

        m.checkups()?;
        assert_eq!(m.get(0, 0).ok_or("ElementMissing")?, 11.0);
        assert_eq!(m.get(7, 0).ok_or("ElementMissing")?, 22.0);
        assert_eq!(m.get(0, 7).ok_or("ElementMissing")?, 33.0);
        assert_eq!(m.get(7, 7).ok_or("ElementMissing")?, 44.0);

        m.state = MatrixState::FACTORING;
        m.swap_rows(0, 7);

        m.checkups()?;
        assert_eq!(m.get(7, 0).ok_or("ElementMissing")?, 11.0);
        assert_eq!(m.get(0, 0).ok_or("ElementMissing")?, 22.0);
        assert_eq!(m.get(7, 7).ok_or("ElementMissing")?, 33.0);
        assert_eq!(m.get(0, 7).ok_or("ElementMissing")?, 44.0);
        Ok(())
    }

    #[test]
    fn test_swap_rows1() -> TestResult {
        let mut m = Matrix::new();

        m.add_element(0, 0, 11.1);
        m.add_element(2, 2, 22.2);

        m.checkups()?;
        assert_eq!(m.get(0, 0).ok_or("ElementMissing")?, 11.1);
        assert_eq!(m.get(2, 2).ok_or("ElementMissing")?, 22.2);
        assert_eq!(m.get(1, 1), None);

        m.state = MatrixState::FACTORING;
        m.swap_rows(0, 2);

        m.checkups()?;
        assert_eq!(m.get(2, 0).ok_or("ElementMissing")?, 11.1);
        assert_eq!(m.get(0, 2).ok_or("ElementMissing")?, 22.2);
        assert_eq!(m.get(1, 1), None);
        Ok(())
    }

    #[test]
    fn test_swap_rows2() -> TestResult {
        let mut m = Matrix::new();

        m.add_element(0, 0, 1.0);
        m.add_element(0, 1, 2.0);
        m.add_element(0, 2, 3.0);
        m.add_element(1, 0, 4.0);
        m.add_element(1, 1, 5.0);
        m.add_element(1, 2, 6.0);
        m.add_element(2, 0, 7.0);
        m.add_element(2, 1, 8.0);
        m.add_element(2, 2, 9.0);

        m.checkups()?;
        m.state = MatrixState::FACTORING;
        m.swap_rows(0, 2);

        m.checkups()?;
        assert_eq!(m.get(0, 0).ok_or("ElementMissing")?, 7.0);
        assert_eq!(m.get(2, 0).ok_or("ElementMissing")?, 1.0);
        // FIXME: check more
        Ok(())
    }

    #[test]
    fn test_swap_rows3() -> TestResult {
        let mut m = Matrix::new();
        m.add_element(1, 0, 71.0);
        m.add_element(2, 0, -11.0);
        m.add_element(2, 2, 99.0);

        m.checkups()?;
        assert_eq!(m.get(1, 0).ok_or("ElementMissing")?, 71.0);
        assert_eq!(m.get(2, 0).ok_or("ElementMissing")?, -11.0);
        assert_eq!(m.get(2, 2).ok_or("ElementMissing")?, 99.0);

        m.state = MatrixState::FACTORING;
        m.swap_rows(0, 2);

        m.checkups()?;
        assert_eq!(m.get(1, 0).ok_or("ElementMissing")?, 71.0);
        assert_eq!(m.get(0, 0).ok_or("ElementMissing")?, -11.0);
        assert_eq!(m.get(0, 2).ok_or("ElementMissing")?, 99.0);
        Ok(())
    }

    #[test]
    fn test_swap_rows4() -> TestResult {
        let mut m = Matrix::new();

        for r in 0..3 {
            for c in 0..3 {
                if r != 0 || c != 1 {
                    m.add_element(r, c, ((r + 1) * (c + 1)) as f64);
                }
            }
        }
        m.checkups()?;

        m.state = MatrixState::FACTORING;
        m.swap_rows(0, 1);

        m.checkups()?;

        // FIXME: add some real checks on this
        Ok(())
    }

    #[test]
    fn test_row_mappings() -> TestResult {
        let mut m = Matrix::<f64>::identity(4);
        m.checkups()?;

        m.state = MatrixState::FACTORING;
        m.axes[ROWS].setup_factoring();
        m.swap_rows(0, 3);

        m.checkups()?;
        assert_eq!(m.axes[Axis::ROWS].mapping.as_ref().unwrap().e2i, vec![3, 1, 2, 0]);
        assert_eq!(m.axes[Axis::ROWS].mapping.as_ref().unwrap().i2e, vec![3, 1, 2, 0]);

        m.swap_rows(0, 2);

        m.checkups()?;
        assert_eq!(m.axes[Axis::ROWS].mapping.as_ref().unwrap().e2i, vec![3, 1, 0, 2]);
        assert_eq!(m.axes[Axis::ROWS].mapping.as_ref().unwrap().i2e, vec![2, 1, 3, 0]);
        Ok(())
    }

    #[test]
    fn test_lu_id3() -> TestResult {
        let mut m = Matrix::<f64>::identity(3);
        m.checkups()?;
        m.lu_factorize()?;
        m.checkups()?;
        assert_eq!(m.get(0, 0).ok_or("ElementMissing")?, 1.0);
        assert_eq!(m.get(1, 1).ok_or("ElementMissing")?, 1.0);
        assert_eq!(m.get(2, 2).ok_or("ElementMissing")?, 1.0);
        return Ok(());
    }

    #[test]
    fn test_lu_lower() -> TestResult {
        // Factors a unit lower-diagonal matrix.  Should leave it unchanged.

        let mut m = Matrix::new();
        m.add_element(0, 0, 1.0);
        m.add_element(1, 0, 1.0);
        m.add_element(2, 0, 1.0);
        m.add_element(1, 1, 1.0);
        m.add_element(2, 1, 1.0);
        m.add_element(2, 2, 1.0);

        m.checkups()?;
        assert_eq!(m.get(0, 0).ok_or("ElementMissing")?, 1.0);
        assert_eq!(m.get(1, 0).ok_or("ElementMissing")?, 1.0);
        assert_eq!(m.get(2, 0).ok_or("ElementMissing")?, 1.0);
        assert_eq!(m.get(1, 1).ok_or("ElementMissing")?, 1.0);
        assert_eq!(m.get(2, 1).ok_or("ElementMissing")?, 1.0);
        assert_eq!(m.get(2, 2).ok_or("ElementMissing")?, 1.0);

        m.lu_factorize()?;

        m.checkups()?;
        assert_eq!(m.get(0, 0).ok_or("ElementMissing")?, 1.0);
        assert_eq!(m.get(1, 0).ok_or("ElementMissing")?, 1.0);
        assert_eq!(m.get(2, 0).ok_or("ElementMissing")?, 1.0);
        assert_eq!(m.get(1, 1).ok_or("ElementMissing")?, 1.0);
        assert_eq!(m.get(2, 1).ok_or("ElementMissing")?, 1.0);
        assert_eq!(m.get(2, 2).ok_or("ElementMissing")?, 1.0);
        Ok(())
    }

    #[test]
    fn test_lu() -> TestResult {
        let mut m = Matrix::from_entries(vec![
            (2, 2, -1.0),
            (2, 1, 5.0),
            (2, 0, 2.0),
            (1, 2, 5.0),
            (1, 1, 2.0),
            (0, 2, 1.0),
            (0, 1, 1.0),
            (0, 0, 1.0),
        ]);
        m.checkups()?;
        m.assert_entries(vec![
            (2, 2, -1.0),
            (2, 1, 5.0),
            (2, 0, 2.0),
            (1, 2, 5.0),
            (1, 1, 2.0),
            (0, 2, 1.0),
            (0, 1, 1.0),
            (0, 0, 1.0),
        ])?;
        m.lu_factorize()?;
        m.checkups()?;
        Ok(())
    }

    #[test]
    fn test_solve() -> TestResult {
        let mut m = Matrix::from_entries(vec![
            (0, 0, 1.0),
            (0, 1, 1.0),
            (0, 2, 1.0),
            (1, 1, 2.0),
            (1, 2, 5.0),
            (2, 0, 2.0),
            (2, 1, 5.0),
            (2, 2, -1.0),
        ]);
        m.checkups()?;
        m.lu_factorize()?;
        m.checkups()?;
        let rhs = vec![6.0, -4.0, 27.0];
        let soln = m.solve(rhs)?;
        let correct = vec![5.0, 3.0, -2.0];
        for k in 0..soln.len() {
            assert!(isclose(soln[k], correct[k]));
        }
        Ok(())
    }

    #[test]
    fn test_solve_id3() -> TestResult {
        let mut m = Matrix::<f64>::identity(3);
        let soln = m.solve(vec![11.1, 30.3, 99.9])?;
        assert_eq!(soln, vec![11.1, 30.3, 99.9]);
        return Ok(());
    }

    #[test]
    fn test_solve_identity() -> TestResult {
        // Test that solutions of Ix=b yield x=b
        for s in 1..10 {
            let mut m = Matrix::<f64>::identity(s);
            let mut rhs: Vec<f64> = vec![];
            for e in 0..s {
                rhs.push(e as f64);
            }

            let soln = m.solve(rhs.clone())?;
            assert_eq!(soln, rhs);
        }
        return Ok(());
    }

    fn isclose(a: f64, b: f64) -> bool {
        return (a - b).abs() < 1e-9;
    }

    #[test]
    fn test_create_complex() -> TestResult {
        let mut m = Matrix::from_entries(vec![(0, 0, Complex::one()), (1, 1, Complex::one())]);
        m.lu_factorize()?;
        Ok(())
    }
    #[test]
    fn test_solve_complex_id2() -> TestResult {
        let mut m = Matrix::from_entries(vec![(0, 0, Complex::one()), (1, 1, Complex::one())]);
        let soln = m.solve(vec![Complex::i(); 2])?;
        assert(soln).eq(vec![Complex::i(); 2])?;
        Ok(())
    }

    #[test]
    fn test_solve_complex() -> TestResult {
        let mut m = Matrix::from_entries(vec![
            (0, 0, Complex::one()),
            (1, 0, Complex::new(-1.0, 0.0)),
            (0, 1, Complex::new(-1.0, 0.0)),
            (1, 1, Complex::new(1.0, 1.0)),
        ]);
        let soln = m.solve(vec![Complex::one(), Complex::zero()])?;
        assert(soln).eq(vec![Complex::new(1.0, -1.0), Complex::new(0.0, -1.0)])?;
        Ok(())
    }
}
