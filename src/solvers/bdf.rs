//! Backward Differentiation Formula (BDF) family of implicit solvers
//!
//! BDF methods are linear multistep methods particularly suited for stiff ODEs.
//! They require a startup solver (DIRK3) for the first few steps to build up
//! the history buffer needed for the multistep formulas.
//!
//! # Fixed-timestep BDF Methods
//!
//! This module implements fixed-timestep BDF solvers of orders 2-6:
//! - BDF2: 2nd order, A-stable
//! - BDF3: 3rd order, A(α)-stable with α ≈ 86°
//! - BDF4: 4th order, A(α)-stable with α ≈ 73°
//! - BDF5: 5th order, A(α)-stable with α ≈ 51°
//! - BDF6: 6th order, A(α)-stable with α ≈ 18° (NOT A-stable)
//!
//! # Startup Method
//!
//! All BDF methods use DIRK3 (a 4-stage, 3rd order L-stable method) for the
//! initial steps until enough history is accumulated for the full BDF method.

use nalgebra::{DMatrix, DVector};
use std::collections::VecDeque;

use super::dirk3::DIRK3;
use super::{ImplicitSolver, Solver, SolverError, SolverStepResult};
use crate::optim::anderson::Anderson;

/// Base BDF solver structure
///
/// Implements the Backward Differentiation Formula family of methods.
#[derive(Debug, Clone)]
struct BDFBase {
    state: DVector<f64>,
    initial: DVector<f64>,
    history: VecDeque<DVector<f64>>,
    order: usize,
    max_history: usize,
    anderson: Anderson,
    max_iter: usize,
    tolerance: f64,
    startup_solver: Option<Box<DIRK3>>,
    needs_startup: bool,
    k_coeffs: Vec<f64>,  // State coefficients
    f_coeff: f64,         // Function coefficient
    current_dt: f64,      // Current timestep
}

impl BDFBase {
    /// Create new BDF solver with specified order
    fn new(initial: DVector<f64>, order: usize) -> Self {
        // BDF coefficients for orders 1-6
        let (k_coeffs, f_coeff) = match order {
            1 => (vec![1.0], 1.0),
            2 => (vec![4.0 / 3.0, -1.0 / 3.0], 2.0 / 3.0),
            3 => (vec![18.0 / 11.0, -9.0 / 11.0, 2.0 / 11.0], 6.0 / 11.0),
            4 => (
                vec![48.0 / 25.0, -36.0 / 25.0, 16.0 / 25.0, -3.0 / 25.0],
                12.0 / 25.0,
            ),
            5 => (
                vec![
                    300.0 / 137.0,
                    -300.0 / 137.0,
                    200.0 / 137.0,
                    -75.0 / 137.0,
                    12.0 / 137.0,
                ],
                60.0 / 137.0,
            ),
            6 => (
                vec![
                    360.0 / 147.0,
                    -450.0 / 147.0,
                    400.0 / 147.0,
                    -225.0 / 147.0,
                    72.0 / 147.0,
                    -10.0 / 147.0,
                ],
                60.0 / 147.0,
            ),
            _ => panic!("BDF order must be 1-6"),
        };

        let startup_solver = Some(Box::new(DIRK3::new(initial.clone())));

        Self {
            state: initial.clone(),
            initial,
            history: VecDeque::with_capacity(order),
            order,
            max_history: order,
            anderson: Anderson::new(3, false),
            max_iter: 50,
            tolerance: 1e-10,
            startup_solver,
            needs_startup: true,
            k_coeffs,
            f_coeff,
            current_dt: 0.0,
        }
    }
}

impl Solver for BDFBase {
    fn state(&self) -> &DVector<f64> {
        &self.state
    }

    fn state_mut(&mut self) -> &mut DVector<f64> {
        &mut self.state
    }

    fn set_state(&mut self, state: DVector<f64>) {
        self.state = state;
    }

    fn buffer(&mut self, dt: f64) {
        self.current_dt = dt;

        // Reset Anderson accelerator
        self.anderson.reset();

        // Add current state to history
        if self.history.len() >= self.max_history {
            self.history.pop_back();
        }
        self.history.push_front(self.state.clone());

        // Check if we still need startup
        self.needs_startup = self.history.len() < self.order;

        // Buffer startup solver if needed
        if self.needs_startup {
            if let Some(ref mut startup) = self.startup_solver {
                startup.buffer(dt);
            }
        }
    }

    fn revert(&mut self) -> Result<(), SolverError> {
        self.state = self.history.pop_front().ok_or(SolverError::EmptyHistory)?;
        Ok(())
    }

    fn reset(&mut self) {
        self.state = self.initial.clone();
        self.history.clear();
        self.anderson.reset();
        self.needs_startup = true;
        if let Some(ref mut startup) = self.startup_solver {
            startup.reset();
        }
    }

    fn order(&self) -> usize {
        self.order
    }

    fn stages(&self) -> usize {
        if self.needs_startup {
            4 // DIRK3 has 4 stages
        } else {
            1 // BDF is single-stage
        }
    }

    fn is_adaptive(&self) -> bool {
        false
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

impl ImplicitSolver for BDFBase {
    fn solve(&mut self, f: &DVector<f64>, _jac: Option<&DMatrix<f64>>, dt: f64) -> Result<f64, SolverError> {
        // Use startup solver if we don't have enough history
        if self.needs_startup {
            if let Some(ref mut startup) = self.startup_solver {
                let err = startup.solve(f, None, dt)?;
                self.state = startup.state().clone();
                return Ok(err);
            }
        }

        // Fixed-point iteration: x_{n+1} = g(f(x_{n+1}))
        // where g(f) = sum(k_i * x_{n-i}) + f_coeff * dt * f
        let mut x = self.state.clone();
        let mut residual = f64::INFINITY;

        for _iter in 0..self.max_iter {
            // Compute g(f) = sum of history terms + dt * f
            let mut g = self.f_coeff * dt * f;

            for (i, k) in self.k_coeffs.iter().enumerate() {
                if i < self.history.len() {
                    g += k * &self.history[i];
                }
            }

            // Anderson acceleration step
            let (x_new, res_norm) = self.anderson.step(&x, &g);
            residual = res_norm;

            // Check convergence
            if residual < self.tolerance {
                self.state = x_new;
                return Ok(residual);
            }

            x = x_new;
        }

        // Failed to converge
        Err(SolverError::ConvergenceFailure(self.max_iter))
    }

    fn step(&mut self, _f: &DVector<f64>, _dt: f64) -> SolverStepResult {
        // For BDF, the step is handled in solve()
        // This is just a placeholder to satisfy the trait
        SolverStepResult::default()
    }
}

// BDF2: 2nd order, A-stable
/// BDF2: 2nd order Backward Differentiation Formula
///
/// # Formula
/// x_{n+1} = 4/3 x_n - 1/3 x_{n-1} + 2/3 h f(x_{n+1})
///
/// # Characteristics
/// - Order: 2
/// - A-stable (excellent for stiff problems)
/// - Requires 1 previous step
///
/// # Stability
/// The workhorse fixed-step stiff solver. A-stability means no eigenvalue
/// in the left half-plane causes instability, regardless of timestep.
#[derive(Debug, Clone)]
pub struct BDF2 {
    base: BDFBase,
}

impl BDF2 {
    pub fn new(initial: DVector<f64>) -> Self {
        Self {
            base: BDFBase::new(initial, 2),
        }
    }
}

impl Solver for BDF2 {
    fn state(&self) -> &DVector<f64> {
        self.base.state()
    }

    fn state_mut(&mut self) -> &mut DVector<f64> {
        self.base.state_mut()
    }

    fn set_state(&mut self, state: DVector<f64>) {
        self.base.set_state(state)
    }

    fn buffer(&mut self, dt: f64) {
        self.base.buffer(dt)
    }

    fn revert(&mut self) -> Result<(), SolverError> {
        self.base.revert()
    }

    fn reset(&mut self) {
        self.base.reset()
    }

    fn order(&self) -> usize {
        2
    }

    fn stages(&self) -> usize {
        self.base.stages()
    }

    fn is_adaptive(&self) -> bool {
        false
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

impl ImplicitSolver for BDF2 {
    fn solve(&mut self, f: &DVector<f64>, jac: Option<&DMatrix<f64>>, dt: f64) -> Result<f64, SolverError> {
        self.base.solve(f, jac, dt)
    }

    fn step(&mut self, f: &DVector<f64>, dt: f64) -> SolverStepResult {
        self.base.step(f, dt)
    }
}

// BDF3: 3rd order
/// BDF3: 3rd order Backward Differentiation Formula
///
/// # Characteristics
/// - Order: 3
/// - A(alpha)-stable with alpha ≈ 86°
/// - Requires 2 previous steps
#[derive(Debug, Clone)]
pub struct BDF3 {
    base: BDFBase,
}

impl BDF3 {
    pub fn new(initial: DVector<f64>) -> Self {
        Self {
            base: BDFBase::new(initial, 3),
        }
    }
}

impl Solver for BDF3 {
    fn state(&self) -> &DVector<f64> {
        self.base.state()
    }

    fn state_mut(&mut self) -> &mut DVector<f64> {
        self.base.state_mut()
    }

    fn set_state(&mut self, state: DVector<f64>) {
        self.base.set_state(state)
    }

    fn buffer(&mut self, dt: f64) {
        self.base.buffer(dt)
    }

    fn revert(&mut self) -> Result<(), SolverError> {
        self.base.revert()
    }

    fn reset(&mut self) {
        self.base.reset()
    }

    fn order(&self) -> usize {
        3
    }

    fn stages(&self) -> usize {
        self.base.stages()
    }

    fn is_adaptive(&self) -> bool {
        false
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

impl ImplicitSolver for BDF3 {
    fn solve(&mut self, f: &DVector<f64>, jac: Option<&DMatrix<f64>>, dt: f64) -> Result<f64, SolverError> {
        self.base.solve(f, jac, dt)
    }

    fn step(&mut self, f: &DVector<f64>, dt: f64) -> SolverStepResult {
        self.base.step(f, dt)
    }
}

// BDF4: 4th order
/// BDF4: 4th order Backward Differentiation Formula
///
/// # Characteristics
/// - Order: 4
/// - A(alpha)-stable with alpha ≈ 73°
/// - Requires 3 previous steps
#[derive(Debug, Clone)]
pub struct BDF4 {
    base: BDFBase,
}

impl BDF4 {
    pub fn new(initial: DVector<f64>) -> Self {
        Self {
            base: BDFBase::new(initial, 4),
        }
    }
}

impl Solver for BDF4 {
    fn state(&self) -> &DVector<f64> {
        self.base.state()
    }

    fn state_mut(&mut self) -> &mut DVector<f64> {
        self.base.state_mut()
    }

    fn set_state(&mut self, state: DVector<f64>) {
        self.base.set_state(state)
    }

    fn buffer(&mut self, dt: f64) {
        self.base.buffer(dt)
    }

    fn revert(&mut self) -> Result<(), SolverError> {
        self.base.revert()
    }

    fn reset(&mut self) {
        self.base.reset()
    }

    fn order(&self) -> usize {
        4
    }

    fn stages(&self) -> usize {
        self.base.stages()
    }

    fn is_adaptive(&self) -> bool {
        false
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

impl ImplicitSolver for BDF4 {
    fn solve(&mut self, f: &DVector<f64>, jac: Option<&DMatrix<f64>>, dt: f64) -> Result<f64, SolverError> {
        self.base.solve(f, jac, dt)
    }

    fn step(&mut self, f: &DVector<f64>, dt: f64) -> SolverStepResult {
        self.base.step(f, dt)
    }
}

// BDF5: 5th order
/// BDF5: 5th order Backward Differentiation Formula
///
/// # Characteristics
/// - Order: 5
/// - A(alpha)-stable with alpha ≈ 51°
/// - Requires 4 previous steps
#[derive(Debug, Clone)]
pub struct BDF5 {
    base: BDFBase,
}

impl BDF5 {
    pub fn new(initial: DVector<f64>) -> Self {
        Self {
            base: BDFBase::new(initial, 5),
        }
    }
}

impl Solver for BDF5 {
    fn state(&self) -> &DVector<f64> {
        self.base.state()
    }

    fn state_mut(&mut self) -> &mut DVector<f64> {
        self.base.state_mut()
    }

    fn set_state(&mut self, state: DVector<f64>) {
        self.base.set_state(state)
    }

    fn buffer(&mut self, dt: f64) {
        self.base.buffer(dt)
    }

    fn revert(&mut self) -> Result<(), SolverError> {
        self.base.revert()
    }

    fn reset(&mut self) {
        self.base.reset()
    }

    fn order(&self) -> usize {
        5
    }

    fn stages(&self) -> usize {
        self.base.stages()
    }

    fn is_adaptive(&self) -> bool {
        false
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

impl ImplicitSolver for BDF5 {
    fn solve(&mut self, f: &DVector<f64>, jac: Option<&DMatrix<f64>>, dt: f64) -> Result<f64, SolverError> {
        self.base.solve(f, jac, dt)
    }

    fn step(&mut self, f: &DVector<f64>, dt: f64) -> SolverStepResult {
        self.base.step(f, dt)
    }
}

// BDF6: 6th order
/// BDF6: 6th order Backward Differentiation Formula
///
/// # Characteristics
/// - Order: 6
/// - A(alpha)-stable with alpha ≈ 18° (NOT A-stable)
/// - Requires 5 previous steps
///
/// # Warning
/// Very narrow stability wedge makes this unsuitable for most stiff problems.
/// Provided mainly for completeness.
#[derive(Debug, Clone)]
pub struct BDF6 {
    base: BDFBase,
}

impl BDF6 {
    pub fn new(initial: DVector<f64>) -> Self {
        Self {
            base: BDFBase::new(initial, 6),
        }
    }
}

impl Solver for BDF6 {
    fn state(&self) -> &DVector<f64> {
        self.base.state()
    }

    fn state_mut(&mut self) -> &mut DVector<f64> {
        self.base.state_mut()
    }

    fn set_state(&mut self, state: DVector<f64>) {
        self.base.set_state(state)
    }

    fn buffer(&mut self, dt: f64) {
        self.base.buffer(dt)
    }

    fn revert(&mut self) -> Result<(), SolverError> {
        self.base.revert()
    }

    fn reset(&mut self) {
        self.base.reset()
    }

    fn order(&self) -> usize {
        6
    }

    fn stages(&self) -> usize {
        self.base.stages()
    }

    fn is_adaptive(&self) -> bool {
        false
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

impl ImplicitSolver for BDF6 {
    fn solve(&mut self, f: &DVector<f64>, jac: Option<&DMatrix<f64>>, dt: f64) -> Result<f64, SolverError> {
        self.base.solve(f, jac, dt)
    }

    fn step(&mut self, f: &DVector<f64>, dt: f64) -> SolverStepResult {
        self.base.step(f, dt)
    }
}
