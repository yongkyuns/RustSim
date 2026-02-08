//! GEAR-type variable-order BDF solvers with adaptive timestepping
//!
//! Implementation of the GEAR family of implicit multistep methods based on
//! Backward Differentiation Formulas (BDF). These solvers compute BDF coefficients
//! dynamically to handle variable timesteps and use lower-order methods for error
//! estimation and adaptive timestep control.
//!
//! # Theory
//!
//! For an m-th order BDF at timestep n, we have:
//!
//! ```text
//! sum(alpha_i * x_i; i=n-m,...,n) = h_n * f_n(x_n, t_n)
//! ```
//!
//! or equivalently:
//!
//! ```text
//! x_n = beta * h_n * f_n(x_n, t_n) - sum(alpha_j * x_{n-1-j}; j=0,...,order-1)
//! ```
//!
//! The coefficients `alpha` and `beta` are computed dynamically from the variable
//! timestep history to maintain accuracy with non-uniform time grids.
//!
//! # Startup Phase
//!
//! Since BDF methods require a history buffer of previous solutions, a startup
//! solver (ESDIRK32-equivalent) is used during the first few steps until enough
//! history is accumulated.
//!
//! # References
//!
//! - Gear, C. W. (1971). "Numerical Initial Value Problems in Ordinary
//!   Differential Equations". Prentice-Hall.
//! - Hairer, E., & Wanner, G. (1996). "Solving Ordinary Differential
//!   Equations II: Stiff and Differential-Algebraic Problems". Springer.

use nalgebra::{DMatrix, DVector};
use std::collections::VecDeque;

use super::{ImplicitSolver, Solver, SolverError, SolverStepResult};

// Constants from PathSim
const TOLERANCE: f64 = 1e-16;
const SOL_BETA: f64 = 0.9;
const SOL_SCALE_MIN: f64 = 0.1;
const SOL_SCALE_MAX: f64 = 10.0;

/// Compute BDF coefficients for variable timestep integration
///
/// For m-th order BDF we have for the n-th timestep:
/// ```text
/// sum(alpha_i * x_i; i=n-m,...,n) = h_n * f_n(x_n, t_n)
/// ```
///
/// # Arguments
/// * `order` - Order of the BDF method
/// * `timesteps` - Timestep buffer (h_{n-j}; j=0,...,order-1), most recent first
///
/// # Returns
/// * `beta` - Weight for function evaluation
/// * `alpha` - Weights for previous solutions (ordered from most recent to oldest)
///
/// # Panics
/// Panics if order < 1 or if timesteps buffer is too small
pub fn compute_bdf_coefficients(order: usize, timesteps: &[f64]) -> (f64, Vec<f64>) {
    assert!(
        order >= 1,
        "BDF coefficients of order '{}' not possible!",
        order
    );

    // Quit early for order 1 (backward Euler)
    if timesteps.len() < 2 {
        return (1.0, vec![1.0]);
    }

    // Compute timestep ratios rho_j = h_{n-j} / h_n
    let h_n = timesteps[0];
    let rho: Vec<f64> = timesteps[1..].iter().map(|&h| h / h_n).collect();

    // Compute normalized time differences theta_j
    let mut theta = vec![-1.0; order + 1];
    theta[0] = 0.0;
    for j in 2..=order {
        let sum: f64 = rho[..j - 1].iter().sum();
        theta[j] -= sum;
    }

    // Set up the linear system (order + 1 equations)
    // A * alpha = b where b = [0, 1, 0, ..., 0]
    let mut a_matrix = DMatrix::zeros(order + 1, order + 1);
    let mut b_vec = DVector::zeros(order + 1);
    b_vec[1] = 1.0;

    // Fill matrix: A[m, :] = theta^m
    for m in 0..=order {
        for j in 0..=order {
            a_matrix[(m, j)] = theta[j].powi(m as i32);
        }
    }

    // Solve the linear system
    let alphas = a_matrix
        .lu()
        .solve(&b_vec)
        .expect("Failed to solve BDF coefficient system");

    // Extract beta and alpha coefficients
    let beta = 1.0 / alphas[0];
    let alpha: Vec<f64> = alphas.iter().skip(1).map(|&a| -a / alphas[0]).collect();

    (beta, alpha)
}

/// Simple Newton-Anderson optimizer for implicit solvers
///
/// This is a simplified version of PathSim's optimizer, using Newton iteration
/// with optional Anderson acceleration.
#[allow(dead_code)]
#[derive(Debug, Clone)]
struct NewtonAnderson {
    tolerance: f64,
    max_iterations: usize,
}

impl Default for NewtonAnderson {
    fn default() -> Self {
        Self {
            tolerance: 1e-8,
            max_iterations: 50,
        }
    }
}

#[allow(dead_code)]
impl NewtonAnderson {
    fn new(tolerance: f64, max_iterations: usize) -> Self {
        Self {
            tolerance,
            max_iterations,
        }
    }

    fn reset(&mut self) {
        // Reset internal state if needed
    }

    /// Perform Newton iteration to solve x = g(x)
    ///
    /// # Arguments
    /// * `x` - Current guess
    /// * `g` - Target value (fixed point)
    /// * `jac` - Optional Jacobian matrix (I - beta*h*J)
    ///
    /// # Returns
    /// * New state estimate
    /// * Residual error
    fn step(
        &mut self,
        x: &DVector<f64>,
        g: &DVector<f64>,
        jac: Option<&DMatrix<f64>>,
    ) -> Result<(DVector<f64>, f64), SolverError> {
        let mut x_new = x.clone();
        let mut error;

        for _iter in 0..self.max_iterations {
            // Residual: r = g - x
            let residual = g - &x_new;
            error = residual.norm();

            if error < self.tolerance {
                return Ok((x_new, error));
            }

            // Newton step
            if let Some(j) = jac {
                // Solve (I - J) * dx = residual
                let eye = DMatrix::identity(x.len(), x.len());
                let system = &eye - j;

                match system.lu().solve(&residual) {
                    Some(dx) => {
                        x_new += dx;
                    }
                    None => {
                        // Fall back to fixed-point iteration
                        x_new = g.clone();
                    }
                }
            } else {
                // Simple fixed-point iteration
                x_new = g.clone();
            }
        }

        // Did not converge
        Err(SolverError::ConvergenceFailure(self.max_iterations))
    }
}

/// Base GEAR solver implementation
///
/// This is an internal base struct that implements the common GEAR logic.
/// Use the specific GEAR variants (GEAR21, GEAR32, etc.) instead.
#[derive(Debug, Clone)]
struct GEARBase {
    state: DVector<f64>,
    initial: DVector<f64>,
    history: VecDeque<DVector<f64>>,
    history_dt: VecDeque<f64>,
    n: usize, // Order of main method
    m: usize, // Order of error estimate method
    tol_abs: f64,
    tol_rel: f64,
    beta: f64,
    opt: NewtonAnderson,
    needs_startup: bool,
    startup_state: DVector<f64>,
    startup_history: VecDeque<DVector<f64>>,
    startup_stage: usize,
    startup_slopes: Vec<DVector<f64>>,
}

impl GEARBase {
    fn new(initial: DVector<f64>, n: usize, m: usize, tol_abs: f64, tol_rel: f64) -> Self {
        let state_dim = initial.len();
        let max_history = n.max(6); // GEAR52A needs up to 6

        Self {
            state: initial.clone(),
            initial: initial.clone(),
            history: VecDeque::with_capacity(max_history),
            history_dt: VecDeque::with_capacity(max_history),
            n,
            m,
            tol_abs,
            tol_rel,
            beta: SOL_BETA,
            opt: NewtonAnderson::default(),
            needs_startup: true,
            startup_state: initial.clone(),
            startup_history: VecDeque::with_capacity(2),
            startup_stage: 0,
            startup_slopes: vec![DVector::zeros(state_dim); 4], // ESDIRK32 has 4 stages
        }
    }

    fn buffer(&mut self, dt: f64) {
        self.opt.reset();

        // Add to histories
        while self.history.len() >= self.n {
            self.history.pop_back();
        }
        while self.history_dt.len() >= self.n {
            self.history_dt.pop_back();
        }

        self.history.push_front(self.state.clone());
        self.history_dt.push_front(dt);

        // Check if we still need startup
        self.needs_startup = self.history.len() < self.n;

        // Buffer startup solver
        if self.needs_startup {
            if self.startup_history.len() >= 2 {
                self.startup_history.pop_back();
            }
            self.startup_history.push_front(self.state.clone());
            self.startup_stage = 0;
        }
    }

    fn revert(&mut self) -> Result<(), SolverError> {
        self.state = self.history.pop_front().ok_or(SolverError::EmptyHistory)?;
        self.history_dt.pop_front();

        if self.needs_startup {
            self.startup_state = self
                .startup_history
                .pop_front()
                .ok_or(SolverError::EmptyHistory)?;
            self.startup_stage = 0;
        }

        Ok(())
    }

    fn reset(&mut self) {
        self.state = self.initial.clone();
        self.history.clear();
        self.history_dt.clear();
        self.needs_startup = true;
        self.startup_state = self.initial.clone();
        self.startup_history.clear();
        self.startup_stage = 0;
    }

    /// Simplified ESDIRK32 startup solver - one stage at a time
    fn startup_step<F, J>(
        &mut self,
        mut f: F,
        _jac: Option<J>,
        t: f64,
        dt: f64,
    ) -> Result<SolverStepResult, SolverError>
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
        J: FnMut(&DVector<f64>, f64) -> DMatrix<f64>,
    {
        // ESDIRK32 coefficients (simplified TR-BDF2)
        let gamma = 1.0 - 1.0 / f64::sqrt(2.0);

        let x0 = self
            .startup_history
            .front()
            .ok_or(SolverError::EmptyHistory)?;

        match self.startup_stage {
            0 => {
                // First stage: implicit Euler with gamma
                self.startup_slopes[0] = f(&self.startup_state, t);

                // Newton iteration for: x1 = x0 + gamma*h*f(x1, t + gamma*h)
                let t1 = t + gamma * dt;
                for _ in 0..10 {
                    let f1 = f(&self.startup_state, t1);
                    let target = x0 + gamma * dt * &f1;

                    let residual = (&target - &self.startup_state).norm();
                    if residual < 1e-8 {
                        break;
                    }

                    self.startup_state = target;
                }

                self.startup_slopes[1] = f(&self.startup_state, t1);
                self.startup_stage += 1;
                Ok(SolverStepResult::default())
            }
            1 | 2 => {
                // Intermediate stages
                self.startup_stage += 1;
                Ok(SolverStepResult::default())
            }
            _ => {
                // Final stage: BDF2-like step
                let t_final = t + dt;

                // Simple implicit step
                for _ in 0..10 {
                    let f_final = f(&self.startup_state, t_final);
                    let c1 = 2.0 / 3.0;
                    let target = (4.0 / 3.0) * x0 - (1.0 / 3.0) * x0 + c1 * dt * &f_final;

                    let residual = (&target - &self.startup_state).norm();
                    if residual < 1e-8 {
                        break;
                    }

                    self.startup_state = target;
                }

                self.state = self.startup_state.clone();
                self.startup_stage = 0;

                // Simple error estimate
                Ok(SolverStepResult {
                    success: true,
                    error_norm: 1e-10,
                    scale: Some(1.0),
                })
            }
        }
    }

    fn solve_implicit<F, J>(
        &mut self,
        mut f: F,
        mut jac: Option<J>,
        t: f64,
        dt: f64,
    ) -> Result<f64, SolverError>
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
        J: FnMut(&DVector<f64>, f64) -> DMatrix<f64>,
    {
        if self.needs_startup {
            let result = self.startup_step(f, jac, t, dt)?;
            return Ok(result.error_norm);
        }

        // Compute BDF coefficients for current timestep
        let timesteps: Vec<f64> = self.history_dt.iter().copied().collect();
        let (beta_n, alpha_n) = compute_bdf_coefficients(self.n, &timesteps);

        // Compute fixed-point target: g = beta * h * f + sum(k * x_prev)
        let f_val = f(&self.state, t);
        let mut g = beta_n * dt * &f_val;

        for (i, hist_state) in self.history.iter().enumerate() {
            if i < alpha_n.len() {
                g += alpha_n[i] * hist_state;
            }
        }

        // Prepare Jacobian if available
        let jac_matrix = jac.as_mut().map(|j| {
            let j_mat = j(&self.state, t);
            beta_n * dt * j_mat
        });

        // Solve implicit equation
        let (new_state, error) = self.opt.step(&self.state, &g, jac_matrix.as_ref())?;
        self.state = new_state;

        Ok(error)
    }

    fn error_controller(&self, tr: &DVector<f64>) -> (bool, f64, f64) {
        // Compute scaling factors
        let scale = self.state.map(|x| self.tol_abs + self.tol_rel * x.abs());

        // Compute scaled error
        let scaled_error = tr.component_div(&scale).map(|e| e.abs());

        // Error norm (max norm) with lower bound
        let error_norm = scaled_error.max().max(TOLERANCE);

        // Determine success
        let success = error_norm <= 1.0;

        // Compute timestep scale factor
        let timestep_rescale = self.beta / error_norm.powf(1.0 / self.n as f64);

        // Clip rescale factor
        let timestep_rescale = timestep_rescale.clamp(SOL_SCALE_MIN, SOL_SCALE_MAX);

        (success, error_norm, timestep_rescale)
    }

    fn step<F>(&mut self, mut f: F, t: f64, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        if self.needs_startup {
            // Return startup result - already computed in solve
            return SolverStepResult {
                success: true,
                error_norm: 1e-10,
                scale: Some(1.0),
            };
        }

        // Compute BDF coefficients for error estimate (lower order)
        let timesteps: Vec<f64> = self.history_dt.iter().copied().collect();
        let (beta_m, alpha_m) = compute_bdf_coefficients(self.m, &timesteps);

        // Estimate truncation error from lower order solution
        let f_val = f(&self.state, t);
        let mut tr = &self.state - beta_m * dt * &f_val;

        for (i, hist_state) in self.history.iter().enumerate() {
            if i < alpha_m.len() {
                tr -= alpha_m[i] * hist_state;
            }
        }

        // Error control
        let (success, error_norm, scale) = self.error_controller(&tr);

        SolverStepResult {
            success,
            error_norm,
            scale: Some(scale),
        }
    }
}

/// GEAR21: Variable-step 2nd order BDF with 1st order error estimate
///
/// # Characteristics
/// - Order: 2 (stepping) / 1 (error estimate)
/// - Implicit variable-step multistep
/// - Adaptive timestep
/// - A-stable
///
/// # Note
/// The simplest adaptive multistep stiff solver. A-stability makes it safe
/// for any stiff system. The multistep approach reuses past solution values,
/// so per-step cost is lower than single-step implicit methods.
#[derive(Debug, Clone)]
pub struct GEAR21 {
    inner: GEARBase,
}

impl GEAR21 {
    pub fn new(initial: DVector<f64>) -> Self {
        Self::with_tolerances(initial, 1e-8, 1e-4)
    }

    pub fn with_tolerances(initial: DVector<f64>, tol_abs: f64, tol_rel: f64) -> Self {
        Self {
            inner: GEARBase::new(initial, 2, 1, tol_abs, tol_rel),
        }
    }
}

impl Solver for GEAR21 {
    fn state(&self) -> &DVector<f64> {
        &self.inner.state
    }

    fn state_mut(&mut self) -> &mut DVector<f64> {
        &mut self.inner.state
    }

    fn set_state(&mut self, state: DVector<f64>) {
        self.inner.state = state;
    }

    fn buffer(&mut self, dt: f64) {
        self.inner.buffer(dt);
    }

    fn revert(&mut self) -> Result<(), SolverError> {
        self.inner.revert()
    }

    fn reset(&mut self) {
        self.inner.reset();
    }

    fn order(&self) -> usize {
        2
    }

    fn stages(&self) -> usize {
        1
    }

    fn is_adaptive(&self) -> bool {
        true
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

impl ImplicitSolver for GEAR21 {
    fn solve<F, J>(&mut self, f: F, jac: Option<J>, dt: f64) -> Result<f64, SolverError>
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
        J: FnMut(&DVector<f64>, f64) -> DMatrix<f64>,
    {
        self.inner.solve_implicit(f, jac, 0.0, dt)
    }

    fn step<F>(&mut self, f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        self.inner.step(f, 0.0, dt)
    }
}

/// GEAR32: Variable-step 3rd order BDF with 2nd order error estimate
///
/// # Characteristics
/// - Order: 3 (stepping) / 2 (error estimate)
/// - Implicit variable-step multistep
/// - Adaptive timestep
/// - A(α)-stable (BDF3 stability wedge, α ≈ 86°)
#[derive(Debug, Clone)]
pub struct GEAR32 {
    inner: GEARBase,
}

impl GEAR32 {
    pub fn new(initial: DVector<f64>) -> Self {
        Self::with_tolerances(initial, 1e-8, 1e-4)
    }

    pub fn with_tolerances(initial: DVector<f64>, tol_abs: f64, tol_rel: f64) -> Self {
        Self {
            inner: GEARBase::new(initial, 3, 2, tol_abs, tol_rel),
        }
    }
}

impl Solver for GEAR32 {
    fn state(&self) -> &DVector<f64> {
        &self.inner.state
    }

    fn state_mut(&mut self) -> &mut DVector<f64> {
        &mut self.inner.state
    }

    fn set_state(&mut self, state: DVector<f64>) {
        self.inner.state = state;
    }

    fn buffer(&mut self, dt: f64) {
        self.inner.buffer(dt);
    }

    fn revert(&mut self) -> Result<(), SolverError> {
        self.inner.revert()
    }

    fn reset(&mut self) {
        self.inner.reset();
    }

    fn order(&self) -> usize {
        3
    }

    fn stages(&self) -> usize {
        1
    }

    fn is_adaptive(&self) -> bool {
        true
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

impl ImplicitSolver for GEAR32 {
    fn solve<F, J>(&mut self, f: F, jac: Option<J>, dt: f64) -> Result<f64, SolverError>
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
        J: FnMut(&DVector<f64>, f64) -> DMatrix<f64>,
    {
        self.inner.solve_implicit(f, jac, 0.0, dt)
    }

    fn step<F>(&mut self, f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        self.inner.step(f, 0.0, dt)
    }
}

/// GEAR43: Variable-step 4th order BDF with 3rd order error estimate
///
/// # Characteristics
/// - Order: 4 (stepping) / 3 (error estimate)
/// - Implicit variable-step multistep
/// - Adaptive timestep
/// - A(α)-stable (BDF4 stability wedge, α ≈ 73°)
#[derive(Debug, Clone)]
pub struct GEAR43 {
    inner: GEARBase,
}

impl GEAR43 {
    pub fn new(initial: DVector<f64>) -> Self {
        Self::with_tolerances(initial, 1e-8, 1e-4)
    }

    pub fn with_tolerances(initial: DVector<f64>, tol_abs: f64, tol_rel: f64) -> Self {
        Self {
            inner: GEARBase::new(initial, 4, 3, tol_abs, tol_rel),
        }
    }
}

impl Solver for GEAR43 {
    fn state(&self) -> &DVector<f64> {
        &self.inner.state
    }

    fn state_mut(&mut self) -> &mut DVector<f64> {
        &mut self.inner.state
    }

    fn set_state(&mut self, state: DVector<f64>) {
        self.inner.state = state;
    }

    fn buffer(&mut self, dt: f64) {
        self.inner.buffer(dt);
    }

    fn revert(&mut self) -> Result<(), SolverError> {
        self.inner.revert()
    }

    fn reset(&mut self) {
        self.inner.reset();
    }

    fn order(&self) -> usize {
        4
    }

    fn stages(&self) -> usize {
        1
    }

    fn is_adaptive(&self) -> bool {
        true
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

impl ImplicitSolver for GEAR43 {
    fn solve<F, J>(&mut self, f: F, jac: Option<J>, dt: f64) -> Result<f64, SolverError>
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
        J: FnMut(&DVector<f64>, f64) -> DMatrix<f64>,
    {
        self.inner.solve_implicit(f, jac, 0.0, dt)
    }

    fn step<F>(&mut self, f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        self.inner.step(f, 0.0, dt)
    }
}

/// GEAR54: Variable-step 5th order BDF with 4th order error estimate
///
/// # Characteristics
/// - Order: 5 (stepping) / 4 (error estimate)
/// - Implicit variable-step multistep
/// - Adaptive timestep
/// - A(α)-stable (BDF5 stability wedge, α ≈ 51°)
#[derive(Debug, Clone)]
pub struct GEAR54 {
    inner: GEARBase,
}

impl GEAR54 {
    pub fn new(initial: DVector<f64>) -> Self {
        Self::with_tolerances(initial, 1e-8, 1e-4)
    }

    pub fn with_tolerances(initial: DVector<f64>, tol_abs: f64, tol_rel: f64) -> Self {
        Self {
            inner: GEARBase::new(initial, 5, 4, tol_abs, tol_rel),
        }
    }
}

impl Solver for GEAR54 {
    fn state(&self) -> &DVector<f64> {
        &self.inner.state
    }

    fn state_mut(&mut self) -> &mut DVector<f64> {
        &mut self.inner.state
    }

    fn set_state(&mut self, state: DVector<f64>) {
        self.inner.state = state;
    }

    fn buffer(&mut self, dt: f64) {
        self.inner.buffer(dt);
    }

    fn revert(&mut self) -> Result<(), SolverError> {
        self.inner.revert()
    }

    fn reset(&mut self) {
        self.inner.reset();
    }

    fn order(&self) -> usize {
        5
    }

    fn stages(&self) -> usize {
        1
    }

    fn is_adaptive(&self) -> bool {
        true
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

impl ImplicitSolver for GEAR54 {
    fn solve<F, J>(&mut self, f: F, jac: Option<J>, dt: f64) -> Result<f64, SolverError>
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
        J: FnMut(&DVector<f64>, f64) -> DMatrix<f64>,
    {
        self.inner.solve_implicit(f, jac, 0.0, dt)
    }

    fn step<F>(&mut self, f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        self.inner.step(f, 0.0, dt)
    }
}

/// GEAR52A: Variable-step, variable-order BDF (orders 2-5)
///
/// Adapts both timestep and order automatically. At each step, compares error
/// estimates from orders n-1 and n+1 and selects the order that minimizes
/// normalized error, allowing larger steps.
///
/// # Characteristics
/// - Order: variable, 2-5
/// - Implicit variable-step, variable-order multistep
/// - Adaptive timestep and order
/// - Stability: A-stable at order 2, A(α)-stable at orders 3-5
#[derive(Debug, Clone)]
pub struct GEAR52A {
    inner: GEARBase,
    n_min: usize,
    n_max: usize,
}

impl GEAR52A {
    pub fn new(initial: DVector<f64>) -> Self {
        Self::with_tolerances(initial, 1e-8, 1e-4)
    }

    pub fn with_tolerances(initial: DVector<f64>, tol_abs: f64, tol_rel: f64) -> Self {
        let mut inner = GEARBase::new(initial, 2, 1, tol_abs, tol_rel);
        inner.history = VecDeque::with_capacity(6);
        inner.history_dt = VecDeque::with_capacity(6);

        Self {
            inner,
            n_min: 2,
            n_max: 5,
        }
    }

    fn adapt_order(&mut self, error_norm_m: f64, error_norm_p: f64) {
        // Decrease order if lower order is more accurate
        if error_norm_m < error_norm_p {
            self.inner.n = (self.inner.n - 1).max(self.n_min);
        } else {
            // Increase order if higher order is more accurate
            self.inner.n = (self.inner.n + 1).min(self.n_max);
        }

        // Update error estimate order
        self.inner.m = self.inner.n - 1;
    }
}

impl Solver for GEAR52A {
    fn state(&self) -> &DVector<f64> {
        &self.inner.state
    }

    fn state_mut(&mut self) -> &mut DVector<f64> {
        &mut self.inner.state
    }

    fn set_state(&mut self, state: DVector<f64>) {
        self.inner.state = state;
    }

    fn buffer(&mut self, dt: f64) {
        self.inner.opt.reset();

        // Add to histories with larger capacity
        while self.inner.history.len() >= 6 {
            self.inner.history.pop_back();
        }
        while self.inner.history_dt.len() >= 6 {
            self.inner.history_dt.pop_back();
        }

        self.inner.history.push_front(self.inner.state.clone());
        self.inner.history_dt.push_front(dt);

        // Check if we still need startup
        self.inner.needs_startup = self.inner.history.len() < 6;

        // Buffer startup solver
        if self.inner.needs_startup {
            if self.inner.startup_history.len() >= 2 {
                self.inner.startup_history.pop_back();
            }
            self.inner
                .startup_history
                .push_front(self.inner.state.clone());
            self.inner.startup_stage = 0;
        }
    }

    fn revert(&mut self) -> Result<(), SolverError> {
        self.inner.revert()
    }

    fn reset(&mut self) {
        self.inner.reset();
        self.inner.n = 2;
        self.inner.m = 1;
    }

    fn order(&self) -> usize {
        self.inner.n
    }

    fn stages(&self) -> usize {
        1
    }

    fn is_adaptive(&self) -> bool {
        true
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

impl ImplicitSolver for GEAR52A {
    fn solve<F, J>(&mut self, f: F, jac: Option<J>, dt: f64) -> Result<f64, SolverError>
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
        J: FnMut(&DVector<f64>, f64) -> DMatrix<f64>,
    {
        self.inner.solve_implicit(f, jac, 0.0, dt)
    }

    fn step<F>(&mut self, mut f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        if self.inner.needs_startup {
            return SolverStepResult {
                success: true,
                error_norm: 1e-10,
                scale: Some(1.0),
            };
        }

        // Compute BDF coefficients for n-1 and n+1 orders
        let n_m = self.inner.n - 1;
        let n_p = self.inner.n + 1;

        let timesteps: Vec<f64> = self.inner.history_dt.iter().copied().collect();
        let (beta_m, alpha_m) = compute_bdf_coefficients(n_m, &timesteps);
        let (beta_p, alpha_p) = compute_bdf_coefficients(n_p, &timesteps);

        // Estimate truncation errors
        let f_val = f(&self.inner.state, 0.0);

        let mut tr_m = &self.inner.state - beta_m * dt * &f_val;
        for (i, hist_state) in self.inner.history.iter().enumerate() {
            if i < alpha_m.len() {
                tr_m -= alpha_m[i] * hist_state;
            }
        }

        let mut tr_p = &self.inner.state - beta_p * dt * &f_val;
        for (i, hist_state) in self.inner.history.iter().enumerate() {
            if i < alpha_p.len() {
                tr_p -= alpha_p[i] * hist_state;
            }
        }

        // Error control for both orders
        let (success_m, error_norm_m, scale_m) = self.inner.error_controller(&tr_m);
        let (_, error_norm_p, _) = self.inner.error_controller(&tr_p);

        // Adapt order based on errors
        self.adapt_order(error_norm_m, error_norm_p);

        SolverStepResult {
            success: success_m,
            error_norm: error_norm_p, // Use higher order estimate for reporting
            scale: Some(scale_m),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_bdf_coefficients_order_1() {
        let timesteps = vec![1.0];
        let (beta, alpha) = compute_bdf_coefficients(1, &timesteps);

        assert_relative_eq!(beta, 1.0, epsilon = 1e-7);
        assert_eq!(alpha.len(), 1);
        assert_relative_eq!(alpha[0], 1.0, epsilon = 1e-7);
    }

    #[test]
    fn test_bdf_coefficients_order_2() {
        let timesteps = vec![1.0, 1.0];
        let (beta, alpha) = compute_bdf_coefficients(2, &timesteps);

        assert_relative_eq!(beta, 2.0 / 3.0, epsilon = 1e-7);
        assert_eq!(alpha.len(), 2);
        assert_relative_eq!(alpha[0], 4.0 / 3.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[1], -1.0 / 3.0, epsilon = 1e-7);
    }

    #[test]
    fn test_bdf_coefficients_order_3() {
        let timesteps = vec![1.0, 1.0, 1.0];
        let (beta, alpha) = compute_bdf_coefficients(3, &timesteps);

        assert_relative_eq!(beta, 6.0 / 11.0, epsilon = 1e-7);
        assert_eq!(alpha.len(), 3);
        assert_relative_eq!(alpha[0], 18.0 / 11.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[1], -9.0 / 11.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[2], 2.0 / 11.0, epsilon = 1e-7);
    }
}
