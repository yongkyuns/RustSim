//! Bogacki-Shampine 3(2) adaptive solver

use nalgebra::DVector;
use std::collections::VecDeque;

use super::{ExplicitSolver, Solver, SolverError, SolverStepResult};

/// Bogacki-Shampine 3(2) pair adaptive solver
///
/// Four-stage, 3rd order with FSAL (First-Same-As-Last) property.
///
/// The underlying method of MATLAB's `ode23`. The First-Same-As-Last
/// (FSAL) property makes the effective cost three stages per accepted step.
///
/// # Characteristics
/// - Order: 3 (propagating) / 2 (embedded)
/// - Stages: 4 (3 effective with FSAL)
/// - Explicit, adaptive timestep
///
/// # Note
/// A good default when moderate accuracy suffices and per-step cost matters
/// more than large step sizes. Fewer stages than 5th order pairs, so faster
/// per step but needs more steps for the same global error. In a RustSim
/// simulation with smooth, non-stiff dynamics and relaxed tolerances this
/// is often the most efficient explicit choice. Switch to RKDP54 when
/// tighter tolerances are required.
///
/// # References
/// - Bogacki, P., & Shampine, L. F. (1989). "A 3(2) pair of
///   Runge-Kutta formulas". Applied Mathematics Letters, 2(4),
///   321-325.
/// - Shampine, L. F., & Reichelt, M. W. (1997). "The MATLAB ODE
///   Suite". SIAM Journal on Scientific Computing, 18(1), 1-22.
#[derive(Debug, Clone)]
pub struct RKBS32 {
    state: DVector<f64>,
    initial: DVector<f64>,
    history: VecDeque<DVector<f64>>,
    slopes: Vec<DVector<f64>>,
    stage: usize,
    tol_abs: f64,
    tol_rel: f64,
    beta: f64,
}

impl RKBS32 {
    /// Create a new RKBS32 solver with the given initial state
    ///
    /// # Arguments
    /// * `initial` - Initial state vector
    pub fn new(initial: DVector<f64>) -> Self {
        Self::with_tolerances(initial, 1e-8, 1e-4)
    }

    /// Create a new RKBS32 solver with custom tolerances
    ///
    /// # Arguments
    /// * `initial` - Initial state vector
    /// * `tol_abs` - Absolute error tolerance
    /// * `tol_rel` - Relative error tolerance
    pub fn with_tolerances(initial: DVector<f64>, tol_abs: f64, tol_rel: f64) -> Self {
        let n = initial.len();
        Self {
            state: initial.clone(),
            initial,
            history: VecDeque::with_capacity(2),
            slopes: vec![DVector::zeros(n); 4],
            stage: 0,
            tol_abs,
            tol_rel,
            beta: 0.9, // Safety factor
        }
    }

    /// Compute error norm and timestep scale factor
    fn error_controller(&self, dt: f64) -> (bool, f64, f64) {
        // Coefficients for local truncation error estimate
        // TR = [-5/72, 1/12, 1/9, -1/8]
        let tr = [-5.0 / 72.0, 1.0 / 12.0, 1.0 / 9.0, -1.0 / 8.0];

        // Compute truncation error slope
        let mut error_slope = DVector::zeros(self.state.len());
        for (i, &coef) in tr.iter().enumerate() {
            error_slope += coef * &self.slopes[i];
        }

        // Compute scaling factors (avoid division by zero)
        let scale = self.state.map(|x| self.tol_abs + self.tol_rel * x.abs());

        // Compute scaled error (element-wise)
        let scaled_error = (dt * &error_slope).component_div(&scale).map(|e| e.abs());

        // Error norm (max norm) with lower bound
        let error_norm = scaled_error.max().max(1e-16);

        // Determine if error is acceptable
        let success = error_norm <= 1.0;

        // Compute timestep scale factor using accuracy order
        // Use minimum of propagating order (3) and embedded order (2)
        let order = 2;
        let mut timestep_scale = self.beta / error_norm.powf(1.0 / (order as f64 + 1.0));

        // Clip rescale factor to reasonable range [0.1, 10.0]
        timestep_scale = timestep_scale.clamp(0.1, 10.0);

        (success, error_norm, timestep_scale)
    }
}

impl Solver for RKBS32 {
    fn state(&self) -> &DVector<f64> {
        &self.state
    }

    fn state_mut(&mut self) -> &mut DVector<f64> {
        &mut self.state
    }

    fn set_state(&mut self, state: DVector<f64>) {
        self.state = state;
    }

    fn buffer(&mut self, _dt: f64) {
        if self.history.len() >= 2 {
            self.history.pop_back();
        }
        self.history.push_front(self.state.clone());
        self.stage = 0;
    }

    fn revert(&mut self) -> Result<(), SolverError> {
        self.state = self.history.pop_front().ok_or(SolverError::EmptyHistory)?;
        self.stage = 0;
        Ok(())
    }

    fn reset(&mut self) {
        self.state = self.initial.clone();
        self.history.clear();
        self.stage = 0;
    }

    fn order(&self) -> usize {
        3
    }

    fn stages(&self) -> usize {
        4
    }

    fn is_adaptive(&self) -> bool {
        true
    }

    fn is_explicit(&self) -> bool {
        true
    }
}

impl ExplicitSolver for RKBS32 {
    fn step<F>(&mut self, mut f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        let x0 = self
            .history
            .front()
            .expect("Must call buffer() before step()");

        // RKBS32 Butcher tableau
        // c (evaluation times) = [0, 1/2, 3/4, 1]
        let c = [0.0, 1.0 / 2.0, 3.0 / 4.0, 1.0];

        // Butcher tableau coefficients (a_ij)
        #[rustfmt::skip]
        let a: [&[f64]; 4] = [
            &[1.0/2.0],
            &[0.0, 3.0/4.0],
            &[2.0/9.0, 1.0/3.0, 4.0/9.0],
            &[2.0/9.0, 1.0/3.0, 4.0/9.0],
        ];

        // Evaluate slope at current stage
        self.slopes[self.stage] = f(&self.state, c[self.stage] * dt);

        // Compute next intermediate state or final state
        if self.stage < 3 {
            // Intermediate stages
            let mut slope_sum = DVector::zeros(x0.len());
            for (i, &coef) in a[self.stage].iter().enumerate() {
                slope_sum += coef * &self.slopes[i];
            }
            self.state = x0 + dt * slope_sum;
            self.stage += 1;

            SolverStepResult::default()
        } else {
            // Final stage - compute error estimate and timestep scale
            let (success, error_norm, scale) = self.error_controller(dt);
            self.stage = 0;

            SolverStepResult {
                success,
                error_norm,
                scale: Some(scale),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_rkbs32_properties() {
        let x0 = DVector::from_vec(vec![1.0]);
        let solver = RKBS32::new(x0);

        assert_eq!(solver.order(), 3);
        assert_eq!(solver.stages(), 4);
        assert!(solver.is_adaptive());
        assert!(solver.is_explicit());
    }

    #[test]
    fn test_rkbs32_exponential_decay() {
        // dx/dt = -x, x(0) = 1
        // Exact solution: x(t) = exp(-t)
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = RKBS32::new(x0);

        let dt = 0.1;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        for _ in 0..n_steps {
            solver.buffer(dt);

            // Perform all 4 stages
            for _ in 0..4 {
                solver.step(|x, _t| -x, dt);
            }
        }

        let exact = (-t_final).exp();
        // RKBS32 is 3rd order
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-4);
    }

    #[test]
    fn test_rkbs32_adaptive_step() {
        // Test that error controller returns valid scale factor
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = RKBS32::new(x0);

        solver.buffer(0.1);

        // Perform all 4 stages
        let mut result = SolverStepResult::default();
        for _ in 0..4 {
            result = solver.step(|x, _t| -x, 0.1);
        }

        // Should have a scale factor
        assert!(result.scale.is_some());
        let scale = result.scale.unwrap();

        // Scale should be in reasonable range [0.1, 10.0]
        assert!(scale >= 0.1 && scale <= 10.0);
    }
}
