//! ESDIRK32 - Four-stage, 3rd order ESDIRK with embedded 2nd order error estimate
//!
//! L-stable and stiffly accurate implicit Runge-Kutta method.

use nalgebra::DVector;
use std::collections::VecDeque;

use super::{ImplicitSolver, Solver, SolverError, SolverStepResult};

/// ESDIRK32 - Four-stage, 3(2) embedded ESDIRK method
///
/// Four-stage, 3rd order ESDIRK method with embedded 2nd order error
/// estimate. L-stable and stiffly accurate.
///
/// # Characteristics
/// - Order: 3 (propagating) / 2 (embedded)
/// - Stages: 4 (1 explicit, 3 implicit)
/// - Adaptive timestep
/// - L-stable, stiffly accurate
/// - Stage order 2 (γ = 1/2)
///
/// # Note
/// The cheapest adaptive implicit Runge-Kutta solver in this library,
/// yet remarkably robust. L-stability and stiff accuracy guarantee that
/// high-frequency parasitic modes are fully damped regardless of
/// timestep, and the optimal stage order of 2 (from γ = 1/2)
/// minimises order reduction on stiff problems. Three implicit stages
/// per step keeps the cost well below ESDIRK43 while still providing
/// adaptive step-size control.
///
/// # References
/// - Kennedy, C. A., & Carpenter, M. H. (2019). "Diagonally implicit
///   Runge-Kutta methods for stiff ODEs". Applied Numerical
///   Mathematics, 146, 221-244.
/// - Hairer, E., & Wanner, G. (1996). "Solving Ordinary Differential
///   Equations II: Stiff and Differential-Algebraic Problems". Springer
///   Series in Computational Mathematics, Vol. 14.
#[derive(Debug, Clone)]
pub struct ESDIRK32 {
    state: DVector<f64>,
    initial: DVector<f64>,
    history: VecDeque<DVector<f64>>,
    slopes: Vec<DVector<f64>>,
    stage: usize,
    tol_abs: f64,
    tol_rel: f64,
    beta: f64,
    max_iterations: usize,
    tolerance_fpi: f64,
}

impl ESDIRK32 {
    /// Create a new ESDIRK32 solver with the given initial state
    pub fn new(initial: DVector<f64>) -> Self {
        Self::with_tolerances(initial, 1e-8, 1e-4)
    }

    /// Create a new ESDIRK32 solver with custom tolerances
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
            beta: 0.9,
            max_iterations: 100,
            tolerance_fpi: 1e-10,
        }
    }

    /// Set fixed-point iteration tolerance
    pub fn with_fpi_tolerance(mut self, tol: f64) -> Self {
        self.tolerance_fpi = tol;
        self
    }

    /// Set maximum iterations for fixed-point iteration
    pub fn with_max_iterations(mut self, max_iter: usize) -> Self {
        self.max_iterations = max_iter;
        self
    }

    /// Get intermediate evaluation times (c coefficients)
    fn eval_stages(&self) -> [f64; 4] {
        [0.0, 1.0, 3.0 / 2.0, 1.0]
    }

    /// Compute error norm and timestep scale factor
    fn error_controller(&self, dt: f64) -> (bool, f64, f64) {
        // Coefficients for truncation error estimate
        let tr = [-1.0 / 9.0, -1.0 / 6.0, -2.0 / 9.0, 1.0 / 2.0];

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
        let order = 2; // min(m, n)
        let mut timestep_scale = self.beta / error_norm.powf(1.0 / (order as f64 + 1.0));

        // Clip rescale factor to reasonable range [0.1, 10.0]
        timestep_scale = timestep_scale.clamp(0.1, 10.0);

        (success, error_norm, timestep_scale)
    }
}

impl Solver for ESDIRK32 {
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
        false
    }
}

impl ImplicitSolver for ESDIRK32 {
    fn step<F>(&mut self, _f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        // For explicit first stage or intermediate stages, no error estimate
        if self.stage < 3 {
            return SolverStepResult::default();
        }

        // Final stage - compute error estimate and timestep scale
        let (success, error_norm, scale) = self.error_controller(dt);
        self.stage = 0;

        SolverStepResult {
            success,
            error_norm,
            scale: Some(scale),
        }
    }

    fn solve<F, J>(&mut self, mut f: F, _jac: Option<J>, dt: f64) -> Result<f64, SolverError>
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
        J: FnMut(&DVector<f64>, f64) -> nalgebra::DMatrix<f64>,
    {
        let x0 = self
            .history
            .front()
            .expect("Must call buffer() before solve()");

        let c = self.eval_stages();

        // Stage 0 is explicit (ESDIRK property)
        if self.stage == 0 {
            self.slopes[0] = f(&self.state, c[0] * dt);
            self.state = x0.clone();
            self.stage += 1;
            return Ok(0.0);
        }

        // Butcher tableau coefficients (lower triangular part)
        let bt: [&[f64]; 4] = [
            &[],                                             // Stage 0: explicit
            &[1.0 / 2.0, 1.0 / 2.0],                         // Stage 1
            &[5.0 / 8.0, 3.0 / 8.0, 1.0 / 2.0],              // Stage 2
            &[7.0 / 18.0, 1.0 / 3.0, -2.0 / 9.0, 1.0 / 2.0], // Stage 3
        ];

        // Diagonal entry (gamma)
        let gamma = 1.0 / 2.0;

        // Compute explicit part of the slope
        let mut slope_explicit = DVector::zeros(x0.len());
        for (i, &coef) in bt[self.stage].iter().enumerate().take(self.stage) {
            slope_explicit += coef * &self.slopes[i];
        }

        // Fixed-point iteration for implicit stage
        let mut residual = f64::INFINITY;
        for iter in 0..self.max_iterations {
            // Evaluate function at current state
            self.slopes[self.stage] = f(&self.state, c[self.stage] * dt);

            // Compute new state
            let x_new = x0 + dt * (&slope_explicit + gamma * &self.slopes[self.stage]);

            // Check convergence
            residual = (&x_new - &self.state).norm();
            self.state = x_new;

            if residual < self.tolerance_fpi {
                break;
            }

            if iter == self.max_iterations - 1 {
                return Err(SolverError::ConvergenceFailure(self.max_iterations));
            }
        }

        // Re-evaluate slope at converged state
        self.slopes[self.stage] = f(&self.state, c[self.stage] * dt);

        // Move to next stage
        self.stage += 1;

        Ok(residual)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_esdirk32_properties() {
        let x0 = DVector::from_vec(vec![1.0]);
        let solver = ESDIRK32::new(x0);

        assert_eq!(solver.order(), 3);
        assert_eq!(solver.stages(), 4);
        assert!(solver.is_adaptive());
        assert!(!solver.is_explicit());
    }

    #[test]
    fn test_esdirk32_exponential_decay() {
        // dx/dt = -x, x(0) = 1, exact: x(t) = exp(-t)
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = ESDIRK32::new(x0);

        let dt = 0.1;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..4 {
                let _ = solver.solve(
                    |x, _t| -x,
                    None::<fn(&DVector<f64>, f64) -> nalgebra::DMatrix<f64>>,
                    dt,
                );
            }
        }

        let exact = (-t_final).exp();
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-4);
    }
}
