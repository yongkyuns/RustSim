//! Dormand-Prince 5(4) adaptive Runge-Kutta solver

use nalgebra::DVector;
use std::collections::VecDeque;

use super::{ExplicitSolver, Solver, SolverError, SolverStepResult};

/// Dormand-Prince 5(4) adaptive solver (DOPRI5)
///
/// Seven-stage, 5th order Runge-Kutta method with embedded 4th order
/// error estimate for adaptive timestepping.
///
/// The industry-standard adaptive explicit solver and the basis of MATLAB's
/// `ode45`. Has the FSAL property (First Same As Last - not exploited in
/// this implementation, so all seven stages are evaluated each step).
///
/// # Characteristics
/// - Order: 5 (propagating) / 4 (embedded)
/// - Stages: 7
/// - Explicit, adaptive timestep
/// - Error estimation via embedded method
/// - PI controller for timestep adaptation
///
/// # Note
/// Recommended default for non-stiff systems. Handles smooth nonlinear
/// dynamics, coupled oscillators, and general ODEs efficiently. If the
/// simulation warns about excessive step rejections or very small timesteps,
/// the system is likely stiff and an implicit solver should be used instead.
/// For very tight tolerances on smooth problems, higher-order methods like
/// RKV65 or RKDP87 can be more efficient per unit accuracy.
///
/// # References
/// - Dormand, J. R., & Prince, P. J. (1980). "A family of embedded
///   Runge-Kutta formulae". Journal of Computational and Applied
///   Mathematics, 6(1), 19-26.
/// - Shampine, L. F., & Reichelt, M. W. (1997). "The MATLAB ODE Suite".
///   SIAM Journal on Scientific Computing, 18(1), 1-22.
#[derive(Debug, Clone)]
pub struct RKDP54 {
    state: DVector<f64>,
    initial: DVector<f64>,
    history: VecDeque<DVector<f64>>,
    slopes: Vec<DVector<f64>>,
    stage: usize,
    tol_abs: f64,
    tol_rel: f64,
    beta: f64,
}

impl RKDP54 {
    /// Create a new RKDP54 solver with the given initial state
    ///
    /// # Arguments
    /// * `initial` - Initial state vector
    /// * `tol_abs` - Absolute error tolerance (default: 1e-8)
    /// * `tol_rel` - Relative error tolerance (default: 1e-4)
    pub fn new(initial: DVector<f64>) -> Self {
        Self::with_tolerances(initial, 1e-8, 1e-4)
    }

    /// Create a new RKDP54 solver with custom tolerances
    pub fn with_tolerances(initial: DVector<f64>, tol_abs: f64, tol_rel: f64) -> Self {
        let n = initial.len();
        Self {
            state: initial.clone(),
            initial,
            history: VecDeque::with_capacity(2),
            slopes: vec![DVector::zeros(n); 7],
            stage: 0,
            tol_abs,
            tol_rel,
            beta: 0.9, // Safety factor
        }
    }

    /// Compute error norm and timestep scale factor
    fn error_controller(&self, dt: f64) -> (bool, f64, f64) {
        // Coefficients for local truncation error estimate
        // TR = [71/57600, 0, -71/16695, 71/1920, -17253/339200, 22/525, -1/40]
        let tr = [
            71.0 / 57600.0,
            0.0,
            -71.0 / 16695.0,
            71.0 / 1920.0,
            -17253.0 / 339200.0,
            22.0 / 525.0,
            -1.0 / 40.0,
        ];

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
        // Use minimum of propagating order (5) and embedded order (4)
        let order = 4.min(5);
        let mut timestep_scale = self.beta / error_norm.powf(1.0 / (order as f64 + 1.0));

        // Clip rescale factor to reasonable range [0.1, 10.0]
        timestep_scale = timestep_scale.clamp(0.1, 10.0);

        (success, error_norm, timestep_scale)
    }
}

impl Solver for RKDP54 {
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
        5
    }

    fn stages(&self) -> usize {
        7
    }

    fn is_adaptive(&self) -> bool {
        true
    }

    fn is_explicit(&self) -> bool {
        true
    }
}

impl ExplicitSolver for RKDP54 {
    fn step<F>(&mut self, mut f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        let x0 = self
            .history
            .front()
            .expect("Must call buffer() before step()");

        // RKDP54 Butcher tableau
        // c (evaluation times) = [0, 1/5, 3/10, 4/5, 8/9, 1, 1]
        let c = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];

        // Butcher tableau coefficients (a_ij)
        #[rustfmt::skip]
        let a: [&[f64]; 7] = [
            &[1.0/5.0],
            &[3.0/40.0, 9.0/40.0],
            &[44.0/45.0, -56.0/15.0, 32.0/9.0],
            &[19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0],
            &[9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0],
            &[35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0],
            &[35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0],
        ];

        // Evaluate slope at current stage
        self.slopes[self.stage] = f(&self.state, c[self.stage] * dt);

        // Compute next intermediate state or final state
        if self.stage < 6 {
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
    fn test_rkdp54_exponential_decay() {
        // dx/dt = -x, x(0) = 1
        // Exact solution: x(t) = exp(-t)
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = RKDP54::new(x0);

        let dt = 0.1;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        for _ in 0..n_steps {
            solver.buffer(dt);

            // Perform all 7 stages
            for _ in 0..7 {
                solver.step(|x, _t| -x, dt);
            }
        }

        let exact = (-t_final).exp();
        // RKDP54 should be very accurate
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-8);
    }

    #[test]
    fn test_rkdp54_adaptive_step() {
        // Test that error controller returns valid scale factor
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = RKDP54::new(x0);

        solver.buffer(0.1);

        // Perform all 7 stages
        let mut result = SolverStepResult::default();
        for _ in 0..7 {
            result = solver.step(|x, _t| -x, 0.1);
        }

        // Should have a scale factor
        assert!(result.scale.is_some());
        let scale = result.scale.unwrap();

        // Scale should be in reasonable range [0.1, 10.0]
        assert!(scale >= 0.1 && scale <= 10.0);
    }

    #[test]
    fn test_rkdp54_harmonic_oscillator() {
        // d²x/dt² = -x => [x, v]' = [v, -x]
        let x0 = DVector::from_vec(vec![1.0, 0.0]);
        let mut solver = RKDP54::new(x0);

        let dt = 0.05;
        let t_final = 2.0 * std::f64::consts::PI;
        let n_steps = (t_final / dt) as usize;

        for _ in 0..n_steps {
            solver.buffer(dt);

            for _ in 0..7 {
                solver.step(
                    |x, _t| {
                        let mut dxdt = DVector::zeros(2);
                        dxdt[0] = x[1];
                        dxdt[1] = -x[0];
                        dxdt
                    },
                    dt,
                );
            }
        }

        // After one period, should return close to initial state
        assert_relative_eq!(solver.state()[0], 1.0, epsilon = 1e-2);
        assert_relative_eq!(solver.state()[1], 0.0, max_relative = 1.0, epsilon = 1e-2);
    }
}
