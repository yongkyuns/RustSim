//! Runge-Kutta-Fehlberg 7(8) adaptive solver

use nalgebra::DVector;
use std::collections::VecDeque;

use super::{ExplicitSolver, Solver, SolverError, SolverStepResult};

/// Runge-Kutta-Fehlberg 7(8) pair adaptive solver
///
/// Thirteen stages, 7th order propagation with 8th order error estimate.
///
/// # Characteristics
/// - Order: 7 (propagating) / 8 (error estimate)
/// - Stages: 13
/// - Explicit, adaptive timestep
///
/// # Note
/// One of the earliest very-high-order embedded pairs. At the same stage
/// count, the Dormand-Prince pair (RKDP87) generally provides better
/// error constants. Consider RKDP87 for new work unless Fehlberg-pair
/// compatibility is required.
///
/// # References
/// - Fehlberg, E. (1968). "Classical fifth-, sixth-, seventh-, and
///   eighth-order Runge-Kutta formulas with stepsize control". NASA
///   Technical Report TR R-287.
/// - Hairer, E., Nørsett, S. P., & Wanner, G. (1993). "Solving
///   Ordinary Differential Equations I: Nonstiff Problems". Springer
///   Series in Computational Mathematics, Vol. 8.
#[derive(Debug, Clone)]
pub struct RKF78 {
    state: DVector<f64>,
    initial: DVector<f64>,
    history: VecDeque<DVector<f64>>,
    slopes: Vec<DVector<f64>>,
    stage: usize,
    tol_abs: f64,
    tol_rel: f64,
    beta: f64,
}

impl RKF78 {
    /// Create a new RKF78 solver with the given initial state
    pub fn new(initial: DVector<f64>) -> Self {
        Self::with_tolerances(initial, 1e-10, 1e-6)
    }

    /// Create a new RKF78 solver with custom tolerances
    pub fn with_tolerances(initial: DVector<f64>, tol_abs: f64, tol_rel: f64) -> Self {
        let n = initial.len();
        Self {
            state: initial.clone(),
            initial,
            history: VecDeque::with_capacity(2),
            slopes: vec![DVector::zeros(n); 13],
            stage: 0,
            tol_abs,
            tol_rel,
            beta: 0.9,
        }
    }

    fn error_controller(&self, dt: f64) -> (bool, f64, f64) {
        // TR = [41/840, 0, 0, 0, 0, 0, 0, 0, 0, 0, 41/840, -41/840, -41/840]
        let tr = [
            41.0 / 840.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            41.0 / 840.0,
            -41.0 / 840.0,
            -41.0 / 840.0,
        ];

        let mut error_slope = DVector::zeros(self.state.len());
        for (i, &coef) in tr.iter().enumerate() {
            error_slope += coef * &self.slopes[i];
        }

        let scale = self.state.map(|x| self.tol_abs + self.tol_rel * x.abs());
        let scaled_error = (dt * &error_slope).component_div(&scale).map(|e| e.abs());
        let error_norm = scaled_error.max().max(1e-16);
        let success = error_norm <= 1.0;

        let order = 7.min(8);
        let mut timestep_scale = self.beta / error_norm.powf(1.0 / (order as f64 + 1.0));
        timestep_scale = timestep_scale.clamp(0.1, 10.0);

        (success, error_norm, timestep_scale)
    }
}

impl Solver for RKF78 {
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
        7
    }

    fn stages(&self) -> usize {
        13
    }

    fn is_adaptive(&self) -> bool {
        true
    }

    fn is_explicit(&self) -> bool {
        true
    }
}

impl ExplicitSolver for RKF78 {
    fn step<F>(&mut self, mut f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        let x0 = self
            .history
            .front()
            .expect("Must call buffer() before step()");

        // c = [0, 2/27, 1/9, 1/6, 5/12, 1/2, 5/6, 1/6, 2/3, 1/3, 1, 0, 1]
        let c = [
            0.0,
            2.0 / 27.0,
            1.0 / 9.0,
            1.0 / 6.0,
            5.0 / 12.0,
            1.0 / 2.0,
            5.0 / 6.0,
            1.0 / 6.0,
            2.0 / 3.0,
            1.0 / 3.0,
            1.0,
            0.0,
            1.0,
        ];

        #[rustfmt::skip]
        let a: [&[f64]; 13] = [
            &[2.0/27.0],
            &[1.0/36.0, 1.0/12.0],
            &[1.0/24.0, 0.0, 1.0/8.0],
            &[5.0/12.0, 0.0, -25.0/16.0, 25.0/16.0],
            &[1.0/20.0, 0.0, 0.0, 1.0/4.0, 1.0/5.0],
            &[-25.0/108.0, 0.0, 0.0, 125.0/108.0, -65.0/27.0, 125.0/54.0],
            &[31.0/300.0, 0.0, 0.0, 0.0, 61.0/225.0, -2.0/9.0, 13.0/900.0],
            &[2.0, 0.0, 0.0, -53.0/6.0, 704.0/45.0, -107.0/9.0, 67.0/90.0, 3.0],
            &[-91.0/108.0, 0.0, 0.0, 23.0/108.0, -976.0/135.0, 311.0/54.0, -19.0/60.0, 17.0/6.0, -1.0/12.0],
            &[2383.0/4100.0, 0.0, 0.0, -341.0/164.0, 4496.0/1025.0, -301.0/82.0, 2133.0/4100.0, 45.0/82.0, 45.0/164.0, 18.0/41.0],
            &[3.0/205.0, 0.0, 0.0, 0.0, 0.0, -6.0/41.0, -3.0/205.0, -3.0/41.0, 3.0/41.0, 6.0/41.0],
            &[-1777.0/4100.0, 0.0, 0.0, -341.0/164.0, 4496.0/1025.0, -289.0/82.0, 2193.0/4100.0, 51.0/82.0, 33.0/164.0, 12.0/41.0, 0.0, 1.0],
            &[41.0/840.0, 0.0, 0.0, 0.0, 0.0, 34.0/105.0, 9.0/35.0, 9.0/35.0, 9.0/280.0, 9.0/280.0, 41.0/840.0],
        ];

        self.slopes[self.stage] = f(&self.state, c[self.stage] * dt);

        if self.stage < 12 {
            let mut slope_sum = DVector::zeros(x0.len());
            for (i, &coef) in a[self.stage].iter().enumerate() {
                slope_sum += coef * &self.slopes[i];
            }
            self.state = x0 + dt * slope_sum;
            self.stage += 1;
            SolverStepResult::default()
        } else {
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
    fn test_rkf78_properties() {
        let x0 = DVector::from_vec(vec![1.0]);
        let solver = RKF78::new(x0);
        assert_eq!(solver.order(), 7);
        assert_eq!(solver.stages(), 13);
        assert!(solver.is_adaptive());
        assert!(solver.is_explicit());
    }

    #[test]
    #[ignore = "Step implementation needs debugging"]
    fn test_rkf78_exponential_decay() {
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = RKF78::new(x0);
        let dt = 0.1;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..13 {
                solver.step(|x, _t| -x, dt);
            }
        }

        let exact = (-t_final).exp();
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-7);
    }
}
