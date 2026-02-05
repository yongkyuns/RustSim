//! Verner 6(5) adaptive solver

use nalgebra::DVector;
use std::collections::VecDeque;

use super::{ExplicitSolver, Solver, SolverError, SolverStepResult};

/// Verner 6(5) "most robust" pair adaptive solver
///
/// Nine stages, 6th order with embedded 5th order error estimate.
///
/// # Characteristics
/// - Order: 6 (propagating) / 5 (embedded)
/// - Stages: 9
/// - Explicit, adaptive timestep
///
/// # Note
/// Fills the gap between 5th order pairs (RKDP54) and the expensive 8th
/// order RKDP87. The extra stages pay off when the dynamics are smooth
/// and tolerances are tight (roughly 1e-8 or below), because the
/// higher order allows much larger steps. For tolerances in the
/// 1e-4 to 1e-6 range, RKDP54 is usually cheaper overall due to fewer stages.
///
/// # References
/// - Verner, J. H. (2010). "Numerically optimal Runge-Kutta pairs
///   with interpolants". Numerical Algorithms, 53(2-3), 383-396.
/// - Hairer, E., Nørsett, S. P., & Wanner, G. (1993). "Solving
///   Ordinary Differential Equations I: Nonstiff Problems". Springer
///   Series in Computational Mathematics, Vol. 8.
#[derive(Debug, Clone)]
pub struct RKV65 {
    state: DVector<f64>,
    initial: DVector<f64>,
    history: VecDeque<DVector<f64>>,
    slopes: Vec<DVector<f64>>,
    stage: usize,
    tol_abs: f64,
    tol_rel: f64,
    beta: f64,
}

impl RKV65 {
    /// Create a new RKV65 solver with the given initial state
    pub fn new(initial: DVector<f64>) -> Self {
        Self::with_tolerances(initial, 1e-8, 1e-4)
    }

    /// Create a new RKV65 solver with custom tolerances
    pub fn with_tolerances(initial: DVector<f64>, tol_abs: f64, tol_rel: f64) -> Self {
        let n = initial.len();
        Self {
            state: initial.clone(),
            initial,
            history: VecDeque::with_capacity(2),
            slopes: vec![DVector::zeros(n); 9],
            stage: 0,
            tol_abs,
            tol_rel,
            beta: 0.9,
        }
    }

    fn error_controller(&self, dt: f64) -> (bool, f64, f64) {
        // Compute coefficients for truncation error
        // _A1 = [11/144, 0, 0, 256/693, 0, 125/504, 125/528, 5/72, 0]
        let a1 = [
            11.0 / 144.0,
            0.0,
            0.0,
            256.0 / 693.0,
            0.0,
            125.0 / 504.0,
            125.0 / 528.0,
            5.0 / 72.0,
            0.0,
        ];
        // _A2 = [28/477, 0, 0, 212/441, -312500/366177, 2125/1764, 0, -2105/35532, 2995/17766]
        let a2 = [
            28.0 / 477.0,
            0.0,
            0.0,
            212.0 / 441.0,
            -312500.0 / 366177.0,
            2125.0 / 1764.0,
            0.0,
            -2105.0 / 35532.0,
            2995.0 / 17766.0,
        ];

        // TR = [a-b for a, b in zip(_A1, _A2)]
        let tr: Vec<f64> = a1
            .iter()
            .zip(a2.iter())
            .map(|(a, b)| a - b)
            .collect();

        let mut error_slope = DVector::zeros(self.state.len());
        for (i, &coef) in tr.iter().enumerate() {
            error_slope += coef * &self.slopes[i];
        }

        let scale = self.state.map(|x| self.tol_abs + self.tol_rel * x.abs());
        let scaled_error = (dt * &error_slope).component_div(&scale).map(|e| e.abs());
        let error_norm = scaled_error.max().max(1e-16);
        let success = error_norm <= 1.0;

        let order = 5.min(6);
        let mut timestep_scale = self.beta / error_norm.powf(1.0 / (order as f64 + 1.0));
        timestep_scale = timestep_scale.clamp(0.1, 10.0);

        (success, error_norm, timestep_scale)
    }
}

impl Solver for RKV65 {
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
        6
    }

    fn stages(&self) -> usize {
        9
    }

    fn is_adaptive(&self) -> bool {
        true
    }

    fn is_explicit(&self) -> bool {
        true
    }
}

impl ExplicitSolver for RKV65 {
    fn step<F>(&mut self, mut f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        let x0 = self
            .history
            .front()
            .expect("Must call buffer() before step()");

        // c = [0, 9/50, 1/6, 1/4, 53/100, 3/5, 4/5, 1, 1]
        let c = [
            0.0,
            9.0 / 50.0,
            1.0 / 6.0,
            1.0 / 4.0,
            53.0 / 100.0,
            3.0 / 5.0,
            4.0 / 5.0,
            1.0,
            1.0,
        ];

        #[rustfmt::skip]
        let a: [&[f64]; 9] = [
            &[9.0/50.0],
            &[29.0/324.0, 25.0/324.0],
            &[1.0/16.0, 0.0, 3.0/16.0],
            &[79129.0/250000.0, 0.0, -261237.0/250000.0, 19663.0/15625.0],
            &[1336883.0/4909125.0, 0.0, -25476.0/30875.0, 194159.0/185250.0, 8225.0/78546.0],
            &[-2459386.0/14727375.0, 0.0, 19504.0/30875.0, 2377474.0/13615875.0, -6157250.0/5773131.0, 902.0/735.0],
            &[2699.0/7410.0, 0.0, -252.0/1235.0, -1393253.0/3993990.0, 236875.0/72618.0, -135.0/49.0, 15.0/22.0],
            &[11.0/144.0, 0.0, 0.0, 256.0/693.0, 0.0, 125.0/504.0, 125.0/528.0, 5.0/72.0],
            &[11.0/144.0, 0.0, 0.0, 256.0/693.0, 0.0, 125.0/504.0, 125.0/528.0, 5.0/72.0],
        ];

        self.slopes[self.stage] = f(&self.state, c[self.stage] * dt);

        if self.stage < 8 {
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
    fn test_rkv65_properties() {
        let x0 = DVector::from_vec(vec![1.0]);
        let solver = RKV65::new(x0);
        assert_eq!(solver.order(), 6);
        assert_eq!(solver.stages(), 9);
        assert!(solver.is_adaptive());
        assert!(solver.is_explicit());
    }

    #[test]
    fn test_rkv65_exponential_decay() {
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = RKV65::new(x0);
        let dt = 0.1;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..9 {
                solver.step(|x, _t| -x, dt);
            }
        }

        let exact = (-t_final).exp();
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-9);
    }
}
