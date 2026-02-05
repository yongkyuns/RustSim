//! ESDIRK4 - Six-stage, 4th order ESDIRK method
//!
//! L-stable and stiffly accurate implicit Runge-Kutta method. Fixed timestep only.

use nalgebra::DVector;
use std::collections::VecDeque;

use super::{ImplicitSolver, Solver, SolverError, SolverStepResult};

/// ESDIRK4 - Six-stage, 4th order ESDIRK method
///
/// Six-stage, 4th order ESDIRK method. L-stable and stiffly accurate.
/// No embedded error estimator; fixed timestep only.
///
/// # Characteristics
/// - Order: 4
/// - Stages: 6 (1 explicit, 5 implicit)
/// - Fixed timestep
/// - L-stable, stiffly accurate
/// - Stage order 2
///
/// # Note
/// Provides 4th order accuracy on stiff systems when the timestep is
/// predetermined (e.g. real-time or hardware-in-the-loop contexts). The
/// explicit first stage reuses the last function evaluation from the
/// previous step, saving one implicit solve per step compared to a fully
/// implicit DIRK. L-stability and stiff accuracy ensure full damping of
/// parasitic high-frequency modes. For adaptive stepping, use ESDIRK43
/// which adds an embedded error estimator at the same stage count.
///
/// # References
/// - Kennedy, C. A., & Carpenter, M. H. (2016). "Diagonally implicit
///   Runge-Kutta methods for ordinary differential equations. A review".
///   NASA/TM-2016-219173.
/// - Hairer, E., & Wanner, G. (1996). "Solving Ordinary Differential
///   Equations II: Stiff and Differential-Algebraic Problems". Springer
///   Series in Computational Mathematics, Vol. 14.
#[derive(Debug, Clone)]
pub struct ESDIRK4 {
    state: DVector<f64>,
    initial: DVector<f64>,
    history: VecDeque<DVector<f64>>,
    slopes: Vec<DVector<f64>>,
    stage: usize,
    max_iterations: usize,
    tolerance_fpi: f64,
}

impl ESDIRK4 {
    /// Create a new ESDIRK4 solver with the given initial state
    pub fn new(initial: DVector<f64>) -> Self {
        let n = initial.len();
        Self {
            state: initial.clone(),
            initial,
            history: VecDeque::with_capacity(2),
            slopes: vec![DVector::zeros(n); 6],
            stage: 0,
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
    fn eval_stages(&self) -> [f64; 6] {
        [0.0, 1.0 / 2.0, 1.0 / 6.0, 37.0 / 40.0, 1.0 / 2.0, 1.0]
    }
}

impl Solver for ESDIRK4 {
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
        4
    }

    fn stages(&self) -> usize {
        6
    }

    fn is_adaptive(&self) -> bool {
        false
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

impl ImplicitSolver for ESDIRK4 {
    fn step<F>(&mut self, mut f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        // No error estimate for fixed-step solver
        SolverStepResult::default()
    }

    fn solve<F, J>(&mut self, mut f: F, mut jac: Option<J>, dt: f64) -> Result<f64, SolverError>
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

        // Butcher tableau coefficients
        let bt: [&[f64]; 6] = [
            &[], // Stage 0: explicit
            &[1.0 / 4.0, 1.0 / 4.0],
            &[-1.0 / 36.0, -1.0 / 18.0, 1.0 / 4.0],
            &[
                -21283.0 / 32000.0,
                -5143.0 / 64000.0,
                90909.0 / 64000.0,
                1.0 / 4.0,
            ],
            &[
                46010759.0 / 749250000.0,
                -737693.0 / 40500000.0,
                10931269.0 / 45500000.0,
                -1140071.0 / 34090875.0,
                1.0 / 4.0,
            ],
            &[
                89.0 / 444.0,
                89.0 / 804756.0,
                -27.0 / 364.0,
                -20000.0 / 171717.0,
                843750.0 / 1140071.0,
                1.0 / 4.0,
            ],
        ];

        // Diagonal entry (gamma)
        let gamma = 1.0 / 4.0;

        // Compute explicit part of the slope
        let mut slope_explicit = DVector::zeros(x0.len());
        for (i, &coef) in bt[self.stage].iter().enumerate().take(self.stage) {
            slope_explicit += coef * &self.slopes[i];
        }

        // Fixed-point iteration for implicit stage
        let mut residual = f64::INFINITY;
        for iter in 0..self.max_iterations {
            self.slopes[self.stage] = f(&self.state, c[self.stage] * dt);

            let x_new = x0 + dt * (&slope_explicit + gamma * &self.slopes[self.stage]);

            residual = (&x_new - &self.state).norm();
            self.state = x_new;

            if residual < self.tolerance_fpi {
                break;
            }

            if iter == self.max_iterations - 1 {
                return Err(SolverError::ConvergenceFailure(self.max_iterations));
            }
        }

        self.slopes[self.stage] = f(&self.state, c[self.stage] * dt);
        self.stage += 1;

        Ok(residual)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_esdirk4_properties() {
        let x0 = DVector::from_vec(vec![1.0]);
        let solver = ESDIRK4::new(x0);

        assert_eq!(solver.order(), 4);
        assert_eq!(solver.stages(), 6);
        assert!(!solver.is_adaptive());
        assert!(!solver.is_explicit());
    }
}
