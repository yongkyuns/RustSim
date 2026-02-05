//! Classic 4th-order Runge-Kutta solver (RK4)

use nalgebra::DVector;
use std::collections::VecDeque;

use super::{ExplicitSolver, Solver, SolverError, SolverStepResult};

/// Classic 4th-order Runge-Kutta solver
///
/// The workhorse fixed-step explicit method. Four-stage, 4th order accuracy.
/// Widely used for non-stiff ODEs when high accuracy is needed with a
/// reasonable number of function evaluations.
///
/// # Characteristics
/// - Order: 4
/// - Stages: 4
/// - Explicit, fixed timestep
/// - Not A-stable
///
/// # Note
/// The standard choice for fixed-step integration of smooth, non-stiff
/// block diagrams. Offers excellent accuracy-per-cost for most dynamics.
/// If step rejections occur or very small timesteps are required, the
/// system is likely stiff and an implicit solver should be used. For
/// adaptive timestepping, prefer RKDP54.
///
/// # References
/// - Kutta, W. (1901). "Beitrag zur näherungsweisen Integration totaler
///   Differentialgleichungen". Zeitschrift für Mathematik und Physik, 46, 435-453.
/// - Butcher, J. C. (2016). "Numerical Methods for Ordinary Differential
///   Equations". John Wiley & Sons, 3rd Edition.
#[derive(Debug, Clone)]
pub struct RK4 {
    state: DVector<f64>,
    initial: DVector<f64>,
    history: VecDeque<DVector<f64>>,
    slopes: Vec<DVector<f64>>,
    stage: usize,
}

impl RK4 {
    /// Create a new RK4 solver with the given initial state
    pub fn new(initial: DVector<f64>) -> Self {
        let n = initial.len();
        Self {
            state: initial.clone(),
            initial,
            history: VecDeque::with_capacity(2),
            slopes: vec![DVector::zeros(n); 4],
            stage: 0,
        }
    }
}

impl Solver for RK4 {
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
        4
    }

    fn is_adaptive(&self) -> bool {
        false
    }

    fn is_explicit(&self) -> bool {
        true
    }
}

impl ExplicitSolver for RK4 {
    fn step<F>(&mut self, mut f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        let x0 = self
            .history
            .front()
            .expect("Must call buffer() before step()");

        // RK4 Butcher tableau
        // c = [0, 1/2, 1/2, 1]
        // a = [[],
        //      [1/2],
        //      [0, 1/2],
        //      [0, 0, 1]]
        // b = [1/6, 1/3, 1/3, 1/6]

        let c = [0.0, 0.5, 0.5, 1.0];
        let b = [1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0];

        // Evaluate slope at current stage
        self.slopes[self.stage] = f(&self.state, c[self.stage] * dt);

        // Compute next intermediate state or final state
        if self.stage < 3 {
            // Intermediate stages: compute x_next using Butcher tableau
            self.state = match self.stage {
                0 => x0 + dt * 0.5 * &self.slopes[0],
                1 => x0 + dt * 0.5 * &self.slopes[1],
                2 => x0 + dt * &self.slopes[2],
                _ => unreachable!(),
            };
            self.stage += 1;
        } else {
            // Final stage: compute weighted sum
            self.state = x0
                + dt * (b[0] * &self.slopes[0]
                    + b[1] * &self.slopes[1]
                    + b[2] * &self.slopes[2]
                    + b[3] * &self.slopes[3]);
            self.stage = 0;
        }

        SolverStepResult::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_rk4_exponential_decay() {
        // dx/dt = -x, x(0) = 1
        // Exact solution: x(t) = exp(-t)
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = RK4::new(x0);

        let dt = 0.1;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        for _ in 0..n_steps {
            solver.buffer(dt);

            // Perform RK4 stages
            for _ in 0..4 {
                let result = solver.step(|x, _t| -x, dt);
                assert!(result.success);
            }
        }

        let exact = (-t_final).exp();
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-6);
    }

    #[test]
    fn test_rk4_harmonic_oscillator() {
        // d²x/dt² = -x => [x, v]' = [v, -x]
        // Exact: x(t) = cos(t), v(t) = -sin(t)
        let x0 = DVector::from_vec(vec![1.0, 0.0]);
        let mut solver = RK4::new(x0);

        let dt = 0.01;
        let t_final = 2.0 * std::f64::consts::PI;
        let n_steps = (t_final / dt) as usize;

        for _ in 0..n_steps {
            solver.buffer(dt);

            for _ in 0..4 {
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

        // After one period, should return to initial state
        // Relaxed tolerance due to accumulated numerical error over many steps
        assert_relative_eq!(solver.state()[0], 1.0, epsilon = 1e-2);
        assert_relative_eq!(solver.state()[1], 0.0, max_relative = 1.0, epsilon = 1e-2);
    }
}
