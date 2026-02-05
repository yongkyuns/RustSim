//! Forward Euler method for numerical integration

use nalgebra::DVector;
use std::collections::VecDeque;

use super::{ExplicitSolver, Solver, SolverError, SolverStepResult};

/// Explicit forward Euler method
///
/// First-order, single-stage explicit integration method.
///
/// # Mathematical Form
/// ```text
/// x_{n+1} = x_n + h * f(x_n, t_n)
/// ```
///
/// # Characteristics
/// - Order: 1
/// - Stages: 1
/// - Explicit, fixed timestep
/// - Not A-stable
///
/// # Note
/// The cheapest solver per step but also the least accurate. Its small stability
/// region requires very small timesteps for moderately dynamic systems, which
/// usually makes higher-order methods more efficient overall. Prefer RK4 for
/// fixed-step or RKDP54 for adaptive integration of non-stiff systems. Only
/// practical when computational cost per step must be absolute minimum and
/// accuracy is secondary.
///
/// # References
/// - Hairer, E., Nørsett, S. P., & Wanner, G. (1993). "Solving Ordinary
///   Differential Equations I: Nonstiff Problems". Springer Series in
///   Computational Mathematics, Vol. 8.
/// - Butcher, J. C. (2016). "Numerical Methods for Ordinary Differential
///   Equations". John Wiley & Sons, 3rd Edition.
#[derive(Debug, Clone)]
pub struct Euler {
    state: DVector<f64>,
    initial: DVector<f64>,
    history: VecDeque<DVector<f64>>,
}

impl Euler {
    /// Create a new Euler solver with the given initial state
    pub fn new(initial: DVector<f64>) -> Self {
        Self {
            state: initial.clone(),
            initial,
            history: VecDeque::with_capacity(2),
        }
    }
}

impl Solver for Euler {
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
    }

    fn revert(&mut self) -> Result<(), SolverError> {
        self.state = self.history.pop_front().ok_or(SolverError::EmptyHistory)?;
        Ok(())
    }

    fn reset(&mut self) {
        self.state = self.initial.clone();
        self.history.clear();
    }

    fn order(&self) -> usize {
        1
    }

    fn stages(&self) -> usize {
        1
    }

    fn is_adaptive(&self) -> bool {
        false
    }

    fn is_explicit(&self) -> bool {
        true
    }
}

impl ExplicitSolver for Euler {
    fn step<F>(&mut self, mut f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        let x0 = self
            .history
            .front()
            .expect("Must call buffer() before step()");

        // Simple forward Euler: x_{n+1} = x_n + dt * f(x_n, t_n)
        self.state = x0 + dt * f(&self.state, 0.0);

        SolverStepResult::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_euler_exponential_decay() {
        // dx/dt = -x, x(0) = 1
        // Exact solution: x(t) = exp(-t)
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = Euler::new(x0);

        let dt = 0.01;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        for _ in 0..n_steps {
            solver.buffer(dt);
            let result = solver.step(|x, _t| -x, dt);
            assert!(result.success);
        }

        let exact = (-t_final).exp();
        // Euler has larger error than RK4
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-2);
    }

    #[test]
    fn test_euler_linear_growth() {
        // dx/dt = 1, x(0) = 0
        // Exact solution: x(t) = t
        let x0 = DVector::from_vec(vec![0.0]);
        let mut solver = Euler::new(x0);

        let dt = 0.1;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        for _ in 0..n_steps {
            solver.buffer(dt);
            solver.step(|_x, _t| DVector::from_vec(vec![1.0]), dt);
        }

        // For constant derivative, Euler should be exact
        assert_relative_eq!(solver.state()[0], t_final, epsilon = 1e-10);
    }
}
