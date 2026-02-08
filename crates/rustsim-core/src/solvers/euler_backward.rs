//! Implicit backward Euler method for numerical integration

use nalgebra::{DMatrix, DVector};
use std::collections::VecDeque;

use super::{ImplicitSolver, NewtonAnderson, Solver, SolverError, SolverStepResult};

/// Implicit backward Euler method
///
/// First-order, A-stable and L-stable implicit integration method.
///
/// # Mathematical Form
/// ```text
/// x_{n+1} = x_n + h * f(x_{n+1}, t_{n+1})
/// ```
///
/// The implicit equation is solved iteratively using fixed-point iteration
/// with Anderson acceleration or Newton's method when Jacobian is provided.
///
/// # Characteristics
/// - Order: 1
/// - Stages: 1 (implicit)
/// - Fixed timestep
/// - A-stable, L-stable
///
/// # Note
/// Maximum stability at the cost of accuracy. L-stability fully damps parasitic
/// high-frequency modes, making this a safe fallback for very stiff systems
/// (e.g., high-gain PID loops or fast electrical dynamics coupled to slow
/// mechanical plant). Because each step requires solving a nonlinear equation,
/// the cost per step is higher than explicit methods.
///
/// # References
/// - Curtiss, C. F., & Hirschfelder, J. O. (1952). "Integration of stiff
///   equations". Proceedings of the National Academy of Sciences, 38(3), 235-243.
/// - Hairer, E., & Wanner, G. (1996). "Solving Ordinary Differential Equations II:
///   Stiff and Differential-Algebraic Problems". Springer Series in Computational
///   Mathematics, Vol. 14.
#[derive(Debug, Clone)]
pub struct EulerBackward {
    state: DVector<f64>,
    initial: DVector<f64>,
    history: VecDeque<DVector<f64>>,
    optimizer: NewtonAnderson,
    stage: usize,
}

impl EulerBackward {
    /// Create a new Euler Backward solver with the given initial state
    pub fn new(initial: DVector<f64>) -> Self {
        Self {
            state: initial.clone(),
            initial,
            history: VecDeque::with_capacity(2),
            optimizer: NewtonAnderson::default(),
            stage: 0,
        }
    }

    /// Get current stage
    pub fn current_stage(&self) -> usize {
        self.stage
    }

    /// Check if this is the first stage
    pub fn is_first_stage(&self) -> bool {
        self.stage == 0
    }

    /// Check if this is the last stage
    pub fn is_last_stage(&self) -> bool {
        self.stage == 0 // Only 1 stage for Euler
    }

    /// Stage iterator yielding evaluation times
    pub fn stages(&mut self, t: f64, dt: f64) -> impl Iterator<Item = f64> + '_ {
        self.stage = 0;
        let eval_stage = t + dt; // Evaluate at t+dt for backward Euler
        std::iter::once(eval_stage).inspect(move |_| {
            self.stage += 1;
        })
    }
}

impl Solver for EulerBackward {
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
        self.optimizer.reset();
    }

    fn revert(&mut self) -> Result<(), SolverError> {
        self.state = self.history.pop_front().ok_or(SolverError::EmptyHistory)?;
        Ok(())
    }

    fn reset(&mut self) {
        self.state = self.initial.clone();
        self.history.clear();
        self.stage = 0;
        self.optimizer.reset();
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
        false
    }
}

impl ImplicitSolver for EulerBackward {
    fn solve<F, J>(&mut self, mut f: F, mut jac: Option<J>, dt: f64) -> Result<f64, SolverError>
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
        J: FnMut(&DVector<f64>, f64) -> DMatrix<f64>,
    {
        // Get current state from history
        let x0 = self
            .history
            .front()
            .ok_or(SolverError::EmptyHistory)?
            .clone();

        // Evaluate f at current state
        let f_val = f(&self.state, 0.0);

        // Update the fixed point equation: x = x0 + dt * f(x, t+dt)
        let g = x0 + dt * f_val;

        // Compute Jacobian if available
        let j_opt = jac.as_mut().map(|j_fn| dt * j_fn(&self.state, 0.0));

        // Use the optimizer to solve for x
        let (x_new, err) = self.optimizer.step(&self.state, &g, j_opt.as_ref());
        self.state = x_new;

        Ok(err)
    }

    fn step<F>(&mut self, _f: F, _dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        // For single-stage backward Euler, the step is already done in solve()
        SolverStepResult::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const MAX_ITER: usize = 100;
    const TOL_FPI: f64 = 1e-10;

    fn integrate_implicit<F, J>(
        solver: &mut EulerBackward,
        f: F,
        jac: J,
        t_start: f64,
        t_end: f64,
        dt: f64,
    ) -> (Vec<f64>, Vec<DVector<f64>>)
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64> + Clone,
        J: FnMut(&DVector<f64>, f64) -> DMatrix<f64> + Clone,
    {
        let mut times = vec![t_start];
        let mut states = vec![solver.state().clone()];
        let mut t = t_start;

        while t < t_end - 1e-10 {
            solver.buffer(dt);

            // Iterate through stages (only 1 for EUB)
            for _t_eval in solver.stages(t, dt) {
                // Fixed-point iteration to solve implicit equation
                for _ in 0..MAX_ITER {
                    let err = solver.solve(f.clone(), Some(jac.clone()), dt).unwrap();

                    if err < TOL_FPI {
                        break;
                    }
                }

                // Perform explicit step (no-op for EUB)
                solver.step(f.clone(), dt);
            }

            t += dt;
            times.push(t);
            states.push(solver.state().clone());
        }

        (times, states)
    }

    #[test]
    fn test_eub_init() {
        let solver = EulerBackward::new(DVector::from_vec(vec![1.0]));
        assert_eq!(solver.order(), 1);
        assert_eq!(solver.stages(), 1);
        assert!(!solver.is_adaptive());
        assert!(!solver.is_explicit());
        assert_eq!(solver.state()[0], 1.0);
    }

    #[test]
    fn test_eub_exponential_decay() {
        // dx/dt = -x, x(0) = 1
        // Exact solution: x(t) = exp(-t)
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = EulerBackward::new(x0);

        let f = |x: &DVector<f64>, _t: f64| -x;
        let jac = |_x: &DVector<f64>, _t: f64| DMatrix::from_element(1, 1, -1.0);

        let dt = 0.1;
        let t_final = 1.0;

        let (_times, states) = integrate_implicit(&mut solver, f, jac, 0.0, t_final, dt);

        let exact = (-t_final).exp();
        let computed = states.last().unwrap()[0];

        // Backward Euler is first order
        assert_relative_eq!(computed, exact, epsilon = 0.05);
    }

    #[test]
    fn test_eub_stiff_decay() {
        // Stiff problem: dx/dt = -1000*x, x(0) = 1
        // Exact solution: x(t) = exp(-1000*t)
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = EulerBackward::new(x0);

        let lambda = -1000.0;
        let f = |x: &DVector<f64>, _t: f64| lambda * x;
        let jac = |_x: &DVector<f64>, _t: f64| DMatrix::from_element(1, 1, lambda);

        // Large timestep that would be unstable for explicit methods
        let dt = 0.01;
        let t_final = 0.1;

        let (_times, states) = integrate_implicit(&mut solver, f, jac, 0.0, t_final, dt);

        let exact = (lambda * t_final).exp();
        let computed = states.last().unwrap()[0];

        // Should remain stable even with stiff problem
        assert!(computed.is_finite());
        assert_relative_eq!(computed, exact, epsilon = 0.1);
    }

    #[test]
    fn test_eub_buffer_and_revert() {
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = EulerBackward::new(x0);

        // Change state
        solver.set_state(DVector::from_vec(vec![2.0]));

        // Buffer current state
        solver.buffer(0.1);

        // Change state again
        solver.set_state(DVector::from_vec(vec![3.0]));

        // Revert should restore buffered state
        solver.revert().unwrap();
        assert_eq!(solver.state()[0], 2.0);
    }

    #[test]
    fn test_eub_reset() {
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = EulerBackward::new(x0);

        // Change state
        solver.set_state(DVector::from_vec(vec![5.0]));

        // Reset should restore initial state
        solver.reset();
        assert_eq!(solver.state()[0], 1.0);
    }

    #[test]
    fn test_eub_stages_iterator() {
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = EulerBackward::new(x0);

        let t = 0.0;
        let dt = 1.0;

        let stages: Vec<f64> = solver.stages(t, dt).collect();

        assert_eq!(stages.len(), 1);
        assert_relative_eq!(stages[0], 1.0, epsilon = 1e-10); // t + dt
    }
}
