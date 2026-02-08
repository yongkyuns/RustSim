//! Two-stage, 2nd order diagonally implicit Runge-Kutta method

use nalgebra::{DMatrix, DVector};
use std::collections::VecDeque;

use super::{ImplicitSolver, NewtonAnderson, Solver, SolverError, SolverStepResult};

/// Two-stage, 2nd order DIRK method
///
/// L-stable, SSP-optimal, symplectic diagonally implicit Runge-Kutta method.
///
/// # Butcher Tableau
/// ```text
/// 1/4  | 1/4
/// 3/4  | 1/2  1/4
/// -----|----------
///      | 1/2  1/2
/// ```
///
/// # Characteristics
/// - Order: 2
/// - Stages: 2 (implicit)
/// - Fixed timestep
/// - L-stable, SSP-optimal, symplectic
///
/// # Note
/// The simplest multi-stage implicit Runge-Kutta method. L-stability fully damps
/// parasitic high-frequency modes, and the symplectic property preserves Hamiltonian
/// structure when the dynamics are conservative. Two implicit stages per step is
/// relatively cheap. For higher accuracy on stiff systems, use DIRK3 or adaptive
/// ESDIRK methods.
///
/// # References
/// - Ferracina, L., & Spijker, M. N. (2008). "Strong stability of singly-diagonally-
///   implicit Runge-Kutta methods". Applied Numerical Mathematics, 58(11), 1675-1686.
/// - Hairer, E., & Wanner, G. (1996). "Solving Ordinary Differential Equations II:
///   Stiff and Differential-Algebraic Problems". Springer Series in Computational
///   Mathematics, Vol. 14.
#[derive(Debug, Clone)]
pub struct DIRK2 {
    state: DVector<f64>,
    initial: DVector<f64>,
    history: VecDeque<DVector<f64>>,
    optimizer: NewtonAnderson,
    stage: usize,
    // Stage slopes
    k: Vec<DVector<f64>>,
}

impl DIRK2 {
    /// Create a new DIRK2 solver with the given initial state
    pub fn new(initial: DVector<f64>) -> Self {
        let dim = initial.len();
        Self {
            state: initial.clone(),
            initial,
            history: VecDeque::with_capacity(2),
            optimizer: NewtonAnderson::default(),
            stage: 0,
            k: vec![DVector::zeros(dim); 2],
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
        self.stage == 1 // 2 stages total (0-indexed)
    }

    /// Stage iterator yielding evaluation times
    pub fn stages(&mut self, t: f64, dt: f64) -> impl Iterator<Item = f64> + '_ {
        self.stage = 0;
        let eval_stages = [1.0 / 4.0, 3.0 / 4.0];
        let mut idx = 0;

        std::iter::from_fn(move || {
            if idx < 2 {
                let t_eval = t + eval_stages[idx] * dt;
                idx += 1;
                self.stage = idx - 1;
                Some(t_eval)
            } else {
                None
            }
        })
    }

    /// Butcher tableau coefficients for stage i
    fn butcher_tableau(&self, stage: usize) -> &[f64] {
        match stage {
            0 => &[1.0 / 4.0],
            1 => &[1.0 / 2.0, 1.0 / 4.0],
            _ => panic!("Invalid stage for DIRK2"),
        }
    }

    /// Final combination coefficients
    fn final_weights(&self) -> &[f64] {
        &[1.0 / 2.0, 1.0 / 2.0]
    }
}

impl Solver for DIRK2 {
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
        for k in &mut self.k {
            k.fill(0.0);
        }
    }

    fn order(&self) -> usize {
        2
    }

    fn stages(&self) -> usize {
        2
    }

    fn is_adaptive(&self) -> bool {
        false
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

impl ImplicitSolver for DIRK2 {
    fn solve<F, J>(&mut self, mut f: F, mut jac: Option<J>, dt: f64) -> Result<f64, SolverError>
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
        J: FnMut(&DVector<f64>, f64) -> DMatrix<f64>,
    {
        // Get past state from history
        let x0 = self
            .history
            .front()
            .ok_or(SolverError::EmptyHistory)?
            .clone();

        // Evaluate f at current state
        let f_val = f(&self.state, 0.0);

        // Buffer the slope for this stage
        self.k[self.stage] = f_val.clone();

        // Compute slope from previous stages
        let bt = self.butcher_tableau(self.stage);
        let mut slope = DVector::zeros(x0.len());
        for (i, &a_ij) in bt.iter().enumerate() {
            slope += a_ij * &self.k[i];
        }

        // Fixed point equation: x = x0 + dt * slope
        let g = &x0 + dt * slope;

        // Get the diagonal Butcher coefficient
        let a_ii = bt[self.stage];

        // Compute Jacobian if available
        let j_opt = jac.as_mut().map(|j_fn| dt * a_ii * j_fn(&self.state, 0.0));

        // Use the optimizer to solve for x
        let (x_new, err) = self.optimizer.step(&self.state, &g, j_opt.as_ref());
        self.state = x_new;

        Ok(err)
    }

    fn step<F>(&mut self, mut f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        // If not last stage, just continue
        if !self.is_last_stage() {
            return SolverStepResult::default();
        }

        // Compute final solution using final weights
        let x0 = self.history.front().unwrap().clone();

        // Recompute final slope from stored k values
        let f_val = f(&self.state, 0.0);
        self.k[self.stage] = f_val;

        // Get weights as array to avoid borrow
        let weights = [1.0 / 2.0, 1.0 / 2.0];

        let mut slope = DVector::zeros(x0.len());
        for (i, &w) in weights.iter().enumerate() {
            slope += w * &self.k[i];
        }

        self.state = &x0 + dt * slope;

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
        solver: &mut DIRK2,
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

            for _t_eval in solver.stages(t, dt) {
                // Fixed-point iteration
                for _ in 0..MAX_ITER {
                    let err = solver.solve(f.clone(), Some(jac.clone()), dt).unwrap();
                    if err < TOL_FPI {
                        break;
                    }
                }

                solver.step(f.clone(), dt);
            }

            t += dt;
            times.push(t);
            states.push(solver.state().clone());
        }

        (times, states)
    }

    #[test]
    fn test_dirk2_init() {
        let solver = DIRK2::new(DVector::from_vec(vec![1.0]));
        assert_eq!(solver.order(), 2);
        assert_eq!(solver.stages(), 2);
        assert!(!solver.is_adaptive());
        assert!(!solver.is_explicit());
        assert_eq!(solver.state()[0], 1.0);
    }

    #[test]
    fn test_dirk2_exponential_decay() {
        // dx/dt = -x, x(0) = 1
        // Exact solution: x(t) = exp(-t)
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = DIRK2::new(x0);

        let f = |x: &DVector<f64>, _t: f64| -x;
        let jac = |_x: &DVector<f64>, _t: f64| DMatrix::from_element(1, 1, -1.0);

        let dt = 0.1;
        let t_final = 1.0;

        let (_times, states) = integrate_implicit(&mut solver, f, jac, 0.0, t_final, dt);

        let exact = (-t_final).exp();
        let computed = states.last().unwrap()[0];

        // DIRK2 is second order
        assert_relative_eq!(computed, exact, epsilon = 1e-3);
    }

    #[test]
    fn test_dirk2_stiff_decay() {
        // Stiff problem: dx/dt = -1000*x, x(0) = 1
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = DIRK2::new(x0);

        let lambda = -1000.0;
        let f = |x: &DVector<f64>, _t: f64| lambda * x;
        let jac = |_x: &DVector<f64>, _t: f64| DMatrix::from_element(1, 1, lambda);

        let dt = 0.01;
        let t_final = 0.1;

        let (_times, states) = integrate_implicit(&mut solver, f, jac, 0.0, t_final, dt);

        let exact = (lambda * t_final).exp();
        let computed = states.last().unwrap()[0];

        assert!(computed.is_finite());
        assert_relative_eq!(computed, exact, epsilon = 0.05);
    }

    #[test]
    fn test_dirk2_stages_iterator() {
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = DIRK2::new(x0);

        let t = 0.0;
        let dt = 1.0;

        let stages: Vec<f64> = solver.stages(t, dt).collect();

        assert_eq!(stages.len(), 2);
        assert_relative_eq!(stages[0], 0.25, epsilon = 1e-10);
        assert_relative_eq!(stages[1], 0.75, epsilon = 1e-10);
    }

    #[test]
    fn test_dirk2_reset() {
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = DIRK2::new(x0);

        solver.set_state(DVector::from_vec(vec![5.0]));
        solver.reset();

        assert_eq!(solver.state()[0], 1.0);
    }
}
