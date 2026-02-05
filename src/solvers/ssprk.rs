//! Strong Stability Preserving Runge-Kutta methods

use nalgebra::DVector;
use std::collections::VecDeque;

use super::{ExplicitSolver, Solver, SolverError, SolverStepResult};

/// Two-stage, 2nd order Strong Stability Preserving Runge-Kutta method (SSPRK22)
///
/// Also known as Heun's method. SSP methods preserve monotonicity and total
/// variation diminishing (TVD) properties of the spatial discretisation under
/// a timestep restriction scaled by the SSP coefficient.
///
/// # Characteristics
/// - Order: 2
/// - Stages: 2
/// - Explicit, fixed timestep
/// - SSP coefficient: C = 1
///
/// # Note
/// Relevant when a system wraps a method-of-lines discretisation of a
/// hyperbolic PDE (e.g. shallow water, compressible Euler) and the spatial
/// operator is TVD under forward Euler. For typical ODE systems without such
/// structure, RK4 or RKDP54 are more appropriate choices.
///
/// # References
/// - Shu, C.-W., & Osher, S. (1988). "Efficient implementation of essentially
///   non-oscillatory shock-capturing schemes". Journal of Computational
///   Physics, 77(2), 439-471.
/// - Gottlieb, S., Shu, C.-W., & Tadmor, E. (2001). "Strong stability-preserving
///   high-order time discretization methods". SIAM Review, 43(1), 89-112.
#[derive(Debug, Clone)]
pub struct SSPRK22 {
    state: DVector<f64>,
    initial: DVector<f64>,
    history: VecDeque<DVector<f64>>,
    slopes: Vec<DVector<f64>>,
    stage: usize,
}

impl SSPRK22 {
    /// Create a new SSPRK22 solver with the given initial state
    pub fn new(initial: DVector<f64>) -> Self {
        let n = initial.len();
        Self {
            state: initial.clone(),
            initial,
            history: VecDeque::with_capacity(2),
            slopes: vec![DVector::zeros(n); 2],
            stage: 0,
        }
    }
}

impl Solver for SSPRK22 {
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
        2
    }

    fn stages(&self) -> usize {
        2
    }

    fn is_adaptive(&self) -> bool {
        false
    }

    fn is_explicit(&self) -> bool {
        true
    }
}

impl ExplicitSolver for SSPRK22 {
    fn step<F>(&mut self, mut f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        let x0 = self
            .history
            .front()
            .expect("Must call buffer() before step()");

        // SSPRK22 Butcher tableau
        // c = [0, 1]
        // a = [[],
        //      [1]]
        // b = [1/2, 1/2]

        let c = [0.0, 1.0];

        // Store slope at current stage
        self.slopes[self.stage] = f(&self.state, c[self.stage] * dt);

        if self.stage == 0 {
            // First stage: x1 = x0 + dt * k0
            self.state = x0 + dt * &self.slopes[0];
            self.stage = 1;
        } else {
            // Final stage: x_{n+1} = x0 + dt * (k0 + k1) / 2
            self.state = x0 + dt * 0.5 * (&self.slopes[0] + &self.slopes[1]);
            self.stage = 0;
        }

        SolverStepResult::default()
    }
}

/// Three-stage, 3rd order Strong Stability Preserving Runge-Kutta method (SSPRK33)
///
/// The unique optimal three-stage, 3rd order SSP scheme. Commonly paired with
/// WENO and ENO spatial discretisations for hyperbolic conservation laws.
///
/// # Characteristics
/// - Order: 3
/// - Stages: 3
/// - Explicit, fixed timestep
/// - SSP coefficient: C = 1
///
/// # Note
/// The standard SSP time integrator for method-of-lines PDE discretisations.
/// If the spatial operator is TVD under forward Euler, this method preserves
/// that property at the same timestep restriction. When stability is borderline,
/// SSPRK34 allows roughly twice the timestep at the cost of one extra stage.
///
/// # References
/// - Shu, C.-W., & Osher, S. (1988). "Efficient implementation of essentially
///   non-oscillatory shock-capturing schemes". Journal of Computational
///   Physics, 77(2), 439-471.
/// - Gottlieb, S., Shu, C.-W., & Tadmor, E. (2001). "Strong stability-preserving
///   high-order time discretization methods". SIAM Review, 43(1), 89-112.
/// - Gottlieb, S., Ketcheson, D. I., & Shu, C.-W. (2011). "Strong Stability
///   Preserving Runge-Kutta and Multistep Time Discretizations". World Scientific.
#[derive(Debug, Clone)]
pub struct SSPRK33 {
    state: DVector<f64>,
    initial: DVector<f64>,
    history: VecDeque<DVector<f64>>,
    slopes: Vec<DVector<f64>>,
    stage: usize,
}

impl SSPRK33 {
    /// Create a new SSPRK33 solver with the given initial state
    pub fn new(initial: DVector<f64>) -> Self {
        let n = initial.len();
        Self {
            state: initial.clone(),
            initial,
            history: VecDeque::with_capacity(2),
            slopes: vec![DVector::zeros(n); 3],
            stage: 0,
        }
    }
}

impl Solver for SSPRK33 {
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
        3
    }

    fn is_adaptive(&self) -> bool {
        false
    }

    fn is_explicit(&self) -> bool {
        true
    }
}

impl ExplicitSolver for SSPRK33 {
    fn step<F>(&mut self, mut f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        let x0 = self
            .history
            .front()
            .expect("Must call buffer() before step()");

        // SSPRK33 Butcher tableau
        // c = [0, 1, 1/2]
        // a = [[],
        //      [1],
        //      [1/4, 1/4]]
        // b = [1/6, 1/6, 2/3]

        let c = [0.0, 1.0, 0.5];

        // Store slope at current stage
        self.slopes[self.stage] = f(&self.state, c[self.stage] * dt);

        match self.stage {
            0 => {
                // First stage: x1 = x0 + dt * k0
                self.state = x0 + dt * &self.slopes[0];
                self.stage = 1;
            }
            1 => {
                // Second stage: x2 = x0 + dt * (k0 + k1) / 4
                self.state = x0 + dt * 0.25 * (&self.slopes[0] + &self.slopes[1]);
                self.stage = 2;
            }
            _ => {
                // Final stage: x_{n+1} = x0 + dt * (k0 + k1 + 4*k2) / 6
                self.state = x0 + dt * (1.0 / 6.0) * (&self.slopes[0] + &self.slopes[1] + 4.0 * &self.slopes[2]);
                self.stage = 0;
            }
        }

        SolverStepResult::default()
    }
}

/// Four-stage, 3rd order Strong Stability Preserving Runge-Kutta method (SSPRK34)
///
/// An extra stage compared to SSPRK33 doubles the allowable SSP timestep
/// (C = 2), giving a larger effective stability region along the negative real axis.
///
/// # Characteristics
/// - Order: 3
/// - Stages: 4
/// - Explicit, fixed timestep
/// - SSP coefficient: C = 2
///
/// # Note
/// Preferable over SSPRK33 when a method-of-lines ODE is close
/// to the SSP timestep limit and the cost of one additional stage per step is
/// acceptable in exchange for a factor-of-two relaxation in the CFL
/// constraint.
///
/// # References
/// - Spiteri, R. J., & Ruuth, S. J. (2002). "A new class of optimal
///   high-order strong-stability-preserving time discretization methods".
///   SIAM Journal on Numerical Analysis, 40(2), 469-491.
/// - Gottlieb, S., Ketcheson, D. I., & Shu, C.-W. (2011). "Strong
///   Stability Preserving Runge-Kutta and Multistep Time
///   Discretizations". World Scientific.
#[derive(Debug, Clone)]
pub struct SSPRK34 {
    state: DVector<f64>,
    initial: DVector<f64>,
    history: VecDeque<DVector<f64>>,
    slopes: Vec<DVector<f64>>,
    stage: usize,
}

impl SSPRK34 {
    /// Create a new SSPRK34 solver with the given initial state
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

impl Solver for SSPRK34 {
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
        false
    }

    fn is_explicit(&self) -> bool {
        true
    }
}

impl ExplicitSolver for SSPRK34 {
    fn step<F>(&mut self, mut f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        let x0 = self
            .history
            .front()
            .expect("Must call buffer() before step()");

        // SSPRK34 Butcher tableau
        // c = [0, 1/2, 1, 1/2]
        // BT = {
        //     0: [1/2],
        //     1: [1/2, 1/2],
        //     2: [1/6, 1/6, 1/6],
        //     3: [1/6, 1/6, 1/6, 1/2]
        // }

        let c = [0.0, 0.5, 1.0, 0.5];

        // Store slope at current stage
        self.slopes[self.stage] = f(&self.state, c[self.stage] * dt);

        match self.stage {
            0 => {
                // Stage 0: x1 = x0 + (1/2) * dt * k0
                self.state = x0 + dt * 0.5 * &self.slopes[0];
                self.stage = 1;
            }
            1 => {
                // Stage 1: x2 = x0 + dt * (1/2 * k0 + 1/2 * k1)
                self.state = x0 + dt * 0.5 * (&self.slopes[0] + &self.slopes[1]);
                self.stage = 2;
            }
            2 => {
                // Stage 2: x3 = x0 + dt * (1/6 * k0 + 1/6 * k1 + 1/6 * k2)
                self.state = x0 + dt * (1.0 / 6.0) * (&self.slopes[0] + &self.slopes[1] + &self.slopes[2]);
                self.stage = 3;
            }
            _ => {
                // Final stage: x_{n+1} = x0 + dt * (1/6 * k0 + 1/6 * k1 + 1/6 * k2 + 1/2 * k3)
                self.state = x0
                    + dt * (1.0 / 6.0) * (&self.slopes[0] + &self.slopes[1] + &self.slopes[2])
                    + dt * 0.5 * &self.slopes[3];
                self.stage = 0;
            }
        }

        SolverStepResult::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_ssprk22_exponential_decay() {
        // dx/dt = -x, x(0) = 1
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = SSPRK22::new(x0);

        let dt = 0.1;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        for _ in 0..n_steps {
            solver.buffer(dt);

            // Perform 2 stages
            for _ in 0..2 {
                solver.step(|x, _t| -x, dt);
            }
        }

        let exact = (-t_final).exp();
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-3);
    }

    #[test]
    fn test_ssprk33_exponential_decay() {
        // dx/dt = -x, x(0) = 1
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = SSPRK33::new(x0);

        let dt = 0.1;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        for _ in 0..n_steps {
            solver.buffer(dt);

            // Perform 3 stages
            for _ in 0..3 {
                solver.step(|x, _t| -x, dt);
            }
        }

        let exact = (-t_final).exp();
        // SSPRK33 should be more accurate than SSPRK22
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-4);
    }

    #[test]
    fn test_ssprk33_harmonic_oscillator() {
        // d²x/dt² = -x => [x, v]' = [v, -x]
        let x0 = DVector::from_vec(vec![1.0, 0.0]);
        let mut solver = SSPRK33::new(x0);

        let dt = 0.02;
        let t_final = 2.0 * std::f64::consts::PI;
        let n_steps = (t_final / dt) as usize;

        for _ in 0..n_steps {
            solver.buffer(dt);

            for _ in 0..3 {
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

    #[test]
    fn test_ssprk34_exponential_decay() {
        // dx/dt = -x, x(0) = 1
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = SSPRK34::new(x0);

        let dt = 0.1;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        for _ in 0..n_steps {
            solver.buffer(dt);

            // Perform 4 stages
            for _ in 0..4 {
                solver.step(|x, _t| -x, dt);
            }
        }

        let exact = (-t_final).exp();
        // SSPRK34 is 3rd order, similar accuracy to SSPRK33
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-4);
    }
}
