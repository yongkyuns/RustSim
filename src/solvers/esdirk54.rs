//! ESDIRK54 - Seven-stage, 5th order ESDIRK with embedded 4th order error estimate
//!
//! L-stable and stiffly accurate implicit Runge-Kutta method.

use nalgebra::DVector;
use std::collections::VecDeque;

use super::{ImplicitSolver, Solver, SolverError, SolverStepResult};

/// ESDIRK54 - Seven-stage, 5(4) embedded ESDIRK method
///
/// Seven-stage, 5th order ESDIRK method with embedded 4th order error
/// estimate. L-stable and stiffly accurate (ESDIRK5(4)7L[2]SA2).
///
/// # Characteristics
/// - Order: 5 (propagating) / 4 (embedded)
/// - Stages: 7 (1 explicit, 6 implicit)
/// - Adaptive timestep
/// - L-stable, stiffly accurate
/// - Stage order 2
///
/// # Note
/// The highest-accuracy L-stable single-step solver in this library before
/// the much more expensive ESDIRK85. Use when tight tolerances are
/// needed on a stiff system (e.g. multi-rate systems combining fast
/// electrical and slow thermal dynamics). At moderate tolerances,
/// ESDIRK43 achieves similar results with fewer implicit solves per
/// step.
///
/// # References
/// - Kennedy, C. A., & Carpenter, M. H. (2019). "Diagonally implicit
///   Runge-Kutta methods for stiff ODEs". Applied Numerical
///   Mathematics, 146, 221-244.
/// - Hairer, E., & Wanner, G. (1996). "Solving Ordinary Differential
///   Equations II: Stiff and Differential-Algebraic Problems". Springer
///   Series in Computational Mathematics, Vol. 14.
#[derive(Debug, Clone)]
pub struct ESDIRK54 {
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

impl ESDIRK54 {
    /// Create a new ESDIRK54 solver with the given initial state
    pub fn new(initial: DVector<f64>) -> Self {
        Self::with_tolerances(initial, 1e-8, 1e-4)
    }

    /// Create a new ESDIRK54 solver with custom tolerances
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
    fn eval_stages(&self) -> [f64; 7] {
        [
            0.0,
            46.0 / 125.0,
            7121331996143.0 / 11335814405378.0,
            49.0 / 353.0,
            3706679970760.0 / 5295570149437.0,
            347.0 / 382.0,
            1.0,
        ]
    }

    /// Compute error norm and timestep scale factor
    fn error_controller(&self, dt: f64) -> (bool, f64, f64) {
        // Coefficients for truncation error estimate (A1 - A2)
        let a1 = [
            -188593204321.0 / 4778616380481.0,
            -188593204321.0 / 4778616380481.0,
            2809310203510.0 / 10304234040467.0,
            1021729336898.0 / 2364210264653.0,
            870612361811.0 / 2470410392208.0,
            -1307970675534.0 / 8059683598661.0,
            23.0 / 125.0,
        ];
        let a2 = [
            -582099335757.0 / 7214068459310.0,
            -582099335757.0 / 7214068459310.0,
            615023338567.0 / 3362626566945.0,
            3192122436311.0 / 6174152374399.0,
            6156034052041.0 / 14430468657929.0,
            -1011318518279.0 / 9693750372484.0,
            1914490192573.0 / 13754262428401.0,
        ];

        let mut tr = [0.0; 7];
        for i in 0..7 {
            tr[i] = a1[i] - a2[i];
        }

        // Compute truncation error slope
        let mut error_slope = DVector::zeros(self.state.len());
        for (i, &coef) in tr.iter().enumerate() {
            error_slope += coef * &self.slopes[i];
        }

        // Compute scaling factors
        let scale = self.state.map(|x| self.tol_abs + self.tol_rel * x.abs());

        // Compute scaled error
        let scaled_error = (dt * &error_slope).component_div(&scale).map(|e| e.abs());

        // Error norm (max norm) with lower bound
        let error_norm = scaled_error.max().max(1e-16);

        // Determine if error is acceptable
        let success = error_norm <= 1.0;

        // Compute timestep scale factor
        let order = 4.min(5); // min(m, n)
        let mut timestep_scale = self.beta / error_norm.powf(1.0 / (order as f64 + 1.0));

        // Clip rescale factor to reasonable range
        timestep_scale = timestep_scale.clamp(0.1, 10.0);

        (success, error_norm, timestep_scale)
    }
}

impl Solver for ESDIRK54 {
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
        false
    }
}

impl ImplicitSolver for ESDIRK54 {
    fn step<F>(&mut self, mut f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        // For explicit first stage or intermediate stages, no error estimate
        if self.stage < 6 {
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
        let bt: [Vec<f64>; 7] = [
            vec![], // Stage 0: explicit
            vec![23.0 / 125.0, 23.0 / 125.0],
            vec![
                791020047304.0 / 3561426431547.0,
                791020047304.0 / 3561426431547.0,
                23.0 / 125.0,
            ],
            vec![
                -158159076358.0 / 11257294102345.0,
                -158159076358.0 / 11257294102345.0,
                -85517644447.0 / 5003708988389.0,
                23.0 / 125.0,
            ],
            vec![
                -1653327111580.0 / 4048416487981.0,
                -1653327111580.0 / 4048416487981.0,
                1514767744496.0 / 9099671765375.0,
                14283835447591.0 / 12247432691556.0,
                23.0 / 125.0,
            ],
            vec![
                -4540011970825.0 / 8418487046959.0,
                -4540011970825.0 / 8418487046959.0,
                -1790937573418.0 / 7393406387169.0,
                10819093665085.0 / 7266595846747.0,
                4109463131231.0 / 7386972500302.0,
                23.0 / 125.0,
            ],
            vec![
                -188593204321.0 / 4778616380481.0,
                -188593204321.0 / 4778616380481.0,
                2809310203510.0 / 10304234040467.0,
                1021729336898.0 / 2364210264653.0,
                870612361811.0 / 2470410392208.0,
                -1307970675534.0 / 8059683598661.0,
                23.0 / 125.0,
            ],
        ];

        // Diagonal entry (gamma)
        let gamma = 23.0 / 125.0;

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
    fn test_esdirk54_properties() {
        let x0 = DVector::from_vec(vec![1.0]);
        let solver = ESDIRK54::new(x0);

        assert_eq!(solver.order(), 5);
        assert_eq!(solver.stages(), 7);
        assert!(solver.is_adaptive());
        assert!(!solver.is_explicit());
    }
}
