//! ESDIRK43 - Six-stage, 4th order ESDIRK with embedded 3rd order error estimate
//!
//! L-stable and stiffly accurate implicit Runge-Kutta method.

use nalgebra::DVector;
use std::collections::VecDeque;

use super::{ImplicitSolver, Solver, SolverError, SolverStepResult};

/// ESDIRK43 - Six-stage, 4(3) embedded ESDIRK method
///
/// Six-stage, 4th order ESDIRK method with embedded 3rd order error
/// estimate. L-stable and stiffly accurate.
///
/// # Characteristics
/// - Order: 4 (propagating) / 3 (embedded)
/// - Stages: 6 (1 explicit, 5 implicit)
/// - Adaptive timestep
/// - L-stable, stiffly accurate
/// - Stage order 2
///
/// # Note
/// Recommended default for stiff systems. L-stability damps
/// high-frequency parasitic modes that arise from stiff subsystems (e.g.
/// PID controllers with large derivative gain, fast electrical or chemical
/// dynamics). The adaptive step-size control concentrates computational
/// effort where the solution changes rapidly. For non-stiff systems,
/// RKDP54 avoids the implicit solve cost and is more efficient. For
/// tighter tolerances on stiff problems, ESDIRK54 provides 5th order
/// accuracy.
///
/// # References
/// - Kennedy, C. A., & Carpenter, M. H. (2019). "Diagonally implicit
///   Runge-Kutta methods for stiff ODEs". Applied Numerical
///   Mathematics, 146, 221-244.
/// - Hairer, E., & Wanner, G. (1996). "Solving Ordinary Differential
///   Equations II: Stiff and Differential-Algebraic Problems". Springer
///   Series in Computational Mathematics, Vol. 14.
#[derive(Debug, Clone)]
pub struct ESDIRK43 {
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

impl ESDIRK43 {
    /// Create a new ESDIRK43 solver with the given initial state
    pub fn new(initial: DVector<f64>) -> Self {
        Self::with_tolerances(initial, 1e-8, 1e-4)
    }

    /// Create a new ESDIRK43 solver with custom tolerances
    pub fn with_tolerances(initial: DVector<f64>, tol_abs: f64, tol_rel: f64) -> Self {
        let n = initial.len();
        Self {
            state: initial.clone(),
            initial,
            history: VecDeque::with_capacity(2),
            slopes: vec![DVector::zeros(n); 6],
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
    fn eval_stages(&self) -> [f64; 6] {
        let sqrt2 = 2.0_f64.sqrt();
        [
            0.0,
            1.0 / 2.0,
            (2.0 - sqrt2) / 4.0,
            2012122486997.0 / 3467029789466.0,
            1.0,
            1.0,
        ]
    }

    /// Compute error norm and timestep scale factor
    fn error_controller(&self, dt: f64) -> (bool, f64, f64) {
        // Coefficients for truncation error estimate (A1 - A2)
        let a1 = [
            657241292721.0 / 9909463049845.0,
            657241292721.0 / 9909463049845.0,
            1290772910128.0 / 5804808736437.0,
            1103522341516.0 / 2197678446715.0,
            -3.0 / 28.0,
            1.0 / 4.0,
        ];
        let a2 = [
            -71925161075.0 / 3900939759889.0,
            -71925161075.0 / 3900939759889.0,
            2973346383745.0 / 8160025745289.0,
            3972464885073.0 / 7694851252693.0,
            -263368882881.0 / 4213126269514.0,
            3295468053953.0 / 15064441987965.0,
        ];

        let mut tr = [0.0; 6];
        for i in 0..6 {
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
        let order = 3.min(4); // min(m, n)
        let mut timestep_scale = self.beta / error_norm.powf(1.0 / (order as f64 + 1.0));

        // Clip rescale factor to reasonable range
        timestep_scale = timestep_scale.clamp(0.1, 10.0);

        (success, error_norm, timestep_scale)
    }
}

impl Solver for ESDIRK43 {
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
        true
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

impl ImplicitSolver for ESDIRK43 {
    fn step<F>(&mut self, mut f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        // For explicit first stage or intermediate stages, no error estimate
        if self.stage < 5 {
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
        let bt: [Vec<f64>; 6] = [
            vec![], // Stage 0: explicit
            vec![1.0 / 4.0, 1.0 / 4.0],
            vec![
                -1356991263433.0 / 26208533697614.0,
                -1356991263433.0 / 26208533697614.0,
                1.0 / 4.0,
            ],
            vec![
                -1778551891173.0 / 14697912885533.0,
                -1778551891173.0 / 14697912885533.0,
                7325038566068.0 / 12797657924939.0,
                1.0 / 4.0,
            ],
            vec![
                -24076725932807.0 / 39344244018142.0,
                -24076725932807.0 / 39344244018142.0,
                9344023789330.0 / 6876721947151.0,
                11302510524611.0 / 18374767399840.0,
                1.0 / 4.0,
            ],
            vec![
                657241292721.0 / 9909463049845.0,
                657241292721.0 / 9909463049845.0,
                1290772910128.0 / 5804808736437.0,
                1103522341516.0 / 2197678446715.0,
                -3.0 / 28.0,
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
    fn test_esdirk43_properties() {
        let x0 = DVector::from_vec(vec![1.0]);
        let solver = ESDIRK43::new(x0);

        assert_eq!(solver.order(), 4);
        assert_eq!(solver.stages(), 6);
        assert!(solver.is_adaptive());
        assert!(!solver.is_explicit());
    }
}
