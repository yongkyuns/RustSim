//! ESDIRK85 - Sixteen-stage, 8th order ESDIRK with embedded 5th order error estimate
//!
//! L-stable and stiffly accurate implicit Runge-Kutta method.

use nalgebra::DVector;
use std::collections::VecDeque;

use super::{ImplicitSolver, Solver, SolverError, SolverStepResult};

/// ESDIRK85 - Sixteen-stage, 8(5) embedded ESDIRK method
///
/// Sixteen-stage, 8th order ESDIRK method with embedded 5th order error
/// estimate. L-stable and stiffly accurate (ESDIRK(16,8)[2]SAL-[(16,5)]).
///
/// # Characteristics
/// - Order: 8 (propagating) / 5 (embedded)
/// - Stages: 16 (1 explicit, 15 implicit)
/// - Adaptive timestep
/// - L-stable, stiffly accurate
/// - Stage order 2
///
/// # Note
/// Fifteen implicit solves per step make this very expensive. It is only
/// justified when the right-hand side evaluation is itself costly (large
/// state dimension, expensive ODE blocks) and very tight tolerances are
/// required so that the 8th order convergence compensates through much
/// larger steps. For generating stiff reference solutions to validate other
/// solvers. In almost all practical simulations, ESDIRK54
/// is the better choice.
///
/// # References
/// - Alamri, Y., & Ketcheson, D. I. (2024). "Very high-order A-stable
///   stiffly accurate diagonally implicit Runge-Kutta methods with
///   error estimators". Journal of Scientific Computing, 100,
///   Article 84.
/// - Kennedy, C. A., & Carpenter, M. H. (2019). "Diagonally implicit
///   Runge-Kutta methods for stiff ODEs". Applied Numerical
///   Mathematics, 146, 221-244.
#[derive(Debug, Clone)]
pub struct ESDIRK85 {
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

impl ESDIRK85 {
    /// Create a new ESDIRK85 solver with the given initial state
    pub fn new(initial: DVector<f64>) -> Self {
        Self::with_tolerances(initial, 1e-8, 1e-4)
    }

    /// Create a new ESDIRK85 solver with custom tolerances
    pub fn with_tolerances(initial: DVector<f64>, tol_abs: f64, tol_rel: f64) -> Self {
        let n = initial.len();
        Self {
            state: initial.clone(),
            initial,
            history: VecDeque::with_capacity(2),
            slopes: vec![DVector::zeros(n); 16],
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
    fn eval_stages(&self) -> [f64; 16] {
        [
            0.0,
            0.234637638717043,
            0.558545926594724,
            0.562667638694992,
            0.697898381329126,
            0.956146958839776,
            0.812903043340468,
            0.148256733818785,
            0.944650387704291,
            0.428471803715736,
            0.984131639774509,
            0.320412672954752,
            0.974077670791771,
            0.852850433853921,
            0.823320301074444,
            1.0,
        ]
    }

    /// Compute error norm and timestep scale factor
    fn error_controller(&self, dt: f64) -> (bool, f64, f64) {
        // Coefficients for truncation error estimate (A1 - A2)
        let a1 = [
            0.0459979286336779,
            0.0780075394482806,
            0.015021874148058,
            0.195180277284195,
            -0.00246643310153235,
            0.0473977117068314,
            -0.0682773558610363,
            0.19568019123878,
            -0.0876765449323747,
            0.177874852409192,
            -0.337519251582222,
            -0.0123255553640736,
            0.311573291192553,
            0.0458604327754991,
            0.278352222645651,
            0.117318819358521,
        ];
        let a2 = [
            0.0603373529853206,
            0.175453809423998,
            0.0537707777611352,
            0.195309248607308,
            0.0135893741970232,
            -0.0221160259296707,
            -0.00726526156430691,
            0.102961059369124,
            0.000900215457460583,
            0.0547959465692338,
            -0.334995726863153,
            0.0464409662093384,
            0.301388101652194,
            0.00524851570622031,
            0.229538601845236,
            0.124643044573514,
        ];

        let mut tr = [0.0; 16];
        for i in 0..16 {
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
        let order = 5; // min(m, n)
        let mut timestep_scale = self.beta / error_norm.powf(1.0 / (order as f64 + 1.0));

        // Clip rescale factor to reasonable range
        timestep_scale = timestep_scale.clamp(0.1, 10.0);

        (success, error_norm, timestep_scale)
    }

    /// Get Butcher tableau for a given stage
    fn butcher_tableau(&self, stage: usize) -> Vec<f64> {
        let gamma = 0.117318819358521;

        match stage {
            0 => vec![],
            1 => vec![gamma, gamma],
            2 => vec![0.0557014605974616, 0.385525646638742, gamma],
            3 => vec![
                0.063493276428895,
                0.373556126263681,
                0.0082994166438953,
                gamma,
            ],
            4 => vec![
                0.0961351856230088,
                0.335558324517178,
                0.207077765910132,
                -0.0581917140797146,
                gamma,
            ],
            5 => vec![
                0.0497669214238319,
                0.384288616546039,
                0.0821728117583936,
                0.120337007107103,
                0.202262782645888,
                gamma,
            ],
            6 => vec![
                0.00626710666809847,
                0.496491452640725,
                -0.111303249827358,
                0.170478821683603,
                0.166517073971103,
                -0.0328669811542241,
                gamma,
            ],
            7 => vec![
                0.0463439767281591,
                0.00306724391019652,
                -0.00816305222386205,
                -0.0353302599538294,
                0.0139313601702569,
                -0.00992014507967429,
                0.0210087909090165,
                gamma,
            ],
            8 => vec![
                0.111574049232048,
                0.467639166482209,
                0.237773114804619,
                0.0798895699267508,
                0.109580615914593,
                0.0307353103825936,
                -0.0404391509541147,
                -0.16942110744293,
                gamma,
            ],
            9 => vec![
                -0.0107072484863877,
                -0.231376703354252,
                0.017541113036611,
                0.144871527682418,
                -0.041855459769806,
                0.0841832168332261,
                -0.0850020937282192,
                0.486170343825899,
                -0.0526717116822739,
                gamma,
            ],
            10 => vec![
                -0.0142238262314935,
                0.14752923682514,
                0.238235830732566,
                0.037950291904103,
                0.252075123381518,
                0.0474266904224567,
                -0.00363139069342027,
                0.274081442388563,
                -0.0599166970745255,
                -0.0527138812389185,
                gamma,
            ],
            11 => vec![
                -0.11837020183211,
                -0.635712481821264,
                0.239738832602538,
                0.330058936651707,
                -0.325784087988237,
                -0.0506514314589253,
                -0.281914404487009,
                0.852596345144291,
                0.651444614298805,
                -0.103476387303591,
                -0.354835880209975,
                gamma,
            ],
            12 => vec![
                -0.00458164025442349,
                0.296219694015248,
                0.322146049419995,
                0.15917778285238,
                0.284864871688843,
                0.185509526463076,
                -0.0784621067883274,
                0.166312223692047,
                -0.284152486083397,
                -0.357125104338944,
                0.078437074055306,
                0.0884129667114481,
                gamma,
            ],
            13 => vec![
                -0.0545561913848106,
                0.675785423442753,
                0.423066443201941,
                -0.000165300126841193,
                0.104252994793763,
                -0.105763019303021,
                -0.15988308809318,
                0.0515050001032011,
                0.56013979290924,
                -0.45781539708603,
                -0.255870699752664,
                0.026960254296416,
                -0.0721245985053681,
                gamma,
            ],
            14 => vec![
                0.0649253995775223,
                -0.0216056457922249,
                -0.073738139377975,
                0.0931033310077225,
                -0.0194339577299149,
                -0.0879623837313009,
                0.057125517179467,
                0.205120850488097,
                0.132576503537441,
                0.489416890627328,
                -0.1106765720501,
                -0.081038793996096,
                0.0606031613503788,
                -0.00241467937442272,
                gamma,
            ],
            15 => vec![
                0.0459979286336779,
                0.0780075394482806,
                0.015021874148058,
                0.195180277284195,
                -0.00246643310153235,
                0.0473977117068314,
                -0.0682773558610363,
                0.19568019123878,
                -0.0876765449323747,
                0.177874852409192,
                -0.337519251582222,
                -0.0123255553640736,
                0.311573291192553,
                0.0458604327754991,
                0.278352222645651,
                gamma,
            ],
            _ => panic!("Invalid stage"),
        }
    }
}

impl Solver for ESDIRK85 {
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
        8
    }

    fn stages(&self) -> usize {
        16
    }

    fn is_adaptive(&self) -> bool {
        true
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

impl ImplicitSolver for ESDIRK85 {
    fn step<F>(&mut self, _f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
    {
        // For explicit first stage or intermediate stages, no error estimate
        if self.stage < 15 {
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

    fn solve<F, J>(&mut self, mut f: F, _jac: Option<J>, dt: f64) -> Result<f64, SolverError>
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

        // Get Butcher tableau for current stage
        let bt = self.butcher_tableau(self.stage);

        // Diagonal entry (gamma)
        let gamma = 0.117318819358521;

        // Compute explicit part of the slope
        let mut slope_explicit = DVector::zeros(x0.len());
        for (i, &coef) in bt.iter().enumerate().take(self.stage) {
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
    fn test_esdirk85_properties() {
        let x0 = DVector::from_vec(vec![1.0]);
        let solver = ESDIRK85::new(x0);

        assert_eq!(solver.order(), 8);
        assert_eq!(solver.stages(), 16);
        assert!(solver.is_adaptive());
        assert!(!solver.is_explicit());
    }
}
