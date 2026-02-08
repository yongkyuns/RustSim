//! Discrete-time Kalman filter for state estimation
//!
//! Uses standard linear algebra notation (uppercase for matrices: F, H, Q, R, P, B, K, S)

#![allow(non_snake_case)]

use crate::block::{Block, DynamicBlock, StepResult};
use nalgebra::{DMatrix, DVector};

/// Discrete-time Kalman filter for state estimation.
///
/// Implements the standard Kalman filter algorithm to estimate the state of a
/// linear dynamic system from noisy measurements. The filter recursively updates
/// state estimates by combining predictions from a system model with incoming
/// measurements, weighted by their respective uncertainties.
///
/// The filter processes measurements at each time step through a two-stage process:
/// prediction (using the system model) and update (incorporating measurements).
///
/// # System Model
///
/// ```text
/// x_{k+1} = F x_k + B u_k + w_k
/// z_k = H x_k + v_k
/// ```
///
/// where `w_k ~ N(0, Q)` is process noise and `v_k ~ N(0, R)` is measurement noise.
///
/// # Algorithm
///
/// At each time step, the filter performs:
///
/// **Prediction:**
/// ```text
/// x_{k|k-1} = F x_{k-1} + B u_k
/// P_{k|k-1} = F P_{k-1} F^T + Q
/// ```
///
/// **Update:**
/// ```text
/// y_k = z_k - H x_{k|k-1}
/// S_k = H P_{k|k-1} H^T + R
/// K_k = P_{k|k-1} H^T S_k^{-1}
/// x_k = x_{k|k-1} + K_k y_k
/// P_k = (I - K_k H) P_{k|k-1}
/// ```
///
/// # Inputs
///
/// The block expects inputs in the following order:
/// - First `m` inputs: measurements `z`
/// - Next `p` inputs (if B is provided): control inputs `u`
///
/// # Outputs
///
/// The block outputs the n-dimensional state estimate `x`.
///
/// # Example
///
/// ```ignore
/// // Constant velocity model: state = [position, velocity]
/// let dt = 0.1;
/// let F = DMatrix::from_row_slice(2, 2, &[1.0, dt, 0.0, 1.0]);
/// let H = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);  // measure position
/// let Q = DMatrix::identity(2, 2) * 0.01;
/// let R = DMatrix::identity(1, 1) * 0.1;
///
/// let mut kf = KalmanFilter::new(F, H, Q, R);
///
/// // Update with measurement
/// kf.set_input(0, measured_position);
/// kf.update(0.0);
/// kf.step(0.0, dt);
///
/// // Get state estimate
/// let position = kf.get_output(0);
/// let velocity = kf.get_output(1);
/// ```
#[derive(Debug, Clone)]
pub struct KalmanFilter {
    // I/O
    inputs: Vec<f64>,
    outputs: Vec<f64>,

    // System matrices
    F: DMatrix<f64>,         // State transition matrix (n x n)
    H: DMatrix<f64>,         // Measurement matrix (m x n)
    Q: DMatrix<f64>,         // Process noise covariance (n x n)
    R: DMatrix<f64>,         // Measurement noise covariance (m x m)
    B: Option<DMatrix<f64>>, // Control input matrix (n x p)

    // State
    x: DVector<f64>, // State estimate (n)
    P: DMatrix<f64>, // Error covariance (n x n)

    // Initial conditions
    x0: DVector<f64>,
    P0: DMatrix<f64>,

    // Dimensions
    n: usize, // State dimension
    m: usize, // Measurement dimension
    p: usize, // Control input dimension

    // For buffer/revert
    buffered_x: DVector<f64>,
    buffered_P: DMatrix<f64>,

    // Discrete time step
    dt: Option<f64>,
}

impl KalmanFilter {
    /// Create a Kalman filter without control input
    ///
    /// # Arguments
    ///
    /// * `F` - State transition matrix (n x n)
    /// * `H` - Measurement matrix (m x n)
    /// * `Q` - Process noise covariance (n x n)
    /// * `R` - Measurement noise covariance (m x m)
    ///
    /// # Panics
    ///
    /// Panics if matrix dimensions are incompatible.
    pub fn new(F: DMatrix<f64>, H: DMatrix<f64>, Q: DMatrix<f64>, R: DMatrix<f64>) -> Self {
        Self::with_options(F, H, Q, R, None, None, None, None)
    }

    /// Create a Kalman filter with control input
    ///
    /// # Arguments
    ///
    /// * `F` - State transition matrix (n x n)
    /// * `H` - Measurement matrix (m x n)
    /// * `Q` - Process noise covariance (n x n)
    /// * `R` - Measurement noise covariance (m x m)
    /// * `B` - Control input matrix (n x p)
    pub fn with_control(
        F: DMatrix<f64>,
        H: DMatrix<f64>,
        Q: DMatrix<f64>,
        R: DMatrix<f64>,
        B: DMatrix<f64>,
    ) -> Self {
        Self::with_options(F, H, Q, R, Some(B), None, None, None)
    }

    /// Create a Kalman filter with all options
    ///
    /// # Arguments
    ///
    /// * `F` - State transition matrix (n x n)
    /// * `H` - Measurement matrix (m x n)
    /// * `Q` - Process noise covariance (n x n)
    /// * `R` - Measurement noise covariance (m x m)
    /// * `B` - Optional control input matrix (n x p)
    /// * `x0` - Optional initial state estimate (defaults to zeros)
    /// * `P0` - Optional initial error covariance (defaults to identity)
    /// * `dt` - Optional discrete time step
    #[allow(clippy::too_many_arguments)]
    pub fn with_options(
        F: DMatrix<f64>,
        H: DMatrix<f64>,
        Q: DMatrix<f64>,
        R: DMatrix<f64>,
        B: Option<DMatrix<f64>>,
        x0: Option<DVector<f64>>,
        P0: Option<DMatrix<f64>>,
        dt: Option<f64>,
    ) -> Self {
        // Validate dimensions
        let (n, n2) = F.shape();
        assert_eq!(n, n2, "F must be square");

        let (m, n3) = H.shape();
        assert_eq!(n, n3, "H must have same number of columns as F has rows");

        let (q1, q2) = Q.shape();
        assert_eq!(n, q1, "Q must have same size as F");
        assert_eq!(n, q2, "Q must be square");

        let (r1, r2) = R.shape();
        assert_eq!(m, r1, "R must have same size as H has rows");
        assert_eq!(m, r2, "R must be square");

        let p = if let Some(ref b) = B {
            let (b1, p_size) = b.shape();
            assert_eq!(n, b1, "B must have same number of rows as F");
            p_size
        } else {
            0
        };

        // Initial conditions
        let x0 = x0.unwrap_or_else(|| DVector::zeros(n));
        let P0 = P0.unwrap_or_else(|| DMatrix::identity(n, n));

        assert_eq!(x0.len(), n, "x0 must have size n");
        assert_eq!(P0.shape(), (n, n), "P0 must be n x n");

        let x = x0.clone();
        let P = P0.clone();

        let num_inputs = m + p;
        let num_outputs = n;

        Self {
            inputs: vec![0.0; num_inputs],
            outputs: vec![0.0; num_outputs],
            F,
            H,
            Q,
            R,
            B,
            x: x.clone(),
            P: P.clone(),
            x0,
            P0,
            n,
            m,
            p,
            buffered_x: x,
            buffered_P: P,
            dt,
        }
    }

    /// Get state dimension
    pub fn state_dim(&self) -> usize {
        self.n
    }

    /// Get measurement dimension
    pub fn measurement_dim(&self) -> usize {
        self.m
    }

    /// Get control input dimension
    pub fn control_dim(&self) -> usize {
        self.p
    }

    /// Get current state estimate
    pub fn state_estimate(&self) -> &DVector<f64> {
        &self.x
    }

    /// Get current error covariance
    pub fn error_covariance(&self) -> &DMatrix<f64> {
        &self.P
    }

    /// Perform Kalman filter update step manually
    ///
    /// This method can be called to trigger a filter update when `dt` is `None`.
    /// When `dt` is set, updates are triggered automatically via `step()`.
    pub fn update_filter(&mut self) {
        self.kf_update();
    }

    /// Perform Kalman filter update step
    fn kf_update(&mut self) {
        // Split inputs into measurements and control
        let z = DVector::from_column_slice(&self.inputs[0..self.m]);
        let u = if self.p > 0 {
            Some(DVector::from_column_slice(
                &self.inputs[self.m..self.m + self.p],
            ))
        } else {
            None
        };

        // Prediction step
        let x_pred = &self.F * &self.x
            + if let Some(ref B) = self.B {
                if let Some(ref u_vec) = u {
                    B * u_vec
                } else {
                    DVector::zeros(self.n)
                }
            } else {
                DVector::zeros(self.n)
            };

        let P_pred = &self.F * &self.P * self.F.transpose() + &self.Q;

        // Update step
        let y = z - &self.H * &x_pred; // Innovation
        let S = &self.H * &P_pred * self.H.transpose() + &self.R; // Innovation covariance

        // Kalman gain: K = P_pred * H^T * S^{-1}
        let K = &P_pred
            * self.H.transpose()
            * S.try_inverse()
                .expect("Innovation covariance S is singular");

        // Update state and covariance
        self.x = x_pred + &K * y;
        self.P = (DMatrix::identity(self.n, self.n) - &K * &self.H) * &P_pred;

        // Update outputs
        for i in 0..self.n {
            self.outputs[i] = self.x[i];
        }
    }
}

impl Block for KalmanFilter {
    const NUM_INPUTS: usize = 0; // Dynamic
    const NUM_OUTPUTS: usize = 0; // Dynamic
    const IS_DYNAMIC: bool = true;

    fn inputs(&self) -> &[f64] {
        &self.inputs
    }

    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut self.inputs
    }

    fn outputs(&self) -> &[f64] {
        &self.outputs
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        &mut self.outputs
    }

    fn update(&mut self, _t: f64) {
        // For discrete-time Kalman filter, update happens in step()
        // or when explicitly triggered if dt is None
        // Just output current state
        for i in 0..self.n {
            self.outputs[i] = self.x[i];
        }
    }

    fn step(&mut self, _t: f64, _dt: f64) -> StepResult {
        // If dt is None, updates are triggered manually via update_filter()
        // Otherwise, update is triggered automatically in step()
        if self.dt.is_some() {
            // Check if this step matches our discrete time step
            // For simplicity, we update on every step when dt is configured
            self.kf_update();
        }
        StepResult::default()
    }

    fn buffer(&mut self) {
        self.buffered_x = self.x.clone();
        self.buffered_P = self.P.clone();
    }

    fn revert(&mut self) {
        self.x = self.buffered_x.clone();
        self.P = self.buffered_P.clone();
        for i in 0..self.n {
            self.outputs[i] = self.x[i];
        }
    }

    fn reset(&mut self) {
        self.x = self.x0.clone();
        self.P = self.P0.clone();
        self.inputs.fill(0.0);
        for i in 0..self.n {
            self.outputs[i] = self.x[i];
        }
        self.buffered_x = self.x.clone();
        self.buffered_P = self.P.clone();
    }
}

impl DynamicBlock for KalmanFilter {
    fn state(&self) -> &[f64] {
        // Kalman filter doesn't integrate ODEs, so this is somewhat abstract
        // We return the state estimate as the "state"
        &self.outputs
    }

    fn state_derivative(&self) -> &[f64] {
        // Discrete-time system - no derivative
        &[]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_init_basic() {
        // Simple 2D system (position and velocity)
        let F = DMatrix::from_row_slice(2, 2, &[1.0, 0.1, 0.0, 1.0]);
        let H = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let Q = DMatrix::identity(2, 2) * 0.01;
        let R = DMatrix::from_row_slice(1, 1, &[0.1]);

        let kf = KalmanFilter::new(F.clone(), H.clone(), Q.clone(), R.clone());

        // Check dimensions
        assert_eq!(kf.n, 2); // state dimension
        assert_eq!(kf.m, 1); // measurement dimension
        assert_eq!(kf.p, 0); // no control input

        // Check default initial conditions
        assert_eq!(kf.x, DVector::zeros(2));
        assert_eq!(kf.P, DMatrix::identity(2, 2));

        // Check matrices
        assert_eq!(kf.F, F);
        assert_eq!(kf.H, H);
        assert_eq!(kf.Q, Q);
        assert_eq!(kf.R, R);
        assert!(kf.B.is_none());

        // Check io dimensions
        assert_eq!(kf.inputs.len(), 1); // m measurements
        assert_eq!(kf.outputs.len(), 2); // n states
    }

    #[test]
    fn test_init_with_initial_conditions() {
        let F = DMatrix::from_row_slice(2, 2, &[1.0, 0.1, 0.0, 1.0]);
        let H = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let Q = DMatrix::identity(2, 2) * 0.01;
        let R = DMatrix::from_row_slice(1, 1, &[0.1]);
        let x0 = DVector::from_row_slice(&[1.0, 2.0]);
        let P0 = DMatrix::identity(2, 2) * 5.0;

        let kf =
            KalmanFilter::with_options(F, H, Q, R, None, Some(x0.clone()), Some(P0.clone()), None);

        assert_eq!(kf.x, x0);
        assert_eq!(kf.P, P0);
    }

    #[test]
    fn test_init_with_control_input() {
        let F = DMatrix::from_row_slice(2, 2, &[1.0, 0.1, 0.0, 1.0]);
        let H = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let Q = DMatrix::identity(2, 2) * 0.01;
        let R = DMatrix::from_row_slice(1, 1, &[0.1]);
        let B = DMatrix::from_row_slice(2, 1, &[0.0, 0.1]);

        let kf = KalmanFilter::with_control(F, H, Q, R, B.clone());

        // Check dimensions
        assert_eq!(kf.n, 2);
        assert_eq!(kf.m, 1);
        assert_eq!(kf.p, 1); // one control input

        assert_eq!(kf.B, Some(B));

        // Check io dimensions (m measurements + p controls)
        assert_eq!(kf.inputs.len(), 2); // 1 measurement + 1 control
        assert_eq!(kf.outputs.len(), 2); // n states
    }

    #[test]
    fn test_init_with_dt() {
        let F = DMatrix::from_row_slice(2, 2, &[1.0, 0.1, 0.0, 1.0]);
        let H = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let Q = DMatrix::identity(2, 2) * 0.01;
        let R = DMatrix::from_row_slice(1, 1, &[0.1]);
        let dt = 0.1;

        let kf = KalmanFilter::with_options(F, H, Q, R, None, None, None, Some(dt));

        assert_eq!(kf.dt, Some(dt));
    }

    #[test]
    fn test_constant_position_estimation() {
        // System: stationary object at position x=5.0
        let F = DMatrix::from_row_slice(1, 1, &[1.0]);
        let H = DMatrix::from_row_slice(1, 1, &[1.0]);
        let Q = DMatrix::from_row_slice(1, 1, &[0.0]);
        let R = DMatrix::from_row_slice(1, 1, &[1.0]);

        // Start with uncertain initial estimate
        let x0 = DVector::from_row_slice(&[0.0]);
        let P0 = DMatrix::from_row_slice(1, 1, &[10.0]);

        let mut kf = KalmanFilter::with_options(F, H, Q, R, None, Some(x0), Some(P0.clone()), None);

        // Simulate measurements of true position = 5.0
        let true_position = 5.0;

        for _ in 0..10 {
            kf.set_input(0, true_position);
            kf.kf_update();
        }

        // After 10 measurements, estimate should be close to true value
        assert_relative_eq!(kf.x[0], true_position, epsilon = 0.5);

        // Covariance should decrease (more certain)
        assert!(kf.P[(0, 0)] < P0[(0, 0)]);
    }

    #[test]
    fn test_constant_velocity_estimation() {
        let dt = 0.1;

        // System: constant velocity model
        let F = DMatrix::from_row_slice(2, 2, &[1.0, dt, 0.0, 1.0]);
        let H = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let Q = DMatrix::from_vec(2, 2, vec![0.001, 0.0, 0.0, 0.001]);
        let R = DMatrix::from_row_slice(1, 1, &[0.5]);

        // True system state
        let mut true_position = 0.0;
        let true_velocity = 2.0;

        // Initial estimate (uncertain)
        let x0 = DVector::from_row_slice(&[0.0, 0.0]);
        let P0 = DMatrix::identity(2, 2) * 2.0;

        let mut kf = KalmanFilter::with_options(F, H, Q, R, None, Some(x0), Some(P0), None);

        // Simulate system over time
        let num_steps = 50;
        for _ in 0..num_steps {
            // True position at this time
            true_position += true_velocity * dt;

            // Noisy measurement
            let measurement = true_position;

            // Update filter
            kf.set_input(0, measurement);
            kf.kf_update();
        }

        // After convergence, velocity estimate should be close to true velocity
        assert_relative_eq!(kf.x[1], true_velocity, epsilon = 0.5);

        // Position should also track
        assert_relative_eq!(kf.x[0], true_position, epsilon = 1.0);
    }

    #[test]
    fn test_with_control_input() {
        let dt = 0.1;

        // System with control input (force applied to mass)
        let F = DMatrix::from_row_slice(2, 2, &[1.0, dt, 0.0, 1.0]);
        let H = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let B = DMatrix::from_row_slice(2, 1, &[0.5 * dt * dt, dt]);
        let Q = DMatrix::from_vec(2, 2, vec![0.01, 0.0, 0.0, 0.01]);
        let R = DMatrix::from_row_slice(1, 1, &[0.1]);

        let x0 = DVector::from_row_slice(&[0.0, 0.0]);
        let P0 = DMatrix::identity(2, 2);

        let mut kf = KalmanFilter::with_options(F, H, Q, R, Some(B), Some(x0), Some(P0), None);

        // Apply constant control input (acceleration = 1.0)
        let control = 1.0;
        let num_steps = 20;

        let mut position = 0.0;
        let mut velocity = 0.0;

        for _ in 0..num_steps {
            // True system evolution with control
            velocity += control * dt;
            position += velocity * dt;

            // Measurement
            let measurement = position;

            // Update filter with measurement and control
            kf.set_input(0, measurement); // measurement
            kf.set_input(1, control); // control input
            kf.kf_update();
        }

        // State estimates should track the true values
        assert_relative_eq!(kf.x[0], position, epsilon = 1.0);
        assert_relative_eq!(kf.x[1], velocity, epsilon = 0.5);
    }

    #[test]
    fn test_output_update() {
        let F = DMatrix::from_row_slice(2, 2, &[1.0, 0.1, 0.0, 1.0]);
        let H = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let Q = DMatrix::identity(2, 2) * 0.01;
        let R = DMatrix::from_row_slice(1, 1, &[0.1]);
        let x0 = DVector::from_row_slice(&[1.0, 2.0]);

        let mut kf = KalmanFilter::with_options(F, H, Q, R, None, Some(x0), None, None);

        // Set measurement
        kf.set_input(0, 3.0);

        // Perform update
        kf.kf_update();

        // Outputs should match state
        for i in 0..kf.n {
            assert_eq!(kf.outputs[i], kf.x[i]);
        }
    }

    #[test]
    fn test_multidimensional_measurement() {
        // 4D state: [x, y, vx, vy] - 2D position and velocity
        let dt = 0.1;
        let F = DMatrix::from_row_slice(
            4,
            4,
            &[
                1.0, 0.0, dt, 0.0, 0.0, 1.0, 0.0, dt, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ],
        );

        // Measure both x and y position
        let H = DMatrix::from_row_slice(2, 4, &[1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]);

        let Q = DMatrix::identity(4, 4) * 0.01;
        let R = DMatrix::identity(2, 2) * 0.1;
        let x0 = DVector::zeros(4);

        let mut kf = KalmanFilter::with_options(F, H, Q, R, None, Some(x0), None, None);

        // Check dimensions
        assert_eq!(kf.n, 4); // state dimension
        assert_eq!(kf.m, 2); // measurement dimension
        assert_eq!(kf.inputs.len(), 2); // 2 measurements
        assert_eq!(kf.outputs.len(), 4); // 4 states

        // Perform a few updates
        for i in 0..10 {
            kf.set_input(0, i as f64 * 0.1); // x measurement
            kf.set_input(1, i as f64 * 0.2); // y measurement
            kf.kf_update();
        }

        // Just check that it runs without error and produces reasonable output
        assert!(kf.x.iter().all(|&x| x.is_finite()));
        assert!(kf.P.iter().all(|&p| p.is_finite()));
    }

    #[test]
    fn test_reset() {
        let F = DMatrix::from_row_slice(2, 2, &[1.0, 0.1, 0.0, 1.0]);
        let H = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let Q = DMatrix::identity(2, 2) * 0.01;
        let R = DMatrix::from_row_slice(1, 1, &[0.1]);
        let x0 = DVector::from_row_slice(&[1.0, 2.0]);
        let P0 = DMatrix::identity(2, 2) * 3.0;

        let mut kf =
            KalmanFilter::with_options(F, H, Q, R, None, Some(x0.clone()), Some(P0.clone()), None);

        // Update filter
        kf.set_input(0, 5.0);
        kf.kf_update();

        // State should have changed
        assert_ne!(kf.x, x0);

        // Reset
        kf.reset();

        // Should be back to initial conditions
        assert_eq!(kf.x, x0);
        assert_eq!(kf.P, P0);
        assert_eq!(kf.inputs[0], 0.0);
    }

    #[test]
    fn test_manual_update_without_dt() {
        // Test that manual update works when dt is None
        let F = DMatrix::from_row_slice(1, 1, &[1.0]);
        let H = DMatrix::from_row_slice(1, 1, &[1.0]);
        let Q = DMatrix::from_row_slice(1, 1, &[0.01]);
        let R = DMatrix::from_row_slice(1, 1, &[0.1]);
        let x0 = DVector::from_row_slice(&[0.0]);

        let mut kf = KalmanFilter::with_options(
            F,
            H,
            Q,
            R,
            None,
            Some(x0),
            None,
            None, // dt = None
        );

        // Set measurement
        kf.set_input(0, 5.0);

        // Initial state
        let initial_state = kf.x[0];

        // Manually trigger update
        kf.update_filter();

        // State should have been updated
        assert_ne!(kf.x[0], initial_state);

        // Output should be updated
        assert_eq!(kf.get_output(0), kf.x[0]);
    }

    #[test]
    fn test_buffer_revert() {
        let F = DMatrix::from_row_slice(2, 2, &[1.0, 0.1, 0.0, 1.0]);
        let H = DMatrix::from_row_slice(1, 2, &[1.0, 0.0]);
        let Q = DMatrix::identity(2, 2) * 0.01;
        let R = DMatrix::from_row_slice(1, 1, &[0.1]);

        let mut kf = KalmanFilter::new(F, H, Q, R);

        // Do one update
        kf.set_input(0, 1.0);
        kf.kf_update();

        // Buffer state
        let x_buffered = kf.x.clone();
        let P_buffered = kf.P.clone();
        kf.buffer();

        // Do another update
        kf.set_input(0, 5.0);
        kf.kf_update();

        // State should have changed
        assert_ne!(kf.x, x_buffered);

        // Revert
        kf.revert();

        // Should be back to buffered state
        assert_eq!(kf.x, x_buffered);
        assert_eq!(kf.P, P_buffered);
    }
}
