//! Transfer Function block: Linear Time-Invariant (LTI) system in polynomial form
//!
//! Implements transfer functions defined by numerator and denominator polynomials:
//!   H(s) = B(s) / A(s) = (b_n s^n + b_{n-1} s^{n-1} + ... + b_0) / (a_m s^m + a_{m-1} s^{m-1} + ... + a_0)
//!
//! The transfer function is converted to state-space form using controllable canonical form
//! (also called controller canonical form), then integrated using RK4.
//!
//! References:
//! - Ogata, K. (2010). Modern Control Engineering (5th ed.). Section 5.6
//! - Chen, C.T. (1999). Linear System Theory and Design (3rd ed.). Section 5.5

use crate::block::{Block, DynamicBlock, StepResult};

/// Transfer Function: SISO LTI system in numerator-denominator polynomial form
///
/// The transfer function H(s) = B(s)/A(s) is converted to state-space form:
///   dx/dt = Ax + Bu
///   y = Cx + Du
///
/// where the realization uses controllable canonical form when n < m,
/// and handles proper (n = m) and improper (n > m) transfer functions appropriately.
///
/// # Type Parameters
///
/// - `N`: Order of the system (number of states = max(len(num)-1, len(den)-1))
/// - `M`: Number of inputs (always 1 for SISO)
/// - `P`: Number of outputs (always 1 for SISO)
///
/// # Polynomial Convention
///
/// Polynomials are specified in **descending powers** of s:
/// - `num = [b_n, b_{n-1}, ..., b_1, b_0]` represents `b_n*s^n + b_{n-1}*s^{n-1} + ... + b_0`
/// - `den = [a_m, a_{m-1}, ..., a_1, a_0]` represents `a_m*s^m + a_{m-1}*s^{m-1} + ... + a_0`
///
/// The leading coefficient of the denominator (a_m) is normalized to 1.
///
/// # Example (First-order system: 1/(s+1))
///
/// ```rust
/// use rustsim::blocks::TransferFunction;
/// use rustsim::Block;
///
/// // H(s) = 1/(s+1)
/// // num = [1.0], den = [1.0, 1.0] represents 1/(1*s + 1)
/// let mut tf = TransferFunction::<1>::new(&[1.0], &[1.0, 1.0]);
///
/// // Apply step input
/// tf.set_input(0, 1.0);
/// tf.update(0.0);
///
/// // Simulate
/// let dt = 0.01;
/// for _ in 0..100 {
///     tf.step(0.0, dt);
/// }
///
/// // After ~5 time constants, output should approach 1.0
/// assert!((tf.get_output(0) - 0.993).abs() < 0.01);
/// ```
///
/// # Example (Second-order system)
///
/// ```rust
/// use rustsim::blocks::TransferFunction;
/// use rustsim::Block;
///
/// // H(s) = 1/(s^2 + 2s + 1) = 1/(s+1)^2
/// // num = [1.0], den = [1.0, 2.0, 1.0]
/// let mut tf = TransferFunction::<2>::new(&[1.0], &[1.0, 2.0, 1.0]);
///
/// tf.set_input(0, 1.0);
/// tf.update(0.0);
///
/// let dt = 0.01;
/// for _ in 0..500 {
///     tf.step(0.0, dt);
/// }
///
/// // Steady state output should be 1.0
/// assert!((tf.get_output(0) - 1.0).abs() < 0.1);
/// ```
///
/// # Example (Transfer function with zeros: (s+1)/(s+2))
///
/// ```rust
/// use rustsim::blocks::TransferFunction;
/// use rustsim::Block;
///
/// // H(s) = (s+1)/(s+2)
/// // num = [1.0, 1.0], den = [1.0, 2.0]
/// let mut tf = TransferFunction::<1>::new(&[1.0, 1.0], &[1.0, 2.0]);
///
/// tf.set_input(0, 1.0);
/// tf.update(0.0);
///
/// // This has direct feedthrough D = b0/a0 = 1/1 = 1.0 when degrees are equal
/// assert_eq!(tf.get_output(0), 0.5); // (1 + 0)/(2 + 0) with initial state 0
/// ```
#[derive(Clone)]
pub struct TransferFunction<const N: usize> {
    /// State matrix (N×N) from controllable canonical form
    a: [[f64; N]; N],
    /// Input matrix (N×1)
    b: [[f64; 1]; N],
    /// Output matrix (1×N)
    c: [[f64; N]; 1],
    /// Feedthrough matrix (1×1)
    d: [[f64; 1]; 1],

    /// External inputs (1-dimensional for SISO)
    inputs: [f64; 1],
    /// System outputs (1-dimensional for SISO)
    outputs: [f64; 1],
    /// Current state (N-dimensional)
    state: [f64; N],
    /// State derivative (computed during update)
    derivative: [f64; N],
    /// Initial state
    initial: [f64; N],

    /// Whether D has non-zero element (direct feedthrough)
    has_passthrough: bool,

    // For buffer/revert
    buffered_state: [f64; N],

    // RK4 intermediate values
    k1: [f64; N],
    k2: [f64; N],
    k3: [f64; N],
    k4: [f64; N],
    temp_state: [f64; N],
}

impl<const N: usize> TransferFunction<N> {
    /// Create Transfer Function block from numerator and denominator polynomials
    ///
    /// # Arguments
    ///
    /// * `num` - Numerator polynomial coefficients in descending powers: [b_n, b_{n-1}, ..., b_0]
    /// * `den` - Denominator polynomial coefficients in descending powers: [a_m, a_{m-1}, ..., a_0]
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - Denominator is empty or all zeros
    /// - Leading coefficient of denominator is zero
    /// - Numerator order exceeds denominator order (improper transfer function)
    ///
    /// # State-Space Realization
    ///
    /// Uses observable canonical form (OCF), matching scipy.signal.TransferFunction.to_ss().
    /// For a transfer function:
    ///   H(s) = (b_n*s^n + ... + b_0) / (s^n + a_{n-1}*s^{n-1} + ... + a_0)
    ///
    /// The state-space realization is:
    ///   A = [-a_{n-1}  -a_{n-2}  ...  -a_1  -a_0 ]
    ///       [   1         0      ...   0     0   ]
    ///       [   0         1      ...   0     0   ]
    ///       [   ⋮         ⋮      ⋱     ⋮     ⋮   ]
    ///       [   0         0      ...   1     0   ]
    ///
    ///   B = [1]
    ///       [0]
    ///       [⋮]
    ///       [0]
    ///
    ///   C = [b_{n-1}  b_{n-2}  ...  b_1  b_0]  (after removing D component)
    ///   D = b_n  (if numerator degree equals denominator degree)
    pub fn new(num: &[f64], den: &[f64]) -> Self {
        assert!(!den.is_empty(), "Denominator cannot be empty");
        assert!(
            den[0] != 0.0,
            "Leading coefficient of denominator cannot be zero"
        );
        assert!(
            num.len() <= den.len(),
            "Improper transfer function (numerator degree > denominator degree) not supported. Got num.len()={}, den.len()={}",
            num.len(),
            den.len()
        );

        // Normalize denominator by leading coefficient
        let leading_coeff = den[0];
        let den_norm: Vec<f64> = den.iter().map(|&x| x / leading_coeff).collect();
        let mut num_norm: Vec<f64> = num.iter().map(|&x| x / leading_coeff).collect();

        // Pad numerator with leading zeros if necessary to match denominator length
        while num_norm.len() < den_norm.len() {
            num_norm.insert(0, 0.0);
        }

        let order = den_norm.len() - 1; // System order is one less than coefficient count

        assert_eq!(
            order,
            N,
            "Transfer function order {} does not match generic parameter N={}. \
             Denominator has {} coefficients (order = len-1 = {}). \
             Use TransferFunction::<{}>",
            order,
            N,
            den_norm.len(),
            order,
            order
        );

        // Extract direct feedthrough (D matrix)
        // If numerator and denominator have same degree, D = b_n/a_n (already normalized, so = b_n)
        let d_value = if num_norm.len() == den_norm.len() {
            num_norm[0]
        } else {
            0.0
        };

        // Compute numerator coefficients for C matrix (after removing D component)
        // If D != 0, we need to compute num - D*den, then take coefficients
        let mut c_coeffs = vec![0.0; N];

        if d_value != 0.0 {
            // Subtract D*den from num to get the strictly proper part
            let mut num_proper = num_norm.clone();
            for i in 0..den_norm.len() {
                num_proper[i] -= d_value * den_norm[i];
            }
            // Take the last N coefficients
            for i in 0..N {
                c_coeffs[i] = num_proper[num_proper.len() - N + i];
            }
        } else {
            // Take last N coefficients directly
            for i in 0..N {
                if i < num_norm.len() {
                    c_coeffs[i] = num_norm[num_norm.len() - N + i];
                }
            }
        }

        // Build state-space matrices in OBSERVABLE canonical form (matching scipy)
        // This is the transpose of the controllable canonical form
        let mut a = [[0.0; N]; N];
        let mut b = [[0.0; 1]; N];
        let mut c = [[0.0; N]; 1];
        let d = [[d_value; 1]; 1];

        if N > 0 {
            // A matrix for observable canonical form:
            // First row: -a_{n-1}, -a_{n-2}, ..., -a_1, -a_0
            // Identity subdiagonal: a[i+1][i] = 1 for i = 0..N-2
            for j in 0..N {
                a[0][j] = -den_norm[j + 1]; // Skip leading 1.0, take coefficients in order
            }
            for i in 0..N - 1 {
                a[i + 1][i] = 1.0;
            }

            // B matrix: [1, 0, 0, ..., 0]^T
            b[0][0] = 1.0;

            // C matrix: matches scipy's observable canonical form
            // C = [c_coeffs[0], c_coeffs[1], ..., c_coeffs[N-1]]
            for i in 0..N {
                c[0][i] = c_coeffs[i];
            }
        }

        let has_passthrough = d_value != 0.0;
        let initial = [0.0; N];

        // Compute initial output: y = Cx + Du (where x = 0, u = 0 initially)
        let outputs = [0.0; 1];

        Self {
            a,
            b,
            c,
            d,
            inputs: [0.0; 1],
            outputs,
            state: initial,
            derivative: [0.0; N],
            initial,
            has_passthrough,
            buffered_state: initial,
            k1: [0.0; N],
            k2: [0.0; N],
            k3: [0.0; N],
            k4: [0.0; N],
            temp_state: [0.0; N],
        }
    }

    /// Get current state value at index
    pub fn state_value(&self, index: usize) -> f64 {
        self.state[index]
    }

    /// Get entire state vector
    pub fn state_vector(&self) -> &[f64; N] {
        &self.state
    }

    /// Get the A matrix (state matrix)
    pub fn a_matrix(&self) -> &[[f64; N]; N] {
        &self.a
    }

    /// Get the B matrix (input matrix)
    pub fn b_matrix(&self) -> &[[f64; 1]; N] {
        &self.b
    }

    /// Get the C matrix (output matrix)
    pub fn c_matrix(&self) -> &[[f64; N]; 1] {
        &self.c
    }

    /// Get the D matrix (feedthrough matrix)
    pub fn d_matrix(&self) -> &[[f64; 1]; 1] {
        &self.d
    }

    /// Returns true if there is direct feedthrough (D != 0)
    pub fn has_passthrough(&self) -> bool {
        self.has_passthrough
    }

    /// Compute state derivative: dx/dt = Ax + Bu
    fn compute_derivative(&self, state: &[f64; N], inputs: &[f64; 1]) -> [f64; N] {
        let mut deriv = [0.0; N];

        // Ax
        for i in 0..N {
            for j in 0..N {
                deriv[i] += self.a[i][j] * state[j];
            }
        }

        // Bu
        for i in 0..N {
            deriv[i] += self.b[i][0] * inputs[0];
        }

        deriv
    }

    /// Compute output: y = Cx + Du
    fn compute_output(&self, state: &[f64; N], inputs: &[f64; 1]) -> [f64; 1] {
        let mut output = [0.0; 1];

        // Cx
        for j in 0..N {
            output[0] += self.c[0][j] * state[j];
        }

        // Du
        output[0] += self.d[0][0] * inputs[0];

        output
    }
}

impl<const N: usize> Block for TransferFunction<N> {
    const NUM_INPUTS: usize = 1;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = true;

    #[inline]
    fn inputs(&self) -> &[f64] {
        &self.inputs
    }

    #[inline]
    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut self.inputs
    }

    #[inline]
    fn outputs(&self) -> &[f64] {
        &self.outputs
    }

    #[inline]
    fn outputs_mut(&mut self) -> &mut [f64] {
        &mut self.outputs
    }

    fn update(&mut self, _t: f64) {
        // Compute output: y = Cx + Du
        self.outputs = self.compute_output(&self.state, &self.inputs);

        // Compute derivative for step: dx/dt = Ax + Bu
        self.derivative = self.compute_derivative(&self.state, &self.inputs);
    }

    fn step(&mut self, _t: f64, dt: f64) -> StepResult {
        // RK4 integration

        // k1 = f(t, y)
        self.k1 = self.compute_derivative(&self.state, &self.inputs);

        // k2 = f(t + dt/2, y + dt/2 * k1)
        for i in 0..N {
            self.temp_state[i] = self.state[i] + 0.5 * dt * self.k1[i];
        }
        self.k2 = self.compute_derivative(&self.temp_state, &self.inputs);

        // k3 = f(t + dt/2, y + dt/2 * k2)
        for i in 0..N {
            self.temp_state[i] = self.state[i] + 0.5 * dt * self.k2[i];
        }
        self.k3 = self.compute_derivative(&self.temp_state, &self.inputs);

        // k4 = f(t + dt, y + dt * k3)
        for i in 0..N {
            self.temp_state[i] = self.state[i] + dt * self.k3[i];
        }
        self.k4 = self.compute_derivative(&self.temp_state, &self.inputs);

        // y_new = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        for i in 0..N {
            self.state[i] +=
                dt / 6.0 * (self.k1[i] + 2.0 * self.k2[i] + 2.0 * self.k3[i] + self.k4[i]);
        }

        // Update output with new state
        self.outputs = self.compute_output(&self.state, &self.inputs);

        StepResult::default()
    }

    fn buffer(&mut self) {
        self.buffered_state.copy_from_slice(&self.state);
    }

    fn revert(&mut self) {
        self.state.copy_from_slice(&self.buffered_state);
        self.outputs = self.compute_output(&self.state, &self.inputs);
    }

    fn reset(&mut self) {
        self.inputs.fill(0.0);
        self.state.copy_from_slice(&self.initial);
        self.outputs = self.compute_output(&self.state, &self.inputs);
        self.derivative.fill(0.0);
        self.buffered_state.copy_from_slice(&self.initial);
    }
}

impl<const N: usize> DynamicBlock for TransferFunction<N> {
    fn state(&self) -> &[f64] {
        &self.state
    }

    fn state_derivative(&self) -> &[f64] {
        &self.derivative
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tf_first_order() {
        // H(s) = 1/(s+1)
        // num = [1], den = [1, 1]
        let tf = TransferFunction::<1>::new(&[1.0], &[1.0, 1.0]);

        // Check state-space matrices
        assert_eq!(tf.a_matrix()[0][0], -1.0); // A = [-1]
        assert_eq!(tf.b_matrix()[0][0], 1.0); // B = [1]
        assert_eq!(tf.c_matrix()[0][0], 1.0); // C = [1]
        assert_eq!(tf.d_matrix()[0][0], 0.0); // D = 0
        assert!(!tf.has_passthrough());
    }

    #[test]
    fn test_tf_second_order() {
        // H(s) = 1/(s^2 + 2s + 1)
        // num = [1], den = [1, 2, 1]
        let tf = TransferFunction::<2>::new(&[1.0], &[1.0, 2.0, 1.0]);

        // Observable canonical form (matching scipy):
        // A = [-2  -1]
        //     [ 1   0]
        assert_eq!(tf.a_matrix()[0][0], -2.0);
        assert_eq!(tf.a_matrix()[0][1], -1.0);
        assert_eq!(tf.a_matrix()[1][0], 1.0);
        assert_eq!(tf.a_matrix()[1][1], 0.0);

        // B = [1]
        //     [0]
        assert_eq!(tf.b_matrix()[0][0], 1.0);
        assert_eq!(tf.b_matrix()[1][0], 0.0);

        // For H(s) = 1/(s^2 + 2s + 1), C = [0, 1]
        assert_eq!(tf.c_matrix()[0][0], 0.0);
        assert_eq!(tf.c_matrix()[0][1], 1.0);

        assert_eq!(tf.d_matrix()[0][0], 0.0);
    }

    #[test]
    fn test_tf_with_zeros() {
        // H(s) = (s+1)/(s+2)
        // num = [1, 1], den = [1, 2]
        let tf = TransferFunction::<1>::new(&[1.0, 1.0], &[1.0, 2.0]);

        // Same degree -> has direct feedthrough
        // D = 1/1 = 1
        // H(s) = 1 + something/(s+2)
        assert!(tf.has_passthrough());
        assert_eq!(tf.d_matrix()[0][0], 1.0);

        // After removing D, we have (s+1) - 1*(s+2) = -1
        // So strictly proper part is -1/(s+2)
        // A = [-2], B = [1], C = [-1], D = 1
        assert_eq!(tf.a_matrix()[0][0], -2.0);
        assert_eq!(tf.b_matrix()[0][0], 1.0);
        assert_eq!(tf.c_matrix()[0][0], -1.0);
    }

    #[test]
    fn test_tf_step_response_first_order() {
        // H(s) = 1/(s+1)
        let mut tf = TransferFunction::<1>::new(&[1.0], &[1.0, 1.0]);

        // Apply step input
        tf.set_input(0, 1.0);
        tf.update(0.0);

        let dt = 0.01;
        let t_final = 5.0; // 5 time constants
        let steps = (t_final / dt) as usize;

        for _ in 0..steps {
            tf.step(0.0, dt);
        }

        // Steady state should be 1.0
        // After 5 time constants, should be ~99.3% of steady state
        let expected = 1.0 * (1.0 - (-5.0_f64).exp());
        assert!(
            (tf.get_output(0) - expected).abs() < 1e-2,
            "Output: {}, Expected: {}",
            tf.get_output(0),
            expected
        );
    }

    #[test]
    fn test_tf_normalization() {
        // H(s) = 2/(2s+2) = 1/(s+1)
        let tf = TransferFunction::<1>::new(&[2.0], &[2.0, 2.0]);

        // After normalization, should be same as 1/(s+1)
        assert_eq!(tf.a_matrix()[0][0], -1.0);
        assert_eq!(tf.b_matrix()[0][0], 1.0);
        assert_eq!(tf.c_matrix()[0][0], 1.0);
        assert_eq!(tf.d_matrix()[0][0], 0.0);
    }

    #[test]
    #[should_panic(expected = "Denominator cannot be empty")]
    fn test_tf_empty_denominator() {
        let _tf = TransferFunction::<1>::new(&[1.0], &[]);
    }

    #[test]
    #[should_panic(expected = "Leading coefficient of denominator cannot be zero")]
    fn test_tf_zero_leading_coeff() {
        let _tf = TransferFunction::<1>::new(&[1.0], &[0.0, 1.0]);
    }

    #[test]
    #[should_panic(expected = "Improper transfer function")]
    fn test_tf_improper() {
        // num has higher degree than den
        let _tf = TransferFunction::<2>::new(&[1.0, 2.0, 3.0], &[1.0, 1.0]);
    }

    #[test]
    fn test_tf_reset() {
        let mut tf = TransferFunction::<2>::new(&[1.0], &[1.0, 2.0, 1.0]);

        tf.set_input(0, 1.0);
        tf.update(0.0);
        tf.step(0.0, 0.1);

        // State should have changed
        assert!(tf.state_value(0) != 0.0 || tf.state_value(1) != 0.0);

        tf.reset();

        // State should be back to zero
        assert_eq!(tf.state_value(0), 0.0);
        assert_eq!(tf.state_value(1), 0.0);
        assert_eq!(tf.get_input(0), 0.0);
        assert_eq!(tf.get_output(0), 0.0);
    }

    #[test]
    fn test_tf_buffer_revert() {
        let mut tf = TransferFunction::<1>::new(&[1.0], &[1.0, 1.0]);

        tf.set_input(0, 1.0);
        tf.update(0.0);
        tf.step(0.0, 0.1);

        let state_after_step = tf.state_value(0);
        tf.buffer();

        // Take another step
        tf.step(0.0, 0.1);
        assert!(tf.state_value(0) != state_after_step);

        // Revert
        tf.revert();
        assert_eq!(tf.state_value(0), state_after_step);
    }
}
