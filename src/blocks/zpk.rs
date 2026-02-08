//! Zero-Pole-Gain (ZPK) Transfer Function Block
//!
//! Implements transfer functions in zeros-poles-gain form:
//!
//! ```text
//! H(s) = K * (s - z₁)(s - z₂)...(s - zₘ) / (s - p₁)(s - p₂)...(s - pₙ)
//! ```

use crate::block::{Block, DynamicBlock, StepResult};
use num_complex::Complex64;
use std::fmt;

/// Zero-Pole-Gain Transfer Function Block
#[derive(Clone)]
pub struct ZerosPoleGain<const N: usize> {
    /// State-space A matrix (N×N)
    a: [[f64; N]; N],
    /// State-space B matrix (N×1)
    b: [[f64; 1]; N],
    /// State-space C matrix (1×N)
    c: [[f64; N]; 1],
    /// State-space D matrix (1×1) - direct feedthrough
    d: [[f64; 1]; 1],

    /// System poles (complex)
    poles: [Complex64; N],
    /// System zeros (complex, can be empty)
    zeros: Vec<Complex64>,
    /// System gain
    gain: f64,

    /// External input
    input: f64,
    /// System output
    output: f64,
    /// Current state (N-dimensional)
    state: [f64; N],
    /// State derivative (computed during update)
    derivative: [f64; N],
    /// Initial state
    initial: [f64; N],

    /// Whether D has non-zero element (direct passthrough)
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

impl<const N: usize> ZerosPoleGain<N> {
    /// Create a new ZPK transfer function block
    pub fn new(zeros: Vec<Complex64>, poles: [Complex64; N], gain: f64) -> Self {
        // Validate that complex poles come in conjugate pairs
        Self::validate_conjugate_pairs(&poles, "poles");
        Self::validate_conjugate_pairs(&zeros, "zeros");

        // Convert ZPK to state-space using controllable canonical form
        let (a, b, c, d) = Self::zpk_to_statespace(&zeros, &poles, gain);

        let has_passthrough = d[0][0].abs() > 1e-15;

        let initial = [0.0; N];

        Self {
            a,
            b,
            c,
            d,
            poles,
            zeros,
            gain,
            input: 0.0,
            output: 0.0,
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

    /// Create from real-valued poles only (convenience constructor)
    pub fn new_real(zeros: Vec<f64>, poles: [f64; N], gain: f64) -> Self {
        let zeros_complex: Vec<Complex64> = zeros.iter().map(|&z| Complex64::new(z, 0.0)).collect();
        let poles_complex: [Complex64; N] =
            poles.map(|p| Complex64::new(p, 0.0));

        Self::new(zeros_complex, poles_complex, gain)
    }

    /// Validate that complex numbers come in conjugate pairs
    fn validate_conjugate_pairs(values: &[Complex64], name: &str) {
        for val in values {
            if val.im.abs() > 1e-15 {
                // This is complex, check if conjugate exists
                let conj = val.conj();
                let has_conjugate = values.iter().any(|v| {
                    (v.re - conj.re).abs() < 1e-12 && (v.im - conj.im).abs() < 1e-12
                });
                assert!(
                    has_conjugate,
                    "Complex {} must come in conjugate pairs: {} has no conjugate",
                    name, val
                );
            }
        }
    }

    /// Convert ZPK to state-space representation
    fn zpk_to_statespace(
        zeros: &[Complex64],
        poles: &[Complex64],
        gain: f64,
    ) -> ([[f64; N]; N], [[f64; 1]; N], [[f64; N]; 1], [[f64; 1]; 1]) {
        if N == 0 {
            // Pure gain system (no dynamics)
            return ([[0.0; N]; N], [[0.0; 1]; N], [[0.0; N]; 1], [[gain; 1]; 1]);
        }

        // Expand denominator polynomial from poles: D(s) = (s - p₁)(s - p₂)...(s - pₙ)
        let den_coeffs = Self::expand_polynomial(poles);

        // Expand numerator polynomial from zeros: N(s) = K * (s - z₁)(s - z₂)...(s - zₘ)
        let mut num_coeffs = if zeros.is_empty() {
            vec![Complex64::new(1.0, 0.0)]
        } else {
            Self::expand_polynomial(zeros)
        };

        // Multiply numerator by gain
        for coeff in &mut num_coeffs {
            *coeff *= gain;
        }

        // Normalize denominator (make leading coefficient 1)
        let leading = den_coeffs[den_coeffs.len() - 1];
        assert!(leading.norm() > 1e-15, "Leading denominator coefficient is zero");

        let den_norm: Vec<f64> = den_coeffs.iter().map(|c| (c / leading).re).collect();
        let num_norm: Vec<f64> = num_coeffs.iter().map(|c| (c / leading).re).collect();

        // Build state-space matrices in controllable canonical form
        let mut a = [[0.0; N]; N];
        let mut b = [[0.0; 1]; N];
        let mut c = [[0.0; N]; 1];
        let mut d = [[0.0; 1]; 1];

        // Build A matrix (companion matrix)
        // Upper diagonal: 1's
        for i in 0..N - 1 {
            a[i][i + 1] = 1.0;
        }
        // Bottom row: -a₀, -a₁, ..., -aₙ₋₁
        for i in 0..N {
            a[N - 1][i] = -den_norm[i];
        }

        // Build B matrix: [0, 0, ..., 0, 1]ᵀ
        b[N - 1][0] = 1.0;

        // Build C and D matrices from numerator
        if num_norm.len() > N {
            // Improper transfer function (more zeros than poles)
            d[0][0] = num_norm[num_norm.len() - 1];

            // Remainder becomes C
            for i in 0..N.min(num_norm.len() - 1) {
                c[0][i] = num_norm[i];
            }
        } else {
            // Proper or strictly proper
            d[0][0] = if num_norm.len() == N + 1 {
                num_norm[N]
            } else {
                0.0
            };

            // C coefficients
            for i in 0..num_norm.len().min(N) {
                c[0][i] = num_norm[i];
            }
        }

        (a, b, c, d)
    }

    /// Expand polynomial from roots: (s - r₁)(s - r₂)...(s - rₙ)
    fn expand_polynomial(roots: &[Complex64]) -> Vec<Complex64> {
        if roots.is_empty() {
            return vec![Complex64::new(1.0, 0.0)];
        }

        // Start with (s - r₀) = -r₀ + s
        let mut coeffs = vec![
            -roots[0],
            Complex64::new(1.0, 0.0),
        ];

        // Multiply by each remaining (s - rᵢ)
        for &root in &roots[1..] {
            let mut new_coeffs = vec![Complex64::new(0.0, 0.0); coeffs.len() + 1];

            // Multiply by s: shifts coefficients up
            for (i, &coeff) in coeffs.iter().enumerate() {
                new_coeffs[i + 1] += coeff;
            }

            // Multiply by -root: multiply all coefficients by -root
            for (i, &coeff) in coeffs.iter().enumerate() {
                new_coeffs[i] -= root * coeff;
            }

            coeffs = new_coeffs;
        }

        // Clean up imaginary parts that should be zero (due to conjugate pairs)
        for coeff in &mut coeffs {
            if coeff.im.abs() < 1e-12 {
                coeff.im = 0.0;
            }
        }

        coeffs
    }

    /// Compute state derivative: dx/dt = Ax + Bu
    fn compute_derivative(&self, state: &[f64; N], input: f64) -> [f64; N] {
        let mut deriv = [0.0; N];

        // Ax
        for i in 0..N {
            for j in 0..N {
                deriv[i] += self.a[i][j] * state[j];
            }
        }

        // Bu
        for i in 0..N {
            deriv[i] += self.b[i][0] * input;
        }

        deriv
    }

    /// Compute output: y = Cx + Du
    fn compute_output(&self, state: &[f64; N], input: f64) -> f64 {
        let mut output = 0.0;

        // Cx
        for j in 0..N {
            output += self.c[0][j] * state[j];
        }

        // Du
        output += self.d[0][0] * input;

        output
    }

    /// Get the system poles
    pub fn poles(&self) -> &[Complex64; N] {
        &self.poles
    }

    /// Get the system zeros
    pub fn zeros(&self) -> &[Complex64] {
        &self.zeros
    }

    /// Get the system gain
    pub fn gain(&self) -> f64 {
        self.gain
    }

    /// Returns true if the system has direct feedthrough (D ≠ 0)
    pub fn has_passthrough(&self) -> bool {
        self.has_passthrough
    }

    /// Get current state value at index
    pub fn state_value(&self, index: usize) -> f64 {
        self.state[index]
    }

    /// Get the state-space matrices (for inspection/testing)
    pub fn state_space_matrices(&self) -> (&[[f64; N]; N], &[[f64; 1]; N], &[[f64; N]; 1], &[[f64; 1]; 1]) {
        (&self.a, &self.b, &self.c, &self.d)
    }
}

impl<const N: usize> Block for ZerosPoleGain<N> {
    const NUM_INPUTS: usize = 1;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = N > 0;

    #[inline]
    fn inputs(&self) -> &[f64] {
        std::slice::from_ref(&self.input)
    }

    #[inline]
    fn inputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.input)
    }

    #[inline]
    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    #[inline]
    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, _t: f64) {
        if N == 0 {
            // Pure gain, no dynamics
            self.output = self.d[0][0] * self.input;
        } else {
            // Compute output: y = Cx + Du
            self.output = self.compute_output(&self.state, self.input);

            // Compute derivative for step: dx/dt = Ax + Bu
            self.derivative = self.compute_derivative(&self.state, self.input);
        }
    }

    fn step(&mut self, _t: f64, dt: f64) -> StepResult {
        if N == 0 {
            // No dynamics to integrate
            return StepResult::default();
        }

        // RK4 integration
        // k1 = f(t, y)
        self.k1 = self.compute_derivative(&self.state, self.input);

        // k2 = f(t + dt/2, y + dt/2 * k1)
        for i in 0..N {
            self.temp_state[i] = self.state[i] + 0.5 * dt * self.k1[i];
        }
        self.k2 = self.compute_derivative(&self.temp_state, self.input);

        // k3 = f(t + dt/2, y + dt/2 * k2)
        for i in 0..N {
            self.temp_state[i] = self.state[i] + 0.5 * dt * self.k2[i];
        }
        self.k3 = self.compute_derivative(&self.temp_state, self.input);

        // k4 = f(t + dt, y + dt * k3)
        for i in 0..N {
            self.temp_state[i] = self.state[i] + dt * self.k3[i];
        }
        self.k4 = self.compute_derivative(&self.temp_state, self.input);

        // y_new = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        for i in 0..N {
            self.state[i] += dt / 6.0 * (self.k1[i] + 2.0 * self.k2[i] + 2.0 * self.k3[i] + self.k4[i]);
        }

        // Update output with new state
        self.output = self.compute_output(&self.state, self.input);

        StepResult::default()
    }

    fn buffer(&mut self) {
        self.buffered_state.copy_from_slice(&self.state);
    }

    fn revert(&mut self) {
        self.state.copy_from_slice(&self.buffered_state);
        self.output = self.compute_output(&self.state, self.input);
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.state.copy_from_slice(&self.initial);
        self.output = self.compute_output(&self.state, self.input);
        self.derivative.fill(0.0);
        self.buffered_state.copy_from_slice(&self.initial);
    }
}

impl<const N: usize> DynamicBlock for ZerosPoleGain<N> {
    fn state(&self) -> &[f64] {
        &self.state
    }

    fn state_derivative(&self) -> &[f64] {
        &self.derivative
    }
}

impl<const N: usize> fmt::Debug for ZerosPoleGain<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct(&format!("ZerosPoleGain<{}>", N))
            .field("poles", &self.poles)
            .field("zeros", &self.zeros)
            .field("gain", &self.gain)
            .field("has_passthrough", &self.has_passthrough)
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_zpk_first_order_real() {
        // H(s) = 1 / (s + 1)
        let zpk = ZerosPoleGain::<1>::new_real(vec![], [-1.0], 1.0);

        assert_eq!(zpk.poles().len(), 1);
        assert_eq!(zpk.zeros().len(), 0);
        assert_eq!(zpk.gain(), 1.0);
        assert!(!zpk.has_passthrough());
    }

    #[test]
    fn test_zpk_expand_polynomial_simple() {
        // (s - 1) = s - 1
        let roots = vec![Complex64::new(1.0, 0.0)];
        let coeffs = ZerosPoleGain::<1>::expand_polynomial(&roots);

        assert_eq!(coeffs.len(), 2);
        assert_relative_eq!(coeffs[0].re, -1.0, epsilon = 1e-10);  // constant term
        assert_relative_eq!(coeffs[1].re, 1.0, epsilon = 1e-10);   // s term
    }

    #[test]
    fn test_zpk_expand_polynomial_quadratic() {
        // (s - 1)(s - 2) = s² - 3s + 2
        let roots = vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(2.0, 0.0),
        ];
        let coeffs = ZerosPoleGain::<2>::expand_polynomial(&roots);

        assert_eq!(coeffs.len(), 3);
        assert_relative_eq!(coeffs[0].re, 2.0, epsilon = 1e-10);   // constant
        assert_relative_eq!(coeffs[1].re, -3.0, epsilon = 1e-10);  // s
        assert_relative_eq!(coeffs[2].re, 1.0, epsilon = 1e-10);   // s²
    }

    #[test]
    fn test_zpk_expand_polynomial_complex_pair() {
        // (s - (1+i))(s - (1-i)) = s² - 2s + 2
        let roots = vec![
            Complex64::new(1.0, 1.0),
            Complex64::new(1.0, -1.0),
        ];
        let coeffs = ZerosPoleGain::<2>::expand_polynomial(&roots);

        assert_eq!(coeffs.len(), 3);
        assert_relative_eq!(coeffs[0].re, 2.0, epsilon = 1e-10);   // constant: (1+i)(1-i) = 2
        assert_relative_eq!(coeffs[0].im, 0.0, epsilon = 1e-10);
        assert_relative_eq!(coeffs[1].re, -2.0, epsilon = 1e-10);  // s: -(1+i + 1-i) = -2
        assert_relative_eq!(coeffs[1].im, 0.0, epsilon = 1e-10);
        assert_relative_eq!(coeffs[2].re, 1.0, epsilon = 1e-10);   // s²
        assert_relative_eq!(coeffs[2].im, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_zpk_state_space_conversion_first_order() {
        // H(s) = 3 / (s + 2) with pole at s=-2
        // Denominator is s+2, so controllable form gives A=[-2]
        let zpk = ZerosPoleGain::<1>::new_real(vec![], [-2.0], 3.0);
        let (a, b, c, d) = zpk.state_space_matrices();

        assert_relative_eq!(a[0][0], -2.0, epsilon = 1e-10);
        assert_relative_eq!(b[0][0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(c[0][0], 3.0, epsilon = 1e-10);
        assert_relative_eq!(d[0][0], 0.0, epsilon = 1e-10);
    }

    #[test]
    #[should_panic(expected = "must come in conjugate pairs")]
    fn test_zpk_validates_complex_poles() {
        // Single complex pole without conjugate should panic
        let _ = ZerosPoleGain::<1>::new(
            vec![],
            [Complex64::new(-1.0, 1.0)],
            1.0,
        );
    }

    #[test]
    fn test_zpk_complex_conjugate_poles() {
        // H(s) = 1 / ((s + 1 + i)(s + 1 - i)) = 1 / (s² + 2s + 2)
        let zpk = ZerosPoleGain::<2>::new(
            vec![],
            [Complex64::new(-1.0, 1.0), Complex64::new(-1.0, -1.0)],
            1.0,
        );

        let (a, b, c, d) = zpk.state_space_matrices();

        // For s² + 2s + 2, controllable form:
        // A = [0   1  ]
        //     [-2 -2 ]
        assert_relative_eq!(a[0][0], 0.0, epsilon = 1e-10);
        assert_relative_eq!(a[0][1], 1.0, epsilon = 1e-10);
        assert_relative_eq!(a[1][0], -2.0, epsilon = 1e-10);
        assert_relative_eq!(a[1][1], -2.0, epsilon = 1e-10);

        assert_relative_eq!(b[0][0], 0.0, epsilon = 1e-10);
        assert_relative_eq!(b[1][0], 1.0, epsilon = 1e-10);

        assert_relative_eq!(d[0][0], 0.0, epsilon = 1e-10);
    }
}
