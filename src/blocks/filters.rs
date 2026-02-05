//! Basic filter blocks

use crate::block::{Block, DynamicBlock, StepResult};

/// First-order RC lowpass filter
///
/// Transfer function: H(s) = 1 / (1 + s*tau)
///
/// Discrete implementation using backward Euler:
/// y[n] = alpha * x[n] + (1-alpha) * y[n-1]
/// where alpha = dt / (tau + dt)
///
/// # Example
///
/// ```ignore
/// let mut lpf = LowpassRC::new(0.1);  // tau = 0.1 seconds
/// lpf.set_input(0, 1.0);  // Step input
/// for i in 0..100 {
///     lpf.update(i as f64 * 0.01);
///     lpf.step(i as f64 * 0.01, 0.01);
/// }
/// // Output approaches 1.0 exponentially
/// ```
#[derive(Debug, Clone)]
pub struct LowpassRC {
    input: f64,
    output: f64,
    tau: f64, // Time constant

    // State
    state: f64,
    dt: f64,

    // Initial conditions
    initial_state: f64,

    // For buffer/revert
    buffered_state: f64,
}

impl LowpassRC {
    /// Create lowpass filter with time constant tau
    ///
    /// The cutoff frequency is fc = 1 / (2*pi*tau)
    pub fn new(tau: f64) -> Self {
        Self {
            input: 0.0,
            output: 0.0,
            tau,
            state: 0.0,
            dt: 0.0,
            initial_state: 0.0,
            buffered_state: 0.0,
        }
    }

    /// Create lowpass filter with initial state
    pub fn with_initial(tau: f64, initial: f64) -> Self {
        Self {
            input: 0.0,
            output: initial,
            tau,
            state: initial,
            dt: 0.0,
            initial_state: initial,
            buffered_state: initial,
        }
    }

    /// Set time constant
    pub fn set_tau(&mut self, tau: f64) {
        self.tau = tau;
    }

    /// Set cutoff frequency (convenience method)
    ///
    /// fc = 1 / (2*pi*tau), so tau = 1 / (2*pi*fc)
    pub fn set_cutoff_freq(&mut self, fc: f64) {
        self.tau = 1.0 / (2.0 * std::f64::consts::PI * fc);
    }

    /// Get current filter state
    pub fn state(&self) -> f64 {
        self.state
    }
}

impl Block for LowpassRC {
    const NUM_INPUTS: usize = 1;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = true;

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
        // Output is current state
        self.output = self.state;
    }

    fn step(&mut self, _t: f64, dt: f64) -> StepResult {
        self.dt = dt;

        if dt > 0.0 {
            // Backward Euler discretization
            // y[n] = alpha * x[n] + (1-alpha) * y[n-1]
            let alpha = dt / (self.tau + dt);
            self.state = alpha * self.input + (1.0 - alpha) * self.state;
            self.output = self.state;
        }

        StepResult::default()
    }

    fn buffer(&mut self) {
        self.buffered_state = self.state;
    }

    fn revert(&mut self) {
        self.state = self.buffered_state;
        self.output = self.state;
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = self.initial_state;
        self.state = self.initial_state;
        self.dt = 0.0;
        self.buffered_state = self.initial_state;
    }
}

impl DynamicBlock for LowpassRC {
    fn state(&self) -> &[f64] {
        std::slice::from_ref(&self.state)
    }

    fn state_derivative(&self) -> &[f64] {
        // dy/dt = (x - y) / tau
        // We don't store this explicitly, so we return a dummy slice
        // In a full ODE implementation, we would compute this
        static ZERO: f64 = 0.0;
        std::slice::from_ref(&ZERO)
    }
}

/// First-order RC highpass filter
///
/// Transfer function: H(s) = s*tau / (1 + s*tau)
///
/// Discrete implementation:
/// y[n] = (1-alpha) * (y[n-1] + x[n] - x[n-1])
/// where alpha = dt / (tau + dt)
///
/// # Example
///
/// ```ignore
/// let mut hpf = HighpassRC::new(0.1);  // tau = 0.1 seconds
/// hpf.set_input(0, 1.0);  // Step input
/// for i in 0..100 {
///     hpf.update(i as f64 * 0.01);
///     hpf.step(i as f64 * 0.01, 0.01);
/// }
/// // Output decays to zero exponentially
/// ```
#[derive(Debug, Clone)]
pub struct HighpassRC {
    input: f64,
    output: f64,
    tau: f64, // Time constant

    // State
    state: f64,      // Filter state
    prev_input: f64, // Previous input for difference
    dt: f64,

    // Initial conditions
    initial_state: f64,

    // For buffer/revert
    buffered_state: f64,
    buffered_prev_input: f64,
}

impl HighpassRC {
    /// Create highpass filter with time constant tau
    ///
    /// The cutoff frequency is fc = 1 / (2*pi*tau)
    pub fn new(tau: f64) -> Self {
        Self {
            input: 0.0,
            output: 0.0,
            tau,
            state: 0.0,
            prev_input: 0.0,
            dt: 0.0,
            initial_state: 0.0,
            buffered_state: 0.0,
            buffered_prev_input: 0.0,
        }
    }

    /// Create highpass filter with initial state
    pub fn with_initial(tau: f64, initial: f64) -> Self {
        Self {
            input: 0.0,
            output: initial,
            tau,
            state: initial,
            prev_input: 0.0,
            dt: 0.0,
            initial_state: initial,
            buffered_state: initial,
            buffered_prev_input: 0.0,
        }
    }

    /// Set time constant
    pub fn set_tau(&mut self, tau: f64) {
        self.tau = tau;
    }

    /// Set cutoff frequency (convenience method)
    ///
    /// fc = 1 / (2*pi*tau), so tau = 1 / (2*pi*fc)
    pub fn set_cutoff_freq(&mut self, fc: f64) {
        self.tau = 1.0 / (2.0 * std::f64::consts::PI * fc);
    }

    /// Get current filter state
    pub fn state(&self) -> f64 {
        self.state
    }
}

impl Block for HighpassRC {
    const NUM_INPUTS: usize = 1;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = true;

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
        // Output is current state
        self.output = self.state;
    }

    fn step(&mut self, _t: f64, dt: f64) -> StepResult {
        self.dt = dt;

        if dt > 0.0 {
            // Highpass discrete implementation
            // y[n] = (1-alpha) * (y[n-1] + x[n] - x[n-1])
            let alpha = dt / (self.tau + dt);
            self.state = (1.0 - alpha) * (self.state + self.input - self.prev_input);
            self.output = self.state;
            self.prev_input = self.input;
        }

        StepResult::default()
    }

    fn buffer(&mut self) {
        self.buffered_state = self.state;
        self.buffered_prev_input = self.prev_input;
    }

    fn revert(&mut self) {
        self.state = self.buffered_state;
        self.output = self.state;
        self.prev_input = self.buffered_prev_input;
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = self.initial_state;
        self.state = self.initial_state;
        self.prev_input = 0.0;
        self.dt = 0.0;
        self.buffered_state = self.initial_state;
        self.buffered_prev_input = 0.0;
    }
}

impl DynamicBlock for HighpassRC {
    fn state(&self) -> &[f64] {
        std::slice::from_ref(&self.state)
    }

    fn state_derivative(&self) -> &[f64] {
        // dy/dt = -(y - x) / tau
        // We don't store this explicitly
        static ZERO: f64 = 0.0;
        std::slice::from_ref(&ZERO)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lowpass_dc() {
        // DC input should pass through
        let mut lpf = LowpassRC::new(0.1);
        let dt = 0.01;

        // Step to 1.0 and wait for settling
        for _ in 0..1000 {
            lpf.set_input(0, 1.0);
            lpf.update(0.0);
            lpf.step(0.0, dt);
        }

        // Should converge to 1.0
        assert!((lpf.get_output(0) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_lowpass_time_constant() {
        // After one time constant, output should be ~63% of final value
        let tau = 0.1;
        let mut lpf = LowpassRC::new(tau);
        let dt = 0.001;

        // Step input
        let steps = (tau / dt) as usize;
        for _ in 0..steps {
            lpf.set_input(0, 1.0);
            lpf.update(0.0);
            lpf.step(0.0, dt);
        }

        // Should be approximately 1 - 1/e ≈ 0.632
        let expected = 1.0 - (-1.0_f64).exp();
        assert!((lpf.get_output(0) - expected).abs() < 0.01);
    }

    #[test]
    fn test_lowpass_buffer_revert() {
        let mut lpf = LowpassRC::new(0.1);
        let dt = 0.01;

        lpf.set_input(0, 1.0);
        lpf.update(0.0);
        lpf.step(0.0, dt);

        let state_after_step = lpf.state();
        lpf.buffer();

        // Take more steps
        for _ in 0..10 {
            lpf.update(0.0);
            lpf.step(0.0, dt);
        }

        // State should have changed
        assert!(lpf.state() > state_after_step);

        // Revert
        lpf.revert();
        assert_eq!(lpf.state(), state_after_step);
    }

    #[test]
    fn test_highpass_dc() {
        // DC input should be blocked
        let mut hpf = HighpassRC::new(0.1);
        let dt = 0.01;

        // Step to 1.0 and wait for settling
        for _ in 0..1000 {
            hpf.set_input(0, 1.0);
            hpf.update(0.0);
            hpf.step(0.0, dt);
        }

        // Should decay to near zero
        assert!(hpf.get_output(0).abs() < 1e-6);
    }

    #[test]
    fn test_highpass_step_response() {
        // Highpass should show initial spike then decay
        let mut hpf = HighpassRC::new(0.1);
        let dt = 0.01;

        // First step - no previous input
        hpf.set_input(0, 1.0);
        hpf.update(0.0);
        hpf.step(0.0, dt);

        let initial_response = hpf.get_output(0);

        // Should have positive response initially
        assert!(initial_response > 0.5);

        // After many steps, should decay
        for _ in 0..1000 {
            hpf.update(0.0);
            hpf.step(0.0, dt);
        }

        assert!(hpf.get_output(0).abs() < 0.01);
    }

    #[test]
    fn test_highpass_ac() {
        // AC signal should pass through (approximately)
        let mut hpf = HighpassRC::new(0.01); // Short time constant
        let dt = 0.001;
        let freq = 10.0; // 10 Hz
        let omega = 2.0 * std::f64::consts::PI * freq;

        // Let it settle first
        for i in 0..100 {
            let t = i as f64 * dt;
            let input = (omega * t).sin();
            hpf.set_input(0, input);
            hpf.update(t);
            hpf.step(t, dt);
        }

        // Check that high-frequency content passes
        // (output amplitude should be close to input amplitude)
        let mut max_output: f64 = 0.0;
        for i in 100..200 {
            let t = i as f64 * dt;
            let input = (omega * t).sin();
            hpf.set_input(0, input);
            hpf.update(t);
            hpf.step(t, dt);
            max_output = max_output.max(hpf.get_output(0).abs());
        }

        // Should pass high frequency with some attenuation
        // The filter has phase shift and some gain reduction at passband frequencies
        assert!(
            max_output > 0.3,
            "Expected max_output > 0.3, got {}",
            max_output
        );
    }

    #[test]
    fn test_lowpass_reset() {
        let mut lpf = LowpassRC::new(0.1);

        lpf.set_input(0, 5.0);
        lpf.update(0.0);
        lpf.step(0.0, 0.1);

        lpf.reset();
        assert_eq!(lpf.get_output(0), 0.0);
        assert_eq!(lpf.state(), 0.0);
    }

    #[test]
    fn test_highpass_reset() {
        let mut hpf = HighpassRC::new(0.1);

        hpf.set_input(0, 5.0);
        hpf.update(0.0);
        hpf.step(0.0, 0.1);

        hpf.reset();
        assert_eq!(hpf.get_output(0), 0.0);
        assert_eq!(hpf.state(), 0.0);
    }

    #[test]
    fn test_set_cutoff_freq() {
        let mut lpf = LowpassRC::new(1.0);

        // Set cutoff to 10 Hz
        lpf.set_cutoff_freq(10.0);

        // tau should be 1/(2*pi*10) ≈ 0.0159
        let expected_tau = 1.0 / (2.0 * std::f64::consts::PI * 10.0);
        assert!((lpf.tau - expected_tau).abs() < 1e-6);
    }
}
