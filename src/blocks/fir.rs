//! Finite Impulse Response (FIR) filter block

use crate::block::{Block, DynamicBlock, StepResult};

/// Finite Impulse Response (FIR) filter with configurable coefficients
///
/// Implements a discrete-time FIR filter that samples the input signal at
/// regular intervals and computes a weighted sum of current and past samples.
///
/// The filter equation is:
/// ```text
/// y[n] = b[0]*x[n] + b[1]*x[n-1] + ... + b[N]*x[n-N]
/// ```
/// where `b` are the filter coefficients and `N` is the filter order.
///
/// # Type Parameters
///
/// - `N`: Number of filter taps (coefficients), determines filter order
///
/// # Example
///
/// ```ignore
/// // Moving average filter: [1/3, 1/3, 1/3]
/// let mut fir = FIR::<3>::new([1.0/3.0, 1.0/3.0, 1.0/3.0], 0.01, 0.0);
/// fir.set_input(0, 1.0);
/// fir.update(0.0);
/// fir.step(0.0, 0.01);
/// ```
#[derive(Debug, Clone)]
pub struct FIR<const N: usize> {
    input: f64,
    output: f64,

    /// Filter coefficients [b0, b1, ..., b_{N-1}]
    coeffs: [f64; N],

    /// Sampling period T
    sample_period: f64,

    /// Initial delay tau before first sample
    initial_delay: f64,

    /// Circular buffer storing input samples
    buffer: [f64; N],

    /// Write index in circular buffer
    write_index: usize,

    /// Current simulation time
    current_time: f64,

    /// Time of next scheduled sample
    next_sample_time: f64,

    /// Flag indicating if we've taken any samples yet
    has_sampled: bool,

    // For buffer/revert
    buffered_buffer: [f64; N],
    buffered_write_index: usize,
    buffered_next_sample_time: f64,
    buffered_has_sampled: bool,
    buffered_output: f64,
}

impl<const N: usize> FIR<N> {
    /// Create a new FIR filter
    ///
    /// # Arguments
    ///
    /// * `coeffs` - Filter coefficients [b0, b1, ..., b_{N-1}]
    /// * `sample_period` - Time between samples (T)
    /// * `initial_delay` - Delay before first sample (tau)
    ///
    /// # Panics
    ///
    /// Panics if N is 0 or if sample_period is negative
    pub fn new(coeffs: [f64; N], sample_period: f64, initial_delay: f64) -> Self {
        assert!(N > 0, "FIR filter must have at least one coefficient");
        assert!(sample_period > 0.0, "Sample period must be positive");
        assert!(initial_delay >= 0.0, "Initial delay must be non-negative");

        Self {
            input: 0.0,
            output: 0.0,
            coeffs,
            sample_period,
            initial_delay,
            buffer: [0.0; N],
            write_index: 0,
            current_time: 0.0,
            next_sample_time: initial_delay,
            has_sampled: false,
            buffered_buffer: [0.0; N],
            buffered_write_index: 0,
            buffered_next_sample_time: initial_delay,
            buffered_has_sampled: false,
            buffered_output: 0.0,
        }
    }

    /// Get filter coefficients
    pub fn coeffs(&self) -> &[f64; N] {
        &self.coeffs
    }

    /// Get sampling period
    pub fn sample_period(&self) -> f64 {
        self.sample_period
    }

    /// Get initial delay
    pub fn initial_delay(&self) -> f64 {
        self.initial_delay
    }

    /// Get current buffer contents (for testing)
    pub fn buffer_state(&self) -> &[f64; N] {
        &self.buffer
    }

    /// Compute FIR output from current buffer state
    fn compute_output(&self) -> f64 {
        let mut sum = 0.0;

        // In PathSim, buffer is a deque with appendleft, so:
        // - buffer[0] is the newest sample (x[n])
        // - buffer[1] is x[n-1], etc.
        //
        // In our implementation:
        // - We write to buffer[write_index], then increment write_index
        // - So buffer[write_index - 1] is the newest sample we just wrote
        // - We need to read backwards from there, wrapping around
        //
        // When compute_output is called, we've already written the current sample
        // and incremented write_index, so:
        // - buffer[(write_index - 1 + N) % N] is newest (x[n])
        // - buffer[(write_index - 2 + N) % N] is x[n-1]
        // - etc.

        for i in 0..N {
            // Buffer index for x[n-i]
            let buffer_idx = (self.write_index + N - 1 - i) % N;
            sum += self.coeffs[i] * self.buffer[buffer_idx];
        }

        sum
    }
}

impl<const N: usize> Block for FIR<N> {
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
        // Output holds constant between samples
        // Already computed in step()
    }

    fn step(&mut self, t: f64, dt: f64) -> StepResult {
        let new_time = t + dt;

        // Check if we should sample at this step
        // We sample if we've crossed the next_sample_time
        if new_time >= self.next_sample_time {
            // Store current input in buffer (at write_index position)
            self.buffer[self.write_index] = self.input;

            // Advance write index (circular)
            self.write_index = (self.write_index + 1) % N;

            // Compute new output
            self.output = self.compute_output();

            // Schedule next sample
            self.next_sample_time += self.sample_period;
            self.has_sampled = true;
        }

        self.current_time = new_time;
        StepResult::default()
    }

    fn buffer(&mut self) {
        self.buffered_buffer.copy_from_slice(&self.buffer);
        self.buffered_write_index = self.write_index;
        self.buffered_next_sample_time = self.next_sample_time;
        self.buffered_has_sampled = self.has_sampled;
        self.buffered_output = self.output;
    }

    fn revert(&mut self) {
        self.buffer.copy_from_slice(&self.buffered_buffer);
        self.write_index = self.buffered_write_index;
        self.next_sample_time = self.buffered_next_sample_time;
        self.has_sampled = self.buffered_has_sampled;
        self.output = self.buffered_output;
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = 0.0;
        self.buffer.fill(0.0);
        self.write_index = 0;
        self.current_time = 0.0;
        self.next_sample_time = self.initial_delay;
        self.has_sampled = false;
        self.buffered_buffer.fill(0.0);
        self.buffered_write_index = 0;
        self.buffered_next_sample_time = self.initial_delay;
        self.buffered_has_sampled = false;
        self.buffered_output = 0.0;
    }
}

impl<const N: usize> DynamicBlock for FIR<N> {
    fn state(&self) -> &[f64] {
        // State is the output value
        std::slice::from_ref(&self.output)
    }

    fn state_derivative(&self) -> &[f64] {
        // FIR has no continuous derivative
        static ZERO: f64 = 0.0;
        std::slice::from_ref(&ZERO)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper to run filter update and step
    fn run_filter_step<const N: usize>(fir: &mut FIR<N>, input: f64, t: f64, dt: f64) {
        fir.set_input(0, input);
        fir.update(t);
        fir.step(t, dt);
    }

    #[test]
    fn test_fir_init() {
        // Default initialization with single coefficient
        let fir = FIR::<1>::new([1.0], 1.0, 0.0);
        assert_eq!(fir.coeffs(), &[1.0]);
        assert_eq!(fir.sample_period(), 1.0);
        assert_eq!(fir.initial_delay(), 0.0);

        // Custom initialization
        let coeffs = [0.5, 0.3, 0.2];
        let fir = FIR::<3>::new(coeffs, 0.1, 0.05);
        assert_eq!(fir.coeffs(), &coeffs);
        assert_eq!(fir.sample_period(), 0.1);
        assert_eq!(fir.initial_delay(), 0.05);
    }

    #[test]
    fn test_fir_len() {
        // FIR filter has no direct passthrough
        let _fir = FIR::<1>::new([1.0], 1.0, 0.0);
        // In RustSim we don't have a len() method on blocks
        // but we verify it's a dynamic block
        assert!(FIR::<1>::IS_DYNAMIC);
    }

    #[test]
    fn test_fir_event_scheduling() {
        // Verify that sampling happens at correct times
        let fir = FIR::<1>::new([1.0], 2.0, 0.5);

        assert_eq!(fir.initial_delay(), 0.5);
        assert_eq!(fir.sample_period(), 2.0);
        assert_eq!(fir.next_sample_time, 0.5);
    }

    #[test]
    fn test_fir_buffer_initialization() {
        // Buffer should be initialized with zeros
        let fir = FIR::<3>::new([0.5, 0.3, 0.2], 1.0, 0.0);

        assert_eq!(fir.buffer.len(), 3);
        for &val in fir.buffer_state() {
            assert_eq!(val, 0.0);
        }
    }

    #[test]
    fn test_fir_passthrough() {
        // Single coefficient [1.0] should pass input through
        let mut fir = FIR::<1>::new([1.0], 1.0, 0.0);

        // Set input and trigger update at t=0
        fir.set_input(0, 5.0);
        fir.update(0.0);

        // Step to t=1.0 (past first sample time at t=0)
        fir.step(0.0, 1.0);

        // Output should equal input
        assert_eq!(fir.get_output(0), 5.0);
    }

    #[test]
    fn test_fir_gain() {
        // Single coefficient for gain
        let mut fir = FIR::<1>::new([2.0], 1.0, 0.0);

        fir.set_input(0, 3.0);
        fir.update(0.0);
        fir.step(0.0, 1.0);

        // Output should be input * gain
        assert_eq!(fir.get_output(0), 6.0);
    }

    #[test]
    fn test_fir_moving_average() {
        // Moving average of 3 samples
        let mut fir = FIR::<3>::new([1.0/3.0, 1.0/3.0, 1.0/3.0], 1.0, 0.0);

        // First sample at t=0
        run_filter_step(&mut fir, 3.0, 0.0, 1.0);
        // Output: 3*1/3 + 0*1/3 + 0*1/3 = 1.0
        assert!((fir.get_output(0) - 1.0).abs() < 1e-10);

        // Second sample at t=1
        run_filter_step(&mut fir, 6.0, 1.0, 1.0);
        // Output: 6*1/3 + 3*1/3 + 0*1/3 = 3.0
        assert!((fir.get_output(0) - 3.0).abs() < 1e-10);

        // Third sample at t=2
        run_filter_step(&mut fir, 9.0, 2.0, 1.0);
        // Output: 9*1/3 + 6*1/3 + 3*1/3 = 6.0
        assert!((fir.get_output(0) - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_fir_multi_tap() {
        // Simple FIR: y[n] = x[n] + 0.5*x[n-1]
        let mut fir = FIR::<2>::new([1.0, 0.5], 1.0, 0.0);

        // First sample
        run_filter_step(&mut fir, 2.0, 0.0, 1.0);
        // Output: 2.0*1.0 + 0.0*0.5 = 2.0
        assert_eq!(fir.get_output(0), 2.0);

        // Second sample
        run_filter_step(&mut fir, 4.0, 1.0, 1.0);
        // Output: 4.0*1.0 + 2.0*0.5 = 5.0
        assert_eq!(fir.get_output(0), 5.0);

        // Third sample
        run_filter_step(&mut fir, 6.0, 2.0, 1.0);
        // Output: 6.0*1.0 + 4.0*0.5 = 8.0
        assert_eq!(fir.get_output(0), 8.0);
    }

    #[test]
    fn test_fir_difference() {
        // Difference filter: y[n] = x[n] - x[n-1]
        let mut fir = FIR::<2>::new([1.0, -1.0], 1.0, 0.0);

        // First sample
        run_filter_step(&mut fir, 5.0, 0.0, 1.0);
        // Output: 5.0*1.0 + 0.0*(-1.0) = 5.0
        assert_eq!(fir.get_output(0), 5.0);

        // Second sample
        run_filter_step(&mut fir, 8.0, 1.0, 1.0);
        // Output: 8.0*1.0 + 5.0*(-1.0) = 3.0
        assert_eq!(fir.get_output(0), 3.0);

        // Third sample
        run_filter_step(&mut fir, 10.0, 2.0, 1.0);
        // Output: 10.0*1.0 + 8.0*(-1.0) = 2.0
        assert_eq!(fir.get_output(0), 2.0);
    }

    #[test]
    fn test_fir_reset() {
        // Test that reset clears the buffer
        let mut fir = FIR::<2>::new([1.0, 0.5], 1.0, 0.0);

        // Add some data
        run_filter_step(&mut fir, 10.0, 0.0, 1.0);
        assert_eq!(fir.get_output(0), 10.0);

        // Reset
        fir.reset();

        // Buffer should be cleared
        for &val in fir.buffer_state() {
            assert_eq!(val, 0.0);
        }
        assert_eq!(fir.get_output(0), 0.0);
    }

    #[test]
    fn test_fir_buffer_revert() {
        let mut fir = FIR::<2>::new([1.0, 0.5], 1.0, 0.0);

        // Take a step
        run_filter_step(&mut fir, 5.0, 0.0, 1.0);

        // Buffer the state
        fir.buffer();
        let buffered_output = fir.get_output(0);

        // Take another step
        run_filter_step(&mut fir, 10.0, 1.0, 1.0);

        // Output should have changed
        assert_ne!(fir.get_output(0), buffered_output);

        // Revert
        fir.revert();

        // Should be back to buffered state
        assert_eq!(fir.get_output(0), buffered_output);
    }
}
