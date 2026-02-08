//! Sample and hold block

use crate::block::{AlgebraicBlock, Block};

/// Sample and hold block
///
/// Samples the input periodically at rate T (sampling period) and holds the value at the output.
/// Optionally supports a delay tau before the first sample.
///
/// # Parameters
///
/// - `T`: Sampling period (time between samples)
/// - `tau`: Initial delay before first sample (default: 0.0)
///
/// # Behavior
///
/// The block samples the input at times: tau, tau + T, tau + 2T, tau + 3T, ...
/// Between samples, the output holds the last sampled value.
///
/// # Example
///
/// ```ignore
/// let mut sh = SampleHold::new(1.0, 0.0);
/// sh.set_input(0, 5.0);
/// sh.update(0.5); // Before first sample at t=1.0
/// assert_eq!(sh.get_output(0), 0.0); // Still holding initial value
///
/// sh.update(1.0); // At sample time
/// assert_eq!(sh.get_output(0), 5.0); // Sampled the input
/// ```
#[derive(Debug, Clone)]
pub struct SampleHold {
    input: f64,
    output: f64,
    period: f64,                   // T - sampling period
    delay: f64,                    // tau - initial delay
    last_sample_time: Option<f64>, // Track when we last sampled
}

impl SampleHold {
    /// Create a new sample-hold with sampling period T and delay tau
    pub fn new(period: f64, delay: f64) -> Self {
        assert!(period > 0.0, "Sampling period must be positive");
        assert!(delay >= 0.0, "Delay must be non-negative");

        Self {
            input: 0.0,
            output: 0.0,
            period,
            delay,
            last_sample_time: None,
        }
    }

    /// Create a sample-hold with default delay (0.0)
    pub fn with_period(period: f64) -> Self {
        Self::new(period, 0.0)
    }

    /// Get sampling period
    pub fn period(&self) -> f64 {
        self.period
    }

    /// Get delay
    pub fn delay(&self) -> f64 {
        self.delay
    }

    /// Get the next scheduled sample time
    pub fn next_sample_time(&self, current_time: f64) -> f64 {
        if current_time < self.delay {
            self.delay
        } else {
            let elapsed = current_time - self.delay;
            let periods = (elapsed / self.period).floor();
            self.delay + (periods + 1.0) * self.period
        }
    }

    /// Check if we should sample at this time
    fn should_sample(&self, t: f64) -> bool {
        if t < self.delay {
            return false;
        }

        match self.last_sample_time {
            None => t >= self.delay,
            Some(last_t) => {
                // Check if enough time has passed for the next sample
                // We use a small epsilon to handle floating point comparison
                let elapsed = t - last_t;
                elapsed >= self.period - 1e-10
            }
        }
    }
}

impl Default for SampleHold {
    fn default() -> Self {
        Self::new(1.0, 0.0)
    }
}

impl Block for SampleHold {
    const NUM_INPUTS: usize = 1;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        std::slice::from_ref(&self.input)
    }

    fn inputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.input)
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, t: f64) {
        // Check if we should sample at this time
        if self.should_sample(t) {
            self.output = self.input;
            self.last_sample_time = Some(t);
        }
        // Otherwise, hold the previous value (do nothing)
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = 0.0;
        self.last_sample_time = None;
    }
}

impl AlgebraicBlock for SampleHold {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_samplehold_init() {
        // Default initialization
        let sh = SampleHold::default();
        assert_eq!(sh.period(), 1.0);
        assert_eq!(sh.delay(), 0.0);

        // Custom initialization
        let sh = SampleHold::new(0.5, 0.1);
        assert_eq!(sh.period(), 0.5);
        assert_eq!(sh.delay(), 0.1);
    }

    #[test]
    fn test_samplehold_event_scheduling() {
        let sh = SampleHold::new(2.0, 0.5);

        // Next sample time should be at delay (0.5)
        assert_eq!(sh.next_sample_time(0.0), 0.5);

        // After first sample, next should be at 0.5 + 2.0 = 2.5
        assert_eq!(sh.next_sample_time(1.0), 2.5);
        assert_eq!(sh.next_sample_time(2.0), 2.5);
        assert_eq!(sh.next_sample_time(2.5), 4.5);
    }

    #[test]
    fn test_samplehold_single_sample() {
        let mut sh = SampleHold::new(1.0, 0.0);

        // Set input
        sh.set_input(0, 5.0);

        // Manually trigger sampling at t=0 (first sample time)
        sh.update(0.0);

        // Output should now match input
        assert_eq!(sh.get_output(0), 5.0);
    }

    #[test]
    fn test_samplehold_multiple_samples() {
        let mut sh = SampleHold::new(1.0, 0.0);

        // First sample
        sh.set_input(0, 3.0);
        sh.update(0.0);
        assert_eq!(sh.get_output(0), 3.0);

        // Second sample with different input
        sh.set_input(0, 7.0);
        sh.update(1.0);
        assert_eq!(sh.get_output(0), 7.0);

        // Third sample
        sh.set_input(0, -2.5);
        sh.update(2.0);
        assert_eq!(sh.get_output(0), -2.5);
    }

    #[test]
    fn test_samplehold_hold_behavior() {
        let mut sh = SampleHold::new(1.0, 0.0);

        // Initial sample at t=0
        sh.set_input(0, 10.0);
        sh.update(0.0);
        assert_eq!(sh.get_output(0), 10.0);

        // Change input but don't sample (at t=0.5, not a sample time)
        sh.set_input(0, 20.0);
        sh.update(0.5);
        assert_eq!(sh.get_output(0), 10.0); // Still holding previous value

        // Now sample at t=1.0
        sh.update(1.0);
        assert_eq!(sh.get_output(0), 20.0); // Now updated
    }

    #[test]
    fn test_samplehold_with_delay() {
        let mut sh = SampleHold::new(1.0, 0.5);

        // Before delay, should not sample
        sh.set_input(0, 5.0);
        sh.update(0.0);
        assert_eq!(sh.get_output(0), 0.0);

        sh.update(0.3);
        assert_eq!(sh.get_output(0), 0.0);

        // At delay time, should sample
        sh.update(0.5);
        assert_eq!(sh.get_output(0), 5.0);

        // Next sample at 0.5 + 1.0 = 1.5
        sh.set_input(0, 10.0);
        sh.update(1.0);
        assert_eq!(sh.get_output(0), 5.0); // Still holding

        sh.update(1.5);
        assert_eq!(sh.get_output(0), 10.0); // Sampled
    }

    #[test]
    fn test_samplehold_reset() {
        let mut sh = SampleHold::new(1.0, 0.0);

        sh.set_input(0, 5.0);
        sh.update(0.0);
        assert_eq!(sh.get_output(0), 5.0);

        sh.reset();
        assert_eq!(sh.get_input(0), 0.0);
        assert_eq!(sh.get_output(0), 0.0);
        assert!(sh.last_sample_time.is_none());
    }

    #[test]
    #[should_panic(expected = "Sampling period must be positive")]
    fn test_samplehold_invalid_period() {
        let _sh = SampleHold::new(0.0, 0.0);
    }

    #[test]
    #[should_panic(expected = "Delay must be non-negative")]
    fn test_samplehold_invalid_delay() {
        let _sh = SampleHold::new(1.0, -0.1);
    }
}
