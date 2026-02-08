//! Time delay block with circular buffer

use crate::block::{Block, DynamicBlock, StepResult};

/// Time delay block: y(t) = u(t - delay)
///
/// Uses a circular buffer to store past values and linear interpolation
/// for arbitrary delay times.
///
/// # Type Parameters
///
/// - `BUFFER_SIZE`: Maximum number of samples to store (const generic)
///
/// # Example
///
/// ```ignore
/// // 100-sample buffer, 0.5 second delay
/// let mut delay = Delay::<100>::new(0.5);
/// delay.set_input(0, 1.0);
/// delay.update(0.0);
/// delay.step(0.0, 0.01);
/// // Output will be 0.0 until 0.5 seconds have elapsed
/// ```
#[derive(Debug, Clone)]
pub struct Delay<const BUFFER_SIZE: usize> {
    input: f64,
    output: f64,
    /// Delay time in seconds
    delay_time: f64,
    /// Circular buffer storing (time, value) pairs
    buffer: [(f64, f64); BUFFER_SIZE],
    /// Current write position in buffer
    write_index: usize,
    /// Number of valid samples in buffer
    count: usize,
    /// Current time
    current_time: f64,

    // For buffer/revert
    buffered_write_index: usize,
    buffered_count: usize,
    buffered_buffer: [(f64, f64); BUFFER_SIZE],
}

impl<const BUFFER_SIZE: usize> Delay<BUFFER_SIZE> {
    /// Create delay block with specified delay time
    ///
    /// # Arguments
    ///
    /// * `delay_time` - Delay in seconds
    ///
    /// # Panics
    ///
    /// Panics if delay_time is negative
    pub fn new(delay_time: f64) -> Self {
        assert!(delay_time >= 0.0, "Delay time must be non-negative");
        assert!(
            BUFFER_SIZE > 1,
            "Buffer size must be at least 2 for interpolation"
        );

        Self {
            input: 0.0,
            output: 0.0,
            delay_time,
            buffer: [(0.0, 0.0); BUFFER_SIZE],
            write_index: 0,
            count: 0,
            current_time: 0.0,
            buffered_write_index: 0,
            buffered_count: 0,
            buffered_buffer: [(0.0, 0.0); BUFFER_SIZE],
        }
    }

    /// Set delay time (in seconds)
    pub fn set_delay(&mut self, delay_time: f64) {
        assert!(delay_time >= 0.0, "Delay time must be non-negative");
        self.delay_time = delay_time;
    }

    /// Get current delay time
    pub fn delay(&self) -> f64 {
        self.delay_time
    }

    /// Get delayed value using linear interpolation
    fn get_delayed_value(&self, target_time: f64) -> f64 {
        if self.count == 0 {
            return 0.0;
        }

        if self.count == 1 {
            // Only return the sample if target_time is at or after the sample time
            if target_time >= self.buffer[0].0 {
                return self.buffer[0].1;
            } else {
                return 0.0;
            }
        }

        // Find the two samples that bracket the target time
        // Samples are ordered from newest (write_index-1) to oldest
        let mut before_idx = None;
        let mut after_idx = None;

        for i in 0..self.count {
            // Calculate index going backwards from write_index
            let idx = if self.write_index > i {
                self.write_index - i - 1
            } else {
                // Wrap around: if write_index=0, i=0 -> BUFFER_SIZE-1
                BUFFER_SIZE - (i - self.write_index) - 1
            };

            let sample_time = self.buffer[idx].0;

            if sample_time <= target_time {
                before_idx = Some(idx);
                break;
            }
            after_idx = Some(idx);
        }

        match (before_idx, after_idx) {
            (Some(before), Some(after)) => {
                // Interpolate between two samples
                let (t0, y0) = self.buffer[before];
                let (t1, y1) = self.buffer[after];

                if (t1 - t0).abs() < 1e-12 {
                    return y0;
                }

                let alpha = (target_time - t0) / (t1 - t0);
                y0 + alpha * (y1 - y0)
            }
            (Some(idx), None) => {
                // Target time is before all samples, return oldest
                self.buffer[idx].1
            }
            (None, Some(idx)) => {
                // Target time is after all samples, return newest
                self.buffer[idx].1
            }
            (None, None) => 0.0,
        }
    }
}

impl<const BUFFER_SIZE: usize> Block for Delay<BUFFER_SIZE> {
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

    fn update(&mut self, t: f64) {
        // Compute delayed time
        let delayed_time = t - self.delay_time;

        // Get interpolated value at delayed time
        self.output = self.get_delayed_value(delayed_time);
    }

    fn step(&mut self, t: f64, dt: f64) -> StepResult {
        // Store current input value with timestamp
        self.buffer[self.write_index] = (self.current_time, self.input);

        // Advance write index (circular)
        self.write_index = (self.write_index + 1) % BUFFER_SIZE;

        // Update count (saturate at buffer size)
        if self.count < BUFFER_SIZE {
            self.count += 1;
        }

        // Advance time
        self.current_time = t + dt;

        StepResult::default()
    }

    fn buffer(&mut self) {
        self.buffered_write_index = self.write_index;
        self.buffered_count = self.count;
        self.buffered_buffer.copy_from_slice(&self.buffer);
    }

    fn revert(&mut self) {
        self.write_index = self.buffered_write_index;
        self.count = self.buffered_count;
        self.buffer.copy_from_slice(&self.buffered_buffer);
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = 0.0;
        self.buffer.fill((0.0, 0.0));
        self.write_index = 0;
        self.count = 0;
        self.current_time = 0.0;
        self.buffered_write_index = 0;
        self.buffered_count = 0;
        self.buffered_buffer.fill((0.0, 0.0));
    }
}

impl<const BUFFER_SIZE: usize> DynamicBlock for Delay<BUFFER_SIZE> {
    fn state(&self) -> &[f64] {
        // State is conceptually the buffer, but we return a slice of the current output
        std::slice::from_ref(&self.output)
    }

    fn state_derivative(&self) -> &[f64] {
        // Delay has no continuous derivative
        std::slice::from_ref(&0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_delay_constant_input() {
        let mut delay = Delay::<100>::new(0.5);
        let dt = 0.01;

        // Simulate: at each timestep, update() computes output, step() advances state
        // Constant input of 1.0 starting from t=0
        let mut t = 0.0;
        for _ in 0..100 {
            delay.set_input(0, 1.0);
            delay.update(t);
            delay.step(t, dt);
            t += dt;
        }

        // After 1 second of simulation with constant input 1.0 and 0.5s delay,
        // the output should be 1.0
        delay.update(t);
        assert!(
            (delay.get_output(0) - 1.0).abs() < 1e-10,
            "At t={}, output={}",
            t,
            delay.get_output(0)
        );
    }

    #[test]
    fn test_delay_step_input() {
        let mut delay = Delay::<200>::new(0.5);
        let dt = 0.01;

        // Step input at t=0.3
        let step_time = 0.3;

        for i in 0..100 {
            let t = i as f64 * dt;

            if t >= step_time {
                delay.set_input(0, 1.0);
            } else {
                delay.set_input(0, 0.0);
            }

            delay.update(t);
            delay.step(t, dt);
        }

        // Check output at t = step_time + delay_time
        let check_time = step_time + delay.delay();
        let steps = (check_time / dt) as usize;

        for i in 0..steps {
            let t = i as f64 * dt;
            if t >= step_time {
                delay.set_input(0, 1.0);
            } else {
                delay.set_input(0, 0.0);
            }
            delay.update(t);
            delay.step(t, dt);
        }

        delay.update(check_time);
        // Should be close to 1.0
        assert!(delay.get_output(0) > 0.9);
    }

    #[test]
    fn test_delay_buffer_revert() {
        let mut delay = Delay::<100>::new(0.1);
        let dt = 0.01;

        delay.set_input(0, 1.0);
        for _ in 0..5 {
            delay.update(0.0);
            delay.step(0.0, dt);
        }

        delay.buffer();

        delay.set_input(0, 2.0);
        for _ in 0..5 {
            delay.update(0.0);
            delay.step(0.0, dt);
        }

        let count_before_revert = delay.count;

        delay.revert();

        assert_ne!(delay.count, count_before_revert);
    }

    #[test]
    fn test_delay_reset() {
        let mut delay = Delay::<100>::new(0.5);
        let dt = 0.01;

        delay.set_input(0, 1.0);
        for _ in 0..100 {
            delay.update(0.0);
            delay.step(0.0, dt);
        }

        delay.reset();

        assert_eq!(delay.count, 0);
        assert_eq!(delay.write_index, 0);
        assert_eq!(delay.get_output(0), 0.0);
    }

    #[test]
    #[should_panic(expected = "Delay time must be non-negative")]
    fn test_delay_negative_time() {
        Delay::<100>::new(-1.0);
    }

    #[test]
    fn test_delay_zero_delay() {
        let mut delay = Delay::<100>::new(0.0);
        let dt = 0.01;

        delay.set_input(0, 5.0);
        delay.update(0.0);
        delay.step(0.0, dt);

        delay.update(dt);
        // With zero delay, output should closely follow input
        assert!((delay.get_output(0) - 5.0).abs() < 0.1);
    }
}
