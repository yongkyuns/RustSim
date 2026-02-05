//! PID controller blocks

use crate::block::{Block, DynamicBlock, StepResult};

/// Basic PID controller
///
/// Implements a PID controller with proportional, integral, and derivative terms.
///
/// # Control Law
///
/// u(t) = Kp * e(t) + Ki * ∫e(τ)dτ + Kd * de(t)/dt
///
/// where e(t) is the error signal (input).
///
/// # Example
///
/// ```ignore
/// let mut pid = PID::new(1.0, 0.5, 0.1);  // Kp=1.0, Ki=0.5, Kd=0.1
/// pid.set_input(0, error);
/// pid.update(t);
/// pid.step(t, dt);
/// let control = pid.get_output(0);
/// ```
#[derive(Debug, Clone)]
pub struct PID {
    // I/O
    input: f64,  // error signal
    output: f64, // control signal

    // Parameters
    kp: f64,
    ki: f64,
    kd: f64,

    // State
    integral: f64,
    prev_error: f64,
    dt: f64, // Stored from step()

    // Initial conditions
    initial_integral: f64,

    // For buffer/revert
    buffered_integral: f64,
    buffered_prev_error: f64,
}

impl PID {
    /// Create PID controller with gains
    pub fn new(kp: f64, ki: f64, kd: f64) -> Self {
        Self {
            input: 0.0,
            output: 0.0,
            kp,
            ki,
            kd,
            integral: 0.0,
            prev_error: 0.0,
            dt: 0.0,
            initial_integral: 0.0,
            buffered_integral: 0.0,
            buffered_prev_error: 0.0,
        }
    }

    /// Create PID with initial integral value
    pub fn with_initial_integral(kp: f64, ki: f64, kd: f64, initial_integral: f64) -> Self {
        Self {
            input: 0.0,
            output: 0.0,
            kp,
            ki,
            kd,
            integral: initial_integral,
            prev_error: 0.0,
            dt: 0.0,
            initial_integral,
            buffered_integral: initial_integral,
            buffered_prev_error: 0.0,
        }
    }

    /// Get current integral state
    pub fn integral(&self) -> f64 {
        self.integral
    }

    /// Set PID gains
    pub fn set_gains(&mut self, kp: f64, ki: f64, kd: f64) {
        self.kp = kp;
        self.ki = ki;
        self.kd = kd;
    }

    /// Reset integral to zero
    pub fn reset_integral(&mut self) {
        self.integral = 0.0;
        self.buffered_integral = 0.0;
    }
}

impl Block for PID {
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
        // Compute derivative term (backward difference)
        let derivative = if self.dt > 0.0 {
            (self.input - self.prev_error) / self.dt
        } else {
            0.0
        };

        // Compute PID output
        self.output = self.kp * self.input + self.ki * self.integral + self.kd * derivative;
    }

    fn step(&mut self, _t: f64, dt: f64) -> StepResult {
        self.dt = dt;

        // Integrate using trapezoidal rule
        self.integral += 0.5 * dt * (self.input + self.prev_error);

        // Update previous error for next derivative calculation
        self.prev_error = self.input;

        StepResult::default()
    }

    fn buffer(&mut self) {
        self.buffered_integral = self.integral;
        self.buffered_prev_error = self.prev_error;
    }

    fn revert(&mut self) {
        self.integral = self.buffered_integral;
        self.prev_error = self.buffered_prev_error;
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = 0.0;
        self.integral = self.initial_integral;
        self.prev_error = 0.0;
        self.dt = 0.0;
        self.buffered_integral = self.initial_integral;
        self.buffered_prev_error = 0.0;
    }
}

impl DynamicBlock for PID {
    fn state(&self) -> &[f64] {
        std::slice::from_ref(&self.integral)
    }

    fn state_derivative(&self) -> &[f64] {
        // Derivative of integral is the error
        std::slice::from_ref(&self.input)
    }
}

/// PID controller with anti-windup
///
/// Prevents integral windup by stopping integration when output is saturated.
///
/// # Anti-Windup Strategy
///
/// When the control output would exceed limits, the integral state is not updated.
/// This prevents the integrator from "winding up" during sustained errors.
///
/// # Example
///
/// ```ignore
/// let mut pid = AntiWindupPID::new(1.0, 0.5, 0.1, -10.0, 10.0);
/// pid.set_input(0, error);
/// pid.update(t);
/// pid.step(t, dt);
/// let control = pid.get_output(0);  // Guaranteed in [-10, 10]
/// ```
#[derive(Debug, Clone)]
pub struct AntiWindupPID {
    pid: PID,
    output_min: f64,
    output_max: f64,

    // Track if we were saturated (for anti-windup)
    is_saturated: bool,
}

impl AntiWindupPID {
    /// Create anti-windup PID with output limits
    pub fn new(kp: f64, ki: f64, kd: f64, output_min: f64, output_max: f64) -> Self {
        Self {
            pid: PID::new(kp, ki, kd),
            output_min,
            output_max,
            is_saturated: false,
        }
    }

    /// Get current integral state
    pub fn integral(&self) -> f64 {
        self.pid.integral()
    }

    /// Set PID gains
    pub fn set_gains(&mut self, kp: f64, ki: f64, kd: f64) {
        self.pid.set_gains(kp, ki, kd);
    }

    /// Set output limits
    pub fn set_limits(&mut self, output_min: f64, output_max: f64) {
        self.output_min = output_min;
        self.output_max = output_max;
    }

    /// Reset integral to zero
    pub fn reset_integral(&mut self) {
        self.pid.reset_integral();
    }
}

impl Block for AntiWindupPID {
    const NUM_INPUTS: usize = 1;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = true;

    #[inline]
    fn inputs(&self) -> &[f64] {
        self.pid.inputs()
    }

    #[inline]
    fn inputs_mut(&mut self) -> &mut [f64] {
        self.pid.inputs_mut()
    }

    #[inline]
    fn outputs(&self) -> &[f64] {
        self.pid.outputs()
    }

    #[inline]
    fn outputs_mut(&mut self) -> &mut [f64] {
        self.pid.outputs_mut()
    }

    fn update(&mut self, t: f64) {
        // Compute unsaturated PID output
        self.pid.update(t);

        // Apply saturation
        let unsaturated = self.pid.output;
        self.pid.output = unsaturated.clamp(self.output_min, self.output_max);

        // Track saturation state for anti-windup
        self.is_saturated = (unsaturated < self.output_min) || (unsaturated > self.output_max);
    }

    fn step(&mut self, t: f64, dt: f64) -> StepResult {
        // Store current integral
        let integral_before = self.pid.integral;

        // Normal PID step
        self.pid.step(t, dt);

        // Anti-windup: if saturated, revert integral update
        if self.is_saturated {
            self.pid.integral = integral_before;
        }

        StepResult::default()
    }

    fn buffer(&mut self) {
        self.pid.buffer();
    }

    fn revert(&mut self) {
        self.pid.revert();
    }

    fn reset(&mut self) {
        self.pid.reset();
        self.is_saturated = false;
    }
}

impl DynamicBlock for AntiWindupPID {
    fn state(&self) -> &[f64] {
        self.pid.state()
    }

    fn state_derivative(&self) -> &[f64] {
        self.pid.state_derivative()
    }
}

/// Rate limiter - limits the rate of change of a signal
///
/// Constrains the slew rate of the output to a maximum value.
///
/// # Example
///
/// ```ignore
/// let mut limiter = RateLimiter::new(1.0);  // Max rate of 1.0 units/sec
/// limiter.set_input(0, 10.0);  // Large step
/// limiter.update(t);
/// limiter.step(t, 0.1);
/// // Output changes by at most 0.1 (= 1.0 * 0.1)
/// ```
#[derive(Debug, Clone)]
pub struct RateLimiter {
    input: f64,
    output: f64,
    max_rate: f64,

    // State
    prev_output: f64,
    dt: f64,

    // Initial conditions
    initial_output: f64,

    // For buffer/revert
    buffered_output: f64,
}

impl RateLimiter {
    /// Create rate limiter with maximum rate
    pub fn new(max_rate: f64) -> Self {
        Self {
            input: 0.0,
            output: 0.0,
            max_rate,
            prev_output: 0.0,
            dt: 0.0,
            initial_output: 0.0,
            buffered_output: 0.0,
        }
    }

    /// Create rate limiter with initial output value
    pub fn with_initial(max_rate: f64, initial: f64) -> Self {
        Self {
            input: 0.0,
            output: initial,
            max_rate,
            prev_output: initial,
            dt: 0.0,
            initial_output: initial,
            buffered_output: initial,
        }
    }

    /// Set maximum rate
    pub fn set_max_rate(&mut self, max_rate: f64) {
        self.max_rate = max_rate;
    }
}

impl Block for RateLimiter {
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
        // Output is current state (computed in step)
        // This is called after step() in the simulation loop
    }

    fn step(&mut self, _t: f64, dt: f64) -> StepResult {
        self.dt = dt;

        if dt > 0.0 {
            // Compute desired change
            let desired_change = self.input - self.prev_output;

            // Limit the change
            let max_change = self.max_rate * dt;
            let limited_change = desired_change.clamp(-max_change, max_change);

            // Update output
            self.output = self.prev_output + limited_change;
            self.prev_output = self.output;
        } else {
            // No time elapsed, output stays the same
            self.output = self.prev_output;
        }

        StepResult::default()
    }

    fn buffer(&mut self) {
        self.buffered_output = self.output;
    }

    fn revert(&mut self) {
        self.output = self.buffered_output;
        self.prev_output = self.buffered_output;
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = self.initial_output;
        self.prev_output = self.initial_output;
        self.dt = 0.0;
        self.buffered_output = self.initial_output;
    }
}

impl DynamicBlock for RateLimiter {
    fn state(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn state_derivative(&self) -> &[f64] {
        // Rate is the derivative
        static ZERO: f64 = 0.0;
        if self.dt > 0.0 {
            // This is a simplification - rate limiter doesn't have a continuous derivative
            std::slice::from_ref(&ZERO)
        } else {
            std::slice::from_ref(&ZERO)
        }
    }
}

/// Saturation - limits output to min/max bounds
///
/// Algebraic block (no state) that clamps input to output limits.
///
/// # Example
///
/// ```ignore
/// let mut sat = Saturation::new(-5.0, 5.0);
/// sat.set_input(0, 10.0);
/// sat.update(t);
/// assert_eq!(sat.get_output(0), 5.0);
/// ```
#[derive(Debug, Clone)]
pub struct Saturation {
    input: f64,
    output: f64,
    min: f64,
    max: f64,
}

impl Saturation {
    /// Create saturation block with limits
    pub fn new(min: f64, max: f64) -> Self {
        Self {
            input: 0.0,
            output: 0.0,
            min,
            max,
        }
    }

    /// Set limits
    pub fn set_limits(&mut self, min: f64, max: f64) {
        self.min = min;
        self.max = max;
    }
}

impl Block for Saturation {
    const NUM_INPUTS: usize = 1;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

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
        self.output = self.input.clamp(self.min, self.max);
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = 0.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pid_proportional_only() {
        let mut pid = PID::new(2.0, 0.0, 0.0);
        let dt = 0.01;

        pid.set_input(0, 1.0);
        pid.update(0.0);

        // Output should be Kp * error = 2.0 * 1.0
        assert_eq!(pid.get_output(0), 2.0);
    }

    #[test]
    fn test_pid_integral() {
        let mut pid = PID::new(0.0, 1.0, 0.0); // Ki = 1.0
        let dt = 0.01;

        // Constant error of 1.0
        for _ in 0..100 {
            pid.set_input(0, 1.0);
            pid.update(0.0);
            pid.step(0.0, dt);
        }

        // Update one more time to compute output with final integral value
        pid.update(0.0);

        // After 1 second, integral should be approximately 1.0
        // Output = Ki * integral = 1.0 * 1.0 = 1.0
        // Note: trapezoidal rule gives 0.995 (first step adds 0.5*dt instead of dt)
        assert!((pid.get_output(0) - 0.995).abs() < 1e-3);
        assert!((pid.integral() - 0.995).abs() < 1e-3);
    }

    #[test]
    fn test_pid_derivative() {
        let mut pid = PID::new(0.0, 0.0, 1.0); // Kd = 1.0
        let dt = 0.01;

        // Step input
        pid.set_input(0, 0.0);
        pid.update(0.0);
        pid.step(0.0, dt);

        pid.set_input(0, 1.0);
        pid.update(0.0);

        // Derivative should be (1.0 - 0.0) / 0.01 = 100.0
        // Output = Kd * derivative = 1.0 * 100.0 = 100.0
        assert!((pid.get_output(0) - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_pid_buffer_revert() {
        let mut pid = PID::new(0.0, 1.0, 0.0);
        let dt = 0.1;

        pid.set_input(0, 1.0);
        pid.update(0.0);
        pid.step(0.0, dt);

        let integral_after_step = pid.integral();
        pid.buffer();

        // Take another step
        pid.set_input(0, 1.0);
        pid.update(0.0);
        pid.step(0.0, dt);

        // Integral should have increased
        assert!(pid.integral() > integral_after_step);

        // Revert
        pid.revert();
        assert_eq!(pid.integral(), integral_after_step);
    }

    #[test]
    fn test_anti_windup_pid_saturation() {
        let mut pid = AntiWindupPID::new(1.0, 0.5, 0.0, -10.0, 10.0);
        let dt = 0.01;

        // Large constant error - would normally windup
        for _ in 0..200 {
            pid.set_input(0, 100.0);
            pid.update(0.0);
            pid.step(0.0, dt);
        }

        // Output should be saturated
        assert_eq!(pid.get_output(0), 10.0);

        // Integral should not have wound up excessively
        // Without anti-windup, integral would be very large
        // With anti-windup, it should be bounded
        assert!(pid.integral() < 100.0);
    }

    #[test]
    fn test_rate_limiter() {
        let mut limiter = RateLimiter::new(1.0); // Max rate 1.0/sec
        let dt = 0.1;

        // Step input from 0 to 10
        limiter.set_input(0, 10.0);
        limiter.update(0.0);
        limiter.step(0.0, dt);

        // Output should change by at most max_rate * dt = 1.0 * 0.1 = 0.1
        assert!((limiter.get_output(0) - 0.1).abs() < 1e-10);

        // Continue stepping
        for _ in 0..99 {
            limiter.update(0.0);
            limiter.step(0.0, dt);
        }

        // After 10 seconds total, output should reach 10.0
        assert!((limiter.get_output(0) - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_rate_limiter_negative() {
        let mut limiter = RateLimiter::with_initial(2.0, 10.0);
        let dt = 0.1;

        // Step down to 0
        limiter.set_input(0, 0.0);
        limiter.update(0.0);
        limiter.step(0.0, dt);

        // Should decrease by max_rate * dt = 2.0 * 0.1 = 0.2
        assert!((limiter.get_output(0) - 9.8).abs() < 1e-10);
    }

    #[test]
    fn test_saturation() {
        let mut sat = Saturation::new(-5.0, 5.0);

        // Test upper saturation
        sat.set_input(0, 10.0);
        sat.update(0.0);
        assert_eq!(sat.get_output(0), 5.0);

        // Test lower saturation
        sat.set_input(0, -10.0);
        sat.update(0.0);
        assert_eq!(sat.get_output(0), -5.0);

        // Test no saturation
        sat.set_input(0, 3.0);
        sat.update(0.0);
        assert_eq!(sat.get_output(0), 3.0);
    }

    #[test]
    fn test_saturation_reset() {
        let mut sat = Saturation::new(-5.0, 5.0);
        sat.set_input(0, 10.0);
        sat.update(0.0);

        sat.reset();
        assert_eq!(sat.get_output(0), 0.0);
    }
}
