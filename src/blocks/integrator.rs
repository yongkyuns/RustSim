//! Integrator block: dy/dt = u

use crate::block::{Block, DynamicBlock, StepResult};

/// Integrator: dy/dt = u
///
/// A single-state integrator using RK4 internally.
///
/// # Example
///
/// ```ignore
/// let mut int = Integrator::new(0.0);  // Initial value 0
/// int.set_input(0, 1.0);               // Constant derivative
/// for _ in 0..100 {
///     int.update(0.0);
///     int.step(0.0, 0.01);
/// }
/// // Output should be approximately 1.0
/// ```
#[derive(Debug, Clone)]
pub struct Integrator {
    input: f64,
    output: f64,
    state: f64,
    derivative: f64,
    initial: f64,
    // For buffer/revert
    buffered_state: f64,
    // RK4 intermediate values
    k1: f64,
    k2: f64,
    k3: f64,
    k4: f64,
}

impl Integrator {
    /// Create integrator with initial value
    pub fn new(initial: f64) -> Self {
        Self {
            input: 0.0,
            output: initial,
            state: initial,
            derivative: 0.0,
            initial,
            buffered_state: initial,
            k1: 0.0,
            k2: 0.0,
            k3: 0.0,
            k4: 0.0,
        }
    }

    /// Get current state value
    pub fn value(&self) -> f64 {
        self.state
    }

    /// Set state value directly (use with caution)
    pub fn set_value(&mut self, value: f64) {
        self.state = value;
        self.output = value;
    }
}

impl Block for Integrator {
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

    #[inline]
    fn update(&mut self, _t: f64) {
        // Output is current state
        self.output = self.state;
        // Store derivative for step()
        self.derivative = self.input;
    }

    fn step(&mut self, _t: f64, dt: f64) -> StepResult {
        // RK4 integration
        // For a simple integrator dy/dt = u, if u is constant over dt:
        // k1 = k2 = k3 = k4 = u
        // y_new = y + dt * u
        //
        // But we use proper RK4 structure for consistency with ODE blocks
        self.k1 = self.derivative;
        self.k2 = self.derivative; // Assuming constant input over dt
        self.k3 = self.derivative;
        self.k4 = self.derivative;

        self.state += dt * (self.k1 + 2.0 * self.k2 + 2.0 * self.k3 + self.k4) / 6.0;
        self.output = self.state;

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
        self.output = self.initial;
        self.state = self.initial;
        self.derivative = 0.0;
        self.buffered_state = self.initial;
    }
}

impl DynamicBlock for Integrator {
    fn state(&self) -> &[f64] {
        std::slice::from_ref(&self.state)
    }

    fn state_derivative(&self) -> &[f64] {
        std::slice::from_ref(&self.derivative)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_integrator_constant_input() {
        let mut int = Integrator::new(0.0);
        let dt = 0.01;

        // Integrate constant input of 1.0 for 1 second
        for _ in 0..100 {
            int.set_input(0, 1.0);
            int.update(0.0);
            int.step(0.0, dt);
        }

        // Should be approximately 1.0
        assert!((int.value() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_integrator_initial_value() {
        let int = Integrator::new(5.0);
        assert_eq!(int.value(), 5.0);
        assert_eq!(int.get_output(0), 5.0);
    }

    #[test]
    fn test_integrator_buffer_revert() {
        let mut int = Integrator::new(0.0);
        int.set_input(0, 1.0);
        int.update(0.0);
        int.step(0.0, 0.1);

        let value_after_step = int.value();
        int.buffer();

        int.set_input(0, 10.0);
        int.update(0.0);
        int.step(0.0, 0.1);

        // Value changed
        assert!(int.value() != value_after_step);

        // Revert
        int.revert();
        assert_eq!(int.value(), value_after_step);
    }

    #[test]
    fn test_integrator_reset() {
        let mut int = Integrator::new(3.0);
        int.set_input(0, 1.0);
        int.update(0.0);
        int.step(0.0, 1.0);

        int.reset();
        assert_eq!(int.value(), 3.0);
    }
}
