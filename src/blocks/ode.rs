//! General ODE block with RK4 solver

use crate::block::{Block, DynamicBlock, StepResult};

/// General ODE block: dx/dt = f(t, x, u)
///
/// Integrates a system of N ordinary differential equations using RK4.
///
/// # Type Parameters
///
/// - `N`: State dimension (const generic)
/// - `F`: Dynamics function with signature `fn(t, state, inputs, derivatives)`
///
/// # Example
///
/// ```ignore
/// // Harmonic oscillator: d²x/dt² = -k*x - c*dx/dt
/// // State: [position, velocity]
/// let ode = ODE::<2, _>::new(
///     [1.0, 0.0],  // Initial state: x=1, v=0
///     |_t, state, _inputs, derivs| {
///         let k = 1.0;  // Spring constant
///         let c = 0.1;  // Damping
///         derivs[0] = state[1];              // dx/dt = v
///         derivs[1] = -k * state[0] - c * state[1];  // dv/dt = -kx - cv
///     }
/// );
/// ```
#[derive(Clone)]
pub struct ODE<const N: usize, F>
where
    F: Fn(f64, &[f64; N], &[f64], &mut [f64; N]),
{
    /// External inputs (variable size)
    inputs: Vec<f64>,
    /// State vector outputs
    outputs: [f64; N],
    /// Current state
    state: [f64; N],
    /// State derivative (computed during update)
    derivative: [f64; N],
    /// Initial state
    initial: [f64; N],
    /// Dynamics function
    dynamics: F,
    /// Number of input ports
    #[allow(dead_code)]
    num_inputs: usize,

    // For buffer/revert
    buffered_state: [f64; N],

    // RK4 intermediate values
    k1: [f64; N],
    k2: [f64; N],
    k3: [f64; N],
    k4: [f64; N],
    temp_state: [f64; N],
}

impl<const N: usize, F> ODE<N, F>
where
    F: Fn(f64, &[f64; N], &[f64], &mut [f64; N]),
{
    /// Create ODE block with initial state and dynamics function
    ///
    /// # Arguments
    ///
    /// * `initial` - Initial state vector
    /// * `dynamics` - Function computing derivatives: f(t, state, inputs, derivatives)
    pub fn new(initial: [f64; N], dynamics: F) -> Self {
        Self::with_inputs(initial, 0, dynamics)
    }

    /// Create ODE block with specified number of inputs
    ///
    /// # Arguments
    ///
    /// * `initial` - Initial state vector
    /// * `num_inputs` - Number of input ports
    /// * `dynamics` - Function computing derivatives: f(t, state, inputs, derivatives)
    pub fn with_inputs(initial: [f64; N], num_inputs: usize, dynamics: F) -> Self {
        Self {
            inputs: vec![0.0; num_inputs],
            outputs: initial,
            state: initial,
            derivative: [0.0; N],
            initial,
            dynamics,
            num_inputs,
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

    /// Set state value at index (use with caution)
    pub fn set_state_value(&mut self, index: usize, value: f64) {
        self.state[index] = value;
        self.outputs[index] = value;
    }

    /// Get entire state vector
    pub fn state_vector(&self) -> &[f64; N] {
        &self.state
    }
}

impl<const N: usize, F> Block for ODE<N, F>
where
    F: Fn(f64, &[f64; N], &[f64], &mut [f64; N]),
{
    const NUM_INPUTS: usize = 0; // Dynamic, overridden by num_inputs field
    const NUM_OUTPUTS: usize = N;
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

    fn update(&mut self, t: f64) {
        // Output current state
        self.outputs.copy_from_slice(&self.state);

        // Compute derivative for current state
        (self.dynamics)(t, &self.state, &self.inputs, &mut self.derivative);
    }

    fn step(&mut self, t: f64, dt: f64) -> StepResult {
        // RK4 integration: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

        // k1 = f(t, y)
        (self.dynamics)(t, &self.state, &self.inputs, &mut self.k1);

        // k2 = f(t + dt/2, y + dt/2 * k1)
        for i in 0..N {
            self.temp_state[i] = self.state[i] + 0.5 * dt * self.k1[i];
        }
        (self.dynamics)(t + 0.5 * dt, &self.temp_state, &self.inputs, &mut self.k2);

        // k3 = f(t + dt/2, y + dt/2 * k2)
        for i in 0..N {
            self.temp_state[i] = self.state[i] + 0.5 * dt * self.k2[i];
        }
        (self.dynamics)(t + 0.5 * dt, &self.temp_state, &self.inputs, &mut self.k3);

        // k4 = f(t + dt, y + dt * k3)
        for i in 0..N {
            self.temp_state[i] = self.state[i] + dt * self.k3[i];
        }
        (self.dynamics)(t + dt, &self.temp_state, &self.inputs, &mut self.k4);

        // y_new = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        for i in 0..N {
            self.state[i] +=
                dt / 6.0 * (self.k1[i] + 2.0 * self.k2[i] + 2.0 * self.k3[i] + self.k4[i]);
        }

        self.outputs.copy_from_slice(&self.state);

        StepResult::default()
    }

    fn buffer(&mut self) {
        self.buffered_state.copy_from_slice(&self.state);
    }

    fn revert(&mut self) {
        self.state.copy_from_slice(&self.buffered_state);
        self.outputs.copy_from_slice(&self.state);
    }

    fn reset(&mut self) {
        self.inputs.fill(0.0);
        self.state.copy_from_slice(&self.initial);
        self.outputs.copy_from_slice(&self.initial);
        self.derivative.fill(0.0);
        self.buffered_state.copy_from_slice(&self.initial);
    }
}

impl<const N: usize, F> DynamicBlock for ODE<N, F>
where
    F: Fn(f64, &[f64; N], &[f64], &mut [f64; N]),
{
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
    fn test_ode_scalar_exponential() {
        // dy/dt = y, y(0) = 1, solution: y(t) = exp(t)
        let mut ode = ODE::<1, _>::new([1.0], |_t, state, _inputs, derivs| {
            derivs[0] = state[0];
        });

        let dt = 0.01;
        let t_final = 1.0;
        let steps = (t_final / dt) as usize;

        for i in 0..steps {
            let t = i as f64 * dt;
            ode.update(t);
            ode.step(t, dt);
        }

        // At t=1, y should be approximately e ≈ 2.71828
        let expected = std::f64::consts::E;
        let error = (ode.state_value(0) - expected).abs();
        assert!(
            error < 1e-4,
            "Error: {}, expected: {}, got: {}",
            error,
            expected,
            ode.state_value(0)
        );
    }

    #[test]
    fn test_ode_harmonic_oscillator() {
        // d²x/dt² = -x, x(0) = 1, v(0) = 0
        // Solution: x(t) = cos(t), v(t) = -sin(t)
        let mut ode = ODE::<2, _>::new([1.0, 0.0], |_t, state, _inputs, derivs| {
            derivs[0] = state[1]; // dx/dt = v
            derivs[1] = -state[0]; // dv/dt = -x
        });

        let dt = 0.01;
        let t_final = 2.0 * std::f64::consts::PI; // One full period

        let mut t = 0.0;
        while t < t_final {
            ode.update(t);
            ode.step(t, dt);
            t += dt;
        }

        // Should return to initial state (allow some error accumulation over full period)
        assert!(
            (ode.state_value(0) - 1.0).abs() < 1e-2,
            "Position error: {}",
            (ode.state_value(0) - 1.0).abs()
        );
        assert!(
            ode.state_value(1).abs() < 1e-2,
            "Velocity error: {}",
            ode.state_value(1).abs()
        );
    }

    #[test]
    fn test_ode_with_inputs() {
        // dy/dt = u (integrator), y(0) = 0
        let mut ode = ODE::<1, _>::with_inputs([0.0], 1, |_t, _state, inputs, derivs| {
            derivs[0] = inputs[0];
        });

        let dt = 0.01;

        // Constant input of 2.0
        ode.set_input(0, 2.0);
        for _ in 0..100 {
            ode.update(0.0);
            ode.step(0.0, dt);
        }

        // After 1 second, should be approximately 2.0
        assert!((ode.state_value(0) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_ode_buffer_revert() {
        let mut ode = ODE::<1, _>::new([0.0], |_t, state, _inputs, derivs| {
            derivs[0] = state[0] + 1.0;
        });

        ode.update(0.0);
        ode.step(0.0, 0.1);
        let value_after_step = ode.state_value(0);

        ode.buffer();

        ode.update(0.0);
        ode.step(0.0, 0.1);

        // Value changed
        assert!(ode.state_value(0) != value_after_step);

        // Revert
        ode.revert();
        assert_eq!(ode.state_value(0), value_after_step);
    }

    #[test]
    fn test_ode_reset() {
        let initial = [1.0, 2.0];
        let mut ode = ODE::<2, _>::new(initial, |_t, _state, _inputs, derivs| {
            derivs[0] = 1.0;
            derivs[1] = 2.0;
        });

        ode.update(0.0);
        ode.step(0.0, 1.0);

        // State changed
        assert!(ode.state_value(0) != initial[0]);
        assert!(ode.state_value(1) != initial[1]);

        ode.reset();
        assert_eq!(ode.state_value(0), initial[0]);
        assert_eq!(ode.state_value(1), initial[1]);
    }

    #[test]
    fn test_ode_lorenz_system() {
        // Lorenz system: chaotic attractor
        // dx/dt = sigma * (y - x)
        // dy/dt = x * (rho - z) - y
        // dz/dt = x * y - beta * z
        let sigma = 10.0;
        let rho = 28.0;
        let beta = 8.0 / 3.0;

        let mut ode = ODE::<3, _>::new([1.0, 1.0, 1.0], move |_t, state, _inputs, derivs| {
            let x = state[0];
            let y = state[1];
            let z = state[2];
            derivs[0] = sigma * (y - x);
            derivs[1] = x * (rho - z) - y;
            derivs[2] = x * y - beta * z;
        });

        let dt = 0.001;

        // Run for a short time
        for _ in 0..1000 {
            ode.update(0.0);
            ode.step(0.0, dt);
        }

        // Just verify it doesn't blow up (Lorenz system is bounded)
        assert!(ode.state_value(0).abs() < 100.0);
        assert!(ode.state_value(1).abs() < 100.0);
        assert!(ode.state_value(2).abs() < 100.0);
    }
}
