//! StateSpace block: Linear Time-Invariant (LTI) system
//!
//! Implements state-space representation:
//!   dx/dt = Ax + Bu
//!   y = Cx + Du
//!
//! where:
//!   - A: n×n state matrix
//!   - B: n×m input matrix
//!   - C: p×n output matrix
//!   - D: p×m feedthrough matrix
//!   - x: n-dimensional state vector
//!   - u: m-dimensional input vector
//!   - y: p-dimensional output vector

use crate::block::{Block, DynamicBlock, StepResult};

/// StateSpace: Linear Time-Invariant state-space model
///
/// # Type Parameters
///
/// - `N`: Number of states
/// - `M`: Number of inputs
/// - `P`: Number of outputs
///
/// # Example (SISO - First order system)
///
/// ```ignore
/// // RC circuit: dy/dt = -y/RC + u/RC, y = x
/// let rc = 1.0;
/// let ss = StateSpace::<1, 1, 1>::new(
///     [[-1.0/rc]],         // A
///     [[1.0/rc]],          // B
///     [[1.0]],             // C
///     [[0.0]],             // D
///     [0.0]                // initial state
/// );
/// ```
///
/// # Example (MIMO - 2 inputs, 2 outputs, 2 states)
///
/// ```ignore
/// let ss = StateSpace::<2, 2, 2>::new(
///     [[-1.0, 0.0],       // A matrix
///      [0.0, -2.0]],
///     [[1.0, 0.0],        // B matrix
///      [0.0, 1.0]],
///     [[1.0, 0.0],        // C matrix
///      [0.0, 1.0]],
///     [[0.0, 0.0],        // D matrix
///      [0.0, 0.0]],
///     [0.0, 0.0]          // initial state
/// );
/// ```
#[derive(Clone)]
pub struct StateSpace<const N: usize, const M: usize, const P: usize> {
    /// State matrix (N×N)
    a: [[f64; N]; N],
    /// Input matrix (N×M)
    b: [[f64; M]; N],
    /// Output matrix (P×N)
    c: [[f64; N]; P],
    /// Feedthrough matrix (P×M)
    d: [[f64; M]; P],

    /// External inputs (M-dimensional)
    inputs: [f64; M],
    /// System outputs (P-dimensional)
    outputs: [f64; P],
    /// Current state (N-dimensional)
    state: [f64; N],
    /// State derivative (computed during update)
    derivative: [f64; N],
    /// Initial state
    initial: [f64; N],

    /// Whether D has any non-zero elements (for passthrough detection)
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

impl<const N: usize, const M: usize, const P: usize> StateSpace<N, M, P> {
    /// Create StateSpace block with matrices and initial state
    ///
    /// # Arguments
    ///
    /// * `a` - State matrix (N×N)
    /// * `b` - Input matrix (N×M)
    /// * `c` - Output matrix (P×N)
    /// * `d` - Feedthrough matrix (P×M)
    /// * `initial` - Initial state vector (N-dimensional)
    pub fn new(
        a: [[f64; N]; N],
        b: [[f64; M]; N],
        c: [[f64; N]; P],
        d: [[f64; M]; P],
        initial: [f64; N],
    ) -> Self {
        // Check if D has any non-zero elements
        let mut has_passthrough = false;
        'outer: for row in &d {
            for &val in row {
                if val != 0.0 {
                    has_passthrough = true;
                    break 'outer;
                }
            }
        }

        // Compute initial output: y = Cx + Du (where u = 0 initially)
        let mut outputs = [0.0; P];
        for i in 0..P {
            for j in 0..N {
                outputs[i] += c[i][j] * initial[j];
            }
        }

        Self {
            a,
            b,
            c,
            d,
            inputs: [0.0; M],
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

    /// Set state value at index (use with caution)
    pub fn set_state_value(&mut self, index: usize, value: f64) {
        self.state[index] = value;
    }

    /// Get entire state vector
    pub fn state_vector(&self) -> &[f64; N] {
        &self.state
    }

    /// Compute state derivative: dx/dt = Ax + Bu
    fn compute_derivative(&self, state: &[f64; N], inputs: &[f64; M]) -> [f64; N] {
        let mut deriv = [0.0; N];

        // Ax
        for i in 0..N {
            for j in 0..N {
                deriv[i] += self.a[i][j] * state[j];
            }
        }

        // Bu
        for i in 0..N {
            for j in 0..M {
                deriv[i] += self.b[i][j] * inputs[j];
            }
        }

        deriv
    }

    /// Compute output: y = Cx + Du
    fn compute_output(&self, state: &[f64; N], inputs: &[f64; M]) -> [f64; P] {
        let mut output = [0.0; P];

        // Cx
        for i in 0..P {
            for j in 0..N {
                output[i] += self.c[i][j] * state[j];
            }
        }

        // Du
        for i in 0..P {
            for j in 0..M {
                output[i] += self.d[i][j] * inputs[j];
            }
        }

        output
    }
}

impl<const N: usize, const M: usize, const P: usize> Block for StateSpace<N, M, P> {
    const NUM_INPUTS: usize = M;
    const NUM_OUTPUTS: usize = P;
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
        // RK4 integration: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

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
            self.state[i] += dt / 6.0 * (self.k1[i] + 2.0 * self.k2[i] + 2.0 * self.k3[i] + self.k4[i]);
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

impl<const N: usize, const M: usize, const P: usize> DynamicBlock for StateSpace<N, M, P> {
    fn state(&self) -> &[f64] {
        &self.state
    }

    fn state_derivative(&self) -> &[f64] {
        &self.derivative
    }
}

/// Helper to detect if StateSpace has passthrough (D != 0)
///
/// In PathSim, this is done via __len__ returning 1 if passthrough exists.
/// We expose this as a method instead.
impl<const N: usize, const M: usize, const P: usize> StateSpace<N, M, P> {
    /// Returns true if the D matrix has any non-zero elements (direct passthrough)
    pub fn has_passthrough(&self) -> bool {
        self.has_passthrough
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_statespace_init_siso() {
        // Single-input single-output
        let ss = StateSpace::<2, 1, 1>::new(
            [[1.0, 0.0], [0.0, 1.0]],  // A = I
            [[1.0], [1.0]],             // B = ones
            [[1.0, 1.0]],               // C = ones
            [[1.0]],                    // D = 1
            [0.0, 0.0],                 // initial = zeros
        );

        assert_eq!(ss.inputs().len(), 1);
        assert_eq!(ss.outputs().len(), 1);
        assert_eq!(ss.state().len(), 2);
    }

    #[test]
    fn test_statespace_init_mimo() {
        // Multi-input multi-output
        let ss = StateSpace::<2, 2, 2>::new(
            [[1.0, 0.0], [0.0, 1.0]],           // A = I
            [[1.0, 1.0], [1.0, 1.0]],           // B = ones
            [[1.0, 1.0], [1.0, 1.0]],           // C = ones
            [[1.0, 1.0], [1.0, 1.0]],           // D = ones
            [1.0, 1.0],                         // initial = ones
        );

        assert_eq!(ss.inputs().len(), 2);
        assert_eq!(ss.outputs().len(), 2);
        assert_eq!(ss.state().len(), 2);
        assert_eq!(ss.state_value(0), 1.0);
        assert_eq!(ss.state_value(1), 1.0);
    }

    #[test]
    fn test_statespace_passthrough() {
        // SISO with no passthrough (D = 0)
        let ss_no_passthrough = StateSpace::<1, 1, 1>::new(
            [[0.0]],
            [[1.0]],
            [[1.0]],
            [[0.0]],  // D = 0
            [0.0],
        );
        assert!(!ss_no_passthrough.has_passthrough());

        // SISO with passthrough (D != 0)
        let ss_with_passthrough = StateSpace::<1, 1, 1>::new(
            [[0.0]],
            [[1.0]],
            [[1.0]],
            [[3.0]],  // D != 0
            [0.0],
        );
        assert!(ss_with_passthrough.has_passthrough());

        // MIMO with no passthrough (D = zeros)
        let ss_mimo_no = StateSpace::<2, 2, 2>::new(
            [[1.0, 0.0], [0.0, 1.0]],
            [[1.0, 0.0], [0.0, 1.0]],
            [[1.0, 0.0], [0.0, 1.0]],
            [[0.0, 0.0], [0.0, 0.0]],  // D = zeros
            [0.0, 0.0],
        );
        assert!(!ss_mimo_no.has_passthrough());

        // MIMO with passthrough (D has non-zero element)
        let ss_mimo_yes = StateSpace::<2, 2, 2>::new(
            [[1.0, 0.0], [0.0, 1.0]],
            [[1.0, 0.0], [0.0, 1.0]],
            [[1.0, 0.0], [0.0, 1.0]],
            [[0.0, 5.0], [0.0, 0.0]],  // D has non-zero element
            [0.0, 0.0],
        );
        assert!(ss_mimo_yes.has_passthrough());
    }

    #[test]
    fn test_statespace_first_order() {
        // First-order system: dy/dt = -ay + bu, y = x
        // Example: RC circuit with R=1, C=1, so tau=1
        // dy/dt = -y + u
        let a = -1.0;
        let b = 1.0;
        let c = 1.0;
        let d = 0.0;

        let mut ss = StateSpace::<1, 1, 1>::new(
            [[a]],
            [[b]],
            [[c]],
            [[d]],
            [0.0],  // y(0) = 0
        );

        // Apply step input u = 1.0
        ss.set_input(0, 1.0);

        let dt = 0.01;
        let t_final = 5.0;  // 5 time constants
        let steps = (t_final / dt) as usize;

        for _ in 0..steps {
            ss.update(0.0);
            ss.step(0.0, dt);
        }

        // Steady state should be y = u/(-a) = 1.0/1.0 = 1.0
        // After 5 time constants, should be ~99.3% of steady state
        let expected = 1.0 * (1.0 - (-5.0_f64).exp());
        assert!(
            (ss.get_output(0) - expected).abs() < 1e-2,
            "Output: {}, Expected: {}",
            ss.get_output(0),
            expected
        );
    }

    #[test]
    fn test_statespace_second_order() {
        // Second-order system: mass-spring-damper
        // m * d²x/dt² + c * dx/dt + k * x = F
        // State: [x, v] where v = dx/dt
        // dx/dt = v
        // dv/dt = (F - k*x - c*v) / m

        let m = 1.0;  // mass
        let c = 0.5;  // damping
        let k = 1.0;  // spring constant

        // State space form:
        // dx/dt = [0, 1; -k/m, -c/m] * x + [0; 1/m] * u
        // y = [1, 0] * x  (output is position)

        let mut ss = StateSpace::<2, 1, 1>::new(
            [[0.0, 1.0], [-k/m, -c/m]],  // A
            [[0.0], [1.0/m]],             // B
            [[1.0, 0.0]],                 // C (output position)
            [[0.0]],                      // D
            [0.0, 0.0],                   // initial: x=0, v=0
        );

        // Apply step input F = 1.0
        ss.set_input(0, 1.0);

        let dt = 0.01;

        // Run for a while
        for _ in 0..1000 {
            ss.update(0.0);
            ss.step(0.0, dt);
        }

        // Steady state: x = F/k = 1.0/1.0 = 1.0, v = 0
        // After enough time, position should settle to 1.0
        assert!(
            (ss.get_output(0) - 1.0).abs() < 0.1,
            "Final position: {}, expected ~1.0",
            ss.get_output(0)
        );
    }

    #[test]
    fn test_statespace_update() {
        // Test that update computes correct output
        let mut ss = StateSpace::<1, 1, 1>::new(
            [[-1.0]],  // A
            [[1.0]],   // B
            [[-1.0]],  // C
            [[1.0]],   // D
            [1.1],     // initial state
        );

        // Initially output should be C*x + D*u = -1.0*1.1 + 1.0*0.0 = -1.1
        assert_eq!(ss.get_output(0), -1.1);

        // Set input
        ss.set_input(0, 3.3);
        ss.update(0.0);

        // Output should be C*x + D*u = -1.0*1.1 + 1.0*3.3 = 2.2
        assert!((ss.get_output(0) - 2.2).abs() < 1e-10);
    }

    #[test]
    fn test_statespace_step() {
        // Test state integration
        let mut ss = StateSpace::<1, 1, 1>::new(
            [[-1.0]],  // A
            [[1.0]],   // B
            [[1.0]],   // C
            [[0.0]],   // D
            [0.0],     // initial state
        );

        ss.set_input(0, 1.0);
        ss.update(0.0);

        let dt = 0.1;
        ss.step(0.0, dt);

        // After one step with dy/dt = -y + u = -0 + 1 = 1
        // y ≈ 0 + 1*0.1 = 0.1 (Euler approximation)
        // RK4 will be slightly different but close
        assert!(ss.state_value(0) > 0.05);
        assert!(ss.state_value(0) < 0.15);
    }

    #[test]
    fn test_statespace_reset() {
        let initial = [1.0, 2.0];
        let mut ss = StateSpace::<2, 1, 2>::new(
            [[0.0, 1.0], [-1.0, 0.0]],
            [[0.0], [1.0]],
            [[1.0, 0.0], [0.0, 1.0]],
            [[0.0], [0.0]],
            initial,
        );

        ss.set_input(0, 5.0);
        ss.update(0.0);
        ss.step(0.0, 0.1);

        // State should have changed
        assert!(ss.state_value(0) != initial[0] || ss.state_value(1) != initial[1]);

        ss.reset();

        // State should be back to initial
        assert_eq!(ss.state_value(0), initial[0]);
        assert_eq!(ss.state_value(1), initial[1]);

        // Inputs should be zero
        assert_eq!(ss.get_input(0), 0.0);
    }

    #[test]
    fn test_statespace_buffer_revert() {
        let mut ss = StateSpace::<2, 1, 1>::new(
            [[0.0, 1.0], [-1.0, 0.0]],
            [[0.0], [1.0]],
            [[1.0, 0.0]],
            [[0.0]],
            [1.0, 0.0],
        );

        ss.set_input(0, 1.0);
        ss.update(0.0);
        ss.step(0.0, 0.1);

        let state_after_step = [ss.state_value(0), ss.state_value(1)];
        ss.buffer();

        // Take another step
        ss.update(0.0);
        ss.step(0.0, 0.1);

        // State should have changed
        assert!(ss.state_value(0) != state_after_step[0] || ss.state_value(1) != state_after_step[1]);

        // Revert
        ss.revert();
        assert_eq!(ss.state_value(0), state_after_step[0]);
        assert_eq!(ss.state_value(1), state_after_step[1]);
    }

    #[test]
    fn test_statespace_mimo_computation() {
        // MIMO system with 2 inputs, 2 outputs, 2 states
        // Test that matrix operations work correctly
        let mut ss = StateSpace::<2, 2, 2>::new(
            [[1.0, 0.0], [0.0, 1.0]],           // A = I
            [[1.0, 1.0], [1.0, 1.0]],           // B = ones
            [[1.0, 1.0], [1.0, 1.0]],           // C = ones
            [[1.0, 1.0], [1.0, 1.0]],           // D = ones
            [0.0, 0.0],                         // initial = zeros
        );

        // Set inputs to sin(t) and cos(t)
        let t: f64 = 1.0;
        ss.set_input(0, t.sin());
        ss.set_input(1, t.cos());
        ss.update(t);

        // With x=[0,0], u=[sin(t), cos(t)]
        // y = Cx + Du = 0 + ones*[sin(t), cos(t)] = [sin(t)+cos(t), sin(t)+cos(t)]
        let expected = t.sin() + t.cos();
        assert!((ss.get_output(0) - expected).abs() < 1e-10);
        assert!((ss.get_output(1) - expected).abs() < 1e-10);
    }
}
