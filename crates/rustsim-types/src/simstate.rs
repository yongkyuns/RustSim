//! Simulation state shared between interpreter and compiled code.
//!
//! This module defines the `SimState` struct which holds the current state
//! of a running simulation. It's designed to be safely passed across dylib
//! boundaries when using the same crate version and compiler.

/// Simulation state shared between interpreter and compiled code.
///
/// This structure contains all the runtime state needed to execute a simulation:
/// - Current time and time step
/// - Internal states (for integrators, delays, memory blocks, etc.)
/// - Input and output values
///
/// # Safety for dylib boundaries
///
/// This struct is designed to be passed across dylib boundaries safely when:
/// - Both sides use the same version of this crate
/// - Both sides are compiled with the same Rust compiler version
/// - The struct uses only `repr(Rust)` types with stable layouts (Vec, f64, etc.)
#[derive(Clone, Debug)]
pub struct SimState {
    /// Current simulation time in seconds
    pub time: f64,

    /// Time step (dt) in seconds
    pub dt: f64,

    /// Internal states for stateful blocks (integrators, delays, etc.)
    ///
    /// Each stateful block is allocated a contiguous range of indices
    /// in this vector to store its internal state between time steps.
    pub states: Vec<f64>,

    /// Input values from external sources
    ///
    /// These represent values coming into the simulation from outside,
    /// such as user inputs, sensor readings, or reference signals.
    pub inputs: Vec<f64>,

    /// Output values computed by the simulation
    ///
    /// These represent the final results of the simulation that should
    /// be logged, displayed, or sent to external systems.
    pub outputs: Vec<f64>,

    /// Runtime-adjustable parameters for blocks
    ///
    /// These allow block parameters (gains, frequencies, thresholds, etc.)
    /// to be modified without recompiling the simulation. Each parameter
    /// is indexed, and the mapping from parameter names to indices is
    /// maintained separately.
    pub params: Vec<f64>,
}

impl SimState {
    /// Create a new simulation state with specified capacities.
    ///
    /// # Arguments
    ///
    /// * `num_states` - Number of internal state variables to allocate
    /// * `num_inputs` - Number of input values to allocate
    /// * `num_outputs` - Number of output values to allocate
    ///
    /// # Examples
    ///
    /// ```
    /// use rustsim_types::SimState;
    ///
    /// // Create a simulation with 10 internal states, 2 inputs, and 3 outputs
    /// let sim_state = SimState::new(10, 2, 3);
    /// assert_eq!(sim_state.states.len(), 10);
    /// assert_eq!(sim_state.inputs.len(), 2);
    /// assert_eq!(sim_state.outputs.len(), 3);
    /// ```
    pub fn new(num_states: usize, num_inputs: usize, num_outputs: usize) -> Self {
        Self::with_params(num_states, num_inputs, num_outputs, 0)
    }

    /// Create a new simulation state with specified capacities including parameters.
    ///
    /// # Arguments
    ///
    /// * `num_states` - Number of internal state variables to allocate
    /// * `num_inputs` - Number of input values to allocate
    /// * `num_outputs` - Number of output values to allocate
    /// * `num_params` - Number of runtime parameters to allocate
    pub fn with_params(num_states: usize, num_inputs: usize, num_outputs: usize, num_params: usize) -> Self {
        Self {
            time: 0.0,
            dt: 0.001, // Default 1ms time step
            states: vec![0.0; num_states],
            inputs: vec![0.0; num_inputs],
            outputs: vec![0.0; num_outputs],
            params: vec![0.0; num_params],
        }
    }

    /// Reset the simulation state to initial conditions.
    ///
    /// This sets time to 0.0 and clears all states, inputs, and outputs
    /// to zero while preserving their allocated capacity. Parameters are
    /// NOT reset as they represent user-configured values.
    pub fn reset(&mut self) {
        self.time = 0.0;
        self.states.fill(0.0);
        self.inputs.fill(0.0);
        self.outputs.fill(0.0);
        // Note: params are intentionally NOT reset - they represent user settings
    }

    /// Set a parameter value by index.
    pub fn set_param(&mut self, index: usize, value: f64) {
        if index < self.params.len() {
            self.params[index] = value;
        }
    }

    /// Get a parameter value by index.
    pub fn get_param(&self, index: usize) -> f64 {
        self.params.get(index).copied().unwrap_or(0.0)
    }

    /// Advance the simulation time by one time step.
    ///
    /// This is typically called at the end of each simulation iteration.
    pub fn step_time(&mut self) {
        self.time += self.dt;
    }
}

impl Default for SimState {
    /// Create a default simulation state with no states, inputs, or outputs.
    ///
    /// The time starts at 0.0 with a default time step of 1ms (0.001s).
    fn default() -> Self {
        Self {
            time: 0.0,
            dt: 0.001,
            states: Vec::new(),
            inputs: Vec::new(),
            outputs: Vec::new(),
            params: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let state = SimState::new(5, 2, 3);
        assert_eq!(state.time, 0.0);
        assert_eq!(state.dt, 0.001);
        assert_eq!(state.states.len(), 5);
        assert_eq!(state.inputs.len(), 2);
        assert_eq!(state.outputs.len(), 3);
        assert!(state.states.iter().all(|&x| x == 0.0));
    }

    #[test]
    fn test_default() {
        let state = SimState::default();
        assert_eq!(state.time, 0.0);
        assert_eq!(state.dt, 0.001);
        assert!(state.states.is_empty());
        assert!(state.inputs.is_empty());
        assert!(state.outputs.is_empty());
    }

    #[test]
    fn test_reset() {
        let mut state = SimState::new(3, 1, 1);
        state.time = 5.0;
        state.states[0] = 1.0;
        state.inputs[0] = 2.0;
        state.outputs[0] = 3.0;

        state.reset();

        assert_eq!(state.time, 0.0);
        assert_eq!(state.states[0], 0.0);
        assert_eq!(state.inputs[0], 0.0);
        assert_eq!(state.outputs[0], 0.0);
    }

    #[test]
    fn test_step_time() {
        let mut state = SimState::new(0, 0, 0);
        assert_eq!(state.time, 0.0);

        state.step_time();
        assert!((state.time - 0.001).abs() < 1e-10);

        state.step_time();
        assert!((state.time - 0.002).abs() < 1e-10);
    }

    #[test]
    fn test_clone() {
        let mut state1 = SimState::new(2, 1, 1);
        state1.time = 1.5;
        state1.states[0] = 10.0;
        state1.inputs[0] = 20.0;

        let state2 = state1.clone();
        assert_eq!(state2.time, 1.5);
        assert_eq!(state2.states[0], 10.0);
        assert_eq!(state2.inputs[0], 20.0);
    }
}
