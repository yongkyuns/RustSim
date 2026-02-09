//! Core Block trait for v2 architecture
//!
//! All I/O sizes are fixed at compile time via associated constants.

/// Result of a simulation step
#[derive(Debug, Clone, Copy, Default)]
pub struct StepResult {
    /// Error norm for adaptive stepping
    pub error_norm: f64,
    /// Suggested timestep scale factor
    pub scale: Option<f64>,
}

/// Core block trait - all sizes known at compile time
///
/// # Design
///
/// - `NUM_INPUTS` and `NUM_OUTPUTS` are compile-time constants
/// - I/O stored as fixed-size arrays, not Vec
/// - No dynamic dispatch - blocks are concrete types
///
/// # Example
///
/// ```ignore
/// pub struct Amplifier {
///     input: f64,
///     output: f64,
///     gain: f64,
/// }
///
/// impl Block for Amplifier {
///     const NUM_INPUTS: usize = 1;
///     const NUM_OUTPUTS: usize = 1;
///     // ...
/// }
/// ```
pub trait Block {
    /// Number of input ports (compile-time constant)
    ///
    /// Note: Also available as a method via `num_inputs()` for trait object compatibility.
    const NUM_INPUTS: usize;

    /// Number of output ports (compile-time constant)
    ///
    /// Note: Also available as a method via `num_outputs()` for trait object compatibility.
    const NUM_OUTPUTS: usize;

    /// Returns true if this block has dynamic state (integrators, ODEs)
    ///
    /// Note: Also available as a method via `is_dynamic()` for trait object compatibility.
    const IS_DYNAMIC: bool = false;

    /// Get number of input ports (runtime method for trait objects)
    fn num_inputs(&self) -> usize {
        Self::NUM_INPUTS
    }

    /// Get number of output ports (runtime method for trait objects)
    fn num_outputs(&self) -> usize {
        Self::NUM_OUTPUTS
    }

    /// Access inputs as slice
    fn inputs(&self) -> &[f64];

    /// Mutable access to inputs
    fn inputs_mut(&mut self) -> &mut [f64];

    /// Access outputs as slice
    fn outputs(&self) -> &[f64];

    /// Mutable access to outputs (rarely needed)
    fn outputs_mut(&mut self) -> &mut [f64];

    /// Evaluate algebraic relationship: outputs = f(inputs, t)
    ///
    /// Called during the update phase of each timestep.
    fn update(&mut self, t: f64);

    /// Advance dynamic state by dt
    ///
    /// Only called for blocks with `IS_DYNAMIC = true`.
    /// Default implementation does nothing.
    fn step(&mut self, _t: f64, _dt: f64) -> StepResult {
        StepResult::default()
    }

    /// Buffer current state for potential revert (adaptive stepping)
    fn buffer(&mut self) {}

    /// Revert to buffered state
    fn revert(&mut self) {}

    /// Reset to initial conditions
    fn reset(&mut self);

    /// Get single input value (convenience)
    #[inline]
    fn get_input(&self, port: usize) -> f64 {
        self.inputs()[port]
    }

    /// Set single input value (convenience)
    #[inline]
    fn set_input(&mut self, port: usize, value: f64) {
        self.inputs_mut()[port] = value;
    }

    /// Get single output value (convenience)
    #[inline]
    fn get_output(&self, port: usize) -> f64 {
        self.outputs()[port]
    }

    /// Check if this block is dynamic (has state that needs integration)
    ///
    /// Default implementation returns IS_DYNAMIC constant.
    /// This method enables runtime checking for trait objects.
    #[inline]
    fn is_dynamic(&self) -> bool {
        Self::IS_DYNAMIC
    }
}

/// Marker trait for algebraic blocks (no internal state)
pub trait AlgebraicBlock: Block {}

/// Marker trait for dynamic blocks (have internal state, need solver)
pub trait DynamicBlock: Block {
    /// Get current state vector
    fn state(&self) -> &[f64];

    /// Get state derivative (computed during update)
    fn state_derivative(&self) -> &[f64];
}
