//! Simulation plugin trait and types.

use thiserror::Error;

/// Errors that can occur during plugin operations
#[derive(Error, Debug)]
pub enum PluginError {
    #[error("Compilation failed: {0}")]
    CompilationFailed(String),

    #[error("Loading failed: {0}")]
    LoadingFailed(String),

    #[error("Simulation step failed: {0}")]
    StepFailed(String),

    #[error("Plugin not loaded")]
    NotLoaded,

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
}

/// Result type for plugin operations
pub type PluginResult<T> = Result<T, PluginError>;

/// Unified interface for simulation plugins (works for both native and WASM)
pub trait SimulationPlugin: Send {
    /// Perform one simulation step
    fn step(&mut self, dt: f64) -> PluginResult<()>;

    /// Get a specific output value
    fn get_output(&self, index: usize) -> f64;

    /// Get all output values
    fn get_outputs(&self) -> Vec<f64> {
        (0..self.output_count()).map(|i| self.get_output(i)).collect()
    }

    /// Set a specific input value
    fn set_input(&mut self, index: usize, value: f64);

    /// Set all input values
    fn set_inputs(&mut self, values: &[f64]) {
        for (i, &value) in values.iter().enumerate() {
            self.set_input(i, value);
        }
    }

    /// Reset the simulation to initial state
    fn reset(&mut self);

    /// Get current simulation time
    fn time(&self) -> f64;

    /// Get number of outputs
    fn output_count(&self) -> usize;

    /// Get number of inputs
    fn input_count(&self) -> usize;
}

/// Simulation result data from a step
#[derive(Debug, Clone, Default)]
pub struct StepResult {
    pub time: f64,
    pub outputs: Vec<f64>,
}
