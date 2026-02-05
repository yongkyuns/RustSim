//! Base solver traits and types

use nalgebra::DVector;
use thiserror::Error;

/// Solver-related errors
#[derive(Error, Debug)]
pub enum SolverError {
    #[error("Solver did not converge after {0} iterations")]
    ConvergenceFailure(usize),

    #[error("Timestep {dt} smaller than minimum {dt_min}")]
    TimestepTooSmall { dt: f64, dt_min: f64 },

    #[error("History buffer is empty")]
    EmptyHistory,

    #[error("Invalid stage index: {stage} >= {max_stages}")]
    InvalidStage { stage: usize, max_stages: usize },
}

/// Result of a solver step
#[derive(Debug, Clone, Copy)]
pub struct SolverStepResult {
    pub success: bool,
    pub error_norm: f64,
    pub scale: Option<f64>,
}

impl Default for SolverStepResult {
    fn default() -> Self {
        Self {
            success: true,
            error_norm: 0.0,
            scale: None,
        }
    }
}

/// Core solver trait for numerical integration
pub trait Solver: Send + Sync {
    /// Get current state vector
    fn state(&self) -> &DVector<f64>;

    /// Get mutable reference to state vector
    fn state_mut(&mut self) -> &mut DVector<f64>;

    /// Set state vector
    fn set_state(&mut self, state: DVector<f64>);

    /// Buffer current state for potential reversion
    fn buffer(&mut self, dt: f64);

    /// Revert to buffered state
    fn revert(&mut self) -> Result<(), SolverError>;

    /// Reset solver to initial state
    fn reset(&mut self);

    /// Order of the method
    fn order(&self) -> usize;

    /// Number of stages
    fn stages(&self) -> usize;

    /// Is this an adaptive solver?
    fn is_adaptive(&self) -> bool;

    /// Is this an explicit solver?
    fn is_explicit(&self) -> bool;
}

/// Explicit solver trait
pub trait ExplicitSolver: Solver {
    /// Perform one step with the given right-hand side function
    fn step<F>(&mut self, f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>;
}

/// Implicit solver trait
pub trait ImplicitSolver: Solver {
    /// Solve the implicit update equation
    fn solve<F, J>(&mut self, f: F, jac: Option<J>, dt: f64) -> Result<f64, SolverError>
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>,
        J: FnMut(&DVector<f64>, f64) -> nalgebra::DMatrix<f64>;

    /// Finalize the timestep and compute error estimate
    fn step<F>(&mut self, f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>;
}
