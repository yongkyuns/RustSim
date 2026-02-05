//! Simulation constants and defaults

/// Default simulation timestep
pub const SIM_TIMESTEP: f64 = 0.001;

/// Minimum timestep for adaptive solvers
pub const SIM_TIMESTEP_MIN: f64 = 1e-12;

/// Maximum timestep for adaptive solvers
pub const SIM_TIMESTEP_MAX: f64 = 1.0;

/// Tolerance for fixed-point iteration convergence
pub const SIM_TOLERANCE_FPI: f64 = 1e-10;

/// Maximum iterations for fixed-point loops
pub const SIM_ITERATIONS_MAX: usize = 100;

/// Small tolerance for numerical comparisons
pub const TOLERANCE: f64 = 1e-16;

/// Minimum scale factor for timestep adjustment
pub const SOL_SCALE_MIN: f64 = 0.1;

/// Maximum scale factor for timestep adjustment
pub const SOL_SCALE_MAX: f64 = 10.0;

/// Safety factor for adaptive error control
pub const SOL_BETA: f64 = 0.9;

/// Default absolute tolerance for local truncation error
pub const SOL_TOLERANCE_LTE_ABS: f64 = 1e-6;

/// Default relative tolerance for local truncation error
pub const SOL_TOLERANCE_LTE_REL: f64 = 1e-3;

/// Default tolerance for event detection
pub const EVT_TOLERANCE: f64 = 1e-4;
