//! Simulation settings types.

use serde::{Deserialize, Serialize};

/// Available ODE solvers
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SolverType {
    /// Forward Euler (1st order)
    Euler,
    /// Heun's method (2nd order)
    Heun,
    /// Classical RK4 (4th order)
    RK4,
    /// Runge-Kutta-Cash-Karp 5(4)
    RKCK54,
    /// Dormand-Prince 5(4)
    DOPRI54,
    /// Bogacki-Shampine 3(2)
    RKBS32,
    /// Strong Stability Preserving RK2
    SSPRK22,
    /// Strong Stability Preserving RK3
    SSPRK33,
}

impl Default for SolverType {
    fn default() -> Self {
        SolverType::RK4
    }
}

impl SolverType {
    pub fn as_str(&self) -> &'static str {
        match self {
            SolverType::Euler => "Euler",
            SolverType::Heun => "Heun",
            SolverType::RK4 => "RK4",
            SolverType::RKCK54 => "RKCK54",
            SolverType::DOPRI54 => "DOPRI54",
            SolverType::RKBS32 => "RKBS32",
            SolverType::SSPRK22 => "SSPRK22",
            SolverType::SSPRK33 => "SSPRK33",
        }
    }
}

/// Simulation settings
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationSettings {
    /// Simulation duration
    pub duration: f64,

    /// Time step
    pub dt: f64,

    /// ODE solver type
    pub solver: SolverType,

    /// Enable adaptive time stepping
    pub adaptive: bool,

    /// Absolute tolerance for adaptive stepping
    pub atol: f64,

    /// Relative tolerance for adaptive stepping
    pub rtol: f64,

    /// Fixed-point iteration tolerance
    pub ftol: f64,

    /// Minimum time step
    pub dt_min: f64,

    /// Maximum time step
    pub dt_max: f64,

    /// Number of ghost traces to show in plots
    pub ghost_traces: usize,

    /// Auto-show plot panel after simulation
    pub plot_results: bool,
}

impl Default for SimulationSettings {
    fn default() -> Self {
        Self {
            duration: 10.0,
            dt: 0.01,
            solver: SolverType::RK4,
            adaptive: false,
            atol: 1e-7,
            rtol: 1e-5,
            ftol: 1e-6,
            dt_min: 1e-8,
            dt_max: 0.1,
            ghost_traces: 3,
            plot_results: true,
        }
    }
}
