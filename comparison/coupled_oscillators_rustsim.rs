//! Coupled Oscillator Benchmark - RustSim Implementation
//!
//! Simulates N=20 coupled harmonic oscillators (chain of masses connected by springs).
//!
//! System equations:
//!     m*x_i'' = -k*(x_i - x_{i-1}) - k*(x_i - x_{i+1}) - c*x_i'
//!
//! Where:
//!     - x_i is the position of mass i
//!     - m is the mass
//!     - k is the spring constant
//!     - c is the damping coefficient
//!     - Fixed boundary conditions: x_0 = x_{N+1} = 0
//!
//! State variables: 40 total (20 positions + 20 velocities)

use nalgebra::DVector;
use rustsim::solvers::{ExplicitSolver, RK4, Solver};
use std::f64::consts::PI;
use std::time::Instant;

/// Coupled oscillator system parameters
struct CoupledOscillators {
    n: usize,      // Number of oscillators
    m: f64,        // Mass
    k: f64,        // Spring constant
    c: f64,        // Damping coefficient
}

impl CoupledOscillators {
    fn new(n: usize, m: f64, k: f64, c: f64) -> Self {
        Self { n, m, k, c }
    }

    /// Initialize state with sinusoidal displacement pattern
    fn initial_state(&self) -> DVector<f64> {
        let mut state = DVector::zeros(2 * self.n);

        // Set initial positions (first N elements)
        for i in 0..self.n {
            state[i] = (PI * (i + 1) as f64 / (self.n + 1) as f64).sin();
        }

        // Velocities (second N elements) are zero
        // Already initialized to zero

        state
    }

    /// Compute derivatives: [dx/dt, dv/dt]
    ///
    /// State vector layout:
    ///     state[0..N] = positions x_i
    ///     state[N..2N] = velocities v_i
    ///
    /// Derivatives:
    ///     dx_i/dt = v_i
    ///     dv_i/dt = (1/m) * Force_i
    ///
    /// Force_i = -k*(x_i - x_{i-1}) - k*(x_i - x_{i+1}) - c*v_i
    ///
    /// With boundary conditions:
    ///     x_{-1} = 0 (left boundary)
    ///     x_N = 0 (right boundary)
    fn derivatives(&self, state: &DVector<f64>, _t: f64) -> DVector<f64> {
        let n = self.n;
        let mut dstate = DVector::zeros(2 * n);

        // Extract positions and velocities
        let positions = state.rows(0, n);
        let velocities = state.rows(n, n);

        // Compute derivatives
        for i in 0..n {
            let x_i = positions[i];
            let v_i = velocities[i];

            // dx_i/dt = v_i
            dstate[i] = v_i;

            // Compute force on oscillator i
            let mut force = 0.0;

            // Left spring force: -k*(x_i - x_{i-1})
            let x_left = if i == 0 { 0.0 } else { positions[i - 1] };
            force -= self.k * (x_i - x_left);

            // Right spring force: -k*(x_i - x_{i+1})
            let x_right = if i == n - 1 { 0.0 } else { positions[i + 1] };
            force -= self.k * (x_i - x_right);

            // Damping force: -c*v_i
            force -= self.c * v_i;

            // dv_i/dt = force / m
            dstate[n + i] = force / self.m;
        }

        dstate
    }

    /// Compute total energy (kinetic + potential)
    fn energy(&self, state: &DVector<f64>) -> f64 {
        let n = self.n;
        let positions = state.rows(0, n);
        let velocities = state.rows(n, n);

        // Kinetic energy: (1/2) * m * sum(v_i^2)
        let mut kinetic = 0.0;
        for i in 0..n {
            kinetic += 0.5 * self.m * velocities[i].powi(2);
        }

        // Potential energy: (1/2) * k * sum((x_i - x_{i-1})^2 + (x_i - x_{i+1})^2)
        // With boundary conditions x_{-1} = 0, x_N = 0
        let mut potential = 0.0;
        for i in 0..n {
            let x_i = positions[i];

            // Left spring
            let x_left = if i == 0 { 0.0 } else { positions[i - 1] };
            potential += 0.5 * self.k * (x_i - x_left).powi(2);

            // Right spring
            let x_right = if i == n - 1 { 0.0 } else { positions[i + 1] };
            potential += 0.5 * self.k * (x_i - x_right).powi(2);
        }

        kinetic + potential
    }
}

fn main() {
    println!("{}", "=".repeat(70));
    println!("Coupled Oscillators Benchmark - RustSim");
    println!("{}", "=".repeat(70));
    println!();

    // Parameters
    let n = 20;
    let m = 1.0;
    let k = 100.0;
    let c = 0.1;

    println!("System parameters:");
    println!("  Number of oscillators: N = {}", n);
    println!("  Mass: m = {}", m);
    println!("  Spring constant: k = {}", k);
    println!("  Damping coefficient: c = {}", c);
    println!("  State dimension: {}", 2 * n);
    println!();

    // Simulation parameters
    let t_start = 0.0;
    let t_end = 100.0;
    let dt = 0.001;

    println!("Simulation parameters:");
    println!("  Time span: [{}, {}] seconds", t_start, t_end);
    println!("  Timestep: dt = {} seconds", dt);
    println!("  Total steps: {}", ((t_end - t_start) / dt) as usize);
    println!("  Solver: RK4");
    println!();

    // Create system
    println!("Creating coupled oscillator system...");
    let system = CoupledOscillators::new(n, m, k, c);

    // Initialize solver
    let initial_state = system.initial_state();
    let mut solver = RK4::new(initial_state.clone());

    // Compute initial energy
    let initial_energy = system.energy(&initial_state);
    println!("  Initial energy: {:.6}", initial_energy);
    println!();

    // Run simulation
    println!("Running simulation...");
    let wall_start = Instant::now();

    let n_steps = ((t_end - t_start) / dt) as usize;

    for _ in 0..n_steps {
        // Buffer state before RK4 stages
        solver.buffer(dt);

        // Perform 4 RK4 stages
        for _ in 0..4 {
            solver.step(|state, time| system.derivatives(state, time), dt);
        }
    }

    let wall_elapsed = wall_start.elapsed();
    let wall_time = wall_elapsed.as_secs_f64();

    println!("Wall-clock time: {:.4} seconds", wall_time);
    println!();

    // Compute final energy
    let final_state = solver.state();
    let final_energy = system.energy(final_state);
    let energy_change = (final_energy - initial_energy).abs();
    let energy_relative_change = if initial_energy > 0.0 {
        energy_change / initial_energy
    } else {
        0.0
    };

    println!("Results:");
    println!("  Final energy: {:.6}", final_energy);
    println!("  Energy change: {:.6e}", energy_change);
    println!("  Relative energy change: {:.6e}", energy_relative_change);
    println!();

    // Performance metrics
    let steps_per_second = n_steps as f64 / wall_time;
    let time_per_step = wall_time / n_steps as f64 * 1e6; // microseconds

    println!("Performance:");
    println!("  Steps per second: {:.1}", steps_per_second);
    println!("  Time per step: {:.3} microseconds", time_per_step);
    println!();

    println!("{}", "=".repeat(70));
}
