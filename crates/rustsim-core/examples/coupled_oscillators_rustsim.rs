//! Coupled Oscillator Example - Connection and Simulation API
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

use rustsim_core::anyblock::AnyBlock;
use rustsim_core::blocks::ODE;
use rustsim_core::Simulation;
use std::f64::consts::PI;
use std::time::Instant;

/// Coupled oscillator system parameters
struct CoupledOscillatorParams {
    n: usize, // Number of oscillators
    m: f64,   // Mass
    k: f64,   // Spring constant
}

impl CoupledOscillatorParams {
    /// Initialize state with sinusoidal displacement pattern
    fn initial_state(&self) -> Vec<f64> {
        let mut state = vec![0.0; 2 * self.n];

        // Set initial positions (first N elements)
        for i in 0..self.n {
            state[i] = (PI * (i + 1) as f64 / (self.n + 1) as f64).sin();
        }

        // Velocities (second N elements) are zero
        state
    }

    /// Compute total energy (kinetic + potential)
    fn energy(&self, state: &[f64]) -> f64 {
        let n = self.n;

        // Kinetic energy: (1/2) * m * sum(v_i^2)
        let mut kinetic = 0.0;
        for i in 0..n {
            kinetic += 0.5 * self.m * state[n + i].powi(2);
        }

        // Potential energy: (1/2) * k * sum((x_i - x_{i-1})^2 + (x_i - x_{i+1})^2)
        let mut potential = 0.0;
        for i in 0..n {
            let x_i = state[i];

            // Left spring
            let x_left = if i == 0 { 0.0 } else { state[i - 1] };
            potential += 0.5 * self.k * (x_i - x_left).powi(2);

            // Right spring
            let x_right = if i == n - 1 { 0.0 } else { state[i + 1] };
            potential += 0.5 * self.k * (x_i - x_right).powi(2);
        }

        kinetic + potential
    }
}

/// Create coupled oscillator dynamics function
///
/// Returns a closure that computes derivatives for the ODE block
fn create_dynamics<const N: usize>(
    n: usize,
    m: f64,
    k: f64,
    c: f64,
) -> impl Fn(f64, &[f64; N], &[f64], &mut [f64; N]) {
    move |_t: f64, state: &[f64; N], _inputs: &[f64], derivs: &mut [f64; N]| {
        // State layout: [positions..., velocities...]
        // Derivatives: [velocities..., accelerations...]

        for i in 0..n {
            let x_i = state[i];
            let v_i = state[n + i];

            // dx_i/dt = v_i
            derivs[i] = v_i;

            // Compute force on oscillator i
            let mut force = 0.0;

            // Left spring force: -k*(x_i - x_{i-1})
            let x_left = if i == 0 { 0.0 } else { state[i - 1] };
            force -= k * (x_i - x_left);

            // Right spring force: -k*(x_i - x_{i+1})
            let x_right = if i == n - 1 { 0.0 } else { state[i + 1] };
            force -= k * (x_i - x_right);

            // Damping force: -c*v_i
            force -= c * v_i;

            // dv_i/dt = force / m
            derivs[n + i] = force / m;
        }
    }
}

fn main() {
    println!("{}", "=".repeat(70));
    println!("Coupled Oscillators - Connection and Simulation API");
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
    println!("  Solver: RK4 (via ODE block)");
    println!();

    // Create system parameters
    println!("Creating coupled oscillator system...");
    let params = CoupledOscillatorParams { n, m, k };

    // Get initial state
    let initial_state = params.initial_state();
    let initial_energy = params.energy(&initial_state);
    println!("  Initial energy: {:.6}", initial_energy);
    println!();

    // Create ODE block with dynamic state size
    // We need to use a fixed-size array for the ODE block, so we'll use a const generic
    // For this example, we'll hardcode N=20 which gives us 40 states
    const STATE_DIM: usize = 40;
    let mut initial_array = [0.0; STATE_DIM];
    initial_array.copy_from_slice(&initial_state);

    let dynamics = create_dynamics::<STATE_DIM>(n, m, k, c);
    let ode = ODE::new(initial_array, dynamics);

    // Create the simulation with just the ODE block
    let blocks: Vec<Box<dyn AnyBlock>> = vec![Box::new(ode)];
    let connections = vec![]; // No connections needed for single block

    let mut sim = Simulation::new(blocks, connections).with_dt(dt);

    // Run simulation
    println!("Running simulation...");
    let wall_start = Instant::now();

    let n_steps = ((t_end - t_start) / dt) as usize;
    sim.run_steps(n_steps);

    let wall_elapsed = wall_start.elapsed();
    let wall_time = wall_elapsed.as_secs_f64();

    println!("Wall-clock time: {:.4} seconds", wall_time);
    println!();

    // Compute final energy
    // Get state from ODE block outputs
    let mut final_state = vec![0.0; STATE_DIM];
    for i in 0..STATE_DIM {
        final_state[i] = sim.get_output(0, i);
    }

    let final_energy = params.energy(&final_state);
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
