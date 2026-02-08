//! Harmonic oscillator benchmark in Rust
//!
//! This implements the same harmonic oscillator simulation as the Python example
//! for performance comparison.
//!
//! System: d²x/dt² = -x
//! Solution: x(t) = x0*cos(t) + v0*sin(t)
//! With x0=1, v0=0: x(t) = cos(t)

use rustsim::{Amplifier, Block, Integrator};
use std::time::Instant;

/// Harmonic oscillator simulation
///
/// Block diagram:
/// ```text
///  ┌──────────────────────────────────┐
///  │                                  │
///  ▼                                  │
/// [gain: -1] ──► [velocity] ──► [position]
/// ```
struct HarmonicOscillator {
    gain: Amplifier,
    velocity: Integrator,
    position: Integrator,
    time: f64,
}

impl HarmonicOscillator {
    fn new(x0: f64, v0: f64) -> Self {
        Self {
            gain: Amplifier::new(-1.0),
            velocity: Integrator::new(v0),
            position: Integrator::new(x0),
            time: 0.0,
        }
    }

    /// Propagate connections: wire outputs to inputs
    #[inline]
    fn propagate(&mut self) {
        self.gain.set_input(0, self.position.get_output(0));
        self.velocity.set_input(0, self.gain.get_output(0));
        self.position.set_input(0, self.velocity.get_output(0));
    }

    /// Evaluate algebraic blocks
    #[inline]
    fn update(&mut self) {
        let t = self.time;

        self.position.update(t);
        self.propagate();

        self.gain.update(t);
        self.propagate();

        self.velocity.update(t);
        self.propagate();
    }

    /// Advance simulation by dt
    #[inline]
    fn step(&mut self, dt: f64) {
        let t = self.time;

        self.update();
        self.velocity.step(t, dt);
        self.position.step(t, dt);

        self.time += dt;
    }

    /// Run for duration
    fn run(&mut self, duration: f64, dt: f64) {
        let end = self.time + duration;
        while self.time < end {
            self.step(dt);
        }
    }

    fn position(&self) -> f64 {
        self.position.get_output(0)
    }

    fn velocity(&self) -> f64 {
        self.velocity.get_output(0)
    }

    fn time(&self) -> f64 {
        self.time
    }
}

fn main() {
    println!("Harmonic Oscillator - Rust Benchmark");
    println!("{}", "=".repeat(50));
    println!();

    // Simulation parameters
    let x0 = 1.0;
    let v0 = 0.0;
    let dt = 0.001;
    let duration = 2.0 * std::f64::consts::PI; // One period

    println!("System: d²x/dt² = -x");
    println!("Initial: x(0) = {}, v(0) = {}", x0, v0);
    println!("Duration: {:.6} seconds", duration);
    println!("Time step: {}", dt);
    println!("Total steps: {}", (duration / dt) as usize);
    println!();

    // Create simulator
    let mut sim = HarmonicOscillator::new(x0, v0);

    // Run benchmark
    let start = Instant::now();
    sim.run(duration, dt);
    let elapsed = start.elapsed();

    let elapsed_secs = elapsed.as_secs_f64();

    // Calculate final values and errors
    let final_position = sim.position();
    let final_velocity = sim.velocity();
    let exact_position = sim.time().cos();
    let exact_velocity = -sim.time().sin();

    let position_error = (final_position - exact_position).abs();
    let velocity_error = (final_velocity - exact_velocity).abs();

    // Energy (should be conserved at 0.5)
    let energy = 0.5 * (final_position.powi(2) + final_velocity.powi(2));
    let energy_error = (energy - 0.5).abs();

    // Print results
    println!("Results:");
    println!("{}", "-".repeat(50));
    println!("  Execution time: {:.6} seconds", elapsed_secs);
    println!(
        "  Steps per second: {:.0}",
        (duration / dt) / elapsed_secs
    );
    println!();
    println!("Final state (t = {:.6}):", sim.time());
    println!(
        "  Position: {:.6} (exact: {:.6})",
        final_position, exact_position
    );
    println!(
        "  Velocity: {:.6} (exact: {:.6})",
        final_velocity, exact_velocity
    );
    println!("  Position error: {:.2e}", position_error);
    println!("  Velocity error: {:.2e}", velocity_error);
    println!("  Energy: {:.6} (should be 0.5)", energy);
    println!("  Energy error: {:.2e}", energy_error);
    println!();

    // Output JSON for benchmark script
    if std::env::args().any(|arg| arg == "--json") {
        let json_output = format!(
            r#"{{
  "language": "Rust",
  "elapsed_time": {},
  "steps": {},
  "steps_per_second": {},
  "final_state": {{
    "time": {},
    "position": {},
    "velocity": {}
  }},
  "exact_state": {{
    "position": {},
    "velocity": {}
  }},
  "errors": {{
    "position": {},
    "velocity": {},
    "energy": {}
  }},
  "energy": {}
}}"#,
            elapsed_secs,
            (duration / dt) as usize,
            (duration / dt) / elapsed_secs,
            sim.time(),
            final_position,
            final_velocity,
            exact_position,
            exact_velocity,
            position_error,
            velocity_error,
            energy_error,
            energy
        );

        std::fs::write("/tmp/rust_benchmark_results.json", json_output)
            .expect("Failed to write results");

        println!("Results written to /tmp/rust_benchmark_results.json");
    }
}
