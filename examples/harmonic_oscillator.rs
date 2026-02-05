//! Harmonic oscillator example using v2 architecture
//!
//! Demonstrates the compile-time static approach without macros.
//!
//! System: d²x/dt² = -x
//! Solution: x(t) = x0*cos(t) + v0*sin(t)
//!
//! With x0=1, v0=0: x(t) = cos(t)

use rustsim::{Amplifier, Block, Integrator};

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
        // position -> gain
        self.gain.set_input(0, self.position.get_output(0));
        // gain -> velocity
        self.velocity.set_input(0, self.gain.get_output(0));
        // velocity -> position
        self.position.set_input(0, self.velocity.get_output(0));
    }

    /// Evaluate algebraic blocks
    #[inline]
    fn update(&mut self) {
        let t = self.time;

        // Update in topological order
        // 1. Position outputs its state
        self.position.update(t);
        self.propagate();

        // 2. Gain computes -position
        self.gain.update(t);
        self.propagate();

        // 3. Velocity outputs its state (acceleration is now set)
        self.velocity.update(t);
        self.propagate();
    }

    /// Advance simulation by dt
    fn step(&mut self, dt: f64) {
        let t = self.time;

        // Update algebraic relationships
        self.update();

        // Step dynamic blocks (integrators)
        self.velocity.step(t, dt);
        self.position.step(t, dt);

        // Advance time
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
    println!("Harmonic Oscillator - v2 Architecture Demo");
    println!("==========================================");
    println!();
    println!("System: d²x/dt² = -x");
    println!("Initial: x(0) = 1.0, v(0) = 0.0");
    println!("Exact solution: x(t) = cos(t)");
    println!();

    let mut sim = HarmonicOscillator::new(1.0, 0.0);
    let dt = 0.001;

    println!(
        "{:>10} {:>12} {:>12} {:>12}",
        "Time", "Position", "Exact", "Error"
    );
    println!("{:-<10} {:-<12} {:-<12} {:-<12}", "", "", "", "");

    // Simulate for 2π (one period)
    let period = 2.0 * std::f64::consts::PI;
    let steps = (period / dt) as usize;

    for i in 0..=steps {
        if i % (steps / 10) == 0 {
            let t = sim.time();
            let x = sim.position();
            let exact = t.cos();
            let error = (x - exact).abs();
            println!("{:10.4} {:12.6} {:12.6} {:12.2e}", t, x, exact, error);
        }
        sim.step(dt);
    }

    println!();
    println!("After one period (t = 2π):");
    println!("  Position: {:.6} (should be ≈ 1.0)", sim.position());
    println!("  Velocity: {:.6} (should be ≈ 0.0)", sim.velocity());
    println!();

    // Energy conservation check
    let energy = 0.5 * (sim.position().powi(2) + sim.velocity().powi(2));
    println!("  Energy:   {:.6} (should be 0.5)", energy);
}
