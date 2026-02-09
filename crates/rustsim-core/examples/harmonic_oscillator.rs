//! Harmonic oscillator example using Connection and Simulation API
//!
//! Demonstrates feedback loops with integrators.
//!
//! System: d²x/dt² = -x
//! Solution: x(t) = x0*cos(t) + v0*sin(t)
//!
//! With x0=1, v0=0: x(t) = cos(t)

use rustsim_core::anyblock::AnyBlock;
use rustsim_core::blocks::{Amplifier, Integrator};
use rustsim_core::{Connection, Simulation};

/// Harmonic oscillator simulation
///
/// Block diagram:
/// ```text
///  ┌──────────────────────────────────┐
///  │                                  │
///  ▼                                  │
/// [gain: -1] ──► [velocity] ──► [position]
/// ```
fn main() {
    println!("Harmonic Oscillator - Connection & Simulation API Demo");
    println!("=======================================================");
    println!();
    println!("System: d²x/dt² = -x");
    println!("Initial: x(0) = 1.0, v(0) = 0.0");
    println!("Exact solution: x(t) = cos(t)");
    println!();

    // 1. Create blocks with descriptive names
    let gain = Amplifier::new(-1.0);     // acceleration = -position
    let velocity = Integrator::new(0.0); // v(0) = 0
    let position = Integrator::new(1.0); // x(0) = 1

    // 2. Collect blocks (indices match connection references)
    // Block 0: gain, Block 1: velocity, Block 2: position
    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(gain),
        Box::new(velocity),
        Box::new(position),
    ];

    // 3. Define connections (feedback loop)
    let connections = vec![
        Connection::from(2, 0).to(0, 0).build(),  // position -> gain
        Connection::from(0, 0).to(1, 0).build(),  // gain -> velocity (acceleration)
        Connection::from(1, 0).to(2, 0).build(),  // velocity -> position
    ];

    // 4. Create and configure simulation
    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.001);

    println!(
        "{:>10} {:>12} {:>12} {:>12}",
        "Time", "Position", "Exact", "Error"
    );
    println!("{:-<10} {:-<12} {:-<12} {:-<12}", "", "", "", "");

    // Simulate for 2π (one period)
    let period = 2.0 * std::f64::consts::PI;
    let dt = sim.dt();
    let steps = (period / dt) as usize;

    for i in 0..=steps {
        if i % (steps / 10) == 0 {
            let t = sim.time();
            let x = sim.get_output(2, 0);  // position
            let exact = t.cos();
            let error = (x - exact).abs();
            println!("{:10.4} {:12.6} {:12.6} {:12.2e}", t, x, exact, error);
        }
        sim.step();
    }

    let final_position = sim.get_output(2, 0);
    let final_velocity = sim.get_output(1, 0);

    println!();
    println!("After one period (t = 2π):");
    println!("  Position: {:.6} (should be ≈ 1.0)", final_position);
    println!("  Velocity: {:.6} (should be ≈ 0.0)", final_velocity);
    println!();

    // Energy conservation check
    let energy = 0.5 * (final_position.powi(2) + final_velocity.powi(2));
    println!("  Energy:   {:.6} (should be 0.5)", energy);
}
