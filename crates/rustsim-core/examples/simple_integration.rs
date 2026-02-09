//! Simple integration example using Connection and Simulation API
//!
//! Demonstrates basic integrator usage: integrate constant input
//!
//! System: dy/dt = 1, y(0) = 0
//! Solution: y(t) = t

use rustsim_core::anyblock::AnyBlock;
use rustsim_core::blocks::{Constant, Integrator};
use rustsim_core::{Connection, Simulation};

fn main() {
    println!("Simple Integration - Connection & Simulation API Demo");
    println!("=====================================================");
    println!();
    println!("System: dy/dt = 1, y(0) = 0");
    println!("Exact solution: y(t) = t");
    println!();

    // 1. Create blocks with descriptive names
    let constant = Constant::new(1.0);       // dy/dt = 1
    let integrator = Integrator::new(0.0);   // y(0) = 0

    // 2. Collect blocks (indices match connection references)
    // Block 0: constant, Block 1: integrator
    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(constant),
        Box::new(integrator),
    ];

    // 3. Define connections
    let connections = vec![
        Connection::from(0, 0).to(1, 0).build(),  // constant -> integrator
    ];

    // 4. Create simulation
    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.01);

    println!(
        "{:>10} {:>12} {:>12} {:>12}",
        "Time", "Output", "Exact", "Error"
    );
    println!("{:-<10} {:-<12} {:-<12} {:-<12}", "", "", "", "");

    // Simulate for 5 seconds
    let duration = 5.0;
    let dt = sim.dt();
    let steps = (duration / dt) as usize;

    for i in 0..=steps {
        if i % (steps / 10) == 0 {
            let t = sim.time();
            let y = sim.get_output(1, 0);
            let exact = t;
            let error = (y - exact).abs();
            println!("{:10.4} {:12.6} {:12.6} {:12.2e}", t, y, exact, error);
        }
        sim.step();
    }

    println!();
    println!("Final value: {:.6} (expected 5.0)", sim.get_output(1, 0));
    println!("Final error: {:.2e}", (sim.get_output(1, 0) - 5.0).abs());
}
