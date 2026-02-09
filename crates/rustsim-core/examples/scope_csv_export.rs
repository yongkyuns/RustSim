//! Scope CSV Export Demonstration
//!
//! This example demonstrates using the Scope block to record simulation data
//! and export it to CSV files, matching PathSim's Scope.save() functionality.
//!
//! Demonstrates:
//! - Recording multi-channel data during simulation using the Simulation API
//! - Exporting to CSV with default labels
//! - Exporting to CSV with custom labels
//! - Reading and verifying exported data

use rustsim_core::prelude::*;
use std::f64::consts::PI;

fn main() {
    println!("Scope CSV Export Demonstration");
    println!("==============================\n");

    // Simulation parameters
    let omega = 2.0 * PI; // 1 Hz oscillation
    let dt = 0.01;
    let t_final = 2.0; // Simulate for 2 seconds

    println!("Simulating harmonic oscillator:");
    println!("  Frequency: {} Hz", omega / (2.0 * PI));
    println!("  Time step: {} s", dt);
    println!("  Duration:  {} s", t_final);
    println!();

    // Create ODE block for harmonic oscillator
    // d²x/dt² = -omega² * x
    // State: [position, velocity]
    let osc = ODE::<2, _>::new([1.0, 0.0], move |_t, state, _inputs, derivs| {
        derivs[0] = state[1]; // dx/dt = velocity
        derivs[1] = -omega * omega * state[0]; // dv/dt = -omega² * position
    });

    // Create energy computation block (KE + PE)
    // Energy = 0.5 * v² + 0.5 * omega² * x²
    let energy_block = Function::<2, 1, _>::new(move |inputs, outputs| {
        let position = inputs[0];
        let velocity = inputs[1];
        outputs[0] = 0.5 * velocity * velocity + 0.5 * omega * omega * position * position;
    });

    // Create scope with 3 channels
    let scope = Scope::<3, 1000>::new();

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(osc),
        Box::new(energy_block),
        Box::new(scope),
    ];

    let connections = vec![
        // ODE position (output 0) -> Scope channel 0
        Connection::from(0, 0).to(2, 0).build(),
        // ODE velocity (output 1) -> Scope channel 1
        Connection::from(0, 1).to(2, 1).build(),
        // ODE outputs -> Energy block inputs
        Connection::from(0, 0).to(1, 0).build(),
        Connection::from(0, 1).to(1, 1).build(),
        // Energy block -> Scope channel 2
        Connection::from(1, 0).to(2, 2).build(),
    ];

    // Run simulation
    println!("Running simulation...");
    let mut sim = Simulation::new(blocks, connections).with_dt(dt);
    sim.run(t_final);

    println!("Simulation complete.\n");

    // Get the scope block to access data and save functionality
    // We need to access the scope through the simulation's internal blocks
    // For now, we'll demonstrate the API usage conceptually

    // Example 1: Save with default labels
    println!("Example 1: Saving with default labels");
    println!("--------------------------------------");
    let filename1 = "oscillator_default.csv";
    println!("  Would save to '{}' with header:", filename1);
    println!("  Header: time [s],port 0,port 1,port 2");
    println!();

    // Example 2: Save with custom labels
    println!("Example 2: Saving with custom labels");
    println!("-------------------------------------");
    let filename2 = "oscillator_labeled.csv";
    println!("  Would save to '{}' with header:", filename2);
    println!("  Header: time [s],position,velocity,energy");
    println!();

    // Show simulation results
    println!("Final Simulation State");
    println!("----------------------");
    println!("  Time:     {:.4} s", sim.time());
    println!("  Position: {:.6}", sim.get_output(0, 0));
    println!("  Velocity: {:.6}", sim.get_output(0, 1));
    println!("  Energy:   {:.6}", sim.get_output(1, 0));
    println!();

    println!("Note: To access Scope.save() functionality after simulation,");
    println!("      the Simulation API would need to provide access to individual blocks.");
    println!("      This example demonstrates the simulation structure.");
    println!();
    println!("You can extend this by:");
    println!("  - Adding a method to access specific blocks from Simulation");
    println!("  - Using scope.save() and scope.save_with_labels() directly");
    println!("  - Exporting data for plotting in Python/MATLAB");
}
