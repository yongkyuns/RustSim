//! Advanced example: PID control system using PathSim-style API
//!
//! Demonstrates a complete closed-loop control system with:
//! - Reference signal (setpoint)
//! - PID controller
//! - Plant (first-order system)
//! - Feedback loop
//!
//! This example showcases the clarity of the PathSim-style API for
//! complex block diagrams.

use rustsim_core::anyblock::AnyBlock;
use rustsim_core::blocks::{Adder, Amplifier, Constant, Integrator, Scope};
use rustsim_core::{BlockExt, ConnectionDef, SimulationBuilder};

fn main() {
    println!("PID Control System - PathSim-Style API");
    println!("=======================================");
    println!();

    // === System Parameters ===
    let setpoint = 10.0; // Desired output
    let kp = 2.0; // Proportional gain
    let ki = 1.0; // Integral gain
    let plant_time_constant = 0.5; // Plant dynamics

    // === Create Blocks ===

    // Reference input
    let reference = Constant::new(setpoint);

    // Error calculation (setpoint - feedback)
    let error_adder = Adder::<2>::with_weights([1.0, -1.0]); // Positive for reference, negative for feedback

    // PID controller components
    let proportional = Amplifier::new(kp);
    let integral = Integrator::new(0.0);
    let integral_gain = Amplifier::new(ki);

    // Control signal combiner
    let control_adder = Adder::<2>::new(); // Add proportional + integral

    // Plant (first-order system: 1/(tau*s + 1))
    // Modeled as: dx/dt = (u - x) / tau
    let plant_feedback = Amplifier::new(-1.0 / plant_time_constant);
    let plant_input = Amplifier::new(1.0 / plant_time_constant);
    let plant_adder = Adder::<2>::new();
    let plant_integrator = Integrator::new(0.0);

    // Scope to record output
    let scope = Scope::<1, 10000>::new();

    // === Define Connections ===
    let connections = vec![
        // Reference to error adder (port 0)
        ConnectionDef::new(&reference, error_adder.port(0)),

        // Error signal to controller
        ConnectionDef::new(&error_adder, &proportional)
            .to(&integral), // Fan-out to P and I terms

        // Integral term processing
        ConnectionDef::new(&integral, &integral_gain),

        // Combine control signals (P + I)
        ConnectionDef::new(&proportional, control_adder.port(0)),
        ConnectionDef::new(&integral_gain, control_adder.port(1)),

        // Control signal to plant
        ConnectionDef::new(&control_adder, &plant_input),
        ConnectionDef::new(&plant_input, plant_adder.port(0)),

        // Plant integrator dynamics
        ConnectionDef::new(&plant_adder, &plant_integrator),

        // Plant feedback
        ConnectionDef::new(&plant_integrator, &plant_feedback)
            .to(scope.port(0))                    // Record output
            .to(error_adder.port(1)),            // Feedback to error calculation

        ConnectionDef::new(&plant_feedback, plant_adder.port(1)),
    ];

    // === Collect Blocks (in order!) ===
    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(reference),
        Box::new(error_adder),
        Box::new(proportional),
        Box::new(integral),
        Box::new(integral_gain),
        Box::new(control_adder),
        Box::new(plant_input),
        Box::new(plant_adder),
        Box::new(plant_integrator),
        Box::new(plant_feedback),
        Box::new(scope),
    ];

    // === Build and Run Simulation ===
    let mut sim = SimulationBuilder::build_from_defs(blocks, connections)
        .with_dt(0.001);

    println!("System Configuration:");
    println!("  Setpoint: {}", setpoint);
    println!("  Kp (proportional gain): {}", kp);
    println!("  Ki (integral gain): {}", ki);
    println!("  Plant time constant: {}", plant_time_constant);
    println!();

    println!("{:>10} {:>12} {:>12} {:>12}", "Time", "Output", "Error", "Error %");
    println!("{:-<10} {:-<12} {:-<12} {:-<12}", "", "", "", "");

    // Run for 10 seconds
    let duration = 10.0;
    let dt = sim.dt();
    let steps = (duration / dt) as usize;
    let print_interval = steps / 20;

    for i in 0..=steps {
        sim.step();

        if i % print_interval == 0 {
            let t = sim.time();

            // Get plant output (from plant_integrator)
            let output = sim.get_output(8, 0);

            let error = setpoint - output;
            let error_pct = (error / setpoint) * 100.0;

            println!(
                "{:10.4} {:12.6} {:12.6} {:12.2}",
                t, output, error, error_pct
            );
        }
    }

    // Get final values
    let final_output = sim.get_output(8, 0);
    let final_error = setpoint - final_output;
    let settling_error_pct = (final_error / setpoint) * 100.0;

    println!();
    println!("Final Results:");
    println!("  Output: {:.6}", final_output);
    println!("  Error: {:.6}", final_error);
    println!("  Settling error: {:.3}%", settling_error_pct);

    // Performance metrics
    if settling_error_pct.abs() < 2.0 {
        println!("  Status: CONVERGED (error < 2%)");
    } else {
        println!("  Status: NOT CONVERGED (error >= 2%)");
    }

    println!();
    println!("This example demonstrates:");
    println!("  1. Complex block diagrams are readable with PathSim-style API");
    println!("  2. Feedback loops work correctly (broken by integrators)");
    println!("  3. Fan-out connections (one output to multiple inputs)");
    println!("  4. Explicit port specification for multi-input blocks");
}
