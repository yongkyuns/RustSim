//! Demonstration of the Simulation runner
//!
//! This example shows how to use the Simulation type to create and run
//! dynamic block-based simulations with automatic execution ordering.
//!
//! The example creates a simple control system:
//! - Step input (reference signal)
//! - Feedback loop with integrator
//! - Amplifier in the feedback path

use rustsim_core::prelude::*;

fn main() {
    println!("Simulation Runner Demo");
    println!("=====================\n");

    // Example 1: Simple integration
    println!("Example 1: Simple Integration (Constant -> Integrator)");
    println!("-------------------------------------------------------");
    simple_integration_example();

    println!("\nExample 2: Feedback Control System");
    println!("-----------------------------------");
    feedback_example();

    println!("\nExample 3: Multi-output Fan-out");
    println!("--------------------------------");
    fanout_example();
}

fn simple_integration_example() {
    // Create blocks
    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Constant::new(2.0)),      // Block 0: constant input
        Box::new(Integrator::new(0.0)),    // Block 1: integrator
    ];

    // Create connections
    let connections = vec![
        Connection::from(0, 0).to(1, 0).build(), // Constant -> Integrator
    ];

    // Create simulation
    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.01); // 10ms time step

    println!("Initial state:");
    println!("  Time: {:.3} s", sim.time());
    println!("  Integrator output: {:.3}", sim.get_output(1, 0));

    // Run for 1 second
    sim.run(1.0);

    println!("\nAfter 1 second:");
    println!("  Time: {:.3} s", sim.time());
    println!("  Integrator output: {:.3}", sim.get_output(1, 0));
    println!("  Expected: 2.0 (integral of constant 2.0 over 1 second)");

    // Reset and run again
    sim.reset();
    println!("\nAfter reset:");
    println!("  Time: {:.3} s", sim.time());
    println!("  Integrator output: {:.3}", sim.get_output(1, 0));
}

fn feedback_example() {
    // Create a feedback system:
    // Step -> (+) -> Amplifier -> Integrator -> (feedback to +)
    //          ^                                      |
    //          |--------------------------------------|

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Step::new(1.0, 0.0)),     // Block 0: step input at t=0
        Box::new(Adder::<2>::new()),       // Block 1: adder (combines reference and feedback)
        Box::new(Amplifier::new(0.5)),     // Block 2: gain in forward path
        Box::new(Integrator::new(0.0)),    // Block 3: integrator (plant)
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).build(), // Step -> Adder input 0
        Connection::from(1, 0).to(2, 0).build(), // Adder -> Amplifier
        Connection::from(2, 0).to(3, 0).build(), // Amplifier -> Integrator
        Connection::from(3, 0).to(1, 1).build(), // Integrator -> Adder input 1 (feedback)
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.01);

    println!("Simulating feedback control system...");
    println!("Execution order: {:?}", sim.execution_order());

    println!("\n{:>8} {:>12} {:>12}", "Time", "Reference", "Output");
    println!("{:-<8} {:-<12} {:-<12}", "", "", "");

    let duration = 5.0;
    let print_interval = 0.5;
    let mut next_print = 0.0;

    while sim.time() < duration {
        if sim.time() >= next_print {
            println!(
                "{:8.2} {:12.4} {:12.4}",
                sim.time(),
                sim.get_output(0, 0), // Reference (step)
                sim.get_output(3, 0), // System output (integrator)
            );
            next_print += print_interval;
        }
        sim.step();
    }
}

fn fanout_example() {
    // Demonstrate fan-out: one signal feeding multiple blocks
    // Sinusoid -> Amplifier(2x)
    //          -> Amplifier(3x)
    //          -> Amplifier(5x)

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Sinusoidal::new(1.0, 1.0, 0.0)), // Block 0: 1Hz sine wave
        Box::new(Amplifier::new(2.0)),             // Block 1: 2x gain
        Box::new(Amplifier::new(3.0)),             // Block 2: 3x gain
        Box::new(Amplifier::new(5.0)),             // Block 3: 5x gain
    ];

    // Single connection with multiple targets (fan-out)
    let connections = vec![
        Connection::from(0, 0)
            .to(1, 0)
            .to(2, 0)
            .to(3, 0)
            .build(),
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.01);

    println!("Fan-out example: one sinusoid feeding three amplifiers");
    println!("\n{:>8} {:>10} {:>10} {:>10} {:>10}",
             "Time", "Source", "Amp 2x", "Amp 3x", "Amp 5x");
    println!("{:-<8} {:-<10} {:-<10} {:-<10} {:-<10}", "", "", "", "", "");

    for i in 0..=20 {
        if i % 5 == 0 {
            println!(
                "{:8.2} {:10.4} {:10.4} {:10.4} {:10.4}",
                sim.time(),
                sim.get_output(0, 0),
                sim.get_output(1, 0),
                sim.get_output(2, 0),
                sim.get_output(3, 0),
            );
        }
        sim.step();
    }

    // Verify fan-out correctness at current time
    let source = sim.get_output(0, 0);
    let amp2 = sim.get_output(1, 0);
    let amp3 = sim.get_output(2, 0);
    let amp5 = sim.get_output(3, 0);

    println!("\nVerification:");
    println!("  Source: {:.6}", source);
    println!("  2x amplified: {:.6} (expected: {:.6})", amp2, 2.0 * source);
    println!("  3x amplified: {:.6} (expected: {:.6})", amp3, 3.0 * source);
    println!("  5x amplified: {:.6} (expected: {:.6})", amp5, 5.0 * source);

    assert!((amp2 - 2.0 * source).abs() < 1e-10);
    assert!((amp3 - 3.0 * source).abs() < 1e-10);
    assert!((amp5 - 5.0 * source).abs() < 1e-10);
    println!("  ✓ All outputs match expected values!");
}
