///! FIR Filter Example
///!
///! Demonstrates using the FIR filter block with different configurations:
///! - Moving average (low-pass)
///! - Differentiator (high-pass)
///! - Custom coefficients
///!
///! Uses the new Connection and Simulation API for cleaner signal flow.

use rustsim_core::prelude::*;

fn main() {
    println!("=== FIR Filter Examples ===\n");

    // Example 1: Moving average filter (3-tap low-pass)
    println!("1. Moving Average Filter [0.25, 0.5, 0.25]");
    println!("   Sinusoidal input with high-frequency noise");

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Sinusoidal::new(1.0, 10.0, 0.0)),    // Block 0: 10 Hz sine
        Box::new(FIR::<3>::new([0.25, 0.5, 0.25], 1.0, 0.0)),  // Block 1: MA filter
        Box::new(Scope::<2, 1000>::new()),             // Block 2: scope (input and output)
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),  // Source -> filter, scope[0]
        Connection::from(1, 0).to(2, 1).build(),           // Filter -> scope[1]
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.001);

    sim.run(0.1);  // Run for 100ms

    println!("   Simulated 100ms with dt=1ms");
    println!("   Input peak: ~1.0, Output smoothed\n");

    // Example 2: Differentiator [1.0, -1.0]
    println!("2. Differentiator FIR [1.0, -1.0]");
    println!("   Ramp input -> constant output");

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Ramp::new(1.0, 0.0)),                 // Block 0: ramp with slope=1
        Box::new(FIR::<2>::new([1.0, -1.0], 1.0, 0.0)), // Block 1: differentiator
        Box::new(Scope::<2, 100>::new()),              // Block 2: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),  // Ramp -> filter, scope[0]
        Connection::from(1, 0).to(2, 1).build(),           // Filter -> scope[1]
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.01);

    sim.run(1.0);  // Run for 1 second

    println!("   Ramp derivative ≈ 1.0 (slope)");
    println!("   Final output: {:.3}\n", sim.get_output(1, 0));

    // Example 3: Custom band-pass FIR
    println!("3. Custom Band-Pass FIR");
    println!("   Multi-frequency input");

    // Create a composite signal: low freq + mid freq + high freq
    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Sinusoidal::new(1.0, 5.0, 0.0)),      // Block 0: 5 Hz
        Box::new(Sinusoidal::new(0.5, 20.0, 0.0)),     // Block 1: 20 Hz
        Box::new(Sinusoidal::new(0.3, 50.0, 0.0)),     // Block 2: 50 Hz
        Box::new(Adder::<3>::new()),                   // Block 3: sum signals
        Box::new(FIR::<5>::new([0.1, 0.2, 0.4, 0.2, 0.1], 1.0, 0.0)), // Block 4: FIR
        Box::new(Scope::<2, 2000>::new()),             // Block 5: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(3, 0).build(),       // 5 Hz -> adder
        Connection::from(1, 0).to(3, 1).build(),       // 20 Hz -> adder
        Connection::from(2, 0).to(3, 2).build(),       // 50 Hz -> adder
        Connection::from(3, 0).to(4, 0).to(5, 0).build(),  // Composite -> filter, scope[0]
        Connection::from(4, 0).to(5, 1).build(),       // Filter -> scope[1]
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.0001);  // High sample rate for filters

    sim.run(0.2);  // 200ms

    println!("   Simulated composite signal through FIR");
    println!("   Input contains 5Hz, 20Hz, 50Hz components");
    println!("   Filter smooths high frequencies\n");

    // Example 4: Step response of low-pass FIR
    println!("4. Step Response - Low-Pass FIR");

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Step::new(1.0, 0.0)),                 // Block 0: unit step at t=0
        Box::new(FIR::<3>::new([0.25, 0.5, 0.25], 1.0, 0.0)), // Block 1: MA filter
        Box::new(Scope::<2, 200>::new()),              // Block 2: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),
        Connection::from(1, 0).to(2, 1).build(),
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.01);

    sim.run(2.0);

    println!("   Step input at t=0");
    println!("   Final filtered output: {:.3}", sim.get_output(1, 0));
    println!("   (Should settle to 1.0)\n");

    println!("=== Example Complete ===");
}
