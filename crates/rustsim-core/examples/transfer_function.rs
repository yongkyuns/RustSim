//! Transfer Function examples demonstrating common use cases
//!
//! Shows how to create and simulate various transfer functions using
//! the Connection and Simulation API for clean signal flow.

use rustsim_core::prelude::*;

fn main() {
    println!("Transfer Function Block Examples");
    println!("================================\n");

    // Example 1: First-order low-pass filter
    example_1_first_order_lowpass();

    // Example 2: Second-order system
    example_2_second_order();

    // Example 3: Lead compensator
    example_3_lead_compensator();

    // Example 4: Pure integrator
    example_4_integrator();

    // Example 5: Transfer function in feedback loop
    example_5_feedback_control();
}

fn example_1_first_order_lowpass() {
    println!("Example 1: First-Order Low-Pass Filter");
    println!("---------------------------------------");
    println!("H(s) = 1/(s+1), Time constant τ = 1");
    println!();

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Step::new(1.0, 0.0)),                 // Block 0: unit step at t=0
        Box::new(TransferFunction::<1>::new(&[1.0], &[1.0, 1.0])),  // Block 1: H(s)
        Box::new(Scope::<2, 500>::new()),              // Block 2: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),  // Step -> TF, scope[0]
        Connection::from(1, 0).to(2, 1).build(),           // TF -> scope[1]
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.01);

    println!(
        "{:>8} {:>12} {:>12} {:>12}",
        "Time", "Output", "Exact", "Error"
    );
    println!("{:-<8} {:-<12} {:-<12} {:-<12}", "", "", "", "");

    // Exact solution: y(t) = 1 - e^(-t)
    for i in 0..=500 {
        if i % 100 == 0 {
            let t = sim.time();
            let output = sim.get_output(1, 0);
            let exact = 1.0 - (-t).exp();
            let error = (output - exact).abs();
            println!("{:8.2} {:12.6} {:12.6} {:12.2e}", t, output, exact, error);
        }
        sim.step();
    }
    println!();
}

fn example_2_second_order() {
    println!("Example 2: Second-Order System");
    println!("-------------------------------");
    println!("H(s) = 1/(s^2 + 2s + 1) = 1/(s+1)^2");
    println!();

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Step::new(1.0, 0.0)),                 // Block 0: unit step
        Box::new(TransferFunction::<2>::new(&[1.0], &[1.0, 2.0, 1.0])),  // Block 1: H(s)
        Box::new(Scope::<2, 500>::new()),              // Block 2: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),
        Connection::from(1, 0).to(2, 1).build(),
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.01);

    println!(
        "{:>8} {:>12} {:>12} {:>12}",
        "Time", "Output", "Exact", "Error"
    );
    println!("{:-<8} {:-<12} {:-<12} {:-<12}", "", "", "", "");

    // Exact solution: y(t) = 1 - (1+t)e^(-t)
    for i in 0..=500 {
        if i % 100 == 0 {
            let t = sim.time();
            let output = sim.get_output(1, 0);
            let exact = 1.0 - (1.0 + t) * (-t).exp();
            let error = (output - exact).abs();
            println!("{:8.2} {:12.6} {:12.6} {:12.2e}", t, output, exact, error);
        }
        sim.step();
    }
    println!();
}

fn example_3_lead_compensator() {
    println!("Example 3: Lead Compensator");
    println!("----------------------------");
    println!("H(s) = (s+2)/(s+10)");
    println!();

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Step::new(1.0, 0.0)),                 // Block 0: unit step
        Box::new(TransferFunction::<1>::new(&[1.0, 2.0], &[1.0, 10.0])),  // Block 1: lead
        Box::new(Scope::<2, 1000>::new()),             // Block 2: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),
        Connection::from(1, 0).to(2, 1).build(),
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.001);

    println!(
        "{:>8} {:>12} {:>12} {:>12}",
        "Time", "Output", "Exact", "Error"
    );
    println!("{:-<8} {:-<12} {:-<12} {:-<12}", "", "", "", "");

    // Exact solution: y(t) = 0.2 + 0.8*e^(-10t)
    for i in 0..=1000 {
        if i % 100 == 0 {
            let t = sim.time();
            let output = sim.get_output(1, 0);
            let exact = 0.2 + 0.8 * (-10.0 * t).exp();
            let error = (output - exact).abs();
            println!("{:8.2} {:12.6} {:12.6} {:12.2e}", t, output, exact, error);
        }
        sim.step();
    }
    println!();
}

fn example_4_integrator() {
    println!("Example 4: Pure Integrator");
    println!("--------------------------");
    println!("H(s) = 1/s");
    println!();

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Constant::new(1.0)),                  // Block 0: constant input
        Box::new(TransferFunction::<1>::new(&[1.0], &[1.0, 0.0])),  // Block 1: integrator
        Box::new(Scope::<2, 500>::new()),              // Block 2: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),
        Connection::from(1, 0).to(2, 1).build(),
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.01);

    println!(
        "{:>8} {:>12} {:>12} {:>12}",
        "Time", "Output", "Exact", "Error"
    );
    println!("{:-<8} {:-<12} {:-<12} {:-<12}", "", "", "", "");

    // Exact solution: y(t) = t
    for i in 0..=500 {
        if i % 100 == 0 {
            let t = sim.time();
            let output = sim.get_output(1, 0);
            let exact = t;
            let error = (output - exact).abs();
            println!("{:8.2} {:12.6} {:12.6} {:12.2e}", t, output, exact, error);
        }
        sim.step();
    }
    println!();
}

fn example_5_feedback_control() {
    println!("Example 5: Transfer Function in Feedback Loop");
    println!("----------------------------------------------");
    println!("Closed-loop system with first-order plant");
    println!("Plant: G(s) = 1/(s+1)");
    println!("Unity feedback with proportional gain K=2");
    println!();

    // Setpoint -> (+) -> Gain -> Plant -> Output
    //              ^                        |
    //              |------------------------|
    //                    (feedback)

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Step::new(1.0, 0.0)),                 // Block 0: setpoint
        Box::new(Adder::<2>::new()),                   // Block 1: error = setpoint - output
        Box::new(Amplifier::new(2.0)),                 // Block 2: controller gain
        Box::new(TransferFunction::<1>::new(&[1.0], &[1.0, 1.0])),  // Block 3: plant
        Box::new(Scope::<3, 500>::new()),              // Block 4: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(4, 0).build(),  // Setpoint -> adder[0], scope[0]
        Connection::from(1, 0).to(2, 0).build(),           // Error -> gain
        Connection::from(2, 0).to(3, 0).build(),           // Gain -> plant
        Connection::from(3, 0).to(1, 1).to(4, 1).build(),  // Output -> adder[1] (feedback), scope[1]
        Connection::from(1, 0).to(4, 2).build(),           // Error -> scope[2]
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.01);

    println!(
        "{:>8} {:>12} {:>12} {:>12}",
        "Time", "Setpoint", "Output", "Error"
    );
    println!("{:-<8} {:-<12} {:-<12} {:-<12}", "", "", "", "");

    for i in 0..=500 {
        if i % 50 == 0 {
            let t = sim.time();
            let setpoint = sim.get_output(0, 0);
            let output = sim.get_output(3, 0);
            let error = sim.get_output(1, 0);
            println!("{:8.2} {:12.6} {:12.6} {:12.6}", t, setpoint, output, error);
        }
        sim.step();
    }

    println!();
    println!("Closed-loop settling:");
    println!("  Final output: {:.6}", sim.get_output(3, 0));
    println!("  Final error:  {:.6}", sim.get_output(1, 0));
    println!();
}
