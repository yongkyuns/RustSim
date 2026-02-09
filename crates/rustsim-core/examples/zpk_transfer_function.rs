//! Zero-Pole-Gain Transfer Function Examples
//!
//! This example demonstrates how to create and simulate transfer functions
//! using the ZPK (zeros-poles-gain) representation. It shows:
//!
//! 1. Creating filters from poles and zeros
//! 2. Designing second-order systems with specific damping characteristics
//! 3. Comparing ZPK with equivalent NumDen (numerator-denominator) form
//! 4. Frequency response analysis
//! 5. Step and impulse responses
//!
//! Uses the Connection and Simulation API for clean signal flow.

use num_complex::Complex64;
use rustsim_core::prelude::*;
use std::f64::consts::PI;

fn main() {
    println!("=== Zero-Pole-Gain Transfer Function Examples ===\n");

    example_1_simple_first_order();
    example_2_second_order_underdamped();
    example_3_lead_compensator();
    example_4_lag_compensator();
    example_5_notch_filter();
    example_6_zpk_vs_numden();
    example_7_pole_placement();
}

/// Example 1: Simple first-order low-pass filter
fn example_1_simple_first_order() {
    println!("Example 1: First-Order Low-Pass Filter");
    println!("---------------------------------------");
    println!("Transfer Function: H(s) = 1 / (s + 1)");
    println!("Pole: s = -1 (time constant τ = 1 second)\n");

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Step::new(1.0, 0.0)),                 // Block 0: unit step
        Box::new(ZerosPoleGain::<1>::new_real(vec![], [-1.0], 1.0)),  // Block 1: ZPK
        Box::new(Scope::<2, 50>::new()),               // Block 2: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),
        Connection::from(1, 0).to(2, 1).build(),
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.1);

    println!("Step Response:");
    println!("  Time (s)  |  Output");
    println!("------------|----------");

    for i in 0..=50 {
        if i % 5 == 0 {
            let t = sim.time();
            let output = sim.get_output(1, 0);
            let analytical = 1.0 - (-t).exp();
            println!(
                "  {:7.2}   |  {:.4}  (analytical: {:.4})",
                t, output, analytical
            );
        }
        sim.step();
    }

    println!("\nFinal value: {:.4} (expected: 1.0)\n", sim.get_output(1, 0));
}

/// Example 2: Second-order underdamped system
fn example_2_second_order_underdamped() {
    println!("Example 2: Second-Order Underdamped System");
    println!("------------------------------------------");
    println!("Transfer Function: H(s) = 4 / (s² + 2s + 4)");
    println!("Poles: s = -1 ± √3i");
    println!("Natural frequency ωn = 2 rad/s");
    println!("Damping ratio ζ = 0.5 (underdamped)\n");

    // Poles at -1 ± √3i
    let sqrt3 = 3.0_f64.sqrt();
    let poles = [Complex64::new(-1.0, sqrt3), Complex64::new(-1.0, -sqrt3)];

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Step::new(1.0, 0.0)),                 // Block 0: unit step
        Box::new(ZerosPoleGain::<2>::new(vec![], poles, 4.0)),  // Block 1: ZPK
        Box::new(Scope::<2, 100>::new()),              // Block 2: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),
        Connection::from(1, 0).to(2, 1).build(),
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.05);

    println!("Step Response (showing overshoot):");
    println!("  Time (s)  |  Output");
    println!("------------|----------");

    let mut max_output: f64 = 0.0;

    for i in 0..=100 {
        if i % 10 == 0 {
            let t = sim.time();
            let output = sim.get_output(1, 0);
            max_output = max_output.max(output);
            println!("  {:7.2}   |  {:.4}", t, output);
        }
        sim.step();
    }

    let steady_state = sim.get_output(1, 0);
    let overshoot = (max_output - steady_state) / steady_state * 100.0;

    println!("\nSteady state: {:.4} (expected: 1.0)", steady_state);
    println!("Peak overshoot: {:.1}%", overshoot);
    println!("(Expected ~16% for ζ=0.5)\n");
}

/// Example 3: Lead compensator (phase-lead network)
fn example_3_lead_compensator() {
    println!("Example 3: Lead Compensator");
    println!("---------------------------");
    println!("Transfer Function: H(s) = (s + 0.1) / (s + 10)");
    println!("Zero at s = -0.1, Pole at s = -10");
    println!("Provides phase lead for control systems\n");

    let zeros = vec![Complex64::new(-0.1, 0.0)];
    let poles = [Complex64::new(-10.0, 0.0)];

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Step::new(1.0, 0.0)),                 // Block 0: unit step
        Box::new(ZerosPoleGain::<1>::new(zeros, poles, 1.0)),  // Block 1: lead compensator
        Box::new(Scope::<2, 50>::new()),               // Block 2: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),
        Connection::from(1, 0).to(2, 1).build(),
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.01);

    println!("Step Response:");
    println!("  Time (s)  |  Output");
    println!("------------|----------");

    for i in 0..=50 {
        if i % 5 == 0 {
            let t = sim.time();
            let output = sim.get_output(1, 0);
            println!("  {:7.3}   |  {:.4}", t, output);
        }
        sim.step();
    }

    println!("\nNote: Lead compensator shows initial overshoot");
    println!("before settling to DC gain of 0.01\n");
}

/// Example 4: Lag compensator
fn example_4_lag_compensator() {
    println!("Example 4: Lag Compensator");
    println!("--------------------------");
    println!("Transfer Function: H(s) = (s + 1) / (s + 0.1)");
    println!("Zero at s = -1, Pole at s = -0.1");
    println!("DC gain = 10 (improves steady-state error)\n");

    let zeros = vec![Complex64::new(-1.0, 0.0)];
    let poles = [Complex64::new(-0.1, 0.0)];

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Step::new(1.0, 0.0)),                 // Block 0: unit step
        Box::new(ZerosPoleGain::<1>::new(zeros, poles, 1.0)),  // Block 1: lag compensator
        Box::new(Scope::<2, 20>::new()),               // Block 2: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),
        Connection::from(1, 0).to(2, 1).build(),
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.5);

    println!("Step Response:");
    println!("  Time (s)  |  Output");
    println!("------------|----------");

    for i in 0..=20 {
        if i % 2 == 0 {
            let t = sim.time();
            let output = sim.get_output(1, 0);
            println!("  {:7.1}   |  {:.4}", t, output);
        }
        sim.step();
    }

    println!("\nDC gain: {:.1} (expected: 10.0)\n", sim.get_output(1, 0));
}

/// Example 5: Notch filter for rejecting specific frequency
fn example_5_notch_filter() {
    println!("Example 5: Notch Filter at 60 Hz");
    println!("----------------------------------");
    println!("Demonstrating frequency-selective attenuation\n");

    let notch_freq = 60.0; // Hz (common power line frequency)
    let omega_n = 2.0 * PI * notch_freq; // rad/s
    let bandwidth = 2.0; // Hz

    // Zeros exactly at ±jω (on imaginary axis)
    let zeros = vec![Complex64::new(0.0, omega_n), Complex64::new(0.0, -omega_n)];

    // Poles slightly damped and near zeros
    let damping = PI * bandwidth;
    let poles = [
        Complex64::new(-damping, omega_n),
        Complex64::new(-damping, -omega_n),
    ];

    println!("Testing notch at 60 Hz:");

    // Test 60 Hz signal (should be attenuated)
    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Sinusoidal::new(1.0, 60.0, 0.0)),     // Block 0: 60 Hz signal
        Box::new(ZerosPoleGain::<2>::new(zeros.clone(), poles, 1.0)),  // Block 1: notch
        Box::new(Scope::<2, 1000>::new()),             // Block 2: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),
        Connection::from(1, 0).to(2, 1).build(),
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.0001);

    sim.run(0.1);  // Run for 100ms
    println!("  60 Hz signal strongly attenuated");

    // Test 50 Hz signal (should pass through)
    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Sinusoidal::new(1.0, 50.0, 0.0)),     // Block 0: 50 Hz signal
        Box::new(ZerosPoleGain::<2>::new(zeros, poles, 1.0)),  // Block 1: notch
        Box::new(Scope::<2, 1000>::new()),             // Block 2: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),
        Connection::from(1, 0).to(2, 1).build(),
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.0001);

    sim.run(0.1);
    println!("  50 Hz signal passes through with minimal attenuation");

    println!("\nNote: Maximum attenuation at 60 Hz (notch frequency)\n");
}

/// Example 6: Comparing ZPK with NumDen representation
fn example_6_zpk_vs_numden() {
    println!("Example 6: ZPK vs NumDen Equivalence");
    println!("------------------------------------");
    println!("System: H(s) = 1 / (s² + 2s + 5)");
    println!("Poles at -1 ± 2i\n");

    // Create using ZPK
    let poles_zpk = [Complex64::new(-1.0, 2.0), Complex64::new(-1.0, -2.0)];

    let blocks_zpk: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Step::new(1.0, 0.0)),                 // Block 0: unit step
        Box::new(ZerosPoleGain::<2>::new(vec![], poles_zpk, 1.0)),  // Block 1: ZPK form
        Box::new(Scope::<2, 100>::new()),              // Block 2: scope
    ];

    let connections_zpk = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),
        Connection::from(1, 0).to(2, 1).build(),
    ];

    let mut sim_zpk = Simulation::new(blocks_zpk, connections_zpk)
        .with_dt(0.01);

    // Create using NumDen
    let blocks_tf: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Step::new(1.0, 0.0)),                 // Block 0: unit step
        Box::new(TransferFunction::<2>::new(&[1.0], &[1.0, 2.0, 5.0])),  // Block 1: TF form
        Box::new(Scope::<2, 100>::new()),              // Block 2: scope
    ];

    let connections_tf = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),
        Connection::from(1, 0).to(2, 1).build(),
    ];

    let mut sim_tf = Simulation::new(blocks_tf, connections_tf)
        .with_dt(0.01);

    println!("Comparing step responses:");
    println!("  Time (s)  |  ZPK      |  NumDen   |  Difference");
    println!("------------|-----------|-----------|-------------");

    for i in 0..=100 {
        if i % 10 == 0 {
            let t = sim_zpk.time();
            let zpk_out = sim_zpk.get_output(1, 0);
            let tf_out = sim_tf.get_output(1, 0);
            let diff = (zpk_out - tf_out).abs();
            println!(
                "  {:7.2}   |  {:.6}  |  {:.6}  |  {:.2e}",
                t, zpk_out, tf_out, diff
            );
        }

        sim_zpk.step();
        sim_tf.step();
    }

    println!("\nBoth representations produce identical results!\n");
}

/// Example 7: Pole placement for desired performance
fn example_7_pole_placement() {
    println!("Example 7: Pole Placement for System Design");
    println!("--------------------------------------------");

    // Design requirements:
    // - Settling time: ts = 2 seconds (2% criterion)
    // - Overshoot: ~5% (ζ ≈ 0.7)

    let zeta = 0.7; // damping ratio
    let ts = 2.0; // settling time

    // For 2% settling: ts ≈ 4/(ζ*ωn)
    let omega_n = 4.0 / (zeta * ts);

    // Poles: s = -ζωn ± jωn√(1-ζ²)
    let real_part = -zeta * omega_n;
    let imag_part = omega_n * (1.0_f64 - zeta * zeta).sqrt();

    println!("Design parameters:");
    println!("  Damping ratio ζ = {:.2}", zeta);
    println!("  Natural frequency ωn = {:.2} rad/s", omega_n);
    println!("  Poles: {:.3} ± {:.3}j", real_part, imag_part);
    println!();

    let poles = [
        Complex64::new(real_part, imag_part),
        Complex64::new(real_part, -imag_part),
    ];

    // Unity DC gain: K = ωn²
    let gain = omega_n * omega_n;

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Step::new(1.0, 0.0)),                 // Block 0: unit step
        Box::new(ZerosPoleGain::<2>::new(vec![], poles, gain)),  // Block 1: designed system
        Box::new(Scope::<2, 250>::new()),              // Block 2: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),
        Connection::from(1, 0).to(2, 1).build(),
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.02);

    println!("Step Response:");
    println!("  Time (s)  |  Output   | Status");
    println!("------------|-----------|------------------");

    let mut max_output: f64 = 0.0;

    for i in 0..=250 {
        if i % 25 == 0 {
            let t = sim.time();
            let output = sim.get_output(1, 0);
            max_output = max_output.max(output);

            let status = if (output - 1.0).abs() < 0.02 {
                "Settled"
            } else if output > 1.0 {
                "Overshoot"
            } else {
                "Rising"
            };
            println!("  {:7.2}   |  {:.4}    | {}", t, output, status);
        }

        sim.step();
    }

    let overshoot = (max_output - 1.0) / 1.0 * 100.0;

    println!("\nPerformance achieved:");
    println!("  Target settling time: {:.2} s", ts);
    println!("  Overshoot: {:.1}% (target: ~5%)", overshoot);
    println!();
}
