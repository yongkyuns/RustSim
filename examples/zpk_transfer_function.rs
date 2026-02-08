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

use num_complex::Complex64;
use rustsim::blocks::{TransferFunction, ZerosPoleGain};
use rustsim::Block;
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

    let mut zpk = ZerosPoleGain::<1>::new_real(vec![], [-1.0], 1.0);

    // Apply step input
    zpk.set_input(0, 1.0);
    zpk.update(0.0);

    println!("Step Response:");
    println!("  Time (s)  |  Output");
    println!("------------|----------");

    let dt = 0.1;
    for i in 0..=50 {
        let t = i as f64 * dt;
        if i % 5 == 0 {
            let analytical = 1.0 - (-t).exp();
            println!(
                "  {:7.2}   |  {:.4}  (analytical: {:.4})",
                t,
                zpk.get_output(0),
                analytical
            );
        }
        zpk.step(t, dt);
    }

    println!("\nFinal value: {:.4} (expected: 1.0)\n", zpk.get_output(0));
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
    let poles = [
        Complex64::new(-1.0, sqrt3),
        Complex64::new(-1.0, -sqrt3),
    ];
    let mut zpk = ZerosPoleGain::<2>::new(vec![], poles, 4.0);

    // Apply step input
    zpk.set_input(0, 1.0);
    zpk.update(0.0);

    println!("Step Response (showing overshoot):");
    println!("  Time (s)  |  Output");
    println!("------------|----------");

    let dt = 0.05;
    let mut max_output: f64 = 0.0;

    for i in 0..=100 {
        let t = i as f64 * dt;
        max_output = max_output.max(zpk.get_output(0));

        if i % 10 == 0 {
            println!("  {:7.2}   |  {:.4}", t, zpk.get_output(0));
        }
        zpk.step(t, dt);
    }

    let steady_state = zpk.get_output(0);
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
    let mut zpk = ZerosPoleGain::<1>::new(zeros, poles, 1.0);

    // Apply step input to see transient behavior
    zpk.set_input(0, 1.0);
    zpk.update(0.0);

    println!("Step Response:");
    println!("  Time (s)  |  Output");
    println!("------------|----------");

    let dt = 0.01;

    for i in 0..=50 {
        let t = i as f64 * dt;
        if i % 5 == 0 {
            println!("  {:7.3}   |  {:.4}", t, zpk.get_output(0));
        }
        zpk.step(t, dt);
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
    let mut zpk = ZerosPoleGain::<1>::new(zeros, poles, 1.0);

    // Apply step input
    zpk.set_input(0, 1.0);
    zpk.update(0.0);

    println!("Step Response:");
    println!("  Time (s)  |  Output");
    println!("------------|----------");

    let dt = 0.5;

    for i in 0..=20 {
        let t = i as f64 * dt;
        if i % 2 == 0 {
            println!("  {:7.1}   |  {:.4}", t, zpk.get_output(0));
        }
        zpk.step(t, dt);
    }

    println!("\nDC gain: {:.1} (expected: 10.0)\n", zpk.get_output(0));
}

/// Example 5: Notch filter for rejecting specific frequency
fn example_5_notch_filter() {
    println!("Example 5: Notch Filter at 60 Hz");
    println!("----------------------------------");

    let notch_freq = 60.0; // Hz (common power line frequency)
    let omega_n = 2.0 * PI * notch_freq; // rad/s
    let bandwidth = 2.0; // Hz

    // Zeros exactly at ±jω (on imaginary axis)
    let zeros = vec![
        Complex64::new(0.0, omega_n),
        Complex64::new(0.0, -omega_n),
    ];

    // Poles slightly damped and near zeros
    let damping = PI * bandwidth;
    let poles = [
        Complex64::new(-damping, omega_n),
        Complex64::new(-damping, -omega_n),
    ];

    let mut zpk = ZerosPoleGain::<2>::new(zeros, poles, 1.0);

    println!("Testing attenuation at different frequencies:\n");

    // Test at different frequencies
    let test_freqs = vec![50.0, 55.0, 60.0, 65.0, 70.0, 100.0];

    for &freq in &test_freqs {
        zpk.reset();

        let dt = 0.0001;
        let duration = 0.1; // Let it settle
        let steps = (duration / dt) as usize;

        // Apply sinusoidal input
        for i in 0..steps {
            let t = i as f64 * dt;
            let input = (2.0 * PI * freq * t).sin();
            zpk.set_input(0, input);
            zpk.update(t);
            zpk.step(t, dt);
        }

        // Measure amplitude in steady state
        let mut max_val = f64::NEG_INFINITY;
        let mut min_val = f64::INFINITY;

        for i in steps..steps + 1000 {
            let t = i as f64 * dt;
            let input = (2.0 * PI * freq * t).sin();
            zpk.set_input(0, input);
            zpk.update(t);
            zpk.step(t, dt);

            max_val = max_val.max(zpk.get_output(0));
            min_val = min_val.min(zpk.get_output(0));
        }

        let amplitude = (max_val - min_val) / 2.0;
        let attenuation_db = if amplitude > 1e-6 {
            20.0 * (amplitude / 1.0).log10()
        } else {
            -100.0
        };

        println!(
            "  {:6.1} Hz: amplitude = {:.4}, attenuation = {:6.1} dB",
            freq, amplitude, attenuation_db
        );
    }

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
    let mut zpk = ZerosPoleGain::<2>::new(vec![], poles_zpk, 1.0);

    // Create using NumDen
    let mut tf = TransferFunction::<2>::new(&[1.0], &[1.0, 2.0, 5.0]);

    // Apply same step input
    zpk.set_input(0, 1.0);
    tf.set_input(0, 1.0);
    zpk.update(0.0);
    tf.update(0.0);

    println!("Comparing step responses:");
    println!("  Time (s)  |  ZPK      |  NumDen   |  Difference");
    println!("------------|-----------|-----------|-------------");

    let dt = 0.01;

    for i in 0..=100 {
        let t = i as f64 * dt;

        if i % 10 == 0 {
            let diff = (zpk.get_output(0) - tf.get_output(0)).abs();
            println!(
                "  {:7.2}   |  {:.6}  |  {:.6}  |  {:.2e}",
                t,
                zpk.get_output(0),
                tf.get_output(0),
                diff
            );
        }

        zpk.step(t, dt);
        tf.step(t, dt);
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

    let mut zpk = ZerosPoleGain::<2>::new(vec![], poles, gain);

    // Apply step input
    zpk.set_input(0, 1.0);
    zpk.update(0.0);

    println!("Step Response:");
    println!("  Time (s)  |  Output   | Status");
    println!("------------|-----------|------------------");

    let dt = 0.02;
    let mut max_output: f64 = 0.0;
    let mut settling_time = 0.0;
    let mut settled = false;

    for i in 0..=250 {
        let t = i as f64 * dt;
        let output = zpk.get_output(0);
        max_output = max_output.max(output);

        // Check if within 2% of final value
        if !settled && (output - 1.0).abs() < 0.02 {
            // Stay within bounds for several samples to confirm
            let mut stays_settled = true;
            for j in 0..20 {
                zpk.step(t, dt);
                if (zpk.get_output(0) - 1.0).abs() >= 0.02 {
                    stays_settled = false;
                    break;
                }
            }
            if stays_settled {
                settling_time = t;
                settled = true;
            }
            // Reset for continued simulation
            for _ in 0..20 {
                zpk.revert();
            }
        }

        if i % 25 == 0 {
            let status = if settled && t >= settling_time {
                "Settled"
            } else if output > 1.0 {
                "Overshoot"
            } else {
                "Rising"
            };
            println!("  {:7.2}   |  {:.4}    | {}", t, output, status);
        }

        zpk.step(t, dt);
    }

    let overshoot = (max_output - 1.0) / 1.0 * 100.0;

    println!("\nPerformance achieved:");
    println!("  Settling time: {:.2} s (target: {:.2} s)", settling_time, ts);
    println!("  Overshoot: {:.1}% (target: ~5%)", overshoot);
    println!();
}
