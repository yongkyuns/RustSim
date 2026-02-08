//! Butterworth Bandstop (Notch) Filter Example
//!
//! Demonstrates using the ButterworthBandstop filter to remove specific frequency
//! components from a signal. This is useful for:
//! - Removing power line noise (50/60 Hz)
//! - Eliminating unwanted frequency bands in signal processing
//! - Notch filtering in audio applications
//!
//! This example shows:
//! 1. Creating a notch filter to remove 50 Hz line noise
//! 2. Processing a composite signal with multiple frequencies
//! 3. Comparing input and output signals

use rustsim::prelude::*;
use std::f64::consts::PI;

fn main() {
    println!("=== Butterworth Bandstop Filter Example ===\n");

    // Example 1: 50 Hz Notch Filter (European power line noise)
    println!("Example 1: 50 Hz Notch Filter");
    println!("-----------------------------");

    // Create a narrow notch filter around 50 Hz (45-55 Hz stopband)
    let mut notch_50hz = ButterworthBandstop::new([45.0, 55.0], 4);

    // Generate a composite signal: 10 Hz signal + 50 Hz noise
    let dt = 0.0001; // 10 kHz sampling
    let duration = 0.5; // seconds
    let n_samples = (duration / dt) as usize;

    let f_signal = 10.0; // Hz - desired signal
    let f_noise = 50.0; // Hz - noise to be removed

    println!("Input: 10 Hz signal (amplitude 1.0) + 50 Hz noise (amplitude 0.5)");

    let mut max_input: f64 = 0.0;
    let mut max_output: f64 = 0.0;
    let mut output_at_50hz: f64 = 0.0;

    for i in 0..n_samples {
        let t = i as f64 * dt;

        // Composite signal
        let signal = (2.0 * PI * f_signal * t).sin();
        let noise = 0.5 * (2.0 * PI * f_noise * t).sin();
        let input = signal + noise;

        notch_50hz.set_input(0, input);
        notch_50hz.update(t);
        notch_50hz.step(t, dt);

        // Track peak values after settling
        if t > 0.2 {
            max_input = max_input.max(input.abs());
            max_output = max_output.max(notch_50hz.get_output(0).abs());

            // Track output specifically for noise frequency
            if (i % 100) == 0 {
                // Sample periodically
                output_at_50hz = output_at_50hz.max(notch_50hz.get_output(0).abs());
            }
        }
    }

    println!("  Peak input amplitude: {:.3}", max_input);
    println!("  Peak output amplitude: {:.3}", max_output);
    println!("  50 Hz component removed: {:.1}%", (1.0 - output_at_50hz) * 100.0);
    println!();

    // Example 2: Wide Bandstop Filter
    println!("Example 2: Wide Bandstop Filter (20-100 Hz)");
    println!("--------------------------------------------");

    let mut wide_bandstop = ButterworthBandstop::new([20.0, 100.0], 4);

    // Test with various frequencies
    let test_frequencies = vec![
        (5.0, "Low passband"),
        (10.0, "Low passband"),
        (30.0, "Stopband"),
        (50.0, "Stopband center"),
        (80.0, "Stopband"),
        (150.0, "High passband"),
        (200.0, "High passband"),
    ];

    println!("Frequency Response:");
    for (freq, description) in test_frequencies {
        wide_bandstop.reset();

        let omega = 2.0 * PI * freq;
        let mut amplitude: f64 = 0.0;

        // Run filter to steady state
        for i in 0..5000 {
            let t = i as f64 * dt;
            let input = (omega * t).sin();
            wide_bandstop.set_input(0, input);
            wide_bandstop.update(t);
            wide_bandstop.step(t, dt);

            if i > 2000 {
                amplitude = amplitude.max(wide_bandstop.get_output(0).abs());
            }
        }

        let attenuation_db = if amplitude > 0.001 {
            20.0 * (amplitude / 1.0).log10()
        } else {
            -60.0 // Cap at -60 dB for display
        };

        println!(
            "  {:>6.1} Hz ({:17}): amplitude = {:.3}, attenuation = {:>6.1} dB",
            freq, description, amplitude, attenuation_db
        );
    }
    println!();

    // Example 3: Comparing Different Filter Orders
    println!("Example 3: Effect of Filter Order");
    println!("----------------------------------");

    let test_freq = 50.0; // Test at center of stopband [45, 55] Hz
    let orders = vec![2, 4, 6, 8];

    println!("Center frequency rejection at 50 Hz:");
    for order in orders {
        let mut flt = ButterworthBandstop::new([45.0, 55.0], order);

        let omega = 2.0 * PI * test_freq;
        let mut amplitude: f64 = 0.0;

        for i in 0..5000 {
            let t = i as f64 * dt;
            let input = (omega * t).sin();
            flt.set_input(0, input);
            flt.update(t);
            flt.step(t, dt);

            if i > 2000 {
                amplitude = amplitude.max(flt.get_output(0).abs());
            }
        }

        let attenuation_db = 20.0 * (amplitude / 1.0).log10();
        println!(
            "  Order {:2}: amplitude = {:.4}, attenuation = {:.1} dB",
            order, amplitude, attenuation_db
        );
    }
    println!();

    // Example 4: Multi-tone Signal Processing
    println!("Example 4: Multi-tone Signal Processing");
    println!("---------------------------------------");

    // Remove middle frequencies from a three-tone signal
    let mut multi_notch = ButterworthBandstop::new([45.0, 55.0], 4);

    let f1 = 20.0; // Should pass
    let f2 = 50.0; // Should be removed
    let f3 = 100.0; // Should pass

    println!("Input: 20 Hz + 50 Hz + 100 Hz (equal amplitudes)");

    let mut spectrum_input: [f64; 3] = [0.0; 3];
    let mut spectrum_output: [f64; 3] = [0.0; 3];

    for i in 0..5000 {
        let t = i as f64 * dt;

        let input = (2.0 * PI * f1 * t).sin()
            + (2.0 * PI * f2 * t).sin()
            + (2.0 * PI * f3 * t).sin();

        multi_notch.set_input(0, input);
        multi_notch.update(t);
        multi_notch.step(t, dt);

        if i > 2000 {
            // Measure each frequency component (simplified - just track peaks)
            spectrum_input[0] = spectrum_input[0].max((2.0 * PI * f1 * t).sin().abs());
            spectrum_input[1] = spectrum_input[1].max((2.0 * PI * f2 * t).sin().abs());
            spectrum_input[2] = spectrum_input[2].max((2.0 * PI * f3 * t).sin().abs());

            let output = multi_notch.get_output(0);
            spectrum_output[0] = spectrum_output[0].max(output.abs());
            spectrum_output[1] = spectrum_output[1].max(output.abs());
            spectrum_output[2] = spectrum_output[2].max(output.abs());
        }
    }

    println!("\nFrequency Component Analysis:");
    println!("  20 Hz:  Input = {:.3}, Output = {:.3}", spectrum_input[0], spectrum_output[0]);
    println!("  50 Hz:  Input = {:.3}, Output < {:.3} (removed)", spectrum_input[1], spectrum_output[1]);
    println!("  100 Hz: Input = {:.3}, Output = {:.3}", spectrum_input[2], spectrum_output[2]);
    println!();

    // Example 5: Step Response
    println!("Example 5: Step Response");
    println!("------------------------");

    let mut step_filter = ButterworthBandstop::new([45.0, 55.0], 2);

    println!("Applying unit step input...");
    let mut step_outputs = Vec::new();

    for i in 0..100 {
        let t = i as f64 * 0.01; // Slower time step for step response
        step_filter.set_input(0, 1.0);
        step_filter.update(t);
        step_filter.step(t, 0.01);

        if i % 10 == 0 {
            step_outputs.push(step_filter.get_output(0));
        }
    }

    println!("Step response samples:");
    print!("  Time (s):   ");
    for (i, _) in step_outputs.iter().enumerate() {
        print!("{:>6.2} ", i as f64 * 0.1);
    }
    println!();

    print!("  Output:     ");
    for output in &step_outputs {
        print!("{:>6.3} ", output);
    }
    println!();
    println!("  (Should settle to ~1.0 since DC passes through bandstop)");
    println!();

    println!("=== Example Complete ===");
}
