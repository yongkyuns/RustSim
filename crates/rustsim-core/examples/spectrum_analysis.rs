//! Spectrum analyzer example
//!
//! Demonstrates the Spectrum block for real-time frequency analysis using FFT.
//! Shows detection of:
//! - Single frequency components
//! - Multiple frequency components
//! - Harmonics
//! - DC offset
//! - Exponential windowing

use rustsim::blocks::Spectrum;
use rustsim::block::Block;
use std::f64::consts::PI;

fn main() {
    println!("=== RustSim Spectrum Analyzer Examples ===\n");

    example_1_single_frequency();
    example_2_multiple_frequencies();
    example_3_harmonic_analysis();
    example_4_dc_component();
    example_5_exponential_windowing();
}

/// Example 1: Detecting a single frequency component
fn example_1_single_frequency() {
    println!("Example 1: Single Frequency Detection");
    println!("--------------------------------------");

    // Signal parameters
    let freq = 10.0; // Hz
    let amplitude = 2.0;

    // Simulation parameters
    let dt = 0.001; // 1 ms timestep
    let sample_rate = 1.0 / dt; // 1000 Hz
    let t_end = 2.0; // 2 seconds

    // Create spectrum analyzer: 1 channel, 1024-sample window
    let mut spectrum = Spectrum::<1, 1024>::new(sample_rate);

    // Run simulation
    let steps = (t_end / dt) as usize;
    let mut t = 0.0;

    for _ in 0..steps {
        let signal = amplitude * (2.0 * PI * freq * t).sin();
        spectrum.set_input(0, signal);
        spectrum.update(t);
        spectrum.step(t, dt);
        t += dt;
    }

    // Get results
    let freqs = spectrum.frequencies();
    let mag = spectrum.magnitude();

    // Find peak
    let (peak_idx, peak_mag) = find_peak(&mag[0]);
    let peak_freq = freqs[peak_idx];

    println!("  Input: {} Hz sinusoid, amplitude {}", freq, amplitude);
    println!("  Detected peak: {:.2} Hz with magnitude {:.3}", peak_freq, peak_mag);
    println!("  Frequency error: {:.3} Hz", (peak_freq - freq).abs());
    println!("  Frequency resolution: {:.3} Hz", spectrum.frequency_resolution());
    println!();
}

/// Example 2: Detecting multiple frequency components
fn example_2_multiple_frequencies() {
    println!("Example 2: Multiple Frequency Components");
    println!("-----------------------------------------");

    // Signal with 3 frequency components
    let freq1 = 5.0;
    let freq2 = 15.0;
    let freq3 = 25.0;
    let amp1 = 1.0;
    let amp2 = 0.5;
    let amp3 = 0.3;

    let dt = 0.001;
    let sample_rate = 1.0 / dt;
    let mut spectrum = Spectrum::<1, 2048>::new(sample_rate); // Larger window for better resolution

    let t_end = 3.0;
    let steps = (t_end / dt) as usize;
    let mut t = 0.0;

    for _ in 0..steps {
        let signal = amp1 * (2.0 * PI * freq1 * t).sin()
            + amp2 * (2.0 * PI * freq2 * t).sin()
            + amp3 * (2.0 * PI * freq3 * t).sin();

        spectrum.set_input(0, signal);
        spectrum.update(t);
        spectrum.step(t, dt);
        t += dt;
    }

    let freqs = spectrum.frequencies();
    let mag = spectrum.magnitude();

    // Find peaks
    let peaks = find_peaks(&mag[0], &freqs, 0.1);

    println!("  Input components:");
    println!("    {:.1} Hz (amplitude {:.1})", freq1, amp1);
    println!("    {:.1} Hz (amplitude {:.1})", freq2, amp2);
    println!("    {:.1} Hz (amplitude {:.1})", freq3, amp3);
    println!("\n  Detected peaks:");
    for (f, m) in peaks.iter().take(5) {
        println!("    {:.2} Hz (magnitude {:.3})", f, m);
    }
    println!();
}

/// Example 3: Harmonic analysis of a rich signal
fn example_3_harmonic_analysis() {
    println!("Example 3: Harmonic Analysis");
    println!("-----------------------------");

    // Fundamental frequency
    let f0 = 10.0;

    let dt = 0.001;
    let sample_rate = 1.0 / dt;
    let mut spectrum = Spectrum::<1, 2048>::new(sample_rate);

    let t_end = 3.0;
    let steps = (t_end / dt) as usize;
    let mut t = 0.0;

    // Create a signal with fundamental + harmonics (like a square wave)
    for _ in 0..steps {
        let w = 2.0 * PI * f0;
        let signal = (w * t).sin()
            + (1.0 / 3.0) * (3.0 * w * t).sin()  // 3rd harmonic
            + (1.0 / 5.0) * (5.0 * w * t).sin()  // 5th harmonic
            + (1.0 / 7.0) * (7.0 * w * t).sin(); // 7th harmonic

        spectrum.set_input(0, signal);
        spectrum.update(t);
        spectrum.step(t, dt);
        t += dt;
    }

    let freqs = spectrum.frequencies();
    let mag = spectrum.magnitude();

    let peaks = find_peaks(&mag[0], &freqs, 0.05);

    println!("  Input: Fundamental at {} Hz plus odd harmonics", f0);
    println!("\n  Detected harmonics:");
    for (f, m) in peaks.iter().take(6) {
        let harmonic = f / f0;
        println!("    {:.1} Hz (magnitude {:.3}) - {:.1}x fundamental", f, m, harmonic);
    }
    println!();
}

/// Example 4: DC component and offset detection
fn example_4_dc_component() {
    println!("Example 4: DC Component Detection");
    println!("----------------------------------");

    let dc_offset = 1.5;
    let ac_freq = 20.0;
    let ac_amp = 0.8;

    let dt = 0.001;
    let sample_rate = 1.0 / dt;
    let mut spectrum = Spectrum::<1, 1024>::new(sample_rate);

    let t_end = 2.0;
    let steps = (t_end / dt) as usize;
    let mut t = 0.0;

    for _ in 0..steps {
        let signal = dc_offset + ac_amp * (2.0 * PI * ac_freq * t).sin();

        spectrum.set_input(0, signal);
        spectrum.update(t);
        spectrum.step(t, dt);
        t += dt;
    }

    let freqs = spectrum.frequencies();
    let mag = spectrum.magnitude();

    println!("  Input: DC offset {:.1} + {} Hz AC (amplitude {:.1})", dc_offset, ac_freq, ac_amp);
    println!("\n  Spectrum:");
    println!("    DC (0 Hz): {:.3}", mag[0][0]);

    let (peak_idx, peak_mag) = find_peak_skip_dc(&mag[0]);
    println!("    AC peak: {:.2} Hz (magnitude {:.3})", freqs[peak_idx], peak_mag);
    println!();
}

/// Example 5: Exponential windowing for time-varying signals
fn example_5_exponential_windowing() {
    println!("Example 5: Exponential Windowing");
    println!("----------------------------------");

    // Compare uniform window vs exponential window

    let test_freq = 10.0;
    let dt = 0.001;
    let sample_rate = 1.0 / dt;

    let mut spectrum_uniform = Spectrum::<1, 1024>::new(sample_rate);
    spectrum_uniform.set_alpha(0.0); // No windowing

    let mut spectrum_exponential = Spectrum::<1, 1024>::new(sample_rate);
    spectrum_exponential.set_alpha(2.0); // Moderate exponential windowing

    // First phase: sinusoid at 10 Hz
    let t_end1 = 1.5;
    let steps1 = (t_end1 / dt) as usize;
    let mut t = 0.0;

    for _ in 0..steps1 {
        let signal = (2.0 * PI * test_freq * t).sin();

        spectrum_uniform.set_input(0, signal);
        spectrum_uniform.update(t);
        spectrum_uniform.step(t, dt);

        spectrum_exponential.set_input(0, signal);
        spectrum_exponential.update(t);
        spectrum_exponential.step(t, dt);

        t += dt;
    }

    println!("  After 1.5s at {} Hz:", test_freq);

    let mag_uniform_1 = spectrum_uniform.magnitude();
    let mag_exp_1 = spectrum_exponential.magnitude();
    let freqs = spectrum_uniform.frequencies();

    let (peak_uniform_1, mag_uniform_1_peak) = find_peak(&mag_uniform_1[0]);
    let (peak_exp_1, mag_exp_1_peak) = find_peak(&mag_exp_1[0]);

    println!("    Uniform window peak: {:.1} Hz (mag {:.3})", freqs[peak_uniform_1], mag_uniform_1_peak);
    println!("    Exponential window peak: {:.1} Hz (mag {:.3})", freqs[peak_exp_1], mag_exp_1_peak);

    // Second phase: zero signal
    let t_end2 = 0.5;
    let steps2 = (t_end2 / dt) as usize;

    for _ in 0..steps2 {
        let signal = 0.0;

        spectrum_uniform.set_input(0, signal);
        spectrum_uniform.update(t);
        spectrum_uniform.step(t, dt);

        spectrum_exponential.set_input(0, signal);
        spectrum_exponential.update(t);
        spectrum_exponential.step(t, dt);

        t += dt;
    }

    println!("\n  After 0.5s of zero signal:");

    let mag_uniform_2 = spectrum_uniform.magnitude();
    let mag_exp_2 = spectrum_exponential.magnitude();

    let (peak_uniform_2, mag_uniform_2_peak) = find_peak(&mag_uniform_2[0]);
    let (peak_exp_2, mag_exp_2_peak) = find_peak(&mag_exp_2[0]);

    println!("    Uniform window peak: {:.1} Hz (mag {:.3})", freqs[peak_uniform_2], mag_uniform_2_peak);
    println!("    Exponential window peak: {:.1} Hz (mag {:.3})", freqs[peak_exp_2], mag_exp_2_peak);

    println!("\n  Exponential window emphasizes recent samples,");
    println!("  so it responds faster to signal changes.");
    println!();
}

// Helper functions

fn find_peak(data: &[f64]) -> (usize, f64) {
    let mut peak_idx = 0;
    let mut peak_val = 0.0;

    for (i, &val) in data.iter().enumerate() {
        if val > peak_val {
            peak_val = val;
            peak_idx = i;
        }
    }

    (peak_idx, peak_val)
}

fn find_peak_skip_dc(data: &[f64]) -> (usize, f64) {
    let mut peak_idx = 1; // Skip DC
    let mut peak_val = 0.0;

    for (i, &val) in data.iter().enumerate().skip(1) {
        if val > peak_val {
            peak_val = val;
            peak_idx = i;
        }
    }

    (peak_idx, peak_val)
}

fn find_peaks(data: &[f64], freqs: &[f64], threshold: f64) -> Vec<(f64, f64)> {
    let mut peaks = vec![];

    for i in 1..data.len() - 1 {
        // Local maximum
        if data[i] > data[i - 1] && data[i] > data[i + 1] && data[i] > threshold {
            peaks.push((freqs[i], data[i]));
        }
    }

    // Sort by magnitude (descending)
    peaks.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    peaks
}
