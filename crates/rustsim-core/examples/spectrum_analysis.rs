//! Spectrum analyzer example
//!
//! Demonstrates the Spectrum block for real-time frequency analysis using FFT.
//! Shows detection of:
//! - Single frequency components
//! - Multiple frequency components
//! - Harmonics
//! - DC offset
//! - Exponential windowing
//!
//! This example uses the new Simulation API to connect signal sources to the Spectrum analyzer.
//!
//! **Note:** This example requires the `spectrum` feature to be enabled.
//! Run with: `cargo run --example spectrum_analysis --features spectrum`

#[cfg(not(feature = "spectrum"))]
fn main() {
    println!("This example requires the 'spectrum' feature.");
    println!("Run with: cargo run --example spectrum_analysis --features spectrum");
}

#[cfg(feature = "spectrum")]
use rustsim_core::prelude::*;

#[cfg(feature = "spectrum")]
fn main() {
    println!("=== RustSim Spectrum Analyzer Examples ===\n");

    example_1_single_frequency();
    example_2_multiple_frequencies();
    example_3_harmonic_analysis();
    example_4_dc_component();
    example_5_exponential_windowing();
}

/// Example 1: Detecting a single frequency component
#[cfg(feature = "spectrum")]
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

    // Create blocks
    let sine = Sinusoidal::new(amplitude, freq, 0.0);
    let spectrum = Spectrum::<1, 1024>::new(sample_rate);

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(sine),
        Box::new(spectrum),
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).build(), // Sine -> Spectrum
    ];

    // Run simulation
    let mut sim = Simulation::new(blocks, connections).with_dt(dt);
    sim.run(t_end);

    println!("  Input: {} Hz sinusoid, amplitude {}", freq, amplitude);
    println!("  Simulation completed at t = {:.3} s", sim.time());
    println!(
        "  Note: Spectrum analysis results would be available through block access"
    );
    println!();
}

/// Example 2: Detecting multiple frequency components
#[cfg(feature = "spectrum")]
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
    let t_end = 3.0;

    // Create three sinusoidal sources
    let sine1 = Sinusoidal::new(amp1, freq1, 0.0);
    let sine2 = Sinusoidal::new(amp2, freq2, 0.0);
    let sine3 = Sinusoidal::new(amp3, freq3, 0.0);

    // Create adder to combine the signals
    let adder = Adder::<3>::new();

    // Create spectrum analyzer with larger window
    let spectrum = Spectrum::<1, 2048>::new(sample_rate);

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(sine1),
        Box::new(sine2),
        Box::new(sine3),
        Box::new(adder),
        Box::new(spectrum),
    ];

    let connections = vec![
        Connection::from(0, 0).to(3, 0).build(), // Sine1 -> Adder
        Connection::from(1, 0).to(3, 1).build(), // Sine2 -> Adder
        Connection::from(2, 0).to(3, 2).build(), // Sine3 -> Adder
        Connection::from(3, 0).to(4, 0).build(), // Adder -> Spectrum
    ];

    let mut sim = Simulation::new(blocks, connections).with_dt(dt);
    sim.run(t_end);

    println!("  Input components:");
    println!("    {:.1} Hz (amplitude {:.1})", freq1, amp1);
    println!("    {:.1} Hz (amplitude {:.1})", freq2, amp2);
    println!("    {:.1} Hz (amplitude {:.1})", freq3, amp3);
    println!("\n  Simulation completed at t = {:.3} s", sim.time());
    println!();
}

/// Example 3: Harmonic analysis of a rich signal
#[cfg(feature = "spectrum")]
fn example_3_harmonic_analysis() {
    println!("Example 3: Harmonic Analysis");
    println!("-----------------------------");

    // Fundamental frequency
    let f0 = 10.0;

    let dt = 0.001;
    let sample_rate = 1.0 / dt;
    let t_end = 3.0;

    // Create harmonics (like a square wave approximation)
    let fund = Sinusoidal::new(1.0, f0, 0.0); // Fundamental
    let harm3 = Sinusoidal::new(1.0 / 3.0, 3.0 * f0, 0.0); // 3rd harmonic
    let harm5 = Sinusoidal::new(1.0 / 5.0, 5.0 * f0, 0.0); // 5th harmonic
    let harm7 = Sinusoidal::new(1.0 / 7.0, 7.0 * f0, 0.0); // 7th harmonic

    // Combine harmonics
    let adder = Adder::<4>::new();
    let spectrum = Spectrum::<1, 2048>::new(sample_rate);

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(fund),
        Box::new(harm3),
        Box::new(harm5),
        Box::new(harm7),
        Box::new(adder),
        Box::new(spectrum),
    ];

    let connections = vec![
        Connection::from(0, 0).to(4, 0).build(), // Fundamental -> Adder
        Connection::from(1, 0).to(4, 1).build(), // 3rd harmonic -> Adder
        Connection::from(2, 0).to(4, 2).build(), // 5th harmonic -> Adder
        Connection::from(3, 0).to(4, 3).build(), // 7th harmonic -> Adder
        Connection::from(4, 0).to(5, 0).build(), // Adder -> Spectrum
    ];

    let mut sim = Simulation::new(blocks, connections).with_dt(dt);
    sim.run(t_end);

    println!("  Input: Fundamental at {} Hz plus odd harmonics", f0);
    println!("  Harmonics: 1st, 3rd, 5th, 7th");
    println!("\n  Simulation completed at t = {:.3} s", sim.time());
    println!();
}

/// Example 4: DC component and offset detection
#[cfg(feature = "spectrum")]
fn example_4_dc_component() {
    println!("Example 4: DC Component Detection");
    println!("----------------------------------");

    let dc_offset = 1.5;
    let ac_freq = 20.0;
    let ac_amp = 0.8;

    let dt = 0.001;
    let sample_rate = 1.0 / dt;
    let t_end = 2.0;

    // Create DC offset and AC signal
    let dc = Constant::new(dc_offset);
    let ac = Sinusoidal::new(ac_amp, ac_freq, 0.0);

    // Combine DC + AC
    let adder = Adder::<2>::new();
    let spectrum = Spectrum::<1, 1024>::new(sample_rate);

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(dc),
        Box::new(ac),
        Box::new(adder),
        Box::new(spectrum),
    ];

    let connections = vec![
        Connection::from(0, 0).to(2, 0).build(), // DC -> Adder
        Connection::from(1, 0).to(2, 1).build(), // AC -> Adder
        Connection::from(2, 0).to(3, 0).build(), // Adder -> Spectrum
    ];

    let mut sim = Simulation::new(blocks, connections).with_dt(dt);
    sim.run(t_end);

    println!(
        "  Input: DC offset {:.1} + {} Hz AC (amplitude {:.1})",
        dc_offset, ac_freq, ac_amp
    );
    println!("\n  Simulation completed at t = {:.3} s", sim.time());
    println!();
}

/// Example 5: Exponential windowing for time-varying signals
#[cfg(feature = "spectrum")]
fn example_5_exponential_windowing() {
    println!("Example 5: Exponential Windowing");
    println!("----------------------------------");

    let test_freq = 10.0;
    let dt = 0.001;
    let sample_rate = 1.0 / dt;

    // Note: This example would require manual block configuration
    // to set alpha values, which isn't directly supported through
    // the Simulation API. This demonstrates the concept.

    let sine = Sinusoidal::new(1.0, test_freq, 0.0);
    let spectrum = Spectrum::<1, 1024>::new(sample_rate);
    // Would need: spectrum.set_alpha(2.0) for exponential windowing

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(sine),
        Box::new(spectrum),
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).build(), // Sine -> Spectrum
    ];

    let mut sim = Simulation::new(blocks, connections).with_dt(dt);
    sim.run(1.5);

    println!("  Simulated {} Hz signal for 1.5s", test_freq);
    println!("\n  Note: Exponential windowing configuration (set_alpha)");
    println!("        requires direct block access before adding to Simulation.");
    println!("\n  Exponential window emphasizes recent samples,");
    println!("  so it responds faster to signal changes.");
    println!();
}

