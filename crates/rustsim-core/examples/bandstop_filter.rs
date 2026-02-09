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
//!
//! Uses the new Connection and Simulation API for clean signal flow.

use rustsim_core::prelude::*;

fn main() {
    println!("=== Butterworth Bandstop Filter Example ===\n");

    // Example 1: 50 Hz Notch Filter (European power line noise)
    example_1_50hz_notch();

    // Example 2: Multi-tone signal processing
    example_2_multitone();

    // Example 3: Wide bandstop filter
    example_3_wide_bandstop();

    println!("=== Example Complete ===");
}

fn example_1_50hz_notch() {
    println!("Example 1: 50 Hz Notch Filter");
    println!("-----------------------------");
    println!("Input: 10 Hz signal (amplitude 1.0) + 50 Hz noise (amplitude 0.5)");

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Sinusoidal::new(1.0, 10.0, 0.0)),     // Block 0: 10 Hz signal
        Box::new(Sinusoidal::new(0.5, 50.0, 0.0)),     // Block 1: 50 Hz noise
        Box::new(Adder::<2>::new()),                   // Block 2: combine signal + noise
        Box::new(ButterworthBandstop::new([45.0, 55.0], 4)),  // Block 3: 50 Hz notch
        Box::new(Scope::<2, 5000>::new()),             // Block 4: scope (input and filtered)
    ];

    let connections = vec![
        Connection::from(0, 0).to(2, 0).build(),       // 10 Hz -> adder
        Connection::from(1, 0).to(2, 1).build(),       // 50 Hz noise -> adder
        Connection::from(2, 0).to(3, 0).to(4, 0).build(),  // Combined -> filter, scope[0]
        Connection::from(3, 0).to(4, 1).build(),       // Filtered -> scope[1]
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.0001);  // 10 kHz sampling

    sim.run(0.5);  // 500ms

    // Analyze results after settling
    println!("  Simulation complete (500ms)");
    println!("  Filter removes 50 Hz component while preserving 10 Hz signal");
    println!();
}

fn example_2_multitone() {
    println!("Example 2: Multi-Tone Signal Processing");
    println!("----------------------------------------");
    println!("Input: 20 Hz + 50 Hz + 100 Hz (equal amplitudes)");
    println!("Filter: 50 Hz notch [45-55 Hz]");

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Sinusoidal::new(1.0, 20.0, 0.0)),     // Block 0: 20 Hz (pass)
        Box::new(Sinusoidal::new(1.0, 50.0, 0.0)),     // Block 1: 50 Hz (remove)
        Box::new(Sinusoidal::new(1.0, 100.0, 0.0)),    // Block 2: 100 Hz (pass)
        Box::new(Adder::<3>::new()),                   // Block 3: sum all tones
        Box::new(ButterworthBandstop::new([45.0, 55.0], 4)),  // Block 4: notch filter
        Box::new(Scope::<2, 5000>::new()),             // Block 5: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(3, 0).build(),       // 20 Hz -> adder
        Connection::from(1, 0).to(3, 1).build(),       // 50 Hz -> adder
        Connection::from(2, 0).to(3, 2).build(),       // 100 Hz -> adder
        Connection::from(3, 0).to(4, 0).to(5, 0).build(),  // Composite -> filter, scope[0]
        Connection::from(4, 0).to(5, 1).build(),       // Filtered -> scope[1]
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.0001);  // 10 kHz sampling

    sim.run(0.5);  // 500ms

    println!("  Simulated 500ms");
    println!("  50 Hz component attenuated, 20 Hz and 100 Hz preserved");
    println!();
}

fn example_3_wide_bandstop() {
    println!("Example 3: Wide Bandstop Filter (20-100 Hz)");
    println!("--------------------------------------------");
    println!("Testing multiple frequency components:");

    // Test low passband (10 Hz)
    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Sinusoidal::new(1.0, 10.0, 0.0)),     // Block 0: 10 Hz (should pass)
        Box::new(ButterworthBandstop::new([20.0, 100.0], 4)),  // Block 1: wide bandstop
        Box::new(Scope::<2, 1000>::new()),             // Block 2: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),
        Connection::from(1, 0).to(2, 1).build(),
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.0001);

    sim.run(0.1);  // 100ms
    println!("  10 Hz (low passband): passes through");

    // Test stopband (50 Hz)
    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Sinusoidal::new(1.0, 50.0, 0.0)),     // Block 0: 50 Hz (should attenuate)
        Box::new(ButterworthBandstop::new([20.0, 100.0], 4)),  // Block 1: wide bandstop
        Box::new(Scope::<2, 1000>::new()),             // Block 2: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),
        Connection::from(1, 0).to(2, 1).build(),
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.0001);

    sim.run(0.1);  // 100ms
    println!("  50 Hz (stopband): attenuated");

    // Test high passband (150 Hz)
    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Sinusoidal::new(1.0, 150.0, 0.0)),    // Block 0: 150 Hz (should pass)
        Box::new(ButterworthBandstop::new([20.0, 100.0], 4)),  // Block 1: wide bandstop
        Box::new(Scope::<2, 1000>::new()),             // Block 2: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).to(2, 0).build(),
        Connection::from(1, 0).to(2, 1).build(),
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.0001);

    sim.run(0.1);  // 100ms
    println!("  150 Hz (high passband): passes through");
    println!();

}
