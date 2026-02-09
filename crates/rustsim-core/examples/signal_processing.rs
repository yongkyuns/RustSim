//! Signal processing example using Connection and Simulation API
//!
//! Demonstrates source blocks, mathematical operations, and signal composition
//!
//! Creates multiple signal processing chains:
//! 1. y = 2*sin(2*pi*t) + 1
//! 2. Complex waveforms with multiple sources
//! 3. Filtering and mixing

use rustsim_core::prelude::*;
use std::f64::consts::PI;

fn main() {
    println!("Signal Processing - Connection API Demo");
    println!("========================================\n");

    // Example 1: y = 2*sin(2*pi*t) + 1
    example_1_amplified_sine_with_offset();

    // Example 2: Multi-frequency composite signal
    example_2_composite_signal();

    // Example 3: Signal mixing and filtering
    example_3_mixing_and_filtering();

    println!("\nDemo complete!");
}

fn example_1_amplified_sine_with_offset() {
    println!("Example 1: Amplified Sine with DC Offset");
    println!("Signal: y = 2*sin(2*pi*t) + 1");
    println!("-----------------------------------------");

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Sinusoidal::new(1.0, 1.0, 0.0)),  // Block 0: sin(2*pi*1*t), amp=1, freq=1Hz
        Box::new(Amplifier::new(2.0)),              // Block 1: gain = 2
        Box::new(Constant::new(1.0)),               // Block 2: DC offset = 1
        Box::new(Adder::<2>::new()),                // Block 3: add amplified sine + DC
        Box::new(Scope::<1, 100>::new()),           // Block 4: record output
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).build(),    // Sine -> Amplifier
        Connection::from(1, 0).to(3, 0).build(),    // Amplifier -> Adder[0]
        Connection::from(2, 0).to(3, 1).build(),    // Constant -> Adder[1]
        Connection::from(3, 0).to(4, 0).build(),    // Adder -> Scope
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.01);

    println!("\n{:>10} {:>12} {:>12} {:>12}",
        "Time", "Output", "Exact", "Error");
    println!("{:-<10} {:-<12} {:-<12} {:-<12}", "", "", "", "");

    // Simulate for one period
    let period = 1.0;
    let steps = (period / 0.01) as usize;

    for i in 0..=steps {
        if i % (steps / 20) == 0 {
            let t = sim.time();
            let y = sim.get_output(3, 0);  // Output from adder
            let exact = 2.0 * (2.0 * PI * t).sin() + 1.0;
            let error = (y - exact).abs();
            println!("{:10.4} {:12.6} {:12.6} {:12.2e}", t, y, exact, error);
        }
        sim.step();
    }
    println!();
}

fn example_2_composite_signal() {
    println!("Example 2: Multi-Frequency Composite Signal");
    println!("Signal: y = sin(2πt) + 0.5*sin(4πt) + 0.25*sin(8πt)");
    println!("----------------------------------------------------");

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Sinusoidal::new(1.0, 1.0, 0.0)),   // Block 0: fundamental (1 Hz)
        Box::new(Sinusoidal::new(1.0, 2.0, 0.0)),   // Block 1: 2nd harmonic (2 Hz)
        Box::new(Amplifier::new(0.5)),              // Block 2: scale 2nd harmonic
        Box::new(Sinusoidal::new(1.0, 4.0, 0.0)),   // Block 3: 4th harmonic (4 Hz)
        Box::new(Amplifier::new(0.25)),             // Block 4: scale 4th harmonic
        Box::new(Adder::<3>::new()),                // Block 5: sum all harmonics
        Box::new(Scope::<4, 200>::new()),           // Block 6: scope (fundamental + harmonics + composite)
    ];

    let connections = vec![
        Connection::from(0, 0).to(5, 0).to(6, 0).build(),  // Fundamental -> adder, scope[0]
        Connection::from(1, 0).to(2, 0).build(),           // 2nd harmonic -> amplifier
        Connection::from(2, 0).to(5, 1).to(6, 1).build(),  // Scaled 2nd -> adder, scope[1]
        Connection::from(3, 0).to(4, 0).build(),           // 4th harmonic -> amplifier
        Connection::from(4, 0).to(5, 2).to(6, 2).build(),  // Scaled 4th -> adder, scope[2]
        Connection::from(5, 0).to(6, 3).build(),           // Composite -> scope[3]
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.001);

    sim.run(1.0);  // Run for one second

    println!("Simulated 1 second of composite waveform");
    println!("Peak amplitude: ~{:.3}", sim.get_output(5, 0));
    println!("Contains harmonics at 1 Hz, 2 Hz, and 4 Hz\n");
}

fn example_3_mixing_and_filtering() {
    println!("Example 3: Signal Mixing and Low-Pass Filtering");
    println!("------------------------------------------------");
    println!("Two signals mixed and filtered:");
    println!("  Signal 1: 5 Hz sine wave");
    println!("  Signal 2: 50 Hz square wave (high frequency)");
    println!("  Low-pass filter smooths the result");

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Sinusoidal::new(1.0, 5.0, 0.0)),     // Block 0: 5 Hz sine
        Box::new(SquareWave::new(0.5, 50.0, 0.0)),    // Block 1: 50 Hz square
        Box::new(Amplifier::new(0.3)),                 // Block 2: attenuate square wave
        Box::new(Adder::<2>::new()),                   // Block 3: mix signals
        Box::new(LowpassRC::new(0.01)),               // Block 4: RC low-pass (tau=0.01, fc ≈ 16 Hz)
        Box::new(Scope::<3, 1000>::new()),            // Block 5: scope
    ];

    let connections = vec![
        Connection::from(0, 0).to(3, 0).to(5, 0).build(),  // 5 Hz sine -> mixer, scope[0]
        Connection::from(1, 0).to(2, 0).build(),           // Square -> attenuator
        Connection::from(2, 0).to(3, 1).build(),           // Attenuated square -> mixer
        Connection::from(3, 0).to(4, 0).to(5, 1).build(),  // Mixed -> filter, scope[1]
        Connection::from(4, 0).to(5, 2).build(),           // Filtered -> scope[2]
    ];

    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.0001);  // High sample rate for 50 Hz signal

    sim.run(0.2);  // 200ms

    println!("\nResults after 200ms:");
    println!("  Raw mixed signal: {:.3}", sim.get_output(3, 0));
    println!("  Filtered output: {:.3}", sim.get_output(4, 0));
    println!("  Filter removes high-frequency square wave component\n");
}
