//! MIMO (Multi-Input Multi-Output) block example
//!
//! This example demonstrates connecting multiple signals to a MIMO block
//! using the PathSim-style API with the `simulation!` macro.

use rustsim_core::blocks::{Adder, Amplifier, Constant};
use rustsim_core::connection::{BlockExt, ConnectionDef};
use rustsim_core::simulation;

fn main() {
    println!("=== MIMO Block Example ===\n");

    // Example 1: Simple 2-input adder
    println!("Example 1: Two constants -> Adder");
    {
        let signal1 = Constant::new(3.0);
        let signal2 = Constant::new(5.0);
        let adder = Adder::<2>::new();

        // Using the simulation! macro - order doesn't matter for blocks
        let mut sim = simulation! {
            connections: [
                ConnectionDef::new(&signal1, adder.port(0)),
                ConnectionDef::new(&signal2, adder.port(1)),
            ],
            blocks: [signal1, signal2, adder]
        };

        sim.step();

        println!("  signal1 = 3.0");
        println!("  signal2 = 5.0");
        println!("  adder output = {} (expected 8.0)", sim.get_output(2, 0));
        println!();
    }

    // Example 2: Weighted adder (subtractor)
    println!("Example 2: Weighted adder as subtractor");
    {
        let a = Constant::new(10.0);
        let b = Constant::new(4.0);
        let sub = Adder::<2>::subtractor(); // weights: [1.0, -1.0]

        let mut sim = simulation! {
            connections: [
                ConnectionDef::new(&a, sub.port(0)),
                ConnectionDef::new(&b, sub.port(1)),
            ],
            blocks: [a, b, sub]
        };

        sim.step();

        println!("  a = 10.0");
        println!("  b = 4.0");
        println!("  a - b = {} (expected 6.0)", sim.get_output(2, 0));
        println!();
    }

    // Example 3: Fan-out to multiple blocks
    println!("Example 3: Fan-out (one source to multiple targets)");
    {
        let source = Constant::new(5.0);
        let gain2x = Amplifier::new(2.0);
        let gain3x = Amplifier::new(3.0);

        let mut sim = simulation! {
            connections: [
                // One source feeding two different amplifiers
                ConnectionDef::new(&source, &gain2x).to(&gain3x),
            ],
            blocks: [source, gain2x, gain3x]
        };

        sim.step();

        println!("  source = 5.0");
        println!("  gain2x output = {} (expected 10.0)", sim.get_output(1, 0));
        println!("  gain3x output = {} (expected 15.0)", sim.get_output(2, 0));
        println!();
    }

    // Example 4: Complex MIMO network
    println!("Example 4: Complex network (3 inputs to adder, then amplified)");
    {
        let s1 = Constant::new(1.0);
        let s2 = Constant::new(2.0);
        let s3 = Constant::new(3.0);
        let adder = Adder::<3>::new();
        let gain = Amplifier::new(2.0);

        let mut sim = simulation! {
            connections: [
                ConnectionDef::new(&s1, adder.port(0)),
                ConnectionDef::new(&s2, adder.port(1)),
                ConnectionDef::new(&s3, adder.port(2)),
                ConnectionDef::new(&adder, &gain),
            ],
            blocks: [s1, s2, s3, adder, gain]
        };

        sim.step();

        println!("  s1 = 1.0, s2 = 2.0, s3 = 3.0");
        println!("  adder output = {} (expected 6.0)", sim.get_output(3, 0));
        println!("  gain output = {} (expected 12.0)", sim.get_output(4, 0));
    }
}
