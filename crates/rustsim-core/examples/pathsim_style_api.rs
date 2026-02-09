//! Example demonstrating the PathSim-style Connection API
//!
//! This example shows the new readable, type-safe connection API that uses
//! block references instead of numeric IDs.

use rustsim_core::anyblock::AnyBlock;
use rustsim_core::blocks::{Amplifier, Constant, Integrator};
use rustsim_core::{BlockExt, ConnectionDef, SimulationBuilder};

fn main() {
    println!("PathSim-Style Connection API Demo");
    println!("===================================");
    println!();

    // === Example 1: Using build_from_defs (recommended) ===
    example_build_from_defs();
    println!();

    // === Example 2: Using explicit ports ===
    example_explicit_ports();
    println!();

    // === Example 3: Using fan-out ===
    example_fanout();
}

fn example_build_from_defs() {
    println!("Example 1: build_from_defs (Recommended Pattern)");
    println!("-------------------------------------------------");

    // 1. Create blocks
    let constant = Constant::new(1.0);
    let gain = Amplifier::new(2.0);
    let integrator = Integrator::new(0.0);

    // 2. Define connections using references (before moving blocks)
    let connections = vec![
        ConnectionDef::new(&constant, &gain),
        ConnectionDef::new(&gain, &integrator),
    ];

    // 3. Collect blocks into Vec (blocks moved here)
    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(constant),
        Box::new(gain),
        Box::new(integrator),
    ];

    // 4. Build simulation
    let mut sim = SimulationBuilder::build_from_defs(blocks, connections).with_dt(0.01);

    println!("Running: Constant(1.0) -> Gain(2.0) -> Integrator");
    println!("Expected: integrator output = 2.0 * t");
    println!();

    // Run for 5 seconds
    sim.run(5.0);

    let final_value = sim.get_output(2, 0); // Integrator is block 2
    println!("Final value: {:.6} (expected 10.0)", final_value);
    println!("Error: {:.2e}", (final_value - 10.0).abs());
}

fn example_explicit_ports() {
    println!("Example 2: Using Explicit Ports");
    println!("--------------------------------");

    // Create blocks
    let constant = Constant::new(5.0);
    let gain = Amplifier::new(3.0);

    // Use explicit port references (port 0 in this case)
    let connections = vec![
        ConnectionDef::new(constant.port(0), gain.port(0)),
    ];

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(constant),
        Box::new(gain),
    ];

    let mut sim = SimulationBuilder::build_from_defs(blocks, connections).with_dt(0.01);
    sim.step();

    println!("Input: 5.0, Gain: 3.0");
    println!("Output: {:.1} (expected 15.0)", sim.get_output(1, 0));
}

fn example_fanout() {
    println!("Example 3: Fan-out (One Source, Multiple Targets)");
    println!("--------------------------------------------------");

    // Create blocks
    let constant = Constant::new(2.0);
    let gain1 = Amplifier::new(2.0);
    let gain2 = Amplifier::new(3.0);
    let gain3 = Amplifier::new(4.0);

    // One source feeding multiple targets
    let connections = vec![
        ConnectionDef::new(&constant, &gain1)
            .to(&gain2)
            .to(&gain3),
    ];

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(constant),
        Box::new(gain1),
        Box::new(gain2),
        Box::new(gain3),
    ];

    let mut sim = SimulationBuilder::build_from_defs(blocks, connections).with_dt(0.01);
    sim.step();

    println!("Input: 2.0");
    println!("Gain1 (x2): {:.1} (expected 4.0)", sim.get_output(1, 0));
    println!("Gain2 (x3): {:.1} (expected 6.0)", sim.get_output(2, 0));
    println!("Gain3 (x4): {:.1} (expected 8.0)", sim.get_output(3, 0));
}
