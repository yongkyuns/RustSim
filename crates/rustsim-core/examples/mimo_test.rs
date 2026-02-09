//! Test MIMO blocks with PathSim-style API

use rustsim_core::prelude::*;
use rustsim_core::connection::{ConnectionDef, BlockExt};
use rustsim_core::simulation::SimulationBuilder;

fn main() {
    println!("MIMO Block Test with PathSim-Style API");
    println!("=======================================\n");

    // Test 1: 2-input Adder
    println!("Test 1: Two signals -> Adder");
    println!("  signal1 (3.0) + signal2 (5.0) = 8.0\n");
    
    let signal1 = Constant::new(3.0);
    let signal2 = Constant::new(5.0);
    let adder = Adder::<2>::new();
    
    let connections = vec![
        ConnectionDef::new(&signal1, adder.port(0)),  // signal1 -> adder input 0
        ConnectionDef::new(&signal2, adder.port(1)),  // signal2 -> adder input 1
    ];
    
    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(signal1),
        Box::new(signal2),
        Box::new(adder),
    ];
    
    let mut sim = SimulationBuilder::build_from_defs(blocks, connections);
    sim.step();
    
    println!("  Result: {} (expected 8.0)\n", sim.get_output(2, 0));

    // Test 2: 3-input Adder with subtraction
    println!("Test 2: Weighted Adder (a + b - c)");
    println!("  5.0 + 3.0 - 2.0 = 6.0\n");
    
    let a = Constant::new(5.0);
    let b = Constant::new(3.0);
    let c = Constant::new(2.0);
    let weighted_adder = Adder::<3>::with_weights([1.0, 1.0, -1.0]);
    
    let connections = vec![
        ConnectionDef::new(&a, weighted_adder.port(0)),
        ConnectionDef::new(&b, weighted_adder.port(1)),
        ConnectionDef::new(&c, weighted_adder.port(2)),
    ];
    
    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(a),
        Box::new(b),
        Box::new(c),
        Box::new(weighted_adder),
    ];
    
    let mut sim = SimulationBuilder::build_from_defs(blocks, connections);
    sim.step();
    
    println!("  Result: {} (expected 6.0)\n", sim.get_output(3, 0));

    // Test 3: Multi-channel Scope
    println!("Test 3: Two signals -> 2-channel Scope");
    
    let ch1 = Constant::new(1.5);
    let ch2 = Constant::new(2.5);
    let scope = Scope::<2, 100>::new();
    
    let connections = vec![
        ConnectionDef::new(&ch1, scope.port(0)),  // ch1 -> scope channel 0
        ConnectionDef::new(&ch2, scope.port(1)),  // ch2 -> scope channel 1
    ];
    
    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(ch1),
        Box::new(ch2),
        Box::new(scope),
    ];
    
    let mut sim = SimulationBuilder::build_from_defs(blocks, connections);
    sim.step();
    
    println!("  Scope channel 0: {} (expected 1.5)", sim.get_output(2, 0));
    println!("  Scope channel 1: {} (expected 2.5)", sim.get_output(2, 1));
    
    println!("\nAll MIMO tests passed!");
}
