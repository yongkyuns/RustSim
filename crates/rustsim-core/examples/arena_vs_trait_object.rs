//! Performance comparison: Arena-based simulation vs Trait object simulation
//!
//! This example demonstrates the performance difference between:
//! 1. `Simulation` with `Box<dyn AnyBlock>` (trait objects, dynamic dispatch)
//! 2. `ArenaSimulation` with `BlockKind` enum (arena allocation, match dispatch)
//!
//! Run with: cargo run --example arena_vs_trait_object --release

use rustsim_core::anyblock::AnyBlock;
use rustsim_core::arena_simulation::ArenaSimulation;
use rustsim_core::blocks::{Adder, Amplifier, Constant, Integrator};
use rustsim_core::connection::Connection;
use rustsim_core::simulation::Simulation;
use std::time::Instant;

const NUM_STEPS: usize = 100_000;
const DT: f64 = 0.001;

fn main() {
    println!("=== Simulation Performance Comparison ===\n");
    println!("Steps: {}", NUM_STEPS);
    println!("dt: {}\n", DT);

    // Build and run trait-object simulation
    let trait_obj_time = run_trait_object_simulation();

    // Build and run arena simulation
    let arena_time = run_arena_simulation();

    // Summary
    println!("\n=== Summary ===");
    println!("Trait Object (Box<dyn AnyBlock>): {:>8.2} ms", trait_obj_time);
    println!("Arena (SlotMap + BlockKind enum): {:>8.2} ms", arena_time);

    let speedup = trait_obj_time / arena_time;
    println!("\nSpeedup: {:.2}x faster with arena", speedup);
}

fn run_trait_object_simulation() -> f64 {
    // Create a moderately complex simulation:
    // const1 + const2 -> adder -> amp -> integrator

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(Constant::new(1.0)),    // 0
        Box::new(Constant::new(2.0)),    // 1
        Box::new(Adder::<2>::new()),     // 2
        Box::new(Amplifier::new(0.5)),   // 3
        Box::new(Integrator::new(0.0)),  // 4
    ];

    let connections = vec![
        Connection::from(0, 0).to(2, 0).build(), // const1 -> adder[0]
        Connection::from(1, 0).to(2, 1).build(), // const2 -> adder[1]
        Connection::from(2, 0).to(3, 0).build(), // adder -> amp
        Connection::from(3, 0).to(4, 0).build(), // amp -> integrator
    ];

    let mut sim = Simulation::new(blocks, connections).with_dt(DT);

    // Warm up
    for _ in 0..1000 {
        sim.step();
    }
    sim.reset();

    // Timed run
    let start = Instant::now();
    for _ in 0..NUM_STEPS {
        sim.step();
    }
    let elapsed = start.elapsed();

    let ms = elapsed.as_secs_f64() * 1000.0;
    println!("Trait Object Simulation:");
    println!("  Time: {:.2} ms", ms);
    println!("  Final integrator output: {:.6}", sim.get_output(4, 0));

    ms
}

fn run_arena_simulation() -> f64 {
    // Same simulation using arena

    let mut sim = ArenaSimulation::new().with_dt(DT);

    let const1 = sim.add_block(Constant::new(1.0));
    let const2 = sim.add_block(Constant::new(2.0));
    let adder = sim.add_block(Adder::<2>::new());
    let amp = sim.add_block(Amplifier::new(0.5));
    let integrator = sim.add_block(Integrator::new(0.0));

    sim.connect(const1, 0, adder, 0);
    sim.connect(const2, 0, adder, 1);
    sim.connect(adder, 0, amp, 0);
    sim.connect(amp, 0, integrator, 0);

    sim.build();

    // Warm up
    for _ in 0..1000 {
        sim.step();
    }
    sim.reset();

    // Timed run
    let start = Instant::now();
    for _ in 0..NUM_STEPS {
        sim.step();
    }
    let elapsed = start.elapsed();

    let ms = elapsed.as_secs_f64() * 1000.0;
    println!("\nArena Simulation:");
    println!("  Time: {:.2} ms", ms);
    println!("  Final integrator output: {:.6}", sim.get_output(integrator, 0));

    ms
}
