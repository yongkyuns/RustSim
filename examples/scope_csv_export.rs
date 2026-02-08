//! Scope CSV Export Demonstration
//!
//! This example demonstrates using the Scope block to record simulation data
//! and export it to CSV files, matching PathSim's Scope.save() functionality.
//!
//! Demonstrates:
//! - Recording multi-channel data during simulation
//! - Exporting to CSV with default labels
//! - Exporting to CSV with custom labels
//! - Reading and verifying exported data

use rustsim::prelude::*;
use std::f64::consts::PI;

/// Simple harmonic oscillator for demonstration
struct SimpleOscillator {
    position: f64,
    velocity: f64,
    omega: f64, // Angular frequency
}

impl SimpleOscillator {
    fn new(x0: f64, v0: f64, omega: f64) -> Self {
        Self {
            position: x0,
            velocity: v0,
            omega,
        }
    }

    fn step(&mut self, dt: f64) {
        // Simple Euler integration for demonstration
        // dx/dt = v
        // dv/dt = -omega^2 * x
        let acceleration = -self.omega * self.omega * self.position;

        self.position += self.velocity * dt;
        self.velocity += acceleration * dt;
    }

    fn position(&self) -> f64 {
        self.position
    }

    fn velocity(&self) -> f64 {
        self.velocity
    }

    fn energy(&self) -> f64 {
        // Total energy: KE + PE
        0.5 * self.velocity * self.velocity + 0.5 * self.omega * self.omega * self.position * self.position
    }
}

fn main() {
    println!("Scope CSV Export Demonstration");
    println!("==============================\n");

    // Simulation parameters
    let omega = 2.0 * PI; // 1 Hz oscillation
    let dt = 0.01;
    let t_final = 2.0; // Simulate for 2 seconds
    let steps = (t_final / dt) as usize;

    println!("Simulating harmonic oscillator:");
    println!("  Frequency: {} Hz", omega / (2.0 * PI));
    println!("  Time step: {} s", dt);
    println!("  Duration:  {} s", t_final);
    println!();

    // Create oscillator and scope
    let mut osc = SimpleOscillator::new(1.0, 0.0, omega);
    let mut scope = Scope::<3, 1000>::new(); // 3 channels: position, velocity, energy

    // Run simulation and record data
    println!("Running simulation...");
    for i in 0..=steps {
        let t = i as f64 * dt;

        // Set scope inputs
        scope.set_input(0, osc.position());
        scope.set_input(1, osc.velocity());
        scope.set_input(2, osc.energy());

        // Update and record
        scope.update(t);
        scope.step(t, dt);

        // Step oscillator
        osc.step(dt);
    }

    println!("Simulation complete. Recorded {} samples.\n", scope.len());

    // Example 1: Save with default labels
    println!("Example 1: Saving with default labels");
    println!("--------------------------------------");
    let filename1 = "oscillator_default.csv";
    match scope.save(filename1) {
        Ok(_) => println!("✓ Saved to '{}'", filename1),
        Err(e) => println!("✗ Error saving: {}", e),
    }

    // Show header format
    println!("  Header: time [s],port 0,port 1,port 2");
    println!();

    // Example 2: Save with custom labels
    println!("Example 2: Saving with custom labels");
    println!("-------------------------------------");
    let filename2 = "oscillator_labeled.csv";
    let labels = ["position", "velocity", "energy"];
    match scope.save_with_labels(filename2, &labels) {
        Ok(_) => println!("✓ Saved to '{}'", filename2),
        Err(e) => println!("✗ Error saving: {}", e),
    }

    // Show header format
    println!("  Header: time [s],position,velocity,energy");
    println!();

    // Example 3: Save to memory buffer (useful for testing)
    println!("Example 3: Saving to in-memory buffer");
    println!("--------------------------------------");
    let mut buffer = Vec::new();
    match scope.save_to_writer(&mut buffer, &labels) {
        Ok(_) => {
            let csv_string = String::from_utf8(buffer).unwrap();
            let line_count = csv_string.lines().count();
            println!("✓ Written to memory buffer");
            println!("  Size: {} bytes", csv_string.len());
            println!("  Lines: {} (1 header + {} data rows)", line_count, line_count - 1);
        }
        Err(e) => println!("✗ Error writing: {}", e),
    }
    println!();

    // Show some statistics from the recorded data
    println!("Data Statistics");
    println!("---------------");
    let data = scope.data();

    if !data.is_empty() {
        let first = &data[0];
        let last = &data[data.len() - 1];

        println!("  First sample:");
        println!("    Time:     {:.4} s", first.0);
        println!("    Position: {:.6}", first.1[0]);
        println!("    Velocity: {:.6}", first.1[1]);
        println!("    Energy:   {:.6}", first.1[2]);
        println!();

        println!("  Last sample:");
        println!("    Time:     {:.4} s", last.0);
        println!("    Position: {:.6}", last.1[0]);
        println!("    Velocity: {:.6}", last.1[1]);
        println!("    Energy:   {:.6}", last.1[2]);
        println!();

        // Find max/min positions
        let mut max_pos = f64::NEG_INFINITY;
        let mut min_pos = f64::INFINITY;
        for (_, values) in &data {
            max_pos = max_pos.max(values[0]);
            min_pos = min_pos.min(values[0]);
        }

        println!("  Position range: [{:.6}, {:.6}]", min_pos, max_pos);
    }

    println!();
    println!("CSV files created:");
    println!("  - {}", filename1);
    println!("  - {}", filename2);
    println!();
    println!("You can now:");
    println!("  - Open these files in Excel, LibreOffice, or any CSV viewer");
    println!("  - Import them into Python/MATLAB for further analysis");
    println!("  - Plot them using your favorite plotting tool");
}
