//! Simple integration example using v2 architecture
//!
//! Demonstrates basic integrator usage: integrate constant input
//!
//! System: dy/dt = 1, y(0) = 0
//! Solution: y(t) = t

use rustsim::{Block, Constant, Integrator};

/// Simple integrator test: integrate a constant
struct SimpleIntegration {
    source: Constant,
    integrator: Integrator,
    time: f64,
}

impl SimpleIntegration {
    fn new() -> Self {
        Self {
            source: Constant::new(1.0),       // dy/dt = 1
            integrator: Integrator::new(0.0), // y(0) = 0
            time: 0.0,
        }
    }

    #[inline]
    fn propagate(&mut self) {
        // Connect source -> integrator
        self.integrator.set_input(0, self.source.get_output(0));
    }

    #[inline]
    fn update(&mut self) {
        let t = self.time;

        // Update blocks in order
        self.source.update(t);
        self.propagate();

        self.integrator.update(t);
        self.propagate();
    }

    fn step(&mut self, dt: f64) {
        // Update algebraic relationships
        self.update();

        // Step dynamic blocks
        self.integrator.step(self.time, dt);

        // Advance time
        self.time += dt;
    }

    fn output(&self) -> f64 {
        self.integrator.get_output(0)
    }

    fn time(&self) -> f64 {
        self.time
    }
}

fn main() {
    println!("Simple Integration - v2 Architecture Demo");
    println!("=========================================");
    println!();
    println!("System: dy/dt = 1, y(0) = 0");
    println!("Exact solution: y(t) = t");
    println!();

    let mut sim = SimpleIntegration::new();
    let dt = 0.01;

    println!(
        "{:>10} {:>12} {:>12} {:>12}",
        "Time", "Output", "Exact", "Error"
    );
    println!("{:-<10} {:-<12} {:-<12} {:-<12}", "", "", "", "");

    // Simulate for 5 seconds
    let duration = 5.0;
    let steps = (duration / dt) as usize;

    for i in 0..=steps {
        if i % (steps / 10) == 0 {
            let t = sim.time();
            let y = sim.output();
            let exact = t;
            let error = (y - exact).abs();
            println!("{:10.4} {:12.6} {:12.6} {:12.2e}", t, y, exact, error);
        }
        sim.step(dt);
    }

    println!();
    println!("Final value: {:.6} (expected 5.0)", sim.output());
    println!("Final error: {:.2e}", (sim.output() - 5.0).abs());
}
