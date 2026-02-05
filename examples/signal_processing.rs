//! Signal processing example using v2 architecture
//!
//! Demonstrates source blocks and mathematical operations
//!
//! Creates: y = 2*sin(t) + 1

use rustsim::{Adder, Amplifier, Block, Constant, Sinusoidal};

/// Signal processing example
struct SignalProcessor {
    sine: Sinusoidal,
    amplifier: Amplifier,
    constant: Constant,
    adder: Adder<2>,
    time: f64,
}

impl SignalProcessor {
    fn new() -> Self {
        Self {
            sine: Sinusoidal::new(1.0, 1.0, 0.0), // sin(2*pi*1*t)
            amplifier: Amplifier::new(2.0),       // gain = 2
            constant: Constant::new(1.0),
            adder: Adder::<2>::new(),
            time: 0.0,
        }
    }

    #[inline]
    fn propagate(&mut self) {
        // sine -> amplifier
        self.amplifier.set_input(0, self.sine.get_output(0));

        // amplifier -> adder[0]
        self.adder.inputs_mut()[0] = self.amplifier.get_output(0);

        // constant -> adder[1]
        self.adder.inputs_mut()[1] = self.constant.get_output(0);
    }

    #[inline]
    fn update(&mut self) {
        let t = self.time;

        // Update in topological order
        self.sine.update(t);
        self.constant.update(t);
        self.propagate();

        self.amplifier.update(t);
        self.propagate();

        self.adder.update(t);
        self.propagate();
    }

    fn step(&mut self, dt: f64) {
        self.update();
        self.time += dt;
    }

    fn output(&self) -> f64 {
        self.adder.get_output(0)
    }

    fn time(&self) -> f64 {
        self.time
    }
}

fn main() {
    println!("Signal Processing - v2 Architecture Demo");
    println!("========================================");
    println!();
    println!("Signal: y = 2*sin(2*pi*t) + 1");
    println!();

    let mut sim = SignalProcessor::new();
    let dt = 0.01;

    println!("{:>10} {:>12} {:>12} {:>12}", "Time", "Output", "Exact", "Error");
    println!("{:-<10} {:-<12} {:-<12} {:-<12}", "", "", "", "");

    // Simulate for one period
    let period = 1.0;
    let steps = (period / dt) as usize;

    for i in 0..=steps {
        if i % (steps / 20) == 0 {
            let t = sim.time();
            let y = sim.output();
            let exact = 2.0 * (2.0 * std::f64::consts::PI * t).sin() + 1.0;
            let error = (y - exact).abs();
            println!("{:10.4} {:12.6} {:12.6} {:12.2e}", t, y, exact, error);
        }
        sim.step(dt);
    }

    println!();
    println!("Demo complete!");
}
