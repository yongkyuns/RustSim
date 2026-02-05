use rustsim::{Adder, Amplifier, Block, Constant, Integrator};

/// Simple v2 demo: constant -> amplifier -> adder
struct SimpleDemo {
    constant: Constant,
    amplifier: Amplifier,
    adder: Adder<2>,
    integrator: Integrator,
    time: f64,
}

impl SimpleDemo {
    fn new() -> Self {
        Self {
            constant: Constant::new(5.0),
            amplifier: Amplifier::new(2.0),
            adder: Adder::<2>::new(),
            integrator: Integrator::new(0.0),
            time: 0.0,
        }
    }

    fn propagate(&mut self) {
        // constant -> amplifier
        self.amplifier.set_input(0, self.constant.get_output(0));

        // amplifier -> adder[0]
        self.adder.inputs_mut()[0] = self.amplifier.get_output(0);

        // constant -> adder[1]
        self.adder.inputs_mut()[1] = self.constant.get_output(0);

        // adder -> integrator
        self.integrator.set_input(0, self.adder.get_output(0));
    }

    fn update(&mut self) {
        let t = self.time;

        // Update blocks in topological order
        self.constant.update(t);
        self.propagate();

        self.amplifier.update(t);
        self.propagate();

        self.adder.update(t);
        self.propagate();

        self.integrator.update(t);
        self.propagate();
    }

    fn step(&mut self, dt: f64) {
        self.update();
        self.integrator.step(self.time, dt);
        self.time += dt;
    }

    fn output(&self) -> f64 {
        self.integrator.get_output(0)
    }
}

fn main() {
    println!("RustSim v2 - Static Compile-Time Block Architecture");
    println!("====================================================\n");

    println!("Creating simple demo simulation:");
    println!("  Constant(5.0) -> Amplifier(2.0) -> Adder (+Constant) -> Integrator");
    println!("  Expected: 2*5 + 5 = 15, integrated over time\n");

    let mut sim = SimpleDemo::new();
    let dt = 0.01;

    println!("Simulating for 1.0 second with dt = {} s\n", dt);

    let duration = 1.0;
    let steps = (duration / dt) as usize;

    // Run simulation
    for _ in 0..steps {
        sim.step(dt);
    }

    println!("Results:");
    println!("  Constant output:    {:.4}", sim.constant.get_output(0));
    println!("  Amplifier output:   {:.4}", sim.amplifier.get_output(0));
    println!("  Adder output:       {:.4}", sim.adder.get_output(0));
    println!("  Integrator output:  {:.4}", sim.integrator.get_output(0));
    println!();

    // Expected: integrator should integrate 15.0 for 1 second = 15.0
    let expected = 15.0;
    let actual = sim.output();
    let error = (actual - expected).abs();

    println!("Expected integrator value: {:.4}", expected);
    println!("Actual integrator value:   {:.4}", actual);
    println!("Error:                     {:.6}", error);
    println!();

    if error < 1e-6 {
        println!("SUCCESS: v2 simulation is working correctly!");
    } else {
        println!("WARNING: Unexpected error of {:.2e}", error);
    }

    println!();
    println!("Key features of v2 architecture:");
    println!("  - Fixed I/O sizes at compile time (no Vec allocations)");
    println!("  - Static dispatch (no dyn Block trait objects)");
    println!("  - Zero-cost abstractions (fully inlined)");
    println!("  - Manual graph assembly (explicit propagate() calls)");
    println!();
    println!("Run examples for more complex simulations:");
    println!("  cargo run --example v2_harmonic_oscillator");
    println!("  cargo run --example v2_simple_integration");
    println!("  cargo run --example v2_signal_processing");
    println!("  cargo run --example v2_van_der_pol");
    println!("  cargo run --example v2_pid_control");
}
