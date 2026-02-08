//! Van der Pol oscillator example using v2 architecture
//!
//! Demonstrates nonlinear oscillator using Function blocks
//!
//! System: d²x/dt² - μ(1 - x²)dx/dt + x = 0
//!
//! Rewritten as two first-order ODEs:
//!   dx/dt = y
//!   dy/dt = μ(1 - x²)y - x

use rustsim::{Adder, Amplifier, Block, Function, Integrator, Pow};

/// Van der Pol oscillator
///
/// Block diagram:
/// ```text
///  x -> [x²] -> [1-x²] -> [μ(1-x²)y] -> [+] -> y_integrator -> y
///                                        ^
///                                        |
///                         x_integrator <-+ <- [-x]
/// ```
struct VanDerPol {
    // State integrators
    x_integrator: Integrator, // dx/dt = y
    y_integrator: Integrator, // dy/dt = μ(1-x²)y - x

    // Nonlinear term: μ(1 - x²)y
    square: Pow,                                           // x²
    one_minus_x2: Adder<2>,                                // 1 - x²
    damping: Function<2, 1, fn(&[f64; 2], &mut [f64; 1])>, // μ(1-x²)y

    // Force term: -x
    neg_x: Amplifier,

    // Acceleration: μ(1-x²)y - x
    acceleration: Adder<2>,

    // Parameters
    mu: f64,
    time: f64,
}

impl VanDerPol {
    fn new(mu: f64, x0: f64, y0: f64) -> Self {
        Self {
            x_integrator: Integrator::new(x0),
            y_integrator: Integrator::new(y0),
            square: Pow::new(2.0),
            one_minus_x2: Adder::<2>::with_weights([1.0, -1.0]),
            damping: Function::new(|inputs, outputs| {
                // inputs[0] = 1 - x²
                // inputs[1] = y
                outputs[0] = inputs[0] * inputs[1];
            }),
            neg_x: Amplifier::new(-1.0),
            acceleration: Adder::<2>::new(),
            mu,
            time: 0.0,
        }
    }

    #[inline]
    fn propagate(&mut self) {
        // x -> square
        self.square.set_input(0, self.x_integrator.get_output(0));

        // 1 - x²
        self.one_minus_x2.inputs_mut()[0] = 1.0; // constant 1
        self.one_minus_x2.inputs_mut()[1] = self.square.get_output(0);

        // μ(1 - x²)y
        self.damping.inputs_mut()[0] = self.one_minus_x2.get_output(0);
        self.damping.inputs_mut()[1] = self.y_integrator.get_output(0);

        // -x
        self.neg_x.set_input(0, self.x_integrator.get_output(0));

        // acceleration = μ(1-x²)y - x
        self.acceleration.inputs_mut()[0] = self.damping.get_output(0);
        self.acceleration.inputs_mut()[1] = self.neg_x.get_output(0);

        // Connect integrators
        self.x_integrator
            .set_input(0, self.y_integrator.get_output(0)); // dx/dt = y
        self.y_integrator
            .set_input(0, self.acceleration.get_output(0)); // dy/dt = acceleration
    }

    #[inline]
    fn update(&mut self) {
        let t = self.time;

        // Update in topological order
        self.x_integrator.update(t);
        self.y_integrator.update(t);
        self.propagate();

        self.square.update(t);
        self.propagate();

        self.one_minus_x2.update(t);
        self.propagate();

        self.damping.update(t);
        self.propagate();

        self.neg_x.update(t);
        self.propagate();

        self.acceleration.update(t);
        self.propagate();
    }

    fn step(&mut self, dt: f64) {
        // Update algebraic relationships
        self.update();

        // Apply damping factor μ
        let accel = self.acceleration.get_output(0) * self.mu;
        self.y_integrator.set_input(0, accel);

        // Step dynamic blocks
        self.x_integrator.step(self.time, dt);
        self.y_integrator.step(self.time, dt);

        self.time += dt;
    }

    fn x(&self) -> f64 {
        self.x_integrator.get_output(0)
    }

    fn y(&self) -> f64 {
        self.y_integrator.get_output(0)
    }

    fn time(&self) -> f64 {
        self.time
    }
}

fn main() {
    println!("Van der Pol Oscillator - v2 Architecture Demo");
    println!("==============================================");
    println!();
    println!("System: d²x/dt² - μ(1 - x²)dx/dt + x = 0");
    println!();

    let mu = 1.0;
    let mut sim = VanDerPol::new(mu, 1.0, 0.0);
    let dt = 0.01;

    println!("Parameters: μ = {}, x(0) = 1.0, y(0) = 0.0", mu);
    println!();
    println!("{:>10} {:>12} {:>12}", "Time", "x", "y");
    println!("{:-<10} {:-<12} {:-<12}", "", "", "");

    // Simulate for several periods
    let duration = 20.0;
    let steps = (duration / dt) as usize;

    for i in 0..=steps {
        if i % (steps / 40) == 0 {
            let t = sim.time();
            let x = sim.x();
            let y = sim.y();
            println!("{:10.4} {:12.6} {:12.6}", t, x, y);
        }
        sim.step(dt);
    }

    println!();
    println!("Final state:");
    println!("  x({:.1}) = {:.6}", sim.time(), sim.x());
    println!("  y({:.1}) = {:.6}", sim.time(), sim.y());
    println!();
    println!("Van der Pol oscillator exhibits a limit cycle.");
    println!("Try different μ values to see stiff (large μ) vs relaxation (small μ) oscillations.");
}
