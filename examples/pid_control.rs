//! PID controller example using v2 architecture
//!
//! Demonstrates a PID controller regulating a first-order plant
//!
//! System:
//!   Plant: dx/dt = -x + u (first-order lag)
//!   Controller: u = PID(setpoint - x)
//!   Goal: x -> setpoint

use rustsim::{Adder, Block, Constant, Integrator, PID};

/// PID-controlled first-order system
///
/// Block diagram:
/// ```text
/// setpoint ─┬─> [error] ─> [PID] ─> [plant] ─> x
///           │                                   │
///           └───────────────────────────────────┘
/// ```
struct PIDSystem {
    setpoint: Constant,
    error_calc: Adder<2>,  // setpoint - measurement
    controller: PID,
    plant: Integrator,     // dx/dt = -x + u
    plant_gain: Adder<2>,  // -x + u
    time: f64,
}

impl PIDSystem {
    fn new(kp: f64, ki: f64, kd: f64, setpoint: f64, x0: f64) -> Self {
        Self {
            setpoint: Constant::new(setpoint),
            error_calc: Adder::<2>::with_weights([1.0, -1.0]), // setpoint - x
            controller: PID::new(kp, ki, kd),
            plant: Integrator::new(x0),
            plant_gain: Adder::<2>::with_weights([-1.0, 1.0]), // -x + u
            time: 0.0,
        }
    }

    #[inline]
    fn propagate(&mut self) {
        // Calculate error: setpoint - x
        self.error_calc.inputs_mut()[0] = self.setpoint.get_output(0);
        self.error_calc.inputs_mut()[1] = self.plant.get_output(0);

        // Controller input is error
        self.controller.set_input(0, self.error_calc.get_output(0));

        // Plant dynamics: dx/dt = -x + u
        self.plant_gain.inputs_mut()[0] = self.plant.get_output(0);
        self.plant_gain.inputs_mut()[1] = self.controller.get_output(0);

        // Plant integrator input
        self.plant.set_input(0, self.plant_gain.get_output(0));
    }

    #[inline]
    fn update(&mut self) {
        let t = self.time;

        // Update in topological order
        self.setpoint.update(t);
        self.plant.update(t);
        self.propagate();

        self.error_calc.update(t);
        self.propagate();

        self.controller.update(t);
        self.propagate();

        self.plant_gain.update(t);
        self.propagate();
    }

    fn step(&mut self, dt: f64) {
        self.update();

        // Step dynamic blocks
        self.controller.step(self.time, dt);
        self.plant.step(self.time, dt);

        self.time += dt;
    }

    fn output(&self) -> f64 {
        self.plant.get_output(0)
    }

    fn error(&self) -> f64 {
        self.error_calc.get_output(0)
    }

    fn control(&self) -> f64 {
        self.controller.get_output(0)
    }

    fn setpoint(&self) -> f64 {
        self.setpoint.get_output(0)
    }

    fn time(&self) -> f64 {
        self.time
    }

    fn set_setpoint(&mut self, value: f64) {
        self.setpoint.set_value(value);
    }
}

fn main() {
    println!("PID Controller - v2 Architecture Demo");
    println!("======================================");
    println!();
    println!("System: First-order plant dx/dt = -x + u");
    println!("Controller: PID with Kp=2.0, Ki=1.0, Kd=0.5");
    println!();

    // Create PID system
    let kp = 2.0;
    let ki = 1.0;
    let kd = 0.5;
    let setpoint = 1.0;
    let x0 = 0.0;

    let mut sim = PIDSystem::new(kp, ki, kd, setpoint, x0);
    let dt = 0.01;

    println!("Initial setpoint: {}", setpoint);
    println!();
    println!("{:>10} {:>12} {:>12} {:>12} {:>12}", "Time", "Setpoint", "Output", "Error", "Control");
    println!("{:-<10} {:-<12} {:-<12} {:-<12} {:-<12}", "", "", "", "", "");

    // Simulate step response
    let duration = 5.0;
    let steps = (duration / dt) as usize;

    for i in 0..=steps {
        // Print every 0.1 seconds
        if i % (steps / 50) == 0 {
            let t = sim.time();
            let sp = sim.setpoint();
            let y = sim.output();
            let e = sim.error();
            let u = sim.control();
            println!("{:10.4} {:12.6} {:12.6} {:12.6} {:12.6}", t, sp, y, e, u);
        }

        sim.step(dt);
    }

    println!();
    println!("Step response complete:");
    println!("  Final output: {:.6} (setpoint: {})", sim.output(), setpoint);
    println!("  Final error:  {:.6}", sim.error());
    println!();

    // Test setpoint change
    println!("Changing setpoint to 2.0...");
    println!();
    sim.set_setpoint(2.0);

    println!("{:>10} {:>12} {:>12} {:>12} {:>12}", "Time", "Setpoint", "Output", "Error", "Control");
    println!("{:-<10} {:-<12} {:-<12} {:-<12} {:-<12}", "", "", "", "", "");

    let duration2 = 5.0;
    let steps2 = (duration2 / dt) as usize;

    for i in 0..=steps2 {
        if i % (steps2 / 50) == 0 {
            let t = sim.time();
            let sp = sim.setpoint();
            let y = sim.output();
            let e = sim.error();
            let u = sim.control();
            println!("{:10.4} {:12.6} {:12.6} {:12.6} {:12.6}", t, sp, y, e, u);
        }

        sim.step(dt);
    }

    println!();
    println!("Final state:");
    println!("  Output: {:.6} (setpoint: 2.0)", sim.output());
    println!("  Error:  {:.6}", sim.error());
    println!();
    println!("PID controller successfully regulates the plant!");
}
