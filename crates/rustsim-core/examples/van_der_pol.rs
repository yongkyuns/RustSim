//! Van der Pol oscillator example using Connection and Simulation API
//!
//! Demonstrates nonlinear oscillator using Function blocks
//!
//! System: d²x/dt² - μ(1 - x²)dx/dt + x = 0
//!
//! Rewritten as two first-order ODEs:
//!   dx/dt = y
//!   dy/dt = μ(1 - x²)y - x

use rustsim_core::anyblock::AnyBlock;
use rustsim_core::blocks::{Adder, Amplifier, Constant, Integrator, Multiplier, Pow};
use rustsim_core::{Connection, Simulation};

/// Van der Pol oscillator
///
/// Block diagram:
/// ```text
///  x -> [x²] -> [1-x²] -> [μ(1-x²)] -> [*y] -> [+] -> y_integrator -> y
///                                                ^
///                                                |
///                                 x_integrator <-+ <- [-x]
/// ```
fn main() {
    println!("Van der Pol Oscillator - Connection & Simulation API Demo");
    println!("==========================================================");
    println!();
    println!("System: d²x/dt² - μ(1 - x²)dx/dt + x = 0");
    println!();

    let mu = 1.0;
    println!("Parameters: μ = {}, x(0) = 1.0, y(0) = 0.0", mu);
    println!();

    // 1. Create blocks with descriptive names
    let x_integrator = Integrator::new(1.0);  // dx/dt = y, x(0) = 1.0
    let y_integrator = Integrator::new(0.0);  // dy/dt = μ(1-x²)y - x, y(0) = 0.0

    let square = Pow::new(2.0);               // x²
    let one_minus_x2 = Adder::<2>::with_weights([1.0, -1.0]); // 1 - x²
    let constant_one = Constant::new(1.0);    // source of 1 for (1 - x²)
    let mu_gain = Amplifier::new(mu);         // μ(1 - x²)
    let damping_mult = Multiplier::<2>::new(); // μ(1-x²) * y
    let neg_x = Amplifier::new(-1.0);         // -x
    let acceleration = Adder::<2>::new();     // μ(1-x²)y - x

    // 2. Collect blocks (indices match connection references)
    // Block 0: x_integrator
    // Block 1: y_integrator
    // Block 2: square (x²)
    // Block 3: one_minus_x2 (1 - x²)
    // Block 4: constant_one
    // Block 5: mu_gain (μ(1-x²))
    // Block 6: damping_mult (μ(1-x²)y)
    // Block 7: neg_x (-x)
    // Block 8: acceleration (sum)
    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(x_integrator),
        Box::new(y_integrator),
        Box::new(square),
        Box::new(one_minus_x2),
        Box::new(constant_one),
        Box::new(mu_gain),
        Box::new(damping_mult),
        Box::new(neg_x),
        Box::new(acceleration),
    ];

    // 3. Define connections
    let connections = vec![
        // x -> square
        Connection::from(0, 0).to(2, 0).build(),
        // 1 -> one_minus_x2 (first input)
        Connection::from(4, 0).to(3, 0).build(),
        // x² -> one_minus_x2 (second input, subtracted)
        Connection::from(2, 0).to(3, 1).build(),
        // (1-x²) -> mu_gain
        Connection::from(3, 0).to(5, 0).build(),
        // μ(1-x²) -> damping_mult (first input)
        Connection::from(5, 0).to(6, 0).build(),
        // y -> damping_mult (second input)
        Connection::from(1, 0).to(6, 1).build(),
        // x -> neg_x
        Connection::from(0, 0).to(7, 0).build(),
        // μ(1-x²)y -> acceleration (first input)
        Connection::from(6, 0).to(8, 0).build(),
        // -x -> acceleration (second input)
        Connection::from(7, 0).to(8, 1).build(),
        // acceleration -> y_integrator
        Connection::from(8, 0).to(1, 0).build(),
        // y -> x_integrator (dx/dt = y)
        Connection::from(1, 0).to(0, 0).build(),
    ];

    // 4. Create and configure simulation
    let mut sim = Simulation::new(blocks, connections)
        .with_dt(0.01);

    println!("{:>10} {:>12} {:>12}", "Time", "x", "y");
    println!("{:-<10} {:-<12} {:-<12}", "", "", "");

    // Simulate for several periods
    let duration = 20.0;
    let dt = sim.dt();
    let steps = (duration / dt) as usize;

    for i in 0..=steps {
        if i % (steps / 40) == 0 {
            let t = sim.time();
            let x = sim.get_output(0, 0);
            let y = sim.get_output(1, 0);
            println!("{:10.4} {:12.6} {:12.6}", t, x, y);
        }
        sim.step();
    }

    let final_x = sim.get_output(0, 0);
    let final_y = sim.get_output(1, 0);

    println!();
    println!("Final state:");
    println!("  x({:.1}) = {:.6}", sim.time(), final_x);
    println!("  y({:.1}) = {:.6}", sim.time(), final_y);
    println!();
    println!("Van der Pol oscillator exhibits a limit cycle.");
    println!("Try different μ values to see stiff (large μ) vs relaxation (small μ) oscillations.");
}
