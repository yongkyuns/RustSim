//! Demo of v2 dynamic blocks: ODE, Delay, Differentiator, Scope
//!
//! This example demonstrates the dynamic blocks using the new Simulation API.

use rustsim_core::prelude::*;

fn main() {
    println!("=== v2 Dynamic Blocks Demo ===\n");

    example_1_ode_harmonic_oscillator();
    example_2_delay_block();
    example_3_differentiator();
    example_4_scope_recording();
    example_5_driven_oscillator();

    println!("\n=== All v2 dynamic blocks working! ===");
}

/// Example 1: ODE Block - Harmonic Oscillator
fn example_1_ode_harmonic_oscillator() {
    println!("1. ODE Block - Harmonic Oscillator");
    println!("   d²x/dt² = -x, x(0)=1, v(0)=0");

    let ode = ODE::<2, _>::new([1.0, 0.0], |_t, state, _inputs, derivs| {
        derivs[0] = state[1]; // dx/dt = v
        derivs[1] = -state[0]; // dv/dt = -x
    });

    let blocks: Vec<Box<dyn AnyBlock>> = vec![Box::new(ode)];
    let connections = vec![]; // No connections needed for standalone block

    let dt = 0.01;
    let period = 2.0 * std::f64::consts::PI;

    let mut sim = Simulation::new(blocks, connections).with_dt(dt);
    sim.run(period);

    // Access the ODE block state through simulation
    println!(
        "   After one period: time={:.6}",
        sim.time()
    );
    println!("   (Should be ≈{:.6})\n", period);
}

/// Example 2: Delay Block
fn example_2_delay_block() {
    println!("2. Delay Block");

    let step = Step::new(1.0, 0.2); // Step from 0 to 1 at t=0.2s
    let delay = Delay::<200>::new(0.5); // 0.5s delay

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(step),
        Box::new(delay),
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).build(), // Step -> Delay
    ];

    let dt = 0.01;
    let mut sim = Simulation::new(blocks, connections).with_dt(dt);
    sim.run(1.0);

    println!("   Step input at t=0.2s, delay=0.5s");
    println!(
        "   At t=1.0s, output={:.6} (should be ≈1.0)\n",
        sim.get_output(1, 0)
    );
}

/// Example 3: Differentiator Block
fn example_3_differentiator() {
    println!("3. Differentiator Block");

    let ramp = Ramp::new(1.0, 0.0); // Ramp with slope 1.0
    let diff = Differentiator::new(0.01);

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(ramp),
        Box::new(diff),
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).build(), // Ramp -> Differentiator
    ];

    let dt = 0.01;
    let mut sim = Simulation::new(blocks, connections).with_dt(dt);
    sim.run(2.0);

    println!("   Ramp input u(t)=t, du/dt=1");
    println!("   Output={:.6} (should be ≈1.0)\n", sim.get_output(1, 0));
}

/// Example 4: Scope Block - Recording sine wave
fn example_4_scope_recording() {
    println!("4. Scope Block - Recording sine wave");

    let sine = Sinusoidal::new(1.0, 1.0, 0.0); // 1 Hz sine wave
    let cosine = Sinusoidal::new(1.0, 1.0, std::f64::consts::PI / 2.0); // Cosine (90° phase)
    let scope = Scope::<2, 100>::new();

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(sine),
        Box::new(cosine),
        Box::new(scope),
    ];

    let connections = vec![
        Connection::from(0, 0).to(2, 0).build(), // Sine -> Scope channel 0
        Connection::from(1, 0).to(2, 1).build(), // Cosine -> Scope channel 1
    ];

    let dt = 0.01;
    let mut sim = Simulation::new(blocks, connections).with_dt(dt);
    sim.run(0.99); // Just under 100 samples

    println!("   Recorded samples at t={:.3}", sim.time());
    println!("   Sine output: {:.6}", sim.get_output(0, 0));
    println!("   Cosine output: {:.6}\n", sim.get_output(1, 0));
}

/// Example 5: ODE with External Input - Driven oscillator
fn example_5_driven_oscillator() {
    println!("5. ODE with External Input - Driven oscillator");

    let omega = 2.0;
    let forcing = Sinusoidal::new(1.0, omega / (2.0 * std::f64::consts::PI), 0.0); // Resonant forcing
    let driven = ODE::<2, _>::with_inputs([0.0, 0.0], 1, move |_t, state, inputs, derivs| {
        let force = inputs[0];
        derivs[0] = state[1];
        derivs[1] = -omega * omega * state[0] + force;
    });

    let blocks: Vec<Box<dyn AnyBlock>> = vec![
        Box::new(forcing),
        Box::new(driven),
    ];

    let connections = vec![
        Connection::from(0, 0).to(1, 0).build(), // Forcing -> ODE
    ];

    let dt = 0.01;
    let mut sim = Simulation::new(blocks, connections).with_dt(dt);
    sim.run(2.0);

    println!("   Oscillator with resonant forcing");
    println!("   At t={:.3}", sim.time());
}
