//! Demo of v2 dynamic blocks: ODE, Delay, Differentiator, Scope

use rustsim::*;

fn main() {
    println!("=== v2 Dynamic Blocks Demo ===\n");

    // 1. ODE Block - Harmonic Oscillator
    println!("1. ODE Block - Harmonic Oscillator");
    println!("   d²x/dt² = -x, x(0)=1, v(0)=0");

    let mut ode = ODE::<2, _>::new([1.0, 0.0], |_t, state, _inputs, derivs| {
        derivs[0] = state[1]; // dx/dt = v
        derivs[1] = -state[0]; // dv/dt = -x
    });

    let dt = 0.01;
    let period = 2.0 * std::f64::consts::PI;
    let mut t = 0.0;

    while t < period {
        ode.update(t);
        ode.step(t, dt);
        t += dt;
    }

    println!(
        "   After one period: x={:.6}, v={:.6}",
        ode.state_value(0),
        ode.state_value(1)
    );
    println!("   (Should return to x≈1, v≈0)\n");

    // 2. Delay Block
    println!("2. Delay Block");
    let mut delay = Delay::<200>::new(0.5);

    t = 0.0;
    for _ in 0..100 {
        delay.set_input(0, if t > 0.2 { 1.0 } else { 0.0 });
        delay.update(t);
        delay.step(t, dt);
        t += dt;
    }

    println!("   Step input at t=0.2s, delay=0.5s");
    println!(
        "   At t=1.0s, output={:.6} (should be ≈1.0)\n",
        delay.get_output(0)
    );

    // 3. Differentiator Block
    println!("3. Differentiator Block");
    let mut diff = Differentiator::new(0.01);

    t = 0.0;
    for i in 0..200 {
        t = i as f64 * dt;
        diff.set_input(0, t); // Ramp input
        diff.update(t);
        diff.step(t, dt);
    }

    println!("   Ramp input u(t)=t, du/dt=1");
    println!("   Output={:.6} (should be ≈1.0)\n", diff.value());

    // 4. Scope Block
    println!("4. Scope Block - Recording sine wave");
    let mut scope = Scope::<2, 100>::new();

    t = 0.0;
    for i in 0..100 {
        t = i as f64 * dt;
        scope.set_input(0, t.sin());
        scope.set_input(1, t.cos());
        scope.update(t);
        scope.step(t, dt);
    }

    println!("   Recorded {} samples", scope.len());

    let data = scope.channel_data(0);
    if let Some((t_last, val_last)) = data.last() {
        println!("   Last sample: t={:.3}, sin={:.6}", t_last, val_last);
    }

    // 5. Complex example: ODE with input
    println!("\n5. ODE with External Input - Driven oscillator");
    let omega = 2.0;
    let mut driven = ODE::<2, _>::with_inputs([0.0, 0.0], 1, move |_t, state, inputs, derivs| {
        let forcing = inputs[0];
        derivs[0] = state[1];
        derivs[1] = -omega * omega * state[0] + forcing;
    });

    t = 0.0;
    for i in 0..200 {
        t = i as f64 * dt;
        driven.set_input(0, (omega * t).sin()); // Resonant forcing
        driven.update(t);
        driven.step(t, dt);
    }

    println!("   Oscillator with resonant forcing");
    println!(
        "   Position={:.6}, Velocity={:.6}",
        driven.state_value(0),
        driven.state_value(1)
    );

    println!("\n=== All v2 dynamic blocks working! ===");
}
