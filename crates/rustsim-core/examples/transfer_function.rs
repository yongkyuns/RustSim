//! Transfer Function examples demonstrating common use cases
//!
//! Shows how to create and simulate various transfer functions

use rustsim::blocks::TransferFunction;
use rustsim::Block;

fn main() {
    println!("Transfer Function Block Examples");
    println!("================================\n");

    // Example 1: First-order low-pass filter
    example_1_first_order_lowpass();

    // Example 2: Second-order system
    example_2_second_order();

    // Example 3: Lead compensator
    example_3_lead_compensator();

    // Example 4: Pure integrator
    example_4_integrator();
}

fn example_1_first_order_lowpass() {
    println!("Example 1: First-Order Low-Pass Filter");
    println!("---------------------------------------");
    println!("H(s) = 1/(s+1), Time constant τ = 1");
    println!();

    let mut tf = TransferFunction::<1>::new(&[1.0], &[1.0, 1.0]);

    // Print state-space matrices
    println!("State-space realization:");
    println!("A = {:?}", tf.a_matrix());
    println!("B = {:?}", tf.b_matrix());
    println!("C = {:?}", tf.c_matrix());
    println!("D = {:?}", tf.d_matrix());
    println!();

    // Apply step input
    tf.set_input(0, 1.0);
    tf.update(0.0);

    let dt = 0.01_f64;
    let mut t = 0.0_f64;

    println!("{:>8} {:>12} {:>12} {:>12}", "Time", "Output", "Exact", "Error");
    println!("{:-<8} {:-<12} {:-<12} {:-<12}", "", "", "", "");

    // Exact solution: y(t) = 1 - e^(-t)
    for _ in 0..500 {
        if t.rem_euclid(1.0) < dt {
            let exact = 1.0 - (-t).exp();
            let output = tf.get_output(0);
            let error = (output - exact).abs();
            println!("{:8.2} {:12.6} {:12.6} {:12.2e}", t, output, exact, error);
        }
        tf.step(t, dt);
        t += dt;
    }
    println!();
}

fn example_2_second_order() {
    println!("Example 2: Second-Order System");
    println!("-------------------------------");
    println!("H(s) = 1/(s^2 + 2s + 1) = 1/(s+1)^2");
    println!();

    let mut tf = TransferFunction::<2>::new(&[1.0], &[1.0, 2.0, 1.0]);

    println!("State-space realization:");
    println!("A = {:?}", tf.a_matrix());
    println!("B = {:?}", tf.b_matrix());
    println!("C = {:?}", tf.c_matrix());
    println!("D = {:?}", tf.d_matrix());
    println!();

    // Apply step input
    tf.set_input(0, 1.0);
    tf.update(0.0);

    let dt = 0.01_f64;
    let mut t = 0.0_f64;

    println!("{:>8} {:>12} {:>12} {:>12}", "Time", "Output", "Exact", "Error");
    println!("{:-<8} {:-<12} {:-<12} {:-<12}", "", "", "", "");

    // Exact solution: y(t) = 1 - (1+t)e^(-t)
    for _ in 0..500 {
        if t.rem_euclid(1.0) < dt {
            let exact = 1.0 - (1.0 + t) * (-t).exp();
            let output = tf.get_output(0);
            let error = (output - exact).abs();
            println!("{:8.2} {:12.6} {:12.6} {:12.2e}", t, output, exact, error);
        }
        tf.step(t, dt);
        t += dt;
    }
    println!();
}

fn example_3_lead_compensator() {
    println!("Example 3: Lead Compensator");
    println!("----------------------------");
    println!("H(s) = (s+2)/(s+10)");
    println!();

    let mut tf = TransferFunction::<1>::new(&[1.0, 2.0], &[1.0, 10.0]);

    println!("Has direct feedthrough: {}", tf.has_passthrough());
    println!("State-space realization:");
    println!("A = {:?}", tf.a_matrix());
    println!("B = {:?}", tf.b_matrix());
    println!("C = {:?}", tf.c_matrix());
    println!("D = {:?}", tf.d_matrix());
    println!();

    // Apply step input
    tf.set_input(0, 1.0);
    tf.update(0.0);

    let dt = 0.001_f64;
    let mut t = 0.0_f64;

    println!("{:>8} {:>12} {:>12} {:>12}", "Time", "Output", "Exact", "Error");
    println!("{:-<8} {:-<12} {:-<12} {:-<12}", "", "", "", "");

    // Exact solution: y(t) = 0.2 + 0.8*e^(-10t)
    for _ in 0..1000 {
        if (t * 10.0).rem_euclid(1.0) < dt * 10.0 {
            let exact = 0.2 + 0.8 * (-10.0 * t).exp();
            let output = tf.get_output(0);
            let error = (output - exact).abs();
            println!("{:8.2} {:12.6} {:12.6} {:12.2e}", t, output, exact, error);
        }
        tf.step(t, dt);
        t += dt;
    }
    println!();
}

fn example_4_integrator() {
    println!("Example 4: Pure Integrator");
    println!("--------------------------");
    println!("H(s) = 1/s");
    println!();

    let mut tf = TransferFunction::<1>::new(&[1.0], &[1.0, 0.0]);

    println!("State-space realization:");
    println!("A = {:?}", tf.a_matrix());
    println!("B = {:?}", tf.b_matrix());
    println!("C = {:?}", tf.c_matrix());
    println!("D = {:?}", tf.d_matrix());
    println!();

    // Apply step input
    tf.set_input(0, 1.0);
    tf.update(0.0);

    let dt = 0.01_f64;
    let mut t = 0.0_f64;

    println!("{:>8} {:>12} {:>12} {:>12}", "Time", "Output", "Exact", "Error");
    println!("{:-<8} {:-<12} {:-<12} {:-<12}", "", "", "", "");

    // Exact solution: y(t) = t
    for _ in 0..500 {
        if t.rem_euclid(1.0) < dt {
            let exact = t;
            let output = tf.get_output(0);
            let error = (output - exact).abs();
            println!("{:8.2} {:12.6} {:12.6} {:12.2e}", t, output, exact, error);
        }
        tf.step(t, dt);
        t += dt;
    }
    println!();
}
