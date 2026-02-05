///! FIR Filter Example
///!
///! Demonstrates using the FIR filter block with different configurations:
///! - Passthrough (single coefficient)
///! - Gain
///! - Moving average
///! - Differentiator

use rustsim::prelude::*;

fn main() {
    println!("=== FIR Filter Examples ===\n");

    // Example 1: Passthrough filter
    println!("1. Passthrough Filter [1.0]");
    let mut passthrough = FIR::<1>::new([1.0], 1.0, 0.0);
    passthrough.set_input(0, 5.0);
    passthrough.step(0.0, 1.0);
    println!("   Input: 5.0 -> Output: {}\n", passthrough.get_output(0));

    // Example 2: Gain filter
    println!("2. Gain Filter [2.0]");
    let mut gain = FIR::<1>::new([2.0], 1.0, 0.0);
    gain.set_input(0, 3.0);
    gain.step(0.0, 1.0);
    println!("   Input: 3.0 -> Output: {}\n", gain.get_output(0));

    // Example 3: Moving average filter
    println!("3. Moving Average Filter [1/3, 1/3, 1/3]");
    let mut ma = FIR::<3>::new([1.0/3.0, 1.0/3.0, 1.0/3.0], 1.0, 0.0);

    let inputs = vec![3.0, 6.0, 9.0, 12.0, 15.0];
    println!("   Inputs: {:?}", inputs);
    print!("   Outputs: [");

    for (i, &val) in inputs.iter().enumerate() {
        ma.set_input(0, val);
        ma.step(i as f64, 1.0);
        print!("{:.2}", ma.get_output(0));
        if i < inputs.len() - 1 {
            print!(", ");
        }
    }
    println!("]\n");

    // Example 4: Differentiator (approximates derivative)
    println!("4. Differentiator [1.0, -1.0]");
    let mut diff = FIR::<2>::new([1.0, -1.0], 1.0, 0.0);

    let ramp = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    println!("   Inputs (ramp): {:?}", ramp);
    print!("   Outputs (derivative): [");

    for (i, &val) in ramp.iter().enumerate() {
        diff.set_input(0, val);
        diff.step(i as f64, 1.0);
        print!("{:.2}", diff.get_output(0));
        if i < ramp.len() - 1 {
            print!(", ");
        }
    }
    println!("]\n");

    // Example 5: Two-tap FIR filter
    println!("5. Two-tap FIR [1.0, 0.5]");
    let mut two_tap = FIR::<2>::new([1.0, 0.5], 1.0, 0.0);

    let step_inputs = vec![0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0];
    println!("   Inputs (step): {:?}", step_inputs);
    print!("   Outputs: [");

    for (i, &val) in step_inputs.iter().enumerate() {
        two_tap.set_input(0, val);
        two_tap.step(i as f64, 1.0);
        print!("{:.2}", two_tap.get_output(0));
        if i < step_inputs.len() - 1 {
            print!(", ");
        }
    }
    println!("]\n");

    println!("=== Example Complete ===");
}
