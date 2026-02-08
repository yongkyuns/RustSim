//! Example demonstrating ADC and DAC converters

use rustsim::block::Block;
use rustsim::blocks::{ADC, DAC};

fn main() {
    println!("=== ADC Example ===");

    // Create a 4-bit ADC with span [-1, 1], period 1.0, delay 0
    let mut adc = ADC::<4>::new([-1.0, 1.0], 1.0, 0.0);

    // Test midrange value
    adc.set_input(0, 0.0);
    adc.update(0.0);
    println!("Input: 0.0");
    println!(
        "Digital code (LSB to MSB): [{}, {}, {}, {}]",
        adc.get_output(0),
        adc.get_output(1),
        adc.get_output(2),
        adc.get_output(3)
    );
    println!("Expected: [0, 0, 0, 1] (code 8)");

    // Test maximum value
    adc.set_input(0, 1.0);
    adc.update(1.0);
    println!("\nInput: 1.0");
    println!(
        "Digital code (LSB to MSB): [{}, {}, {}, {}]",
        adc.get_output(0),
        adc.get_output(1),
        adc.get_output(2),
        adc.get_output(3)
    );
    println!("Expected: [1, 1, 1, 1] (code 15)");

    // Test minimum value
    adc.set_input(0, -1.0);
    adc.update(2.0);
    println!("\nInput: -1.0");
    println!(
        "Digital code (LSB to MSB): [{}, {}, {}, {}]",
        adc.get_output(0),
        adc.get_output(1),
        adc.get_output(2),
        adc.get_output(3)
    );
    println!("Expected: [0, 0, 0, 0] (code 0)");

    println!("\n=== DAC Example ===");

    // Create a 4-bit DAC with span [-1, 1], period 1.0, delay 0
    let mut dac = DAC::<4>::new([-1.0, 1.0], 1.0, 0.0);

    // Test code 0
    dac.set_input(0, 0.0);
    dac.set_input(1, 0.0);
    dac.set_input(2, 0.0);
    dac.set_input(3, 0.0);
    dac.update(0.0);
    println!("Digital code: 0");
    println!("Analog output: {}", dac.get_output(0));
    println!("Expected: -1.0");

    // Test code 15 (all 1s)
    dac.set_input(0, 1.0);
    dac.set_input(1, 1.0);
    dac.set_input(2, 1.0);
    dac.set_input(3, 1.0);
    dac.update(1.0);
    println!("\nDigital code: 15");
    println!("Analog output: {}", dac.get_output(0));
    println!("Expected: 1.0");

    // Test code 8 (midrange)
    dac.set_input(0, 0.0);
    dac.set_input(1, 0.0);
    dac.set_input(2, 0.0);
    dac.set_input(3, 1.0);
    dac.update(2.0);
    println!("\nDigital code: 8");
    println!("Analog output: {}", dac.get_output(0));
    let expected = -1.0 + 2.0 * (8.0 / 15.0);
    println!("Expected: {} (approx 0.0667)", expected);

    println!("\n=== Round-trip Example ===");

    // ADC -> DAC round-trip
    let mut adc = ADC::<8>::new([0.0, 5.0], 1.0, 0.0);
    let mut dac = DAC::<8>::new([0.0, 5.0], 1.0, 0.0);

    let test_values = vec![0.0, 1.25, 2.5, 3.75, 5.0];

    for val in test_values {
        // Sample with ADC
        adc.set_input(0, val);
        adc.update(0.0);

        // Transfer bits to DAC
        for i in 0..8 {
            dac.set_input(i, adc.get_output(i));
        }

        // Convert back to analog
        dac.update(0.0);

        println!(
            "Original: {:.2}, Reconstructed: {:.4}",
            val,
            dac.get_output(0)
        );
    }
}
