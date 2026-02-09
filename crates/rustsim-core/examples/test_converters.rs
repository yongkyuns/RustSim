//! Example demonstrating ADC and DAC converters
//!
//! This example demonstrates the ADC (Analog-to-Digital) and DAC (Digital-to-Analog)
//! converter blocks using the new Simulation API.

use rustsim_core::prelude::*;

fn main() {
    println!("=== ADC and DAC Converter Examples ===\n");

    example_1_adc_basic();
    example_2_dac_basic();
    example_3_round_trip();
}

/// Example 1: Basic ADC operation
fn example_1_adc_basic() {
    println!("Example 1: ADC (Analog-to-Digital Converter)");
    println!("---------------------------------------------");

    // Test with different constant inputs
    let test_values = vec![-1.0, 0.0, 1.0];

    for val in test_values {
        // Create a 4-bit ADC with span [-1, 1], period 1.0, delay 0
        let constant = Constant::new(val);
        let adc = ADC::<4>::new([-1.0, 1.0], 1.0, 0.0);

        let blocks: Vec<Box<dyn AnyBlock>> = vec![
            Box::new(constant),
            Box::new(adc),
        ];

        let connections = vec![
            Connection::from(0, 0).to(1, 0).build(), // Constant -> ADC
        ];

        let mut sim = Simulation::new(blocks, connections).with_dt(0.1);
        sim.run(0.1);

        println!("\nInput: {:.1}", val);
        println!(
            "Digital code (LSB to MSB): [{:.0}, {:.0}, {:.0}, {:.0}]",
            sim.get_output(1, 0),
            sim.get_output(1, 1),
            sim.get_output(1, 2),
            sim.get_output(1, 3)
        );
    }
    println!();
}

/// Example 2: Basic DAC operation
fn example_2_dac_basic() {
    println!("Example 2: DAC (Digital-to-Analog Converter)");
    println!("---------------------------------------------");

    // Test with different digital codes
    let test_codes = vec![
        (0, 0, 0, 0, "code 0"),
        (1, 1, 1, 1, "code 15"),
        (0, 0, 0, 1, "code 8"),
    ];

    for (b0, b1, b2, b3, description) in test_codes {
        // Create constant blocks for each bit
        let bit0 = Constant::new(b0 as f64);
        let bit1 = Constant::new(b1 as f64);
        let bit2 = Constant::new(b2 as f64);
        let bit3 = Constant::new(b3 as f64);
        let dac = DAC::<4>::new([-1.0, 1.0], 1.0, 0.0);

        let blocks: Vec<Box<dyn AnyBlock>> = vec![
            Box::new(bit0),
            Box::new(bit1),
            Box::new(bit2),
            Box::new(bit3),
            Box::new(dac),
        ];

        let connections = vec![
            Connection::from(0, 0).to(4, 0).build(), // Bit 0 -> DAC
            Connection::from(1, 0).to(4, 1).build(), // Bit 1 -> DAC
            Connection::from(2, 0).to(4, 2).build(), // Bit 2 -> DAC
            Connection::from(3, 0).to(4, 3).build(), // Bit 3 -> DAC
        ];

        let mut sim = Simulation::new(blocks, connections).with_dt(0.1);
        sim.run(0.1);

        println!("\nDigital {} [{}, {}, {}, {}]:", description, b0, b1, b2, b3);
        println!("Analog output: {:.6}", sim.get_output(4, 0));
    }
    println!();
}

/// Example 3: Round-trip conversion (ADC -> DAC)
fn example_3_round_trip() {
    println!("Example 3: Round-trip ADC -> DAC Conversion");
    println!("--------------------------------------------");

    let test_values = vec![0.0, 1.25, 2.5, 3.75, 5.0];

    println!("8-bit converters with span [0.0, 5.0]\n");

    for val in test_values {
        // Create signal source
        let source = Constant::new(val);
        let adc = ADC::<8>::new([0.0, 5.0], 1.0, 0.0);
        let dac = DAC::<8>::new([0.0, 5.0], 1.0, 0.0);

        let blocks: Vec<Box<dyn AnyBlock>> = vec![
            Box::new(source),
            Box::new(adc),
            Box::new(dac),
        ];

        // Connect source to ADC, then all 8 ADC outputs to DAC inputs
        let mut connections = vec![
            Connection::from(0, 0).to(1, 0).build(), // Source -> ADC
        ];

        // Connect all 8 bits from ADC to DAC
        for i in 0..8 {
            connections.push(Connection::from(1, i).to(2, i).build());
        }

        let mut sim = Simulation::new(blocks, connections).with_dt(0.1);
        sim.run(0.1);

        println!(
            "Original: {:.2}, Reconstructed: {:.4}",
            val,
            sim.get_output(2, 0)
        );
    }
    println!();
}
