//! Integration tests for ADC and DAC converters

use rustsim::blocks::{ADC, DAC};
use rustsim::block::Block;

#[test]
fn test_adc_dac_round_trip() {
    // Test round-trip conversion with 8-bit converters
    let mut adc = ADC::<8>::new([0.0, 5.0], 1.0, 0.0);
    let mut dac = DAC::<8>::new([0.0, 5.0], 1.0, 0.0);

    let test_values = vec![0.0, 1.0, 2.5, 4.0, 5.0];

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

        // Should be very close (within quantization error)
        let reconstructed = dac.get_output(0);
        let error = (val - reconstructed).abs();
        let max_error = 5.0 / 255.0; // Maximum quantization error for 8-bit

        assert!(error <= max_error,
                "Round-trip error too large: input={}, output={}, error={}",
                val, reconstructed, error);
    }
}

#[test]
fn test_adc_periodic_sampling() {
    let mut adc = ADC::<4>::new([-1.0, 1.0], 1.0, 0.0);

    // Set input
    adc.set_input(0, 0.5);

    // First sample at t=0
    adc.update(0.0);
    let first_code: Vec<f64> = (0..4).map(|i| adc.get_output(i)).collect();

    // Change input but don't trigger sample (t=0.5 is not a sample time)
    adc.set_input(0, -0.5);
    adc.update(0.5);
    let same_code: Vec<f64> = (0..4).map(|i| adc.get_output(i)).collect();

    // Output should still be the same (held from last sample)
    assert_eq!(first_code, same_code, "ADC should hold previous value between samples");

    // Now sample at t=1.0
    adc.update(1.0);
    let new_code: Vec<f64> = (0..4).map(|i| adc.get_output(i)).collect();

    // Output should be different now
    assert_ne!(first_code, new_code, "ADC should update at sample time");
}

#[test]
fn test_dac_periodic_update() {
    let mut dac = DAC::<4>::new([-1.0, 1.0], 1.0, 0.0);

    // Set code to 5 (binary: 0101)
    dac.set_input(0, 1.0);
    dac.set_input(1, 0.0);
    dac.set_input(2, 1.0);
    dac.set_input(3, 0.0);

    // First update at t=0
    dac.update(0.0);
    let first_output = dac.get_output(0);

    // Change input but don't trigger update (t=0.5 is not an update time)
    dac.set_input(0, 0.0);
    dac.set_input(1, 0.0);
    dac.set_input(2, 0.0);
    dac.set_input(3, 0.0);
    dac.update(0.5);
    let same_output = dac.get_output(0);

    // Output should still be the same (held from last update)
    assert_eq!(first_output, same_output, "DAC should hold previous value between updates");

    // Now update at t=1.0
    dac.update(1.0);
    let new_output = dac.get_output(0);

    // Output should be different now
    assert_ne!(first_output, new_output, "DAC should update at scheduled time");
}

#[test]
fn test_adc_with_delay() {
    let mut adc = ADC::<4>::new([-1.0, 1.0], 1.0, 0.5);

    adc.set_input(0, 0.5);

    // Before delay, should not sample
    adc.update(0.0);
    let initial_code: Vec<f64> = (0..4).map(|i| adc.get_output(i)).collect();

    // At delay time, should sample
    adc.update(0.5);
    let sampled_code: Vec<f64> = (0..4).map(|i| adc.get_output(i)).collect();

    // Should have sampled
    assert_ne!(initial_code, sampled_code, "ADC should sample after delay");
}

#[test]
fn test_dac_with_delay() {
    let mut dac = DAC::<4>::new([-1.0, 1.0], 1.0, 0.5);

    // Set some code
    dac.set_input(0, 1.0);
    dac.set_input(1, 1.0);
    dac.set_input(2, 0.0);
    dac.set_input(3, 0.0);

    // Before delay, should not update
    dac.update(0.0);
    assert_eq!(dac.get_output(0), 0.0, "DAC should not update before delay");

    // At delay time, should update
    dac.update(0.5);
    assert_ne!(dac.get_output(0), 0.0, "DAC should update after delay");
}

#[test]
fn test_adc_all_codes() {
    // Test that a 2-bit ADC produces all 4 possible codes
    let mut adc = ADC::<2>::new([0.0, 1.0], 1.0, 0.0);

    let test_inputs = vec![0.0, 0.3, 0.6, 0.9];
    let expected_codes = vec![
        vec![0.0, 0.0], // 0
        vec![1.0, 0.0], // 1
        vec![0.0, 1.0], // 2
        vec![1.0, 1.0], // 3
    ];

    for (input, expected) in test_inputs.iter().zip(expected_codes.iter()) {
        adc.set_input(0, *input);
        adc.update(0.0);

        let actual: Vec<f64> = (0..2).map(|i| adc.get_output(i)).collect();
        assert_eq!(actual, *expected, "Input {} should produce code {:?}", input, expected);
    }
}

#[test]
fn test_dac_all_codes() {
    // Test that a 2-bit DAC produces all 4 possible outputs
    let mut dac = DAC::<2>::new([0.0, 3.0], 1.0, 0.0);

    let codes = vec![
        (0.0, 0.0), // 0
        (1.0, 0.0), // 1
        (0.0, 1.0), // 2
        (1.0, 1.0), // 3
    ];

    let expected_outputs = vec![0.0, 1.0, 2.0, 3.0];

    for ((bit0, bit1), expected) in codes.iter().zip(expected_outputs.iter()) {
        dac.set_input(0, *bit0);
        dac.set_input(1, *bit1);
        dac.update(0.0);

        let actual = dac.get_output(0);
        assert!((actual - expected).abs() < 1e-10,
                "Code ({}, {}) should produce output {}, got {}",
                bit0, bit1, expected, actual);
    }
}
