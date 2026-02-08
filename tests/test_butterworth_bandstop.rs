//! Comprehensive tests for ButterworthBandstop filter
//!
//! Tests the bandstop (notch) filter implementation against various conditions
//! including frequency response, different orders, and comparison with expected behavior.

use rustsim::prelude::*;
use std::f64::consts::PI;

/// Helper function to compute steady-state amplitude at a given frequency
fn measure_steady_state_amplitude(
    filter: &mut ButterworthBandstop,
    frequency: f64,
    dt: f64,
    settle_steps: usize,
    measure_steps: usize,
) -> f64 {
    let omega = 2.0 * PI * frequency;
    let mut max_output: f64 = 0.0;

    // Run for settling + measurement time
    for i in 0..(settle_steps + measure_steps) {
        let t = i as f64 * dt;
        let input = (omega * t).sin();
        filter.set_input(0, input);
        filter.update(t);
        filter.step(t, dt);

        // Only measure after settling
        if i >= settle_steps {
            max_output = max_output.max(filter.get_output(0).abs());
        }
    }

    max_output
}

#[test]
fn test_bandstop_initialization() {
    let flt = ButterworthBandstop::new([45.0, 55.0], 2);
    assert_eq!(flt.get_output(0), 0.0);
    assert_eq!(flt.get_input(0), 0.0);
}

#[test]
fn test_bandstop_dc_passthrough() {
    // DC (0 Hz) should pass through a bandstop filter
    let mut flt = ButterworthBandstop::new([40.0, 60.0], 2);
    let dt = 0.001;

    // Apply constant DC input
    for _ in 0..2000 {
        flt.set_input(0, 1.5);
        flt.update(0.0);
        flt.step(0.0, dt);
    }

    // Output should approach the DC value
    assert!(
        (flt.get_output(0) - 1.5).abs() < 0.2,
        "DC not passing through: expected ~1.5, got {}",
        flt.get_output(0)
    );
}

#[test]
fn test_bandstop_center_frequency_attenuation() {
    // Test that center frequency is attenuated
    let fc_low: f64 = 45.0;
    let fc_high: f64 = 55.0;
    let fc_center = (fc_low * fc_high).sqrt(); // Geometric mean

    let mut flt = ButterworthBandstop::new([fc_low, fc_high], 4);
    let dt = 0.0001;

    let amplitude = measure_steady_state_amplitude(&mut flt, fc_center, dt, 2000, 3000);

    // Center frequency should be significantly attenuated
    assert!(
        amplitude < 0.2,
        "Center frequency not sufficiently attenuated: {}",
        amplitude
    );
}

#[test]
fn test_bandstop_stopband_rejection() {
    // Test multiple frequencies within the stopband
    let fc_low = 45.0;
    let fc_high = 55.0;

    let mut flt = ButterworthBandstop::new([fc_low, fc_high], 4);
    let dt = 0.0001;

    // Test frequencies within the stopband
    let test_frequencies = vec![46.0, 48.0, 50.0, 52.0, 54.0];

    for freq in test_frequencies {
        flt.reset();
        let amplitude = measure_steady_state_amplitude(&mut flt, freq, dt, 2000, 3000);

        assert!(
            amplitude < 0.5,
            "Frequency {} Hz not sufficiently attenuated: amplitude = {}",
            freq,
            amplitude
        );
    }
}

#[test]
fn test_bandstop_lower_passband() {
    // Test that frequencies well below the stopband pass through
    let mut flt = ButterworthBandstop::new([45.0, 55.0], 2);
    let dt = 0.0001;

    // Test low-frequency passband
    let test_frequencies = vec![5.0, 10.0, 20.0, 30.0];

    for freq in test_frequencies {
        flt.reset();
        let amplitude = measure_steady_state_amplitude(&mut flt, freq, dt, 2000, 3000);

        assert!(
            amplitude > 0.6,
            "Low frequency {} Hz excessively attenuated: amplitude = {}",
            freq,
            amplitude
        );
    }
}

#[test]
fn test_bandstop_upper_passband() {
    // Test that frequencies well above the stopband pass through
    let mut flt = ButterworthBandstop::new([45.0, 55.0], 2);
    let dt = 0.0001;

    // Test high-frequency passband
    let test_frequencies = vec![80.0, 100.0, 150.0, 200.0];

    for freq in test_frequencies {
        flt.reset();
        let amplitude = measure_steady_state_amplitude(&mut flt, freq, dt, 2000, 3000);

        assert!(
            amplitude > 0.5,
            "High frequency {} Hz excessively attenuated: amplitude = {}",
            freq,
            amplitude
        );
    }
}

#[test]
fn test_bandstop_different_orders() {
    // Test that higher-order filters provide better attenuation
    let fc_low = 45.0;
    let fc_high = 55.0;
    let fc_center = 50.0;
    let dt = 0.0001;

    let orders = vec![2, 4, 6];
    let mut attenuations = Vec::new();

    for order in orders {
        let mut flt = ButterworthBandstop::new([fc_low, fc_high], order);
        let amplitude = measure_steady_state_amplitude(&mut flt, fc_center, dt, 2000, 3000);
        attenuations.push(amplitude);
    }

    // Higher order should give better attenuation (lower amplitude)
    assert!(
        attenuations[1] < attenuations[0],
        "4th order ({}) should attenuate more than 2nd order ({})",
        attenuations[1],
        attenuations[0]
    );

    assert!(
        attenuations[2] < attenuations[1],
        "6th order ({}) should attenuate more than 4th order ({})",
        attenuations[2],
        attenuations[1]
    );
}

#[test]
fn test_bandstop_narrow_notch() {
    // Test a narrow notch filter (60 Hz line noise rejection)
    let mut flt = ButterworthBandstop::new([58.0, 62.0], 4);
    let dt = 0.0001;

    // 60 Hz should be strongly rejected
    let amplitude_60hz = measure_steady_state_amplitude(&mut flt, 60.0, dt, 2000, 3000);
    assert!(
        amplitude_60hz < 0.2,
        "60 Hz not sufficiently attenuated: {}",
        amplitude_60hz
    );

    // 50 Hz should pass (outside notch)
    flt.reset();
    let amplitude_50hz = measure_steady_state_amplitude(&mut flt, 50.0, dt, 2000, 3000);
    assert!(
        amplitude_50hz > 0.6,
        "50 Hz excessively attenuated: {}",
        amplitude_50hz
    );

    // 70 Hz should pass (outside notch)
    flt.reset();
    let amplitude_70hz = measure_steady_state_amplitude(&mut flt, 70.0, dt, 2000, 3000);
    assert!(
        amplitude_70hz > 0.6,
        "70 Hz excessively attenuated: {}",
        amplitude_70hz
    );
}

#[test]
fn test_bandstop_composite_signal() {
    // Test with a composite signal: sum of frequencies inside and outside stopband
    let mut flt = ButterworthBandstop::new([45.0, 55.0], 4);
    let dt = 0.0001;

    let f_pass = 20.0; // In passband
    let f_stop = 50.0; // In stopband
    let omega_pass = 2.0 * PI * f_pass;
    let omega_stop = 2.0 * PI * f_stop;

    let mut max_output: f64 = 0.0;
    let mut input_50hz_component: f64 = 0.0;

    for i in 0..5000 {
        let t = i as f64 * dt;

        // Composite input: equal amplitude sinusoids
        let input = (omega_pass * t).sin() + (omega_stop * t).sin();
        flt.set_input(0, input);
        flt.update(t);
        flt.step(t, dt);

        if i > 2000 {
            max_output = max_output.max(flt.get_output(0).abs());

            // Track what the 50 Hz component would be without filtering
            input_50hz_component = input_50hz_component.max((omega_stop * t).sin().abs());
        }
    }

    // Output should be dominated by 20 Hz component (~1.0 amplitude)
    // and the 50 Hz component should be largely removed
    // Since we have two equal-amplitude inputs, without filtering max would be ~2.0
    // With filtering, 50 Hz is removed, so max should be closer to 1.0-1.5
    assert!(
        max_output < 1.8,
        "Composite signal: stopband component not sufficiently removed: {}",
        max_output
    );
    assert!(
        max_output > 0.7,
        "Composite signal: passband component excessively attenuated: {}",
        max_output
    );
}

#[test]
fn test_bandstop_reset() {
    let mut flt = ButterworthBandstop::new([45.0, 55.0], 2);

    // Run for a while
    for i in 0..100 {
        let t = i as f64 * 0.001;
        flt.set_input(0, (100.0 * t).sin());
        flt.update(t);
        flt.step(t, 0.001);
    }

    // State should be non-zero
    assert!(flt.state().iter().any(|&x| x != 0.0));

    // Reset
    flt.reset();

    // All states should be zero
    assert!(flt.state().iter().all(|&x| x == 0.0));
    assert_eq!(flt.get_output(0), 0.0);
}

#[test]
fn test_bandstop_buffer_revert() {
    let mut flt = ButterworthBandstop::new([45.0, 55.0], 2);

    // Run and buffer
    flt.set_input(0, 1.0);
    flt.update(0.0);
    flt.step(0.0, 0.01);
    let state_after_step: Vec<f64> = flt.state().to_vec();
    flt.buffer();

    // Continue running
    for _ in 0..10 {
        flt.update(0.0);
        flt.step(0.0, 0.01);
    }

    // State should have changed
    assert!(flt.state() != &state_after_step[..]);

    // Revert
    flt.revert();

    // Should match buffered state
    for (i, &val) in state_after_step.iter().enumerate() {
        assert_eq!(flt.state()[i], val);
    }
}

#[test]
fn test_bandstop_wide_stopband() {
    // Test a wide stopband filter
    let mut flt = ButterworthBandstop::new([20.0, 100.0], 4);
    let dt = 0.0001;

    // Frequencies in the middle of the stopband should be well attenuated
    let stopband_freqs_center = vec![35.0, 45.0, 55.0, 65.0];
    for freq in stopband_freqs_center {
        flt.reset();
        let amplitude = measure_steady_state_amplitude(&mut flt, freq, dt, 2000, 3000);
        assert!(
            amplitude < 0.4,
            "Wide stopband: {} Hz not sufficiently attenuated: {}",
            freq,
            amplitude
        );
    }

    // Frequencies near the edges have less attenuation
    let stopband_freqs_edge = vec![25.0, 85.0];
    for freq in stopband_freqs_edge {
        flt.reset();
        let amplitude = measure_steady_state_amplitude(&mut flt, freq, dt, 2000, 3000);
        assert!(
            amplitude < 0.7,
            "Wide stopband edge: {} Hz not sufficiently attenuated: {}",
            freq,
            amplitude
        );
    }

    // Frequencies outside should pass
    let passband_freqs = vec![5.0, 10.0, 150.0, 200.0];
    for freq in passband_freqs {
        flt.reset();
        let amplitude = measure_steady_state_amplitude(&mut flt, freq, dt, 2000, 3000);
        assert!(
            amplitude > 0.5,
            "Wide stopband: {} Hz passband excessively attenuated: {}",
            freq,
            amplitude
        );
    }
}

#[test]
fn test_bandstop_step_response() {
    // Test step response - should eventually settle to step value
    let mut flt = ButterworthBandstop::new([40.0, 60.0], 2);
    let dt = 0.001;

    // Apply step input
    for _ in 0..3000 {
        flt.set_input(0, 1.0);
        flt.update(0.0);
        flt.step(0.0, dt);
    }

    // Should settle to approximately 1.0 (DC passes through)
    assert!(
        (flt.get_output(0) - 1.0).abs() < 0.1,
        "Step response did not settle correctly: {}",
        flt.get_output(0)
    );
}

#[test]
#[should_panic(expected = "Low cutoff must be less than high cutoff")]
fn test_bandstop_invalid_cutoff_order() {
    ButterworthBandstop::new([100.0, 50.0], 2);
}

#[test]
#[should_panic(expected = "Filter order must be positive")]
fn test_bandstop_zero_order() {
    ButterworthBandstop::new([45.0, 55.0], 0);
}
