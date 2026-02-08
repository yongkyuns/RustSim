//! Comprehensive tests for TransferFunction block
//!
//! Tests various transfer function configurations and validates against analytical solutions

use rustsim::blocks::TransferFunction;
use rustsim::Block;

const TOL: f64 = 1e-3; // Tolerance for numerical integration errors

/// Test 1: First-order low-pass filter: H(s) = 1/(s+1)
#[test]
fn test_tf_01_first_order_lowpass() {
    // H(s) = 1/(s+1), time constant tau = 1
    let mut tf = TransferFunction::<1>::new(&[1.0], &[1.0, 1.0]);

    // Step response: y(t) = 1 - e^(-t)
    tf.set_input(0, 1.0);
    tf.update(0.0);

    let dt = 0.001;
    let mut t = 0.0;

    // Test at specific time points
    let test_points = [
        (1.0, 1.0 - (-1.0_f64).exp()),  // t = 1 tau
        (2.0, 1.0 - (-2.0_f64).exp()),  // t = 2 tau
        (3.0, 1.0 - (-3.0_f64).exp()),  // t = 3 tau
        (5.0, 1.0 - (-5.0_f64).exp()),  // t = 5 tau (99.3%)
    ];

    for (t_target, y_expected) in test_points {
        while t < t_target {
            tf.step(t, dt);
            t += dt;
        }
        let y_actual = tf.get_output(0);
        assert!(
            (y_actual - y_expected).abs() < TOL,
            "At t={}: expected {}, got {}, error={}",
            t_target,
            y_expected,
            y_actual,
            (y_actual - y_expected).abs()
        );
    }
}

/// Test 2: First-order high-pass filter: H(s) = s/(s+1)
#[test]
fn test_tf_02_first_order_highpass() {
    // H(s) = s/(s+1) = 1 - 1/(s+1)
    let mut tf = TransferFunction::<1>::new(&[1.0, 0.0], &[1.0, 1.0]);

    // Step response: y(t) = e^(-t)
    // This has direct feedthrough D=1
    tf.set_input(0, 1.0);
    tf.update(0.0);

    // Initial output should be 1 (direct feedthrough)
    assert!((tf.get_output(0) - 1.0).abs() < TOL);

    let dt = 0.01;
    let mut t = 0.0;

    // Test at t = 1, 2, 3
    for t_target in [1.0, 2.0, 3.0] {
        while t < t_target {
            tf.step(t, dt);
            t += dt;
        }
        let y_expected = (-t).exp();
        let y_actual = tf.get_output(0);
        assert!(
            (y_actual - y_expected).abs() < TOL,
            "At t={}: expected {}, got {}",
            t_target,
            y_expected,
            y_actual
        );
    }
}

/// Test 3: Second-order underdamped: H(s) = ω_n^2 / (s^2 + 2ζω_n s + ω_n^2)
#[test]
fn test_tf_03_second_order_underdamped() {
    // Natural frequency ωn = 1, damping ratio ζ = 0.2 (underdamped)
    let wn = 1.0;
    let zeta = 0.2;

    // H(s) = 1 / (s^2 + 0.4s + 1)
    let mut tf = TransferFunction::<2>::new(&[1.0], &[1.0, 2.0 * zeta * wn, wn * wn]);

    // Step response for underdamped second-order system:
    // y(t) = 1 - e^(-ζω_n t) * (cos(ω_d t) + (ζ/√(1-ζ²)) sin(ω_d t))
    // where ω_d = ω_n √(1-ζ²) is the damped natural frequency

    tf.set_input(0, 1.0);
    tf.update(0.0);

    let dt = 0.01;
    let mut t = 0.0;

    let wd = wn * (1.0 - zeta * zeta).sqrt(); // damped frequency

    let test_times = [1.0, 5.0, 10.0, 20.0];
    for &t_target in &test_times {
        while t < t_target {
            tf.step(t, dt);
            t += dt;
        }

        let decay = (-zeta * wn * t).exp();
        let y_expected =
            1.0 - decay * ((wd * t).cos() + (zeta / (1.0 - zeta * zeta).sqrt()) * (wd * t).sin());
        let y_actual = tf.get_output(0);

        assert!(
            (y_actual - y_expected).abs() < 0.01, // Slightly relaxed tolerance for 2nd order
            "At t={}: expected {}, got {}",
            t_target,
            y_expected,
            y_actual
        );
    }
}

/// Test 4: Second-order critically damped: H(s) = 1 / (s+1)^2
#[test]
fn test_tf_04_second_order_critically_damped() {
    // H(s) = 1 / (s^2 + 2s + 1) = 1/(s+1)^2
    let mut tf = TransferFunction::<2>::new(&[1.0], &[1.0, 2.0, 1.0]);

    // Step response: y(t) = 1 - (1 + t) e^(-t)
    tf.set_input(0, 1.0);
    tf.update(0.0);

    let dt = 0.01;
    let mut t = 0.0;

    for t_target in [1.0, 2.0, 3.0, 5.0] {
        while t < t_target {
            tf.step(t, dt);
            t += dt;
        }

        let y_expected = 1.0 - (1.0 + t) * (-t).exp();
        let y_actual = tf.get_output(0);

        assert!(
            (y_actual - y_expected).abs() < 0.01,
            "At t={}: expected {}, got {}",
            t_target,
            y_expected,
            y_actual
        );
    }
}

/// Test 5: Second-order overdamped: H(s) = 1 / (s^2 + 3s + 2) = 1/((s+1)(s+2))
#[test]
fn test_tf_05_second_order_overdamped() {
    // H(s) = 1 / ((s+1)(s+2))
    let mut tf = TransferFunction::<2>::new(&[1.0], &[1.0, 3.0, 2.0]);

    // Partial fractions: H(s) = 1/(s+1) - 1/(s+2)
    // Step response: Y(s) = H(s)/s = 1/(s(s+1)(s+2))
    // By partial fractions: Y(s) = 1/2 * 1/s - 1/(s+1) + 1/2 * 1/(s+2)
    // Inverse Laplace: y(t) = 1/2 - e^(-t) + (1/2)e^(-2t)

    tf.set_input(0, 1.0);
    tf.update(0.0);

    let dt = 0.01;
    let mut t = 0.0;

    for t_target in [0.5, 1.0, 2.0, 5.0] {
        while t < t_target {
            tf.step(t, dt);
            t += dt;
        }

        // Correct step response for 1/((s+1)(s+2))
        let y_expected = 0.5 - (-t).exp() + 0.5 * (-2.0 * t).exp();
        let y_actual = tf.get_output(0);

        assert!(
            (y_actual - y_expected).abs() < 0.01,
            "At t={}: expected {}, got {}",
            t_target,
            y_expected,
            y_actual
        );
    }
}

/// Test 6: Lead compensator: H(s) = (s+2)/(s+10)
#[test]
fn test_tf_06_lead_compensator() {
    // H(s) = (s+2)/(s+10)
    let mut tf = TransferFunction::<1>::new(&[1.0, 2.0], &[1.0, 10.0]);

    // This has direct feedthrough since num and den have same degree
    assert!(tf.has_passthrough());

    // Rewrite: H(s) = 1 - 8/(s+10)
    // Step response: y(t) = 1 - 0.8(1 - e^(-10t)) = 0.2 + 0.8 e^(-10t)

    tf.set_input(0, 1.0);
    tf.update(0.0);

    let dt = 0.001;
    let mut t = 0.0;

    for t_target in [0.1, 0.2, 0.5, 1.0] {
        while t < t_target {
            tf.step(t, dt);
            t += dt;
        }

        let y_expected = 0.2 + 0.8 * (-10.0 * t).exp();
        let y_actual = tf.get_output(0);

        assert!(
            (y_actual - y_expected).abs() < 0.01,
            "At t={}: expected {}, got {}",
            t_target,
            y_expected,
            y_actual
        );
    }
}

/// Test 7: Lag compensator: H(s) = (s+1)/(s+0.1)
#[test]
fn test_tf_07_lag_compensator() {
    // H(s) = (s+1)/(s+0.1)
    let mut tf = TransferFunction::<1>::new(&[1.0, 1.0], &[1.0, 0.1]);

    assert!(tf.has_passthrough());

    tf.set_input(0, 1.0);
    tf.update(0.0);

    let dt = 0.01;
    let mut t = 0.0;

    // Steady state should be 1/0.1 = 10
    while t < 50.0 {
        tf.step(t, dt);
        t += dt;
    }

    assert!(
        (tf.get_output(0) - 10.0).abs() < 0.1,
        "Steady state expected 10.0, got {}",
        tf.get_output(0)
    );
}

/// Test 8: Pure integrator: H(s) = 1/s
#[test]
fn test_tf_08_pure_integrator() {
    // H(s) = 1/s
    let mut tf = TransferFunction::<1>::new(&[1.0], &[1.0, 0.0]);

    // Step response: y(t) = t
    tf.set_input(0, 1.0);
    tf.update(0.0);

    let dt = 0.01;
    let mut t = 0.0;

    for t_target in [1.0, 2.0, 5.0] {
        while t < t_target {
            tf.step(t, dt);
            t += dt;
        }

        let y_expected = t;
        let y_actual = tf.get_output(0);

        assert!(
            (y_actual - y_expected).abs() < 0.05,
            "At t={}: expected {}, got {}",
            t_target,
            y_expected,
            y_actual
        );
    }
}

/// Test 9: Double integrator: H(s) = 1/s^2
#[test]
fn test_tf_09_double_integrator() {
    // H(s) = 1/s^2
    let mut tf = TransferFunction::<2>::new(&[1.0], &[1.0, 0.0, 0.0]);

    // Step response: y(t) = t^2/2
    tf.set_input(0, 1.0);
    tf.update(0.0);

    let dt = 0.01;
    let mut t = 0.0;

    for t_target in [1.0, 2.0, 3.0] {
        while t < t_target {
            tf.step(t, dt);
            t += dt;
        }

        let y_expected = t * t / 2.0;
        let y_actual = tf.get_output(0);

        assert!(
            (y_actual - y_expected).abs() < 0.05,
            "At t={}: expected {}, got {}",
            t_target,
            y_expected,
            y_actual
        );
    }
}

/// Test 10: Notch filter: H(s) = (s^2 + 1) / (s^2 + 0.1s + 1)
#[test]
fn test_tf_10_notch_filter() {
    // H(s) = (s^2 + 1) / (s^2 + 0.1s + 1)
    // This is a notch at ω = 1 rad/s with light damping
    let mut tf = TransferFunction::<2>::new(&[1.0, 0.0, 1.0], &[1.0, 0.1, 1.0]);

    // Has direct feedthrough (same degree)
    assert!(tf.has_passthrough());

    tf.set_input(0, 1.0);
    tf.update(0.0);

    let dt = 0.001;

    // Run for a while
    for _ in 0..10000 {
        tf.step(0.0, dt);
    }

    // Steady state: H(0) = 1/1 = 1
    assert!(
        (tf.get_output(0) - 1.0).abs() < 0.05,
        "Steady state expected 1.0, got {}",
        tf.get_output(0)
    );
}

/// Test 11: Third-order system
#[test]
fn test_tf_11_third_order() {
    // H(s) = 1 / (s^3 + 6s^2 + 11s + 6) = 1/((s+1)(s+2)(s+3))
    let mut tf = TransferFunction::<3>::new(&[1.0], &[1.0, 6.0, 11.0, 6.0]);

    tf.set_input(0, 1.0);
    tf.update(0.0);

    let dt = 0.01;

    // Run to steady state
    for _ in 0..1000 {
        tf.step(0.0, dt);
    }

    // DC gain: H(0) = 1/6
    assert!(
        (tf.get_output(0) - 1.0 / 6.0).abs() < 0.01,
        "Steady state expected {}, got {}",
        1.0 / 6.0,
        tf.get_output(0)
    );
}

/// Test 12: Fourth-order Butterworth filter
#[test]
fn test_tf_12_fourth_order_butterworth() {
    // 4th order Butterworth with ωc = 1 rad/s
    // H(s) = 1 / (s^4 + 2.613s^3 + 3.414s^2 + 2.613s + 1)
    let mut tf = TransferFunction::<4>::new(&[1.0], &[1.0, 2.613, 3.414, 2.613, 1.0]);

    tf.set_input(0, 1.0);
    tf.update(0.0);

    let dt = 0.01;

    // Run to steady state
    for _ in 0..2000 {
        tf.step(0.0, dt);
    }

    // DC gain: H(0) = 1/1 = 1
    assert!(
        (tf.get_output(0) - 1.0).abs() < 0.05,
        "Steady state expected 1.0, got {}",
        tf.get_output(0)
    );
}

/// Test 13: Resonant system: H(s) = 1/(s^2 + 0.1s + 1)
#[test]
fn test_tf_13_resonant_system() {
    // Very lightly damped second-order (ζ = 0.05)
    let wn = 1.0;
    let zeta = 0.05;

    let mut tf = TransferFunction::<2>::new(&[1.0], &[1.0, 2.0 * zeta * wn, wn * wn]);

    tf.set_input(0, 1.0);
    tf.update(0.0);

    let dt = 0.001;

    // Should show significant overshoot
    let mut max_output: f64 = 0.0;
    for _ in 0..5000 {
        tf.step(0.0, dt);
        max_output = max_output.max(tf.get_output(0));
    }

    // Overshoot for ζ=0.05 should be approximately 85%
    // Peak value ≈ 1 + e^(-πζ/√(1-ζ²)) ≈ 1.85
    assert!(
        max_output > 1.5,
        "Expected significant overshoot, got max = {}",
        max_output
    );

    // Eventually oscillates around 1.0 (very slow decay due to low damping)
    // After a long time, should be closer to 1.0
    for _ in 0..50000 {
        tf.step(0.0, dt);
    }
    // With very light damping (ζ=0.05), settling takes a very long time
    // The system oscillates with slow decay - just verify it's in reasonable range
    assert!(
        (tf.get_output(0) - 1.0).abs() < 0.3,
        "Should be oscillating around 1.0 but got {}",
        tf.get_output(0)
    );
}

/// Test 14: Impulse response approximation (narrow pulse)
#[test]
fn test_tf_14_impulse_response() {
    // H(s) = 1/(s+1)
    let mut tf = TransferFunction::<1>::new(&[1.0], &[1.0, 1.0]);

    // Approximate impulse with narrow pulse: area = 1, duration = 0.01
    let dt = 0.001;
    let pulse_duration = 0.01;
    let pulse_amplitude = 1.0 / pulse_duration;

    tf.set_input(0, pulse_amplitude);
    tf.update(0.0);

    let mut t = 0.0;
    for _ in 0..10 {
        // 10 steps = 0.01s
        tf.step(t, dt);
        t += dt;
    }

    // Switch to zero input
    tf.set_input(0, 0.0);

    // Impulse response of 1/(s+1) is e^(-t)
    // After pulse at t ≈ 0.01, output should be close to e^(-t)
    for t_target in [0.5, 1.0, 2.0] {
        while t < t_target {
            tf.step(t, dt);
            t += dt;
        }

        let y_expected = (-t).exp();
        let y_actual = tf.get_output(0);

        assert!(
            (y_actual - y_expected).abs() < 0.1, // Relaxed tolerance for impulse approximation
            "At t={}: expected {}, got {}",
            t_target,
            y_expected,
            y_actual
        );
    }
}

/// Test 15: Frequency response at DC (steady state gain)
#[test]
fn test_tf_15_dc_gain() {
    let test_cases = vec![
        // (num, den, expected_dc_gain)
        (vec![1.0], vec![1.0, 1.0], 1.0), // 1/(s+1) -> DC gain = 1
        (vec![2.0], vec![1.0, 1.0], 2.0), // 2/(s+1) -> DC gain = 2
        (vec![1.0], vec![1.0, 2.0, 1.0], 1.0), // 1/(s^2+2s+1) -> DC gain = 1
        (vec![5.0], vec![1.0, 2.0, 1.0], 5.0), // 5/(s^2+2s+1) -> DC gain = 5
        (vec![1.0, 2.0], vec![1.0, 4.0], 0.5), // (s+2)/(s+4) -> DC gain = 2/4 = 0.5
    ];

    for (num, den, expected_gain) in test_cases {
        let order = den.len() - 1;

        match order {
            1 => {
                let mut tf = TransferFunction::<1>::new(&num, &den);
                tf.set_input(0, 1.0);
                tf.update(0.0);
                let dt = 0.01;
                for _ in 0..1000 {
                    tf.step(0.0, dt);
                }
                assert!(
                    (tf.get_output(0) - expected_gain).abs() < 0.01,
                    "DC gain test failed for {:?}/{:?}: expected {}, got {}",
                    num,
                    den,
                    expected_gain,
                    tf.get_output(0)
                );
            }
            2 => {
                let mut tf = TransferFunction::<2>::new(&num, &den);
                tf.set_input(0, 1.0);
                tf.update(0.0);
                let dt = 0.01;
                for _ in 0..1000 {
                    tf.step(0.0, dt);
                }
                assert!(
                    (tf.get_output(0) - expected_gain).abs() < 0.01,
                    "DC gain test failed for {:?}/{:?}: expected {}, got {}",
                    num,
                    den,
                    expected_gain,
                    tf.get_output(0)
                );
            }
            _ => {}
        }
    }
}
