//! Bouncing ball system evaluation tests
//!
//! Tests event detection and handling on a bouncing ball system:
//! dy/dt = v
//! dv/dt = -g
//!
//! Event: y = 0 (ground collision) -> v = -e*v (coefficient of restitution e)
//!
//! This tests the ability to detect zero crossings and handle
//! discontinuous state changes (velocity reversal)

use approx::assert_relative_eq;
use nalgebra::DVector;
use rustsim::solvers::*;

const G: f64 = 9.81; // gravitational acceleration

/// Compute theoretical bounce times for a bouncing ball
/// with initial height h and coefficient of restitution e=1 (elastic)
fn bounce_times_analytical(h: f64, n_bounces: usize) -> Vec<f64> {
    let t_fall = (2.0 * h / G).sqrt();
    let mut times = Vec::new();

    for i in 0..n_bounces {
        times.push(t_fall + (i as f64) * 2.0 * t_fall);
    }

    times
}

#[test]
fn test_bouncingball_event_detection_rkf21() {
    let h = 1.0; // initial height
    let e = 1.0; // coefficient of restitution (elastic)

    let mut state = DVector::from_vec(vec![h, 0.0]); // [position, velocity]
    let mut solver = RKF21::new(state.clone());

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            x[1],      // dy/dt = v
            -G,        // dv/dt = -g
        ])
    };

    let mut t = 0.0;
    let mut dt = 0.001;
    let t_final = 10.0;
    let mut bounce_times = Vec::new();

    while t < t_final {
        solver.buffer(dt);

        // RKF21 has 3 stages
        let mut result = SolverStepResult::default();
        for _ in 0..3 {
            result = solver.step(&f, dt);
        }

        if result.success {
            // Check for ground collision (y <= 0)
            if solver.state()[0] <= 0.0 {
                // Record bounce time
                bounce_times.push(t);

                // Reverse velocity (elastic collision)
                let mut new_state = solver.state().clone();
                new_state[0] = new_state[0].abs();
                new_state[1] = -e * new_state[1];

                // Reset solver with new state
                solver = RKF21::new(new_state);
            }

            t += dt;
            if let Some(scale) = result.scale {
                dt = (dt * scale.min(2.0).max(0.5)).min(0.01);
            }
        } else {
            solver.revert().unwrap();
            dt *= 0.5;
        }
    }

    // Check that we detected the bounces at approximately the right times
    let expected_times = bounce_times_analytical(h, 5);

    assert!(
        bounce_times.len() >= 5,
        "Should detect at least 5 bounces, got {}",
        bounce_times.len()
    );

    for (i, &expected) in expected_times.iter().enumerate() {
        let actual = bounce_times[i];
        let error = (actual - expected).abs();
        assert!(
            error < 0.01,
            "Bounce {} time error: expected {}, got {}, error {}",
            i,
            expected,
            actual,
            error
        );
    }
}

#[test]
fn test_bouncingball_event_detection_rkbs32() {
    let h = 1.0;
    let e = 1.0;

    let mut state = DVector::from_vec(vec![h, 0.0]);
    let mut solver = RKBS32::new(state.clone());

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![x[1], -G])
    };

    let mut t = 0.0;
    let mut dt = 0.001;
    let t_final = 10.0;
    let mut bounce_times = Vec::new();

    while t < t_final {
        solver.buffer(dt);

        // RKBS32 has 4 stages
        let mut result = SolverStepResult::default();
        for _ in 0..4 {
            result = solver.step(&f, dt);
        }

        if result.success {
            if solver.state()[0] <= 0.0 {
                bounce_times.push(t);
                let mut new_state = solver.state().clone();
                new_state[0] = new_state[0].abs();
                new_state[1] = -e * new_state[1];
                solver = RKBS32::new(new_state);
            }

            t += dt;
            if let Some(scale) = result.scale {
                dt = (dt * scale.min(2.0).max(0.5)).min(0.01);
            }
        } else {
            solver.revert().unwrap();
            dt *= 0.5;
        }
    }

    let expected_times = bounce_times_analytical(h, 5);
    assert!(bounce_times.len() >= 5);

    for (i, &expected) in expected_times.iter().enumerate() {
        let actual = bounce_times[i];
        assert!(
            (actual - expected).abs() < 0.01,
            "Bounce {} time error too large",
            i
        );
    }
}

#[test]
fn test_bouncingball_event_detection_rkf45() {
    let h = 1.0;
    let e = 1.0;

    let mut state = DVector::from_vec(vec![h, 0.0]);
    let mut solver = RKF45::new(state.clone());

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![x[1], -G])
    };

    let mut t = 0.0;
    let mut dt = 0.001;
    let t_final = 10.0;
    let mut bounce_times = Vec::new();

    while t < t_final {
        solver.buffer(dt);

        // RKF45 has 6 stages
        let mut result = SolverStepResult::default();
        for _ in 0..6 {
            result = solver.step(&f, dt);
        }

        if result.success {
            if solver.state()[0] <= 0.0 {
                bounce_times.push(t);
                let mut new_state = solver.state().clone();
                new_state[0] = new_state[0].abs();
                new_state[1] = -e * new_state[1];
                solver = RKF45::new(new_state);
            }

            t += dt;
            if let Some(scale) = result.scale {
                dt = (dt * scale.min(2.0).max(0.5)).min(0.01);
            }
        } else {
            solver.revert().unwrap();
            dt *= 0.5;
        }
    }

    let expected_times = bounce_times_analytical(h, 5);
    assert!(bounce_times.len() >= 5);

    for (i, &expected) in expected_times.iter().enumerate() {
        let actual = bounce_times[i];
        assert!(
            (actual - expected).abs() < 0.005,
            "Bounce {} time error too large",
            i
        );
    }
}

#[test]
fn test_bouncingball_event_detection_rkck54() {
    let h = 1.0;
    let e = 1.0;

    let mut state = DVector::from_vec(vec![h, 0.0]);
    let mut solver = RKCK54::new(state.clone());

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![x[1], -G])
    };

    let mut t = 0.0;
    let mut dt = 0.001;
    let t_final = 10.0;
    let mut bounce_times = Vec::new();

    while t < t_final {
        solver.buffer(dt);

        // RKCK54 has 6 stages
        let mut result = SolverStepResult::default();
        for _ in 0..6 {
            result = solver.step(&f, dt);
        }

        if result.success {
            if solver.state()[0] <= 0.0 {
                bounce_times.push(t);
                let mut new_state = solver.state().clone();
                new_state[0] = new_state[0].abs();
                new_state[1] = -e * new_state[1];
                solver = RKCK54::new(new_state);
            }

            t += dt;
            if let Some(scale) = result.scale {
                dt = (dt * scale.min(2.0).max(0.5)).min(0.01);
            }
        } else {
            solver.revert().unwrap();
            dt *= 0.5;
        }
    }

    let expected_times = bounce_times_analytical(h, 5);
    assert!(bounce_times.len() >= 5);

    for (i, &expected) in expected_times.iter().enumerate() {
        let actual = bounce_times[i];
        assert!(
            (actual - expected).abs() < 0.005,
            "Bounce {} time error too large",
            i
        );
    }
}

#[test]
fn test_bouncingball_event_detection_rkdp54() {
    let h = 1.0;
    let e = 1.0;

    let mut state = DVector::from_vec(vec![h, 0.0]);
    let mut solver = RKDP54::new(state.clone());

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![x[1], -G])
    };

    let mut t = 0.0;
    let mut dt = 0.001;
    let t_final = 10.0;
    let mut bounce_times = Vec::new();

    while t < t_final {
        solver.buffer(dt);

        // RKDP54 has 7 stages
        let mut result = SolverStepResult::default();
        for _ in 0..7 {
            result = solver.step(&f, dt);
        }

        if result.success {
            if solver.state()[0] <= 0.0 {
                bounce_times.push(t);
                let mut new_state = solver.state().clone();
                new_state[0] = new_state[0].abs();
                new_state[1] = -e * new_state[1];
                solver = RKDP54::new(new_state);
            }

            t += dt;
            if let Some(scale) = result.scale {
                dt = (dt * scale.min(2.0).max(0.5)).min(0.01);
            }
        } else {
            solver.revert().unwrap();
            dt *= 0.5;
        }
    }

    let expected_times = bounce_times_analytical(h, 5);
    assert!(bounce_times.len() >= 5);

    for (i, &expected) in expected_times.iter().enumerate() {
        let actual = bounce_times[i];
        assert!(
            (actual - expected).abs() < 0.002,
            "Bounce {} time error too large",
            i
        );
    }
}

#[test]
fn test_bouncingball_event_detection_rkv65() {
    let h = 1.0;
    let e = 1.0;

    let mut state = DVector::from_vec(vec![h, 0.0]);
    let mut solver = RKV65::new(state.clone());

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![x[1], -G])
    };

    let mut t = 0.0;
    let mut dt = 0.001;
    let t_final = 10.0;
    let mut bounce_times = Vec::new();

    while t < t_final {
        solver.buffer(dt);

        // RKV65 has 8 stages
        let mut result = SolverStepResult::default();
        for _ in 0..8 {
            result = solver.step(&f, dt);
        }

        if result.success {
            if solver.state()[0] <= 0.0 {
                bounce_times.push(t);
                let mut new_state = solver.state().clone();
                new_state[0] = new_state[0].abs();
                new_state[1] = -e * new_state[1];
                solver = RKV65::new(new_state);
            }

            t += dt;
            if let Some(scale) = result.scale {
                dt = (dt * scale.min(2.0).max(0.5)).min(0.01);
            }
        } else {
            solver.revert().unwrap();
            dt *= 0.5;
        }
    }

    let expected_times = bounce_times_analytical(h, 5);
    assert!(bounce_times.len() >= 5);

    for (i, &expected) in expected_times.iter().enumerate() {
        let actual = bounce_times[i];
        assert!(
            (actual - expected).abs() < 0.002,
            "Bounce {} time error too large",
            i
        );
    }
}

#[test]
fn test_bouncingball_event_detection_rkdp87() {
    let h = 1.0;
    let e = 1.0;

    let mut state = DVector::from_vec(vec![h, 0.0]);
    let mut solver = RKDP87::new(state.clone());

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![x[1], -G])
    };

    let mut t = 0.0;
    let mut dt = 0.001;
    let t_final = 10.0;
    let mut bounce_times = Vec::new();

    while t < t_final {
        solver.buffer(dt);

        // RKDP87 has 13 stages
        let mut result = SolverStepResult::default();
        for _ in 0..13 {
            result = solver.step(&f, dt);
        }

        if result.success {
            if solver.state()[0] <= 0.0 {
                bounce_times.push(t);
                let mut new_state = solver.state().clone();
                new_state[0] = new_state[0].abs();
                new_state[1] = -e * new_state[1];
                solver = RKDP87::new(new_state);
            }

            t += dt;
            if let Some(scale) = result.scale {
                dt = (dt * scale.min(2.0).max(0.5)).min(0.01);
            }
        } else {
            solver.revert().unwrap();
            dt *= 0.5;
        }
    }

    let expected_times = bounce_times_analytical(h, 5);
    assert!(bounce_times.len() >= 5);

    for (i, &expected) in expected_times.iter().enumerate() {
        let actual = bounce_times[i];
        assert!(
            (actual - expected).abs() < 0.001,
            "Bounce {} time error too large",
            i
        );
    }
}

#[test]
fn test_bouncingball_energy_decay() {
    // Test inelastic collision (e < 1)
    let h = 1.0;
    let e = 0.8; // coefficient of restitution < 1

    let mut state = DVector::from_vec(vec![h, 0.0]);
    let mut solver = RKDP54::new(state.clone());

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![x[1], -G])
    };

    let mut t = 0.0;
    let mut dt = 0.001;
    let t_final = 10.0;
    let mut max_heights = Vec::new();

    let mut last_y = h;
    let mut max_y = h;

    while t < t_final {
        solver.buffer(dt);

        // RKDP54 has 7 stages
        let mut result = SolverStepResult::default();
        for _ in 0..7 {
            result = solver.step(&f, dt);
        }

        if result.success {
            let y = solver.state()[0];

            // Track maximum height
            if y > max_y {
                max_y = y;
            }

            // Detect bounce
            if solver.state()[0] <= 0.0 {
                max_heights.push(max_y);
                max_y = 0.0;

                let mut new_state = solver.state().clone();
                new_state[0] = new_state[0].abs();
                new_state[1] = -e * new_state[1];
                solver = RKDP54::new(new_state);
            }

            t += dt;
            if let Some(scale) = result.scale {
                dt = (dt * scale.min(2.0).max(0.5)).min(0.01);
            }

            last_y = y;
        } else {
            solver.revert().unwrap();
            dt *= 0.5;
        }
    }

    // Verify that maximum heights decay exponentially
    assert!(max_heights.len() >= 5, "Should have at least 5 bounces");

    for i in 1..max_heights.len().min(5) {
        let ratio = max_heights[i] / max_heights[i - 1];
        // For coefficient of restitution e, the ratio of successive max heights
        // should be approximately e^2
        let expected_ratio = e * e;
        let error = (ratio - expected_ratio).abs();

        assert!(
            error < 0.1,
            "Height ratio {} should be close to e^2={}, got ratio={}, error={}",
            i,
            expected_ratio,
            ratio,
            error
        );
    }
}
