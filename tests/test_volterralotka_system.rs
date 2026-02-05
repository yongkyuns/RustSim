//! Volterra-Lotka (Predator-Prey) system evaluation tests
//!
//! Tests solvers on the classic predator-prey system:
//! dx/dt = αx - βxy
//! dy/dt = δxy - γy
//!
//! Parameters: α=1.0, β=0.1, δ=0.5, γ=1.2
//! Initial conditions: x=10 (predator), y=5 (prey)
//!
//! This system exhibits periodic oscillations and is conservative

use approx::assert_relative_eq;
use nalgebra::DVector;
use rustsim::solvers::*;

/// High-accuracy reference solution for Volterra-Lotka system
fn volterralotka_reference(
    t: f64,
    x0: &[f64; 2],
    alpha: f64,
    beta: f64,
    delta: f64,
    gamma: f64,
) -> [f64; 2] {
    let mut state = DVector::from_vec(vec![x0[0], x0[1]]);
    let mut solver = RKDP54::new(state.clone());

    let dt_base = 0.00005;
    let steps = (t / dt_base) as usize;

    for _ in 0..steps {
        solver.buffer(dt_base);

        let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
            DVector::from_vec(vec![
                alpha * x[0] - beta * x[0] * x[1],
                -gamma * x[1] + delta * x[0] * x[1],
            ])
        };

        // RKDP54 has 7 stages
        for _ in 0..7 {
            solver.step(&f, dt_base);
        }
    }

    [solver.state()[0], solver.state()[1]]
}

#[test]
fn test_volterralotka_rkf21() {
    let alpha = 1.0;
    let beta = 0.1;
    let delta = 0.5;
    let gamma = 1.2;
    let x0 = [10.0, 5.0];
    let t_final = 5.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            alpha * x[0] - beta * x[0] * x[1],
            -gamma * x[1] + delta * x[0] * x[1],
        ])
    };

    let tolerances = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = RKF21::new(initial);

        let mut t = 0.0;
        let mut dt = 0.001;

        while t < t_final {
            solver.buffer(dt);

            // RKF21 has 3 stages
            let mut result = SolverStepResult::default();
            for _ in 0..3 {
                result = solver.step(&f, dt);
            }

            if result.success {
                t += dt;
                if let Some(scale) = result.scale {
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = volterralotka_reference(t_final, &x0, alpha, beta, delta, gamma);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(error < tol * 1000.0, "RKF21 tol={} error={}", tol, error);
    }
}

#[test]
fn test_volterralotka_rkbs32() {
    let alpha = 1.0;
    let beta = 0.1;
    let delta = 0.5;
    let gamma = 1.2;
    let x0 = [10.0, 5.0];
    let t_final = 5.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            alpha * x[0] - beta * x[0] * x[1],
            -gamma * x[1] + delta * x[0] * x[1],
        ])
    };

    let tolerances = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = RKBS32::new(initial);

        let mut t = 0.0;
        let mut dt = 0.001;

        while t < t_final {
            solver.buffer(dt);

            // RKBS32 has 4 stages
            let mut result = SolverStepResult::default();
            for _ in 0..4 {
                result = solver.step(&f, dt);
            }

            if result.success {
                t += dt;
                if let Some(scale) = result.scale {
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = volterralotka_reference(t_final, &x0, alpha, beta, delta, gamma);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(error < tol * 100.0, "RKBS32 tol={} error={}", tol, error);
    }
}

#[test]
fn test_volterralotka_rkf45() {
    let alpha = 1.0;
    let beta = 0.1;
    let delta = 0.5;
    let gamma = 1.2;
    let x0 = [10.0, 5.0];
    let t_final = 5.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            alpha * x[0] - beta * x[0] * x[1],
            -gamma * x[1] + delta * x[0] * x[1],
        ])
    };

    let tolerances = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = RKF45::new(initial);

        let mut t = 0.0;
        let mut dt = 0.001;

        while t < t_final {
            solver.buffer(dt);

            // RKF45 has 6 stages
            let mut result = SolverStepResult::default();
            for _ in 0..6 {
                result = solver.step(&f, dt);
            }

            if result.success {
                t += dt;
                if let Some(scale) = result.scale {
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = volterralotka_reference(t_final, &x0, alpha, beta, delta, gamma);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(error < tol * 500.0, "RKF45 tol={} error={}", tol, error);
    }
}

#[test]
fn test_volterralotka_rkck54() {
    let alpha = 1.0;
    let beta = 0.1;
    let delta = 0.5;
    let gamma = 1.2;
    let x0 = [10.0, 5.0];
    let t_final = 5.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            alpha * x[0] - beta * x[0] * x[1],
            -gamma * x[1] + delta * x[0] * x[1],
        ])
    };

    let tolerances = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = RKCK54::new(initial);

        let mut t = 0.0;
        let mut dt = 0.001;

        while t < t_final {
            solver.buffer(dt);

            // RKCK54 has 6 stages
            let mut result = SolverStepResult::default();
            for _ in 0..6 {
                result = solver.step(&f, dt);
            }

            if result.success {
                t += dt;
                if let Some(scale) = result.scale {
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = volterralotka_reference(t_final, &x0, alpha, beta, delta, gamma);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(error < tol * 50.0, "RKCK54 tol={} error={}", tol, error);
    }
}

#[test]
fn test_volterralotka_rkdp54() {
    let alpha = 1.0;
    let beta = 0.1;
    let delta = 0.5;
    let gamma = 1.2;
    let x0 = [10.0, 5.0];
    let t_final = 5.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            alpha * x[0] - beta * x[0] * x[1],
            -gamma * x[1] + delta * x[0] * x[1],
        ])
    };

    let tolerances = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = RKDP54::new(initial);

        let mut t = 0.0;
        let mut dt = 0.001;

        while t < t_final {
            solver.buffer(dt);

            // RKDP54 has 7 stages
            let mut result = SolverStepResult::default();
            for _ in 0..7 {
                result = solver.step(&f, dt);
            }

            if result.success {
                t += dt;
                if let Some(scale) = result.scale {
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = volterralotka_reference(t_final, &x0, alpha, beta, delta, gamma);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(error < tol * 500.0, "RKDP54 tol={} error={}", tol, error);
    }
}

#[test]
fn test_volterralotka_rkv65() {
    let alpha = 1.0;
    let beta = 0.1;
    let delta = 0.5;
    let gamma = 1.2;
    let x0 = [10.0, 5.0];
    let t_final = 5.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            alpha * x[0] - beta * x[0] * x[1],
            -gamma * x[1] + delta * x[0] * x[1],
        ])
    };

    let tolerances = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = RKV65::new(initial);

        let mut t = 0.0;
        let mut dt = 0.001;

        while t < t_final {
            solver.buffer(dt);

            // RKV65 has 8 stages
            let mut result = SolverStepResult::default();
            for _ in 0..8 {
                result = solver.step(&f, dt);
            }

            if result.success {
                t += dt;
                if let Some(scale) = result.scale {
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = volterralotka_reference(t_final, &x0, alpha, beta, delta, gamma);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(error < tol * 20.0, "RKV65 tol={} error={}", tol, error);
    }
}

#[test]
fn test_volterralotka_rkdp87() {
    let alpha = 1.0;
    let beta = 0.1;
    let delta = 0.5;
    let gamma = 1.2;
    let x0 = [10.0, 5.0];
    let t_final = 5.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            alpha * x[0] - beta * x[0] * x[1],
            -gamma * x[1] + delta * x[0] * x[1],
        ])
    };

    let tolerances = [1e-4, 1e-5, 1e-6, 1e-7, 1e-8];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = RKDP87::new(initial);

        let mut t = 0.0;
        let mut dt = 0.001;

        while t < t_final {
            solver.buffer(dt);

            // RKDP87 has 13 stages
            let mut result = SolverStepResult::default();
            for _ in 0..13 {
                result = solver.step(&f, dt);
            }

            if result.success {
                t += dt;
                if let Some(scale) = result.scale {
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = volterralotka_reference(t_final, &x0, alpha, beta, delta, gamma);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(error < tol * 200.0, "RKDP87 tol={} error={}", tol, error);
    }
}
