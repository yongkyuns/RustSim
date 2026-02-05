//! Brusselator system evaluation tests
//!
//! Tests solvers on the Brusselator autocatalytic reaction system:
//! dx/dt = a - x - bx + x²y
//! dy/dt = bx - x²y
//!
//! Parameters: a=0.4, b=1.2
//! Initial conditions: x=0, y=0
//!
//! The Brusselator exhibits oscillatory behavior and limit cycles

use approx::assert_relative_eq;
use nalgebra::DVector;
use rustsim::solvers::*;

/// High-accuracy reference solution for Brusselator system
fn brusselator_reference(t: f64, x0: &[f64; 2], a: f64, b: f64) -> [f64; 2] {
    let mut state = DVector::from_vec(vec![x0[0], x0[1]]);
    let mut solver = RKDP54::new(state.clone());

    let dt_base = 0.0001;
    let steps = (t / dt_base) as usize;

    for _ in 0..steps {
        solver.buffer(dt_base);

        let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
            DVector::from_vec(vec![
                a - x[0] - b * x[0] + x[0] * x[0] * x[1],
                b * x[0] - x[0] * x[0] * x[1],
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
fn test_brusselator_rkbs32() {
    let a = 0.4;
    let b = 1.2;
    let x0 = [0.0, 0.0];
    let t_final = 50.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            a - x[0] - b * x[0] + x[0] * x[0] * x[1],
            b * x[0] - x[0] * x[0] * x[1],
        ])
    };

    let tolerances = [1e-5, 1e-6, 1e-7];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = RKBS32::new(initial);

        let mut t = 0.0;
        let mut dt = 0.01;

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
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.5);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = brusselator_reference(t_final, &x0, a, b);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(
            error < tol * 1000.0,
            "RKBS32 tol={} error={}",
            tol,
            error
        );
    }
}

#[test]
fn test_brusselator_rkf45() {
    let a = 0.4;
    let b = 1.2;
    let x0 = [0.0, 0.0];
    let t_final = 50.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            a - x[0] - b * x[0] + x[0] * x[0] * x[1],
            b * x[0] - x[0] * x[0] * x[1],
        ])
    };

    let tolerances = [1e-5, 1e-6, 1e-7];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = RKF45::new(initial);

        let mut t = 0.0;
        let mut dt = 0.01;

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
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.5);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = brusselator_reference(t_final, &x0, a, b);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(
            error < tol * 500.0,
            "RKF45 tol={} error={}",
            tol,
            error
        );
    }
}

#[test]
fn test_brusselator_rkck54() {
    let a = 0.4;
    let b = 1.2;
    let x0 = [0.0, 0.0];
    let t_final = 50.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            a - x[0] - b * x[0] + x[0] * x[0] * x[1],
            b * x[0] - x[0] * x[0] * x[1],
        ])
    };

    let tolerances = [1e-5, 1e-6, 1e-7];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = RKCK54::new(initial);

        let mut t = 0.0;
        let mut dt = 0.01;

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
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.5);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = brusselator_reference(t_final, &x0, a, b);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(
            error < tol * 500.0,
            "RKCK54 tol={} error={}",
            tol,
            error
        );
    }
}

#[test]
fn test_brusselator_rkdp54() {
    let a = 0.4;
    let b = 1.2;
    let x0 = [0.0, 0.0];
    let t_final = 50.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            a - x[0] - b * x[0] + x[0] * x[0] * x[1],
            b * x[0] - x[0] * x[0] * x[1],
        ])
    };

    let tolerances = [1e-5, 1e-6, 1e-7];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = RKDP54::new(initial);

        let mut t = 0.0;
        let mut dt = 0.01;

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
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.5);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = brusselator_reference(t_final, &x0, a, b);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(
            error < tol * 200.0,
            "RKDP54 tol={} error={}",
            tol,
            error
        );
    }
}

#[test]
fn test_brusselator_rkv65() {
    let a = 0.4;
    let b = 1.2;
    let x0 = [0.0, 0.0];
    let t_final = 50.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            a - x[0] - b * x[0] + x[0] * x[0] * x[1],
            b * x[0] - x[0] * x[0] * x[1],
        ])
    };

    let tolerances = [1e-5, 1e-6, 1e-7];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = RKV65::new(initial);

        let mut t = 0.0;
        let mut dt = 0.01;

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
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.5);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = brusselator_reference(t_final, &x0, a, b);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(
            error < tol * 200.0,
            "RKV65 tol={} error={}",
            tol,
            error
        );
    }
}

#[test]
fn test_brusselator_rkdp87() {
    let a = 0.4;
    let b = 1.2;
    let x0 = [0.0, 0.0];
    let t_final = 50.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            a - x[0] - b * x[0] + x[0] * x[0] * x[1],
            b * x[0] - x[0] * x[0] * x[1],
        ])
    };

    let tolerances = [1e-5, 1e-6, 1e-7];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = RKDP87::new(initial);

        let mut t = 0.0;
        let mut dt = 0.01;

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
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.5);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = brusselator_reference(t_final, &x0, a, b);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(
            error < tol * 100.0,
            "RKDP87 tol={} error={}",
            tol,
            error
        );
    }
}
