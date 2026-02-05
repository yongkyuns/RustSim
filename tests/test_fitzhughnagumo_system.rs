//! FitzHugh-Nagumo neuron model evaluation tests
//!
//! Tests solvers on the FitzHugh-Nagumo system:
//! dv/dt = v - v³/3 - w + R·Iext
//! dw/dt = (v + a - bw)/τ
//!
//! Parameters: a=0.7, b=0.8, τ=12.5, R=1.0, Iext=0.5
//! Initial conditions: v=0, w=0
//!
//! This is a simplified neuron model exhibiting excitability and oscillations

use approx::assert_relative_eq;
use nalgebra::DVector;
use rustsim::solvers::*;

/// High-accuracy reference solution for FitzHugh-Nagumo system
fn fitzhughnagumo_reference(
    t: f64,
    x0: &[f64; 2],
    a: f64,
    b: f64,
    tau: f64,
    r: f64,
    iext: f64,
) -> [f64; 2] {
    let mut state = DVector::from_vec(vec![x0[0], x0[1]]);
    let mut solver = RKDP54::new(state.clone());

    let dt_base = 0.0001;
    let steps = (t / dt_base) as usize;

    for _ in 0..steps {
        solver.buffer(dt_base);

        let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
            let v = x[0];
            let w = x[1];
            DVector::from_vec(vec![
                v - (1.0 / 3.0) * v * v * v - w + r * iext,
                (v + a - b * w) / tau,
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
fn test_fitzhughnagumo_rkbs32() {
    let a = 0.7;
    let b = 0.8;
    let tau = 12.5;
    let r = 1.0;
    let iext = 0.5;
    let x0 = [0.0, 0.0];
    let t_final = 50.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        let v = x[0];
        let w = x[1];
        DVector::from_vec(vec![
            v - (1.0 / 3.0) * v * v * v - w + r * iext,
            (v + a - b * w) / tau,
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

        let reference = fitzhughnagumo_reference(t_final, &x0, a, b, tau, r, iext);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(error < tol * 1000.0, "RKBS32 tol={} error={}", tol, error);
    }
}

#[test]
fn test_fitzhughnagumo_rkf45() {
    let a = 0.7;
    let b = 0.8;
    let tau = 12.5;
    let r = 1.0;
    let iext = 0.5;
    let x0 = [0.0, 0.0];
    let t_final = 50.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        let v = x[0];
        let w = x[1];
        DVector::from_vec(vec![
            v - (1.0 / 3.0) * v * v * v - w + r * iext,
            (v + a - b * w) / tau,
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

        let reference = fitzhughnagumo_reference(t_final, &x0, a, b, tau, r, iext);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(error < tol * 500.0, "RKF45 tol={} error={}", tol, error);
    }
}

#[test]
fn test_fitzhughnagumo_rkck54() {
    let a = 0.7;
    let b = 0.8;
    let tau = 12.5;
    let r = 1.0;
    let iext = 0.5;
    let x0 = [0.0, 0.0];
    let t_final = 50.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        let v = x[0];
        let w = x[1];
        DVector::from_vec(vec![
            v - (1.0 / 3.0) * v * v * v - w + r * iext,
            (v + a - b * w) / tau,
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

        let reference = fitzhughnagumo_reference(t_final, &x0, a, b, tau, r, iext);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(error < tol * 500.0, "RKCK54 tol={} error={}", tol, error);
    }
}

#[test]
fn test_fitzhughnagumo_rkdp54() {
    let a = 0.7;
    let b = 0.8;
    let tau = 12.5;
    let r = 1.0;
    let iext = 0.5;
    let x0 = [0.0, 0.0];
    let t_final = 50.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        let v = x[0];
        let w = x[1];
        DVector::from_vec(vec![
            v - (1.0 / 3.0) * v * v * v - w + r * iext,
            (v + a - b * w) / tau,
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

        let reference = fitzhughnagumo_reference(t_final, &x0, a, b, tau, r, iext);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(error < tol * 200.0, "RKDP54 tol={} error={}", tol, error);
    }
}

#[test]
fn test_fitzhughnagumo_rkv65() {
    let a = 0.7;
    let b = 0.8;
    let tau = 12.5;
    let r = 1.0;
    let iext = 0.5;
    let x0 = [0.0, 0.0];
    let t_final = 50.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        let v = x[0];
        let w = x[1];
        DVector::from_vec(vec![
            v - (1.0 / 3.0) * v * v * v - w + r * iext,
            (v + a - b * w) / tau,
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

        let reference = fitzhughnagumo_reference(t_final, &x0, a, b, tau, r, iext);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(error < tol * 200.0, "RKV65 tol={} error={}", tol, error);
    }
}

#[test]
fn test_fitzhughnagumo_rkdp87() {
    let a = 0.7;
    let b = 0.8;
    let tau = 12.5;
    let r = 1.0;
    let iext = 0.5;
    let x0 = [0.0, 0.0];
    let t_final = 50.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        let v = x[0];
        let w = x[1];
        DVector::from_vec(vec![
            v - (1.0 / 3.0) * v * v * v - w + r * iext,
            (v + a - b * w) / tau,
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

        let reference = fitzhughnagumo_reference(t_final, &x0, a, b, tau, r, iext);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(error < tol * 100.0, "RKDP87 tol={} error={}", tol, error);
    }
}
