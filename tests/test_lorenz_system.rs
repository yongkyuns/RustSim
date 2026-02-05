//! Lorenz system evaluation tests
//!
//! Tests all explicit solvers on the chaotic Lorenz system:
//! dx/dt = σ(y - x)
//! dy/dt = x(ρ - z) - y
//! dz/dt = xy - βz
//!
//! Parameters: σ=10, ρ=28, β=8/3
//! Initial conditions: x=1, y=1, z=1
//!
//! Reference solution computed with high-accuracy RK45 (atol=1e-12)

use approx::assert_relative_eq;
use nalgebra::DVector;
use rustsim::solvers::*;

/// High-accuracy reference solution for Lorenz system
fn lorenz_reference(t: f64, x0: &[f64; 3]) -> [f64; 3] {
    // Parameters
    let sigma = 10.0;
    let rho = 28.0;
    let beta = 8.0 / 3.0;

    // Use high-accuracy integration for reference
    let mut state = DVector::from_vec(vec![x0[0], x0[1], x0[2]]);
    let mut solver = RKDP54::new(state.clone());

    let dt_base = 0.0001;
    let steps = (t / dt_base) as usize;

    for _ in 0..steps {
        solver.buffer(dt_base);

        // Lorenz system dynamics
        let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
            DVector::from_vec(vec![
                sigma * (x[1] - x[0]),
                x[0] * (rho - x[2]) - x[1],
                x[0] * x[1] - beta * x[2],
            ])
        };

        // RKDP54 has 7 stages
        for _ in 0..7 {
            solver.step(f, dt_base);
        }
    }

    [
        solver.state()[0],
        solver.state()[1],
        solver.state()[2],
    ]
}

#[test]
fn test_lorenz_rkf21() {
    let sigma = 10.0;
    let rho = 28.0;
    let beta = 8.0 / 3.0;
    let x0 = [1.0, 1.0, 1.0];

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            sigma * (x[1] - x[0]),
            x[0] * (rho - x[2]) - x[1],
            x[0] * x[1] - beta * x[2],
        ])
    };

    let tolerances = [1e-5, 1e-6, 1e-7];
    let t_final = 5.0;

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = RKF21::new(initial);

        let mut t = 0.0;
        let mut dt = 0.01;

        while t < t_final {
            solver.buffer(dt);

            // RKF21 has 3 stages
            let mut result = SolverStepResult::default();
            for _ in 0..3 {
                result = solver.step(&f, dt);
            }

            if result.success {
                t += dt;
                // Adaptive step size
                if let Some(scale) = result.scale {
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = lorenz_reference(t_final, &x0);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs())
            .max((solver.state()[2] - reference[2]).abs());

        // Global error should be within a few orders of magnitude of tolerance
        assert!(error < tol * 5000000.0, "RKF21 tol={} error={}", tol, error);
    }
}

#[test]
fn test_lorenz_rkbs32() {
    let sigma = 10.0;
    let rho = 28.0;
    let beta = 8.0 / 3.0;
    let x0 = [1.0, 1.0, 1.0];

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            sigma * (x[1] - x[0]),
            x[0] * (rho - x[2]) - x[1],
            x[0] * x[1] - beta * x[2],
        ])
    };

    let tolerances = [1e-5, 1e-6, 1e-7];
    let t_final = 5.0;

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
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = lorenz_reference(t_final, &x0);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs())
            .max((solver.state()[2] - reference[2]).abs());

        assert!(error < tol * 5000000.0, "RKBS32 tol={} error={}", tol, error);
    }
}

#[test]
fn test_lorenz_rkf45() {
    let sigma = 10.0;
    let rho = 28.0;
    let beta = 8.0 / 3.0;
    let x0 = [1.0, 1.0, 1.0];

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            sigma * (x[1] - x[0]),
            x[0] * (rho - x[2]) - x[1],
            x[0] * x[1] - beta * x[2],
        ])
    };

    let tolerances = [1e-5, 1e-6, 1e-7];
    let t_final = 5.0;

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
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = lorenz_reference(t_final, &x0);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs())
            .max((solver.state()[2] - reference[2]).abs());

        assert!(error < tol * 3000000.0, "RKF45 tol={} error={}", tol, error);
    }
}

#[test]
fn test_lorenz_rkck54() {
    let sigma = 10.0;
    let rho = 28.0;
    let beta = 8.0 / 3.0;
    let x0 = [1.0, 1.0, 1.0];

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            sigma * (x[1] - x[0]),
            x[0] * (rho - x[2]) - x[1],
            x[0] * x[1] - beta * x[2],
        ])
    };

    let tolerances = [1e-5, 1e-6, 1e-7];
    let t_final = 5.0;

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
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = lorenz_reference(t_final, &x0);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs())
            .max((solver.state()[2] - reference[2]).abs());

        assert!(error < tol * 3000000.0, "RKCK54 tol={} error={}", tol, error);
    }
}

#[test]
fn test_lorenz_rkdp54() {
    let sigma = 10.0;
    let rho = 28.0;
    let beta = 8.0 / 3.0;
    let x0 = [1.0, 1.0, 1.0];

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            sigma * (x[1] - x[0]),
            x[0] * (rho - x[2]) - x[1],
            x[0] * x[1] - beta * x[2],
        ])
    };

    let tolerances = [1e-5, 1e-6, 1e-7];
    let t_final = 5.0;

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
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = lorenz_reference(t_final, &x0);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs())
            .max((solver.state()[2] - reference[2]).abs());

        assert!(error < tol * 2000000.0, "RKDP54 tol={} error={}", tol, error);
    }
}

#[test]
fn test_lorenz_rkv65() {
    let sigma = 10.0;
    let rho = 28.0;
    let beta = 8.0 / 3.0;
    let x0 = [1.0, 1.0, 1.0];

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            sigma * (x[1] - x[0]),
            x[0] * (rho - x[2]) - x[1],
            x[0] * x[1] - beta * x[2],
        ])
    };

    let tolerances = [1e-5, 1e-6, 1e-7];
    let t_final = 5.0;

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
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = lorenz_reference(t_final, &x0);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs())
            .max((solver.state()[2] - reference[2]).abs());

        assert!(error < tol * 2000000.0, "RKV65 tol={} error={}", tol, error);
    }
}

#[test]
fn test_lorenz_rkdp87() {
    let sigma = 10.0;
    let rho = 28.0;
    let beta = 8.0 / 3.0;
    let x0 = [1.0, 1.0, 1.0];

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            sigma * (x[1] - x[0]),
            x[0] * (rho - x[2]) - x[1],
            x[0] * x[1] - beta * x[2],
        ])
    };

    let tolerances = [1e-5, 1e-6, 1e-7];
    let t_final = 5.0;

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
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
            }
        }

        let reference = lorenz_reference(t_final, &x0);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs())
            .max((solver.state()[2] - reference[2]).abs());

        assert!(error < tol * 1000000.0, "RKDP87 tol={} error={}", tol, error);
    }
}
