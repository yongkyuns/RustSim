//! Van der Pol oscillator evaluation tests
//!
//! Tests implicit solvers on the stiff Van der Pol system:
//! dx/dt = y
//! dy/dt = μ(1 - x²)y - x
//!
//! Parameters: μ=10 (stiff case)
//! Initial conditions: x=2, y=0
//!
//! This is a classic test for stiff solvers due to the μ(1-x²)y term

use approx::assert_relative_eq;
use nalgebra::DVector;
use rustsim::solvers::*;

/// High-accuracy reference solution for Van der Pol oscillator
fn vanderpol_reference(t: f64, x0: &[f64; 2], mu: f64) -> [f64; 2] {
    // Use high-accuracy implicit solver for reference
    let mut state = DVector::from_vec(vec![x0[0], x0[1]]);
    let mut solver = ESDIRK54::new(state.clone());

    let dt_base = 0.0001;
    let steps = (t / dt_base) as usize;

    for _ in 0..steps {
        solver.buffer(dt_base);

        let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
            DVector::from_vec(vec![
                x[1],
                mu * (1.0 - x[0] * x[0]) * x[1] - x[0],
            ])
        };

        // ESDIRK54 has 7 stages
        for _ in 0..7 {
            solver.step(&f, dt_base);
        }
    }

    [solver.state()[0], solver.state()[1]]
}

#[test]
fn test_vanderpol_esdirk32() {
    let mu = 10.0;
    let x0 = [2.0, 0.0];
    let t_final = mu;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            x[1],
            mu * (1.0 - x[0] * x[0]) * x[1] - x[0],
        ])
    };

    let tolerances = [1e-6];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = ESDIRK32::new(initial);

        let mut t = 0.0;
        let mut dt = 0.01;

        while t < t_final {
            solver.buffer(dt);

            // ESDIRK32 has 3 stages
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
                if dt < 1e-8 {
                    panic!("Step size too small");
                }
            }
        }

        let reference = vanderpol_reference(t_final, &x0, mu);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(
            error < tol * 100.0,
            "ESDIRK32 tol={} error={}",
            tol,
            error
        );
    }
}

#[test]
fn test_vanderpol_esdirk43() {
    let mu = 10.0;
    let x0 = [2.0, 0.0];
    let t_final = mu;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            x[1],
            mu * (1.0 - x[0] * x[0]) * x[1] - x[0],
        ])
    };

    let tolerances = [1e-6];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = ESDIRK43::new(initial);

        let mut t = 0.0;
        let mut dt = 0.01;

        while t < t_final {
            solver.buffer(dt);

            // ESDIRK43 has 5 stages
            let mut result = SolverStepResult::default();
            for _ in 0..5 {
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
                if dt < 1e-8 {
                    panic!("Step size too small");
                }
            }
        }

        let reference = vanderpol_reference(t_final, &x0, mu);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(
            error < tol * 100.0,
            "ESDIRK43 tol={} error={}",
            tol,
            error
        );
    }
}

#[test]
fn test_vanderpol_gear21() {
    let mu = 10.0;
    let x0 = [2.0, 0.0];
    let t_final = mu;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            x[1],
            mu * (1.0 - x[0] * x[0]) * x[1] - x[0],
        ])
    };

    let tolerances = [1e-6];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = GEAR21::new(initial);

        let mut t = 0.0;
        let mut dt = 0.01;

        while t < t_final {
            solver.buffer(dt);

            // GEAR21 is implicit, single stage
            let result = solver.step(&f, dt);

            if result.success {
                t += dt;
                if let Some(scale) = result.scale {
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
                if dt < 1e-8 {
                    panic!("Step size too small");
                }
            }
        }

        let reference = vanderpol_reference(t_final, &x0, mu);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(
            error < tol * 200.0,
            "GEAR21 tol={} error={}",
            tol,
            error
        );
    }
}

#[test]
fn test_vanderpol_gear32() {
    let mu = 10.0;
    let x0 = [2.0, 0.0];
    let t_final = mu;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            x[1],
            mu * (1.0 - x[0] * x[0]) * x[1] - x[0],
        ])
    };

    let tolerances = [1e-6];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = GEAR32::new(initial);

        let mut t = 0.0;
        let mut dt = 0.01;

        while t < t_final {
            solver.buffer(dt);

            let result = solver.step(&f, dt);

            if result.success {
                t += dt;
                if let Some(scale) = result.scale {
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
                if dt < 1e-8 {
                    panic!("Step size too small");
                }
            }
        }

        let reference = vanderpol_reference(t_final, &x0, mu);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(
            error < tol * 150.0,
            "GEAR32 tol={} error={}",
            tol,
            error
        );
    }
}

#[test]
fn test_vanderpol_gear43() {
    let mu = 10.0;
    let x0 = [2.0, 0.0];
    let t_final = mu;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            x[1],
            mu * (1.0 - x[0] * x[0]) * x[1] - x[0],
        ])
    };

    let tolerances = [1e-6];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = GEAR43::new(initial);

        let mut t = 0.0;
        let mut dt = 0.01;

        while t < t_final {
            solver.buffer(dt);

            let result = solver.step(&f, dt);

            if result.success {
                t += dt;
                if let Some(scale) = result.scale {
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
                if dt < 1e-8 {
                    panic!("Step size too small");
                }
            }
        }

        let reference = vanderpol_reference(t_final, &x0, mu);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(
            error < tol * 100.0,
            "GEAR43 tol={} error={}",
            tol,
            error
        );
    }
}

#[test]
fn test_vanderpol_gear54() {
    let mu = 10.0;
    let x0 = [2.0, 0.0];
    let t_final = mu;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            x[1],
            mu * (1.0 - x[0] * x[0]) * x[1] - x[0],
        ])
    };

    let tolerances = [1e-6];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = GEAR54::new(initial);

        let mut t = 0.0;
        let mut dt = 0.01;

        while t < t_final {
            solver.buffer(dt);

            let result = solver.step(&f, dt);

            if result.success {
                t += dt;
                if let Some(scale) = result.scale {
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
                if dt < 1e-8 {
                    panic!("Step size too small");
                }
            }
        }

        let reference = vanderpol_reference(t_final, &x0, mu);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(
            error < tol * 100.0,
            "GEAR54 tol={} error={}",
            tol,
            error
        );
    }
}

#[test]
fn test_vanderpol_gear52a() {
    let mu = 10.0;
    let x0 = [2.0, 0.0];
    let t_final = mu;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            x[1],
            mu * (1.0 - x[0] * x[0]) * x[1] - x[0],
        ])
    };

    let tolerances = [1e-6];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = GEAR52A::new(initial);

        let mut t = 0.0;
        let mut dt = 0.01;

        while t < t_final {
            solver.buffer(dt);

            let result = solver.step(&f, dt);

            if result.success {
                t += dt;
                if let Some(scale) = result.scale {
                    dt = (dt * scale.min(2.0).max(0.5)).min(0.1);
                }
            } else {
                solver.revert().unwrap();
                dt *= 0.5;
                if dt < 1e-8 {
                    panic!("Step size too small");
                }
            }
        }

        let reference = vanderpol_reference(t_final, &x0, mu);
        let error = (solver.state()[0] - reference[0]).abs()
            .max((solver.state()[1] - reference[1]).abs());

        assert!(
            error < tol * 100.0,
            "GEAR52A tol={} error={}",
            tol,
            error
        );
    }
}
