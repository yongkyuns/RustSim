//! Robertson chemical kinetics system evaluation tests
//!
//! Tests implicit solvers on the highly stiff Robertson system:
//! dx/dt = -0.04x + 10⁴yz
//! dy/dt = 0.04x - 10⁴yz - 3×10⁷y²
//! dz/dt = 3×10⁷y²
//!
//! Initial conditions: x=1, y=0, z=0
//!
//! This is one of the most challenging stiff test problems with
//! time scales spanning many orders of magnitude

use approx::assert_relative_eq;
use nalgebra::DVector;
use rustsim::solvers::*;

/// High-accuracy reference solution for Robertson system
fn robertson_reference(t: f64, x0: &[f64; 3]) -> [f64; 3] {
    let a = 0.04;
    let b = 1e4;
    let c = 3e7;

    // Use high-accuracy implicit solver for reference
    let mut state = DVector::from_vec(vec![x0[0], x0[1], x0[2]]);
    let mut solver = ESDIRK85::new(state.clone());

    let dt_base = 0.00001;
    let steps = (t / dt_base) as usize;

    for _ in 0..steps {
        solver.buffer(dt_base);

        let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
            DVector::from_vec(vec![
                -a * x[0] + b * x[1] * x[2],
                a * x[0] - b * x[1] * x[2] - c * x[1] * x[1],
                c * x[1] * x[1],
            ])
        };

        // ESDIRK85 has 9 stages
        for _ in 0..9 {
            solver.step(&f, dt_base);
        }
    }

    [solver.state()[0], solver.state()[1], solver.state()[2]]
}

#[test]
fn test_robertson_esdirk32() {
    let a = 0.04;
    let b = 1e4;
    let c = 3e7;
    let x0 = [1.0, 0.0, 0.0];
    let t_final = 10.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            -a * x[0] + b * x[1] * x[2],
            a * x[0] - b * x[1] * x[2] - c * x[1] * x[1],
            c * x[1] * x[1],
        ])
    };

    let tolerances = [1e-6];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = ESDIRK32::new(initial);

        let mut t = 0.0;
        let mut dt = 0.001;

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
                if dt < 1e-10 {
                    panic!("Step size too small at t={}", t);
                }
            }
        }

        let reference = robertson_reference(t_final, &x0);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs())
            .max((solver.state()[2] - reference[2]).abs());

        assert!(error < tol * 100.0, "ESDIRK32 tol={} error={}", tol, error);
    }
}

#[test]
fn test_robertson_esdirk43() {
    let a = 0.04;
    let b = 1e4;
    let c = 3e7;
    let x0 = [1.0, 0.0, 0.0];
    let t_final = 10.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            -a * x[0] + b * x[1] * x[2],
            a * x[0] - b * x[1] * x[2] - c * x[1] * x[1],
            c * x[1] * x[1],
        ])
    };

    let tolerances = [1e-6];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = ESDIRK43::new(initial);

        let mut t = 0.0;
        let mut dt = 0.001;

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
                if dt < 1e-10 {
                    panic!("Step size too small at t={}", t);
                }
            }
        }

        let reference = robertson_reference(t_final, &x0);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs())
            .max((solver.state()[2] - reference[2]).abs());

        assert!(error < tol * 100.0, "ESDIRK43 tol={} error={}", tol, error);
    }
}

#[test]
fn test_robertson_esdirk54() {
    let a = 0.04;
    let b = 1e4;
    let c = 3e7;
    let x0 = [1.0, 0.0, 0.0];
    let t_final = 10.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            -a * x[0] + b * x[1] * x[2],
            a * x[0] - b * x[1] * x[2] - c * x[1] * x[1],
            c * x[1] * x[1],
        ])
    };

    let tolerances = [1e-6];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = ESDIRK54::new(initial);

        let mut t = 0.0;
        let mut dt = 0.001;

        while t < t_final {
            solver.buffer(dt);

            // ESDIRK54 has 7 stages
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
                if dt < 1e-10 {
                    panic!("Step size too small at t={}", t);
                }
            }
        }

        let reference = robertson_reference(t_final, &x0);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs())
            .max((solver.state()[2] - reference[2]).abs());

        assert!(error < tol * 50.0, "ESDIRK54 tol={} error={}", tol, error);
    }
}

#[test]
fn test_robertson_esdirk85() {
    let a = 0.04;
    let b = 1e4;
    let c = 3e7;
    let x0 = [1.0, 0.0, 0.0];
    let t_final = 10.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            -a * x[0] + b * x[1] * x[2],
            a * x[0] - b * x[1] * x[2] - c * x[1] * x[1],
            c * x[1] * x[1],
        ])
    };

    let tolerances = [1e-6];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = ESDIRK85::new(initial);

        let mut t = 0.0;
        let mut dt = 0.001;

        while t < t_final {
            solver.buffer(dt);

            // ESDIRK85 has 9 stages
            let mut result = SolverStepResult::default();
            for _ in 0..9 {
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
                if dt < 1e-10 {
                    panic!("Step size too small at t={}", t);
                }
            }
        }

        let reference = robertson_reference(t_final, &x0);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs())
            .max((solver.state()[2] - reference[2]).abs());

        assert!(error < tol * 20.0, "ESDIRK85 tol={} error={}", tol, error);
    }
}

#[test]
fn test_robertson_gear21() {
    let a = 0.04;
    let b = 1e4;
    let c = 3e7;
    let x0 = [1.0, 0.0, 0.0];
    let t_final = 10.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            -a * x[0] + b * x[1] * x[2],
            a * x[0] - b * x[1] * x[2] - c * x[1] * x[1],
            c * x[1] * x[1],
        ])
    };

    let tolerances = [1e-6];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = GEAR21::new(initial);

        let mut t = 0.0;
        let mut dt = 0.001;

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
                if dt < 1e-10 {
                    panic!("Step size too small at t={}", t);
                }
            }
        }

        let reference = robertson_reference(t_final, &x0);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs())
            .max((solver.state()[2] - reference[2]).abs());

        assert!(error < tol * 200.0, "GEAR21 tol={} error={}", tol, error);
    }
}

#[test]
fn test_robertson_gear32() {
    let a = 0.04;
    let b = 1e4;
    let c = 3e7;
    let x0 = [1.0, 0.0, 0.0];
    let t_final = 10.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            -a * x[0] + b * x[1] * x[2],
            a * x[0] - b * x[1] * x[2] - c * x[1] * x[1],
            c * x[1] * x[1],
        ])
    };

    let tolerances = [1e-6];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = GEAR32::new(initial);

        let mut t = 0.0;
        let mut dt = 0.001;

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
                if dt < 1e-10 {
                    panic!("Step size too small at t={}", t);
                }
            }
        }

        let reference = robertson_reference(t_final, &x0);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs())
            .max((solver.state()[2] - reference[2]).abs());

        assert!(error < tol * 150.0, "GEAR32 tol={} error={}", tol, error);
    }
}

#[test]
fn test_robertson_gear43() {
    let a = 0.04;
    let b = 1e4;
    let c = 3e7;
    let x0 = [1.0, 0.0, 0.0];
    let t_final = 10.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            -a * x[0] + b * x[1] * x[2],
            a * x[0] - b * x[1] * x[2] - c * x[1] * x[1],
            c * x[1] * x[1],
        ])
    };

    let tolerances = [1e-6];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = GEAR43::new(initial);

        let mut t = 0.0;
        let mut dt = 0.001;

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
                if dt < 1e-10 {
                    panic!("Step size too small at t={}", t);
                }
            }
        }

        let reference = robertson_reference(t_final, &x0);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs())
            .max((solver.state()[2] - reference[2]).abs());

        assert!(error < tol * 100.0, "GEAR43 tol={} error={}", tol, error);
    }
}

#[test]
fn test_robertson_gear54() {
    let a = 0.04;
    let b = 1e4;
    let c = 3e7;
    let x0 = [1.0, 0.0, 0.0];
    let t_final = 10.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            -a * x[0] + b * x[1] * x[2],
            a * x[0] - b * x[1] * x[2] - c * x[1] * x[1],
            c * x[1] * x[1],
        ])
    };

    let tolerances = [1e-6];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = GEAR54::new(initial);

        let mut t = 0.0;
        let mut dt = 0.001;

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
                if dt < 1e-10 {
                    panic!("Step size too small at t={}", t);
                }
            }
        }

        let reference = robertson_reference(t_final, &x0);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs())
            .max((solver.state()[2] - reference[2]).abs());

        assert!(error < tol * 100.0, "GEAR54 tol={} error={}", tol, error);
    }
}

#[test]
fn test_robertson_gear52a() {
    let a = 0.04;
    let b = 1e4;
    let c = 3e7;
    let x0 = [1.0, 0.0, 0.0];
    let t_final = 10.0;

    let f = |x: &DVector<f64>, _t: f64| -> DVector<f64> {
        DVector::from_vec(vec![
            -a * x[0] + b * x[1] * x[2],
            a * x[0] - b * x[1] * x[2] - c * x[1] * x[1],
            c * x[1] * x[1],
        ])
    };

    let tolerances = [1e-6];

    for &tol in &tolerances {
        let initial = DVector::from_vec(x0.to_vec());
        let mut solver = GEAR52A::new(initial);

        let mut t = 0.0;
        let mut dt = 0.001;

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
                if dt < 1e-10 {
                    panic!("Step size too small at t={}", t);
                }
            }
        }

        let reference = robertson_reference(t_final, &x0);
        let error = (solver.state()[0] - reference[0])
            .abs()
            .max((solver.state()[1] - reference[1]).abs())
            .max((solver.state()[2] - reference[2]).abs());

        assert!(error < tol * 100.0, "GEAR52A tol={} error={}", tol, error);
    }
}
