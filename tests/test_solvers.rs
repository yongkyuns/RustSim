//! Integration tests for numerical solvers
//!
//! Tests all solver implementations against known ODE solutions

use approx::assert_relative_eq;
use nalgebra::DVector;

// Import solver types and traits from the solvers module
mod solvers {
    pub use rustsim::solvers::*;
}

use solvers::{Euler, ExplicitSolver, RK4, RKDP54, Solver, SSPRK22, SSPRK33};

#[test]
fn test_all_solvers_exponential_decay() {
    // Test equation: dx/dt = -x, x(0) = 1
    // Exact solution: x(t) = exp(-t)

    let x0 = DVector::from_vec(vec![1.0]);
    let dt = 0.1;
    let t_final = 1.0;
    let n_steps = (t_final / dt) as usize;
    let exact = (-t_final as f64).exp();

    // Test Euler
    {
        let mut solver = Euler::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            solver.step(|x, _t| -x, dt);
        }
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-2);
    }

    // Test SSPRK22
    {
        let mut solver = SSPRK22::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..2 {
                solver.step(|x, _t| -x, dt);
            }
        }
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-3);
    }

    // Test SSPRK33
    {
        let mut solver = SSPRK33::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..3 {
                solver.step(|x, _t| -x, dt);
            }
        }
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-4);
    }

    // Test RK4
    {
        let mut solver = RK4::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..4 {
                solver.step(|x, _t| -x, dt);
            }
        }
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-6);
    }

    // Test RKDP54
    {
        let mut solver = RKDP54::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..7 {
                solver.step(|x, _t| -x, dt);
            }
        }
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-8);
    }
}

#[test]
fn test_solver_properties() {
    let x0 = DVector::from_vec(vec![1.0]);

    // Check solver properties
    let euler = Euler::new(x0.clone());
    assert_eq!(euler.order(), 1);
    assert_eq!(euler.stages(), 1);
    assert!(!euler.is_adaptive());
    assert!(euler.is_explicit());

    let ssprk22 = SSPRK22::new(x0.clone());
    assert_eq!(ssprk22.order(), 2);
    assert_eq!(ssprk22.stages(), 2);
    assert!(!ssprk22.is_adaptive());

    let ssprk33 = SSPRK33::new(x0.clone());
    assert_eq!(ssprk33.order(), 3);
    assert_eq!(ssprk33.stages(), 3);
    assert!(!ssprk33.is_adaptive());

    let rk4 = RK4::new(x0.clone());
    assert_eq!(rk4.order(), 4);
    assert_eq!(rk4.stages(), 4);
    assert!(!rk4.is_adaptive());

    let rkdp54 = RKDP54::new(x0.clone());
    assert_eq!(rkdp54.order(), 5);
    assert_eq!(rkdp54.stages(), 7);
    assert!(rkdp54.is_adaptive());
}

#[test]
fn test_rkdp54_adaptive_error_control() {
    // Verify that RKDP54 provides error estimates
    let x0 = DVector::from_vec(vec![1.0]);
    let mut solver = RKDP54::new(x0);

    solver.buffer(0.1);

    // Step through all 7 stages
    let mut result = None;
    for _ in 0..7 {
        result = Some(solver.step(|x, _t| -x, 0.1));
    }

    let result = result.unwrap();

    // Should have a scale factor on the last stage
    assert!(result.scale.is_some());

    // Scale should be in reasonable range [0.1, 10.0]
    let scale = result.scale.unwrap();
    assert!(scale >= 0.1 && scale <= 10.0);

    // For this smooth problem, step should succeed
    assert!(result.success);
}

#[test]
fn test_solver_reset() {
    let x0 = DVector::from_vec(vec![1.0]);
    let mut solver = RK4::new(x0.clone());

    // Integrate for some steps
    for _ in 0..5 {
        solver.buffer(0.1);
        for _ in 0..4 {
            solver.step(|x, _t| -x, 0.1);
        }
    }

    // State should have changed
    assert!(solver.state()[0] < 1.0);

    // Reset should restore initial state
    solver.reset();
    assert_eq!(solver.state()[0], 1.0);
}

#[test]
fn test_solver_revert() {
    let x0 = DVector::from_vec(vec![1.0]);
    let mut solver = RK4::new(x0);

    // Buffer initial state
    solver.buffer(0.1);
    let buffered_state = solver.state()[0];

    // Take a step
    for _ in 0..4 {
        solver.step(|x, _t| -x, 0.1);
    }

    // State should have changed
    assert!(solver.state()[0] < buffered_state);

    // Revert should restore buffered state
    solver.revert().unwrap();
    assert_eq!(solver.state()[0], buffered_state);
}

#[test]
fn test_harmonic_oscillator() {
    // d²x/dt² = -x => [x, v]' = [v, -x]
    // Exact: x(t) = cos(t), v(t) = -sin(t)

    let x0 = DVector::from_vec(vec![1.0, 0.0]);
    let mut solver = RK4::new(x0);

    let dt = 0.01;
    let t_final = 2.0 * std::f64::consts::PI;
    let n_steps = (t_final / dt) as usize;

    for _ in 0..n_steps {
        solver.buffer(dt);
        for _ in 0..4 {
            solver.step(
                |x, _t| {
                    let mut dxdt = DVector::zeros(2);
                    dxdt[0] = x[1];
                    dxdt[1] = -x[0];
                    dxdt
                },
                dt,
            );
        }
    }

    // After one period, should return to initial state
    assert_relative_eq!(solver.state()[0], 1.0, epsilon = 1e-4);
    assert_relative_eq!(solver.state()[1], 0.0, epsilon = 1e-4);
}
