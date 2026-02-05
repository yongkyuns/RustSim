//! Comprehensive tests for all Runge-Kutta solvers
//!
//! Tests match PathSim test structure exactly with:
//! - Convergence order verification
//! - Adaptive timestepping tests
//! - Stage iterator tests
//! - Comparison across all solvers

use approx::assert_relative_eq;
use nalgebra::DVector;
use rustsim::solvers::{
    ExplicitSolver, RK4, RKBS32, RKCK54, RKDP54, RKDP87, RKF21, RKF45, RKF78, RKV65, Solver,
    SSPRK22, SSPRK33, SSPRK34,
};

// Reference problems for testing
struct Problem {
    name: &'static str,
    x0: DVector<f64>,
    t_span: (f64, f64),
}

fn exponential_decay() -> Problem {
    Problem {
        name: "exponential_decay",
        x0: DVector::from_vec(vec![1.0]),
        t_span: (0.0, 5.0),
    }
}

fn logistic() -> Problem {
    Problem {
        name: "logistic",
        x0: DVector::from_vec(vec![0.5]),
        t_span: (0.0, 10.0),
    }
}

fn quadratic() -> Problem {
    Problem {
        name: "quadratic",
        x0: DVector::from_vec(vec![1.0]),
        t_span: (0.0, 0.6),
    }
}

fn sin_decay() -> Problem {
    Problem {
        name: "sin_decay",
        x0: DVector::from_vec(vec![1.0]),
        t_span: (0.0, 10.0),
    }
}

fn polynomial() -> Problem {
    Problem {
        name: "polynomial",
        x0: DVector::from_vec(vec![0.0]),
        t_span: (0.0, 5.0),
    }
}

fn get_rhs(problem: &Problem) -> impl Fn(&DVector<f64>, f64) -> DVector<f64> + '_ {
    move |x: &DVector<f64>, t: f64| -> DVector<f64> {
        match problem.name {
            "exponential_decay" => -x,
            "logistic" => x.map(|xi| xi * (1.0 - xi)),
            "quadratic" => x.map(|xi| xi * xi),
            "sin_decay" => x.map(|xi| -xi * t.sin()),
            "polynomial" => {
                let mut result = DVector::zeros(1);
                result[0] = t * t - x[0];
                result
            }
            _ => panic!("Unknown problem"),
        }
    }
}

fn get_exact_solution(problem: &Problem, t: f64) -> DVector<f64> {
    match problem.name {
        "exponential_decay" => DVector::from_vec(vec![(-t).exp()]),
        "logistic" => DVector::from_vec(vec![1.0 / (1.0 + (-t).exp())]),
        "quadratic" => DVector::from_vec(vec![1.0 / (1.0 - t)]),
        "sin_decay" => DVector::from_vec(vec![(t.cos() - 1.0).exp()]),
        "polynomial" => {
            DVector::from_vec(vec![t * t - 2.0 * t + 2.0 - 2.0 * (-t).exp()])
        }
        _ => panic!("Unknown problem"),
    }
}

// ============================================================================
// RKF21 Tests
// ============================================================================

#[test]
fn test_rkf21_properties() {
    let x0 = DVector::from_vec(vec![1.0]);
    let solver = RKF21::new(x0);

    assert_eq!(solver.order(), 2);
    assert_eq!(solver.stages(), 3);
    assert!(solver.is_adaptive());
    assert!(solver.is_explicit());
}

#[test]
fn test_rkf21_stage_iterator() {
    let problem = exponential_decay();
    let mut solver = RKF21::new(problem.x0);
    let dt = 0.1;

    solver.buffer(dt);

    for stage in 0..solver.stages() {
        let result = solver.step(get_rhs(&problem), dt);

        if stage < solver.stages() - 1 {
            // Intermediate stages
            assert!(result.success);
            assert_eq!(result.error_norm, 0.0);
            assert!(result.scale.is_none());
        } else {
            // Final stage
            assert!(result.scale.is_some());
            let scale = result.scale.unwrap();
            assert!(scale >= 0.1 && scale <= 10.0);
        }
    }
}

#[test]
fn test_rkf21_convergence() {
    let problem = exponential_decay();
    let duration = problem.t_span.1 - problem.t_span.0;

    // Test with different timesteps
    let divisions = vec![10.0, 20.0, 40.0, 80.0, 160.0];
    let mut errors = Vec::new();

    for div in &divisions {
        let mut solver = RKF21::new(problem.x0.clone());
        let dt = duration / div;
        let n_steps = *div as usize;

        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..solver.stages() {
                solver.step(get_rhs(&problem), dt);
            }
        }

        let exact = get_exact_solution(&problem, problem.t_span.1);
        let error = (solver.state() - exact).norm();
        errors.push(error);
    }

    // Errors should be monotonically decreasing
    for i in 1..errors.len() {
        assert!(errors[i] < errors[i - 1]);
    }

    // Check convergence order (should be at least 1st order globally)
    let log_h: Vec<f64> = divisions.iter().map(|d| (duration / d).ln()).collect();
    let log_e: Vec<f64> = errors.iter().map(|e: &f64| e.ln()).collect();

    // Simple linear regression to estimate order
    let n = log_h.len() as f64;
    let sum_x: f64 = log_h.iter().sum();
    let sum_y: f64 = log_e.iter().sum();
    let sum_xx: f64 = log_h.iter().map(|x| x * x).sum();
    let sum_xy: f64 = log_h.iter().zip(log_e.iter()).map(|(x, y)| x * y).sum();

    let slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);

    // For RKF21, expect convergence order > 1 (2nd order method)
    assert!(slope > 1.0, "Convergence order {} is too low", slope);
}

#[test]
fn test_rkf21_adaptive() {
    let problem = exponential_decay();
    let mut solver = RKF21::with_tolerances(problem.x0.clone(), 1e-4, 0.0);
    let duration = problem.t_span.1 - problem.t_span.0;
    let dt = duration / 100.0;

    let mut t = problem.t_span.0;
    let mut current_dt = dt;

    while t < problem.t_span.1 {
        solver.buffer(current_dt);

        let mut result = None;
        for _ in 0..solver.stages() {
            result = Some(solver.step(get_rhs(&problem), current_dt));
        }

        let result = result.unwrap();

        if result.success {
            t += current_dt;
            if let Some(scale) = result.scale {
                let new_dt = current_dt * scale;
                let remaining = problem.t_span.1 - t;
                current_dt = new_dt.min(remaining);
            }
        } else {
            solver.revert().unwrap();
            if let Some(scale) = result.scale {
                current_dt *= scale;
            }
        }
    }

    let exact = get_exact_solution(&problem, problem.t_span.1);
    let error = (solver.state() - exact).norm();

    // Error should be within tolerance range
    assert!(
        error < solver.tol_abs * 50.0,
        "Error {} exceeds tolerance",
        error
    );
}

// ============================================================================
// RKBS32 Tests
// ============================================================================

#[test]
fn test_rkbs32_properties() {
    let x0 = DVector::from_vec(vec![1.0]);
    let solver = RKBS32::new(x0);

    assert_eq!(solver.order(), 3);
    assert_eq!(solver.stages(), 4);
    assert!(solver.is_adaptive());
    assert!(solver.is_explicit());
}

#[test]
fn test_rkbs32_convergence() {
    let problem = exponential_decay();
    let duration = problem.t_span.1 - problem.t_span.0;

    let divisions = vec![10.0, 20.0, 40.0, 80.0];
    let mut errors = Vec::new();

    for div in &divisions {
        let mut solver = RKBS32::new(problem.x0.clone());
        let dt = duration / div;
        let n_steps = *div as usize;

        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..solver.stages() {
                solver.step(get_rhs(&problem), dt);
            }
        }

        let exact = get_exact_solution(&problem, problem.t_span.1);
        let error = (solver.state() - exact).norm();
        errors.push(error);
    }

    // Check monotonicity
    for i in 1..errors.len() {
        assert!(errors[i] < errors[i - 1]);
    }

    // Check convergence order (should be at least 2nd order)
    let log_h: Vec<f64> = divisions.iter().map(|d| (duration / d).ln()).collect();
    let log_e: Vec<f64> = errors.iter().map(|e: &f64| e.ln()).collect();

    let n = log_h.len() as f64;
    let sum_x: f64 = log_h.iter().sum();
    let sum_y: f64 = log_e.iter().sum();
    let sum_xx: f64 = log_h.iter().map(|x| x * x).sum();
    let sum_xy: f64 = log_h.iter().zip(log_e.iter()).map(|(x, y)| x * y).sum();

    let slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);

    assert!(slope > 2.0, "Convergence order {} is too low", slope);
}

#[test]
fn test_rkbs32_adaptive() {
    let problem = exponential_decay();
    let mut solver = RKBS32::with_tolerances(problem.x0.clone(), 1e-3, 0.0);
    let duration = problem.t_span.1 - problem.t_span.0;
    let dt = duration / 100.0;

    let mut t = problem.t_span.0;
    let mut current_dt = dt;

    while t < problem.t_span.1 {
        solver.buffer(current_dt);

        let mut result = None;
        for _ in 0..solver.stages() {
            result = Some(solver.step(get_rhs(&problem), current_dt));
        }

        let result = result.unwrap();

        if result.success {
            t += current_dt;
            if let Some(scale) = result.scale {
                let new_dt = current_dt * scale;
                let remaining = problem.t_span.1 - t;
                current_dt = new_dt.min(remaining);
            }
        } else {
            solver.revert().unwrap();
            if let Some(scale) = result.scale {
                current_dt *= scale;
            }
        }
    }

    let exact = get_exact_solution(&problem, problem.t_span.1);
    let error = (solver.state() - exact).norm();

    assert!(
        error < solver.tol_abs * 10.0,
        "Error {} exceeds tolerance",
        error
    );
}

// ============================================================================
// RKF45, RKCK54, RKV65, RKF78, RKDP87 Tests
// ============================================================================

#[test]
fn test_rkf45_properties() {
    let x0 = DVector::from_vec(vec![1.0]);
    let solver = RKF45::new(x0);
    assert_eq!(solver.order(), 5);
    assert_eq!(solver.stages(), 6);
    assert!(solver.is_adaptive());
    assert!(solver.is_explicit());
}

#[test]
fn test_rkck54_properties() {
    let x0 = DVector::from_vec(vec![1.0]);
    let solver = RKCK54::new(x0);
    assert_eq!(solver.order(), 5);
    assert_eq!(solver.stages(), 6);
    assert!(solver.is_adaptive());
    assert!(solver.is_explicit());
}

#[test]
fn test_rkv65_properties() {
    let x0 = DVector::from_vec(vec![1.0]);
    let solver = RKV65::new(x0);
    assert_eq!(solver.order(), 6);
    assert_eq!(solver.stages(), 9);
    assert!(solver.is_adaptive());
    assert!(solver.is_explicit());
}

#[test]
fn test_rkf78_properties() {
    let x0 = DVector::from_vec(vec![1.0]);
    let solver = RKF78::new(x0);
    assert_eq!(solver.order(), 7);
    assert_eq!(solver.stages(), 13);
    assert!(solver.is_adaptive());
    assert!(solver.is_explicit());
}

#[test]
fn test_rkdp87_properties() {
    let x0 = DVector::from_vec(vec![1.0]);
    let solver = RKDP87::new(x0);
    assert_eq!(solver.order(), 8);
    assert_eq!(solver.stages(), 13);
    assert!(solver.is_adaptive());
    assert!(solver.is_explicit());
}

#[test]
fn test_ssprk34_properties() {
    let x0 = DVector::from_vec(vec![1.0]);
    let solver = SSPRK34::new(x0);
    assert_eq!(solver.order(), 3);
    assert_eq!(solver.stages(), 4);
    assert!(!solver.is_adaptive());
    assert!(solver.is_explicit());
}

// ============================================================================
// Comparison Tests
// ============================================================================

#[test]
fn test_all_rk_exponential_decay() {
    // Test equation: dx/dt = -x, x(0) = 1
    // Exact solution: x(t) = exp(-t)

    let x0 = DVector::from_vec(vec![1.0]);
    let dt = 0.1;
    let t_final = 1.0;
    let n_steps = (t_final / dt) as usize;
    let exact = (-t_final as f64).exp();

    // Test RKF21
    {
        let mut solver = RKF21::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..3 {
                solver.step(|x, _t| -x, dt);
            }
        }
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-3);
    }

    // Test RKBS32
    {
        let mut solver = RKBS32::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..4 {
                solver.step(|x, _t| -x, dt);
            }
        }
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-4);
    }

    // Test RKF45
    {
        let mut solver = RKF45::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..6 {
                solver.step(|x, _t| -x, dt);
            }
        }
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-6);
    }

    // Test RKCK54
    {
        let mut solver = RKCK54::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..6 {
                solver.step(|x, _t| -x, dt);
            }
        }
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-7);
    }

    // Test RKV65
    {
        let mut solver = RKV65::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..9 {
                solver.step(|x, _t| -x, dt);
            }
        }
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-9);
    }

    // Test RKF78
    {
        let mut solver = RKF78::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..13 {
                solver.step(|x, _t| -x, dt);
            }
        }
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-11);
    }

    // Test RKDP87
    {
        let mut solver = RKDP87::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..13 {
                solver.step(|x, _t| -x, dt);
            }
        }
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-12);
    }

    // Test SSPRK34
    {
        let mut solver = SSPRK34::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..4 {
                solver.step(|x, _t| -x, dt);
            }
        }
        assert_relative_eq!(solver.state()[0], exact, epsilon = 1e-4);
    }
}

#[test]
fn test_all_rk_harmonic_oscillator() {
    // d²x/dt² = -x => [x, v]' = [v, -x]
    // Exact: x(t) = cos(t), v(t) = -sin(t)

    let x0 = DVector::from_vec(vec![1.0, 0.0]);
    let dt = 0.05;
    let t_final = 2.0 * std::f64::consts::PI;
    let n_steps = (t_final / dt) as usize;

    let rhs = |x: &DVector<f64>, _t: f64| {
        let mut dxdt = DVector::zeros(2);
        dxdt[0] = x[1];
        dxdt[1] = -x[0];
        dxdt
    };

    // Test RK4
    {
        let mut solver = RK4::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..4 {
                solver.step(&rhs, dt);
            }
        }
        assert_relative_eq!(solver.state()[0], 1.0, epsilon = 1e-3);
        assert_relative_eq!(solver.state()[1], 0.0, epsilon = 1e-3);
    }

    // Test RKDP54
    {
        let mut solver = RKDP54::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..7 {
                solver.step(&rhs, dt);
            }
        }
        assert_relative_eq!(solver.state()[0], 1.0, epsilon = 1e-4);
        assert_relative_eq!(solver.state()[1], 0.0, epsilon = 1e-4);
    }

    // Test RKF45
    {
        let mut solver = RKF45::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..6 {
                solver.step(&rhs, dt);
            }
        }
        assert_relative_eq!(solver.state()[0], 1.0, epsilon = 1e-4);
        assert_relative_eq!(solver.state()[1], 0.0, epsilon = 1e-4);
    }

    // Test RKCK54
    {
        let mut solver = RKCK54::new(x0.clone());
        for _ in 0..n_steps {
            solver.buffer(dt);
            for _ in 0..6 {
                solver.step(&rhs, dt);
            }
        }
        assert_relative_eq!(solver.state()[0], 1.0, epsilon = 1e-4);
        assert_relative_eq!(solver.state()[1], 0.0, epsilon = 1e-4);
    }
}
