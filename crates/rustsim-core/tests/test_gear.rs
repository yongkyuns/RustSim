//! Comprehensive tests for GEAR family of solvers
//!
//! Tests match PathSim's test_gear.py exactly to ensure test parity

use approx::assert_relative_eq;
use nalgebra::DVector;
use rustsim::solvers::{
    compute_bdf_coefficients, ImplicitSolver, Solver, GEAR21, GEAR32, GEAR43, GEAR52A, GEAR54,
};

// ===========================================================================================
// REFERENCE PROBLEMS (matching PathSim)
// ===========================================================================================

struct Problem {
    name: &'static str,
    x0: f64,
    t_span: (f64, f64),
}

impl Problem {
    fn func(&self, x: &DVector<f64>, _t: f64) -> DVector<f64> {
        match self.name {
            "exponential_decay" => -x,
            "logistic" => x.map(|xi| xi * (1.0 - xi)),
            "quadratic" => x.map(|xi| xi * xi),
            "sin_decay" => x.map(|xi| -xi * _t.sin()),
            "polynomial" => {
                let mut result = DVector::zeros(x.len());
                result[0] = _t * _t - x[0];
                result
            }
            _ => panic!("Unknown problem"),
        }
    }

    fn jac(&self, x: &DVector<f64>, _t: f64) -> nalgebra::DMatrix<f64> {
        let mut j = nalgebra::DMatrix::zeros(x.len(), x.len());
        match self.name {
            "exponential_decay" => {
                j[(0, 0)] = -1.0;
            }
            "logistic" => {
                j[(0, 0)] = 1.0 - 2.0 * x[0];
            }
            "quadratic" => {
                j[(0, 0)] = 2.0 * x[0];
            }
            "sin_decay" => {
                j[(0, 0)] = -_t.sin();
            }
            "polynomial" => {
                j[(0, 0)] = -1.0;
            }
            _ => panic!("Unknown problem"),
        }
        j
    }

    fn solution(&self, t: f64) -> f64 {
        match self.name {
            "exponential_decay" => (-t).exp(),
            "logistic" => 1.0 / (1.0 + (-t).exp()),
            "quadratic" => 1.0 / (1.0 - t),
            "sin_decay" => (t.cos() - 1.0).exp(),
            "polynomial" => t * t - 2.0 * t + 2.0 - 2.0 * (-t).exp(),
            _ => panic!("Unknown problem"),
        }
    }
}

const PROBLEMS: &[Problem] = &[
    Problem {
        name: "exponential_decay",
        x0: 1.0,
        t_span: (0.0, 5.0),
    },
    Problem {
        name: "logistic",
        x0: 0.5,
        t_span: (0.0, 10.0),
    },
    Problem {
        name: "quadratic",
        x0: 1.0,
        t_span: (0.0, 0.6),
    },
    Problem {
        name: "sin_decay",
        x0: 1.0,
        t_span: (0.0, 10.0),
    },
    Problem {
        name: "polynomial",
        x0: 0.0,
        t_span: (0.0, 5.0),
    },
];

// ===========================================================================================
// HELPER FUNCTIONS
// ===========================================================================================

fn integrate_fixed<S: ImplicitSolver>(
    solver: &mut S,
    problem: &Problem,
    dt: f64,
) -> (Vec<f64>, Vec<f64>) {
    let mut time = problem.t_span.0;
    let t_end = problem.t_span.1;

    let mut times = vec![time];
    let mut states = vec![problem.x0];

    while time < t_end {
        solver.buffer(dt);

        // Solve implicit equation
        let problem_name = problem.name;
        let _ = solver.solve(
            |x: &DVector<f64>, t: f64| -> DVector<f64> {
                match problem_name {
                    "exponential_decay" => -x,
                    "logistic" => x.map(|xi| xi * (1.0 - xi)),
                    "quadratic" => x.map(|xi| xi * xi),
                    "sin_decay" => x.map(|xi| -xi * t.sin()),
                    "polynomial" => {
                        let mut result = DVector::zeros(x.len());
                        result[0] = t * t - x[0];
                        result
                    }
                    _ => panic!("Unknown problem"),
                }
            },
            Some(|x: &DVector<f64>, t: f64| -> nalgebra::DMatrix<f64> {
                let mut j = nalgebra::DMatrix::zeros(x.len(), x.len());
                match problem_name {
                    "exponential_decay" => {
                        j[(0, 0)] = -1.0;
                    }
                    "logistic" => {
                        j[(0, 0)] = 1.0 - 2.0 * x[0];
                    }
                    "quadratic" => {
                        j[(0, 0)] = 2.0 * x[0];
                    }
                    "sin_decay" => {
                        j[(0, 0)] = -t.sin();
                    }
                    "polynomial" => {
                        j[(0, 0)] = -1.0;
                    }
                    _ => panic!("Unknown problem"),
                }
                j
            }),
            dt,
        );

        // Finalize step (but don't use adaptive logic)
        let _ = solver.step(
            |x: &DVector<f64>, t: f64| -> DVector<f64> {
                match problem_name {
                    "exponential_decay" => -x,
                    "logistic" => x.map(|xi| xi * (1.0 - xi)),
                    "quadratic" => x.map(|xi| xi * xi),
                    "sin_decay" => x.map(|xi| -xi * t.sin()),
                    "polynomial" => {
                        let mut result = DVector::zeros(x.len());
                        result[0] = t * t - x[0];
                        result
                    }
                    _ => panic!("Unknown problem"),
                }
            },
            dt,
        );

        time += dt;
        times.push(time);
        states.push(solver.state()[0]);
    }

    (times, states)
}

fn integrate_adaptive<S: ImplicitSolver>(
    solver: &mut S,
    problem: &Problem,
    dt_max: f64,
    _tol_abs: f64,
) -> (Vec<f64>, Vec<f64>) {
    let mut time = problem.t_span.0;
    let t_end = problem.t_span.1;
    let dt_min = 1e-10;

    let mut times = vec![time];
    let mut states = vec![problem.x0];
    let mut dt = dt_max.min(t_end - time);

    let problem_name = problem.name;

    while time < t_end {
        solver.buffer(dt);

        // Solve implicit equation
        let solve_result = solver.solve(
            |x: &DVector<f64>, t: f64| -> DVector<f64> {
                match problem_name {
                    "exponential_decay" => -x,
                    "logistic" => x.map(|xi| xi * (1.0 - xi)),
                    "quadratic" => x.map(|xi| xi * xi),
                    "sin_decay" => x.map(|xi| -xi * t.sin()),
                    "polynomial" => {
                        let mut result = DVector::zeros(x.len());
                        result[0] = t * t - x[0];
                        result
                    }
                    _ => panic!("Unknown problem"),
                }
            },
            Some(|x: &DVector<f64>, t: f64| -> nalgebra::DMatrix<f64> {
                let mut j = nalgebra::DMatrix::zeros(x.len(), x.len());
                match problem_name {
                    "exponential_decay" => {
                        j[(0, 0)] = -1.0;
                    }
                    "logistic" => {
                        j[(0, 0)] = 1.0 - 2.0 * x[0];
                    }
                    "quadratic" => {
                        j[(0, 0)] = 2.0 * x[0];
                    }
                    "sin_decay" => {
                        j[(0, 0)] = -t.sin();
                    }
                    "polynomial" => {
                        j[(0, 0)] = -1.0;
                    }
                    _ => panic!("Unknown problem"),
                }
                j
            }),
            dt,
        );

        if solve_result.is_err() {
            // Convergence failed, reduce timestep and retry
            let _ = solver.revert();
            dt *= 0.5;
            continue;
        }

        // Finalize step and get error estimate
        let result = solver.step(
            |x: &DVector<f64>, t: f64| -> DVector<f64> {
                match problem_name {
                    "exponential_decay" => -x,
                    "logistic" => x.map(|xi| xi * (1.0 - xi)),
                    "quadratic" => x.map(|xi| xi * xi),
                    "sin_decay" => x.map(|xi| -xi * t.sin()),
                    "polynomial" => {
                        let mut result = DVector::zeros(x.len());
                        result[0] = t * t - x[0];
                        result
                    }
                    _ => panic!("Unknown problem"),
                }
            },
            dt,
        );

        if result.success {
            // Step accepted
            time += dt;
            times.push(time);
            states.push(solver.state()[0]);

            // Adjust timestep for next step
            if let Some(scale) = result.scale {
                dt = (scale * dt).clamp(dt_min, dt_max);
            }

            // Don't overshoot
            if time + dt > t_end {
                dt = t_end - time;
            }
        } else {
            // Step rejected, reduce timestep and retry
            let _ = solver.revert();
            if let Some(scale) = result.scale {
                dt = (scale * dt).max(dt_min);
            } else {
                dt *= 0.5;
            }
        }
    }

    (times, states)
}

// ===========================================================================================
// TEST: BDF COEFFICIENTS
// ===========================================================================================

#[test]
fn test_gear_bdf_coefficients() {
    // Test BDF coefficient computation for constant timesteps

    // Order 1
    {
        let (beta, alpha) = compute_bdf_coefficients(1, &[1.0]);
        assert_relative_eq!(beta, 1.0, epsilon = 1e-7);
        assert_eq!(alpha.len(), 1);
        assert_relative_eq!(alpha[0], 1.0, epsilon = 1e-7);
    }

    // Order 2
    {
        let (beta, alpha) = compute_bdf_coefficients(2, &[1.0, 1.0]);
        assert_relative_eq!(beta, 2.0 / 3.0, epsilon = 1e-7);
        assert_eq!(alpha.len(), 2);
        assert_relative_eq!(alpha[0], 4.0 / 3.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[1], -1.0 / 3.0, epsilon = 1e-7);
    }

    // Order 3
    {
        let (beta, alpha) = compute_bdf_coefficients(3, &[1.0, 1.0, 1.0]);
        assert_relative_eq!(beta, 6.0 / 11.0, epsilon = 1e-7);
        assert_eq!(alpha.len(), 3);
        assert_relative_eq!(alpha[0], 18.0 / 11.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[1], -9.0 / 11.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[2], 2.0 / 11.0, epsilon = 1e-7);
    }

    // Order 4
    {
        let (beta, alpha) = compute_bdf_coefficients(4, &[1.0, 1.0, 1.0, 1.0]);
        assert_relative_eq!(beta, 12.0 / 25.0, epsilon = 1e-7);
        assert_eq!(alpha.len(), 4);
        assert_relative_eq!(alpha[0], 48.0 / 25.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[1], -36.0 / 25.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[2], 16.0 / 25.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[3], -3.0 / 25.0, epsilon = 1e-7);
    }

    // Order 5
    {
        let (beta, alpha) = compute_bdf_coefficients(5, &[1.0, 1.0, 1.0, 1.0, 1.0]);
        assert_relative_eq!(beta, 60.0 / 137.0, epsilon = 1e-7);
        assert_eq!(alpha.len(), 5);
        assert_relative_eq!(alpha[0], 300.0 / 137.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[1], -300.0 / 137.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[2], 200.0 / 137.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[3], -75.0 / 137.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[4], 12.0 / 137.0, epsilon = 1e-7);
    }

    // Order 6
    {
        let (beta, alpha) = compute_bdf_coefficients(6, &[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);
        assert_relative_eq!(beta, 60.0 / 147.0, epsilon = 1e-7);
        assert_eq!(alpha.len(), 6);
        assert_relative_eq!(alpha[0], 360.0 / 147.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[1], -450.0 / 147.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[2], 400.0 / 147.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[3], -225.0 / 147.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[4], 72.0 / 147.0, epsilon = 1e-7);
        assert_relative_eq!(alpha[5], -10.0 / 147.0, epsilon = 1e-7);
    }
}

#[test]
fn test_gear_history_buffer() {
    // Test that GEAR solvers maintain correct history buffer length

    // GEAR21 should maintain history of 2
    {
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = GEAR21::new(x0);

        for _k in 0..10 {
            solver.buffer(1.0);
            // Note: history is internal, can't directly test length
            // But we can verify no panic occurs
        }

        assert_eq!(solver.order(), 2);
    }

    // GEAR32 should maintain history of 3
    {
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = GEAR32::new(x0);

        for _k in 0..10 {
            solver.buffer(1.0);
        }

        assert_eq!(solver.order(), 3);
    }

    // GEAR54 should maintain history of 5
    {
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = GEAR54::new(x0);

        for _k in 0..10 {
            solver.buffer(1.0);
        }

        assert_eq!(solver.order(), 5);
    }
}

// ===========================================================================================
// TEST: GEAR21 CONVERGENCE
// ===========================================================================================

#[test]
fn test_gear21_convergence() {
    // Test that GEAR21 shows monotonic error decrease and correct convergence order

    for problem in PROBLEMS {
        let divisions = [16.0, 32.0, 64.0, 128.0, 256.0];
        let duration = problem.t_span.1 - problem.t_span.0;
        let mut errors = Vec::new();
        let mut timesteps = Vec::new();

        for &div in &divisions {
            let dt = duration / div;
            let x0 = DVector::from_vec(vec![problem.x0]);
            let mut solver = GEAR21::new(x0);

            let (times, states) = integrate_fixed(&mut solver, problem, dt);

            // Compute error
            let mut err_sum = 0.0;
            for (t, &x_num) in times.iter().zip(states.iter()) {
                let x_exact = problem.solution(*t);
                err_sum += (x_num - x_exact).abs();
            }
            let mean_error = err_sum / states.len() as f64;

            errors.push(mean_error);
            timesteps.push(dt);
        }

        // Check monotonic decrease
        for i in 1..errors.len() {
            assert!(
                errors[i] < errors[i - 1],
                "Error not decreasing for {} at step {}",
                problem.name,
                i
            );
        }

        // Check convergence order (expect > 1.5 due to startup phase)
        let log_h: Vec<f64> = timesteps.iter().map(|&h| h.ln()).collect();
        let log_e: Vec<f64> = errors.iter().map(|&e| e.ln()).collect();

        // Linear regression
        let n = log_h.len() as f64;
        let sum_h: f64 = log_h.iter().sum();
        let sum_e: f64 = log_e.iter().sum();
        let sum_hh: f64 = log_h.iter().map(|&h| h * h).sum();
        let sum_he: f64 = log_h.iter().zip(&log_e).map(|(&h, &e)| h * e).sum();

        let slope = (n * sum_he - sum_h * sum_e) / (n * sum_hh - sum_h * sum_h);

        assert!(
            slope > 1.5,
            "Convergence order {} < 1.5 for problem {}",
            slope,
            problem.name
        );
    }
}

// ===========================================================================================
// TEST: GEAR21 ADAPTIVE
// ===========================================================================================

#[test]
fn test_gear21_adaptive() {
    // Test adaptive timestepping maintains error within tolerance

    for problem in PROBLEMS {
        let x0 = DVector::from_vec(vec![problem.x0]);
        let tol_abs = 1e-5;
        let mut solver = GEAR21::with_tolerances(x0, tol_abs, 0.0);

        let duration = problem.t_span.1 - problem.t_span.0;
        let (times, states) = integrate_adaptive(&mut solver, problem, duration, tol_abs);

        // Compute mean error
        let mut err_sum = 0.0;
        for (t, &x_num) in times.iter().zip(states.iter()) {
            let x_exact = problem.solution(*t);
            err_sum += x_num - x_exact;
        }
        let mean_error = err_sum / states.len() as f64;

        // Error should be controlled to within 10x tolerance (global error)
        assert!(
            mean_error.abs() < tol_abs * 10.0,
            "Mean error {} exceeds tolerance for {}",
            mean_error,
            problem.name
        );
    }
}

// ===========================================================================================
// TEST: GEAR32 CONVERGENCE
// ===========================================================================================

#[test]
fn test_gear32_convergence() {
    // Test GEAR32 convergence order (should be > 2)

    for problem in PROBLEMS {
        let divisions = [16.0, 32.0, 64.0, 128.0, 256.0];
        let duration = problem.t_span.1 - problem.t_span.0;
        let mut errors = Vec::new();
        let mut timesteps = Vec::new();

        for &div in &divisions {
            let dt = duration / div;
            let x0 = DVector::from_vec(vec![problem.x0]);
            let mut solver = GEAR32::new(x0);

            let (times, states) = integrate_fixed(&mut solver, problem, dt);

            let mut err_sum = 0.0;
            for (t, &x_num) in times.iter().zip(states.iter()) {
                let x_exact = problem.solution(*t);
                err_sum += (x_num - x_exact).abs();
            }
            let mean_error = err_sum / states.len() as f64;

            errors.push(mean_error);
            timesteps.push(dt);
        }

        // Check monotonic decrease
        for i in 1..errors.len() {
            assert!(
                errors[i] < errors[i - 1],
                "Error not decreasing for {}",
                problem.name
            );
        }

        // Check convergence order > 2
        let log_h: Vec<f64> = timesteps.iter().map(|&h| h.ln()).collect();
        let log_e: Vec<f64> = errors.iter().map(|&e| e.ln()).collect();

        let n = log_h.len() as f64;
        let sum_h: f64 = log_h.iter().sum();
        let sum_e: f64 = log_e.iter().sum();
        let sum_hh: f64 = log_h.iter().map(|&h| h * h).sum();
        let sum_he: f64 = log_h.iter().zip(&log_e).map(|(&h, &e)| h * e).sum();

        let slope = (n * sum_he - sum_h * sum_e) / (n * sum_hh - sum_h * sum_h);

        assert!(
            slope > 2.0,
            "Convergence order {} < 2.0 for problem {}",
            slope,
            problem.name
        );
    }
}

// ===========================================================================================
// TEST: GEAR32 ADAPTIVE
// ===========================================================================================

#[test]
fn test_gear32_adaptive() {
    for problem in PROBLEMS {
        let x0 = DVector::from_vec(vec![problem.x0]);
        let tol_abs = 1e-5;
        let mut solver = GEAR32::with_tolerances(x0, tol_abs, 0.0);

        let duration = problem.t_span.1 - problem.t_span.0;
        let (times, states) = integrate_adaptive(&mut solver, problem, duration, tol_abs);

        let mut err_sum = 0.0;
        for (t, &x_num) in times.iter().zip(states.iter()) {
            let x_exact = problem.solution(*t);
            err_sum += x_num - x_exact;
        }
        let mean_error = err_sum / states.len() as f64;

        assert!(
            mean_error.abs() < tol_abs * 10.0,
            "Mean error {} exceeds tolerance for {}",
            mean_error,
            problem.name
        );
    }
}

// ===========================================================================================
// TEST: GEAR43 CONVERGENCE
// ===========================================================================================

#[test]
fn test_gear43_convergence() {
    // Test GEAR43 convergence order (should be > 3)

    for problem in PROBLEMS {
        let divisions = [10.0, 20.0, 40.0, 80.0, 160.0];
        let duration = problem.t_span.1 - problem.t_span.0;
        let mut errors = Vec::new();
        let mut timesteps = Vec::new();

        for &div in &divisions {
            let dt = duration / div;
            let x0 = DVector::from_vec(vec![problem.x0]);
            let mut solver = GEAR43::new(x0);

            let (times, states) = integrate_fixed(&mut solver, problem, dt);

            let mut err_sum = 0.0;
            for (t, &x_num) in times.iter().zip(states.iter()) {
                let x_exact = problem.solution(*t);
                err_sum += (x_num - x_exact).abs();
            }
            let mean_error = err_sum / states.len() as f64;

            errors.push(mean_error);
            timesteps.push(dt);
        }

        // Check monotonic decrease
        for i in 1..errors.len() {
            assert!(
                errors[i] < errors[i - 1],
                "Error not decreasing for {}",
                problem.name
            );
        }

        // Check convergence order > 3
        let log_h: Vec<f64> = timesteps.iter().map(|&h| h.ln()).collect();
        let log_e: Vec<f64> = errors.iter().map(|&e| e.ln()).collect();

        let n = log_h.len() as f64;
        let sum_h: f64 = log_h.iter().sum();
        let sum_e: f64 = log_e.iter().sum();
        let sum_hh: f64 = log_h.iter().map(|&h| h * h).sum();
        let sum_he: f64 = log_h.iter().zip(&log_e).map(|(&h, &e)| h * e).sum();

        let slope = (n * sum_he - sum_h * sum_e) / (n * sum_hh - sum_h * sum_h);

        assert!(
            slope > 3.0,
            "Convergence order {} < 3.0 for problem {}",
            slope,
            problem.name
        );
    }
}

// ===========================================================================================
// TEST: GEAR43 ADAPTIVE
// ===========================================================================================

#[test]
fn test_gear43_adaptive() {
    for problem in PROBLEMS {
        let x0 = DVector::from_vec(vec![problem.x0]);
        let tol_abs = 1e-5;
        let mut solver = GEAR43::with_tolerances(x0, tol_abs, 0.0);

        let duration = problem.t_span.1 - problem.t_span.0;
        let (times, states) = integrate_adaptive(&mut solver, problem, duration, tol_abs);

        let mut err_sum = 0.0;
        for (t, &x_num) in times.iter().zip(states.iter()) {
            let x_exact = problem.solution(*t);
            err_sum += x_num - x_exact;
        }
        let mean_error = err_sum / states.len() as f64;

        assert!(
            mean_error.abs() < tol_abs * 10.0,
            "Mean error {} exceeds tolerance for {}",
            mean_error,
            problem.name
        );
    }
}

// ===========================================================================================
// TEST: GEAR54 CONVERGENCE
// ===========================================================================================

#[test]
fn test_gear54_convergence() {
    // Test GEAR54 convergence order (should be > 3 due to ESDIRK32 startup)

    for problem in PROBLEMS {
        let divisions = [10.0, 20.0, 40.0, 80.0, 160.0];
        let duration = problem.t_span.1 - problem.t_span.0;
        let mut errors = Vec::new();
        let mut timesteps = Vec::new();

        for &div in &divisions {
            let dt = duration / div;
            let x0 = DVector::from_vec(vec![problem.x0]);
            let mut solver = GEAR54::new(x0);

            let (times, states) = integrate_fixed(&mut solver, problem, dt);

            let mut err_sum = 0.0;
            for (t, &x_num) in times.iter().zip(states.iter()) {
                let x_exact = problem.solution(*t);
                err_sum += (x_num - x_exact).abs();
            }
            let mean_error = err_sum / states.len() as f64;

            errors.push(mean_error);
            timesteps.push(dt);
        }

        // Check monotonic decrease
        for i in 1..errors.len() {
            assert!(
                errors[i] < errors[i - 1],
                "Error not decreasing for {}",
                problem.name
            );
        }

        // Check convergence order > 3 (limited by startup)
        let log_h: Vec<f64> = timesteps.iter().map(|&h| h.ln()).collect();
        let log_e: Vec<f64> = errors.iter().map(|&e| e.ln()).collect();

        let n = log_h.len() as f64;
        let sum_h: f64 = log_h.iter().sum();
        let sum_e: f64 = log_e.iter().sum();
        let sum_hh: f64 = log_h.iter().map(|&h| h * h).sum();
        let sum_he: f64 = log_h.iter().zip(&log_e).map(|(&h, &e)| h * e).sum();

        let slope = (n * sum_he - sum_h * sum_e) / (n * sum_hh - sum_h * sum_h);

        assert!(
            slope > 3.0,
            "Convergence order {} < 3.0 for problem {}",
            slope,
            problem.name
        );
    }
}

// ===========================================================================================
// TEST: GEAR54 ADAPTIVE
// ===========================================================================================

#[test]
fn test_gear54_adaptive() {
    for problem in PROBLEMS {
        let x0 = DVector::from_vec(vec![problem.x0]);
        let tol_abs = 1e-5;
        let mut solver = GEAR54::with_tolerances(x0, tol_abs, 0.0);

        let duration = problem.t_span.1 - problem.t_span.0;
        let (times, states) = integrate_adaptive(&mut solver, problem, duration, tol_abs);

        let mut err_sum = 0.0;
        for (t, &x_num) in times.iter().zip(states.iter()) {
            let x_exact = problem.solution(*t);
            err_sum += x_num - x_exact;
        }
        let mean_error = err_sum / states.len() as f64;

        assert!(
            mean_error.abs() < tol_abs * 10.0,
            "Mean error {} exceeds tolerance for {}",
            mean_error,
            problem.name
        );
    }
}

// ===========================================================================================
// TEST: GEAR52A CONVERGENCE
// ===========================================================================================

#[test]
fn test_gear52a_convergence() {
    // Test GEAR52A (variable order) convergence (should be > 2)

    for problem in PROBLEMS {
        let divisions = [10.0, 20.0, 40.0, 80.0, 160.0];
        let duration = problem.t_span.1 - problem.t_span.0;
        let mut errors = Vec::new();
        let mut timesteps = Vec::new();

        for &div in &divisions {
            let dt = duration / div;
            let x0 = DVector::from_vec(vec![problem.x0]);
            let mut solver = GEAR52A::new(x0);

            let (times, states) = integrate_fixed(&mut solver, problem, dt);

            let mut err_sum = 0.0;
            for (t, &x_num) in times.iter().zip(states.iter()) {
                let x_exact = problem.solution(*t);
                err_sum += (x_num - x_exact).abs();
            }
            let mean_error = err_sum / states.len() as f64;

            errors.push(mean_error);
            timesteps.push(dt);
        }

        // Check convergence order > 2 (variable order should achieve this)
        let log_h: Vec<f64> = timesteps.iter().map(|&h| h.ln()).collect();
        let log_e: Vec<f64> = errors.iter().map(|&e| e.ln()).collect();

        let n = log_h.len() as f64;
        let sum_h: f64 = log_h.iter().sum();
        let sum_e: f64 = log_e.iter().sum();
        let sum_hh: f64 = log_h.iter().map(|&h| h * h).sum();
        let sum_he: f64 = log_h.iter().zip(&log_e).map(|(&h, &e)| h * e).sum();

        let slope = (n * sum_he - sum_h * sum_e) / (n * sum_hh - sum_h * sum_h);

        assert!(
            slope > 2.0,
            "Convergence order {} < 2.0 for problem {}",
            slope,
            problem.name
        );
    }
}

// ===========================================================================================
// TEST: GEAR52A ADAPTIVE
// ===========================================================================================

#[test]
fn test_gear52a_adaptive() {
    for problem in PROBLEMS {
        let x0 = DVector::from_vec(vec![problem.x0]);
        let tol_abs = 1e-5;
        let mut solver = GEAR52A::with_tolerances(x0, tol_abs, 0.0);

        let duration = problem.t_span.1 - problem.t_span.0;
        let (times, states) = integrate_adaptive(&mut solver, problem, duration, tol_abs);

        let mut err_sum = 0.0;
        for (t, &x_num) in times.iter().zip(states.iter()) {
            let x_exact = problem.solution(*t);
            err_sum += x_num - x_exact;
        }
        let mean_error = err_sum / states.len() as f64;

        assert!(
            mean_error.abs() < tol_abs * 10.0,
            "Mean error {} exceeds tolerance for {}",
            mean_error,
            problem.name
        );
    }
}
