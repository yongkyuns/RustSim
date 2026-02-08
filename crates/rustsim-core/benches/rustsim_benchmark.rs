//! Comprehensive RustSim Benchmark Suite
//!
//! This benchmark tests various ODE solvers on computationally challenging systems:
//! 1. Robertson's stiff chemical kinetics problem
//! 2. Van der Pol oscillator with large mu (stiff)
//! 3. Large coupled oscillator system (20 masses)
//! 4. Lorenz chaotic attractor (long-time integration)
//!
//! Each system is solved with multiple solvers to compare performance and accuracy.

use nalgebra::DVector;
use rustsim::solvers::{ExplicitSolver, RK4, RKDP54, RKF45, Euler, RKCK54};
use std::time::Instant;

// ============================================================================
// PROBLEM 1: Robertson's Stiff Chemical Kinetics
// ============================================================================
// dy1/dt = -0.04*y1 + 1e4*y2*y3
// dy2/dt = 0.04*y1 - 1e4*y2*y3 - 3e7*y2^2
// dy3/dt = 3e7*y2^2
//
// Initial: y1=1, y2=0, y3=0
// This is a classic stiff problem with vastly different time scales
// ============================================================================

fn robertson_rhs(state: &DVector<f64>, _t: f64) -> DVector<f64> {
    let y1 = state[0];
    let y2 = state[1];
    let y3 = state[2];

    let mut dydt = DVector::zeros(3);
    dydt[0] = -0.04 * y1 + 1e4 * y2 * y3;
    dydt[1] = 0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2 * y2;
    dydt[2] = 3e7 * y2 * y2;

    dydt
}

// ============================================================================
// PROBLEM 2: Van der Pol Oscillator with Large Mu (Stiff)
// ============================================================================
// dx/dt = y
// dy/dt = mu*(1 - x^2)*y - x
//
// With mu=1000, this becomes very stiff
// ============================================================================

fn van_der_pol_rhs(state: &DVector<f64>, _t: f64, mu: f64) -> DVector<f64> {
    let x = state[0];
    let y = state[1];

    let mut dydt = DVector::zeros(2);
    dydt[0] = y;
    dydt[1] = mu * (1.0 - x * x) * y - x;

    dydt
}

// ============================================================================
// PROBLEM 3: Large Coupled Oscillator System
// ============================================================================
// A chain of 20 masses connected by springs
// Each mass has position x_i and velocity v_i
// State dimension: 40 (20 positions + 20 velocities)
//
// Equations:
// dx_i/dt = v_i
// dv_i/dt = k*(x_{i-1} - 2*x_i + x_{i+1}) - c*v_i
//
// Boundary conditions: x_0 = x_{N+1} = 0 (fixed ends)
// ============================================================================

fn coupled_oscillators_rhs(state: &DVector<f64>, _t: f64) -> DVector<f64> {
    let n = 20; // Number of masses
    let k = 1.0; // Spring constant
    let c = 0.1; // Damping coefficient

    let mut dydt = DVector::zeros(2 * n);

    // First half: velocities (dx_i/dt = v_i)
    for i in 0..n {
        dydt[i] = state[n + i]; // v_i
    }

    // Second half: accelerations (dv_i/dt)
    for i in 0..n {
        let x_left = if i == 0 { 0.0 } else { state[i - 1] };
        let x_curr = state[i];
        let x_right = if i == n - 1 { 0.0 } else { state[i + 1] };
        let v_curr = state[n + i];

        // Spring force: k*(x_{i-1} - 2*x_i + x_{i+1})
        // Damping force: -c*v_i
        dydt[n + i] = k * (x_left - 2.0 * x_curr + x_right) - c * v_curr;
    }

    dydt
}

// ============================================================================
// PROBLEM 4: Lorenz Chaotic Attractor
// ============================================================================
// dx/dt = sigma*(y - x)
// dy/dt = x*(rho - z) - y
// dz/dt = x*y - beta*z
//
// Standard parameters: sigma=10, rho=28, beta=8/3
// ============================================================================

fn lorenz_rhs(state: &DVector<f64>, _t: f64) -> DVector<f64> {
    let sigma = 10.0;
    let rho = 28.0;
    let beta = 8.0 / 3.0;

    let x = state[0];
    let y = state[1];
    let z = state[2];

    let mut dydt = DVector::zeros(3);
    dydt[0] = sigma * (y - x);
    dydt[1] = x * (rho - z) - y;
    dydt[2] = x * y - beta * z;

    dydt
}

// ============================================================================
// Benchmark Runner Functions
// ============================================================================

fn run_solver<S>(
    solver: &mut S,
    rhs: impl Fn(&DVector<f64>, f64) -> DVector<f64>,
    dt: f64,
    t_final: f64,
    stages: usize,
) -> (f64, DVector<f64>)
where
    S: ExplicitSolver,
{
    let n_steps = (t_final / dt) as usize;
    let start = Instant::now();

    for _ in 0..n_steps {
        solver.buffer(dt);

        // Execute all stages of the solver
        for _ in 0..stages {
            solver.step(&rhs, dt);
        }
    }

    let elapsed = start.elapsed().as_secs_f64();
    let final_state = solver.state().clone();

    (elapsed, final_state)
}

// ============================================================================
// BENCHMARK 1: Robertson's Stiff Problem
// ============================================================================

fn benchmark_robertson() {
    println!("\n{:=<70}", "");
    println!("BENCHMARK 1: Robertson's Stiff Chemical Kinetics");
    println!("{:=<70}", "");
    println!("Problem: 3 equations with vastly different time scales");
    println!("Stiffness: Very high (ratio ~1e7)");
    println!();

    let initial = DVector::from_vec(vec![1.0, 0.0, 0.0]);
    let dt = 1e-4; // Small timestep required for stability
    let t_final = 0.1; // Even short integration is challenging

    println!("Time step: {:.6}, Final time: {}", dt, t_final);
    println!("Number of steps: {}", (t_final / dt) as usize);
    println!();

    // Test with different solvers
    let solvers: Vec<(&str, usize)> = vec![
        ("Euler", 1),
        ("RK4", 4),
        ("RKF45", 6),
        ("RKCK54", 6),
        ("RKDP54", 7),
    ];

    for (solver_name, stages) in solvers {
        match solver_name {
            "Euler" => {
                let mut solver = Euler::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, robertson_rhs, dt, t_final, stages);
                println!("{:<12} Time: {:.4}s | Final state: [{:.6e}, {:.6e}, {:.6e}]",
                    solver_name, time, final_state[0], final_state[1], final_state[2]);
            }
            "RK4" => {
                let mut solver = RK4::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, robertson_rhs, dt, t_final, stages);
                println!("{:<12} Time: {:.4}s | Final state: [{:.6e}, {:.6e}, {:.6e}]",
                    solver_name, time, final_state[0], final_state[1], final_state[2]);
            }
            "RKF45" => {
                let mut solver = RKF45::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, robertson_rhs, dt, t_final, stages);
                println!("{:<12} Time: {:.4}s | Final state: [{:.6e}, {:.6e}, {:.6e}]",
                    solver_name, time, final_state[0], final_state[1], final_state[2]);
            }
            "RKCK54" => {
                let mut solver = RKCK54::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, robertson_rhs, dt, t_final, stages);
                println!("{:<12} Time: {:.4}s | Final state: [{:.6e}, {:.6e}, {:.6e}]",
                    solver_name, time, final_state[0], final_state[1], final_state[2]);
            }
            "RKDP54" => {
                let mut solver = RKDP54::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, robertson_rhs, dt, t_final, stages);
                println!("{:<12} Time: {:.4}s | Final state: [{:.6e}, {:.6e}, {:.6e}]",
                    solver_name, time, final_state[0], final_state[1], final_state[2]);
            }
            _ => {}
        }
    }
}

// ============================================================================
// BENCHMARK 2: Van der Pol Oscillator (Stiff)
// ============================================================================

fn benchmark_van_der_pol() {
    println!("\n{:=<70}", "");
    println!("BENCHMARK 2: Van der Pol Oscillator (mu=100, stiff)");
    println!("{:=<70}", "");
    println!("Problem: 2 equations with high stiffness parameter");
    println!("Stiffness: High (mu=100)");
    println!();

    let initial = DVector::from_vec(vec![2.0, 0.0]);
    let mu = 100.0;
    let dt = 1e-3;
    let t_final = 1.0;

    println!("Time step: {:.6}, Final time: {}", dt, t_final);
    println!("Number of steps: {}", (t_final / dt) as usize);
    println!();

    let solvers: Vec<(&str, usize)> = vec![
        ("Euler", 1),
        ("RK4", 4),
        ("RKF45", 6),
        ("RKCK54", 6),
        ("RKDP54", 7),
    ];

    for (solver_name, stages) in solvers {
        let rhs = |state: &DVector<f64>, t: f64| van_der_pol_rhs(state, t, mu);

        match solver_name {
            "Euler" => {
                let mut solver = Euler::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, rhs, dt, t_final, stages);
                println!("{:<12} Time: {:.4}s | Final state: [{:.6}, {:.6}]",
                    solver_name, time, final_state[0], final_state[1]);
            }
            "RK4" => {
                let mut solver = RK4::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, rhs, dt, t_final, stages);
                println!("{:<12} Time: {:.4}s | Final state: [{:.6}, {:.6}]",
                    solver_name, time, final_state[0], final_state[1]);
            }
            "RKF45" => {
                let mut solver = RKF45::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, rhs, dt, t_final, stages);
                println!("{:<12} Time: {:.4}s | Final state: [{:.6}, {:.6}]",
                    solver_name, time, final_state[0], final_state[1]);
            }
            "RKCK54" => {
                let mut solver = RKCK54::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, rhs, dt, t_final, stages);
                println!("{:<12} Time: {:.4}s | Final state: [{:.6}, {:.6}]",
                    solver_name, time, final_state[0], final_state[1]);
            }
            "RKDP54" => {
                let mut solver = RKDP54::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, rhs, dt, t_final, stages);
                println!("{:<12} Time: {:.4}s | Final state: [{:.6}, {:.6}]",
                    solver_name, time, final_state[0], final_state[1]);
            }
            _ => {}
        }
    }
}

// ============================================================================
// BENCHMARK 3: Large Coupled Oscillator System
// ============================================================================

fn benchmark_coupled_oscillators() {
    println!("\n{:=<70}", "");
    println!("BENCHMARK 3: Coupled Oscillator System (20 masses)");
    println!("{:=<70}", "");
    println!("Problem: 40 equations (20 positions + 20 velocities)");
    println!("Stiffness: Moderate");
    println!();

    // Initial condition: sinusoidal displacement
    let n = 20;
    let mut initial = DVector::zeros(2 * n);
    for i in 0..n {
        // Initial position: sine wave
        initial[i] = (std::f64::consts::PI * (i as f64) / (n as f64)).sin();
        // Initial velocity: zero
        initial[n + i] = 0.0;
    }

    let dt = 0.01;
    let t_final = 10.0; // Long integration

    println!("Time step: {:.6}, Final time: {}", dt, t_final);
    println!("Number of steps: {}", (t_final / dt) as usize);
    println!();

    let solvers: Vec<(&str, usize)> = vec![
        ("Euler", 1),
        ("RK4", 4),
        ("RKF45", 6),
        ("RKCK54", 6),
        ("RKDP54", 7),
    ];

    for (solver_name, stages) in solvers {
        match solver_name {
            "Euler" => {
                let mut solver = Euler::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, coupled_oscillators_rhs, dt, t_final, stages);
                let energy = calculate_oscillator_energy(&final_state);
                println!("{:<12} Time: {:.4}s | Energy: {:.6} | x[0]={:.6}, x[10]={:.6}",
                    solver_name, time, energy, final_state[0], final_state[10]);
            }
            "RK4" => {
                let mut solver = RK4::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, coupled_oscillators_rhs, dt, t_final, stages);
                let energy = calculate_oscillator_energy(&final_state);
                println!("{:<12} Time: {:.4}s | Energy: {:.6} | x[0]={:.6}, x[10]={:.6}",
                    solver_name, time, energy, final_state[0], final_state[10]);
            }
            "RKF45" => {
                let mut solver = RKF45::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, coupled_oscillators_rhs, dt, t_final, stages);
                let energy = calculate_oscillator_energy(&final_state);
                println!("{:<12} Time: {:.4}s | Energy: {:.6} | x[0]={:.6}, x[10]={:.6}",
                    solver_name, time, energy, final_state[0], final_state[10]);
            }
            "RKCK54" => {
                let mut solver = RKCK54::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, coupled_oscillators_rhs, dt, t_final, stages);
                let energy = calculate_oscillator_energy(&final_state);
                println!("{:<12} Time: {:.4}s | Energy: {:.6} | x[0]={:.6}, x[10]={:.6}",
                    solver_name, time, energy, final_state[0], final_state[10]);
            }
            "RKDP54" => {
                let mut solver = RKDP54::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, coupled_oscillators_rhs, dt, t_final, stages);
                let energy = calculate_oscillator_energy(&final_state);
                println!("{:<12} Time: {:.4}s | Energy: {:.6} | x[0]={:.6}, x[10]={:.6}",
                    solver_name, time, energy, final_state[0], final_state[10]);
            }
            _ => {}
        }
    }
}

fn calculate_oscillator_energy(state: &DVector<f64>) -> f64 {
    let n = state.len() / 2;
    let k = 1.0;
    let m = 1.0;

    let mut energy = 0.0;

    // Kinetic energy
    for i in 0..n {
        energy += 0.5 * m * state[n + i].powi(2);
    }

    // Potential energy
    for i in 0..n {
        let x_left = if i == 0 { 0.0 } else { state[i - 1] };
        let x_curr = state[i];
        energy += 0.5 * k * (x_curr - x_left).powi(2);
    }

    energy
}

// ============================================================================
// BENCHMARK 4: Lorenz Chaotic System (Long Integration)
// ============================================================================

fn benchmark_lorenz() {
    println!("\n{:=<70}", "");
    println!("BENCHMARK 4: Lorenz Chaotic Attractor (Long Integration)");
    println!("{:=<70}", "");
    println!("Problem: 3 equations, chaotic dynamics");
    println!("Stiffness: Low (but sensitive to errors)");
    println!();

    let initial = DVector::from_vec(vec![1.0, 1.0, 1.0]);
    let dt = 0.001;
    let t_final = 50.0; // Long integration to accumulate errors

    println!("Time step: {:.6}, Final time: {}", dt, t_final);
    println!("Number of steps: {}", (t_final / dt) as usize);
    println!();

    let solvers: Vec<(&str, usize)> = vec![
        ("Euler", 1),
        ("RK4", 4),
        ("RKF45", 6),
        ("RKCK54", 6),
        ("RKDP54", 7),
    ];

    for (solver_name, stages) in solvers {
        match solver_name {
            "Euler" => {
                let mut solver = Euler::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, lorenz_rhs, dt, t_final, stages);
                println!("{:<12} Time: {:.4}s | Final state: [{:.6}, {:.6}, {:.6}]",
                    solver_name, time, final_state[0], final_state[1], final_state[2]);
            }
            "RK4" => {
                let mut solver = RK4::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, lorenz_rhs, dt, t_final, stages);
                println!("{:<12} Time: {:.4}s | Final state: [{:.6}, {:.6}, {:.6}]",
                    solver_name, time, final_state[0], final_state[1], final_state[2]);
            }
            "RKF45" => {
                let mut solver = RKF45::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, lorenz_rhs, dt, t_final, stages);
                println!("{:<12} Time: {:.4}s | Final state: [{:.6}, {:.6}, {:.6}]",
                    solver_name, time, final_state[0], final_state[1], final_state[2]);
            }
            "RKCK54" => {
                let mut solver = RKCK54::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, lorenz_rhs, dt, t_final, stages);
                println!("{:<12} Time: {:.4}s | Final state: [{:.6}, {:.6}, {:.6}]",
                    solver_name, time, final_state[0], final_state[1], final_state[2]);
            }
            "RKDP54" => {
                let mut solver = RKDP54::new(initial.clone());
                let (time, final_state) = run_solver(&mut solver, lorenz_rhs, dt, t_final, stages);
                println!("{:<12} Time: {:.4}s | Final state: [{:.6}, {:.6}, {:.6}]",
                    solver_name, time, final_state[0], final_state[1], final_state[2]);
            }
            _ => {}
        }
    }
}

// ============================================================================
// Main Benchmark Runner
// ============================================================================

fn main() {
    println!("\n{:#<70}", "");
    println!("# RustSim Comprehensive Benchmark Suite");
    println!("{:#<70}", "");
    println!();
    println!("Testing ODE solvers on challenging dynamical systems:");
    println!("  1. Robertson's stiff chemical kinetics");
    println!("  2. Van der Pol oscillator (stiff, mu=100)");
    println!("  3. Coupled oscillator system (20 masses, 40 equations)");
    println!("  4. Lorenz chaotic attractor (long integration)");
    println!();
    println!("Solvers tested:");
    println!("  - Euler:  Order 1, 1 stage  (baseline)");
    println!("  - RK4:    Order 4, 4 stages (classic)");
    println!("  - RKF45:  Order 5, 6 stages (adaptive)");
    println!("  - RKCK54: Order 5, 6 stages (adaptive)");
    println!("  - RKDP54: Order 5, 7 stages (adaptive, MATLAB's ode45)");
    println!();

    let total_start = Instant::now();

    // Run all benchmarks
    benchmark_robertson();
    benchmark_van_der_pol();
    benchmark_coupled_oscillators();
    benchmark_lorenz();

    let total_time = total_start.elapsed().as_secs_f64();

    println!("\n{:=<70}", "");
    println!("BENCHMARK COMPLETE");
    println!("{:=<70}", "");
    println!("Total benchmark time: {:.2}s", total_time);
    println!();
}
