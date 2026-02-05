//! Basic integration tests for RustSim using only available blocks

use approx::assert_relative_eq;
use nalgebra::DVector;
use rustsim::prelude::*;
use rustsim::solvers::{ExplicitSolver, Solver, RK4};

#[test]
fn test_adder_basic() {
    let mut adder = Adder::<3>::new();

    // Set inputs
    adder.inputs_mut()[0] = 1.0;
    adder.inputs_mut()[1] = 2.0;
    adder.inputs_mut()[2] = 3.0;

    // Update
    adder.update(0.0);

    // Check output
    assert_relative_eq!(*adder.outputs().get(0).unwrap(), 6.0, epsilon = 1e-10);
}

#[test]
fn test_adder_weighted() {
    let mut adder = Adder::<3>::with_weights([1.0, -1.0, 2.0]);

    // Set inputs
    adder.inputs_mut()[0] = 10.0;
    adder.inputs_mut()[1] = 3.0;
    adder.inputs_mut()[2] = 2.0;

    // Update
    adder.update(0.0);

    // Check output: 10*1 - 3*1 + 2*2 = 10 - 3 + 4 = 11
    assert_relative_eq!(*adder.outputs().get(0).unwrap(), 11.0, epsilon = 1e-10);
}

#[test]
fn test_rk4_exponential_decay() {
    // Test dx/dt = -k*x, x(0) = 1
    // Analytical solution: x(t) = exp(-k*t)
    let k = 0.5;
    let initial = DVector::from_vec(vec![1.0]);
    let mut solver = RK4::new(initial);

    let dt = 0.01;
    let t_final = 5.0;
    let steps = (t_final / dt) as usize;

    for _ in 0..steps {
        solver.buffer(dt);

        // dx/dt = -k*x
        let f = |x: &DVector<f64>, _t: f64| DVector::from_vec(vec![-k * x[0]]);

        // RK4 requires 4 stage evaluations
        for _ in 0..4 {
            solver.step(f, dt);
        }
    }

    // Check final value
    let expected = (-k * t_final).exp();
    let actual = solver.state()[0];
    assert_relative_eq!(actual, expected, epsilon = 1e-6);
}

#[test]
fn test_rk4_harmonic_oscillator() {
    // Test d²x/dt² = -ω²x (simple harmonic oscillator)
    // With x(0) = 1, dx/dt(0) = 0
    // Solution: x(t) = cos(ωt), dx/dt = -ω*sin(ωt)

    let omega = 2.0 * std::f64::consts::PI;
    let initial = DVector::from_vec(vec![1.0, 0.0]); // [x, dx/dt]
    let mut solver = RK4::new(initial);

    let dt = 0.001;
    let t_final = 1.0; // One period
    let steps = (t_final / dt) as usize;

    for _ in 0..steps {
        solver.buffer(dt);

        // [dx/dt, d²x/dt²] = [v, -ω²x]
        let f = |state: &DVector<f64>, _t: f64| {
            DVector::from_vec(vec![state[1], -omega * omega * state[0]])
        };

        for _ in 0..4 {
            solver.step(f, dt);
        }
    }

    // After one period, should return to initial state
    let final_state = solver.state();
    assert_relative_eq!(final_state[0], 1.0, epsilon = 1e-3);
    assert_relative_eq!(final_state[1], 0.0, epsilon = 1e-2);
}

// NOTE: The following tests use the v1 Simulation API which has been removed
// in favor of the new v2 manual struct approach. These tests are disabled
// as part of the architectural migration.

// #[test]
// fn test_simulation_basic() {
//     let mut sim = Simulation::new();
//     // Create blocks
//     let adder1 = sim.add_block(Adder::new(2));
//     let adder2 = sim.add_block(Adder::new(2));
//     // Check block count
//     assert_eq!(sim.block_count(), 2);
//     // Connect blocks
//     sim.connect(adder1, adder2);
//     assert_eq!(sim.connection_count(), 1);
// }

// #[test]
// fn test_simulation_run() {
//     let mut sim = Simulation::new();
//     let adder = sim.add_block(Adder::new(2));
//     // Set initial inputs
//     if let Some(block) = sim.get_block_mut(adder) {
//         block.inputs_mut().set(0, 1.0);
//         block.inputs_mut().set(1, 2.0);
//     }
//     // Run simulation
//     let stats = sim.run(0.1).expect("Simulation failed");
//     assert!(stats.steps > 0);
//     assert_eq!(stats.successful_steps, stats.steps);
//     assert_eq!(stats.failed_steps, 0);
// }

// #[test]
// fn test_simulation_reset() {
//     let mut sim = Simulation::new();
//     let adder = sim.add_block(Adder::new(2));
//     // Run simulation
//     let _stats = sim.run(1.0).expect("Simulation failed");
//     assert!(sim.time > 0.0);
//     // Reset
//     sim.reset();
//     assert_eq!(sim.time, 0.0);
// }

// #[test]
// fn test_simulation_statistics() {
//     let mut sim = Simulation::new();
//     let _adder = sim.add_block(Adder::new(2));
//     sim.dt = 0.1;
//     let duration = 1.0;
//     let stats = sim.run(duration).expect("Simulation failed");
//     // Should have run approximately 10 steps (allowing for off-by-one)
//     assert!(stats.steps >= 10 && stats.steps <= 11);
//     assert_eq!(stats.successful_steps, stats.steps);
//     assert_eq!(stats.failed_steps, 0);
//     assert!(stats.runtime_ms > 0.0);
// }
