//! Integration tests for RustSim

use approx::assert_relative_eq;
use rustsim::prelude::*;
use rustsim::blocks::{Amplifier, Adder, Constant, Integrator, Scope};
use rustsim::solvers::{ExplicitSolver, RK4, Solver};
use nalgebra::DVector;

#[test]
fn test_amplifier_basic() {
    let mut amp = Amplifier::new(2.5);

    // Set input
    amp.inputs_mut()[0] = 4.0;

    // Update
    amp.update(0.0);

    // Check output
    assert_relative_eq!(amp*.outputs().get(0).unwrap(), 10.0, epsilon = 1e-10);
}

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
    assert_relative_eq!(adder*.outputs().get(0).unwrap(), 6.0, epsilon = 1e-10);
}

#[test]
fn test_subtractor() {
    let mut sub = Adder::with_weights([1.0, -1.0]);

    // Set inputs
    sub.inputs_mut()[0] = 5.0;
    sub.inputs_mut()[1] = 3.0;

    // Update
    sub.update(0.0);

    // Check output
    assert_relative_eq!(sub*.outputs().get(0).unwrap(), 2.0, epsilon = 1e-10);
}

#[test]
fn test_constant_source() {
    let mut constant = Constant::new(42.0);

    // Update
    constant.update(0.0);

    // Check output
    assert_relative_eq!(constant*.outputs().get(0).unwrap(), 42.0, epsilon = 1e-10);

    // Should stay constant
    constant.update(100.0);
    assert_relative_eq!(constant*.outputs().get(0).unwrap(), 42.0, epsilon = 1e-10);
}

// Note: FunctionSource is not in the public API yet
// #[test]
// fn test_function_source() {
//     let mut source = FunctionSource::new(|t| t * 2.0);
//
//     // Update at different times
//     source.update(0.0);
//     assert_relative_eq!(source*.outputs().get(0).unwrap(), 0.0, epsilon = 1e-10);
//
//     source.update(5.0);
//     assert_relative_eq!(source*.outputs().get(0).unwrap(), 10.0, epsilon = 1e-10);
//
//     source.update(10.0);
//     assert_relative_eq!(source*.outputs().get(0).unwrap(), 20.0, epsilon = 1e-10);
// }

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

    for i in 0..steps {
        let _t = i as f64 * dt;
        solver.buffer(dt);

        // dx/dt = -k*x
        let f = |x: &DVector<f64>, _t: f64| {
            DVector::from_vec(vec![-k * x[0]])
        };

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
fn test_rk4_sine_integration() {
    // Test d²x/dt² = -x (simple harmonic oscillator)
    // With x(0) = 0, dx/dt(0) = 1
    // Solution: x(t) = sin(t), dx/dt = cos(t)

    let initial = DVector::from_vec(vec![0.0, 1.0]); // [x, dx/dt]
    let mut solver = RK4::new(initial);

    let dt = 0.01;
    let t_final = 2.0 * std::f64::consts::PI; // One period
    let steps = (t_final / dt) as usize;

    for i in 0..steps {
        let _t = i as f64 * dt;
        solver.buffer(dt);

        // [dx/dt, d²x/dt²] = [v, -x]
        let f = |state: &DVector<f64>, _t: f64| {
            DVector::from_vec(vec![state[1], -state[0]])
        };

        for _ in 0..4 {
            solver.step(f, dt);
        }
    }

    // After one period, should return to initial state
    let final_state = solver.state();
    assert_relative_eq!(final_state[0], 0.0, epsilon = 1e-4);
    assert_relative_eq!(final_state[1], 1.0, epsilon = 1e-4);
}

// NOTE: The following tests use the v1 Simulation API which has been removed
// in favor of the new v2 manual struct approach. These tests are disabled
// as part of the architectural migration.

// #[test]
// fn test_simulation_basic() {
//     let mut sim = Simulation::new();
//     // Create blocks
//     let amp = sim.add_block(Amplifier::new(2.0));
//     let _constant = sim.add_block(Constant::new(3.0));
//     // Check block count
//     assert_eq!(sim.block_count(), 2);
//     // Connect blocks
//     sim.connect(amp, amp);
//     assert_eq!(sim.connection_count(), 1);
// }

// #[test]
// fn test_simulation_loop() {
//     let mut sim = Simulation::new();
//     // Create a simple amplifier chain
//     let amp1 = sim.add_block(Amplifier::new(2.0));
//     let amp2 = sim.add_block(Amplifier::new(3.0));
//     // Set initial input
//     if let Some(block) = sim.get_block_mut(amp1) {
//         block.inputs_mut()[0] = 1.0;
//     }
//     // Connect amp1 -> amp2
//     sim.connect(amp1, amp2);
//     // Run simulation
//     let stats = sim.run(0.1).expect("Simulation failed");
//     assert!(stats.steps > 0);
//     assert_eq!(stats.successful_steps, stats.steps);
//     assert_eq!(stats.failed_steps, 0);
// }

// #[test]
// fn test_integrator_constant_input() {
//     // Integrate a constant: ∫1 dt = t
//     let mut sim = Simulation::new();
//     let constant = sim.add_block(Constant::new(1.0));
//     let mut integrator = Integrator::new(0.0);
//     // Initialize the integrator with an RK4 solver
//     integrator.set_solver(RK4::new(DVector::from_vec(vec![0.0])));
//     let integrator_id = sim.add_block(integrator);
//     sim.connect(constant, integrator_id);
//     // Run for 1 second with dt=0.01
//     sim.dt = 0.01;
//     let duration = 1.0;
//     let _stats = sim.run(duration).expect("Simulation failed");
//     // Check integrator output
//     if let Some(block) = sim.get_block(integrator_id) {
//         let output = block*.outputs().get(0).unwrap();
//         assert_relative_eq!(output, duration, epsilon = 1e-2);
//     }
// }

// Note: ZeroCrossing closure needs to capture by move
// #[test]
// fn test_zero_crossing_event() {
//     let event = ZeroCrossing::new(
//         |t| t - 0.5, // Cross zero at t=0.5
//         None,  // No action function for now
//         1e-6,
//     );
//
//     let mut sim = Simulation::new();
//     sim.add_event(event);
//
//     sim.dt = 0.01;
//     let _stats = sim.run(1.0).expect("Simulation failed");
//
//     // Event should have been detected (though not resolved in current implementation)
//     // This test mainly verifies the API works
// }

// NOTE: Scope API has changed in v2 - this test needs to be rewritten for new API
// #[test]
// fn test_scope_capture() {
//     let mut scope = Scope::new();
//     // Simulate some values
//     for i in 0..10 {
//         let t = i as f64 * 0.1;
//         scope.inputs_mut()[0] = t * 2.0;
//         scope.update(t);
//         scope.sample(t, 0.1);
//     }
//     // Check history
//     assert_eq!(scope.history().len(), 10);
//     for (i, &(t, val)) in scope.history().iter().enumerate() {
//         let expected_t = i as f64 * 0.1;
//         let expected_val = expected_t * 2.0;
//         assert_relative_eq!(t, expected_t, epsilon = 1e-10);
//         assert_relative_eq!(val, expected_val, epsilon = 1e-10);
//     }
// }

// #[test]
// fn test_simulation_reset() {
//     let mut sim = Simulation::new();
//     let amp = sim.add_block(Amplifier::new(2.0));
//     // Set input and run
//     if let Some(block) = sim.get_block_mut(amp) {
//         block.inputs_mut()[0] = 5.0;
//     }
//     let _stats = sim.run(1.0).expect("Simulation failed");
//     assert!(sim.time > 0.0);
//     // Reset
//     sim.reset();
//     assert_eq!(sim.time, 0.0);
// }

// #[test]
// fn test_simulation_statistics() {
//     let mut sim = Simulation::new();
//     let _amp = sim.add_block(Amplifier::new(1.5));
//     sim.dt = 0.1;
//     let duration = 1.0;
//     let stats = sim.run(duration).expect("Simulation failed");
//     // Should have run 10 steps
//     assert_eq!(stats.steps, 10);
//     assert_eq!(stats.successful_steps, 10);
//     assert_eq!(stats.failed_steps, 0);
//     assert!(stats.runtime_ms > 0.0);
// }
