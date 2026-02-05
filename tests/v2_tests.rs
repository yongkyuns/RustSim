//! Comprehensive tests for v2 architecture
//!
//! Tests all implemented blocks and simulation patterns

use rustsim::*;

// ============================================================================
// Block Tests
// ============================================================================

#[test]
fn test_amplifier() {
    let mut amp = Amplifier::new(3.0);
    amp.set_input(0, 4.0);
    amp.update(0.0);
    assert_eq!(amp.get_output(0), 12.0);
}

#[test]
fn test_adder() {
    let mut adder = Adder::<3>::new();
    adder.inputs_mut()[0] = 1.0;
    adder.inputs_mut()[1] = 2.0;
    adder.inputs_mut()[2] = 3.0;
    adder.update(0.0);
    assert_eq!(adder.get_output(0), 6.0);
}

#[test]
fn test_adder_weighted() {
    let mut adder = Adder::<2>::with_weights([2.0, -1.0]);
    adder.inputs_mut()[0] = 5.0;
    adder.inputs_mut()[1] = 3.0;
    adder.update(0.0);
    assert_eq!(adder.get_output(0), 7.0); // 2*5 - 1*3
}

#[test]
fn test_subtractor() {
    let mut sub = Adder::<2>::subtractor();
    sub.inputs_mut()[0] = 10.0;
    sub.inputs_mut()[1] = 3.0;
    sub.update(0.0);
    assert_eq!(sub.get_output(0), 7.0);
}

#[test]
fn test_constant() {
    let mut constant = Constant::new(42.0);
    constant.update(0.0);
    assert_eq!(constant.get_output(0), 42.0);

    constant.update(100.0);
    assert_eq!(constant.get_output(0), 42.0);
}

#[test]
fn test_sinusoidal() {
    let mut sine = Sinusoidal::new(2.0, 1.0, 0.0);

    sine.update(0.0);
    assert!((sine.get_output(0) - 0.0).abs() < 1e-10);

    sine.update(0.25);
    assert!((sine.get_output(0) - 2.0).abs() < 1e-10);

    sine.update(0.5);
    assert!((sine.get_output(0) - 0.0).abs() < 1e-10);
}

#[test]
fn test_step() {
    let mut step = Step::new(5.0, 1.0);

    step.update(0.5);
    assert_eq!(step.get_output(0), 0.0);

    step.update(1.0);
    assert_eq!(step.get_output(0), 5.0);

    step.update(2.0);
    assert_eq!(step.get_output(0), 5.0);
}

#[test]
fn test_integrator_constant_input() {
    let mut integrator = Integrator::new(0.0);
    let dt = 0.01;

    // Integrate constant 1.0 for 1 second
    for _ in 0..100 {
        integrator.set_input(0, 1.0);
        integrator.update(0.0);
        integrator.step(0.0, dt);
    }

    assert!((integrator.value() - 1.0).abs() < 1e-10);
}

#[test]
fn test_integrator_initial_value() {
    let integrator = Integrator::new(5.0);
    assert_eq!(integrator.value(), 5.0);
    assert_eq!(integrator.get_output(0), 5.0);
}

#[test]
fn test_function_block() {
    let mut square = Function::<1, 1, _>::new(|inputs, outputs| {
        outputs[0] = inputs[0] * inputs[0];
    });

    square.set_input(0, 3.0);
    square.update(0.0);
    assert_eq!(square.get_output(0), 9.0);
}

#[test]
fn test_function_multi_input() {
    let mut mult = Function::<2, 1, _>::new(|inputs, outputs| {
        outputs[0] = inputs[0] * inputs[1];
    });

    mult.inputs_mut()[0] = 3.0;
    mult.inputs_mut()[1] = 4.0;
    mult.update(0.0);
    assert_eq!(mult.get_output(0), 12.0);
}

// ============================================================================
// Math Block Tests
// ============================================================================

#[test]
fn test_math_sin() {
    let mut sin = Sin::new();
    sin.set_input(0, std::f64::consts::PI / 2.0);
    sin.update(0.0);
    assert!((sin.get_output(0) - 1.0).abs() < 1e-10);
}

#[test]
fn test_math_cos() {
    let mut cos = Cos::new();
    cos.set_input(0, std::f64::consts::PI);
    cos.update(0.0);
    assert!((cos.get_output(0) + 1.0).abs() < 1e-10);
}

#[test]
fn test_math_exp() {
    let mut exp = Exp::new();
    exp.set_input(0, 1.0);
    exp.update(0.0);
    assert!((exp.get_output(0) - std::f64::consts::E).abs() < 1e-10);
}

#[test]
fn test_math_sqrt() {
    let mut sqrt = Sqrt::new();
    sqrt.set_input(0, 16.0);
    sqrt.update(0.0);
    assert_eq!(sqrt.get_output(0), 4.0);
}

#[test]
fn test_math_abs() {
    let mut abs = Abs::new();
    abs.set_input(0, -5.5);
    abs.update(0.0);
    assert_eq!(abs.get_output(0), 5.5);
}

#[test]
fn test_math_pow() {
    let mut pow = Pow::new(3.0);
    pow.set_input(0, 2.0);
    pow.update(0.0);
    assert_eq!(pow.get_output(0), 8.0);
}

#[test]
fn test_math_clip() {
    let mut clip = Clip::new(-1.0, 1.0);

    clip.set_input(0, 0.5);
    clip.update(0.0);
    assert_eq!(clip.get_output(0), 0.5);

    clip.set_input(0, 2.0);
    clip.update(0.0);
    assert_eq!(clip.get_output(0), 1.0);

    clip.set_input(0, -3.0);
    clip.update(0.0);
    assert_eq!(clip.get_output(0), -1.0);
}

// ============================================================================
// Reset Tests
// ============================================================================

#[test]
fn test_amplifier_reset() {
    let mut amp = Amplifier::new(2.0);
    amp.set_input(0, 5.0);
    amp.update(0.0);
    amp.reset();
    assert_eq!(amp.get_input(0), 0.0);
    assert_eq!(amp.get_output(0), 0.0);
}

#[test]
fn test_integrator_reset() {
    let mut integrator = Integrator::new(3.0);
    integrator.set_input(0, 1.0);
    integrator.update(0.0);
    integrator.step(0.0, 1.0);

    integrator.reset();
    assert_eq!(integrator.value(), 3.0);
    assert_eq!(integrator.get_input(0), 0.0);
}

// ============================================================================
// Buffer/Revert Tests
// ============================================================================

#[test]
fn test_integrator_buffer_revert() {
    let mut integrator = Integrator::new(0.0);

    // Take a step
    integrator.set_input(0, 1.0);
    integrator.update(0.0);
    integrator.step(0.0, 0.1);
    let value1 = integrator.value();

    // Buffer state
    integrator.buffer();

    // Take another step
    integrator.set_input(0, 10.0);
    integrator.update(0.0);
    integrator.step(0.0, 0.1);
    let value2 = integrator.value();

    // Values should be different
    assert!(value2 != value1);

    // Revert to buffered state
    integrator.revert();
    assert_eq!(integrator.value(), value1);
}

// ============================================================================
// Simulation Pattern Tests
// ============================================================================

#[test]
fn test_simple_chain() {
    // Test: constant -> amplifier
    struct SimpleChain {
        constant: Constant,
        amplifier: Amplifier,
    }

    impl SimpleChain {
        fn new() -> Self {
            Self {
                constant: Constant::new(5.0),
                amplifier: Amplifier::new(2.0),
            }
        }

        fn update(&mut self) {
            self.constant.update(0.0);
            self.amplifier.set_input(0, self.constant.get_output(0));
            self.amplifier.update(0.0);
        }

        fn output(&self) -> f64 {
            self.amplifier.get_output(0)
        }
    }

    let mut sim = SimpleChain::new();
    sim.update();
    assert_eq!(sim.output(), 10.0);
}

#[test]
fn test_harmonic_oscillator() {
    // Simple harmonic oscillator: d²x/dt² = -x
    struct HarmonicOscillator {
        gain: Amplifier,
        velocity: Integrator,
        position: Integrator,
        time: f64,
    }

    impl HarmonicOscillator {
        fn new(x0: f64, v0: f64) -> Self {
            Self {
                gain: Amplifier::new(-1.0),
                velocity: Integrator::new(v0),
                position: Integrator::new(x0),
                time: 0.0,
            }
        }

        fn propagate(&mut self) {
            self.gain.set_input(0, self.position.get_output(0));
            self.velocity.set_input(0, self.gain.get_output(0));
            self.position.set_input(0, self.velocity.get_output(0));
        }

        fn update(&mut self) {
            let t = self.time;
            self.position.update(t);
            self.propagate();
            self.gain.update(t);
            self.propagate();
            self.velocity.update(t);
            self.propagate();
        }

        fn step(&mut self, dt: f64) {
            self.update();
            self.velocity.step(self.time, dt);
            self.position.step(self.time, dt);
            self.time += dt;
        }

        fn position(&self) -> f64 {
            self.position.get_output(0)
        }

        fn velocity(&self) -> f64 {
            self.velocity.get_output(0)
        }
    }

    let mut sim = HarmonicOscillator::new(1.0, 0.0);
    let dt = 0.001;

    // Simulate for one period
    let period = 2.0 * std::f64::consts::PI;
    let steps = (period / dt) as usize;

    for _ in 0..steps {
        sim.step(dt);
    }

    // Should return to initial state (relaxed tolerance for RK4 integration)
    assert!((sim.position() - 1.0).abs() < 1e-2);
    assert!((sim.velocity() - 0.0).abs() < 1e-2);

    // Check energy conservation
    let energy = 0.5 * (sim.position().powi(2) + sim.velocity().powi(2));
    assert!((energy - 0.5).abs() < 1e-2);
}

#[test]
fn test_exponential_decay() {
    // Test: dx/dt = -k*x, x(0) = 1
    // Exact solution: x(t) = exp(-k*t)
    struct ExponentialDecay {
        gain: Amplifier,
        integrator: Integrator,
        time: f64,
    }

    impl ExponentialDecay {
        fn new(k: f64, x0: f64) -> Self {
            Self {
                gain: Amplifier::new(-k),
                integrator: Integrator::new(x0),
                time: 0.0,
            }
        }

        fn propagate(&mut self) {
            self.gain.set_input(0, self.integrator.get_output(0));
            self.integrator.set_input(0, self.gain.get_output(0));
        }

        fn update(&mut self) {
            let t = self.time;
            self.integrator.update(t);
            self.propagate();
            self.gain.update(t);
            self.propagate();
        }

        fn step(&mut self, dt: f64) {
            self.update();
            self.integrator.step(self.time, dt);
            self.time += dt;
        }

        fn value(&self) -> f64 {
            self.integrator.get_output(0)
        }
    }

    let k = 0.5;
    let mut sim = ExponentialDecay::new(k, 1.0);
    let dt = 0.01;
    let t_final = 5.0;
    let steps = (t_final / dt) as usize;

    for _ in 0..steps {
        sim.step(dt);
    }

    let expected = (-k * t_final).exp();
    let error = (sim.value() - expected).abs();
    assert!(error < 2e-3, "Error {} exceeds tolerance", error);
}

// ============================================================================
// Block Properties Tests
// ============================================================================

#[test]
fn test_block_num_ports() {
    let amp = Amplifier::new(1.0);
    assert_eq!(Amplifier::NUM_INPUTS, 1);
    assert_eq!(Amplifier::NUM_OUTPUTS, 1);
    assert_eq!(amp.inputs().len(), 1);
    assert_eq!(amp.outputs().len(), 1);

    let adder = Adder::<3>::new();
    assert_eq!(Adder::<3>::NUM_INPUTS, 3);
    assert_eq!(Adder::<3>::NUM_OUTPUTS, 1);
    assert_eq!(adder.inputs().len(), 3);
    assert_eq!(adder.outputs().len(), 1);

    let constant = Constant::new(1.0);
    assert_eq!(Constant::NUM_INPUTS, 0);
    assert_eq!(Constant::NUM_OUTPUTS, 1);
    assert_eq!(constant.inputs().len(), 0);
    assert_eq!(constant.outputs().len(), 1);
}

#[test]
fn test_block_is_dynamic() {
    assert!(!Amplifier::IS_DYNAMIC);
    assert!(!Constant::IS_DYNAMIC);
    assert!(Integrator::IS_DYNAMIC);
}

// ============================================================================
// Complex Simulation Tests
// ============================================================================

#[test]
fn test_feedback_loop() {
    // Test a simple feedback: x[n+1] = 0.9 * x[n]
    struct FeedbackLoop {
        gain: Amplifier,
        integrator: Integrator,
        time: f64,
    }

    impl FeedbackLoop {
        fn new(x0: f64) -> Self {
            Self {
                gain: Amplifier::new(0.9),
                integrator: Integrator::new(x0),
                time: 0.0,
            }
        }

        fn step(&mut self, dt: f64) {
            // Feedback: integrator output -> gain -> integrator input
            self.integrator.update(self.time);
            self.gain.set_input(0, self.integrator.get_output(0));
            self.gain.update(self.time);
            self.integrator.set_input(0, self.gain.get_output(0) - self.integrator.get_output(0));
            self.integrator.step(self.time, dt);
            self.time += dt;
        }

        fn value(&self) -> f64 {
            self.integrator.get_output(0)
        }
    }

    let mut sim = FeedbackLoop::new(10.0);

    // Should decay toward zero
    for _ in 0..100 {
        sim.step(0.01);
    }

    assert!(sim.value() < 10.0);
    assert!(sim.value() > 0.0);
}

#[test]
fn test_multi_rate_simulation() {
    // Test that we can have different dt values
    let mut integrator = Integrator::new(0.0);

    // Step with dt = 0.1
    integrator.set_input(0, 1.0);
    integrator.update(0.0);
    integrator.step(0.0, 0.1);

    let value1 = integrator.value();

    // Step with dt = 0.01
    integrator.update(0.0);
    integrator.step(0.0, 0.01);

    let value2 = integrator.value();

    assert!(value2 > value1);
    assert!((value1 - 0.1).abs() < 1e-10);
    assert!((value2 - 0.11).abs() < 1e-10);
}
