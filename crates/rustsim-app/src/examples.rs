//! Example simulation graphs for RustSim.
//!
//! This module provides pre-built example simulations that demonstrate
//! various control systems and dynamic modeling concepts.

use rustsim_types::{Connection, NodeInstance, Position, SimulationGraph};
use serde_json::json;

/// Helper to create a connection with auto-generated ID
fn conn(
    id: usize,
    source_node_id: &str,
    source_port_index: usize,
    target_node_id: &str,
    target_port_index: usize,
) -> Connection {
    Connection::new(
        format!("conn-{}", id),
        source_node_id.to_string(),
        source_port_index,
        target_node_id.to_string(),
        target_port_index,
    )
}

/// Create a PID controller example.
///
/// This example demonstrates a PID controller regulating a first-order plant.
/// The setpoint is a step function, and the controller drives the plant output
/// to track the setpoint.
///
/// System structure:
/// ```text
/// Step -> Adder(+,-) -> PID -> Integrator -> Gain -> Scope
///             ^                                 |
///             |_________________________________|
/// ```
pub fn pid_controller() -> SimulationGraph {
    let mut graph = SimulationGraph::new();

    // Step source (setpoint)
    let mut step = NodeInstance::new(
        "step".to_string(),
        "Step".to_string(),
        Position::new(50.0, 150.0),
    );
    step.add_output("out");
    step.set_param("step_time", json!(1.0));
    step.set_param("final_value", json!(1.0));

    // Adder for error signal
    let mut adder = NodeInstance::new(
        "adder".to_string(),
        "Adder".to_string(),
        Position::new(200.0, 150.0),
    );
    adder.add_input("+");
    adder.add_input("-");
    adder.add_output("sum");

    // PID controller
    let mut pid = NodeInstance::new(
        "pid".to_string(),
        "PID".to_string(),
        Position::new(350.0, 150.0),
    );
    pid.add_input("error");
    pid.add_output("out");
    pid.set_param("kp", json!(1.5));
    pid.set_param("ki", json!(0.5));
    pid.set_param("kd", json!(0.1));

    // Plant integrator
    let mut integrator = NodeInstance::new(
        "integrator".to_string(),
        "Integrator".to_string(),
        Position::new(500.0, 150.0),
    );
    integrator.add_input("in");
    integrator.add_output("out");
    integrator.set_param("initial_value", json!(0.0));

    // Plant gain
    let mut gain = NodeInstance::new(
        "gain".to_string(),
        "Amplifier".to_string(),
        Position::new(650.0, 150.0),
    );
    gain.add_input("in");
    gain.add_output("out");
    gain.set_param("gain", json!(0.4));

    // Scope for recording
    let mut scope = NodeInstance::new(
        "scope".to_string(),
        "Scope".to_string(),
        Position::new(800.0, 150.0),
    );
    scope.add_input("setpoint");
    scope.add_input("output");
    scope.add_input("error");

    // Add nodes
    graph.add_node(step);
    graph.add_node(adder);
    graph.add_node(pid);
    graph.add_node(integrator);
    graph.add_node(gain);
    graph.add_node(scope);

    // Connections
    graph.add_connection(conn(1, "step", 0, "adder", 0)); // Step -> Adder[0]
    graph.add_connection(conn(2, "step", 0, "scope", 0)); // Step -> Scope[0]
    graph.add_connection(conn(3, "adder", 0, "pid", 0)); // Adder -> PID
    graph.add_connection(conn(4, "adder", 0, "scope", 2)); // Adder -> Scope[2]
    graph.add_connection(conn(5, "pid", 0, "integrator", 0)); // PID -> Integrator
    graph.add_connection(conn(6, "integrator", 0, "gain", 0)); // Integrator -> Gain
    graph.add_connection(conn(7, "gain", 0, "adder", 1)); // Gain -> Adder[1] (feedback)
    graph.add_connection(conn(8, "gain", 0, "scope", 1)); // Gain -> Scope[1]

    graph
}

/// Create a harmonic oscillator (spring-mass-damper) example.
///
/// This example simulates a damped spring-mass system with the equation:
/// m*x'' + c*x' + k*x = 0
///
/// Rearranged: x'' = (-c*x' - k*x) / m
pub fn harmonic_oscillator() -> SimulationGraph {
    let mut graph = SimulationGraph::new();

    // Initial conditions
    let x0 = 2.0;
    let v0 = 5.0;

    // Parameters: mass=0.8, damping=0.2, spring=1.5
    let m = 0.8;
    let c = 0.2;
    let k = 1.5;

    // Integrator for velocity
    let mut int_vel = NodeInstance::new(
        "int_vel".to_string(),
        "Integrator".to_string(),
        Position::new(500.0, 100.0),
    );
    int_vel.add_input("accel");
    int_vel.add_output("velocity");
    int_vel.set_param("initial_value", json!(v0));

    // Integrator for position
    let mut int_pos = NodeInstance::new(
        "int_pos".to_string(),
        "Integrator".to_string(),
        Position::new(650.0, 100.0),
    );
    int_pos.add_input("velocity");
    int_pos.add_output("position");
    int_pos.set_param("initial_value", json!(x0));

    // Damping gain (c)
    let mut amp_c = NodeInstance::new(
        "amp_c".to_string(),
        "Amplifier".to_string(),
        Position::new(200.0, 50.0),
    );
    amp_c.add_input("in");
    amp_c.add_output("out");
    amp_c.set_param("gain", json!(c));

    // Spring gain (k)
    let mut amp_k = NodeInstance::new(
        "amp_k".to_string(),
        "Amplifier".to_string(),
        Position::new(200.0, 200.0),
    );
    amp_k.add_input("in");
    amp_k.add_output("out");
    amp_k.set_param("gain", json!(k));

    // Adder for c*v + k*x
    let mut adder = NodeInstance::new(
        "adder".to_string(),
        "Adder".to_string(),
        Position::new(300.0, 100.0),
    );
    adder.add_input("in0");
    adder.add_input("in1");
    adder.add_output("sum");

    // Final gain (-1/m)
    let mut amp_m = NodeInstance::new(
        "amp_m".to_string(),
        "Amplifier".to_string(),
        Position::new(400.0, 100.0),
    );
    amp_m.add_input("in");
    amp_m.add_output("out");
    amp_m.set_param("gain", json!(-1.0 / m));

    // Scope to record velocity and position
    let mut scope = NodeInstance::new(
        "scope".to_string(),
        "Scope".to_string(),
        Position::new(800.0, 100.0),
    );
    scope.add_input("velocity");
    scope.add_input("position");

    // Add nodes
    graph.add_node(int_vel);
    graph.add_node(int_pos);
    graph.add_node(amp_c);
    graph.add_node(amp_k);
    graph.add_node(adder);
    graph.add_node(amp_m);
    graph.add_node(scope);

    // Connections
    graph.add_connection(conn(1, "int_vel", 0, "int_pos", 0)); // velocity -> position integrator
    graph.add_connection(conn(2, "int_vel", 0, "amp_c", 0)); // velocity -> damping
    graph.add_connection(conn(3, "int_vel", 0, "scope", 0)); // velocity -> scope
    graph.add_connection(conn(4, "int_pos", 0, "amp_k", 0)); // position -> spring
    graph.add_connection(conn(5, "int_pos", 0, "scope", 1)); // position -> scope
    graph.add_connection(conn(6, "amp_c", 0, "adder", 0)); // damping -> adder
    graph.add_connection(conn(7, "amp_k", 0, "adder", 1)); // spring -> adder
    graph.add_connection(conn(8, "adder", 0, "amp_m", 0)); // sum -> mass gain
    graph.add_connection(conn(9, "amp_m", 0, "int_vel", 0)); // acceleration -> velocity (feedback)

    graph
}

/// Create a simple feedback system example.
///
/// This demonstrates a first-order system with negative feedback:
/// x' = input - x
pub fn simple_feedback() -> SimulationGraph {
    let mut graph = SimulationGraph::new();

    // Step source
    let mut step = NodeInstance::new(
        "step".to_string(),
        "Step".to_string(),
        Position::new(50.0, 100.0),
    );
    step.add_output("out");
    step.set_param("step_time", json!(3.0));
    step.set_param("final_value", json!(1.0));

    // Adder
    let mut adder = NodeInstance::new(
        "adder".to_string(),
        "Adder".to_string(),
        Position::new(200.0, 100.0),
    );
    adder.add_input("+");
    adder.add_input("-");
    adder.add_output("sum");

    // Integrator
    let mut integrator = NodeInstance::new(
        "integrator".to_string(),
        "Integrator".to_string(),
        Position::new(350.0, 100.0),
    );
    integrator.add_input("in");
    integrator.add_output("out");
    integrator.set_param("initial_value", json!(2.0));

    // Negative feedback gain
    let mut amp = NodeInstance::new(
        "amp".to_string(),
        "Amplifier".to_string(),
        Position::new(500.0, 100.0),
    );
    amp.add_input("in");
    amp.add_output("out");
    amp.set_param("gain", json!(-1.0));

    // Scope
    let mut scope = NodeInstance::new(
        "scope".to_string(),
        "Scope".to_string(),
        Position::new(650.0, 100.0),
    );
    scope.add_input("step");
    scope.add_input("response");

    graph.add_node(step);
    graph.add_node(adder);
    graph.add_node(integrator);
    graph.add_node(amp);
    graph.add_node(scope);

    graph.add_connection(conn(1, "step", 0, "adder", 0)); // Step -> Adder[0]
    graph.add_connection(conn(2, "step", 0, "scope", 0)); // Step -> Scope[0]
    graph.add_connection(conn(3, "amp", 0, "adder", 1)); // Amp -> Adder[1] (feedback)
    graph.add_connection(conn(4, "adder", 0, "integrator", 0)); // Adder -> Integrator
    graph.add_connection(conn(5, "integrator", 0, "amp", 0)); // Integrator -> Amp
    graph.add_connection(conn(6, "integrator", 0, "scope", 1)); // Integrator -> Scope[1]

    graph
}

/// Create a sinusoidal source demo.
pub fn sinusoidal_demo() -> SimulationGraph {
    let mut graph = SimulationGraph::new();

    let mut sin = NodeInstance::new(
        "sin".to_string(),
        "Sinusoidal".to_string(),
        Position::new(100.0, 100.0),
    );
    sin.add_output("out");
    sin.set_param("amplitude", json!(1.0));
    sin.set_param("frequency", json!(1.0));
    sin.set_param("phase", json!(0.0));

    let mut scope = NodeInstance::new(
        "scope".to_string(),
        "Scope".to_string(),
        Position::new(300.0, 100.0),
    );
    scope.add_input("signal");

    graph.add_node(sin);
    graph.add_node(scope);

    graph.add_connection(conn(1, "sin", 0, "scope", 0));

    graph
}

/// Get all available example names and descriptions.
pub fn list_examples() -> Vec<(&'static str, &'static str)> {
    vec![
        ("PID Controller", "PID controller with first-order plant"),
        ("Harmonic Oscillator", "Spring-mass-damper system"),
        ("Simple Feedback", "First-order negative feedback"),
        ("Sinusoidal Demo", "Basic sinusoidal source"),
    ]
}

/// Load an example by name.
pub fn load_example(name: &str) -> Option<SimulationGraph> {
    match name {
        "PID Controller" => Some(pid_controller()),
        "Harmonic Oscillator" => Some(harmonic_oscillator()),
        "Simple Feedback" => Some(simple_feedback()),
        "Sinusoidal Demo" => Some(sinusoidal_demo()),
        _ => None,
    }
}
