//! Unit tests for AppState logic (no UI dependencies)

use rustsim_types::Position;

// Import the state module - we need to make it public or use a test helper
// For now, we'll test the types directly

#[test]
fn test_position_creation() {
    let pos = Position::new(100.0, 200.0);
    assert_eq!(pos.x, 100.0);
    assert_eq!(pos.y, 200.0);
}

#[test]
fn test_position_zero() {
    let pos = Position::zero();
    assert_eq!(pos.x, 0.0);
    assert_eq!(pos.y, 0.0);
}

#[test]
fn test_snap_to_grid() {
    use rustsim_types::snap_to_grid;

    assert_eq!(snap_to_grid(0.0), 0.0);
    assert_eq!(snap_to_grid(5.0), 10.0); // Rounds up
    assert_eq!(snap_to_grid(4.9), 0.0); // Rounds down
    assert_eq!(snap_to_grid(15.0), 20.0);
    assert_eq!(snap_to_grid(25.0), 30.0);
    assert_eq!(snap_to_grid(-5.0), -10.0);
}

#[test]
fn test_node_instance_creation() {
    use rustsim_types::NodeInstance;

    let node = NodeInstance::new(
        "node-1".to_string(),
        "Integrator".to_string(),
        Position::new(100.0, 100.0),
    );

    assert_eq!(node.id, "node-1");
    assert_eq!(node.block_type, "Integrator");
    assert_eq!(node.name, "Integrator"); // Default name is block type
    assert_eq!(node.position.x, 100.0);
    assert_eq!(node.inputs.len(), 0);
    assert_eq!(node.outputs.len(), 0);
}

#[test]
fn test_node_add_ports() {
    use rustsim_types::NodeInstance;

    let mut node = NodeInstance::new("node-1".to_string(), "Adder".to_string(), Position::zero());

    node.add_input("in 0");
    node.add_input("in 1");
    node.add_output("out");

    assert_eq!(node.inputs.len(), 2);
    assert_eq!(node.outputs.len(), 1);
    assert_eq!(node.inputs[0].name, "in 0");
    assert_eq!(node.inputs[1].name, "in 1");
    assert_eq!(node.outputs[0].name, "out");
}

#[test]
fn test_node_parameters() {
    use rustsim_types::NodeInstance;
    use serde_json::json;

    let mut node = NodeInstance::new(
        "node-1".to_string(),
        "Constant".to_string(),
        Position::zero(),
    );

    node.set_param("value", json!(42.0));

    let value = node.get_param("value");
    assert!(value.is_some());
    assert_eq!(value.unwrap().as_f64(), Some(42.0));

    // Non-existent param
    assert!(node.get_param("nonexistent").is_none());
}

#[test]
fn test_connection_creation() {
    use rustsim_types::Connection;

    let conn = Connection::new(
        "conn-1".to_string(),
        "node-1".to_string(),
        0,
        "node-2".to_string(),
        0,
    );

    assert_eq!(conn.id, "conn-1");
    assert_eq!(conn.source_node_id, "node-1");
    assert_eq!(conn.source_port_index, 0);
    assert_eq!(conn.target_node_id, "node-2");
    assert_eq!(conn.target_port_index, 0);
    assert!(conn.waypoints.is_empty());
}

#[test]
fn test_simulation_graph() {
    use rustsim_types::{Connection, NodeInstance, SimulationGraph};

    let mut graph = SimulationGraph::new();

    // Add nodes
    let node1 = NodeInstance::new("n1".to_string(), "Constant".to_string(), Position::zero());
    let node2 = NodeInstance::new(
        "n2".to_string(),
        "Scope".to_string(),
        Position::new(200.0, 0.0),
    );

    graph.add_node(node1);
    graph.add_node(node2);

    assert_eq!(graph.nodes.len(), 2);
    assert!(graph.get_node("n1").is_some());
    assert!(graph.get_node("n2").is_some());

    // Add connection
    let conn = Connection::new("c1".to_string(), "n1".to_string(), 0, "n2".to_string(), 0);
    graph.add_connection(conn);

    assert_eq!(graph.connections.len(), 1);

    // Query connections
    let from_n1 = graph.get_connections_from("n1");
    assert_eq!(from_n1.len(), 1);

    let to_n2 = graph.get_connections_to("n2");
    assert_eq!(to_n2.len(), 1);

    // Remove node (should also remove connections)
    graph.remove_node("n1");
    assert_eq!(graph.nodes.len(), 1);
    assert_eq!(graph.connections.len(), 0);
}

#[test]
fn test_graph_file_serialization() {
    use rustsim_types::{GraphFile, NodeInstance, Position};

    let mut file = GraphFile::new("Test Model");

    let node = NodeInstance::new(
        "n1".to_string(),
        "Constant".to_string(),
        Position::new(50.0, 50.0),
    );
    file.graph.add_node(node);

    // Serialize
    let json = file.to_json().expect("Serialization failed");
    assert!(json.contains("Test Model"));
    assert!(json.contains("Constant"));

    // Deserialize
    let loaded = GraphFile::from_json(&json).expect("Deserialization failed");
    assert_eq!(loaded.metadata.name, "Test Model");
    assert_eq!(loaded.graph.nodes.len(), 1);
}

#[test]
fn test_simulation_settings_defaults() {
    use rustsim_types::SimulationSettings;

    let settings = SimulationSettings::default();

    assert_eq!(settings.duration, 10.0);
    assert_eq!(settings.dt, 0.01);
    assert!(!settings.adaptive);
    assert_eq!(settings.ghost_traces, 3);
}

#[test]
fn test_block_category_shape_mapping() {
    use rustsim_types::{BlockCategory, NodeShape};

    assert_eq!(BlockCategory::Sources.default_shape(), NodeShape::Pill);
    assert_eq!(BlockCategory::Dynamic.default_shape(), NodeShape::Rect);
    assert_eq!(BlockCategory::Algebraic.default_shape(), NodeShape::Rect);
    assert_eq!(BlockCategory::Mixed.default_shape(), NodeShape::Mixed);
    assert_eq!(BlockCategory::Recording.default_shape(), NodeShape::Pill);
}

#[test]
fn test_port_instance_creation() {
    use rustsim_types::{PortDirection, PortInstance};

    let input = PortInstance::new_input("node-1", 0, "in");
    assert_eq!(input.id, "node-1-input-0");
    assert_eq!(input.direction, PortDirection::Input);
    assert_eq!(input.index, 0);

    let output = PortInstance::new_output("node-1", 0, "out");
    assert_eq!(output.id, "node-1-output-0");
    assert_eq!(output.direction, PortDirection::Output);
}
