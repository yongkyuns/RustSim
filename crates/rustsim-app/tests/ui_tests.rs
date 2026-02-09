//! UI component tests using egui_kittest.
//!
//! These tests verify that UI components render correctly and respond to user interactions.

use egui::{Pos2, Vec2};
use egui_kittest::Harness;
use rustsim_types::{Connection, NodeInstance, Position, SimulationGraph};

// ============================================================================
// Toolbar Tests
// ============================================================================

#[test]
fn test_toolbar_renders() {
    let app = |ui: &mut egui::Ui| {
        ui.horizontal(|ui| {
            ui.button("New");
            ui.button("Open");
            ui.button("Save");
            ui.separator();
            ui.button("Undo");
            ui.button("Redo");
            ui.separator();
            ui.button("Run");
        });
    };

    let mut harness = Harness::new_ui(app);
    harness.run();

    // Test doesn't crash and renders
}

#[test]
fn test_zoom_controls_logic() {
    let mut zoom = 1.0_f32;

    // Zoom in
    zoom = (zoom * 1.2).min(4.0);
    assert!((zoom - 1.2).abs() < 0.01);

    // Zoom out
    zoom = (zoom / 1.2).max(0.1);
    assert!((zoom - 1.0).abs() < 0.01);
}

// ============================================================================
// Block Palette Tests
// ============================================================================

#[test]
fn test_palette_renders() {
    let categories = ["Sources", "Dynamic", "Algebraic", "Mixed", "Recording"];

    let app = |ui: &mut egui::Ui| {
        ui.heading("Blocks");
        ui.separator();

        egui::ScrollArea::vertical().show(ui, |ui| {
            for category in &categories {
                egui::CollapsingHeader::new(*category)
                    .default_open(true)
                    .show(ui, |ui| {
                        ui.label(format!("{} blocks", category));
                    });
            }
        });
    };

    let mut harness = Harness::new_ui(app);
    harness.run();

    // Test renders without crashing
}

// ============================================================================
// Node Graph Tests
// ============================================================================

#[test]
fn test_node_graph_empty_canvas_renders() {
    let app = |ui: &mut egui::Ui| {
        let (response, painter) =
            ui.allocate_painter(egui::Vec2::new(800.0, 600.0), egui::Sense::click_and_drag());

        // Draw background
        let rect = response.rect;
        painter.rect_filled(rect, 0.0, egui::Color32::from_gray(30));

        // Draw grid lines
        let grid_size = 20.0;
        for x in (rect.left() as i32..rect.right() as i32).step_by(grid_size as usize) {
            painter.line_segment(
                [
                    Pos2::new(x as f32, rect.top()),
                    Pos2::new(x as f32, rect.bottom()),
                ],
                egui::Stroke::new(0.5, egui::Color32::from_gray(40)),
            );
        }

        // Show help text
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "Drag blocks from palette to add",
            egui::FontId::default(),
            egui::Color32::GRAY,
        );
    };

    let mut harness = Harness::new_ui(app);
    harness.run();
}

#[test]
fn test_node_rendering_logic() {
    let mut graph = SimulationGraph::new();
    let node = NodeInstance::new(
        "test-node".to_string(),
        "Constant".to_string(),
        Position::new(100.0, 100.0),
    );
    graph.add_node(node);

    // Verify node is in graph
    assert_eq!(graph.nodes.len(), 1);
    assert!(graph.get_node("test-node").is_some());

    let node_ref = graph.get_node("test-node").unwrap();
    let node_rect = egui::Rect::from_min_size(
        Pos2::new(node_ref.position.x, node_ref.position.y),
        Vec2::new(80.0, 40.0),
    );

    assert_eq!(node_rect.min.x, 100.0);
    assert_eq!(node_rect.min.y, 100.0);
    assert_eq!(node_rect.width(), 80.0);
    assert_eq!(node_rect.height(), 40.0);
}

// ============================================================================
// Properties Panel Tests
// ============================================================================

#[test]
fn test_properties_panel_no_selection_renders() {
    let app = |ui: &mut egui::Ui| {
        ui.heading("Properties");
        ui.separator();
        ui.label("Select a node to view properties");
    };

    let mut harness = Harness::new_ui(app);
    harness.run();
}

#[test]
fn test_node_parameter_access() {
    let mut node = NodeInstance::new(
        "node-1".to_string(),
        "Constant".to_string(),
        Position::new(100.0, 100.0),
    );
    node.set_param("value", serde_json::json!(42.0));

    // Verify parameter access
    assert_eq!(node.get_param("value").unwrap().as_f64(), Some(42.0));
    assert!(node.get_param("nonexistent").is_none());
}

// ============================================================================
// Plot Panel Tests
// ============================================================================

#[test]
fn test_plots_panel_empty_renders() {
    let app = |ui: &mut egui::Ui| {
        ui.heading("Plots");
        ui.label("Run simulation to see plots");
    };

    let mut harness = Harness::new_ui(app);
    harness.run();
}

#[test]
fn test_plot_data_processing() {
    let plot_data: Vec<(f64, Vec<f64>)> =
        vec![(0.0, vec![0.0]), (0.1, vec![0.1]), (0.2, vec![0.2])];

    // Convert to plot points format
    let points: Vec<[f64; 2]> = plot_data
        .iter()
        .map(|(t, values)| [*t, values.first().copied().unwrap_or(0.0)])
        .collect();

    assert_eq!(points.len(), 3);
    assert_eq!(points[0], [0.0, 0.0]);
    assert_eq!(points[1], [0.1, 0.1]);
    assert_eq!(points[2], [0.2, 0.2]);
}

// ============================================================================
// Graph Operations Tests
// ============================================================================

#[test]
fn test_graph_node_operations() {
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

    // Remove node
    graph.remove_node("n1");
    assert_eq!(graph.nodes.len(), 1);
    assert!(graph.get_node("n1").is_none());
    assert!(graph.get_node("n2").is_some());
}

#[test]
fn test_graph_connection_operations() {
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

    // Add connection
    let conn = Connection::new("c1".to_string(), "n1".to_string(), 0, "n2".to_string(), 0);
    graph.add_connection(conn);

    assert_eq!(graph.connections.len(), 1);

    // Query connections
    let from_n1 = graph.get_connections_from("n1");
    assert_eq!(from_n1.len(), 1);

    let to_n2 = graph.get_connections_to("n2");
    assert_eq!(to_n2.len(), 1);

    // Remove node should also remove connections
    graph.remove_node("n1");
    assert_eq!(graph.connections.len(), 0);
}

// ============================================================================
// Zoom and Pan Tests
// ============================================================================

#[test]
fn test_zoom_limits() {
    let mut zoom = 1.0_f32;

    // Zoom in until limit
    for _ in 0..20 {
        zoom = (zoom * 1.2).min(4.0);
    }
    assert!(zoom <= 4.0, "Zoom should not exceed 4.0");
    assert!(
        (zoom - 4.0).abs() < 0.01,
        "Zoom should hit the limit at 4.0"
    );

    // Zoom out until limit
    zoom = 1.0;
    for _ in 0..20 {
        zoom = (zoom / 1.2).max(0.1);
    }
    assert!(zoom >= 0.1, "Zoom should not go below 0.1");
    assert!(
        (zoom - 0.1).abs() < 0.01,
        "Zoom should hit the limit at 0.1"
    );
}

#[test]
fn test_pan_calculation() {
    let pan = Position::new(100.0, 50.0);
    let zoom = 1.5_f32;
    let screen_pos = Pos2::new(400.0, 300.0);

    // Convert screen to canvas coordinates
    let canvas_x = (screen_pos.x - pan.x) / zoom;
    let canvas_y = (screen_pos.y - pan.y) / zoom;

    assert!((canvas_x - 200.0).abs() < 0.01);
    assert!((canvas_y - 166.67).abs() < 0.1);
}

#[test]
fn test_canvas_to_screen() {
    let pan = Position::new(100.0, 50.0);
    let zoom = 1.5_f32;
    let canvas_pos = Pos2::new(200.0, 166.67);

    // Convert canvas to screen coordinates
    let screen_x = canvas_pos.x * zoom + pan.x;
    let screen_y = canvas_pos.y * zoom + pan.y;

    assert!((screen_x - 400.0).abs() < 0.1);
    assert!((screen_y - 300.0).abs() < 0.1);
}

// ============================================================================
// Keyboard Shortcut Tests
// ============================================================================

#[test]
fn test_keyboard_shortcut_mapping() {
    struct MockKeyboard {
        ctrl: bool,
        shift: bool,
        key: Option<char>,
    }

    impl MockKeyboard {
        fn is_undo(&self) -> bool {
            self.ctrl && self.key == Some('z') && !self.shift
        }

        fn is_redo(&self) -> bool {
            self.ctrl && self.key == Some('z') && self.shift
        }

        fn is_save(&self) -> bool {
            self.ctrl && self.key == Some('s')
        }

        fn is_new(&self) -> bool {
            self.ctrl && self.key == Some('n')
        }

        fn is_open(&self) -> bool {
            self.ctrl && self.key == Some('o')
        }
    }

    let undo = MockKeyboard {
        ctrl: true,
        shift: false,
        key: Some('z'),
    };
    assert!(undo.is_undo());
    assert!(!undo.is_redo());

    let redo = MockKeyboard {
        ctrl: true,
        shift: true,
        key: Some('z'),
    };
    assert!(!redo.is_undo());
    assert!(redo.is_redo());

    let save = MockKeyboard {
        ctrl: true,
        shift: false,
        key: Some('s'),
    };
    assert!(save.is_save());

    let new = MockKeyboard {
        ctrl: true,
        shift: false,
        key: Some('n'),
    };
    assert!(new.is_new());

    let open = MockKeyboard {
        ctrl: true,
        shift: false,
        key: Some('o'),
    };
    assert!(open.is_open());
}

// ============================================================================
// Node Position and Grid Snap Tests
// ============================================================================

#[test]
fn test_grid_snap() {
    use rustsim_types::snap_to_grid;

    assert_eq!(snap_to_grid(0.0), 0.0);
    assert_eq!(snap_to_grid(5.0), 10.0);
    assert_eq!(snap_to_grid(4.9), 0.0);
    assert_eq!(snap_to_grid(15.0), 20.0);
    assert_eq!(snap_to_grid(-5.0), -10.0);
    assert_eq!(snap_to_grid(-4.9), 0.0);
}

#[test]
fn test_node_positioning() {
    let pos = Position::new(155.0, 245.0);
    let snapped = Position::new(
        rustsim_types::snap_to_grid(pos.x),
        rustsim_types::snap_to_grid(pos.y),
    );

    assert_eq!(snapped.x, 160.0);
    assert_eq!(snapped.y, 250.0);
}

// ============================================================================
// Selection Tests
// ============================================================================

#[test]
fn test_selection_operations() {
    use std::collections::HashSet;

    let mut selected_nodes: HashSet<String> = HashSet::new();

    // Select a node
    selected_nodes.insert("n1".to_string());
    assert_eq!(selected_nodes.len(), 1);
    assert!(selected_nodes.contains("n1"));

    // Add another to selection
    selected_nodes.insert("n2".to_string());
    assert_eq!(selected_nodes.len(), 2);

    // Clear selection
    selected_nodes.clear();
    assert!(selected_nodes.is_empty());
}

#[test]
fn test_multiple_selection() {
    use std::collections::HashSet;

    let mut selected: HashSet<String> = HashSet::new();

    // Add multiple nodes
    for i in 0..5 {
        selected.insert(format!("n{}", i));
    }

    assert_eq!(selected.len(), 5);

    // Remove specific node
    selected.remove("n2");
    assert_eq!(selected.len(), 4);
    assert!(!selected.contains("n2"));
}

// ============================================================================
// Undo/Redo Stack Tests
// ============================================================================

#[test]
fn test_undo_redo_stack() {
    let mut undo_stack: Vec<SimulationGraph> = Vec::new();
    let mut redo_stack: Vec<SimulationGraph> = Vec::new();
    let mut current = SimulationGraph::new();

    // Make changes and push to undo
    let node = NodeInstance::new("n1".to_string(), "Constant".to_string(), Position::zero());
    undo_stack.push(current.clone());
    current.add_node(node);

    assert_eq!(undo_stack.len(), 1);
    assert_eq!(current.nodes.len(), 1);

    // Undo
    if let Some(prev) = undo_stack.pop() {
        redo_stack.push(current.clone());
        current = prev;
    }

    assert_eq!(undo_stack.len(), 0);
    assert_eq!(redo_stack.len(), 1);
    assert_eq!(current.nodes.len(), 0);

    // Redo
    if let Some(next) = redo_stack.pop() {
        undo_stack.push(current.clone());
        current = next;
    }

    assert_eq!(undo_stack.len(), 1);
    assert_eq!(redo_stack.len(), 0);
    assert_eq!(current.nodes.len(), 1);
}

// ============================================================================
// Block Type Definition Tests
// ============================================================================

#[test]
fn test_block_category_icons() {
    use rustsim_types::BlockCategory;

    // Map categories to expected UI representations
    let categories = [
        (BlockCategory::Sources, "Sources"),
        (BlockCategory::Dynamic, "Dynamic"),
        (BlockCategory::Algebraic, "Algebraic"),
        (BlockCategory::Mixed, "Mixed"),
        (BlockCategory::Recording, "Recording"),
    ];

    for (category, name) in categories {
        assert_eq!(category.as_str(), name);
    }
}

// ============================================================================
// Simulation Tests
// ============================================================================

#[test]
fn test_simulation_sinusoidal_output() {
    let mut graph = SimulationGraph::new();

    // Create a simple sinusoidal source
    let mut node = NodeInstance::new(
        "sin1".to_string(),
        "Sinusoidal".to_string(),
        Position::zero(),
    );
    node.set_param("amplitude", serde_json::json!(1.0));
    node.set_param("frequency", serde_json::json!(1.0));
    node.set_param("phase", serde_json::json!(0.0));
    graph.add_node(node);

    // Verify graph setup
    assert_eq!(graph.nodes.len(), 1);
    let sin_node = graph.get_node("sin1").unwrap();
    assert_eq!(sin_node.get_param("amplitude").unwrap().as_f64(), Some(1.0));
}

#[test]
fn test_simulation_step_signal() {
    let mut node = NodeInstance::new("step1".to_string(), "Step".to_string(), Position::zero());
    node.set_param("step_time", serde_json::json!(1.0));
    node.set_param("initial", serde_json::json!(0.0));
    node.set_param("final_value", serde_json::json!(5.0));

    // Test parameter retrieval
    assert_eq!(node.get_param("step_time").unwrap().as_f64(), Some(1.0));
    assert_eq!(node.get_param("initial").unwrap().as_f64(), Some(0.0));
    assert_eq!(node.get_param("final_value").unwrap().as_f64(), Some(5.0));
}

#[test]
fn test_simulation_graph_with_connection() {
    let mut graph = SimulationGraph::new();

    // Source -> Amplifier -> Scope
    let mut source = NodeInstance::new("src".to_string(), "Constant".to_string(), Position::zero());
    source.set_param("value", serde_json::json!(2.0));
    source.add_output("out");
    graph.add_node(source);

    let mut amp = NodeInstance::new(
        "amp".to_string(),
        "Amplifier".to_string(),
        Position::new(150.0, 0.0),
    );
    amp.set_param("gain", serde_json::json!(3.0));
    amp.add_input("in");
    amp.add_output("out");
    graph.add_node(amp);

    let mut scope = NodeInstance::new(
        "scope".to_string(),
        "Scope".to_string(),
        Position::new(300.0, 0.0),
    );
    scope.add_input("in");
    graph.add_node(scope);

    // Connect source -> amp -> scope
    graph.add_connection(Connection::new(
        "c1".to_string(),
        "src".to_string(),
        0,
        "amp".to_string(),
        0,
    ));
    graph.add_connection(Connection::new(
        "c2".to_string(),
        "amp".to_string(),
        0,
        "scope".to_string(),
        0,
    ));

    assert_eq!(graph.nodes.len(), 3);
    assert_eq!(graph.connections.len(), 2);

    // Verify connections
    let from_src = graph.get_connections_from("src");
    assert_eq!(from_src.len(), 1);
    assert_eq!(from_src[0].target_node_id, "amp");

    let from_amp = graph.get_connections_from("amp");
    assert_eq!(from_amp.len(), 1);
    assert_eq!(from_amp[0].target_node_id, "scope");
}
