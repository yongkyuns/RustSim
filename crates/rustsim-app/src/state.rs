//! Application state management.

use std::collections::{HashMap, HashSet};

use rustsim_types::{
    BlockCategory, BlockTypeDefinition, Connection, NodeInstance, ParamDefinition, PortConfig,
    PortDefinition, Position, SimulationGraph, SimulationSettings,
};

/// Application state
pub struct AppState {
    /// The simulation graph
    pub graph: SimulationGraph,

    /// Simulation settings
    pub settings: SimulationSettings,

    /// Block type definitions
    pub block_types: Vec<BlockTypeDefinition>,

    /// Selected node IDs
    pub selected_nodes: HashSet<String>,

    /// Selected connection IDs
    pub selected_connections: HashSet<String>,

    /// Canvas pan offset
    pub pan: Position,

    /// Canvas zoom level
    pub zoom: f32,

    /// Simulation running state
    running: bool,

    /// Simulation time
    pub sim_time: f64,

    /// Recorded plot data: (time, outputs)
    pub plot_data: Vec<(f64, Vec<f64>)>,

    /// Undo stack
    undo_stack: Vec<UndoState>,

    /// Redo stack
    redo_stack: Vec<UndoState>,

    /// ID counter for generating unique IDs
    id_counter: u64,

    /// Currently dragging node
    pub dragging_node: Option<String>,

    /// Connection being created
    pub pending_connection: Option<PendingConnection>,
}

/// State snapshot for undo/redo
#[derive(Clone)]
struct UndoState {
    graph: SimulationGraph,
    settings: SimulationSettings,
}

/// Pending connection being drawn
pub struct PendingConnection {
    pub source_node: String,
    pub source_port: usize,
    pub current_pos: Position,
}

impl AppState {
    pub fn new() -> Self {
        Self {
            graph: SimulationGraph::new(),
            settings: SimulationSettings::default(),
            block_types: Self::create_block_types(),
            selected_nodes: HashSet::new(),
            selected_connections: HashSet::new(),
            pan: Position::zero(),
            zoom: 1.0,
            running: false,
            sim_time: 0.0,
            plot_data: Vec::new(),
            undo_stack: Vec::new(),
            redo_stack: Vec::new(),
            id_counter: 0,
            dragging_node: None,
            pending_connection: None,
        }
    }

    /// Generate a unique ID
    pub fn generate_id(&mut self) -> String {
        self.id_counter += 1;
        format!("node-{}", self.id_counter)
    }

    /// Create block type definitions for all RustSim blocks
    fn create_block_types() -> Vec<BlockTypeDefinition> {
        vec![
            // Sources
            BlockTypeDefinition::new("Constant", "Constant output value", BlockCategory::Sources)
                .with_param("value", ParamDefinition::number("Output value", 1.0))
                .with_ports(PortConfig {
                    inputs: vec![],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 0,
                    min_outputs: 1,
                    max_inputs: Some(0),
                    max_outputs: Some(1),
                    sync_ports: false,
                }),
            BlockTypeDefinition::new(
                "Sinusoidal",
                "Sinusoidal signal generator",
                BlockCategory::Sources,
            )
            .with_param("amplitude", ParamDefinition::number("Amplitude", 1.0))
            .with_param("frequency", ParamDefinition::number("Frequency (Hz)", 1.0))
            .with_param("phase", ParamDefinition::number("Phase (rad)", 0.0))
            .with_ports(PortConfig {
                inputs: vec![],
                outputs: vec![PortDefinition::output("out")],
                min_inputs: 0,
                min_outputs: 1,
                max_inputs: Some(0),
                max_outputs: Some(1),
                sync_ports: false,
            }),
            BlockTypeDefinition::new("Step", "Step signal", BlockCategory::Sources)
                .with_param("step_time", ParamDefinition::number("Step time", 1.0))
                .with_param("initial", ParamDefinition::number("Initial value", 0.0))
                .with_param("final_value", ParamDefinition::number("Final value", 1.0)),
            BlockTypeDefinition::new("Ramp", "Ramp signal", BlockCategory::Sources)
                .with_param("slope", ParamDefinition::number("Slope", 1.0))
                .with_param("start_time", ParamDefinition::number("Start time", 0.0)),
            // Dynamic blocks
            BlockTypeDefinition::new("Integrator", "Integrates input signal", BlockCategory::Dynamic)
                .with_param("initial_value", ParamDefinition::number("Initial value", 0.0))
                .with_ports(PortConfig {
                    inputs: vec![PortDefinition::input("in")],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 1,
                    min_outputs: 1,
                    max_inputs: None,
                    max_outputs: None,
                    sync_ports: true,
                }),
            BlockTypeDefinition::new(
                "Differentiator",
                "Differentiates input signal",
                BlockCategory::Dynamic,
            )
            .with_ports(PortConfig {
                inputs: vec![PortDefinition::input("in")],
                outputs: vec![PortDefinition::output("out")],
                min_inputs: 1,
                min_outputs: 1,
                max_inputs: None,
                max_outputs: None,
                sync_ports: true,
            }),
            BlockTypeDefinition::new("Delay", "Time delay", BlockCategory::Dynamic)
                .with_param("delay", ParamDefinition::number("Delay (s)", 0.1)),
            BlockTypeDefinition::new("PID", "PID controller", BlockCategory::Dynamic)
                .with_param("kp", ParamDefinition::number("Proportional gain", 1.0))
                .with_param("ki", ParamDefinition::number("Integral gain", 0.0))
                .with_param("kd", ParamDefinition::number("Derivative gain", 0.0)),
            // Algebraic blocks
            BlockTypeDefinition::new("Adder", "Adds input signals", BlockCategory::Algebraic)
                .with_ports(PortConfig {
                    inputs: vec![
                        PortDefinition::input("in 0"),
                        PortDefinition::input("in 1"),
                    ],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 2,
                    min_outputs: 1,
                    max_inputs: None,
                    max_outputs: Some(1),
                    sync_ports: false,
                }),
            BlockTypeDefinition::new(
                "Multiplier",
                "Multiplies input signals",
                BlockCategory::Algebraic,
            )
            .with_ports(PortConfig {
                inputs: vec![
                    PortDefinition::input("in 0"),
                    PortDefinition::input("in 1"),
                ],
                outputs: vec![PortDefinition::output("out")],
                min_inputs: 2,
                min_outputs: 1,
                max_inputs: None,
                max_outputs: Some(1),
                sync_ports: false,
            }),
            BlockTypeDefinition::new("Amplifier", "Amplifies signal by gain", BlockCategory::Algebraic)
                .with_param("gain", ParamDefinition::number("Gain", 1.0))
                .with_ports(PortConfig {
                    inputs: vec![PortDefinition::input("in")],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 1,
                    min_outputs: 1,
                    max_inputs: None,
                    max_outputs: None,
                    sync_ports: true,
                }),
            BlockTypeDefinition::new("Sin", "Sine function", BlockCategory::Algebraic),
            BlockTypeDefinition::new("Cos", "Cosine function", BlockCategory::Algebraic),
            BlockTypeDefinition::new("Abs", "Absolute value", BlockCategory::Algebraic),
            BlockTypeDefinition::new("Sqrt", "Square root", BlockCategory::Algebraic),
            // Mixed blocks
            BlockTypeDefinition::new("SampleHold", "Sample and hold", BlockCategory::Mixed),
            BlockTypeDefinition::new("Counter", "Event counter", BlockCategory::Mixed),
            // Recording blocks
            BlockTypeDefinition::new("Scope", "Records time-domain signals", BlockCategory::Recording)
                .with_ports(PortConfig {
                    inputs: vec![PortDefinition::input("in 0")],
                    outputs: vec![],
                    min_inputs: 1,
                    min_outputs: 0,
                    max_inputs: None,
                    max_outputs: Some(0),
                    sync_ports: false,
                }),
            BlockTypeDefinition::new("Spectrum", "Frequency analysis", BlockCategory::Recording)
                .with_ports(PortConfig {
                    inputs: vec![PortDefinition::input("in")],
                    outputs: vec![],
                    min_inputs: 1,
                    min_outputs: 0,
                    max_inputs: None,
                    max_outputs: Some(0),
                    sync_ports: false,
                }),
        ]
    }

    /// Get block type by name
    pub fn get_block_type(&self, name: &str) -> Option<&BlockTypeDefinition> {
        self.block_types.iter().find(|b| b.name == name)
    }

    /// Add a node to the graph
    pub fn add_node(&mut self, block_type: &str, position: Position) -> String {
        self.save_undo_state();

        let id = self.generate_id();
        let mut node = NodeInstance::new(id.clone(), block_type.to_string(), position);

        // Set up default ports based on block type
        if let Some(block_def) = self.get_block_type(block_type) {
            for port_def in &block_def.ports.inputs {
                node.add_input(&port_def.name);
            }
            for port_def in &block_def.ports.outputs {
                node.add_output(&port_def.name);
            }

            // Set default parameters
            for (name, param_def) in &block_def.params {
                if let Some(default) = &param_def.default {
                    if let Ok(value) = serde_json::from_str(default) {
                        node.set_param(name, value);
                    }
                }
            }
        }

        self.graph.add_node(node);
        id
    }

    /// Add a connection between nodes
    pub fn add_connection(
        &mut self,
        source_node: &str,
        source_port: usize,
        target_node: &str,
        target_port: usize,
    ) -> String {
        self.save_undo_state();

        let id = format!("conn-{}", self.id_counter);
        self.id_counter += 1;

        let connection = Connection::new(
            id.clone(),
            source_node.to_string(),
            source_port,
            target_node.to_string(),
            target_port,
        );

        self.graph.add_connection(connection);
        id
    }

    /// Delete selected items
    pub fn delete_selected(&mut self) {
        if self.selected_nodes.is_empty() && self.selected_connections.is_empty() {
            return;
        }

        self.save_undo_state();

        for node_id in self.selected_nodes.drain() {
            self.graph.remove_node(&node_id);
        }

        for conn_id in self.selected_connections.drain() {
            self.graph.remove_connection(&conn_id);
        }
    }

    /// Rotate selected nodes by 90 degrees
    pub fn rotate_selected(&mut self) {
        if self.selected_nodes.is_empty() {
            return;
        }

        self.save_undo_state();

        for node_id in &self.selected_nodes {
            if let Some(node) = self.graph.get_node_mut(node_id) {
                node.rotation = (node.rotation + 1) % 4;
            }
        }
    }

    /// Clear selection
    pub fn clear_selection(&mut self) {
        self.selected_nodes.clear();
        self.selected_connections.clear();
    }

    /// Zoom in
    pub fn zoom_in(&mut self) {
        self.zoom = (self.zoom * 1.2).min(4.0);
    }

    /// Zoom out
    pub fn zoom_out(&mut self) {
        self.zoom = (self.zoom / 1.2).max(0.1);
    }

    /// Fit all nodes in view
    pub fn fit_view(&mut self) {
        if self.graph.nodes.is_empty() {
            self.pan = Position::zero();
            self.zoom = 1.0;
            return;
        }

        // Calculate bounding box
        let mut min_x = f32::MAX;
        let mut min_y = f32::MAX;
        let mut max_x = f32::MIN;
        let mut max_y = f32::MIN;

        for node in self.graph.nodes.values() {
            min_x = min_x.min(node.position.x);
            min_y = min_y.min(node.position.y);
            max_x = max_x.max(node.position.x + 80.0);
            max_y = max_y.max(node.position.y + 40.0);
        }

        // Center the view
        self.pan = Position::new(-(min_x + max_x) / 2.0, -(min_y + max_y) / 2.0);
        self.zoom = 1.0;
    }

    // File operations

    pub fn new_file(&mut self) {
        self.graph = SimulationGraph::new();
        self.settings = SimulationSettings::default();
        self.selected_nodes.clear();
        self.selected_connections.clear();
        self.plot_data.clear();
        self.undo_stack.clear();
        self.redo_stack.clear();
        self.id_counter = 0;
    }

    pub fn open_file(&mut self) {
        // TODO: Implement file dialog
    }

    pub fn save_file(&mut self) {
        // TODO: Implement file save
    }

    // Undo/Redo

    fn save_undo_state(&mut self) {
        self.undo_stack.push(UndoState {
            graph: self.graph.clone(),
            settings: self.settings.clone(),
        });
        self.redo_stack.clear();

        // Limit undo stack size
        if self.undo_stack.len() > 50 {
            self.undo_stack.remove(0);
        }
    }

    pub fn undo(&mut self) {
        if let Some(state) = self.undo_stack.pop() {
            self.redo_stack.push(UndoState {
                graph: self.graph.clone(),
                settings: self.settings.clone(),
            });
            self.graph = state.graph;
            self.settings = state.settings;
        }
    }

    pub fn redo(&mut self) {
        if let Some(state) = self.redo_stack.pop() {
            self.undo_stack.push(UndoState {
                graph: self.graph.clone(),
                settings: self.settings.clone(),
            });
            self.graph = state.graph;
            self.settings = state.settings;
        }
    }

    // Simulation control

    pub fn is_running(&self) -> bool {
        self.running
    }

    pub fn run_simulation(&mut self) {
        self.running = true;
        self.sim_time = 0.0;
        self.plot_data.clear();
    }

    pub fn stop_simulation(&mut self) {
        self.running = false;
    }

    pub fn step_simulation(&mut self) {
        // TODO: Step the compiled simulation
        self.sim_time += self.settings.dt;

        // Dummy data for now
        self.plot_data
            .push((self.sim_time, vec![self.sim_time.sin()]));
    }
}

impl Default for AppState {
    fn default() -> Self {
        Self::new()
    }
}
