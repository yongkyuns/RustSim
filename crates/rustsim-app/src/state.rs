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

    /// Integrator states for simulation
    integrator_states: HashMap<String, f64>,

    /// Node being edited (double-clicked)
    pub editing_node: Option<String>,

    /// Whether to generate C ABI exports in code viewer
    pub generate_c_abi: bool,
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
            integrator_states: HashMap::new(),
            editing_node: None,
            generate_c_abi: false,
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
                .with_param("final_value", ParamDefinition::number("Final value", 1.0))
                .with_ports(PortConfig {
                    inputs: vec![],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 0,
                    min_outputs: 1,
                    max_inputs: Some(0),
                    max_outputs: Some(1),
                    sync_ports: false,
                }),
            BlockTypeDefinition::new("Ramp", "Ramp signal", BlockCategory::Sources)
                .with_param("slope", ParamDefinition::number("Slope", 1.0))
                .with_param("start_time", ParamDefinition::number("Start time", 0.0))
                .with_ports(PortConfig {
                    inputs: vec![],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 0,
                    min_outputs: 1,
                    max_inputs: Some(0),
                    max_outputs: Some(1),
                    sync_ports: false,
                }),
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
                .with_param("delay", ParamDefinition::number("Delay (s)", 0.1))
                .with_ports(PortConfig {
                    inputs: vec![PortDefinition::input("in")],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 1,
                    min_outputs: 1,
                    max_inputs: Some(1),
                    max_outputs: Some(1),
                    sync_ports: false,
                }),
            BlockTypeDefinition::new("PID", "PID controller", BlockCategory::Dynamic)
                .with_param("kp", ParamDefinition::number("Proportional gain", 1.0))
                .with_param("ki", ParamDefinition::number("Integral gain", 0.0))
                .with_param("kd", ParamDefinition::number("Derivative gain", 0.0))
                .with_ports(PortConfig {
                    inputs: vec![PortDefinition::input("error")],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 1,
                    min_outputs: 1,
                    max_inputs: Some(1),
                    max_outputs: Some(1),
                    sync_ports: false,
                }),
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
            BlockTypeDefinition::new("Sin", "Sine function", BlockCategory::Algebraic)
                .with_ports(PortConfig {
                    inputs: vec![PortDefinition::input("in")],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 1,
                    min_outputs: 1,
                    max_inputs: Some(1),
                    max_outputs: Some(1),
                    sync_ports: false,
                }),
            BlockTypeDefinition::new("Cos", "Cosine function", BlockCategory::Algebraic)
                .with_ports(PortConfig {
                    inputs: vec![PortDefinition::input("in")],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 1,
                    min_outputs: 1,
                    max_inputs: Some(1),
                    max_outputs: Some(1),
                    sync_ports: false,
                }),
            BlockTypeDefinition::new("Abs", "Absolute value", BlockCategory::Algebraic)
                .with_ports(PortConfig {
                    inputs: vec![PortDefinition::input("in")],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 1,
                    min_outputs: 1,
                    max_inputs: Some(1),
                    max_outputs: Some(1),
                    sync_ports: false,
                }),
            BlockTypeDefinition::new("Sqrt", "Square root", BlockCategory::Algebraic)
                .with_ports(PortConfig {
                    inputs: vec![PortDefinition::input("in")],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 1,
                    min_outputs: 1,
                    max_inputs: Some(1),
                    max_outputs: Some(1),
                    sync_ports: false,
                }),
            // Mixed blocks
            BlockTypeDefinition::new("SampleHold", "Sample and hold", BlockCategory::Mixed)
                .with_ports(PortConfig {
                    inputs: vec![PortDefinition::input("in"), PortDefinition::input("trigger")],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 2,
                    min_outputs: 1,
                    max_inputs: Some(2),
                    max_outputs: Some(1),
                    sync_ports: false,
                }),
            BlockTypeDefinition::new("Counter", "Event counter", BlockCategory::Mixed)
                .with_ports(PortConfig {
                    inputs: vec![PortDefinition::input("trigger")],
                    outputs: vec![PortDefinition::output("count")],
                    min_inputs: 1,
                    min_outputs: 1,
                    max_inputs: Some(1),
                    max_outputs: Some(1),
                    sync_ports: false,
                }),
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

    // Port management

    /// Check if an input port can be added to a node
    pub fn can_add_input(&self, node_id: &str) -> bool {
        let node = match self.graph.get_node(node_id) {
            Some(n) => n,
            None => return false,
        };
        let block_def = match self.get_block_type(&node.block_type) {
            Some(b) => b,
            None => return false,
        };
        match block_def.ports.max_inputs {
            Some(max) => node.inputs.len() < max,
            None => true, // Unlimited
        }
    }

    /// Check if an input port can be removed from a node
    pub fn can_remove_input(&self, node_id: &str) -> bool {
        let node = match self.graph.get_node(node_id) {
            Some(n) => n,
            None => return false,
        };
        let block_def = match self.get_block_type(&node.block_type) {
            Some(b) => b,
            None => return false,
        };
        node.inputs.len() > block_def.ports.min_inputs
    }

    /// Check if an output port can be added to a node
    pub fn can_add_output(&self, node_id: &str) -> bool {
        let node = match self.graph.get_node(node_id) {
            Some(n) => n,
            None => return false,
        };
        let block_def = match self.get_block_type(&node.block_type) {
            Some(b) => b,
            None => return false,
        };
        // For sync_ports blocks, outputs are controlled by inputs
        if block_def.ports.sync_ports {
            return false;
        }
        match block_def.ports.max_outputs {
            Some(max) => node.outputs.len() < max,
            None => true, // Unlimited
        }
    }

    /// Check if an output port can be removed from a node
    pub fn can_remove_output(&self, node_id: &str) -> bool {
        let node = match self.graph.get_node(node_id) {
            Some(n) => n,
            None => return false,
        };
        let block_def = match self.get_block_type(&node.block_type) {
            Some(b) => b,
            None => return false,
        };
        // For sync_ports blocks, outputs are controlled by inputs
        if block_def.ports.sync_ports {
            return false;
        }
        node.outputs.len() > block_def.ports.min_outputs
    }

    /// Check if a node has sync_ports enabled
    pub fn has_sync_ports(&self, node_id: &str) -> bool {
        let node = match self.graph.get_node(node_id) {
            Some(n) => n,
            None => return false,
        };
        self.get_block_type(&node.block_type)
            .map(|b| b.ports.sync_ports)
            .unwrap_or(false)
    }

    /// Add an input port to a node
    pub fn add_input_port(&mut self, node_id: &str) -> bool {
        if !self.can_add_input(node_id) {
            return false;
        }

        self.save_undo_state();

        let sync_ports = self.has_sync_ports(node_id);

        if let Some(node) = self.graph.get_node_mut(node_id) {
            let port_name = format!("in {}", node.inputs.len());
            node.add_input(&port_name);

            // For sync_ports blocks, also add an output
            if sync_ports {
                let out_name = format!("out {}", node.outputs.len());
                node.add_output(&out_name);
            }
            true
        } else {
            false
        }
    }

    /// Remove the last input port from a node
    pub fn remove_input_port(&mut self, node_id: &str) -> bool {
        if !self.can_remove_input(node_id) {
            return false;
        }

        self.save_undo_state();

        let sync_ports = self.has_sync_ports(node_id);

        // Get the port index to remove
        let port_index = if let Some(node) = self.graph.get_node(node_id) {
            node.inputs.len().saturating_sub(1)
        } else {
            return false;
        };

        // Remove connections to this port
        self.graph.connections.retain(|conn| {
            !(conn.target_node_id == node_id && conn.target_port_index == port_index)
        });

        // For sync_ports, also remove connections from the corresponding output
        if sync_ports {
            self.graph.connections.retain(|conn| {
                !(conn.source_node_id == node_id && conn.source_port_index == port_index)
            });
        }

        if let Some(node) = self.graph.get_node_mut(node_id) {
            node.remove_input();
            if sync_ports {
                node.remove_output();
            }
            true
        } else {
            false
        }
    }

    /// Add an output port to a node
    pub fn add_output_port(&mut self, node_id: &str) -> bool {
        if !self.can_add_output(node_id) {
            return false;
        }

        self.save_undo_state();

        if let Some(node) = self.graph.get_node_mut(node_id) {
            let port_name = format!("out {}", node.outputs.len());
            node.add_output(&port_name);
            true
        } else {
            false
        }
    }

    /// Remove the last output port from a node
    pub fn remove_output_port(&mut self, node_id: &str) -> bool {
        if !self.can_remove_output(node_id) {
            return false;
        }

        self.save_undo_state();

        // Get the port index to remove
        let port_index = if let Some(node) = self.graph.get_node(node_id) {
            node.outputs.len().saturating_sub(1)
        } else {
            return false;
        };

        // Remove connections from this port
        self.graph.connections.retain(|conn| {
            !(conn.source_node_id == node_id && conn.source_port_index == port_index)
        });

        if let Some(node) = self.graph.get_node_mut(node_id) {
            node.remove_output();
            true
        } else {
            false
        }
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
        self.integrator_states.clear();

        // Initialize integrator states from node parameters
        for (id, node) in &self.graph.nodes {
            if node.block_type == "Integrator" {
                let initial = node.get_param("initial_value")
                    .and_then(|v| v.as_f64())
                    .unwrap_or(0.0);
                self.integrator_states.insert(id.clone(), initial);
            }
        }
    }

    pub fn stop_simulation(&mut self) {
        self.running = false;
    }

    pub fn step_simulation(&mut self) {
        let dt = self.settings.dt;
        let mut outputs: Vec<f64> = Vec::new();

        // Get execution order (simple topological sort)
        let order = self.get_execution_order();

        // Store computed values for wiring
        let mut values: std::collections::HashMap<String, f64> = std::collections::HashMap::new();

        for node_id in &order {
            if let Some(node) = self.graph.get_node(node_id) {
                let output = match node.block_type.as_str() {
                    "Constant" => {
                        node.get_param("value")
                            .and_then(|v| v.as_f64())
                            .unwrap_or(1.0)
                    }
                    "Sinusoidal" => {
                        let amp = node.get_param("amplitude").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        let freq = node.get_param("frequency").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        let phase = node.get_param("phase").and_then(|v| v.as_f64()).unwrap_or(0.0);
                        amp * (2.0 * std::f64::consts::PI * freq * self.sim_time + phase).sin()
                    }
                    "Step" => {
                        let step_time = node.get_param("step_time").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        let initial = node.get_param("initial").and_then(|v| v.as_f64()).unwrap_or(0.0);
                        let final_val = node.get_param("final_value").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        if self.sim_time >= step_time { final_val } else { initial }
                    }
                    "Ramp" => {
                        let slope = node.get_param("slope").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        let start_time = node.get_param("start_time").and_then(|v| v.as_f64()).unwrap_or(0.0);
                        if self.sim_time >= start_time {
                            slope * (self.sim_time - start_time)
                        } else {
                            0.0
                        }
                    }
                    "Amplifier" | "Gain" => {
                        let gain = node.get_param("gain")
                            .or_else(|| node.get_param("k"))
                            .and_then(|v| v.as_f64())
                            .unwrap_or(1.0);
                        let input = self.get_input_value(node_id, 0, &values);
                        input * gain
                    }
                    "Adder" => {
                        let mut sum = 0.0;
                        for i in 0..node.inputs.len().max(2) {
                            sum += self.get_input_value(node_id, i, &values);
                        }
                        sum
                    }
                    "Multiplier" => {
                        let mut product = 1.0;
                        for i in 0..node.inputs.len().max(2) {
                            product *= self.get_input_value(node_id, i, &values);
                        }
                        product
                    }
                    "Sin" => {
                        let input = self.get_input_value(node_id, 0, &values);
                        input.sin()
                    }
                    "Cos" => {
                        let input = self.get_input_value(node_id, 0, &values);
                        input.cos()
                    }
                    "Abs" => {
                        let input = self.get_input_value(node_id, 0, &values);
                        input.abs()
                    }
                    "Sqrt" => {
                        let input = self.get_input_value(node_id, 0, &values);
                        input.abs().sqrt()
                    }
                    "Integrator" => {
                        let input = self.get_input_value(node_id, 0, &values);
                        // Get current state value and update it
                        let current = *self.integrator_states.get(node_id).unwrap_or(&0.0);
                        let new_val = current + input * dt;
                        self.integrator_states.insert(node_id.clone(), new_val);
                        new_val
                    }
                    "Scope" | "Spectrum" => {
                        // Recording blocks just pass through their input
                        let input = self.get_input_value(node_id, 0, &values);
                        outputs.push(input);
                        input
                    }
                    _ => 0.0,
                };

                values.insert(node_id.clone(), output);
            }
        }

        // If no scope blocks, record all source block outputs
        if outputs.is_empty() {
            for node_id in &order {
                if let Some(node) = self.graph.get_node(node_id) {
                    if matches!(node.block_type.as_str(), "Constant" | "Sinusoidal" | "Step" | "Ramp" | "Integrator") {
                        if let Some(&val) = values.get(node_id) {
                            outputs.push(val);
                        }
                    }
                }
            }
        }

        self.plot_data.push((self.sim_time, outputs));
        self.sim_time += dt;
    }

    /// Get execution order using topological sort (stable ordering)
    fn get_execution_order(&self) -> Vec<String> {
        let mut in_degree: HashMap<String, usize> = HashMap::new();
        let mut adjacency: HashMap<String, Vec<String>> = HashMap::new();

        // Initialize
        for id in self.graph.nodes.keys() {
            in_degree.insert(id.clone(), 0);
            adjacency.insert(id.clone(), Vec::new());
        }

        // Build graph
        for conn in &self.graph.connections {
            if let Some(adj) = adjacency.get_mut(&conn.source_node_id) {
                adj.push(conn.target_node_id.clone());
            }
            if let Some(deg) = in_degree.get_mut(&conn.target_node_id) {
                *deg += 1;
            }
        }

        // Kahn's algorithm with stable ordering (sort by ID for determinism)
        let mut ready: Vec<String> = in_degree
            .iter()
            .filter(|(_, &deg)| deg == 0)
            .map(|(id, _)| id.clone())
            .collect();
        ready.sort(); // Alphabetical order for stability

        let mut result = Vec::new();
        while let Some(node) = ready.pop() {
            result.push(node.clone());

            if let Some(neighbors) = adjacency.get(&node) {
                let mut new_ready = Vec::new();
                for neighbor in neighbors {
                    if let Some(deg) = in_degree.get_mut(neighbor) {
                        *deg -= 1;
                        if *deg == 0 {
                            new_ready.push(neighbor.clone());
                        }
                    }
                }
                // Sort new ready nodes and add them
                new_ready.sort();
                new_ready.reverse(); // Reverse so pop() gives smallest first
                ready.extend(new_ready);
            }
        }

        // Add any remaining nodes in sorted order (in case of cycles)
        let mut remaining: Vec<_> = self.graph.nodes.keys()
            .filter(|id| !result.contains(id))
            .cloned()
            .collect();
        remaining.sort();
        result.extend(remaining);

        result
    }

    /// Get input value from connected node
    fn get_input_value(&self, target_node_id: &str, target_port: usize, values: &std::collections::HashMap<String, f64>) -> f64 {
        for conn in &self.graph.connections {
            if conn.target_node_id == target_node_id && conn.target_port_index == target_port {
                if let Some(&val) = values.get(&conn.source_node_id) {
                    return val;
                }
            }
        }
        0.0
    }
}

impl Default for AppState {
    fn default() -> Self {
        Self::new()
    }
}
