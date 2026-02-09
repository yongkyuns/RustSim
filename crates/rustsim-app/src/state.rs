//! Application state management.

use std::collections::{HashMap, HashSet, VecDeque};
use std::time::{Duration, Instant};

use rustsim_core::block::Block;
use rustsim_core::block_kind::BlockKind;
use rustsim_core::blocks::*;
use rustsim_types::{
    BlockCategory, BlockTypeDefinition, Connection, NodeInstance, ParamDefinition, PortConfig,
    PortDefinition, Position, SimulationGraph, SimulationSettings,
};

use crate::examples;
use crate::layout;
use crate::ui::PlotSettings;

#[cfg(not(target_arch = "wasm32"))]
use crate::compiler::CompiledSimulation;

#[cfg(not(target_arch = "wasm32"))]
use std::sync::mpsc::{self, Receiver};

/// Simulation result snapshot for ghost traces
#[derive(Clone)]
pub struct SimulationResult {
    pub plot_data: Vec<(f64, Vec<f64>)>,
    pub plot_labels: Vec<String>,
}

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

    /// Whether to show grid on canvas
    pub show_canvas_grid: bool,

    /// Simulation running state
    running: bool,

    /// Simulation time
    pub sim_time: f64,

    /// Recorded plot data: (time, outputs)
    pub plot_data: Vec<(f64, Vec<f64>)>,

    /// Plot signal labels
    pub plot_labels: Vec<String>,

    /// Per-node plot data for Scope nodes: node_id -> Vec<(time, values)>
    pub scope_data: HashMap<String, Vec<(f64, Vec<f64>)>>,

    /// Per-node spectrum data for Spectrum nodes: node_id -> Vec<(freq, magnitude)>
    pub spectrum_data: HashMap<String, Vec<(f64, Vec<f64>)>>,

    /// Plot visualization settings
    pub plot_settings: PlotSettings,

    /// Result history for ghost traces (max 6)
    pub result_history: VecDeque<SimulationResult>,

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

    /// Block instances for interpreter mode (keyed by node_id)
    blocks: HashMap<String, BlockKind>,

    /// Execution order for interpreter (cached topological sort)
    execution_order: Vec<String>,

    /// Node being edited (double-clicked)
    pub editing_node: Option<String>,

    /// Whether to use compiled mode instead of interpreter
    pub use_compiled_mode: bool,

    /// Compiled simulation instance (only on native platforms)
    #[cfg(not(target_arch = "wasm32"))]
    compiled_sim: Option<CompiledSimulation>,

    /// Compilation status message
    pub compilation_status: CompilationStatus,

    /// Compilation log messages
    pub compilation_log: Vec<String>,

    /// Receiver for compilation messages (from background thread)
    #[cfg(not(target_arch = "wasm32"))]
    compilation_receiver: Option<Receiver<CompilationMessage>>,

    /// Total time spent in step_simulation calls
    total_step_time: Duration,

    /// Number of steps executed
    step_count: u64,
}

/// Status of compilation
#[derive(Clone, Debug, PartialEq)]
pub enum CompilationStatus {
    /// Not compiled yet
    NotCompiled,
    /// Currently compiling
    Compiling,
    /// Compilation succeeded
    Ready,
    /// Compilation failed with error message
    Error(String),
}

/// Message from compilation thread
#[cfg(not(target_arch = "wasm32"))]
pub enum CompilationMessage {
    /// Log message
    Log(String),
    /// Compilation completed successfully
    Success(CompiledSimulation),
    /// Compilation failed
    Error(String),
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
            // Load PID Controller as default example
            graph: examples::pid_controller(),
            settings: SimulationSettings::default(),
            block_types: Self::create_block_types(),
            selected_nodes: HashSet::new(),
            selected_connections: HashSet::new(),
            pan: Position::zero(),
            zoom: 1.0,
            show_canvas_grid: true,
            running: false,
            sim_time: 0.0,
            plot_data: Vec::new(),
            plot_labels: Vec::new(),
            scope_data: HashMap::new(),
            spectrum_data: HashMap::new(),
            plot_settings: PlotSettings::default(),
            result_history: VecDeque::new(),
            undo_stack: Vec::new(),
            redo_stack: Vec::new(),
            id_counter: 0,
            dragging_node: None,
            pending_connection: None,
            blocks: HashMap::new(),
            execution_order: Vec::new(),
            editing_node: None,
            use_compiled_mode: false,
            #[cfg(not(target_arch = "wasm32"))]
            compiled_sim: None,
            compilation_status: CompilationStatus::NotCompiled,
            compilation_log: Vec::new(),
            #[cfg(not(target_arch = "wasm32"))]
            compilation_receiver: None,
            total_step_time: Duration::ZERO,
            step_count: 0,
        }
    }

    /// Load an example simulation by name
    pub fn load_example(&mut self, name: &str) {
        if let Some(graph) = examples::load_example(name) {
            // Stop any running simulation
            self.running = false;

            // Reset state
            self.graph = graph;
            self.sim_time = 0.0;
            self.plot_data.clear();
            self.plot_labels.clear();
            self.scope_data.clear();
            self.spectrum_data.clear();
            self.blocks.clear();
            self.execution_order.clear();
            self.selected_nodes.clear();
            self.selected_connections.clear();
            self.result_history.clear();
            self.undo_stack.clear();
            self.redo_stack.clear();

            // Reset compilation state
            self.use_compiled_mode = false;
            #[cfg(not(target_arch = "wasm32"))]
            {
                self.compiled_sim = None;
            }
            self.compilation_status = CompilationStatus::NotCompiled;
            self.compilation_log.clear();
        }
    }

    /// Get list of available examples
    pub fn available_examples() -> Vec<(&'static str, &'static str)> {
        examples::list_examples()
    }

    /// Auto-arrange blocks using layered layout algorithm
    pub fn auto_layout(&mut self) {
        layout::auto_layout(&mut self.graph);
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
            BlockTypeDefinition::new(
                "Integrator",
                "Integrates input signal",
                BlockCategory::Dynamic,
            )
            .with_param(
                "initial_value",
                ParamDefinition::number("Initial value", 0.0),
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
                    inputs: vec![PortDefinition::input("in 0"), PortDefinition::input("in 1")],
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
                inputs: vec![PortDefinition::input("in 0"), PortDefinition::input("in 1")],
                outputs: vec![PortDefinition::output("out")],
                min_inputs: 2,
                min_outputs: 1,
                max_inputs: None,
                max_outputs: Some(1),
                sync_ports: false,
            }),
            BlockTypeDefinition::new(
                "Amplifier",
                "Amplifies signal by gain",
                BlockCategory::Algebraic,
            )
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
            BlockTypeDefinition::new("Sin", "Sine function", BlockCategory::Algebraic).with_ports(
                PortConfig {
                    inputs: vec![PortDefinition::input("in")],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 1,
                    min_outputs: 1,
                    max_inputs: Some(1),
                    max_outputs: Some(1),
                    sync_ports: false,
                },
            ),
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
            BlockTypeDefinition::new("Abs", "Absolute value", BlockCategory::Algebraic).with_ports(
                PortConfig {
                    inputs: vec![PortDefinition::input("in")],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 1,
                    min_outputs: 1,
                    max_inputs: Some(1),
                    max_outputs: Some(1),
                    sync_ports: false,
                },
            ),
            BlockTypeDefinition::new("Sqrt", "Square root", BlockCategory::Algebraic).with_ports(
                PortConfig {
                    inputs: vec![PortDefinition::input("in")],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 1,
                    min_outputs: 1,
                    max_inputs: Some(1),
                    max_outputs: Some(1),
                    sync_ports: false,
                },
            ),
            // Mixed blocks
            BlockTypeDefinition::new("SampleHold", "Sample and hold", BlockCategory::Mixed)
                .with_ports(PortConfig {
                    inputs: vec![
                        PortDefinition::input("in"),
                        PortDefinition::input("trigger"),
                    ],
                    outputs: vec![PortDefinition::output("out")],
                    min_inputs: 2,
                    min_outputs: 1,
                    max_inputs: Some(2),
                    max_outputs: Some(1),
                    sync_ports: false,
                }),
            BlockTypeDefinition::new("Counter", "Event counter", BlockCategory::Mixed).with_ports(
                PortConfig {
                    inputs: vec![PortDefinition::input("trigger")],
                    outputs: vec![PortDefinition::output("count")],
                    min_inputs: 1,
                    min_outputs: 1,
                    max_inputs: Some(1),
                    max_outputs: Some(1),
                    sync_ports: false,
                },
            ),
            // Recording blocks
            BlockTypeDefinition::new(
                "Scope",
                "Records time-domain signals",
                BlockCategory::Recording,
            )
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

    /// Invalidate compiled simulation when graph changes
    fn invalidate_compiled_sim(&mut self) {
        #[cfg(not(target_arch = "wasm32"))]
        {
            self.compiled_sim = None;
        }
        if self.compilation_status != CompilationStatus::NotCompiled {
            self.compilation_status = CompilationStatus::NotCompiled;
        }
    }

    /// Add a node to the graph
    pub fn add_node(&mut self, block_type: &str, position: Position) -> String {
        self.save_undo_state();
        self.invalidate_compiled_sim();

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
        self.invalidate_compiled_sim();

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
        self.invalidate_compiled_sim();

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
        self.plot_labels.clear();
        self.scope_data.clear();
        self.spectrum_data.clear();
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

    /// Start async compilation of the current simulation graph
    #[cfg(not(target_arch = "wasm32"))]
    pub fn compile_simulation(&mut self) -> Result<(), String> {
        use rustsim_codegen::CodeGenerator;

        // Clear previous compilation log
        self.compilation_log.clear();

        // Set status to compiling
        self.compilation_status = CompilationStatus::Compiling;
        self.compilation_log.push("Starting compilation...".to_string());

        // Generate code with parameter mapping
        self.compilation_log.push("Generating Rust code from simulation graph...".to_string());
        let generator = CodeGenerator::new()
            .with_c_abi(false)
            .with_library_mode(false);
        let (code, param_mapping) = generator.generate_with_params(&self.graph, &self.settings);

        self.compilation_log.push(format!(
            "Code generation complete. {} runtime parameters.",
            param_mapping.params.len()
        ));

        // Count states, inputs, and outputs
        let (num_inputs, num_outputs, num_states) = self.count_io_for_compilation();
        self.compilation_log.push(format!(
            "Simulation dimensions: {} inputs, {} outputs, {} states",
            num_inputs, num_outputs, num_states
        ));

        // Get workspace root
        self.compilation_log.push("Locating workspace root...".to_string());
        let workspace_root = std::env::current_dir()
            .ok()
            .and_then(|mut p| {
                // Try to find workspace root by looking for Cargo.toml with [workspace]
                loop {
                    let cargo_toml = p.join("Cargo.toml");
                    if cargo_toml.exists() {
                        if let Ok(contents) = std::fs::read_to_string(&cargo_toml) {
                            if contents.contains("[workspace]") {
                                return Some(p);
                            }
                        }
                    }
                    if !p.pop() {
                        break;
                    }
                }
                None
            })
            .ok_or_else(|| "Could not find workspace root".to_string())?;

        let workspace_root_str = workspace_root
            .to_str()
            .ok_or_else(|| "Invalid workspace path".to_string())?
            .to_string();

        self.compilation_log.push(format!("Workspace root: {}", workspace_root_str));

        // Create channel for communication with compilation thread
        let (tx, rx) = mpsc::channel();
        self.compilation_receiver = Some(rx);

        // Spawn compilation in background thread
        std::thread::spawn(move || {
            // Send log messages
            let _ = tx.send(CompilationMessage::Log("Creating temporary build directory...".to_string()));
            let _ = tx.send(CompilationMessage::Log("Writing generated code to src/lib.rs...".to_string()));
            let _ = tx.send(CompilationMessage::Log("Creating Cargo.toml manifest...".to_string()));
            let _ = tx.send(CompilationMessage::Log("Running 'cargo build --release --lib'...".to_string()));
            let _ = tx.send(CompilationMessage::Log("This may take a while on first compilation...".to_string()));

            // Perform compilation with parameter mapping
            match CompiledSimulation::with_params(&code, &workspace_root_str, num_states, num_inputs, num_outputs, param_mapping) {
                Ok(mut sim) => {
                    let _ = tx.send(CompilationMessage::Log("Compilation successful!".to_string()));
                    let _ = tx.send(CompilationMessage::Log("Loading compiled library...".to_string()));
                    sim.init();
                    let _ = tx.send(CompilationMessage::Log("Initializing simulation state...".to_string()));
                    let _ = tx.send(CompilationMessage::Success(sim));
                }
                Err(e) => {
                    let _ = tx.send(CompilationMessage::Log(format!("Compilation failed: {}", e)));
                    let _ = tx.send(CompilationMessage::Error(e.to_string()));
                }
            }
        });

        self.compilation_log.push("Background compilation thread started.".to_string());
        Ok(())
    }

    /// Poll compilation progress (call this regularly from UI update)
    #[cfg(not(target_arch = "wasm32"))]
    pub fn poll_compilation(&mut self) {
        let mut should_close = false;

        if let Some(ref rx) = self.compilation_receiver {
            // Process all available messages without blocking
            while let Ok(msg) = rx.try_recv() {
                match msg {
                    CompilationMessage::Log(text) => {
                        self.compilation_log.push(text);
                    }
                    CompilationMessage::Success(sim) => {
                        self.compiled_sim = Some(sim);
                        self.compilation_status = CompilationStatus::Ready;
                        self.compilation_log.push("Compilation complete and ready!".to_string());
                        should_close = true;
                    }
                    CompilationMessage::Error(err) => {
                        self.compilation_status = CompilationStatus::Error(err.clone());
                        self.compilation_log.push(format!("Error: {}", err));
                        should_close = true;
                    }
                }
            }
        }

        if should_close {
            self.compilation_receiver = None; // Close channel
        }
    }

    /// Count inputs, outputs, and states for compiled simulation
    #[cfg(not(target_arch = "wasm32"))]
    fn count_io_for_compilation(&self) -> (usize, usize, usize) {
        let num_inputs = 0;
        let mut num_outputs = 0;
        let mut num_states = 0;

        for node in self.graph.nodes.values() {
            match node.block_type.as_str() {
                "Scope" | "Spectrum" => {
                    // Recording blocks contribute to outputs
                    num_outputs += node.inputs.len().max(1);
                }
                "Integrator" => {
                    num_states += 1;
                    num_outputs += 1;
                }
                "Differentiator" => {
                    num_states += 1;
                    num_outputs += 1;
                }
                "Constant" | "Sinusoidal" | "Step" | "Ramp" => {
                    num_outputs += 1;
                }
                _ => {
                    num_outputs += node.outputs.len().max(1);
                }
            }
        }

        (num_inputs.max(1), num_outputs.max(1), num_states.max(1))
    }

    /// On WASM, compilation is not supported yet
    #[cfg(target_arch = "wasm32")]
    pub fn compile_simulation(&mut self) -> Result<(), String> {
        Err("Compilation not supported on WASM yet".to_string())
    }

    pub fn run_simulation(&mut self) {
        // If compiled mode is enabled and not yet compiled, compile first
        #[cfg(not(target_arch = "wasm32"))]
        if self.use_compiled_mode && self.compiled_sim.is_none() {
            if let Err(e) = self.compile_simulation() {
                eprintln!("Compilation failed: {}. Falling back to interpreter.", e);
                self.use_compiled_mode = false;
            }
        }

        self.running = true;
        self.sim_time = 0.0;
        self.plot_data.clear();
        self.plot_labels.clear();
        self.scope_data.clear();
        self.spectrum_data.clear();

        // Clear result history when starting a new simulation to prevent duplicate traces
        self.result_history.clear();

        // Reset timing statistics
        self.total_step_time = Duration::ZERO;
        self.step_count = 0;

        // Reset compiled simulation and sync parameters if in compiled mode
        #[cfg(not(target_arch = "wasm32"))]
        if self.use_compiled_mode {
            // Sync parameters from graph before running
            self.sync_compiled_params();

            if let Some(ref mut sim) = self.compiled_sim {
                // Reset time and internal states (but not params)
                sim.state_mut().time = 0.0;
                sim.state_mut().states.fill(0.0);
                sim.state_mut().outputs.fill(0.0);
            }
        }

        // Create block instances from graph
        self.create_blocks_from_graph();

        // Collect signal labels from Scope blocks
        for (_id, node) in &self.graph.nodes {
            if node.block_type == "Scope" {
                // Add labels from Scope input ports
                for input in &node.inputs {
                    self.plot_labels.push(input.name.clone());
                }
            }
        }

        // If no scope blocks, create default labels for source blocks
        if self.plot_labels.is_empty() {
            let order = self.get_execution_order();
            for node_id in &order {
                if let Some(node) = self.graph.get_node(node_id) {
                    if matches!(
                        node.block_type.as_str(),
                        "Constant" | "Sinusoidal" | "Step" | "Ramp" | "Integrator"
                    ) {
                        self.plot_labels
                            .push(format!("{} ({})", node.name, node.block_type));
                    }
                }
            }
        }
    }

    pub fn stop_simulation(&mut self) {
        self.running = false;

        // Save result to history when simulation completes
        if !self.plot_data.is_empty() {
            self.save_result_to_history();
        }
    }

    /// Sync parameters from graph to compiled simulation.
    /// Call this when parameters are changed in the UI before running the simulation.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn sync_compiled_params(&mut self) {
        if let Some(ref mut sim) = self.compiled_sim {
            let mapping = sim.param_mapping().clone();

            // Update each parameter from the graph
            for param in &mapping.params {
                if let Some(node) = self.graph.get_node(&param.node_id) {
                    if let Some(value) = node.get_param(&param.name).and_then(|v| v.as_f64()) {
                        sim.set_param(&param.node_id, &param.name, value);
                    }
                }
            }

            // Reinitialize blocks with new parameters
            sim.reinit();
        }
    }

    #[cfg(target_arch = "wasm32")]
    pub fn sync_compiled_params(&mut self) {
        // No-op on WASM
    }

    /// Save current simulation result to history for ghost traces
    fn save_result_to_history(&mut self) {
        let result = SimulationResult {
            plot_data: self.plot_data.clone(),
            plot_labels: self.plot_labels.clone(),
        };

        self.result_history.push_front(result);

        // Keep max 6 results
        while self.result_history.len() > 6 {
            self.result_history.pop_back();
        }
    }

    /// Clear result history
    pub fn clear_result_history(&mut self) {
        self.result_history.clear();
    }

    /// Create block instances from the current graph
    fn create_blocks_from_graph(&mut self) {
        self.blocks.clear();
        self.execution_order = self.get_execution_order();

        for node_id in &self.execution_order.clone() {
            if let Some(node) = self.graph.get_node(node_id) {
                let block = self.create_block_from_node(node);
                if let Some(b) = block {
                    self.blocks.insert(node_id.clone(), b);
                }
            }
        }
    }

    /// Create a BlockKind instance from a node
    fn create_block_from_node(&self, node: &NodeInstance) -> Option<BlockKind> {
        let block_type = node.block_type.as_str();

        Some(match block_type {
            "Constant" => {
                let value = node.get_param("value").and_then(|v| v.as_f64()).unwrap_or(1.0);
                BlockKind::Constant(Constant::new(value))
            }
            "Sinusoidal" => {
                let amp = node.get_param("amplitude").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let freq = node.get_param("frequency").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let phase = node.get_param("phase").and_then(|v| v.as_f64()).unwrap_or(0.0);
                BlockKind::Sinusoidal(Sinusoidal::new(amp, freq, phase))
            }
            "Step" => {
                let step_time = node.get_param("step_time").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let final_val = node.get_param("final_value").and_then(|v| v.as_f64()).unwrap_or(1.0);
                BlockKind::Step(Step::new(final_val, step_time))
            }
            "Ramp" => {
                let slope = node.get_param("slope").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let start = node.get_param("start_time").and_then(|v| v.as_f64()).unwrap_or(0.0);
                BlockKind::Ramp(Ramp::new(slope, start))
            }
            "Amplifier" | "Gain" => {
                let gain = node.get_param("gain")
                    .or_else(|| node.get_param("k"))
                    .and_then(|v| v.as_f64())
                    .unwrap_or(1.0);
                BlockKind::Amplifier(Amplifier::new(gain))
            }
            "Integrator" => {
                let initial = node.get_param("initial_value").and_then(|v| v.as_f64()).unwrap_or(0.0);
                BlockKind::Integrator(Integrator::new(initial))
            }
            "Differentiator" => {
                let tau = node.get_param("tau").and_then(|v| v.as_f64()).unwrap_or(0.01);
                BlockKind::Differentiator(Differentiator::new(tau))
            }
            "PID" => {
                let kp = node.get_param("kp").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let ki = node.get_param("ki").and_then(|v| v.as_f64()).unwrap_or(0.0);
                let kd = node.get_param("kd").and_then(|v| v.as_f64()).unwrap_or(0.0);
                BlockKind::PID(PID::new(kp, ki, kd))
            }
            "Adder" => {
                // Determine weights based on port names (+ or -)
                let n = node.inputs.len().max(2);
                let mut weights = vec![1.0; n];
                for (i, port) in node.inputs.iter().enumerate() {
                    if port.name.contains('-') {
                        weights[i] = -1.0;
                    }
                }
                match n {
                    2 => BlockKind::Adder2(Adder::<2>::with_weights([weights[0], weights[1]])),
                    3 => BlockKind::Adder3(Adder::<3>::with_weights([weights[0], weights[1], weights[2]])),
                    _ => BlockKind::Adder4(Adder::<4>::with_weights([weights[0], weights[1], weights[2], weights.get(3).copied().unwrap_or(1.0)])),
                }
            }
            "Multiplier" => {
                let n = node.inputs.len().max(2);
                match n {
                    2 => BlockKind::Multiplier2(Multiplier::<2>::new()),
                    3 => BlockKind::Multiplier3(Multiplier::<3>::new()),
                    _ => BlockKind::Multiplier4(Multiplier::<4>::new()),
                }
            }
            "Saturation" => {
                let min = node.get_param("min").and_then(|v| v.as_f64()).unwrap_or(-1.0);
                let max = node.get_param("max").and_then(|v| v.as_f64()).unwrap_or(1.0);
                BlockKind::Saturation(Saturation::new(min, max))
            }
            "Sin" => BlockKind::Sin(Sin::new()),
            "Cos" => BlockKind::Cos(Cos::new()),
            "Tan" => BlockKind::Tan(Tan::new()),
            "Abs" => BlockKind::Abs(Abs::new()),
            "Sqrt" => BlockKind::Sqrt(Sqrt::new()),
            "Exp" => BlockKind::Exp(Exp::new()),
            "Log" => BlockKind::Log(Log::new()),
            "Sign" => BlockKind::Sign(Sign::new()),
            "Pow" => {
                let exp = node.get_param("exponent").and_then(|v| v.as_f64()).unwrap_or(2.0);
                BlockKind::Pow(Pow::new(exp))
            }
            "Scope" => {
                let n = node.inputs.len().max(1);
                match n {
                    1 => BlockKind::Scope1(Scope::<1, 1000>::new()),
                    2 => BlockKind::Scope2(Scope::<2, 1000>::new()),
                    _ => BlockKind::Scope4(Scope::<4, 1000>::new()),
                }
            }
            "Pulse" => {
                let amp = node.get_param("amplitude").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let period = node.get_param("period").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let width = node.get_param("pulse_width").and_then(|v| v.as_f64()).unwrap_or(0.5);
                let rise = node.get_param("rise_time").and_then(|v| v.as_f64()).unwrap_or(0.0);
                let fall = node.get_param("fall_time").and_then(|v| v.as_f64()).unwrap_or(0.0);
                BlockKind::Pulse(Pulse::new(amp, period, width, rise, fall))
            }
            "SquareWave" => {
                let amp = node.get_param("amplitude").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let freq = node.get_param("frequency").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let duty = node.get_param("duty_cycle").and_then(|v| v.as_f64()).unwrap_or(0.5);
                BlockKind::SquareWave(SquareWave::new(amp, freq, duty))
            }
            "TriangleWave" => {
                let amp = node.get_param("amplitude").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let freq = node.get_param("frequency").and_then(|v| v.as_f64()).unwrap_or(1.0);
                BlockKind::TriangleWave(TriangleWave::new(amp, freq))
            }
            "Comparator" => {
                let threshold = node.get_param("threshold").and_then(|v| v.as_f64()).unwrap_or(0.0);
                BlockKind::Comparator(Comparator::new(threshold))
            }
            "Relay" => {
                let off_thresh = node.get_param("off_threshold").and_then(|v| v.as_f64()).unwrap_or(0.0);
                let on_thresh = node.get_param("on_threshold").and_then(|v| v.as_f64()).unwrap_or(0.5);
                BlockKind::Relay(Relay::new(off_thresh, on_thresh))
            }
            "LowpassRC" => {
                let cutoff = node.get_param("cutoff").and_then(|v| v.as_f64()).unwrap_or(10.0);
                BlockKind::LowpassRC(LowpassRC::new(cutoff))
            }
            "HighpassRC" => {
                let cutoff = node.get_param("cutoff").and_then(|v| v.as_f64()).unwrap_or(10.0);
                BlockKind::HighpassRC(HighpassRC::new(cutoff))
            }
            "RateLimiter" => {
                let rate = node.get_param("rate").and_then(|v| v.as_f64()).unwrap_or(1.0);
                BlockKind::RateLimiter(RateLimiter::new(rate))
            }
            "Min" => BlockKind::Min2(Min::<2>::new()),
            "Max" => BlockKind::Max2(Max::<2>::new()),
            _ => return None, // Unknown block type
        })
    }

    pub fn step_simulation(&mut self) {
        let step_start = Instant::now();
        let dt = self.settings.dt;

        // Use compiled simulation if available
        #[cfg(not(target_arch = "wasm32"))]
        if self.use_compiled_mode {
            if let Some(ref mut sim) = self.compiled_sim {
                // Record time BEFORE stepping (to match interpreter behavior)
                let t = self.sim_time;

                // Step the compiled simulation
                sim.step(dt);

                // Get all outputs from compiled simulation
                let all_outputs: Vec<f64> = (0..sim.state().outputs.len())
                    .map(|i| sim.get_output(i))
                    .collect();

                // Extract only Scope outputs (matching interpreter behavior)
                // The compiled sim outputs all blocks, but we only want Scope inputs for plotting
                // IMPORTANT: Use execution_order to match code generator's output order
                let mut scope_outputs: Vec<f64> = Vec::new();
                let mut output_idx = 0;

                // Build execution order to match code generator
                let order = self.execution_order.clone();

                // Iterate through nodes in topological order (matching code generator)
                for node_id in &order {
                    if let Some(node) = self.graph.get_node(node_id) {
                        match node.block_type.as_str() {
                            "Scope" | "Spectrum" => {
                                // Recording blocks - collect their outputs for plotting
                                let input_count = node.inputs.len().max(1);
                                let mut node_data = Vec::new();
                                for _ in 0..input_count {
                                    if output_idx < all_outputs.len() {
                                        let val = all_outputs[output_idx];
                                        scope_outputs.push(val);
                                        node_data.push(val);
                                        output_idx += 1;
                                    }
                                }
                                // Store in scope_data for preview
                                if !node_data.is_empty() {
                                    self.scope_data
                                        .entry(node_id.clone())
                                        .or_insert_with(Vec::new)
                                        .push((t, node_data));
                                }
                            }
                            _ => {
                                // Other blocks - skip their outputs (they're in all_outputs but we don't plot them)
                                output_idx += node.outputs.len().max(1);
                            }
                        }
                    }
                }

                // If no scope blocks, fall back to source block outputs
                if scope_outputs.is_empty() {
                    output_idx = 0;
                    for node_id in &order {
                        if let Some(node) = self.graph.get_node(node_id) {
                            if matches!(
                                node.block_type.as_str(),
                                "Constant" | "Sinusoidal" | "Step" | "Ramp" | "Integrator"
                            ) {
                                if output_idx < all_outputs.len() {
                                    scope_outputs.push(all_outputs[output_idx]);
                                }
                            }
                            // Advance output_idx for all blocks
                            match node.block_type.as_str() {
                                "Scope" | "Spectrum" => output_idx += node.inputs.len().max(1),
                                _ => output_idx += node.outputs.len().max(1),
                            }
                        }
                    }
                }

                self.plot_data.push((t, scope_outputs));
                self.sim_time += dt;

                // Record timing
                self.total_step_time += step_start.elapsed();
                self.step_count += 1;
                return;
            } else {
                // No compiled sim available, fall back to interpreter
                self.use_compiled_mode = false;
            }
        }

        // Interpreter mode using real block instances
        let t = self.sim_time;
        let mut outputs: Vec<f64> = Vec::new();

        // Process blocks in execution order
        let order = self.execution_order.clone();
        for node_id in &order {
            // Wire inputs from connected source blocks
            let connections: Vec<_> = self.graph.get_connections_to(node_id)
                .iter()
                .map(|c| (c.source_node_id.clone(), c.source_port_index, c.target_port_index))
                .collect();

            for (source_id, source_port, target_port) in connections {
                if let Some(source_block) = self.blocks.get(&source_id) {
                    let value = source_block.get_output(source_port);
                    if let Some(target_block) = self.blocks.get_mut(node_id) {
                        target_block.set_input(target_port, value);
                    }
                }
            }

            // Update and step the block
            if let Some(block) = self.blocks.get_mut(node_id) {
                block.update(t);
                block.step(t, dt);
            }
        }

        // Collect outputs from Scope blocks
        for node_id in &order {
            if let Some(node) = self.graph.get_node(node_id) {
                if node.block_type == "Scope" || node.block_type == "Spectrum" {
                    if let Some(block) = self.blocks.get(node_id) {
                        let mut node_data = Vec::new();
                        for i in 0..node.inputs.len().max(1) {
                            let val = block.get_output(i);
                            outputs.push(val);
                            node_data.push(val);
                        }
                        self.scope_data
                            .entry(node_id.clone())
                            .or_insert_with(Vec::new)
                            .push((t, node_data));
                    }
                }
            }
        }

        // If no scope blocks, record outputs from source blocks
        if outputs.is_empty() {
            for node_id in &order {
                if let Some(node) = self.graph.get_node(node_id) {
                    if matches!(
                        node.block_type.as_str(),
                        "Constant" | "Sinusoidal" | "Step" | "Ramp" | "Integrator"
                    ) {
                        if let Some(block) = self.blocks.get(node_id) {
                            outputs.push(block.get_output(0));
                        }
                    }
                }
            }
        }

        self.plot_data.push((t, outputs));
        self.sim_time += dt;

        // Record timing
        self.total_step_time += step_start.elapsed();
        self.step_count += 1;
    }

    /// Get average step execution time
    pub fn average_step_time(&self) -> Option<Duration> {
        if self.step_count > 0 {
            Some(self.total_step_time / self.step_count as u32)
        } else {
            None
        }
    }

    /// Get total step count
    pub fn get_step_count(&self) -> u64 {
        self.step_count
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
        let mut remaining: Vec<_> = self
            .graph
            .nodes
            .keys()
            .filter(|id| !result.contains(id))
            .cloned()
            .collect();
        remaining.sort();
        result.extend(remaining);

        result
    }
}

impl Default for AppState {
    fn default() -> Self {
        Self::new()
    }
}
