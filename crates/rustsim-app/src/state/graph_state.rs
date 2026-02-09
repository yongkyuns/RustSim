//! Graph state management - contains the simulation graph and block type definitions.

use rustsim_types::{
    BlockCategory, BlockTypeDefinition, Connection, NodeInstance, NodeShape, ParamDefinition,
    PortConfig, PortDefinition, Position, SimulationGraph,
};

use crate::subsystem::SubsystemNavigator;

/// State for managing the simulation graph
pub struct GraphState {
    /// The simulation graph
    pub graph: SimulationGraph,

    /// Block type definitions
    pub block_types: Vec<BlockTypeDefinition>,

    /// ID counter for generating unique IDs
    id_counter: u64,

    /// Subsystem navigator for hierarchical graphs
    pub subsystem_navigator: SubsystemNavigator,
}

impl GraphState {
    pub fn new() -> Self {
        Self {
            graph: SimulationGraph::new(),
            block_types: Self::create_block_types(),
            id_counter: 0,
            subsystem_navigator: SubsystemNavigator::new(),
        }
    }

    pub fn with_graph(graph: SimulationGraph) -> Self {
        Self {
            graph,
            block_types: Self::create_block_types(),
            id_counter: 0,
            subsystem_navigator: SubsystemNavigator::new(),
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
            })
            .with_shape(NodeShape::Triangle),
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
            // Subsystem blocks
            BlockTypeDefinition::new(
                "Subsystem",
                "Hierarchical subsystem containing nested blocks",
                BlockCategory::Subsystem,
            )
            .with_ports(PortConfig {
                inputs: vec![],
                outputs: vec![],
                min_inputs: 0,
                min_outputs: 0,
                max_inputs: None,
                max_outputs: None,
                sync_ports: false,
            }),
            BlockTypeDefinition::new(
                "Interface",
                "Subsystem boundary interface (internal use)",
                BlockCategory::Subsystem,
            )
            .with_ports(PortConfig {
                inputs: vec![],
                outputs: vec![],
                min_inputs: 0,
                min_outputs: 0,
                max_inputs: None,
                max_outputs: None,
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

    /// Remove a node from the graph
    pub fn remove_node(&mut self, node_id: &str) {
        self.graph.remove_node(node_id);
    }

    /// Add a connection between nodes
    pub fn add_connection(
        &mut self,
        source_node: &str,
        source_port: usize,
        target_node: &str,
        target_port: usize,
    ) -> String {
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

    /// Remove a connection from the graph
    pub fn remove_connection(&mut self, conn_id: &str) {
        self.graph.remove_connection(conn_id);
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

    /// Reset the graph to an empty state
    pub fn reset(&mut self) {
        self.graph = SimulationGraph::new();
        self.id_counter = 0;
    }
}

impl Default for GraphState {
    fn default() -> Self {
        Self::new()
    }
}
