//! Simulation graph types.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::{Annotation, Connection, NodeInstance, SimulationSettings};

/// A simulation graph containing nodes, connections, and annotations
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct SimulationGraph {
    /// All nodes in the graph, keyed by ID
    pub nodes: HashMap<String, NodeInstance>,

    /// All connections between nodes
    pub connections: Vec<Connection>,

    /// Annotations (text labels)
    #[serde(default)]
    pub annotations: HashMap<String, Annotation>,
}

impl SimulationGraph {
    /// Create an empty graph
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a node to the graph
    pub fn add_node(&mut self, node: NodeInstance) {
        self.nodes.insert(node.id.clone(), node);
    }

    /// Remove a node and all its connections
    pub fn remove_node(&mut self, node_id: &str) {
        self.nodes.remove(node_id);
        self.connections
            .retain(|c| c.source_node_id != node_id && c.target_node_id != node_id);
    }

    /// Add a connection
    pub fn add_connection(&mut self, connection: Connection) {
        self.connections.push(connection);
    }

    /// Remove a connection by ID
    pub fn remove_connection(&mut self, connection_id: &str) {
        self.connections.retain(|c| c.id != connection_id);
    }

    /// Get a node by ID
    pub fn get_node(&self, id: &str) -> Option<&NodeInstance> {
        self.nodes.get(id)
    }

    /// Get a mutable node by ID
    pub fn get_node_mut(&mut self, id: &str) -> Option<&mut NodeInstance> {
        self.nodes.get_mut(id)
    }

    /// Get connections from a specific node
    pub fn get_connections_from(&self, node_id: &str) -> Vec<&Connection> {
        self.connections
            .iter()
            .filter(|c| c.source_node_id == node_id)
            .collect()
    }

    /// Get connections to a specific node
    pub fn get_connections_to(&self, node_id: &str) -> Vec<&Connection> {
        self.connections
            .iter()
            .filter(|c| c.target_node_id == node_id)
            .collect()
    }
}

/// Nested graph for subsystems
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct SubsystemGraph {
    /// Nodes inside the subsystem
    pub nodes: Vec<NodeInstance>,

    /// Connections inside the subsystem
    pub connections: Vec<Connection>,

    /// Annotations
    #[serde(default)]
    pub annotations: Vec<Annotation>,
}

/// Complete file format for saving/loading simulations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GraphFile {
    /// File format version
    pub version: String,

    /// File metadata
    pub metadata: FileMetadata,

    /// The simulation graph
    pub graph: SimulationGraph,

    /// User-defined code context (expressions, functions)
    #[serde(default)]
    pub code_context: String,

    /// Simulation settings
    #[serde(default)]
    pub settings: SimulationSettings,
}

impl GraphFile {
    /// Current file format version
    pub const VERSION: &'static str = "1.0.0";

    /// Create a new graph file
    pub fn new(name: &str) -> Self {
        Self {
            version: Self::VERSION.to_string(),
            metadata: FileMetadata::new(name),
            graph: SimulationGraph::new(),
            code_context: String::new(),
            settings: SimulationSettings::default(),
        }
    }

    /// Load from JSON string
    pub fn from_json(json: &str) -> Result<Self, serde_json::Error> {
        serde_json::from_str(json)
    }

    /// Serialize to JSON string
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }
}

/// File metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileMetadata {
    /// File name
    pub name: String,

    /// Description
    #[serde(default)]
    pub description: String,

    /// Creation timestamp (ISO 8601)
    pub created: String,

    /// Last modified timestamp (ISO 8601)
    pub modified: String,
}

impl FileMetadata {
    pub fn new(name: &str) -> Self {
        let now = chrono_lite_now();
        Self {
            name: name.to_string(),
            description: String::new(),
            created: now.clone(),
            modified: now,
        }
    }
}

/// Simple timestamp function (avoids chrono dependency)
fn chrono_lite_now() -> String {
    // Return a placeholder - in real usage, use proper datetime
    "2024-01-01T00:00:00Z".to_string()
}
