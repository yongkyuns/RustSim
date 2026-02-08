//! Port types for block inputs and outputs.

use serde::{Deserialize, Serialize};

/// Direction of a port (input or output)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum PortDirection {
    Input,
    Output,
}

impl PortDirection {
    pub fn as_str(&self) -> &'static str {
        match self {
            PortDirection::Input => "input",
            PortDirection::Output => "output",
        }
    }
}

/// Default port colors
pub mod port_colors {
    pub const DEFAULT: &str = "#969696";
    pub const SIGNAL: &str = "#64c8ff";
    pub const CONTROL: &str = "#ffc864";
    pub const DATA: &str = "#c864ff";
}

/// A port instance on a node
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PortInstance {
    /// Unique identifier (format: "{node_id}-{direction}-{index}")
    pub id: String,

    /// Parent node ID
    pub node_id: String,

    /// Display name
    pub name: String,

    /// Port direction (input or output)
    pub direction: PortDirection,

    /// Index within the node's inputs or outputs
    pub index: usize,

    /// Visual color (hex string)
    pub color: String,
}

impl PortInstance {
    /// Create a new input port
    pub fn new_input(node_id: &str, index: usize, name: &str) -> Self {
        Self {
            id: format!("{}-input-{}", node_id, index),
            node_id: node_id.to_string(),
            name: name.to_string(),
            direction: PortDirection::Input,
            index,
            color: port_colors::DEFAULT.to_string(),
        }
    }

    /// Create a new output port
    pub fn new_output(node_id: &str, index: usize, name: &str) -> Self {
        Self {
            id: format!("{}-output-{}", node_id, index),
            node_id: node_id.to_string(),
            name: name.to_string(),
            direction: PortDirection::Output,
            index,
            color: port_colors::DEFAULT.to_string(),
        }
    }

    /// Create a port with a specific color
    pub fn with_color(mut self, color: &str) -> Self {
        self.color = color.to_string();
        self
    }
}

/// Port definition in a block type (defines default port configuration)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PortDefinition {
    /// Default port name
    pub name: String,

    /// Direction
    pub direction: PortDirection,

    /// Default color
    #[serde(default = "default_port_color")]
    pub color: String,
}

fn default_port_color() -> String {
    port_colors::DEFAULT.to_string()
}

impl PortDefinition {
    pub fn input(name: &str) -> Self {
        Self {
            name: name.to_string(),
            direction: PortDirection::Input,
            color: port_colors::DEFAULT.to_string(),
        }
    }

    pub fn output(name: &str) -> Self {
        Self {
            name: name.to_string(),
            direction: PortDirection::Output,
            color: port_colors::DEFAULT.to_string(),
        }
    }
}
