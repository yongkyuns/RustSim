//! Node instance types for simulation graphs.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::{PortInstance, SubsystemGraph};

/// Position in 2D space
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Position {
    pub x: f32,
    pub y: f32,
}

impl Position {
    pub fn new(x: f32, y: f32) -> Self {
        Self { x, y }
    }

    pub fn zero() -> Self {
        Self { x: 0.0, y: 0.0 }
    }
}

impl Default for Position {
    fn default() -> Self {
        Self::zero()
    }
}

/// A node instance in the simulation graph
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeInstance {
    /// Unique identifier for this node
    pub id: String,

    /// Block type (e.g., "Integrator", "Constant", "Adder")
    pub block_type: String,

    /// User-editable display name
    pub name: String,

    /// Position on the canvas
    pub position: Position,

    /// Input ports
    pub inputs: Vec<PortInstance>,

    /// Output ports
    pub outputs: Vec<PortInstance>,

    /// Parameter values (user-edited)
    pub params: HashMap<String, serde_json::Value>,

    /// Parameters to display inline on the node
    #[serde(default)]
    pub pinned_params: Vec<String>,

    /// Custom color override (hex string, e.g., "#ff0000")
    #[serde(skip_serializing_if = "Option::is_none")]
    pub color: Option<String>,

    /// Rotation in 90-degree increments (0, 1, 2, 3)
    #[serde(default)]
    pub rotation: u8,

    /// Nested graph for Subsystem nodes
    #[serde(skip_serializing_if = "Option::is_none")]
    pub graph: Option<Box<SubsystemGraph>>,
}

impl NodeInstance {
    /// Create a new node instance with default settings
    pub fn new(id: String, block_type: String, position: Position) -> Self {
        Self {
            id,
            block_type: block_type.clone(),
            name: block_type,
            position,
            inputs: Vec::new(),
            outputs: Vec::new(),
            params: HashMap::new(),
            pinned_params: Vec::new(),
            color: None,
            rotation: 0,
            graph: None,
        }
    }

    /// Add an input port
    pub fn add_input(&mut self, name: &str) {
        let index = self.inputs.len();
        self.inputs.push(PortInstance::new_input(&self.id, index, name));
    }

    /// Add an output port
    pub fn add_output(&mut self, name: &str) {
        let index = self.outputs.len();
        self.outputs
            .push(PortInstance::new_output(&self.id, index, name));
    }

    /// Remove the last input port, returns true if removed
    pub fn remove_input(&mut self) -> bool {
        if self.inputs.is_empty() {
            return false;
        }
        self.inputs.pop();
        true
    }

    /// Remove the last output port, returns true if removed
    pub fn remove_output(&mut self) -> bool {
        if self.outputs.is_empty() {
            return false;
        }
        self.outputs.pop();
        true
    }

    /// Set a parameter value
    pub fn set_param(&mut self, name: &str, value: serde_json::Value) {
        self.params.insert(name.to_string(), value);
    }

    /// Get a parameter value
    pub fn get_param(&self, name: &str) -> Option<&serde_json::Value> {
        self.params.get(name)
    }
}

/// An annotation (text label) on the canvas
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Annotation {
    /// Unique identifier
    pub id: String,

    /// Position on the canvas
    pub position: Position,

    /// Markdown content (supports LaTeX with $...$ and $$...$$)
    pub content: String,

    /// Width in pixels
    pub width: f32,

    /// Height in pixels
    pub height: f32,

    /// Custom color
    #[serde(skip_serializing_if = "Option::is_none")]
    pub color: Option<String>,

    /// Font size in pixels
    #[serde(skip_serializing_if = "Option::is_none")]
    pub font_size: Option<f32>,
}

impl Annotation {
    pub fn new(id: String, position: Position, content: String) -> Self {
        Self {
            id,
            position,
            content,
            width: 200.0,
            height: 100.0,
            color: None,
            font_size: None,
        }
    }
}
