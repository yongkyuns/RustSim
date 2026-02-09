//! Block type definitions and metadata.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::PortDefinition;

/// Category of a block type
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BlockCategory {
    Sources,
    Dynamic,
    Algebraic,
    Mixed,
    Recording,
    Subsystem,
}

impl BlockCategory {
    pub fn as_str(&self) -> &'static str {
        match self {
            BlockCategory::Sources => "Sources",
            BlockCategory::Dynamic => "Dynamic",
            BlockCategory::Algebraic => "Algebraic",
            BlockCategory::Mixed => "Mixed",
            BlockCategory::Recording => "Recording",
            BlockCategory::Subsystem => "Subsystem",
        }
    }

    /// Get the default shape for this category
    pub fn default_shape(&self) -> NodeShape {
        match self {
            BlockCategory::Sources => NodeShape::Pill,
            BlockCategory::Dynamic => NodeShape::Rect,
            BlockCategory::Algebraic => NodeShape::Rect,
            BlockCategory::Mixed => NodeShape::Mixed,
            BlockCategory::Recording => NodeShape::Pill,
            BlockCategory::Subsystem => NodeShape::Rect,
        }
    }
}

/// Visual shape of a node
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum NodeShape {
    /// Rounded ends (border-radius: 20px)
    Pill,
    /// Rectangle (border-radius: 4px)
    Rect,
    /// Circular (border-radius: 16px)
    Circle,
    /// Diamond shape
    Diamond,
    /// Asymmetric corners (12px 4px 12px 4px)
    Mixed,
    /// Right triangle pointing right (for amplifiers/gains)
    Triangle,
}

impl Default for NodeShape {
    fn default() -> Self {
        NodeShape::Rect
    }
}

/// Parameter type for block configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ParamType {
    Number,
    Integer,
    Boolean,
    String,
    Array,
    Callable,
    Any,
}

/// Parameter definition for a block type
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParamDefinition {
    /// Parameter type
    #[serde(rename = "type")]
    pub param_type: ParamType,

    /// Default value (serialized as string)
    pub default: Option<String>,

    /// Description of the parameter
    pub description: String,

    /// Minimum value (for numeric types)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min: Option<f64>,

    /// Maximum value (for numeric types)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max: Option<f64>,

    /// Allowed values (for enum-like strings)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub options: Option<Vec<String>>,
}

impl ParamDefinition {
    pub fn number(description: &str, default: f64) -> Self {
        Self {
            param_type: ParamType::Number,
            default: Some(default.to_string()),
            description: description.to_string(),
            min: None,
            max: None,
            options: None,
        }
    }

    pub fn integer(description: &str, default: i64) -> Self {
        Self {
            param_type: ParamType::Integer,
            default: Some(default.to_string()),
            description: description.to_string(),
            min: None,
            max: None,
            options: None,
        }
    }

    pub fn boolean(description: &str, default: bool) -> Self {
        Self {
            param_type: ParamType::Boolean,
            default: Some(default.to_string()),
            description: description.to_string(),
            min: None,
            max: None,
            options: None,
        }
    }

    pub fn with_range(mut self, min: f64, max: f64) -> Self {
        self.min = Some(min);
        self.max = Some(max);
        self
    }
}

/// Port configuration for a block type
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PortConfig {
    /// Input port definitions
    pub inputs: Vec<PortDefinition>,

    /// Output port definitions
    pub outputs: Vec<PortDefinition>,

    /// Minimum number of inputs
    #[serde(default = "default_one")]
    pub min_inputs: usize,

    /// Minimum number of outputs
    #[serde(default = "default_one")]
    pub min_outputs: usize,

    /// Maximum number of inputs (None = unlimited)
    pub max_inputs: Option<usize>,

    /// Maximum number of outputs (None = unlimited)
    pub max_outputs: Option<usize>,

    /// Whether output count should match input count
    #[serde(default)]
    pub sync_ports: bool,
}

fn default_one() -> usize {
    1
}

impl Default for PortConfig {
    fn default() -> Self {
        Self {
            inputs: vec![PortDefinition::input("in")],
            outputs: vec![PortDefinition::output("out")],
            min_inputs: 1,
            min_outputs: 1,
            max_inputs: None,
            max_outputs: None,
            sync_ports: false,
        }
    }
}

/// Definition of a block type
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlockTypeDefinition {
    /// Block type name (e.g., "Integrator")
    pub name: String,

    /// Short description
    pub description: String,

    /// Full documentation (HTML)
    #[serde(default)]
    pub doc_html: String,

    /// Category
    pub category: BlockCategory,

    /// Visual shape
    #[serde(default)]
    pub shape: NodeShape,

    /// Parameter definitions
    #[serde(default)]
    pub params: HashMap<String, ParamDefinition>,

    /// Port configuration
    #[serde(default)]
    pub ports: PortConfig,
}

impl BlockTypeDefinition {
    /// Create a simple block type with default ports
    pub fn new(name: &str, description: &str, category: BlockCategory) -> Self {
        Self {
            name: name.to_string(),
            description: description.to_string(),
            doc_html: String::new(),
            category,
            shape: category.default_shape(),
            params: HashMap::new(),
            ports: PortConfig::default(),
        }
    }

    /// Add a parameter to this block type
    pub fn with_param(mut self, name: &str, param: ParamDefinition) -> Self {
        self.params.insert(name.to_string(), param);
        self
    }

    /// Set the port configuration
    pub fn with_ports(mut self, ports: PortConfig) -> Self {
        self.ports = ports;
        self
    }

    /// Set the visual shape
    pub fn with_shape(mut self, shape: NodeShape) -> Self {
        self.shape = shape;
        self
    }
}
