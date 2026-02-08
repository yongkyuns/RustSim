//! Connection types for wiring blocks together.

use serde::{Deserialize, Serialize};

use crate::Position;

/// A connection between two nodes
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Connection {
    /// Unique identifier
    pub id: String,

    /// Source node ID
    pub source_node_id: String,

    /// Source port index
    pub source_port_index: usize,

    /// Target node ID
    pub target_node_id: String,

    /// Target port index
    pub target_port_index: usize,

    /// Optional user-placed waypoints for custom routing
    #[serde(default)]
    pub waypoints: Vec<Waypoint>,
}

impl Connection {
    /// Create a new connection
    pub fn new(
        id: String,
        source_node_id: String,
        source_port_index: usize,
        target_node_id: String,
        target_port_index: usize,
    ) -> Self {
        Self {
            id,
            source_node_id,
            source_port_index,
            target_node_id,
            target_port_index,
            waypoints: Vec::new(),
        }
    }

    /// Add a waypoint
    pub fn add_waypoint(&mut self, position: Position, is_user_waypoint: bool) {
        let id = format!("{}-wp-{}", self.id, self.waypoints.len());
        self.waypoints.push(Waypoint {
            id,
            position,
            is_user_waypoint,
        });
    }
}

/// A waypoint in a connection path
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Waypoint {
    /// Unique identifier
    pub id: String,

    /// Position on the canvas
    pub position: Position,

    /// Whether this waypoint was placed by the user (vs auto-calculated)
    pub is_user_waypoint: bool,
}
