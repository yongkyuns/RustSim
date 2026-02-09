//! Canvas state management - contains canvas view and interaction state.

use std::collections::HashSet;
use rustsim_types::Position;

/// Pending connection being drawn
pub struct PendingConnection {
    pub source_node: String,
    pub source_port: usize,
    pub current_pos: Position,
}

/// State for managing canvas view and interaction
pub struct CanvasState {
    /// Canvas pan offset
    pub pan: Position,

    /// Canvas zoom level
    pub zoom: f32,

    /// Whether to show grid on canvas
    pub show_grid: bool,

    /// Selected node IDs
    pub selected_nodes: HashSet<String>,

    /// Selected connection IDs
    pub selected_connections: HashSet<String>,

    /// Currently dragging node
    pub dragging_node: Option<String>,

    /// Connection being created
    pub pending_connection: Option<PendingConnection>,
}

impl CanvasState {
    pub fn new() -> Self {
        Self {
            pan: Position::zero(),
            zoom: 1.0,
            show_grid: true,
            selected_nodes: HashSet::new(),
            selected_connections: HashSet::new(),
            dragging_node: None,
            pending_connection: None,
        }
    }

    /// Clear all selections
    pub fn clear_selection(&mut self) {
        self.selected_nodes.clear();
        self.selected_connections.clear();
    }

    /// Select a node
    pub fn select_node(&mut self, node_id: String) {
        self.selected_nodes.insert(node_id);
    }

    /// Deselect a node
    pub fn deselect_node(&mut self, node_id: &str) {
        self.selected_nodes.remove(node_id);
    }

    /// Select a connection
    pub fn select_connection(&mut self, conn_id: String) {
        self.selected_connections.insert(conn_id);
    }

    /// Deselect a connection
    pub fn deselect_connection(&mut self, conn_id: &str) {
        self.selected_connections.remove(conn_id);
    }

    /// Check if a node is selected
    pub fn is_node_selected(&self, node_id: &str) -> bool {
        self.selected_nodes.contains(node_id)
    }

    /// Check if a connection is selected
    pub fn is_connection_selected(&self, conn_id: &str) -> bool {
        self.selected_connections.contains(conn_id)
    }

    /// Zoom in
    pub fn zoom_in(&mut self) {
        self.zoom = (self.zoom * 1.2).min(4.0);
    }

    /// Zoom out
    pub fn zoom_out(&mut self) {
        self.zoom = (self.zoom / 1.2).max(0.1);
    }

    /// Reset zoom
    pub fn reset_zoom(&mut self) {
        self.zoom = 1.0;
    }

    /// Pan the canvas
    pub fn pan_by(&mut self, dx: f32, dy: f32) {
        self.pan.x += dx;
        self.pan.y += dy;
    }

    /// Reset pan
    pub fn reset_pan(&mut self) {
        self.pan = Position::zero();
    }

    /// Start dragging a node
    pub fn start_dragging(&mut self, node_id: String) {
        self.dragging_node = Some(node_id);
    }

    /// Stop dragging
    pub fn stop_dragging(&mut self) {
        self.dragging_node = None;
    }

    /// Start creating a connection
    pub fn start_connection(&mut self, source_node: String, source_port: usize, pos: Position) {
        self.pending_connection = Some(PendingConnection {
            source_node,
            source_port,
            current_pos: pos,
        });
    }

    /// Update pending connection position
    pub fn update_connection(&mut self, pos: Position) {
        if let Some(ref mut conn) = self.pending_connection {
            conn.current_pos = pos;
        }
    }

    /// Cancel pending connection
    pub fn cancel_connection(&mut self) {
        self.pending_connection = None;
    }

    /// Reset canvas state
    pub fn reset(&mut self) {
        self.pan = Position::zero();
        self.zoom = 1.0;
        self.selected_nodes.clear();
        self.selected_connections.clear();
        self.dragging_node = None;
        self.pending_connection = None;
    }
}

impl Default for CanvasState {
    fn default() -> Self {
        Self::new()
    }
}
