//! Shared types for RustSim simulation graphs.
//!
//! This crate defines the core data structures used across all RustSim components:
//! - Node instances and their properties
//! - Connections between nodes
//! - Block type definitions and metadata
//! - Simulation settings and configuration

mod block;
mod connection;
mod graph;
mod node;
mod port;
mod settings;
mod simstate;

pub use block::*;
pub use connection::*;
pub use graph::*;
pub use node::*;
pub use port::*;
pub use settings::*;
pub use simstate::*;

/// Grid size for snapping positions (in pixels)
pub const GRID_SIZE: f32 = 10.0;

/// Port spacing (2 grid units)
pub const PORT_SPACING: f32 = 20.0;

/// Base node width
pub const NODE_BASE_WIDTH: f32 = 80.0;

/// Base node height
pub const NODE_BASE_HEIGHT: f32 = 40.0;

/// Snap a value to the grid
pub fn snap_to_grid(value: f32) -> f32 {
    (value / GRID_SIZE).round() * GRID_SIZE
}
