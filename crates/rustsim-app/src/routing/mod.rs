//! Wire routing module for orthogonal connections between nodes.

mod astar;
mod orthogonal;

#[cfg(test)]
mod tests;

pub use astar::{AStarRouter, Direction, Rect};
pub use orthogonal::{OrthogonalRouter, RotationOptimizer};
