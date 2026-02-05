//! Optimization utilities
//!
//! Provides Anderson acceleration, operators with linearization,
//! and numerical differentiation.

pub mod anderson;
mod operator;

pub use anderson::*;
pub use operator::*;
