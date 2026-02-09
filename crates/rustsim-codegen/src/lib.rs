//! Code generation for RustSim.
//!
//! Converts simulation graphs to executable Rust code.
//!
//! Three generators are available:
//! - `EmbeddedGenerator`: Zero-cost abstraction for embedded/production deployment
//! - `PathSimGenerator`: Clean runtime API using ArenaSimulation
//! - `CodeGenerator`: Full-featured with C ABI exports for dynamic loading

mod embedded_generator;
mod generator;
mod pathsim_generator;

pub use embedded_generator::EmbeddedGenerator;
pub use generator::*;
pub use pathsim_generator::PathSimGenerator;
