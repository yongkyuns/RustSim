//! Code generation for RustSim.
//!
//! Converts simulation graphs to executable Rust code with C ABI exports
//! for dynamic loading.

mod generator;

pub use generator::*;
