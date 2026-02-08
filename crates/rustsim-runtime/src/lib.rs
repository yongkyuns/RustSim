//! Platform-abstracted simulation runtime for RustSim.
//!
//! This crate provides a unified interface for loading and running simulations
//! across different platforms:
//! - Native: Compiles to dynamic library and loads with libloading
//! - WASM: Sends to server for compilation, loads result as WASM module

mod platform;
mod plugin;

pub use platform::*;
pub use plugin::*;
