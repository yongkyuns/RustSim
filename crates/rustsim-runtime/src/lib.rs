//! Platform-abstracted simulation runtime for RustSim.
//!
//! This crate provides a unified interface for loading and running simulations
//! across different platforms:
//! - Native: Compiles to dynamic library and loads with libloading
//! - WASM: Sends to server for compilation, loads result as WASM module
//!
//! It also provides backend abstractions for different execution strategies:
//! - Interpreter: Direct block-by-block execution
//! - Compiled: JIT-compiled native code execution

mod backend;
mod interpreter_backend;
mod platform;
mod plugin;

pub use backend::*;
pub use interpreter_backend::*;
pub use platform::*;
pub use plugin::*;
