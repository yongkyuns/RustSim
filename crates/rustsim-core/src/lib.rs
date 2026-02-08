//! RustSim - Block-based time-domain simulation framework
//!
//! A Rust port of PathSim, providing high-performance simulation of dynamical systems
//! using block diagrams and numerical integration.
//!
//! # Architecture
//!
//! RustSim uses a compile-time static architecture with:
//! - Fixed I/O sizes at compile time (const generics)
//! - Graph assembled at compile time (manual struct or macro)
//! - Fully static dispatch - zero overhead
//!
//! # Example
//!
//! ```rust,ignore
//! use rustsim::*;
//!
//! // Define simulation as a struct
//! struct MySimulation {
//!     source: Constant,
//!     amp: Amplifier,
//!     int: Integrator,
//!     time: f64,
//! }
//!
//! impl MySimulation {
//!     fn step(&mut self, dt: f64) {
//!         self.source.update(self.time);
//!         self.amp.set_input(0, self.source.get_output(0));
//!         self.amp.update(self.time);
//!         self.int.set_input(0, self.amp.get_output(0));
//!         self.int.update(self.time);
//!         self.int.step(self.time, dt);
//!         self.time += dt;
//!     }
//! }
//! ```

// Core block trait and types
pub mod block;
pub mod blocks;
pub mod events;
pub mod optim;
pub mod solvers;
pub mod utils;

pub use block::{Block, DynamicBlock, StepResult};
pub use blocks::*;

/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::block::{Block, DynamicBlock, StepResult};
    pub use crate::blocks::*;
    pub use crate::solvers::*;
}
