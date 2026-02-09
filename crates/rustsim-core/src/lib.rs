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
pub mod anyblock;
pub mod arena_simulation;
#[cfg(test)]
mod arena_tests;
pub mod block;
pub mod block_kind;
pub mod blocks;
pub mod connection;
pub mod events;
#[macro_use]
pub mod macros;
pub mod optim;
pub mod simulation;
pub mod solvers;
pub mod utils;

pub use anyblock::AnyBlock;
pub use arena_simulation::{ArenaConnection, ArenaSimulation, BlockId};
pub use block::{Block, DynamicBlock, StepResult};
pub use block_kind::BlockKind;
pub use blocks::*;
pub use connection::{AsSource, AsTarget, BlockExt, Connection, ConnectionBuilder, ConnectionDef, Port, PortRef};
pub use simulation::{Simulation, SimulationBuilder};

/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::anyblock::AnyBlock;
    pub use crate::arena_simulation::{ArenaConnection, ArenaSimulation, BlockId};
    pub use crate::block::{Block, DynamicBlock, StepResult};
    pub use crate::block_kind::BlockKind;
    pub use crate::blocks::*;
    pub use crate::connection::{AsSource, AsTarget, BlockExt, Connection, ConnectionBuilder, ConnectionDef, Port, PortRef};
    pub use crate::simulation::{Simulation, SimulationBuilder};
    pub use crate::solvers::*;
}
