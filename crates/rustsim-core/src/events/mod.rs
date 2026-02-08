//! Event handling for hybrid systems
//!
//! Provides zero-crossing detection, condition-based events,
//! and scheduled events for discontinuous dynamics.

mod base;
mod condition;
mod schedule;
mod zerocrossing;

pub use base::*;
pub use condition::*;
pub use schedule::*;
pub use zerocrossing::*;
