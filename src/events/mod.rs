//! Event handling for hybrid systems
//!
//! Provides zero-crossing detection, condition-based events,
//! and scheduled events for discontinuous dynamics.

mod base;
mod zerocrossing;
mod condition;
mod schedule;

pub use base::*;
pub use zerocrossing::*;
pub use condition::*;
pub use schedule::*;
