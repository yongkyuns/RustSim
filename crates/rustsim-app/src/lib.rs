//! RustSim application library

pub mod app;
pub mod compiler;
pub mod examples;
pub mod layout;
pub mod plotting;
pub mod spectrum_utils;
pub mod state;
pub mod ui;

// Re-export commonly used items
pub use app::RustSimApp;
pub use examples::{list_examples, load_example};
pub use state::AppState;
