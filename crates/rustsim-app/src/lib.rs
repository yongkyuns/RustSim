//! RustSim application library

pub mod app;
pub mod compiler;
pub mod examples;
pub mod icons;
pub mod layout;
pub mod plotting;
pub mod routing;
pub mod spectrum_utils;
pub mod state;
pub mod subsystem;
pub mod ui;

// Re-export commonly used items
pub use app::RustSimApp;
pub use examples::{list_examples, load_example};
pub use routing::{AStarRouter, Direction, OrthogonalRouter, Rect};
pub use state::AppState;
pub use subsystem::{SubsystemEditor, SubsystemNavigator};
