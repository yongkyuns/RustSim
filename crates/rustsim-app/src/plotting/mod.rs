//! Plot utilities module.

pub mod decimation;
pub mod export;

pub use decimation::decimate_minmax;
pub use export::{export_csv, PlotData};
