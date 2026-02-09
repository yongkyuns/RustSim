//! UI components.

mod code_viewer;
mod compilation_log;
mod node_graph;
pub mod node_preview;
mod palette;
mod plot_settings;
mod plots;
mod properties;
mod spectrum_plot;
mod toolbar;

pub use code_viewer::render_code_viewer;
pub use compilation_log::render_compilation_log;
pub use node_graph::render_node_graph;
pub use node_preview::{
    decimate_minmax as preview_decimate, draw_plot_preview, draw_spectrum_preview,
};
pub use palette::render_block_palette;
pub use plot_settings::{
    get_default_color, LineStyle, MarkerStyle, PlotSettings, PlotViewMode, TraceStyle,
    DEFAULT_COLORS,
};
pub use plots::render_plots;
pub use properties::render_properties_panel;
pub use spectrum_plot::{render_spectrum_plot, render_spectrum_plot_linear};
pub use toolbar::render_toolbar;
