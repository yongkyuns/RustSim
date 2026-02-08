//! Plot panel UI component.

use egui::Ui;
use egui_plot::{Line, Plot, PlotPoints};

use crate::state::AppState;

/// Render the plots panel
pub fn render_plots(ui: &mut Ui, state: &mut AppState) {
    ui.horizontal(|ui| {
        ui.heading("Plots");
        ui.separator();

        if ui.button("Clear").clicked() {
            state.plot_data.clear();
        }

        ui.label(format!("{} points", state.plot_data.len()));
    });

    ui.separator();

    if state.plot_data.is_empty() {
        ui.centered_and_justified(|ui| {
            ui.label("Run simulation to see results");
        });
        return;
    }

    // Create plot
    Plot::new("simulation_plot")
        .height(ui.available_height() - 10.0)
        .show_axes(true)
        .show_grid(true)
        .legend(egui_plot::Legend::default())
        .show(ui, |plot_ui| {
            // Assume first output channel for now
            let points: PlotPoints = state
                .plot_data
                .iter()
                .filter_map(|(t, outputs)| {
                    outputs.first().map(|&y| [*t, y])
                })
                .collect();

            let line = Line::new(points)
                .name("Output 0")
                .width(2.0);

            plot_ui.line(line);
        });
}
