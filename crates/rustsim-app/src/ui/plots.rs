//! Plot panel UI component.

use egui::{Color32, Ui};
use egui_plot::{Line, LineStyle as EguiLineStyle, MarkerShape, Plot, PlotPoints, Points};

use crate::plotting::{decimate_minmax, export_csv, PlotData};
use crate::state::AppState;
use crate::ui::{get_default_color, LineStyle, MarkerStyle};

const DECIMATION_THRESHOLD: usize = 1000;

/// Render the plots panel
pub fn render_plots(ui: &mut Ui, state: &mut AppState) {
    ui.horizontal(|ui| {
        ui.heading("Plots");
        ui.separator();

        if ui.button("Clear").clicked() {
            state.plot_data.clear();
        }

        if ui.button("Export CSV").clicked() {
            export_plot_data(state);
        }

        if ui.button("⚙ Settings").clicked() {
            state.editing_node = Some("plot_settings".to_string());
        }

        ui.label(format!("{} points", state.plot_data.len()));
    });

    ui.separator();

    // Show plot settings dialog if requested
    if state.editing_node == Some("plot_settings".to_string()) {
        render_plot_settings_dialog(ui, state);
    }

    if state.plot_data.is_empty() {
        ui.centered_and_justified(|ui| {
            ui.label("Run simulation to see results");
        });
        return;
    }

    // Determine the number of signals from the first data point
    let num_signals = state
        .plot_data
        .first()
        .map(|(_, outputs)| outputs.len())
        .unwrap_or(0);

    // Prepare plot data with log scale transformation if needed
    let plot_data = prepare_plot_data(state, num_signals);

    // Create plot with axis labels
    // Use available size to fill the panel properly
    let available = ui.available_size();
    let mut plot = Plot::new("simulation_plot")
        .width(available.x)
        .height(available.y.max(50.0))
        .show_axes(true)
        .x_axis_label(if state.plot_settings.x_axis_log {
            "Time [s] (log)"
        } else {
            "Time [s]"
        })
        .y_axis_label(if state.plot_settings.y_axis_log {
            "Amplitude (log)"
        } else {
            "Amplitude"
        });

    if state.plot_settings.show_grid {
        plot = plot.show_grid(true);
    }

    if state.plot_settings.show_legend {
        plot = plot.legend(egui_plot::Legend::default());
    }

    plot.show(ui, |plot_ui| {
        // Plot each signal channel
        for signal_idx in 0..num_signals {
            let trace_key = format!("signal-{}", signal_idx);

            // Get or create trace style
            let default_color = get_default_color(signal_idx);
            let style = state
                .plot_settings
                .get_or_create_trace_style(&trace_key, Some(default_color))
                .clone();

            // Skip if not visible
            if !style.visible {
                continue;
            }

            // Extract time and values for this signal
            let mut time_data: Vec<f64> = Vec::new();
            let mut signal_data: Vec<f64> = Vec::new();

            for (t, outputs) in &plot_data {
                if let Some(&y) = outputs.get(signal_idx) {
                    time_data.push(*t);
                    signal_data.push(y);
                }
            }

            // Apply decimation if we have too many points
            let (plot_time, plot_signal) = if time_data.len() > DECIMATION_THRESHOLD {
                decimate_minmax(&time_data, &signal_data, DECIMATION_THRESHOLD / 2)
            } else {
                (time_data, signal_data)
            };

            // Convert to plot points
            let points: PlotPoints = plot_time
                .iter()
                .zip(plot_signal.iter())
                .map(|(&t, &y)| [t, y])
                .collect();

            // Get signal name from labels or use default
            let signal_name = state
                .plot_labels
                .get(signal_idx)
                .cloned()
                .unwrap_or_else(|| format!("Signal {}", signal_idx));

            // Convert color to egui Color32
            let color = if let Some([r, g, b]) = style.color {
                Color32::from_rgb(r, g, b)
            } else {
                let [r, g, b] = default_color;
                Color32::from_rgb(r, g, b)
            };

            // Create line with style
            let mut line = Line::new(points).name(signal_name).color(color).width(2.0);

            // Apply line style
            line = match style.line_style {
                LineStyle::Solid => line.style(EguiLineStyle::Solid),
                LineStyle::Dash => line.style(EguiLineStyle::dashed_dense()),
                LineStyle::Dot => line.style(EguiLineStyle::dotted_dense()),
                LineStyle::None => line.style(EguiLineStyle::Solid).width(0.0),
            };

            // Apply marker style if not None
            if style.marker_style != MarkerStyle::None {
                let marker_shape = match style.marker_style {
                    MarkerStyle::Circle => MarkerShape::Circle,
                    MarkerStyle::Square => MarkerShape::Square,
                    MarkerStyle::Triangle => MarkerShape::Up,
                    MarkerStyle::None => MarkerShape::Circle, // fallback
                };
                // TODO: markers not available in this version of egui_plot
            }

            plot_ui.line(line);
        }
    });
}

/// Prepare plot data with log scale transformation if needed
fn prepare_plot_data(state: &AppState, _num_signals: usize) -> Vec<(f64, Vec<f64>)> {
    if !state.plot_settings.x_axis_log && !state.plot_settings.y_axis_log {
        return state.plot_data.clone();
    }

    state
        .plot_data
        .iter()
        .filter_map(|(t, outputs)| {
            let x = if state.plot_settings.x_axis_log {
                if *t > 0.0 {
                    t.log10()
                } else {
                    return None; // Skip non-positive values in log scale
                }
            } else {
                *t
            };

            let y_values: Vec<f64> = if state.plot_settings.y_axis_log {
                outputs
                    .iter()
                    .map(|&y| {
                        if y > 0.0 {
                            y.log10()
                        } else {
                            f64::NEG_INFINITY
                        }
                    })
                    .collect()
            } else {
                outputs.clone()
            };

            Some((x, y_values))
        })
        .collect()
}

/// Render plot settings dialog
fn render_plot_settings_dialog(ui: &mut Ui, state: &mut AppState) {
    egui::Window::new("Plot Settings")
        .collapsible(false)
        .resizable(true)
        .default_width(400.0)
        .show(ui.ctx(), |ui| {
            ui.heading("Global Settings");
            ui.separator();

            ui.horizontal(|ui| {
                ui.label("Show Legend:");
                ui.checkbox(&mut state.plot_settings.show_legend, "");
            });

            ui.horizontal(|ui| {
                ui.label("Show Grid:");
                ui.checkbox(&mut state.plot_settings.show_grid, "");
            });

            ui.horizontal(|ui| {
                ui.label("X-Axis Log Scale:");
                ui.checkbox(&mut state.plot_settings.x_axis_log, "");
            });

            ui.horizontal(|ui| {
                ui.label("Y-Axis Log Scale:");
                ui.checkbox(&mut state.plot_settings.y_axis_log, "");
            });

            ui.separator();
            ui.heading("Per-Signal Settings");
            ui.separator();

            // Determine the number of signals
            let num_signals = state
                .plot_data
                .first()
                .map(|(_, outputs)| outputs.len())
                .unwrap_or(0);

            if num_signals == 0 {
                ui.label("No signals to configure. Run a simulation first.");
            } else {
                egui::ScrollArea::vertical()
                    .max_height(300.0)
                    .show(ui, |ui| {
                        for signal_idx in 0..num_signals {
                            let trace_key = format!("signal-{}", signal_idx);

                            // Get or create trace style
                            let default_color = get_default_color(signal_idx);
                            let style = state
                                .plot_settings
                                .get_or_create_trace_style(&trace_key, Some(default_color));

                            // Get signal name
                            let signal_name = state
                                .plot_labels
                                .get(signal_idx)
                                .cloned()
                                .unwrap_or_else(|| format!("Signal {}", signal_idx));

                            ui.group(|ui| {
                                ui.horizontal(|ui| {
                                    ui.checkbox(&mut style.visible, "");
                                    ui.label(&signal_name);
                                });

                                ui.horizontal(|ui| {
                                    ui.label("Line Style:");
                                    egui::ComboBox::from_id_salt(format!(
                                        "line_style_{}",
                                        signal_idx
                                    ))
                                    .selected_text(style.line_style.name())
                                    .show_ui(ui, |ui| {
                                        for line_style in LineStyle::all() {
                                            ui.selectable_value(
                                                &mut style.line_style,
                                                line_style.clone(),
                                                line_style.name(),
                                            );
                                        }
                                    });
                                });

                                ui.horizontal(|ui| {
                                    ui.label("Marker Style:");
                                    egui::ComboBox::from_id_salt(format!(
                                        "marker_style_{}",
                                        signal_idx
                                    ))
                                    .selected_text(style.marker_style.name())
                                    .show_ui(ui, |ui| {
                                        for marker_style in MarkerStyle::all() {
                                            ui.selectable_value(
                                                &mut style.marker_style,
                                                marker_style.clone(),
                                                marker_style.name(),
                                            );
                                        }
                                    });
                                });

                                ui.horizontal(|ui| {
                                    ui.label("Color:");
                                    let mut color = if let Some([r, g, b]) = style.color {
                                        Color32::from_rgb(r, g, b)
                                    } else {
                                        let [r, g, b] = default_color;
                                        Color32::from_rgb(r, g, b)
                                    };

                                    let mut color_array = [color.r(), color.g(), color.b()];
                                    if ui.color_edit_button_srgb(&mut color_array).changed() {
                                        color = Color32::from_rgb(
                                            color_array[0],
                                            color_array[1],
                                            color_array[2],
                                        );
                                        style.color = Some(color_array);
                                    }

                                    if ui.button("Reset").clicked() {
                                        style.color = None;
                                    }
                                });
                            });

                            ui.add_space(5.0);
                        }
                    });
            }

            ui.separator();
            if ui.button("Close").clicked() {
                state.editing_node = None;
            }
        });
}

/// Export plot data to CSV
fn export_plot_data(state: &AppState) {
    if state.plot_data.is_empty() {
        return;
    }

    // Extract time and signals
    let time: Vec<f64> = state.plot_data.iter().map(|(t, _)| *t).collect();

    // Determine number of signals
    let num_signals = state
        .plot_data
        .first()
        .map(|(_, outputs)| outputs.len())
        .unwrap_or(0);

    // Extract each signal
    let mut signals: Vec<Vec<f64>> = vec![Vec::with_capacity(time.len()); num_signals];
    let mut labels: Vec<String> = Vec::with_capacity(num_signals);

    for i in 0..num_signals {
        // Use custom label if available, otherwise default
        labels.push(
            state
                .plot_labels
                .get(i)
                .cloned()
                .unwrap_or_else(|| format!("Signal {}", i)),
        );

        for (_, outputs) in &state.plot_data {
            if let Some(&value) = outputs.get(i) {
                signals[i].push(value);
            }
        }
    }

    // Create plot data
    let plot_data = PlotData::TimeDomain {
        time,
        signals,
        labels,
    };

    // Generate CSV
    let _csv = export_csv(&plot_data);

// TODO:     // Save or download the CSV
// TODO:     #[cfg(not(target_arch = "wasm32"))]
// TODO:     {
// TODO:         // Native: use file dialog
// TODO:         if let Some(path) = rfd::FileDialog::new()
// TODO:             .set_file_name("plot_data.csv")
// TODO:             .add_filter("CSV", &["csv"])
// TODO:             .save_file()
// TODO:         {
// TODO:             if let Err(e) = std::fs::write(&path, csv) {
// TODO:                 eprintln!("Failed to save CSV: {}", e);
// TODO:             }
// TODO:         }
// TODO:     }
// TODO: 
// TODO:     #[cfg(target_arch = "wasm32")]
// TODO:     {
// TODO:         // WASM: trigger browser download
// TODO:         use crate::plotting::export::download_csv_browser;
// TODO:         download_csv_browser(&csv, "plot_data.csv");
// TODO:     }
// TODO: Export function end
}
