//! Spectrum plot UI component for FFT visualization.

use egui::{Color32, Ui};
use egui_plot::{Line, Plot, PlotPoints};

use crate::spectrum_utils::{compute_fft, format_frequency};
use crate::state::AppState;
use crate::ui::get_default_color;

/// Render spectrum plot panel
pub fn render_spectrum_plot(ui: &mut Ui, state: &AppState) {
    ui.horizontal(|ui| {
        ui.heading("Spectrum Analysis");
        ui.separator();
        ui.label(format!("{} points", state.plot_data().len()));
    });

    ui.separator();

    if state.plot_data().is_empty() {
        ui.centered_and_justified(|ui| {
            ui.label("Run simulation to see spectrum");
        });
        return;
    }

    // Compute FFT on the plot data
    let sample_rate = 1.0 / state.settings().dt;
    let window_size = 1024; // Default FFT window size

    let (frequencies, magnitudes) = compute_fft(state.plot_data(), sample_rate, window_size);

    if frequencies.is_empty() {
        ui.centered_and_justified(|ui| {
            ui.label("No spectrum data available");
        });
        return;
    }

    // Determine the number of signals
    let num_signals = magnitudes.len();

    // Convert magnitudes to dB for better visualization
    let magnitudes_db: Vec<Vec<f64>> = magnitudes
        .iter()
        .map(|channel| {
            channel
                .iter()
                .map(|&mag| {
                    if mag > 0.0 {
                        20.0 * mag.log10()
                    } else {
                        -120.0 // Floor at -120 dB
                    }
                })
                .collect()
        })
        .collect();

    // Create spectrum plot
    let mut plot = Plot::new("spectrum_plot")
        .height(ui.available_height() - 10.0)
        .show_axes(true)
        .x_axis_label("Frequency [Hz]")
        .y_axis_label("Magnitude [dB]");

    if state.plot_settings().show_grid {
        plot = plot.show_grid(true);
    }

    if state.plot_settings().show_legend {
        plot = plot.legend(egui_plot::Legend::default());
    }

    plot.show(ui, |plot_ui| {
        // Plot each signal channel
        for signal_idx in 0..num_signals {
            let trace_key = format!("signal-{}", signal_idx);

            // Get or create trace style
            let default_color = get_default_color(signal_idx);
            let style = state
                .plot_settings()
                .get_trace_style(&trace_key)
                .cloned()
                .unwrap_or_else(|| {
                    // Create default style
                    crate::ui::TraceStyle {
                        visible: true,
                        color: Some(default_color),
                        line_style: crate::ui::LineStyle::Solid,
                        marker_style: crate::ui::MarkerStyle::None,
                    }
                });

            // Skip if not visible
            if !style.visible {
                continue;
            }

            // Extract spectrum points for this signal
            let points: PlotPoints = frequencies
                .iter()
                .zip(magnitudes_db[signal_idx].iter())
                .map(|(&freq, &mag_db)| [freq, mag_db])
                .collect();

            // Get signal name from labels or use default
            let signal_name = state
                .plot_labels()
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
            let line = Line::new(points)
                .name(format!("{} (Spectrum)", signal_name))
                .color(color)
                .width(2.0);

            plot_ui.line(line);
        }
    });

    // Display spectrum info
    ui.separator();
    ui.horizontal(|ui| {
        ui.label(format!("Sample Rate: {:.1} Hz", sample_rate));
        ui.label(format!("FFT Size: {}", window_size));
        let freq_resolution = sample_rate / window_size as f64;
        ui.label(format!("Frequency Resolution: {:.2} Hz", freq_resolution));
        let nyquist = sample_rate / 2.0;
        ui.label(format!(
            "Nyquist Frequency: {} Hz",
            format_frequency(nyquist)
        ));
    });
}

/// Render spectrum plot with linear magnitude scale
pub fn render_spectrum_plot_linear(ui: &mut Ui, state: &AppState) {
    ui.horizontal(|ui| {
        ui.heading("Spectrum Analysis (Linear)");
        ui.separator();
        ui.label(format!("{} points", state.plot_data().len()));
    });

    ui.separator();

    if state.plot_data().is_empty() {
        ui.centered_and_justified(|ui| {
            ui.label("Run simulation to see spectrum");
        });
        return;
    }

    // Compute FFT on the plot data
    let sample_rate = 1.0 / state.settings().dt;
    let window_size = 1024;

    let (frequencies, magnitudes) = compute_fft(state.plot_data(), sample_rate, window_size);

    if frequencies.is_empty() {
        ui.centered_and_justified(|ui| {
            ui.label("No spectrum data available");
        });
        return;
    }

    let num_signals = magnitudes.len();

    // Create spectrum plot with linear scale
    let mut plot = Plot::new("spectrum_plot_linear")
        .height(ui.available_height() - 10.0)
        .show_axes(true)
        .x_axis_label("Frequency [Hz]")
        .y_axis_label("Magnitude");

    if state.plot_settings().show_grid {
        plot = plot.show_grid(true);
    }

    if state.plot_settings().show_legend {
        plot = plot.legend(egui_plot::Legend::default());
    }

    plot.show(ui, |plot_ui| {
        for signal_idx in 0..num_signals {
            let trace_key = format!("signal-{}", signal_idx);
            let default_color = get_default_color(signal_idx);
            let style = state
                .plot_settings()
                .get_trace_style(&trace_key)
                .cloned()
                .unwrap_or_else(|| crate::ui::TraceStyle {
                    visible: true,
                    color: Some(default_color),
                    line_style: crate::ui::LineStyle::Solid,
                    marker_style: crate::ui::MarkerStyle::None,
                });

            if !style.visible {
                continue;
            }

            let points: PlotPoints = frequencies
                .iter()
                .zip(magnitudes[signal_idx].iter())
                .map(|(&freq, &mag)| [freq, mag])
                .collect();

            let signal_name = state
                .plot_labels()
                .get(signal_idx)
                .cloned()
                .unwrap_or_else(|| format!("Signal {}", signal_idx));

            let color = if let Some([r, g, b]) = style.color {
                Color32::from_rgb(r, g, b)
            } else {
                let [r, g, b] = default_color;
                Color32::from_rgb(r, g, b)
            };

            let line = Line::new(points)
                .name(format!("{} (Spectrum)", signal_name))
                .color(color)
                .width(2.0);

            plot_ui.line(line);
        }
    });
}
