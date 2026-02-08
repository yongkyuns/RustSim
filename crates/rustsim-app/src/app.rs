//! Main application structure.

use eframe::egui;

use crate::state::AppState;
use crate::ui::{
    render_block_palette, render_node_graph, render_plots, render_properties_panel,
    render_toolbar,
};

/// Main application
pub struct RustSimApp {
    /// Application state
    state: AppState,

    /// UI visibility flags
    show_palette: bool,
    show_properties: bool,
    show_plots: bool,

    /// Panel sizes (for resizing)
    palette_width: f32,
    properties_width: f32,
    plots_height: f32,
}

impl RustSimApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        // Configure fonts and style
        Self::configure_style(&cc.egui_ctx);

        Self {
            state: AppState::new(),
            show_palette: true,
            show_properties: false,
            show_plots: true,
            palette_width: 200.0,
            properties_width: 300.0,
            plots_height: 200.0,
        }
    }

    fn configure_style(ctx: &egui::Context) {
        let mut style = (*ctx.style()).clone();

        // Use a more modern look
        style.visuals.window_rounding = egui::Rounding::same(8.0);
        style.visuals.menu_rounding = egui::Rounding::same(4.0);
        style.visuals.popup_shadow = egui::epaint::Shadow::NONE;

        ctx.set_style(style);
    }

    fn handle_keyboard(&mut self, ctx: &egui::Context) {
        let input = ctx.input(|i| i.clone());

        // Ctrl/Cmd shortcuts
        if input.modifiers.command {
            if input.key_pressed(egui::Key::S) {
                self.state.save_file();
            }
            if input.key_pressed(egui::Key::O) {
                self.state.open_file();
            }
            if input.key_pressed(egui::Key::N) {
                self.state.new_file();
            }
            if input.key_pressed(egui::Key::Z) {
                if input.modifiers.shift {
                    self.state.redo();
                } else {
                    self.state.undo();
                }
            }
            if input.key_pressed(egui::Key::Enter) {
                self.state.run_simulation();
            }
        }

        // Single key shortcuts (when not typing)
        if !ctx.wants_keyboard_input() {
            if input.key_pressed(egui::Key::Escape) {
                if self.state.is_running() {
                    self.state.stop_simulation();
                } else {
                    self.state.clear_selection();
                }
            }
            if input.key_pressed(egui::Key::Delete) || input.key_pressed(egui::Key::Backspace) {
                self.state.delete_selected();
            }
            if input.key_pressed(egui::Key::R) {
                self.state.rotate_selected();
            }
            if input.key_pressed(egui::Key::F) {
                self.state.fit_view();
            }
        }
    }
}

impl eframe::App for RustSimApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Handle keyboard shortcuts
        self.handle_keyboard(ctx);

        // Top toolbar
        egui::TopBottomPanel::top("toolbar")
            .exact_height(48.0)
            .show(ctx, |ui| {
                render_toolbar(ui, &mut self.state);
            });

        // Left panel - Block palette
        if self.show_palette {
            egui::SidePanel::left("palette")
                .resizable(true)
                .default_width(self.palette_width)
                .min_width(150.0)
                .max_width(350.0)
                .show(ctx, |ui| {
                    self.palette_width = ui.available_width();
                    render_block_palette(ui, &mut self.state);
                });
        }

        // Right panel - Properties
        if self.show_properties {
            egui::SidePanel::right("properties")
                .resizable(true)
                .default_width(self.properties_width)
                .min_width(200.0)
                .max_width(500.0)
                .show(ctx, |ui| {
                    self.properties_width = ui.available_width();
                    render_properties_panel(ui, &mut self.state);
                });
        }

        // Bottom panel - Plots
        if self.show_plots {
            egui::TopBottomPanel::bottom("plots")
                .resizable(true)
                .default_height(self.plots_height)
                .min_height(100.0)
                .max_height(500.0)
                .show(ctx, |ui| {
                    self.plots_height = ui.available_height();
                    render_plots(ui, &mut self.state);
                });
        }

        // Central panel - Node graph
        egui::CentralPanel::default().show(ctx, |ui| {
            // Panel toggle buttons at the edge
            ui.horizontal(|ui| {
                if ui
                    .selectable_label(self.show_palette, "📦")
                    .on_hover_text("Block Palette")
                    .clicked()
                {
                    self.show_palette = !self.show_palette;
                }
                if ui
                    .selectable_label(self.show_properties, "⚙")
                    .on_hover_text("Properties")
                    .clicked()
                {
                    self.show_properties = !self.show_properties;
                }
                if ui
                    .selectable_label(self.show_plots, "📊")
                    .on_hover_text("Plots")
                    .clicked()
                {
                    self.show_plots = !self.show_plots;
                }

                ui.separator();

                // Zoom controls
                ui.label(format!("Zoom: {:.0}%", self.state.zoom * 100.0));
                if ui.button("−").clicked() {
                    self.state.zoom_out();
                }
                if ui.button("+").clicked() {
                    self.state.zoom_in();
                }
                if ui.button("Fit").clicked() {
                    self.state.fit_view();
                }
            });

            ui.separator();

            // Node graph canvas
            render_node_graph(ui, &mut self.state);
        });

        // Continuous repaint during simulation
        if self.state.is_running() {
            ctx.request_repaint();
        }
    }
}
