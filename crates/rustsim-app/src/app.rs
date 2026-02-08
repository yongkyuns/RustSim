//! Main application structure.

use eframe::egui;

use crate::state::AppState;
use crate::ui::{
    render_block_palette, render_code_viewer, render_node_graph, render_plots,
    render_properties_panel, render_toolbar,
};

/// Main application
pub struct RustSimApp {
    /// Application state
    state: AppState,

    /// UI visibility flags
    show_palette: bool,
    show_properties: bool,
    show_plots: bool,
    show_code: bool,

    /// Panel sizes (for resizing)
    palette_width: f32,
    properties_width: f32,
    plots_height: f32,
    code_width: f32,
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
            show_code: false,
            palette_width: 200.0,
            properties_width: 300.0,
            plots_height: 200.0,
            code_width: 400.0,
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

    fn show_parameter_editor(&mut self, ctx: &egui::Context) {
        let editing_node = self.state.editing_node.clone();

        if let Some(node_id) = editing_node {
            let mut open = true;

            // Get node info for the title
            let title = if let Some(node) = self.state.graph.get_node(&node_id) {
                format!("{} - {}", node.name, node.block_type)
            } else {
                "Edit Parameters".to_string()
            };

            egui::Window::new(title)
                .id(egui::Id::new("param_editor"))
                .open(&mut open)
                .resizable(true)
                .default_width(300.0)
                .show(ctx, |ui| {
                    // Get block definition for parameter info
                    let block_def = if let Some(node) = self.state.graph.get_node(&node_id) {
                        self.state.get_block_type(&node.block_type).cloned()
                    } else {
                        None
                    };

                    if let Some(node) = self.state.graph.get_node_mut(&node_id) {
                        // Node name
                        ui.horizontal(|ui| {
                            ui.label("Name:");
                            ui.text_edit_singleline(&mut node.name);
                        });

                        ui.separator();
                        ui.heading("Parameters");

                        if let Some(ref block_def) = block_def {
                            for (param_name, param_def) in &block_def.params {
                                ui.horizontal(|ui| {
                                    // Label with description as tooltip
                                    let label = ui.label(format!("{}:", param_name));
                                    if !param_def.description.is_empty() {
                                        label.on_hover_text(&param_def.description);
                                    }

                                    // Get current value or default
                                    let current_val = node
                                        .get_param(param_name)
                                        .and_then(|v| v.as_f64())
                                        .or_else(|| {
                                            param_def.default.as_ref().and_then(|d| d.parse().ok())
                                        })
                                        .unwrap_or(0.0);

                                    let mut val = current_val;
                                    if ui
                                        .add(egui::DragValue::new(&mut val).speed(0.01))
                                        .changed()
                                    {
                                        node.set_param(param_name, serde_json::json!(val));
                                    }
                                });
                            }
                        }

                        if block_def.is_none() || block_def.as_ref().map_or(true, |b| b.params.is_empty()) {
                            ui.label("No parameters for this block");
                        }
                    }
                });

            if !open {
                self.state.editing_node = None;
            }
        }
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

        // Right panel - Code viewer (outermost right)
        if self.show_code {
            egui::SidePanel::right("code_viewer")
                .resizable(true)
                .default_width(self.code_width)
                .min_width(300.0)
                .max_width(800.0)
                .show(ctx, |ui| {
                    self.code_width = ui.available_width();
                    render_code_viewer(ui, &mut self.state);
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
                if ui
                    .selectable_label(self.show_code, "</>" )
                    .on_hover_text("Code Viewer")
                    .clicked()
                {
                    self.show_code = !self.show_code;
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

        // Parameter edit window (shown on double-click)
        self.show_parameter_editor(ctx);

        // Step simulation when running
        if self.state.is_running() {
            // Step the simulation
            self.state.step_simulation();

            // Stop when duration is reached
            if self.state.sim_time >= self.state.settings.duration {
                self.state.stop_simulation();
            }

            // Request continuous repaint
            ctx.request_repaint();
        }
    }
}
