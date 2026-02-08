//! Toolbar UI component.

use egui::Ui;

use crate::state::AppState;

/// Render the top toolbar
pub fn render_toolbar(ui: &mut Ui, state: &mut AppState) {
    ui.horizontal_centered(|ui| {
        ui.add_space(8.0);

        // File operations group
        ui.horizontal(|ui| {
            if ui.button("📄 New").clicked() {
                state.new_file();
            }
            if ui.button("📂 Open").clicked() {
                state.open_file();
            }
            if ui.button("💾 Save").clicked() {
                state.save_file();
            }
        });

        ui.separator();

        // Simulation controls group
        ui.horizontal(|ui| {
            if state.is_running() {
                if ui
                    .button("⏹ Stop")
                    .on_hover_text("Stop simulation (Esc)")
                    .clicked()
                {
                    state.stop_simulation();
                }
            } else {
                if ui
                    .button("▶ Run")
                    .on_hover_text("Run simulation (Ctrl+Enter)")
                    .clicked()
                {
                    state.run_simulation();
                }
            }

            if ui
                .button("⏭ Step")
                .on_hover_text("Single step")
                .clicked()
            {
                state.step_simulation();
            }

            if ui.button("↺ Reset").on_hover_text("Reset to t=0").clicked() {
                state.sim_time = 0.0;
                state.plot_data.clear();
            }
        });

        ui.separator();

        // Status display
        ui.horizontal(|ui| {
            if state.is_running() {
                ui.spinner();
                ui.label(format!("t = {:.3}s", state.sim_time));
            } else {
                ui.label("Ready");
            }
        });

        ui.separator();

        // Edit operations
        ui.horizontal(|ui| {
            ui.add_enabled_ui(!state.selected_nodes.is_empty(), |ui| {
                if ui.button("🗑 Delete").clicked() {
                    state.delete_selected();
                }
                if ui.button("↻ Rotate").clicked() {
                    state.rotate_selected();
                }
            });
        });

        // Settings (right side)
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.add_space(8.0);

            // Solver selection
            egui::ComboBox::from_id_salt("solver")
                .selected_text(state.settings.solver.as_str())
                .show_ui(ui, |ui| {
                    use rustsim_types::SolverType;
                    ui.selectable_value(&mut state.settings.solver, SolverType::Euler, "Euler");
                    ui.selectable_value(&mut state.settings.solver, SolverType::Heun, "Heun");
                    ui.selectable_value(&mut state.settings.solver, SolverType::RK4, "RK4");
                    ui.selectable_value(&mut state.settings.solver, SolverType::RKCK54, "RKCK54");
                    ui.selectable_value(&mut state.settings.solver, SolverType::DOPRI54, "DOPRI54");
                });

            ui.label("Solver:");

            ui.separator();

            // Duration and dt
            ui.add(
                egui::DragValue::new(&mut state.settings.dt)
                    .speed(0.001)
                    .range(1e-6..=1.0)
                    .prefix("dt: "),
            );

            ui.add(
                egui::DragValue::new(&mut state.settings.duration)
                    .speed(0.1)
                    .range(0.1..=1000.0)
                    .prefix("T: ")
                    .suffix("s"),
            );
        });
    });
}
