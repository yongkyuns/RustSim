//! Toolbar UI component.

use egui::Ui;

use crate::state::AppState;

/// Render the top toolbar
pub fn render_toolbar(ui: &mut Ui, state: &mut AppState, show_compilation_log: &mut bool) {
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

        // Examples menu
        egui::menu::menu_button(ui, "📚 Examples", |ui| {
            for (name, description) in AppState::available_examples() {
                if ui.button(name).on_hover_text(description).clicked() {
                    state.load_example(name);
                    ui.close_menu();
                }
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

            if ui.button("⏭ Step").on_hover_text("Single step").clicked() {
                state.step_simulation();
            }

            if ui.button("↺ Reset").on_hover_text("Reset to t=0").clicked() {
                state.reset_to_zero();
            }
        });

        ui.separator();

        // Status display
        ui.horizontal(|ui| {
            if state.is_running() {
                ui.spinner();
                ui.label(format!("t = {:.3}s", state.sim_time()));
            } else if state.get_step_count() > 0 {
                ui.label(format!("t = {:.3}s", state.sim_time()));
            } else {
                ui.label("Ready");
            }

            // Show step timing if available
            if let Some(avg_time) = state.average_step_time() {
                let steps = state.get_step_count();
                let avg_us = avg_time.as_nanos() as f64 / 1000.0;
                #[cfg(not(target_arch = "wasm32"))]
                let mode = if state.use_compiled_mode() { "compiled" } else { "interp" };
                #[cfg(target_arch = "wasm32")]
                let mode = "interp";

                if avg_us < 1000.0 {
                    ui.label(format!("| {} steps @ {:.1}µs/step ({})", steps, avg_us, mode))
                        .on_hover_text("Average step execution time");
                } else {
                    let avg_ms = avg_us / 1000.0;
                    ui.label(format!("| {} steps @ {:.2}ms/step ({})", steps, avg_ms, mode))
                        .on_hover_text("Average step execution time");
                }
            }
        });

        ui.separator();

        // Edit operations
        ui.horizontal(|ui| {
            if ui.button("⚡ Auto Layout")
                .on_hover_text("Arrange blocks automatically")
                .clicked()
            {
                state.auto_layout();
            }

            ui.add_enabled_ui(!state.selected_nodes().is_empty(), |ui| {
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
                .selected_text(state.settings().solver.as_str())
                .show_ui(ui, |ui| {
                    use rustsim_types::SolverType;
                    ui.selectable_value(&mut state.settings_mut().solver, SolverType::Euler, "Euler");
                    ui.selectable_value(&mut state.settings_mut().solver, SolverType::Heun, "Heun");
                    ui.selectable_value(&mut state.settings_mut().solver, SolverType::RK4, "RK4");
                    ui.selectable_value(&mut state.settings_mut().solver, SolverType::RKCK54, "RKCK54");
                    ui.selectable_value(&mut state.settings_mut().solver, SolverType::DOPRI54, "DOPRI54");
                });

            ui.label("Solver:");

            ui.separator();

            // Compiled mode toggle (only on native platforms)
            #[cfg(not(target_arch = "wasm32"))]
            {
                use crate::state::CompilationStatus;

                // Show compilation status indicator
                match state.compilation_status() {
                    CompilationStatus::NotCompiled => {}
                    CompilationStatus::Compiling => {
                        ui.spinner();
                        ui.label("[...]");
                        // Show log window button when compiling
                        if ui.button("Show Log").clicked() {
                            *show_compilation_log = true;
                        }
                    }
                    CompilationStatus::Ready => {
                        if ui.button("[OK]")
                            .on_hover_text("Compilation successful - click to view log")
                            .clicked()
                        {
                            *show_compilation_log = true;
                        }
                    }
                    CompilationStatus::Error(err) => {
                        if ui.button("[ERR]")
                            .on_hover_text(format!("Compilation error: {} - click to view log", err))
                            .clicked()
                        {
                            *show_compilation_log = true;
                        }
                    }
                }

                // Compile button (shown when compiled mode is enabled but not compiled)
                if state.use_compiled_mode() && matches!(state.compilation_status(), CompilationStatus::NotCompiled | CompilationStatus::Error(_)) {
                    if ui.button("Compile")
                        .on_hover_text("Compile simulation to native code")
                        .clicked()
                    {
                        if let Err(e) = state.compile_simulation() {
                            eprintln!("Compilation failed: {}", e);
                        }
                        *show_compilation_log = true; // Auto-show log when starting compilation
                    }
                }

                // Compiled mode checkbox
                let was_compiled = state.use_compiled_mode();
                ui.checkbox(state.use_compiled_mode_mut(), "Compiled Mode")
                    .on_hover_text("Use compiled native code instead of interpreter (faster)");

                // If toggled on, trigger compilation
                if state.use_compiled_mode() && !was_compiled {
                    if let Err(e) = state.compile_simulation() {
                        eprintln!("Auto-compilation failed: {}", e);
                    }
                    *show_compilation_log = true; // Auto-show log when starting auto-compilation
                }
            }

            ui.separator();

            // Duration and dt
            ui.add(
                egui::DragValue::new(&mut state.settings_mut().dt)
                    .speed(0.001)
                    .range(1e-6..=1.0)
                    .prefix("dt: "),
            );

            ui.add(
                egui::DragValue::new(&mut state.settings_mut().duration)
                    .speed(0.1)
                    .range(0.1..=1000.0)
                    .prefix("T: ")
                    .suffix("s"),
            );
        });
    });
}
