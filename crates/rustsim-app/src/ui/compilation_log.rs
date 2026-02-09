//! Compilation log viewer UI component.

use egui::ScrollArea;

use crate::state::AppState;

/// Render the compilation log window
pub fn render_compilation_log(ctx: &egui::Context, state: &mut AppState, open: &mut bool) {
    egui::Window::new("Compilation Log")
        .id(egui::Id::new("compilation_log_window"))
        .open(open)
        .resizable(true)
        .default_width(600.0)
        .default_height(400.0)
        .show(ctx, |ui| {
            ui.heading("Compilation Progress");

            ui.separator();

            // Show compilation status
            ui.horizontal(|ui| {
                ui.label("Status:");
                match &state.compilation_status {
                    crate::state::CompilationStatus::NotCompiled => {
                        ui.label("Not compiled");
                    }
                    crate::state::CompilationStatus::Compiling => {
                        ui.spinner();
                        ui.label("Compiling...");
                    }
                    crate::state::CompilationStatus::Ready => {
                        ui.label("[OK] Ready");
                    }
                    crate::state::CompilationStatus::Error(_) => {
                        ui.label("[ERR] Failed");
                    }
                }
            });

            ui.separator();

            // Log messages in a scrollable area
            ScrollArea::vertical()
                .auto_shrink([false; 2])
                .stick_to_bottom(true)
                .show(ui, |ui| {
                    if state.compilation_log.is_empty() {
                        ui.label("No compilation activity yet.");
                    } else {
                        for (i, msg) in state.compilation_log.iter().enumerate() {
                            ui.horizontal(|ui| {
                                ui.monospace(format!("[{}] {}", i + 1, msg));
                            });
                        }
                    }
                });

            ui.separator();

            // Action buttons
            ui.horizontal(|ui| {
                if ui.button("Clear Log").clicked() {
                    state.compilation_log.clear();
                }

                #[cfg(not(target_arch = "wasm32"))]
                {
                    use crate::state::CompilationStatus;

                    if matches!(
                        state.compilation_status,
                        CompilationStatus::NotCompiled | CompilationStatus::Error(_)
                    ) {
                        if ui.button("Compile").clicked() {
                            if let Err(e) = state.compile_simulation() {
                                state.compilation_log.push(format!("Failed to start compilation: {}", e));
                            }
                        }
                    }
                }
            });
        });
}
