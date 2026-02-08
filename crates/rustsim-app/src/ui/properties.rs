//! Properties panel UI component.

use egui::Ui;

use crate::state::AppState;

/// Render the properties panel
pub fn render_properties_panel(ui: &mut Ui, state: &mut AppState) {
    ui.heading("Properties");
    ui.separator();

    if state.selected_nodes.is_empty() {
        ui.label("Select a node to view properties");
        return;
    }

    // Show properties for first selected node
    let node_id = state.selected_nodes.iter().next().cloned();

    if let Some(node_id) = node_id {
        if let Some(node) = state.graph.get_node(&node_id) {
            ui.label(format!("Type: {}", node.block_type));
            ui.separator();

            // Name field
            ui.horizontal(|ui| {
                ui.label("Name:");
                // Note: Can't mutably borrow here, would need to restructure
                ui.label(&node.name);
            });

            ui.separator();
            ui.label("Parameters:");

            // Get block type definition for parameter info
            if let Some(block_def) = state.get_block_type(&node.block_type) {
                for (param_name, param_def) in &block_def.params {
                    ui.horizontal(|ui| {
                        ui.label(format!("{}:", param_name));

                        if let Some(value) = node.get_param(param_name) {
                            ui.label(value.to_string());
                        } else if let Some(default) = &param_def.default {
                            ui.label(format!("{} (default)", default));
                        }
                    });
                }
            }

            ui.separator();
            ui.label(format!("Inputs: {}", node.inputs.len()));
            ui.label(format!("Outputs: {}", node.outputs.len()));
            ui.label(format!(
                "Position: ({:.0}, {:.0})",
                node.position.x, node.position.y
            ));
        }
    }
}
