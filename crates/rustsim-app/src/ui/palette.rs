//! Block palette UI component.

use egui::Ui;

use crate::state::AppState;
use rustsim_types::{BlockCategory, Position};

/// Render the block palette sidebar
pub fn render_block_palette(ui: &mut Ui, state: &mut AppState) {
    ui.heading("Blocks");
    ui.separator();

    // Search bar
    // TODO: Add search functionality

    // Collect block info first to avoid borrow issues
    let categories = [
        BlockCategory::Sources,
        BlockCategory::Dynamic,
        BlockCategory::Algebraic,
        BlockCategory::Mixed,
        BlockCategory::Recording,
    ];

    let blocks_by_category: Vec<_> = categories
        .iter()
        .map(|&category| {
            let blocks: Vec<_> = state
                .block_types()
                .iter()
                .filter(|b| b.category == category)
                .map(|b| (b.name.clone(), b.description.clone(), b.shape))
                .collect();
            (category, blocks)
        })
        .collect();

    let pan = *state.pan();

    egui::ScrollArea::vertical().show(ui, |ui| {
        for (category, blocks) in &blocks_by_category {
            if blocks.is_empty() {
                continue;
            }

            egui::CollapsingHeader::new(category.as_str())
                .default_open(true)
                .show(ui, |ui| {
                    ui.horizontal_wrapped(|ui| {
                        for (name, description, shape) in blocks {
                            let response = ui.add(BlockTile::new(name, description, *shape));

                            if response.clicked() {
                                // Add block at center of visible canvas area
                                // Account for node size to truly center it
                                let position = Position::new(
                                    -pan.x - rustsim_types::NODE_BASE_WIDTH / 2.0,
                                    -pan.y - rustsim_types::NODE_BASE_HEIGHT / 2.0,
                                );
                                state.add_node(name, position);
                            }
                        }
                    });
                });
        }
    });
}

/// Block tile widget for the palette
struct BlockTile<'a> {
    name: &'a str,
    description: &'a str,
    shape: rustsim_types::NodeShape,
}

impl<'a> BlockTile<'a> {
    fn new(name: &'a str, description: &'a str, shape: rustsim_types::NodeShape) -> Self {
        Self {
            name,
            description,
            shape,
        }
    }
}

impl<'a> egui::Widget for BlockTile<'a> {
    fn ui(self, ui: &mut Ui) -> egui::Response {
        let desired_size = egui::vec2(70.0, 50.0);

        let (rect, response) = ui.allocate_exact_size(desired_size, egui::Sense::click_and_drag());

        if ui.is_rect_visible(rect) {
            let visuals = ui.style().interact(&response);

            // Draw tile background
            let rounding = match self.shape {
                rustsim_types::NodeShape::Pill => egui::Rounding::same(10.0),
                rustsim_types::NodeShape::Rect => egui::Rounding::same(4.0),
                rustsim_types::NodeShape::Circle => egui::Rounding::same(8.0),
                rustsim_types::NodeShape::Mixed => egui::Rounding {
                    nw: 8.0,
                    ne: 4.0,
                    sw: 8.0,
                    se: 4.0,
                },
                _ => egui::Rounding::same(4.0),
            };

            let bg_color = if response.hovered() {
                visuals.bg_fill
            } else {
                ui.style().visuals.widgets.inactive.bg_fill
            };

            ui.painter().rect(
                rect,
                rounding,
                bg_color,
                egui::Stroke::new(1.0, visuals.fg_stroke.color),
            );

            // Draw block name
            ui.painter().text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                self.name,
                egui::FontId::proportional(11.0),
                visuals.text_color(),
            );
        }

        response.on_hover_text(self.description)
    }
}
