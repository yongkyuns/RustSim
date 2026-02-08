//! Node graph canvas UI component.

use egui::{Color32, Pos2, Rect, Sense, Stroke, Ui, Vec2};

use crate::state::AppState;
use rustsim_types::{snap_to_grid, NodeShape, Position, NODE_BASE_HEIGHT, NODE_BASE_WIDTH};

/// Render the node graph canvas
pub fn render_node_graph(ui: &mut Ui, state: &mut AppState) {
    let (response, painter) =
        ui.allocate_painter(ui.available_size(), Sense::click_and_drag());

    let canvas_rect = response.rect;
    let canvas_center = canvas_rect.center();

    // Handle pan
    if response.dragged_by(egui::PointerButton::Middle)
        || (response.dragged_by(egui::PointerButton::Primary) && ui.input(|i| i.modifiers.shift))
    {
        let delta = response.drag_delta();
        state.pan.x += delta.x / state.zoom;
        state.pan.y += delta.y / state.zoom;
    }

    // Handle zoom with scroll
    if response.hovered() {
        let scroll_delta = ui.input(|i| i.raw_scroll_delta.y);
        if scroll_delta != 0.0 {
            let zoom_factor = 1.0 + scroll_delta * 0.001;
            state.zoom = (state.zoom * zoom_factor).clamp(0.1, 4.0);
        }
    }

    // Transform functions
    let to_screen = |pos: Position| -> Pos2 {
        Pos2::new(
            canvas_center.x + (pos.x + state.pan.x) * state.zoom,
            canvas_center.y + (pos.y + state.pan.y) * state.zoom,
        )
    };

    let to_canvas = |pos: Pos2| -> Position {
        Position::new(
            (pos.x - canvas_center.x) / state.zoom - state.pan.x,
            (pos.y - canvas_center.y) / state.zoom - state.pan.y,
        )
    };

    // Draw grid
    draw_grid(&painter, canvas_rect, state.pan, state.zoom);

    // Draw connections
    for conn in &state.graph.connections {
        let source_node = state.graph.get_node(&conn.source_node_id);
        let target_node = state.graph.get_node(&conn.target_node_id);

        if let (Some(source), Some(target)) = (source_node, target_node) {
            let source_pos = get_output_port_position(source, conn.source_port_index);
            let target_pos = get_input_port_position(target, conn.target_port_index);

            let source_screen = to_screen(source_pos);
            let target_screen = to_screen(target_pos);

            let is_selected = state.selected_connections.contains(&conn.id);
            let color = if is_selected {
                Color32::from_rgb(100, 200, 255)
            } else {
                Color32::from_rgb(150, 150, 150)
            };

            // Draw bezier curve
            draw_connection(&painter, source_screen, target_screen, color);
        }
    }

    // Draw pending connection
    if let Some(pending) = &state.pending_connection {
        if let Some(source_node) = state.graph.get_node(&pending.source_node) {
            let source_pos = get_output_port_position(source_node, pending.source_port);
            let source_screen = to_screen(source_pos);
            let target_screen = to_screen(pending.current_pos);

            draw_connection(&painter, source_screen, target_screen, Color32::from_rgb(100, 200, 255));
        }
    }

    // Draw nodes
    let mut clicked_node: Option<String> = None;
    let mut drag_delta = Vec2::ZERO;

    for (id, node) in &state.graph.nodes {
        let node_pos = to_screen(node.position);
        let node_size = Vec2::new(NODE_BASE_WIDTH * state.zoom, NODE_BASE_HEIGHT * state.zoom);
        let node_rect = Rect::from_min_size(node_pos, node_size);

        let is_selected = state.selected_nodes.contains(id);

        // Get block shape
        let shape = state
            .get_block_type(&node.block_type)
            .map(|b| b.shape)
            .unwrap_or(NodeShape::Rect);

        // Draw node
        draw_node(
            &painter,
            node_rect,
            &node.name,
            shape,
            is_selected,
            node.rotation,
        );

        // Draw ports
        let input_count = node.inputs.len();
        let output_count = node.outputs.len();

        for i in 0..input_count {
            let port_pos = to_screen(get_input_port_position(node, i));
            draw_port(&painter, port_pos, true, false);
        }

        for i in 0..output_count {
            let port_pos = to_screen(get_output_port_position(node, i));
            draw_port(&painter, port_pos, false, false);
        }

        // Handle node interaction
        let node_response = ui.interact(node_rect, egui::Id::new(id), Sense::click_and_drag());

        if node_response.clicked() {
            clicked_node = Some(id.clone());
        }

        if node_response.dragged() && state.dragging_node.as_ref() == Some(id) {
            drag_delta = node_response.drag_delta();
        }

        if node_response.drag_started() {
            state.dragging_node = Some(id.clone());
        }

        if node_response.drag_stopped() {
            state.dragging_node = None;
        }
    }

    // Handle node selection
    if let Some(node_id) = clicked_node {
        let ctrl = ui.input(|i| i.modifiers.ctrl || i.modifiers.command);
        if ctrl {
            // Toggle selection
            if state.selected_nodes.contains(&node_id) {
                state.selected_nodes.remove(&node_id);
            } else {
                state.selected_nodes.insert(node_id);
            }
        } else {
            // Single selection
            state.selected_nodes.clear();
            state.selected_nodes.insert(node_id);
        }
    } else if response.clicked() {
        // Clicked on empty space - clear selection
        state.clear_selection();
    }

    // Update dragged node positions
    if drag_delta != Vec2::ZERO {
        for node_id in &state.selected_nodes {
            if let Some(node) = state.graph.get_node_mut(node_id) {
                node.position.x += drag_delta.x / state.zoom;
                node.position.y += drag_delta.y / state.zoom;
                // Snap to grid on release
                if state.dragging_node.is_none() {
                    node.position.x = snap_to_grid(node.position.x);
                    node.position.y = snap_to_grid(node.position.y);
                }
            }
        }
    }
}

fn draw_grid(painter: &egui::Painter, rect: Rect, pan: Position, zoom: f32) {
    let grid_size = 20.0 * zoom;
    let grid_color = Color32::from_rgba_unmultiplied(128, 128, 128, 30);

    let offset_x = (pan.x * zoom) % grid_size;
    let offset_y = (pan.y * zoom) % grid_size;

    // Vertical lines
    let mut x = rect.min.x + offset_x;
    while x < rect.max.x {
        painter.line_segment(
            [Pos2::new(x, rect.min.y), Pos2::new(x, rect.max.y)],
            Stroke::new(1.0, grid_color),
        );
        x += grid_size;
    }

    // Horizontal lines
    let mut y = rect.min.y + offset_y;
    while y < rect.max.y {
        painter.line_segment(
            [Pos2::new(rect.min.x, y), Pos2::new(rect.max.x, y)],
            Stroke::new(1.0, grid_color),
        );
        y += grid_size;
    }
}

fn draw_node(
    painter: &egui::Painter,
    rect: Rect,
    name: &str,
    shape: NodeShape,
    selected: bool,
    _rotation: u8,
) {
    let rounding = match shape {
        NodeShape::Pill => egui::Rounding::same(rect.height() / 2.0),
        NodeShape::Rect => egui::Rounding::same(4.0),
        NodeShape::Circle => egui::Rounding::same(rect.height() / 2.0),
        NodeShape::Diamond => egui::Rounding::same(4.0),
        NodeShape::Mixed => egui::Rounding {
            nw: 12.0,
            ne: 4.0,
            sw: 12.0,
            se: 4.0,
        },
    };

    let bg_color = Color32::from_rgb(45, 45, 48);
    let border_color = if selected {
        Color32::from_rgb(100, 200, 255)
    } else {
        Color32::from_rgb(80, 80, 80)
    };

    // Selection glow
    if selected {
        painter.rect(
            rect.expand(3.0),
            rounding,
            Color32::from_rgba_unmultiplied(100, 200, 255, 50),
            Stroke::NONE,
        );
    }

    // Node background
    painter.rect(rect, rounding, bg_color, Stroke::new(1.0, border_color));

    // Node name
    painter.text(
        rect.center(),
        egui::Align2::CENTER_CENTER,
        name,
        egui::FontId::proportional(12.0),
        Color32::WHITE,
    );
}

fn draw_port(painter: &egui::Painter, pos: Pos2, is_input: bool, is_hovered: bool) {
    let radius = 4.0;
    let color = if is_hovered {
        Color32::from_rgb(100, 200, 255)
    } else {
        Color32::from_rgb(150, 150, 150)
    };

    // Hollow circle for ports
    painter.circle_stroke(pos, radius, Stroke::new(2.0, color));
}

fn draw_connection(painter: &egui::Painter, from: Pos2, to: Pos2, color: Color32) {
    // Simple bezier curve
    let ctrl_dist = ((to.x - from.x).abs() / 2.0).max(30.0);

    let ctrl1 = Pos2::new(from.x + ctrl_dist, from.y);
    let ctrl2 = Pos2::new(to.x - ctrl_dist, to.y);

    // Draw cubic bezier
    let points: Vec<Pos2> = (0..=20)
        .map(|i| {
            let t = i as f32 / 20.0;
            let t2 = t * t;
            let t3 = t2 * t;
            let mt = 1.0 - t;
            let mt2 = mt * mt;
            let mt3 = mt2 * mt;

            Pos2::new(
                mt3 * from.x + 3.0 * mt2 * t * ctrl1.x + 3.0 * mt * t2 * ctrl2.x + t3 * to.x,
                mt3 * from.y + 3.0 * mt2 * t * ctrl1.y + 3.0 * mt * t2 * ctrl2.y + t3 * to.y,
            )
        })
        .collect();

    for window in points.windows(2) {
        painter.line_segment([window[0], window[1]], Stroke::new(2.0, color));
    }
}

fn get_input_port_position(node: &rustsim_types::NodeInstance, port_index: usize) -> Position {
    let port_count = node.inputs.len();
    let spacing = 20.0;
    let total_height = (port_count.saturating_sub(1)) as f32 * spacing;
    let start_y = node.position.y + NODE_BASE_HEIGHT / 2.0 - total_height / 2.0;

    Position::new(node.position.x, start_y + port_index as f32 * spacing)
}

fn get_output_port_position(node: &rustsim_types::NodeInstance, port_index: usize) -> Position {
    let port_count = node.outputs.len();
    let spacing = 20.0;
    let total_height = (port_count.saturating_sub(1)) as f32 * spacing;
    let start_y = node.position.y + NODE_BASE_HEIGHT / 2.0 - total_height / 2.0;

    Position::new(
        node.position.x + NODE_BASE_WIDTH,
        start_y + port_index as f32 * spacing,
    )
}
