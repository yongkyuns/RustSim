//! Node graph canvas UI component.

use egui::{Color32, Pos2, Rect, Sense, Stroke, Ui, Vec2};

use crate::icons;
use crate::state::AppState;
use crate::ui::node_preview;
use rustsim_types::{snap_to_grid, NodeShape, Position, NODE_BASE_HEIGHT, NODE_BASE_WIDTH};

/// Render the node graph canvas
pub fn render_node_graph(ui: &mut Ui, state: &mut AppState) {
    let (response, painter) = ui.allocate_painter(ui.available_size(), Sense::click_and_drag());

    let canvas_rect = response.rect;
    let canvas_center = canvas_rect.center();

    // Handle pan
    if response.dragged_by(egui::PointerButton::Middle)
        || (response.dragged_by(egui::PointerButton::Primary) && ui.input(|i| i.modifiers.shift))
    {
        let delta = response.drag_delta();
        let zoom = state.zoom();
        state.pan_mut().x += delta.x / zoom;
        state.pan_mut().y += delta.y / zoom;
    }

    // Handle zoom with scroll
    if response.hovered() {
        let scroll_delta = ui.input(|i| i.raw_scroll_delta.y);
        if scroll_delta != 0.0 {
            let zoom_factor = 1.0 + scroll_delta * 0.001;
            let new_zoom = (state.zoom() * zoom_factor).clamp(0.1, 4.0);
            *state.zoom_mut() = new_zoom;
        }
    }

    // Cache values for closures
    let pan = *state.pan();
    let zoom = state.zoom();

    // Transform functions
    let to_screen = |pos: Position| -> Pos2 {
        Pos2::new(
            canvas_center.x + (pos.x + pan.x) * zoom,
            canvas_center.y + (pos.y + pan.y) * zoom,
        )
    };

    let _to_canvas = |pos: Pos2| -> Position {
        Position::new(
            (pos.x - canvas_center.x) / zoom - pan.x,
            (pos.y - canvas_center.y) / zoom - pan.y,
        )
    };

    // Draw grid (if enabled)
    if state.show_canvas_grid() {
        draw_grid(&painter, canvas_rect, pan, zoom);
    }

    // Draw connections
    for conn in &state.graph().connections {
        let source_node = state.graph().get_node(&conn.source_node_id);
        let target_node = state.graph().get_node(&conn.target_node_id);

        if let (Some(source), Some(target)) = (source_node, target_node) {
            let source_pos = get_output_port_position(source, conn.source_port_index);
            let target_pos = get_input_port_position(target, conn.target_port_index);

            let source_screen = to_screen(source_pos);
            let target_screen = to_screen(target_pos);

            let is_selected = state.selected_connections().contains(&conn.id);
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
    if let Some(pending) = state.pending_connection() {
        if let Some(source_node) = state.graph().get_node(&pending.source_node) {
            let source_pos = get_output_port_position(source_node, pending.source_port);
            let source_screen = to_screen(source_pos);
            let target_screen = to_screen(pending.current_pos);

            draw_connection(
                &painter,
                source_screen,
                target_screen,
                Color32::from_rgb(100, 200, 255),
            );
        }
    }

    // Draw nodes and collect port hit areas
    let mut clicked_node: Option<String> = None;
    let mut drag_delta = Vec2::ZERO;
    let mut port_interactions: Vec<PortInteraction> = Vec::new();
    let mut port_actions: Vec<PortAction> = Vec::new();

    // Collect node data first to avoid borrow conflicts when mutating state
    struct NodeRenderData {
        id: String,
        position: Position,
        name: String,
        block_type: String,
        rotation: u8,
        input_count: usize,
        output_count: usize,
        shape: NodeShape,
        is_selected: bool,
    }

    let nodes_data: Vec<NodeRenderData> = state.graph().nodes.iter().map(|(id, node)| {
        let shape = state
            .get_block_type(&node.block_type)
            .map(|b| b.shape)
            .unwrap_or(NodeShape::Rect);
        let is_selected = state.selected_nodes().contains(id);

        NodeRenderData {
            id: id.clone(),
            position: node.position,
            name: node.name.clone(),
            block_type: node.block_type.clone(),
            rotation: node.rotation,
            input_count: node.inputs.len(),
            output_count: node.outputs.len(),
            shape,
            is_selected,
        }
    }).collect();


    // Now iterate over collected data without holding state borrow
    for node_data in &nodes_data {
        let node_pos = to_screen(node_data.position);
        let node_size = Vec2::new(NODE_BASE_WIDTH * zoom, NODE_BASE_HEIGHT * zoom);
        let node_rect = Rect::from_min_size(node_pos, node_size);

        // Draw node with icon
        draw_node(
            ui,
            &painter,
            node_rect,
            &node_data.name,
            &node_data.block_type,
            node_data.shape,
            node_data.is_selected,
            node_data.rotation,
        );

        // Draw and handle ports
        let port_radius = 6.0 * zoom;

        // Create a temporary NodeInstance-like struct for port position calculation
        let temp_node = rustsim_types::NodeInstance {
            id: node_data.id.clone(),
            name: node_data.name.clone(),
            block_type: node_data.block_type.clone(),
            position: node_data.position,
            rotation: node_data.rotation,
            inputs: (0..node_data.input_count).map(|i| rustsim_types::PortInstance::new_input(&node_data.id, i, &format!("in_{}", i))).collect(),
            outputs: (0..node_data.output_count).map(|i| rustsim_types::PortInstance::new_output(&node_data.id, i, &format!("out_{}", i))).collect(),
            params: Default::default(),
            pinned_params: Vec::new(),
            color: None,
            graph: None,
        };

        // Input ports
        for i in 0..node_data.input_count {
            let port_canvas_pos = get_input_port_position(&temp_node, i);
            let port_pos = to_screen(port_canvas_pos);
            let port_id = egui::Id::new((node_data.id.as_str(), "input", i));
            let port_rect = Rect::from_center_size(port_pos, Vec2::splat(port_radius * 2.0));
            let port_response = ui.interact(port_rect, port_id, Sense::click_and_drag());

            let is_hovered = port_response.hovered();
            draw_port(&painter, port_pos, true, is_hovered);

            port_interactions.push(PortInteraction {
                node_id: node_data.id.clone(),
                port_index: i,
                is_output: false,
                screen_pos: port_pos,
                canvas_pos: port_canvas_pos,
                response: port_response,
            });
        }

        // Output ports
        for i in 0..node_data.output_count {
            let port_canvas_pos = get_output_port_position(&temp_node, i);
            let port_pos = to_screen(port_canvas_pos);
            let port_id = egui::Id::new((node_data.id.as_str(), "output", i));
            let port_rect = Rect::from_center_size(port_pos, Vec2::splat(port_radius * 2.0));
            let port_response = ui.interact(port_rect, port_id, Sense::click_and_drag());

            let is_hovered = port_response.hovered();
            draw_port(&painter, port_pos, false, is_hovered);

            port_interactions.push(PortInteraction {
                node_id: node_data.id.clone(),
                port_index: i,
                is_output: true,
                screen_pos: port_pos,
                canvas_pos: port_canvas_pos,
                response: port_response,
            });
        }

        // Draw plot preview for Scope and Spectrum nodes
        if node_data.block_type == "Scope" || node_data.block_type == "Spectrum" {
            draw_node_preview(&painter, node_rect, &node_data.id, &node_data.block_type, state);
        }

        // Handle node interaction
        let node_response = ui.interact(node_rect, egui::Id::new(&node_data.id), Sense::click_and_drag());

        // Start dragging immediately on mouse down - select and prepare to drag
        if node_response.drag_started() {
            *state.dragging_node_mut() = Some(node_data.id.clone());
            // Select this node for dragging (unless Ctrl is held for multi-select)
            let ctrl = ui.input(|i| i.modifiers.ctrl || i.modifiers.command);
            if !ctrl && !state.selected_nodes().contains(&node_data.id) {
                state.selected_nodes_mut().clear();
            }
            state.selected_nodes_mut().insert(node_data.id.clone());
        }

        // Apply drag delta while dragging
        if node_response.dragged() && state.dragging_node().as_ref() == Some(&node_data.id) {
            drag_delta = node_response.drag_delta();
        }

        // On drag stop, clear selection if it was just a click (no movement)
        if node_response.drag_stopped() {
            *state.dragging_node_mut() = None;
        }

        // Simple click without drag - handle selection toggle
        if node_response.clicked() {
            clicked_node = Some(node_data.id.clone());
        }

        // Double-click to open parameter editor
        if node_response.double_clicked() {
            state.editing_node = Some(node_data.id.clone());
        }
    }

    // Draw port controls for selected nodes (after all nodes are drawn)
    let selected_nodes: Vec<String> = state.selected_nodes().iter().cloned().collect();
    for node_id in &selected_nodes {
        if let Some(node) = state.graph().get_node(node_id) {
            let node_pos = to_screen(node.position);
            let node_size = Vec2::new(NODE_BASE_WIDTH * zoom, NODE_BASE_HEIGHT * zoom);
            let node_rect = Rect::from_min_size(node_pos, node_size);

            // Check port capabilities
            let can_add_input = state.can_add_input(node_id);
            let can_remove_input = state.can_remove_input(node_id);
            let can_add_output = state.can_add_output(node_id);
            let can_remove_output = state.can_remove_output(node_id);
            let sync_ports = state.has_sync_ports(node_id);

            let btn_size = Vec2::new(16.0 * zoom, 16.0 * zoom);

            // Input port controls (left side)
            if can_add_input || can_remove_input {
                let input_ctrl_center =
                    Pos2::new(node_rect.left() - 12.0 * zoom, node_rect.center().y);

                // Plus button
                let plus_rect = Rect::from_center_size(
                    Pos2::new(input_ctrl_center.x, input_ctrl_center.y - 10.0 * zoom),
                    btn_size,
                );
                let plus_response = ui.interact(
                    plus_rect,
                    egui::Id::new((node_id.as_str(), "add_input")),
                    Sense::click(),
                );

                if can_add_input {
                    draw_port_button(&painter, plus_rect, "+", plus_response.hovered());
                    if plus_response.clicked() {
                        port_actions.push(PortAction::AddInput(node_id.clone()));
                    }
                }

                // Minus button
                let minus_rect = Rect::from_center_size(
                    Pos2::new(input_ctrl_center.x, input_ctrl_center.y + 10.0 * zoom),
                    btn_size,
                );
                let minus_response = ui.interact(
                    minus_rect,
                    egui::Id::new((node_id.as_str(), "remove_input")),
                    Sense::click(),
                );

                if can_remove_input {
                    draw_port_button(&painter, minus_rect, "-", minus_response.hovered());
                    if minus_response.clicked() {
                        port_actions.push(PortAction::RemoveInput(node_id.clone()));
                    }
                }
            }

            // Output port controls (right side) - hidden for sync_ports blocks
            if !sync_ports && (can_add_output || can_remove_output) {
                let output_ctrl_center =
                    Pos2::new(node_rect.right() + 12.0 * zoom, node_rect.center().y);

                // Plus button
                let plus_rect = Rect::from_center_size(
                    Pos2::new(
                        output_ctrl_center.x,
                        output_ctrl_center.y - 10.0 * zoom,
                    ),
                    btn_size,
                );
                let plus_response = ui.interact(
                    plus_rect,
                    egui::Id::new((node_id.as_str(), "add_output")),
                    Sense::click(),
                );

                if can_add_output {
                    draw_port_button(&painter, plus_rect, "+", plus_response.hovered());
                    if plus_response.clicked() {
                        port_actions.push(PortAction::AddOutput(node_id.clone()));
                    }
                }

                // Minus button
                let minus_rect = Rect::from_center_size(
                    Pos2::new(
                        output_ctrl_center.x,
                        output_ctrl_center.y + 10.0 * zoom,
                    ),
                    btn_size,
                );
                let minus_response = ui.interact(
                    minus_rect,
                    egui::Id::new((node_id.as_str(), "remove_output")),
                    Sense::click(),
                );

                if can_remove_output {
                    draw_port_button(&painter, minus_rect, "-", minus_response.hovered());
                    if minus_response.clicked() {
                        port_actions.push(PortAction::RemoveOutput(node_id.clone()));
                    }
                }
            }
        }
    }

    // Apply port actions
    for action in port_actions {
        match action {
            PortAction::AddInput(id) => {
                state.add_input_port(&id);
            }
            PortAction::RemoveInput(id) => {
                state.remove_input_port(&id);
            }
            PortAction::AddOutput(id) => {
                state.add_output_port(&id);
            }
            PortAction::RemoveOutput(id) => {
                state.remove_output_port(&id);
            }
        }
    }

    // Handle port interactions for connection creation
    handle_port_connections(ui, state, &port_interactions, canvas_center);

    // Handle node selection
    if let Some(node_id) = clicked_node {
        let ctrl = ui.input(|i| i.modifiers.ctrl || i.modifiers.command);
        if ctrl {
            // Toggle selection
            if state.selected_nodes().contains(&node_id) {
                state.selected_nodes_mut().remove(&node_id);
            } else {
                state.selected_nodes_mut().insert(node_id);
            }
        } else {
            // Single selection
            state.selected_nodes_mut().clear();
            state.selected_nodes_mut().insert(node_id);
        }
    } else if response.clicked() {
        // Clicked on empty space - clear selection
        state.clear_selection();
    }

    // Update dragged node positions
    if drag_delta != Vec2::ZERO {
        let nodes_to_update: Vec<String> = state.selected_nodes().iter().cloned().collect();
        let should_snap = state.dragging_node().is_none();
        for node_id in &nodes_to_update {
            if let Some(node) = state.graph_mut().get_node_mut(node_id) {
                node.position.x += drag_delta.x / zoom;
                node.position.y += drag_delta.y / zoom;
                // Snap to grid on release
                if should_snap {
                    node.position.x = snap_to_grid(node.position.x);
                    node.position.y = snap_to_grid(node.position.y);
                }
            }
        }
    }
}

fn draw_grid(painter: &egui::Painter, rect: Rect, pan: Position, zoom: f32) {
    let grid_size = 20.0 * zoom;
    let grid_color = Color32::from_rgba_unmultiplied(128, 128, 128, 5);

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
    ui: &mut Ui,
    painter: &egui::Painter,
    rect: Rect,
    name: &str,
    block_type: &str,
    shape: NodeShape,
    selected: bool,
    _rotation: u8,
) {
    let bg_color = Color32::from_rgb(45, 45, 48);
    let border_color = if selected {
        Color32::from_rgb(100, 200, 255)
    } else {
        Color32::from_rgb(80, 80, 80)
    };

    // Handle triangle shape separately
    if matches!(shape, NodeShape::Triangle) {
        let corner_radius = 4.0;

        // Selection glow
        if selected {
            let glow_rect = rect.expand(3.0);
            let glow_path = create_rounded_triangle_path(glow_rect, corner_radius + 1.0);
            painter.add(egui::Shape::Path(egui::epaint::PathShape {
                points: glow_path,
                closed: true,
                fill: Color32::from_rgba_unmultiplied(100, 200, 255, 50),
                stroke: Stroke::NONE.into(),
            }));
        }

        // Triangle background with rounded corners
        let triangle_path = create_rounded_triangle_path(rect, corner_radius);
        painter.add(egui::Shape::Path(egui::epaint::PathShape {
            points: triangle_path,
            closed: true,
            fill: bg_color,
            stroke: Stroke::new(1.0, border_color).into(),
        }));

        // For triangle, draw icon and text in the left 2/3 of the shape
        draw_node_content_triangle(ui, painter, rect, name, block_type);
        return;
    }

    // Regular shapes with rounding
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
        NodeShape::Triangle => egui::Rounding::ZERO, // Handled above
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

    // Calculate layout: icon on top, text below, vertically centered
    let icon_size = NODE_ICON_SIZE.min(rect.width() - NODE_PADDING * 2.0).max(NODE_MIN_ICON_SIZE);

    // Get the icon URI for this block type
    let icon_uri = icons::get_icon_uri(block_type);

    // Total content height (icon + gap + text)
    let total_content_height = icon_size + NODE_GAP + NODE_TEXT_HEIGHT;

    // Draw icon and name vertically stacked and centered
    if rect.height() > total_content_height + NODE_PADDING * 2.0 {
        // Tall enough for icon + name below, center the whole block
        let content_top = rect.center().y - total_content_height / 2.0;
        let icon_center_y = content_top + icon_size / 2.0;
        let icon_rect = Rect::from_center_size(
            Pos2::new(rect.center().x, icon_center_y),
            Vec2::splat(icon_size),
        );

        // Draw icon using egui Image
        ui.allocate_new_ui(egui::UiBuilder::new().max_rect(icon_rect), |ui| {
            ui.add(
                egui::Image::new(&icon_uri)
                    .fit_to_exact_size(Vec2::splat(icon_size))
                    .tint(Color32::from_rgb(200, 200, 200)),
            );
        });

        // Draw name below icon
        let text_y = icon_rect.bottom() + NODE_GAP + NODE_TEXT_HEIGHT / 2.0;
        painter.text(
            Pos2::new(rect.center().x, text_y),
            egui::Align2::CENTER_CENTER,
            name,
            egui::FontId::proportional(NODE_FONT_SIZE),
            Color32::WHITE,
        );
    } else if rect.height() > icon_size + NODE_PADDING * 2.0 {
        // Just icon, no text (too short for both)
        let icon_rect = Rect::from_center_size(rect.center(), Vec2::splat(icon_size));
        ui.allocate_new_ui(egui::UiBuilder::new().max_rect(icon_rect), |ui| {
            ui.add(
                egui::Image::new(&icon_uri)
                    .fit_to_exact_size(Vec2::splat(icon_size))
                    .tint(Color32::from_rgb(200, 200, 200)),
            );
        });
    } else {
        // Too small, just draw name
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            name,
            egui::FontId::proportional(NODE_FONT_SIZE),
            Color32::WHITE,
        );
    }
}

/// Create a rounded triangle path (right-pointing triangle with rounded corners)
fn create_rounded_triangle_path(rect: Rect, radius: f32) -> Vec<Pos2> {
    let mut points = Vec::new();

    // Three vertices of the triangle (clockwise order)
    let top_left = Pos2::new(rect.left(), rect.top());
    let right_tip = Pos2::new(rect.right(), rect.center().y);
    let bottom_left = Pos2::new(rect.left(), rect.bottom());

    // Number of points per corner arc
    let segments = 4;

    // Helper to add rounded corner points
    // prev -> corner -> next defines the path direction
    let add_rounded_corner = |points: &mut Vec<Pos2>, prev: Pos2, corner: Pos2, next: Pos2, r: f32| {
        // Direction from corner to previous point
        let to_prev = Vec2::new(prev.x - corner.x, prev.y - corner.y).normalized();
        // Direction from corner to next point
        let to_next = Vec2::new(next.x - corner.x, next.y - corner.y).normalized();

        // Points where the arc starts and ends (offset from corner)
        let arc_start = Pos2::new(corner.x + to_prev.x * r, corner.y + to_prev.y * r);
        let arc_end = Pos2::new(corner.x + to_next.x * r, corner.y + to_next.y * r);

        // Generate arc points using quadratic bezier approximation
        for i in 0..=segments {
            let t = i as f32 / segments as f32;
            // Quadratic bezier: (1-t)^2 * start + 2*(1-t)*t * corner + t^2 * end
            let mt = 1.0 - t;
            let x = mt * mt * arc_start.x + 2.0 * mt * t * corner.x + t * t * arc_end.x;
            let y = mt * mt * arc_start.y + 2.0 * mt * t * corner.y + t * t * arc_end.y;
            points.push(Pos2::new(x, y));
        }
    };

    // Traverse clockwise: top-left → right-tip → bottom-left

    // Top-left corner: coming from bottom-left (via left edge), going to right-tip
    add_rounded_corner(&mut points, bottom_left, top_left, right_tip, radius);

    // Right tip: coming from top-left, going to bottom-left
    add_rounded_corner(&mut points, top_left, right_tip, bottom_left, radius);

    // Bottom-left corner: coming from right-tip, going to top-left (via left edge)
    add_rounded_corner(&mut points, right_tip, bottom_left, top_left, radius);

    points
}

// Shared constants for node content layout
const NODE_ICON_SIZE: f32 = 24.0;
const NODE_TEXT_HEIGHT: f32 = 14.0;
const NODE_PADDING: f32 = 4.0;
const NODE_GAP: f32 = 2.0;
const NODE_FONT_SIZE: f32 = 10.0;
const NODE_MIN_ICON_SIZE: f32 = 12.0;

/// Draw content (icon + label) inside a triangular node
fn draw_node_content_triangle(
    ui: &mut Ui,
    painter: &egui::Painter,
    rect: Rect,
    name: &str,
    block_type: &str,
) {
    // For a right-pointing triangle, the usable area is roughly the left 60%
    // The centroid of the triangle is at (left + width/3, center_y)
    let content_center_x = rect.left() + rect.width() * 0.35;
    let content_center_y = rect.center().y;

    // Use same sizing as regular nodes for consistency
    let icon_size = NODE_ICON_SIZE.min(rect.width() * 0.5 - NODE_PADDING * 2.0).max(NODE_MIN_ICON_SIZE);

    let icon_uri = icons::get_icon_uri(block_type);

    // Total content height
    let total_content_height = icon_size + NODE_GAP + NODE_TEXT_HEIGHT;

    if rect.height() > total_content_height + NODE_PADDING * 2.0 {
        // Icon + text
        let content_top = content_center_y - total_content_height / 2.0;
        let icon_center_y = content_top + icon_size / 2.0;
        let icon_rect = Rect::from_center_size(
            Pos2::new(content_center_x, icon_center_y),
            Vec2::splat(icon_size),
        );

        ui.allocate_new_ui(egui::UiBuilder::new().max_rect(icon_rect), |ui| {
            ui.add(
                egui::Image::new(&icon_uri)
                    .fit_to_exact_size(Vec2::splat(icon_size))
                    .tint(Color32::from_rgb(200, 200, 200)),
            );
        });

        let text_y = icon_rect.bottom() + NODE_GAP + NODE_TEXT_HEIGHT / 2.0;
        painter.text(
            Pos2::new(content_center_x, text_y),
            egui::Align2::CENTER_CENTER,
            name,
            egui::FontId::proportional(NODE_FONT_SIZE),
            Color32::WHITE,
        );
    } else if rect.height() > icon_size + NODE_PADDING * 2.0 {
        // Just icon
        let icon_rect = Rect::from_center_size(
            Pos2::new(content_center_x, content_center_y),
            Vec2::splat(icon_size),
        );
        ui.allocate_new_ui(egui::UiBuilder::new().max_rect(icon_rect), |ui| {
            ui.add(
                egui::Image::new(&icon_uri)
                    .fit_to_exact_size(Vec2::splat(icon_size))
                    .tint(Color32::from_rgb(200, 200, 200)),
            );
        });
    } else {
        // Just name
        painter.text(
            Pos2::new(content_center_x, content_center_y),
            egui::Align2::CENTER_CENTER,
            name,
            egui::FontId::proportional(NODE_FONT_SIZE),
            Color32::WHITE,
        );
    }
}

fn draw_port(painter: &egui::Painter, pos: Pos2, _is_input: bool, is_hovered: bool) {
    let radius = 4.0;
    let color = if is_hovered {
        Color32::from_rgb(100, 200, 255)
    } else {
        Color32::from_rgb(150, 150, 150)
    };

    // Hollow circle for ports
    painter.circle_stroke(pos, radius, Stroke::new(2.0, color));
}

fn draw_port_button(painter: &egui::Painter, rect: Rect, label: &str, is_hovered: bool) {
    let bg_color = if is_hovered {
        Color32::from_rgb(80, 80, 85)
    } else {
        Color32::from_rgb(55, 55, 60)
    };
    let text_color = if is_hovered {
        Color32::from_rgb(100, 200, 255)
    } else {
        Color32::from_rgb(180, 180, 180)
    };

    painter.rect(
        rect,
        egui::Rounding::same(3.0),
        bg_color,
        Stroke::new(1.0, Color32::from_rgb(80, 80, 80)),
    );

    painter.text(
        rect.center(),
        egui::Align2::CENTER_CENTER,
        label,
        egui::FontId::proportional(12.0),
        text_color,
    );
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

/// Port action to apply after rendering
enum PortAction {
    AddInput(String),
    RemoveInput(String),
    AddOutput(String),
    RemoveOutput(String),
}

/// Port interaction data
struct PortInteraction {
    node_id: String,
    port_index: usize,
    is_output: bool,
    #[allow(dead_code)]
    screen_pos: Pos2,
    canvas_pos: Position,
    response: egui::Response,
}

/// Handle port interactions for connection creation
fn handle_port_connections(
    ui: &mut Ui,
    state: &mut AppState,
    ports: &[PortInteraction],
    canvas_center: Pos2,
) {
    let mouse_pos = ui.input(|i| i.pointer.hover_pos());

    // Check if any output port started a drag (begin connection)
    for port in ports.iter().filter(|p| p.is_output) {
        if port.response.drag_started() {
            *state.pending_connection_mut() = Some(crate::state::PendingConnection {
                source_node: port.node_id.clone(),
                source_port: port.port_index,
                current_pos: port.canvas_pos,
            });
        }
    }

    // Update pending connection position
    let zoom = state.zoom();
    let pan = *state.pan();
    if let Some(pending) = state.pending_connection_mut() {
        if let Some(mouse) = mouse_pos {
            // Convert screen position to canvas position
            let canvas_x = (mouse.x - canvas_center.x) / zoom - pan.x;
            let canvas_y = (mouse.y - canvas_center.y) / zoom - pan.y;
            pending.current_pos = Position::new(canvas_x, canvas_y);
        }
    }

    // Check if pending connection was released over an input port
    if state.pending_connection().is_some() {
        let released = ui.input(|i| i.pointer.any_released());

        if released {
            // Find if we're hovering over an input port
            let mut target_port: Option<(String, usize)> = None;

            for port in ports.iter().filter(|p| !p.is_output) {
                if port.response.hovered() {
                    target_port = Some((port.node_id.clone(), port.port_index));
                    break;
                }
            }

            if let Some((target_node, target_port_idx)) = target_port {
                if let Some(pending) = state.pending_connection_mut().take() {
                    // Create the connection (only if not connecting to same node)
                    if pending.source_node != target_node {
                        state.add_connection(
                            &pending.source_node,
                            pending.source_port,
                            &target_node,
                            target_port_idx,
                        );
                    }
                }
            } else {
                // Released without target - cancel connection
                *state.pending_connection_mut() = None;
            }
        }
    }
}

/// Draw plot preview for Scope or Spectrum nodes
fn draw_node_preview(
    painter: &egui::Painter,
    node_rect: Rect,
    node_id: &str,
    block_type: &str,
    state: &AppState,
) {
    // Define preview area below the node
    let preview_width = node_rect.width() * 0.9;
    let preview_height = 40.0;
    let preview_margin = 5.0;

    let preview_rect = Rect::from_min_size(
        Pos2::new(
            node_rect.min.x + (node_rect.width() - preview_width) / 2.0,
            node_rect.max.y + preview_margin,
        ),
        Vec2::new(preview_width, preview_height),
    );

    // Color palette for multiple channels
    const CHANNEL_COLORS: [Color32; 6] = [
        Color32::from_rgb(100, 200, 255), // Cyan
        Color32::from_rgb(255, 150, 100), // Orange
        Color32::from_rgb(150, 255, 150), // Green
        Color32::from_rgb(255, 150, 255), // Pink
        Color32::from_rgb(255, 255, 100), // Yellow
        Color32::from_rgb(150, 150, 255), // Purple
    ];

    match block_type {
        "Scope" => {
            // Get scope data for this node
            if let Some(data) = state.scope_data().get(node_id) {
                if !data.is_empty() {
                    // Determine number of channels from first data point
                    let num_channels = data.first().map(|(_, v)| v.len()).unwrap_or(0);

                    // Draw each channel with a different color
                    for channel_idx in 0..num_channels {
                        let channel_data: Vec<(f64, f64)> = data
                            .iter()
                            .filter_map(|(t, values)| {
                                values.get(channel_idx).map(|&v| (*t, v))
                            })
                            .collect();

                        let color = CHANNEL_COLORS[channel_idx % CHANNEL_COLORS.len()];
                        node_preview::draw_plot_preview(painter, preview_rect, &channel_data, color);
                    }
                }
            }
        }
        "Spectrum" => {
            // Get spectrum data for this node
            if let Some(data) = state.spectrum_data().get(node_id) {
                if !data.is_empty() {
                    // Use first channel for preview
                    let channel_data: Vec<(f64, f64)> = data
                        .iter()
                        .map(|(f, values)| (*f, values.first().copied().unwrap_or(0.0)))
                        .collect();

                    node_preview::draw_plot_preview(
                        painter,
                        preview_rect,
                        &channel_data,
                        Color32::from_rgb(255, 200, 100),
                    );
                }
            }
        }
        _ => {}
    }
}
