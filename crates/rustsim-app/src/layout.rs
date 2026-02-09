//! Auto-layout for simulation graphs using layout-rs (Sugiyama-style layered layout).

use crate::routing::RotationOptimizer;
use layout::backends::svg::SVGWriter;
use layout::core::base::Orientation;
use layout::core::geometry::Point;
use layout::core::style::StyleAttr;
use layout::std_shapes::shapes::{Arrow, Element, ShapeKind};
use layout::topo::layout::VisualGraph;
use rustsim_types::{Position, SimulationGraph};
use std::collections::HashMap;

/// Node size for layout calculation
const NODE_WIDTH: f64 = 120.0;
const NODE_HEIGHT: f64 = 60.0;

/// Spacing multipliers
const HORIZONTAL_SPACING: f32 = 1.5;
const VERTICAL_SPACING: f32 = 1.3;

/// Compute auto-layout positions for a simulation graph.
///
/// Uses the Sugiyama algorithm (layered/hierarchical layout) which is ideal
/// for directed dataflow graphs like block diagrams.
///
/// Returns a map of node_id -> new Position
pub fn compute_layout(graph: &SimulationGraph) -> HashMap<String, Position> {
    let mut positions = HashMap::new();

    if graph.nodes.is_empty() {
        return positions;
    }

    // Create a VisualGraph for layout-rs
    let mut vg = VisualGraph::new(Orientation::LeftToRight);

    // Map from our node IDs to layout-rs handles
    let mut id_to_handle = HashMap::new();
    let mut handle_to_id = HashMap::new();

    // Add nodes to the visual graph
    for (node_id, _node) in &graph.nodes {
        let shape = ShapeKind::new_box(node_id);
        let style = StyleAttr::simple();
        let size = Point::new(NODE_WIDTH, NODE_HEIGHT);
        let element = Element::create(shape, style, Orientation::LeftToRight, size);

        let handle = vg.add_node(element);
        id_to_handle.insert(node_id.clone(), handle);
        handle_to_id.insert(handle.get_index(), node_id.clone());
    }

    // Add edges to the visual graph
    for conn in &graph.connections {
        if let (Some(&from_handle), Some(&to_handle)) = (
            id_to_handle.get(&conn.source_node_id),
            id_to_handle.get(&conn.target_node_id),
        ) {
            let arrow = Arrow::simple("");
            vg.add_edge(arrow, from_handle, to_handle);
        }
    }

    // Run the layout algorithm
    // We need a render backend even if we don't use the output
    let mut svg = SVGWriter::new();
    vg.do_it(false, false, false, &mut svg);

    // Extract positions from the laid-out graph
    for node_handle in vg.iter_nodes() {
        let pos = vg.pos(node_handle);
        if let Some(node_id) = handle_to_id.get(&node_handle.get_index()) {
            // Convert layout-rs Position to our Position
            // pos.middle() returns the center point of the node
            let middle = pos.middle();
            // Apply spacing multipliers and offset for better visual appearance
            let x = middle.x as f32 * HORIZONTAL_SPACING + 50.0;
            let y = middle.y as f32 * VERTICAL_SPACING + 50.0;
            positions.insert(node_id.clone(), Position::new(x, y));
        }
    }

    positions
}

/// Apply computed layout positions to a simulation graph
pub fn apply_layout(graph: &mut SimulationGraph, positions: &HashMap<String, Position>) {
    for (node_id, pos) in positions {
        if let Some(node) = graph.nodes.get_mut(node_id) {
            node.position = *pos;
        }
    }
}

/// Perform auto-layout on a simulation graph (compute and apply in one step)
/// Also optimizes node rotations to minimize wire crossings.
pub fn auto_layout(graph: &mut SimulationGraph) {
    // Step 1: Compute and apply positions using Sugiyama algorithm
    let positions = compute_layout(graph);
    apply_layout(graph, &positions);

    // Step 2: Optimize rotations to minimize wire crossings
    optimize_rotations(graph);
}

/// Optimize node rotations to minimize wire crossings
pub fn optimize_rotations(graph: &mut SimulationGraph) {
    if graph.nodes.is_empty() || graph.connections.is_empty() {
        return;
    }

    let optimizer = RotationOptimizer::new();
    let optimal_rotations = optimizer.optimize_rotations(&graph.nodes, &graph.connections);
    RotationOptimizer::apply_rotations(&mut graph.nodes, &optimal_rotations);
}

/// Perform position layout only (without rotation optimization)
pub fn auto_layout_positions_only(graph: &mut SimulationGraph) {
    let positions = compute_layout(graph);
    apply_layout(graph, &positions);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::examples;
    use crate::routing::OrthogonalRouter;

    #[test]
    fn test_auto_layout_pid() {
        let mut graph = examples::pid_controller();
        let original_positions: HashMap<_, _> = graph
            .nodes
            .iter()
            .map(|(id, n)| (id.clone(), n.position))
            .collect();

        auto_layout(&mut graph);

        // Check that positions changed (layout was applied)
        let mut any_changed = false;
        for (id, node) in &graph.nodes {
            if let Some(orig) = original_positions.get(id) {
                if node.position.x != orig.x || node.position.y != orig.y {
                    any_changed = true;
                    break;
                }
            }
        }
        assert!(any_changed, "Auto-layout should change at least one position");
    }

    #[test]
    fn test_auto_layout_empty_graph() {
        let mut graph = SimulationGraph::new();
        auto_layout(&mut graph);
        // Should not panic on empty graph
        assert!(graph.nodes.is_empty());
    }

    #[test]
    fn test_rotation_optimization() {
        let mut graph = examples::pid_controller();

        // Get crossings before optimization
        let router = OrthogonalRouter::new();
        let crossings_before = router.count_crossings(&graph.connections, &graph.nodes);

        // Optimize rotations
        optimize_rotations(&mut graph);

        // Crossings should not increase
        let crossings_after = router.count_crossings(&graph.connections, &graph.nodes);
        assert!(
            crossings_after <= crossings_before,
            "Rotation optimization should not increase crossings"
        );
    }

    #[test]
    fn test_auto_layout_includes_rotation() {
        let mut graph = examples::pid_controller();

        // Set all rotations to non-zero to check if optimization changes them
        for node in graph.nodes.values_mut() {
            node.rotation = 2; // 180 degrees
        }

        auto_layout(&mut graph);

        // After auto_layout, rotations should be optimized
        // (may or may not change depending on the layout)
        // Just verify no panic and graph is still valid
        assert!(!graph.nodes.is_empty());
    }
}
