//! Orthogonal router for connections with waypoint support and rotation awareness.

use rustsim_types::{Connection, NodeInstance, Position, NODE_BASE_HEIGHT, NODE_BASE_WIDTH};
use std::collections::HashMap;

use super::astar::{AStarRouter, Direction, Rect};

/// Get the direction a port faces based on node rotation
/// Rotation: 0 = normal, 1 = 90° CW, 2 = 180°, 3 = 270° CW
fn rotate_direction(base_direction: Direction, rotation: u8) -> Direction {
    let rotations = rotation % 4;
    let mut dir = base_direction;
    for _ in 0..rotations {
        dir = match dir {
            Direction::Right => Direction::Down,
            Direction::Down => Direction::Left,
            Direction::Left => Direction::Up,
            Direction::Up => Direction::Right,
        };
    }
    dir
}

/// Get the position offset for a port based on rotation
fn get_port_offset(
    is_output: bool,
    port_index: usize,
    port_count: usize,
    rotation: u8,
) -> (f32, f32) {
    let spacing = 20.0;
    let total_span = (port_count.saturating_sub(1)) as f32 * spacing;
    let port_offset = port_index as f32 * spacing - total_span / 2.0;

    // Base positions: outputs on right, inputs on left
    // After rotation, positions change
    match rotation % 4 {
        0 => {
            // Normal: outputs right, inputs left
            if is_output {
                (NODE_BASE_WIDTH, NODE_BASE_HEIGHT / 2.0 + port_offset)
            } else {
                (0.0, NODE_BASE_HEIGHT / 2.0 + port_offset)
            }
        }
        1 => {
            // 90° CW: outputs bottom, inputs top
            if is_output {
                (NODE_BASE_WIDTH / 2.0 + port_offset, NODE_BASE_HEIGHT)
            } else {
                (NODE_BASE_WIDTH / 2.0 + port_offset, 0.0)
            }
        }
        2 => {
            // 180°: outputs left, inputs right
            if is_output {
                (0.0, NODE_BASE_HEIGHT / 2.0 - port_offset)
            } else {
                (NODE_BASE_WIDTH, NODE_BASE_HEIGHT / 2.0 - port_offset)
            }
        }
        3 => {
            // 270° CW: outputs top, inputs bottom
            if is_output {
                (NODE_BASE_WIDTH / 2.0 - port_offset, 0.0)
            } else {
                (NODE_BASE_WIDTH / 2.0 - port_offset, NODE_BASE_HEIGHT)
            }
        }
        _ => unreachable!(),
    }
}

/// Orthogonal router that uses A* for path segments
pub struct OrthogonalRouter {
    astar: AStarRouter,
}

impl OrthogonalRouter {
    /// Create a new orthogonal router with default grid size
    pub fn new() -> Self {
        Self {
            astar: AStarRouter::default(),
        }
    }

    /// Create a new orthogonal router with custom grid size
    pub fn with_grid_size(grid_size: f32) -> Self {
        Self {
            astar: AStarRouter::new(grid_size),
        }
    }

    /// Route a connection between nodes (rotation-aware)
    pub fn route(
        &self,
        connection: &Connection,
        nodes: &HashMap<String, NodeInstance>,
    ) -> Vec<Position> {
        // Get source and target nodes
        let source_node = match nodes.get(&connection.source_node_id) {
            Some(n) => n,
            None => return Vec::new(),
        };

        let target_node = match nodes.get(&connection.target_node_id) {
            Some(n) => n,
            None => return Vec::new(),
        };

        // Calculate port positions with rotation
        let source_pos =
            self.get_output_port_position(source_node, connection.source_port_index);
        let target_pos = self.get_input_port_position(target_node, connection.target_port_index);

        // Get port directions based on rotation
        let source_dir = rotate_direction(Direction::Right, source_node.rotation);
        let target_dir = rotate_direction(Direction::Left, target_node.rotation);

        // Collect obstacles (all node bounding boxes except source and target)
        let obstacles = self.collect_obstacles(
            nodes,
            &connection.source_node_id,
            &connection.target_node_id,
        );

        // Get user waypoints
        let user_waypoints: Vec<Position> = connection
            .waypoints
            .iter()
            .filter(|wp| wp.is_user_waypoint)
            .map(|wp| wp.position)
            .collect();

        // Route with or without waypoints
        if user_waypoints.is_empty() {
            // Direct routing from source to target
            self.astar.route_connection(
                source_pos,
                source_dir,
                target_pos,
                target_dir,
                &obstacles,
            )
        } else {
            // Route through waypoints
            self.route_with_waypoints(source_pos, target_pos, &user_waypoints, &obstacles)
        }
    }

    /// Route from start to end through a series of waypoints
    pub fn route_with_waypoints(
        &self,
        start: Position,
        end: Position,
        waypoints: &[Position],
        obstacles: &[Rect],
    ) -> Vec<Position> {
        if waypoints.is_empty() {
            return self.astar.find_path(start, end, obstacles);
        }

        let mut full_path = Vec::new();

        // Route from start to first waypoint
        let mut segment = self.astar.find_path(start, waypoints[0], obstacles);
        full_path.extend_from_slice(&segment);

        // Route between waypoints
        for i in 0..waypoints.len() - 1 {
            segment = self.astar.find_path(waypoints[i], waypoints[i + 1], obstacles);
            // Skip the first point to avoid duplicates
            if !segment.is_empty() {
                full_path.extend_from_slice(&segment[1..]);
            }
        }

        // Route from last waypoint to end
        segment = self.astar.find_path(
            waypoints[waypoints.len() - 1],
            end,
            obstacles,
        );
        if !segment.is_empty() {
            full_path.extend_from_slice(&segment[1..]);
        }

        full_path
    }

    /// Get the position of an output port on a node (rotation-aware)
    fn get_output_port_position(&self, node: &NodeInstance, port_index: usize) -> Position {
        let (dx, dy) = get_port_offset(true, port_index, node.outputs.len(), node.rotation);
        Position::new(node.position.x + dx, node.position.y + dy)
    }

    /// Get the position of an input port on a node (rotation-aware)
    fn get_input_port_position(&self, node: &NodeInstance, port_index: usize) -> Position {
        let (dx, dy) = get_port_offset(false, port_index, node.inputs.len(), node.rotation);
        Position::new(node.position.x + dx, node.position.y + dy)
    }

    /// Collect obstacle rectangles from all nodes except source and target
    fn collect_obstacles(
        &self,
        nodes: &HashMap<String, NodeInstance>,
        source_id: &str,
        target_id: &str,
    ) -> Vec<Rect> {
        nodes
            .iter()
            .filter(|(id, _)| *id != source_id && *id != target_id)
            .map(|(_, node)| {
                // Add a small margin around nodes
                let margin = 5.0;
                Rect::new(
                    node.position.x - margin,
                    node.position.y - margin,
                    node.position.x + NODE_BASE_WIDTH + margin,
                    node.position.y + NODE_BASE_HEIGHT + margin,
                )
            })
            .collect()
    }

    /// Calculate total wire length for all connections
    pub fn calculate_total_wire_length(
        &self,
        connections: &[Connection],
        nodes: &HashMap<String, NodeInstance>,
    ) -> f32 {
        let mut total = 0.0;
        for conn in connections {
            let path = self.route(conn, nodes);
            for i in 0..path.len().saturating_sub(1) {
                let dx = path[i + 1].x - path[i].x;
                let dy = path[i + 1].y - path[i].y;
                total += (dx * dx + dy * dy).sqrt();
            }
        }
        total
    }

    /// Count wire crossings between all connections
    pub fn count_crossings(
        &self,
        connections: &[Connection],
        nodes: &HashMap<String, NodeInstance>,
    ) -> usize {
        let paths: Vec<Vec<Position>> = connections
            .iter()
            .map(|conn| self.route(conn, nodes))
            .collect();

        let mut crossings = 0;

        // Check each pair of paths for crossings
        for i in 0..paths.len() {
            for j in (i + 1)..paths.len() {
                crossings += self.count_path_crossings(&paths[i], &paths[j]);
            }
        }

        crossings
    }

    /// Count crossings between two paths
    fn count_path_crossings(&self, path1: &[Position], path2: &[Position]) -> usize {
        let mut crossings = 0;

        for i in 0..path1.len().saturating_sub(1) {
            let a1 = &path1[i];
            let a2 = &path1[i + 1];

            for j in 0..path2.len().saturating_sub(1) {
                let b1 = &path2[j];
                let b2 = &path2[j + 1];

                if self.segments_cross(a1, a2, b1, b2) {
                    crossings += 1;
                }
            }
        }

        crossings
    }

    /// Check if two line segments cross (orthogonal segments only)
    fn segments_cross(&self, a1: &Position, a2: &Position, b1: &Position, b2: &Position) -> bool {
        let a_horizontal = (a1.y - a2.y).abs() < 0.01;
        let b_horizontal = (b1.y - b2.y).abs() < 0.01;

        // Only count crossings between perpendicular segments
        if a_horizontal == b_horizontal {
            return false;
        }

        if a_horizontal {
            // A is horizontal, B is vertical
            let ax_min = a1.x.min(a2.x);
            let ax_max = a1.x.max(a2.x);
            let by_min = b1.y.min(b2.y);
            let by_max = b1.y.max(b2.y);

            b1.x > ax_min && b1.x < ax_max && a1.y > by_min && a1.y < by_max
        } else {
            // A is vertical, B is horizontal
            let ay_min = a1.y.min(a2.y);
            let ay_max = a1.y.max(a2.y);
            let bx_min = b1.x.min(b2.x);
            let bx_max = b1.x.max(b2.x);

            b1.y > ay_min && b1.y < ay_max && a1.x > bx_min && a1.x < bx_max
        }
    }
}

impl Default for OrthogonalRouter {
    fn default() -> Self {
        Self::new()
    }
}

/// Auto-rotation optimizer to minimize wire crossings
pub struct RotationOptimizer {
    router: OrthogonalRouter,
}

impl RotationOptimizer {
    pub fn new() -> Self {
        Self {
            router: OrthogonalRouter::new(),
        }
    }

    /// Find optimal rotations for all nodes to minimize crossings
    /// Returns a map of node_id -> optimal_rotation
    pub fn optimize_rotations(
        &self,
        nodes: &HashMap<String, NodeInstance>,
        connections: &[Connection],
    ) -> HashMap<String, u8> {
        let node_ids: Vec<String> = nodes.keys().cloned().collect();
        let mut best_rotations: HashMap<String, u8> = nodes
            .iter()
            .map(|(id, node)| (id.clone(), node.rotation))
            .collect();

        let mut best_crossings = self.evaluate_crossings(nodes, connections, &best_rotations);
        let mut best_length = self.evaluate_wire_length(nodes, connections, &best_rotations);

        // Greedy optimization: try rotating each node and keep improvements
        let mut improved = true;
        let mut iterations = 0;
        const MAX_ITERATIONS: usize = 100;

        while improved && iterations < MAX_ITERATIONS {
            improved = false;
            iterations += 1;

            for node_id in &node_ids {
                let current_rotation = best_rotations[node_id];

                // Try all 4 rotations
                for rotation in 0..4u8 {
                    if rotation == current_rotation {
                        continue;
                    }

                    let mut test_rotations = best_rotations.clone();
                    test_rotations.insert(node_id.clone(), rotation);

                    let crossings = self.evaluate_crossings(nodes, connections, &test_rotations);
                    let length = self.evaluate_wire_length(nodes, connections, &test_rotations);

                    // Prefer fewer crossings, then shorter wires
                    if crossings < best_crossings
                        || (crossings == best_crossings && length < best_length - 1.0)
                    {
                        best_crossings = crossings;
                        best_length = length;
                        best_rotations.insert(node_id.clone(), rotation);
                        improved = true;
                    }
                }
            }
        }

        best_rotations
    }

    /// Apply optimal rotations to nodes
    pub fn apply_rotations(
        nodes: &mut HashMap<String, NodeInstance>,
        rotations: &HashMap<String, u8>,
    ) {
        for (node_id, rotation) in rotations {
            if let Some(node) = nodes.get_mut(node_id) {
                node.rotation = *rotation;
            }
        }
    }

    fn evaluate_crossings(
        &self,
        nodes: &HashMap<String, NodeInstance>,
        connections: &[Connection],
        rotations: &HashMap<String, u8>,
    ) -> usize {
        // Create temporary nodes with test rotations
        let mut test_nodes = nodes.clone();
        for (id, rotation) in rotations {
            if let Some(node) = test_nodes.get_mut(id) {
                node.rotation = *rotation;
            }
        }
        self.router.count_crossings(connections, &test_nodes)
    }

    fn evaluate_wire_length(
        &self,
        nodes: &HashMap<String, NodeInstance>,
        connections: &[Connection],
        rotations: &HashMap<String, u8>,
    ) -> f32 {
        let mut test_nodes = nodes.clone();
        for (id, rotation) in rotations {
            if let Some(node) = test_nodes.get_mut(id) {
                node.rotation = *rotation;
            }
        }
        self.router.calculate_total_wire_length(connections, &test_nodes)
    }
}

impl Default for RotationOptimizer {
    fn default() -> Self {
        Self::new()
    }
}
