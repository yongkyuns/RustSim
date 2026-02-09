//! Tests for wire routing algorithms.

use super::*;
use rustsim_types::{Connection, NodeInstance, Position};
use std::collections::HashMap;

#[test]
fn test_astar_finds_direct_path() {
    let router = AStarRouter::new(10.0);
    let start = Position::new(0.0, 0.0);
    let end = Position::new(100.0, 0.0);
    let obstacles = vec![];

    let path = router.find_path(start, end, &obstacles);

    // Should have at least start and end
    assert!(path.len() >= 2);
    assert_eq!(path[0], start);
    assert_eq!(path[path.len() - 1], end);

    // All points should be orthogonal (horizontal or vertical moves)
    for i in 0..path.len() - 1 {
        let dx = (path[i + 1].x - path[i].x).abs();
        let dy = (path[i + 1].y - path[i].y).abs();
        assert!(
            dx < 0.01 || dy < 0.01,
            "Path should be orthogonal at segment {}",
            i
        );
    }
}

#[test]
fn test_astar_avoids_obstacles() {
    let router = AStarRouter::new(10.0);
    let start = Position::new(0.0, 0.0);
    let end = Position::new(100.0, 0.0);

    // Place obstacle in the middle
    let obstacles = vec![Rect::new(40.0, -20.0, 60.0, 20.0)];

    let path = router.find_path(start, end, &obstacles);

    // Should find a path
    assert!(path.len() >= 2);

    // Path should not pass through obstacle
    for point in &path {
        let inside_obstacle = point.x >= 40.0
            && point.x <= 60.0
            && point.y >= -20.0
            && point.y <= 20.0;
        assert!(
            !inside_obstacle,
            "Path should not go through obstacle at ({}, {})",
            point.x,
            point.y
        );
    }
}

#[test]
fn test_astar_orthogonal_only() {
    let router = AStarRouter::new(10.0);
    let start = Position::new(0.0, 0.0);
    let end = Position::new(50.0, 50.0);
    let obstacles = vec![];

    let path = router.find_path(start, end, &obstacles);

    // Check that all segments are orthogonal (horizontal or vertical)
    for i in 0..path.len() - 1 {
        let dx = (path[i + 1].x - path[i].x).abs();
        let dy = (path[i + 1].y - path[i].y).abs();

        // Either horizontal or vertical, not diagonal
        assert!(
            dx < 0.01 || dy < 0.01,
            "Segment {} should be orthogonal: ({}, {}) -> ({}, {})",
            i,
            path[i].x,
            path[i].y,
            path[i + 1].x,
            path[i + 1].y
        );
    }
}

#[test]
fn test_astar_start_equals_end() {
    let router = AStarRouter::new(10.0);
    let pos = Position::new(50.0, 50.0);
    let obstacles = vec![];

    let path = router.find_path(pos, pos, &obstacles);

    // Should return a single point or very short path
    assert!(!path.is_empty());
    assert_eq!(path[0], pos);
}

#[test]
fn test_astar_no_path_possible() {
    let router = AStarRouter::new(10.0);
    let start = Position::new(0.0, 0.0);
    let end = Position::new(100.0, 0.0);

    // Create a wall that completely blocks the path
    let obstacles = vec![
        Rect::new(40.0, -100.0, 60.0, 100.0), // Vertical wall
    ];

    let path = router.find_path(start, end, &obstacles);

    // Should still return some path (fallback to direct line)
    assert!(!path.is_empty());
}

#[test]
fn test_route_connection_with_directions() {
    let router = AStarRouter::new(10.0);
    let source = Position::new(0.0, 0.0);
    let target = Position::new(100.0, 50.0);
    let obstacles = vec![];

    let path = router.route_connection(
        source,
        Direction::Right,
        target,
        Direction::Left,
        &obstacles,
    );

    // Should start at source and end at target
    assert_eq!(path[0], source);
    assert_eq!(path[path.len() - 1], target);

    // Should have intermediate points
    assert!(path.len() >= 4);

    // Check orthogonality
    for i in 0..path.len() - 1 {
        let dx = (path[i + 1].x - path[i].x).abs();
        let dy = (path[i + 1].y - path[i].y).abs();
        assert!(dx < 0.01 || dy < 0.01, "Path should be orthogonal");
    }
}

#[test]
fn test_orthogonal_router_waypoints() {
    let router = OrthogonalRouter::new();
    let start = Position::new(0.0, 0.0);
    let end = Position::new(100.0, 100.0);
    let waypoints = vec![Position::new(50.0, 0.0), Position::new(50.0, 100.0)];
    let obstacles = vec![];

    let path = router.route_with_waypoints(start, end, &waypoints, &obstacles);

    // Should pass through waypoints area
    assert!(!path.is_empty());
    assert_eq!(path[0], start);
    assert_eq!(path[path.len() - 1], end);

    // Path should be reasonably close to waypoints
    // (exact match depends on grid alignment)
    let has_point_near_wp1 = path
        .iter()
        .any(|p| (p.x - 50.0).abs() < 15.0 && (p.y - 0.0).abs() < 15.0);
    let has_point_near_wp2 = path
        .iter()
        .any(|p| (p.x - 50.0).abs() < 15.0 && (p.y - 100.0).abs() < 15.0);

    assert!(
        has_point_near_wp1,
        "Path should go near first waypoint (50, 0)"
    );
    assert!(
        has_point_near_wp2,
        "Path should go near second waypoint (50, 100)"
    );
}

#[test]
fn test_orthogonal_router_no_waypoints() {
    let router = OrthogonalRouter::new();
    let start = Position::new(0.0, 0.0);
    let end = Position::new(100.0, 50.0);
    let waypoints = vec![];
    let obstacles = vec![];

    let path = router.route_with_waypoints(start, end, &waypoints, &obstacles);

    assert!(!path.is_empty());
    assert_eq!(path[0], start);
    assert_eq!(path[path.len() - 1], end);
}

#[test]
fn test_orthogonal_router_with_nodes() {
    let router = OrthogonalRouter::new();

    // Create test nodes
    let mut nodes = HashMap::new();

    let mut source = NodeInstance::new(
        "source".to_string(),
        "Constant".to_string(),
        Position::new(0.0, 0.0),
    );
    source.add_output("out");

    let mut target = NodeInstance::new(
        "target".to_string(),
        "Scope".to_string(),
        Position::new(200.0, 0.0),
    );
    target.add_input("in");

    // Add an obstacle node in the middle
    let obstacle = NodeInstance::new(
        "obstacle".to_string(),
        "Gain".to_string(),
        Position::new(80.0, -10.0),
    );

    nodes.insert("source".to_string(), source);
    nodes.insert("target".to_string(), target);
    nodes.insert("obstacle".to_string(), obstacle);

    // Create connection
    let connection = Connection::new(
        "conn1".to_string(),
        "source".to_string(),
        0,
        "target".to_string(),
        0,
    );

    let path = router.route(&connection, &nodes);

    // Should find a path
    assert!(!path.is_empty());

    // Path should be orthogonal
    for i in 0..path.len() - 1 {
        let dx = (path[i + 1].x - path[i].x).abs();
        let dy = (path[i + 1].y - path[i].y).abs();
        assert!(dx < 0.01 || dy < 0.01, "Path should be orthogonal");
    }
}

#[test]
fn test_path_simplification() {
    let router = AStarRouter::new(10.0);

    // Create a path with redundant collinear points
    let path = vec![
        Position::new(0.0, 0.0),
        Position::new(10.0, 0.0),
        Position::new(20.0, 0.0), // Collinear with prev and next
        Position::new(30.0, 0.0),
        Position::new(30.0, 10.0),
        Position::new(30.0, 20.0), // Collinear with prev and next
        Position::new(30.0, 30.0),
    ];

    let simplified = router.simplify_path(&path);

    // Should remove the collinear points
    assert!(simplified.len() < path.len());

    // Should keep start and end
    assert_eq!(simplified[0], path[0]);
    assert_eq!(simplified[simplified.len() - 1], path[path.len() - 1]);

    // Should only have corner points
    // Expected: (0,0) -> (30,0) -> (30,30)
    assert_eq!(simplified.len(), 3);
}

#[test]
fn test_grid_alignment() {
    let grid_size = 10.0;
    let router = AStarRouter::new(grid_size);

    let start = Position::new(5.0, 5.0); // Not aligned
    let end = Position::new(95.0, 95.0); // Not aligned
    let obstacles = vec![];

    let path = router.find_path(start, end, &obstacles);

    // All intermediate points (except start/end) should be grid-aligned
    for i in 0..path.len() {
        let x_aligned = (path[i].x / grid_size).round() * grid_size;
        let y_aligned = (path[i].y / grid_size).round() * grid_size;

        // Points should be close to grid positions
        assert!(
            (path[i].x - x_aligned).abs() < 0.1,
            "Point x={} should be grid-aligned",
            path[i].x
        );
        assert!(
            (path[i].y - y_aligned).abs() < 0.1,
            "Point y={} should be grid-aligned",
            path[i].y
        );
    }
}

#[test]
fn test_multiple_obstacles() {
    let router = AStarRouter::new(10.0);
    let start = Position::new(0.0, 0.0);
    let end = Position::new(100.0, 0.0);

    // Create multiple obstacles forming a maze
    let obstacles = vec![
        Rect::new(20.0, -10.0, 30.0, 10.0),
        Rect::new(40.0, -10.0, 50.0, 10.0),
        Rect::new(60.0, -10.0, 70.0, 10.0),
    ];

    let path = router.find_path(start, end, &obstacles);

    // Should find a path
    assert!(!path.is_empty());

    // Verify no collisions with obstacles
    for point in &path {
        for obs in &obstacles {
            let inside = point.x >= obs.min_x - 1.0
                && point.x <= obs.max_x + 1.0
                && point.y >= obs.min_y - 1.0
                && point.y <= obs.max_y + 1.0;
            assert!(
                !inside,
                "Path point ({}, {}) should not be inside obstacle",
                point.x,
                point.y
            );
        }
    }
}
