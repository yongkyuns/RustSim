//! A* pathfinding for grid-based wire routing.

use rustsim_types::Position;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};

/// Direction a wire exits from a port
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    Left,
    Right,
    Up,
    Down,
}

/// Rectangle for obstacle representation
#[derive(Debug, Clone, Copy)]
pub struct Rect {
    pub min_x: f32,
    pub min_y: f32,
    pub max_x: f32,
    pub max_y: f32,
}

impl Rect {
    pub fn new(min_x: f32, min_y: f32, max_x: f32, max_y: f32) -> Self {
        Self {
            min_x,
            min_y,
            max_x,
            max_y,
        }
    }

    /// Check if a point is inside the rectangle (with small margin)
    pub fn contains(&self, pos: &Position, margin: f32) -> bool {
        pos.x >= self.min_x - margin
            && pos.x <= self.max_x + margin
            && pos.y >= self.min_y - margin
            && pos.y <= self.max_y + margin
    }
}

/// Grid position for pathfinding
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct GridPos {
    x: i32,
    y: i32,
}

impl GridPos {
    fn new(x: i32, y: i32) -> Self {
        Self { x, y }
    }

    fn from_position(pos: Position, grid_size: f32) -> Self {
        Self {
            x: (pos.x / grid_size).round() as i32,
            y: (pos.y / grid_size).round() as i32,
        }
    }

    fn to_position(&self, grid_size: f32) -> Position {
        Position::new(self.x as f32 * grid_size, self.y as f32 * grid_size)
    }

    fn neighbors(&self) -> [GridPos; 4] {
        [
            GridPos::new(self.x + 1, self.y), // Right
            GridPos::new(self.x - 1, self.y), // Left
            GridPos::new(self.x, self.y + 1), // Down
            GridPos::new(self.x, self.y - 1), // Up
        ]
    }

    fn manhattan_distance(&self, other: &GridPos) -> u32 {
        ((self.x - other.x).abs() + (self.y - other.y).abs()) as u32
    }
}

/// Node in the A* search
#[derive(Debug, Clone)]
struct Node {
    pos: GridPos,
    g_score: u32, // Cost from start
    f_score: u32, // g_score + heuristic
}

impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.f_score == other.f_score
    }
}

impl Eq for Node {}

impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Node {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reverse ordering for min-heap behavior
        other.f_score.cmp(&self.f_score)
    }
}

/// A* router for orthogonal wire routing
pub struct AStarRouter {
    grid_size: f32,
}

impl AStarRouter {
    /// Create a new A* router with the specified grid size
    pub fn new(grid_size: f32) -> Self {
        Self { grid_size }
    }

    /// Find a path from start to end avoiding obstacles
    pub fn find_path(
        &self,
        start: Position,
        end: Position,
        obstacles: &[Rect],
    ) -> Vec<Position> {
        let start_grid = GridPos::from_position(start, self.grid_size);
        let end_grid = GridPos::from_position(end, self.grid_size);

        // If start and end are the same, return empty path
        if start_grid == end_grid {
            return vec![start];
        }

        // Initialize data structures
        let mut open_set = BinaryHeap::new();
        let mut came_from: HashMap<GridPos, GridPos> = HashMap::new();
        let mut g_score: HashMap<GridPos, u32> = HashMap::new();
        let mut closed_set: HashSet<GridPos> = HashSet::new();

        // Add start node
        g_score.insert(start_grid, 0);
        let h = start_grid.manhattan_distance(&end_grid);
        open_set.push(Node {
            pos: start_grid,
            g_score: 0,
            f_score: h,
        });

        while let Some(current) = open_set.pop() {
            // Check if we reached the goal
            if current.pos == end_grid {
                return self.reconstruct_path(&came_from, current.pos, start_grid);
            }

            // Skip if already processed
            if closed_set.contains(&current.pos) {
                continue;
            }

            closed_set.insert(current.pos);

            // Explore neighbors
            for neighbor_pos in current.pos.neighbors() {
                // Skip if already processed
                if closed_set.contains(&neighbor_pos) {
                    continue;
                }

                // Check if neighbor is blocked by obstacle
                let neighbor_world = neighbor_pos.to_position(self.grid_size);
                if self.is_blocked(&neighbor_world, obstacles) {
                    continue;
                }

                // Calculate tentative g_score
                let tentative_g = current.g_score + 1;

                // Check if this path is better
                if tentative_g < *g_score.get(&neighbor_pos).unwrap_or(&u32::MAX) {
                    came_from.insert(neighbor_pos, current.pos);
                    g_score.insert(neighbor_pos, tentative_g);

                    let h = neighbor_pos.manhattan_distance(&end_grid);
                    open_set.push(Node {
                        pos: neighbor_pos,
                        g_score: tentative_g,
                        f_score: tentative_g + h,
                    });
                }
            }
        }

        // No path found, return direct line
        vec![start, end]
    }

    /// Route a connection with exit directions
    pub fn route_connection(
        &self,
        source_pos: Position,
        source_dir: Direction,
        target_pos: Position,
        target_dir: Direction,
        obstacles: &[Rect],
    ) -> Vec<Position> {
        // Create intermediate points based on exit directions
        let exit_offset = self.grid_size * 2.0;

        let source_exit = match source_dir {
            Direction::Right => Position::new(source_pos.x + exit_offset, source_pos.y),
            Direction::Left => Position::new(source_pos.x - exit_offset, source_pos.y),
            Direction::Up => Position::new(source_pos.x, source_pos.y - exit_offset),
            Direction::Down => Position::new(source_pos.x, source_pos.y + exit_offset),
        };

        let target_entry = match target_dir {
            Direction::Right => Position::new(target_pos.x + exit_offset, target_pos.y),
            Direction::Left => Position::new(target_pos.x - exit_offset, target_pos.y),
            Direction::Up => Position::new(target_pos.x, target_pos.y - exit_offset),
            Direction::Down => Position::new(target_pos.x, target_pos.y + exit_offset),
        };

        // Find path from exit point to entry point
        let middle_path = self.find_path(source_exit, target_entry, obstacles);

        // Combine: source -> exit -> middle path -> entry -> target
        let mut path = vec![source_pos];
        path.push(source_exit);
        path.extend(middle_path);
        path.push(target_entry);
        path.push(target_pos);

        // Simplify the path to remove redundant points
        self.simplify_path(&path)
    }

    /// Check if a position is blocked by any obstacle
    fn is_blocked(&self, pos: &Position, obstacles: &[Rect]) -> bool {
        // Use a small margin to avoid edge cases
        let margin = self.grid_size * 0.1;
        obstacles.iter().any(|obs| obs.contains(pos, margin))
    }

    /// Reconstruct the path from the came_from map
    fn reconstruct_path(
        &self,
        came_from: &HashMap<GridPos, GridPos>,
        mut current: GridPos,
        start: GridPos,
    ) -> Vec<Position> {
        let mut path = vec![current.to_position(self.grid_size)];

        while current != start {
            if let Some(&prev) = came_from.get(&current) {
                current = prev;
                path.push(current.to_position(self.grid_size));
            } else {
                break;
            }
        }

        path.reverse();
        self.simplify_path(&path)
    }

    /// Simplify path by removing collinear points
    pub fn simplify_path(&self, path: &[Position]) -> Vec<Position> {
        if path.len() <= 2 {
            return path.to_vec();
        }

        let mut simplified = vec![path[0]];

        for i in 1..path.len() - 1 {
            let prev = simplified.last().unwrap();
            let current = &path[i];
            let next = &path[i + 1];

            // Check if current point is collinear with prev and next
            let dx1 = current.x - prev.x;
            let dy1 = current.y - prev.y;
            let dx2 = next.x - current.x;
            let dy2 = next.y - current.y;

            // If not on the same line (either horizontal or vertical), keep the point
            let is_horizontal_1 = dy1.abs() < 0.01;
            let is_horizontal_2 = dy2.abs() < 0.01;
            let is_vertical_1 = dx1.abs() < 0.01;
            let is_vertical_2 = dx2.abs() < 0.01;

            if !(is_horizontal_1 && is_horizontal_2) && !(is_vertical_1 && is_vertical_2) {
                simplified.push(*current);
            }
        }

        // Always add the last point
        simplified.push(path[path.len() - 1]);
        simplified
    }
}

impl Default for AStarRouter {
    fn default() -> Self {
        Self::new(10.0)
    }
}
