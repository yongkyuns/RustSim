//! History state management - contains undo/redo functionality.

use rustsim_types::{SimulationGraph, SimulationSettings};

/// State snapshot for undo/redo
#[derive(Clone)]
struct UndoState {
    graph: SimulationGraph,
    settings: SimulationSettings,
}

/// State for managing undo/redo history
pub struct HistoryState {
    /// Undo stack
    undo_stack: Vec<UndoState>,

    /// Redo stack
    redo_stack: Vec<UndoState>,

    /// Maximum undo stack size
    max_undo_size: usize,
}

impl HistoryState {
    pub fn new() -> Self {
        Self {
            undo_stack: Vec::new(),
            redo_stack: Vec::new(),
            max_undo_size: 50,
        }
    }

    /// Save current state to undo stack
    pub fn save_state(&mut self, graph: &SimulationGraph, settings: &SimulationSettings) {
        self.undo_stack.push(UndoState {
            graph: graph.clone(),
            settings: settings.clone(),
        });
        self.redo_stack.clear();

        // Limit undo stack size
        if self.undo_stack.len() > self.max_undo_size {
            self.undo_stack.remove(0);
        }
    }

    /// Undo the last operation
    pub fn undo(
        &mut self,
        current_graph: &SimulationGraph,
        current_settings: &SimulationSettings,
    ) -> Option<(SimulationGraph, SimulationSettings)> {
        if let Some(state) = self.undo_stack.pop() {
            // Save current state to redo stack
            self.redo_stack.push(UndoState {
                graph: current_graph.clone(),
                settings: current_settings.clone(),
            });

            Some((state.graph, state.settings))
        } else {
            None
        }
    }

    /// Redo the last undone operation
    pub fn redo(
        &mut self,
        current_graph: &SimulationGraph,
        current_settings: &SimulationSettings,
    ) -> Option<(SimulationGraph, SimulationSettings)> {
        if let Some(state) = self.redo_stack.pop() {
            // Save current state to undo stack
            self.undo_stack.push(UndoState {
                graph: current_graph.clone(),
                settings: current_settings.clone(),
            });

            Some((state.graph, state.settings))
        } else {
            None
        }
    }

    /// Check if undo is available
    pub fn can_undo(&self) -> bool {
        !self.undo_stack.is_empty()
    }

    /// Check if redo is available
    pub fn can_redo(&self) -> bool {
        !self.redo_stack.is_empty()
    }

    /// Clear all history
    pub fn clear(&mut self) {
        self.undo_stack.clear();
        self.redo_stack.clear();
    }

    /// Get undo stack size
    pub fn undo_size(&self) -> usize {
        self.undo_stack.len()
    }

    /// Get redo stack size
    pub fn redo_size(&self) -> usize {
        self.redo_stack.len()
    }
}

impl Default for HistoryState {
    fn default() -> Self {
        Self::new()
    }
}
