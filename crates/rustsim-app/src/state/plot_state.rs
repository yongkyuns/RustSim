//! Plot state management - contains plot data and visualization settings.

use std::collections::{HashMap, VecDeque};
use crate::ui::PlotSettings;

/// Simulation result snapshot for ghost traces
#[derive(Clone)]
pub struct SimulationResult {
    pub plot_data: Vec<(f64, Vec<f64>)>,
    pub plot_labels: Vec<String>,
}

/// State for managing plot data and settings
pub struct PlotState {
    /// Recorded plot data: (time, outputs)
    pub plot_data: Vec<(f64, Vec<f64>)>,

    /// Plot signal labels
    pub plot_labels: Vec<String>,

    /// Per-node plot data for Scope nodes: node_id -> Vec<(time, values)>
    pub scope_data: HashMap<String, Vec<(f64, Vec<f64>)>>,

    /// Per-node spectrum data for Spectrum nodes: node_id -> Vec<(freq, magnitude)>
    pub spectrum_data: HashMap<String, Vec<(f64, Vec<f64>)>>,

    /// Plot visualization settings
    pub plot_settings: PlotSettings,

    /// Result history for ghost traces (max 6)
    pub result_history: VecDeque<SimulationResult>,
}

impl PlotState {
    pub fn new() -> Self {
        Self {
            plot_data: Vec::new(),
            plot_labels: Vec::new(),
            scope_data: HashMap::new(),
            spectrum_data: HashMap::new(),
            plot_settings: PlotSettings::default(),
            result_history: VecDeque::new(),
        }
    }

    /// Clear all plot data
    pub fn clear(&mut self) {
        self.plot_data.clear();
        self.plot_labels.clear();
        self.scope_data.clear();
        self.spectrum_data.clear();
    }

    /// Add a data point to the main plot
    pub fn add_plot_point(&mut self, time: f64, values: Vec<f64>) {
        self.plot_data.push((time, values));
    }

    /// Set plot labels
    pub fn set_labels(&mut self, labels: Vec<String>) {
        self.plot_labels = labels;
    }

    /// Add plot labels
    pub fn add_labels(&mut self, labels: impl IntoIterator<Item = String>) {
        self.plot_labels.extend(labels);
    }

    /// Add data point to a scope node
    pub fn add_scope_data(&mut self, node_id: String, time: f64, values: Vec<f64>) {
        self.scope_data
            .entry(node_id)
            .or_insert_with(Vec::new)
            .push((time, values));
    }

    /// Add spectrum data for a node
    pub fn add_spectrum_data(&mut self, node_id: String, freq: f64, magnitudes: Vec<f64>) {
        self.spectrum_data
            .entry(node_id)
            .or_insert_with(Vec::new)
            .push((freq, magnitudes));
    }

    /// Save current simulation result to history for ghost traces
    pub fn save_to_history(&mut self) {
        let result = SimulationResult {
            plot_data: self.plot_data.clone(),
            plot_labels: self.plot_labels.clone(),
        };

        self.result_history.push_front(result);

        // Keep max 6 results
        while self.result_history.len() > 6 {
            self.result_history.pop_back();
        }
    }

    /// Clear result history
    pub fn clear_history(&mut self) {
        self.result_history.clear();
    }

    /// Get scope data for a specific node
    pub fn get_scope_data(&self, node_id: &str) -> Option<&Vec<(f64, Vec<f64>)>> {
        self.scope_data.get(node_id)
    }

    /// Get spectrum data for a specific node
    pub fn get_spectrum_data(&self, node_id: &str) -> Option<&Vec<(f64, Vec<f64>)>> {
        self.spectrum_data.get(node_id)
    }

    /// Reset plot state
    pub fn reset(&mut self) {
        self.plot_data.clear();
        self.plot_labels.clear();
        self.scope_data.clear();
        self.spectrum_data.clear();
        self.result_history.clear();
    }
}

impl Default for PlotState {
    fn default() -> Self {
        Self::new()
    }
}
