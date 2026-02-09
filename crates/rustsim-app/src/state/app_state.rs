//! Main application state combining all state modules.

use super::{CanvasState, CompilationStatus, GraphState, HistoryState, PlotState, SimulationState};
use crate::examples;
use crate::layout;
use rustsim_types::Position;

#[cfg(not(target_arch = "wasm32"))]
use rustsim_codegen::CodeGenerator;

#[cfg(not(target_arch = "wasm32"))]
use std::sync::mpsc;

#[cfg(not(target_arch = "wasm32"))]
use super::CompilationMessage;

/// Unified application state
pub struct AppState {
    /// Graph state (nodes, connections, block types)
    pub graph_state: GraphState,

    /// Canvas state (pan, zoom, selection)
    pub canvas_state: CanvasState,

    /// Simulation state (running, data, compiled sim)
    pub simulation_state: SimulationState,

    /// Plot state (settings, ghost traces)
    pub plot_state: PlotState,

    /// History state (undo/redo)
    pub history_state: HistoryState,

    /// Node being edited (double-clicked)
    pub editing_node: Option<String>,
}

impl AppState {
    pub fn new() -> Self {
        // Load PID Controller as default example
        let graph_state = GraphState::with_graph(examples::pid_controller());

        Self {
            graph_state,
            canvas_state: CanvasState::new(),
            simulation_state: SimulationState::new(),
            plot_state: PlotState::new(),
            history_state: HistoryState::new(),
            editing_node: None,
        }
    }

    /// Delegate method to access graph
    pub fn graph(&self) -> &rustsim_types::SimulationGraph {
        &self.graph_state.graph
    }

    /// Delegate method to access graph mutably
    pub fn graph_mut(&mut self) -> &mut rustsim_types::SimulationGraph {
        &mut self.graph_state.graph
    }

    /// Delegate method to get block type
    pub fn get_block_type(&self, name: &str) -> Option<&rustsim_types::BlockTypeDefinition> {
        self.graph_state.get_block_type(name)
    }

    // File operations

    pub fn new_file(&mut self) {
        self.graph_state.reset();
        self.simulation_state.settings = rustsim_types::SimulationSettings::default();
        self.canvas_state.clear_selection();
        self.plot_state.clear();
        self.history_state.clear();
    }

    pub fn open_file(&mut self) {
        // TODO: Implement file dialog
    }

    pub fn save_file(&mut self) {
        // TODO: Implement file save
    }

    // Example loading

    /// Load an example simulation by name
    pub fn load_example(&mut self, name: &str) {
        if let Some(example_graph) = examples::load_example(name) {
            // Stop any running simulation
            self.simulation_state.stop();

            // Reset state
            self.graph_state = GraphState::with_graph(example_graph);
            self.simulation_state.reset();
            self.canvas_state.clear_selection();
            self.plot_state.clear();
            self.history_state.clear();

            // Reset compilation state
            self.simulation_state.use_compiled_mode = false;
            self.simulation_state.invalidate_compiled_sim();
            self.simulation_state.compilation_log.clear();
        }
    }

    /// Get list of available examples
    pub fn available_examples() -> Vec<(&'static str, &'static str)> {
        examples::list_examples()
    }

    // Layout operations

    /// Auto-arrange blocks using layered layout algorithm
    pub fn auto_layout(&mut self) {
        layout::auto_layout(&mut self.graph_state.graph);
    }

    // Graph operations

    /// Add a node to the graph
    pub fn add_node(&mut self, block_type: &str, position: Position) -> String {
        self.save_undo_state();
        self.simulation_state.invalidate_compiled_sim();
        self.graph_state.add_node(block_type, position)
    }

    /// Add a connection between nodes
    pub fn add_connection(
        &mut self,
        source_node: &str,
        source_port: usize,
        target_node: &str,
        target_port: usize,
    ) -> String {
        self.save_undo_state();
        self.simulation_state.invalidate_compiled_sim();
        self.graph_state.add_connection(source_node, source_port, target_node, target_port)
    }

    /// Delete selected items
    pub fn delete_selected(&mut self) {
        if self.canvas_state.selected_nodes.is_empty() && self.canvas_state.selected_connections.is_empty() {
            return;
        }

        self.save_undo_state();
        self.simulation_state.invalidate_compiled_sim();

        for node_id in self.canvas_state.selected_nodes.drain() {
            self.graph_state.remove_node(&node_id);
        }

        for conn_id in self.canvas_state.selected_connections.drain() {
            self.graph_state.remove_connection(&conn_id);
        }
    }

    /// Rotate selected nodes by 90 degrees
    pub fn rotate_selected(&mut self) {
        if self.canvas_state.selected_nodes.is_empty() {
            return;
        }

        self.save_undo_state();

        for node_id in &self.canvas_state.selected_nodes.clone() {
            if let Some(node) = self.graph_state.graph.get_node_mut(node_id) {
                node.rotation = (node.rotation + 1) % 4;
            }
        }
    }

    /// Clear selection
    pub fn clear_selection(&mut self) {
        self.canvas_state.clear_selection();
    }

    // Port management

    /// Check if an input port can be added
    pub fn can_add_input(&self, node_id: &str) -> bool {
        self.graph_state.can_add_input(node_id)
    }

    /// Check if an input port can be removed
    pub fn can_remove_input(&self, node_id: &str) -> bool {
        self.graph_state.can_remove_input(node_id)
    }

    /// Check if an output port can be added
    pub fn can_add_output(&self, node_id: &str) -> bool {
        self.graph_state.can_add_output(node_id)
    }

    /// Check if an output port can be removed
    pub fn can_remove_output(&self, node_id: &str) -> bool {
        self.graph_state.can_remove_output(node_id)
    }

    /// Check if a node has sync_ports
    pub fn has_sync_ports(&self, node_id: &str) -> bool {
        self.graph_state.has_sync_ports(node_id)
    }

    /// Add an input port
    pub fn add_input_port(&mut self, node_id: &str) -> bool {
        self.save_undo_state();
        self.graph_state.add_input_port(node_id)
    }

    /// Remove an input port
    pub fn remove_input_port(&mut self, node_id: &str) -> bool {
        self.save_undo_state();
        self.graph_state.remove_input_port(node_id)
    }

    /// Add an output port
    pub fn add_output_port(&mut self, node_id: &str) -> bool {
        self.save_undo_state();
        self.graph_state.add_output_port(node_id)
    }

    /// Remove an output port
    pub fn remove_output_port(&mut self, node_id: &str) -> bool {
        self.save_undo_state();
        self.graph_state.remove_output_port(node_id)
    }

    // View operations

    /// Zoom in
    pub fn zoom_in(&mut self) {
        self.canvas_state.zoom_in();
    }

    /// Zoom out
    pub fn zoom_out(&mut self) {
        self.canvas_state.zoom_out();
    }

    /// Fit all nodes in view
    pub fn fit_view(&mut self) {
        if self.graph_state.graph.nodes.is_empty() {
            self.canvas_state.reset_pan();
            self.canvas_state.reset_zoom();
            return;
        }

        // Calculate bounding box
        let mut min_x = f32::MAX;
        let mut min_y = f32::MAX;
        let mut max_x = f32::MIN;
        let mut max_y = f32::MIN;

        for node in self.graph_state.graph.nodes.values() {
            min_x = min_x.min(node.position.x);
            min_y = min_y.min(node.position.y);
            max_x = max_x.max(node.position.x + 80.0);
            max_y = max_y.max(node.position.y + 40.0);
        }

        // Center the view
        self.canvas_state.pan = Position::new(-(min_x + max_x) / 2.0, -(min_y + max_y) / 2.0);
        self.canvas_state.zoom = 1.0;
    }

    // Undo/Redo

    fn save_undo_state(&mut self) {
        self.history_state.save_state(&self.graph_state.graph, &self.simulation_state.settings);
    }

    pub fn undo(&mut self) {
        if let Some((graph, settings)) = self.history_state.undo(&self.graph_state.graph, &self.simulation_state.settings) {
            self.graph_state.graph = graph;
            self.simulation_state.settings = settings;
        }
    }

    pub fn redo(&mut self) {
        if let Some((graph, settings)) = self.history_state.redo(&self.graph_state.graph, &self.simulation_state.settings) {
            self.graph_state.graph = graph;
            self.simulation_state.settings = settings;
        }
    }

    // Simulation control

    /// Check if simulation is running
    pub fn is_running(&self) -> bool {
        self.simulation_state.is_running()
    }

    /// Get average step time
    pub fn average_step_time(&self) -> Option<std::time::Duration> {
        self.simulation_state.average_step_time()
    }

    /// Get step count
    pub fn get_step_count(&self) -> u64 {
        self.simulation_state.get_step_count()
    }

    /// Start async compilation of the current simulation graph
    #[cfg(not(target_arch = "wasm32"))]
    pub fn compile_simulation(&mut self) -> Result<(), String> {
        use crate::compiler::CompiledSimulation;

        // Clear previous compilation log
        self.simulation_state.compilation_log.clear();

        // Set status to compiling
        self.simulation_state.compilation_status = CompilationStatus::Compiling;
        self.simulation_state.compilation_log.push("Starting compilation...".to_string());

        // Generate code with parameter mapping
        self.simulation_state.compilation_log.push("Generating Rust code from simulation graph...".to_string());
        let generator = CodeGenerator::new()
            .with_c_abi(false)
            .with_library_mode(false);
        let (code, param_mapping) = generator.generate_with_params(&self.graph_state.graph, &self.simulation_state.settings);

        self.simulation_state.compilation_log.push(format!(
            "Code generation complete. {} runtime parameters.",
            param_mapping.params.len()
        ));

        // Count states, inputs, and outputs
        let (num_inputs, num_outputs, num_states) = self.count_io_for_compilation();
        self.simulation_state.compilation_log.push(format!(
            "Simulation dimensions: {} inputs, {} outputs, {} states",
            num_inputs, num_outputs, num_states
        ));

        // Get workspace root
        self.simulation_state.compilation_log.push("Locating workspace root...".to_string());
        let workspace_root = std::env::current_dir()
            .ok()
            .and_then(|mut p| {
                // Try to find workspace root by looking for Cargo.toml with [workspace]
                loop {
                    let cargo_toml = p.join("Cargo.toml");
                    if cargo_toml.exists() {
                        if let Ok(contents) = std::fs::read_to_string(&cargo_toml) {
                            if contents.contains("[workspace]") {
                                return Some(p);
                            }
                        }
                    }
                    if !p.pop() {
                        break;
                    }
                }
                None
            })
            .ok_or_else(|| "Could not find workspace root".to_string())?;

        let workspace_root_str = workspace_root
            .to_str()
            .ok_or_else(|| "Invalid workspace path".to_string())?
            .to_string();

        self.simulation_state.compilation_log.push(format!("Workspace root: {}", workspace_root_str));

        // Create channel for communication with compilation thread
        let (tx, rx) = mpsc::channel();
        self.simulation_state.set_compilation_receiver(rx);

        // Spawn compilation in background thread
        std::thread::spawn(move || {
            // Send log messages
            let _ = tx.send(CompilationMessage::Log("Creating temporary build directory...".to_string()));
            let _ = tx.send(CompilationMessage::Log("Writing generated code to src/lib.rs...".to_string()));
            let _ = tx.send(CompilationMessage::Log("Creating Cargo.toml manifest...".to_string()));
            let _ = tx.send(CompilationMessage::Log("Running 'cargo build --release --lib'...".to_string()));
            let _ = tx.send(CompilationMessage::Log("This may take a while on first compilation...".to_string()));

            // Perform compilation with parameter mapping
            match CompiledSimulation::with_params(&code, &workspace_root_str, num_states, num_inputs, num_outputs, param_mapping) {
                Ok(mut sim) => {
                    let _ = tx.send(CompilationMessage::Log("Compilation successful!".to_string()));
                    let _ = tx.send(CompilationMessage::Log("Loading compiled library...".to_string()));
                    sim.init();
                    let _ = tx.send(CompilationMessage::Log("Initializing simulation state...".to_string()));
                    let _ = tx.send(CompilationMessage::Success(sim));
                }
                Err(e) => {
                    let _ = tx.send(CompilationMessage::Log(format!("Compilation failed: {}", e)));
                    let _ = tx.send(CompilationMessage::Error(e.to_string()));
                }
            }
        });

        self.simulation_state.compilation_log.push("Background compilation thread started.".to_string());
        Ok(())
    }

    /// Poll compilation progress (call this regularly from UI update)
    #[cfg(not(target_arch = "wasm32"))]
    pub fn poll_compilation(&mut self) {
        use crate::compiler::CompiledSimulation;

        let mut should_close = false;
        let mut compiled_sim: Option<CompiledSimulation> = None;
        let mut status: Option<CompilationStatus> = None;

        if let Some(rx) = self.simulation_state.take_compilation_receiver() {
            // Process all available messages without blocking
            while let Ok(msg) = rx.try_recv() {
                match msg {
                    CompilationMessage::Log(text) => {
                        self.simulation_state.compilation_log.push(text);
                    }
                    CompilationMessage::Success(sim) => {
                        compiled_sim = Some(sim);
                        status = Some(CompilationStatus::Ready);
                        self.simulation_state.compilation_log.push("Compilation complete and ready!".to_string());
                        should_close = true;
                    }
                    CompilationMessage::Error(err) => {
                        status = Some(CompilationStatus::Error(err.clone()));
                        self.simulation_state.compilation_log.push(format!("Error: {}", err));
                        should_close = true;
                    }
                }
            }

            if !should_close {
                // Put the receiver back
                self.simulation_state.set_compilation_receiver(rx);
            }
        }

        // Update state after processing messages
        if let Some(sim) = compiled_sim {
            self.simulation_state.set_compiled_sim(sim);
        }
        if let Some(s) = status {
            self.simulation_state.compilation_status = s;
        }
    }

    /// Count inputs, outputs, and states for compiled simulation
    #[cfg(not(target_arch = "wasm32"))]
    fn count_io_for_compilation(&self) -> (usize, usize, usize) {
        let num_inputs = 0;
        let mut num_outputs = 0;
        let mut num_states = 0;

        for node in self.graph_state.graph.nodes.values() {
            match node.block_type.as_str() {
                "Scope" | "Spectrum" => {
                    // Recording blocks contribute to outputs
                    num_outputs += node.inputs.len().max(1);
                }
                "Integrator" => {
                    num_states += 1;
                    num_outputs += 1;
                }
                "Differentiator" => {
                    num_states += 1;
                    num_outputs += 1;
                }
                "Constant" | "Sinusoidal" | "Step" | "Ramp" => {
                    num_outputs += 1;
                }
                _ => {
                    num_outputs += node.outputs.len().max(1);
                }
            }
        }

        (num_inputs.max(1), num_outputs.max(1), num_states.max(1))
    }

    /// On WASM, compilation is not supported yet
    #[cfg(target_arch = "wasm32")]
    pub fn compile_simulation(&mut self) -> Result<(), String> {
        Err("Compilation not supported on WASM yet".to_string())
    }

    #[cfg(target_arch = "wasm32")]
    pub fn poll_compilation(&mut self) {
        // No-op on WASM
    }

    pub fn run_simulation(&mut self) {
        // If compiled mode is enabled and not yet compiled, compile first
        #[cfg(not(target_arch = "wasm32"))]
        if self.simulation_state.use_compiled_mode && self.simulation_state.compiled_sim().is_none() {
            if let Err(e) = self.compile_simulation() {
                eprintln!("Compilation failed: {}. Falling back to interpreter.", e);
                self.simulation_state.use_compiled_mode = false;
            }
        }

        self.simulation_state.start();
        self.plot_state.clear();

        // Reset compiled simulation and sync parameters if in compiled mode
        #[cfg(not(target_arch = "wasm32"))]
        if self.simulation_state.use_compiled_mode {
            // Sync parameters from graph before running
            self.sync_compiled_params();

            if let Some(sim) = self.simulation_state.compiled_sim_mut() {
                // Reset time and internal states (but not params)
                sim.state_mut().time = 0.0;
                sim.state_mut().states.fill(0.0);
                sim.state_mut().outputs.fill(0.0);
            }
        }

        // Create block instances from graph
        self.simulation_state.create_blocks_from_graph(&self.graph_state.graph);

        // Collect signal labels from Scope blocks
        for (_id, node) in &self.graph_state.graph.nodes {
            if node.block_type == "Scope" {
                // Add labels from Scope input ports
                for input in &node.inputs {
                    self.plot_state.add_labels(std::iter::once(input.name.clone()));
                }
            }
        }

        // If no scope blocks, create default labels for source blocks
        if self.plot_state.plot_labels.is_empty() {
            let order = self.simulation_state.execution_order().to_vec();
            for node_id in &order {
                if let Some(node) = self.graph_state.graph.get_node(node_id) {
                    if matches!(
                        node.block_type.as_str(),
                        "Constant" | "Sinusoidal" | "Step" | "Ramp" | "Integrator"
                    ) {
                        self.plot_state.add_labels(std::iter::once(format!("{} ({})", node.name, node.block_type)));
                    }
                }
            }
        }
    }

    pub fn stop_simulation(&mut self) {
        self.simulation_state.stop();

        // Save result to history when simulation completes
        if !self.plot_state.plot_data.is_empty() {
            self.plot_state.save_to_history();
        }
    }

    /// Reset simulation to t=0 without stopping
    pub fn reset_to_zero(&mut self) {
        self.simulation_state.sim_time = 0.0;
        self.plot_state.clear();
    }

    /// Sync parameters from graph to compiled simulation.
    /// Call this when parameters are changed in the UI before running the simulation.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn sync_compiled_params(&mut self) {
        if let Some(sim) = self.simulation_state.compiled_sim_mut() {
            let mapping = sim.param_mapping().clone();

            // Update each parameter from the graph
            for param in &mapping.params {
                if let Some(node) = self.graph_state.graph.get_node(&param.node_id) {
                    if let Some(value) = node.get_param(&param.name).and_then(|v| v.as_f64()) {
                        sim.set_param(&param.node_id, &param.name, value);
                    }
                }
            }

            // Reinitialize blocks with new parameters
            sim.reinit();
        }
    }

    #[cfg(target_arch = "wasm32")]
    pub fn sync_compiled_params(&mut self) {
        // No-op on WASM
    }

    pub fn step_simulation(&mut self) {
        #[cfg(not(target_arch = "wasm32"))]
        if self.simulation_state.use_compiled_mode {
            // Extract values before getting mutable borrow on simulation_state
            let t = self.simulation_state.sim_time;
            let dt = self.simulation_state.settings.dt;

            if let Some(sim) = self.simulation_state.compiled_sim_mut() {

                // Step the compiled simulation
                sim.step(dt);

                // Get all outputs from compiled simulation
                let all_outputs: Vec<f64> = (0..sim.state().outputs.len())
                    .map(|i| sim.get_output(i))
                    .collect();

                // Extract only Scope outputs (matching interpreter behavior)
                let mut scope_outputs: Vec<f64> = Vec::new();
                let mut output_idx = 0;

                // Build execution order to match code generator
                let order = self.simulation_state.execution_order().to_vec();

                // Iterate through nodes in topological order (matching code generator)
                for node_id in &order {
                    if let Some(node) = self.graph_state.graph.get_node(node_id) {
                        match node.block_type.as_str() {
                            "Scope" | "Spectrum" => {
                                // Recording blocks - collect their outputs for plotting
                                let input_count = node.inputs.len().max(1);
                                let mut node_data = Vec::new();
                                for _ in 0..input_count {
                                    if output_idx < all_outputs.len() {
                                        let val = all_outputs[output_idx];
                                        scope_outputs.push(val);
                                        node_data.push(val);
                                        output_idx += 1;
                                    }
                                }
                                // Store in scope_data for preview
                                if !node_data.is_empty() {
                                    self.plot_state.add_scope_data(node_id.clone(), t, node_data);
                                }
                            }
                            _ => {
                                // Other blocks - skip their outputs
                                output_idx += node.outputs.len().max(1);
                            }
                        }
                    }
                }

                // If no scope blocks, fall back to source block outputs
                if scope_outputs.is_empty() {
                    output_idx = 0;
                    for node_id in &order {
                        if let Some(node) = self.graph_state.graph.get_node(node_id) {
                            if matches!(
                                node.block_type.as_str(),
                                "Constant" | "Sinusoidal" | "Step" | "Ramp" | "Integrator"
                            ) {
                                if output_idx < all_outputs.len() {
                                    scope_outputs.push(all_outputs[output_idx]);
                                }
                            }
                            // Advance output_idx for all blocks
                            match node.block_type.as_str() {
                                "Scope" | "Spectrum" => output_idx += node.inputs.len().max(1),
                                _ => output_idx += node.outputs.len().max(1),
                            }
                        }
                    }
                }

                self.plot_state.add_plot_point(t, scope_outputs);
                self.simulation_state.sim_time += dt;
                return;
            } else {
                // No compiled sim available, fall back to interpreter
                self.simulation_state.use_compiled_mode = false;
            }
        }

        // Interpreter mode
        let t = self.simulation_state.sim_time;
        self.simulation_state.step(&self.graph_state.graph);

        // Collect outputs from Scope blocks
        let mut outputs: Vec<f64> = Vec::new();
        let order = self.simulation_state.execution_order().to_vec();

        for node_id in &order {
            if let Some(node) = self.graph_state.graph.get_node(node_id) {
                if node.block_type == "Scope" || node.block_type == "Spectrum" {
                    let mut node_data = Vec::new();
                    for i in 0..node.inputs.len().max(1) {
                        if let Some(val) = self.simulation_state.get_block_output(node_id, i) {
                            outputs.push(val);
                            node_data.push(val);
                        }
                    }
                    if !node_data.is_empty() {
                        self.plot_state.add_scope_data(node_id.clone(), t, node_data);
                    }
                }
            }
        }

        // If no scope blocks, record outputs from source blocks
        if outputs.is_empty() {
            for node_id in &order {
                if let Some(node) = self.graph_state.graph.get_node(node_id) {
                    if matches!(
                        node.block_type.as_str(),
                        "Constant" | "Sinusoidal" | "Step" | "Ramp" | "Integrator"
                    ) {
                        if let Some(val) = self.simulation_state.get_block_output(node_id, 0) {
                            outputs.push(val);
                        }
                    }
                }
            }
        }

        self.plot_state.add_plot_point(t, outputs);
    }

    /// Clear result history
    pub fn clear_result_history(&mut self) {
        self.plot_state.clear_history();
    }

    // Backward compatibility accessors for direct field access
    // These allow existing UI code to continue working with minimal changes

    /// Access canvas zoom
    pub fn zoom(&self) -> f32 {
        self.canvas_state.zoom
    }

    /// Access canvas zoom mutably
    pub fn zoom_mut(&mut self) -> &mut f32 {
        &mut self.canvas_state.zoom
    }

    /// Access canvas pan
    pub fn pan(&self) -> &Position {
        &self.canvas_state.pan
    }

    /// Access canvas pan mutably
    pub fn pan_mut(&mut self) -> &mut Position {
        &mut self.canvas_state.pan
    }

    /// Access show_canvas_grid flag
    pub fn show_canvas_grid(&self) -> bool {
        self.canvas_state.show_grid
    }

    /// Access show_canvas_grid flag mutably
    pub fn show_canvas_grid_mut(&mut self) -> &mut bool {
        &mut self.canvas_state.show_grid
    }

    /// Access simulation time
    pub fn sim_time(&self) -> f64 {
        self.simulation_state.sim_time
    }

    /// Access simulation settings
    pub fn settings(&self) -> &rustsim_types::SimulationSettings {
        &self.simulation_state.settings
    }

    /// Access simulation settings mutably
    pub fn settings_mut(&mut self) -> &mut rustsim_types::SimulationSettings {
        &mut self.simulation_state.settings
    }

    /// Access compilation status
    pub fn compilation_status(&self) -> &CompilationStatus {
        &self.simulation_state.compilation_status
    }

    /// Access compilation log
    pub fn compilation_log(&self) -> &[String] {
        &self.simulation_state.compilation_log
    }

    /// Clear compilation log
    pub fn clear_compilation_log(&mut self) {
        self.simulation_state.compilation_log.clear();
    }

    /// Add to compilation log
    pub fn push_compilation_log(&mut self, msg: String) {
        self.simulation_state.compilation_log.push(msg);
    }

    /// Access selected nodes
    pub fn selected_nodes(&self) -> &std::collections::HashSet<String> {
        &self.canvas_state.selected_nodes
    }

    /// Access selected nodes mutably
    pub fn selected_nodes_mut(&mut self) -> &mut std::collections::HashSet<String> {
        &mut self.canvas_state.selected_nodes
    }

    /// Access selected connections
    pub fn selected_connections(&self) -> &std::collections::HashSet<String> {
        &self.canvas_state.selected_connections
    }

    /// Access selected connections mutably
    pub fn selected_connections_mut(&mut self) -> &mut std::collections::HashSet<String> {
        &mut self.canvas_state.selected_connections
    }

    /// Access dragging node
    pub fn dragging_node(&self) -> &Option<String> {
        &self.canvas_state.dragging_node
    }

    /// Access dragging node mutably
    pub fn dragging_node_mut(&mut self) -> &mut Option<String> {
        &mut self.canvas_state.dragging_node
    }

    /// Access pending connection
    pub fn pending_connection(&self) -> &Option<super::PendingConnection> {
        &self.canvas_state.pending_connection
    }

    /// Access pending connection mutably
    pub fn pending_connection_mut(&mut self) -> &mut Option<super::PendingConnection> {
        &mut self.canvas_state.pending_connection
    }

    /// Access plot data
    pub fn plot_data(&self) -> &[(f64, Vec<f64>)] {
        &self.plot_state.plot_data
    }

    /// Clear plot data
    pub fn clear_plot_data(&mut self) {
        self.plot_state.clear();
    }

    /// Access plot labels
    pub fn plot_labels(&self) -> &[String] {
        &self.plot_state.plot_labels
    }

    /// Access scope data
    pub fn scope_data(&self) -> &std::collections::HashMap<String, Vec<(f64, Vec<f64>)>> {
        &self.plot_state.scope_data
    }

    /// Access spectrum data
    pub fn spectrum_data(&self) -> &std::collections::HashMap<String, Vec<(f64, Vec<f64>)>> {
        &self.plot_state.spectrum_data
    }

    /// Access plot settings
    pub fn plot_settings(&self) -> &crate::ui::PlotSettings {
        &self.plot_state.plot_settings
    }

    /// Access plot settings mutably
    pub fn plot_settings_mut(&mut self) -> &mut crate::ui::PlotSettings {
        &mut self.plot_state.plot_settings
    }

    /// Access result history
    pub fn result_history(&self) -> &std::collections::VecDeque<super::SimulationResult> {
        &self.plot_state.result_history
    }

    /// Access result history mutably
    pub fn result_history_mut(&mut self) -> &mut std::collections::VecDeque<super::SimulationResult> {
        &mut self.plot_state.result_history
    }

    /// Access block types
    pub fn block_types(&self) -> &[rustsim_types::BlockTypeDefinition] {
        &self.graph_state.block_types
    }

    /// Access use_compiled_mode flag
    #[cfg(not(target_arch = "wasm32"))]
    pub fn use_compiled_mode(&self) -> bool {
        self.simulation_state.use_compiled_mode
    }

    /// Access use_compiled_mode flag mutably
    #[cfg(not(target_arch = "wasm32"))]
    pub fn use_compiled_mode_mut(&mut self) -> &mut bool {
        &mut self.simulation_state.use_compiled_mode
    }

    /// Generate a unique ID
    pub fn generate_id(&mut self) -> String {
        self.graph_state.generate_id()
    }
}

impl Default for AppState {
    fn default() -> Self {
        Self::new()
    }
}
