//! Simulation state management - contains runtime simulation state and control.

use std::collections::HashMap;
use std::time::{Duration, Instant};

use rustsim_core::block::Block;
use rustsim_core::block_kind::BlockKind;
use rustsim_core::blocks::*;
use rustsim_types::{NodeInstance, SimulationGraph, SimulationSettings};

#[cfg(not(target_arch = "wasm32"))]
use crate::compiler::CompiledSimulation;

#[cfg(not(target_arch = "wasm32"))]
use std::sync::mpsc::Receiver;

/// Status of compilation
#[derive(Clone, Debug, PartialEq)]
pub enum CompilationStatus {
    /// Not compiled yet
    NotCompiled,
    /// Currently compiling
    Compiling,
    /// Compilation succeeded
    Ready,
    /// Compilation failed with error message
    Error(String),
}

/// Message from compilation thread
#[cfg(not(target_arch = "wasm32"))]
pub enum CompilationMessage {
    /// Log message
    Log(String),
    /// Compilation completed successfully
    Success(CompiledSimulation),
    /// Compilation failed
    Error(String),
}

/// State for managing simulation runtime
pub struct SimulationState {
    /// Simulation running state
    running: bool,

    /// Current simulation time
    pub sim_time: f64,

    /// Simulation settings
    pub settings: SimulationSettings,

    /// Whether to use compiled mode instead of interpreter
    pub use_compiled_mode: bool,

    /// Compilation status message
    pub compilation_status: CompilationStatus,

    /// Compilation log messages
    pub compilation_log: Vec<String>,

    /// Receiver for compilation messages (from background thread)
    #[cfg(not(target_arch = "wasm32"))]
    compilation_receiver: Option<Receiver<CompilationMessage>>,

    /// Compiled simulation instance (only on native platforms)
    #[cfg(not(target_arch = "wasm32"))]
    compiled_sim: Option<CompiledSimulation>,

    /// Total time spent in step_simulation calls
    total_step_time: Duration,

    /// Number of steps executed
    step_count: u64,

    /// Block instances for interpreter mode (keyed by node_id)
    blocks: HashMap<String, BlockKind>,

    /// Execution order for interpreter (cached topological sort)
    execution_order: Vec<String>,
}

impl SimulationState {
    pub fn new() -> Self {
        Self {
            running: false,
            sim_time: 0.0,
            settings: SimulationSettings::default(),
            use_compiled_mode: false,
            compilation_status: CompilationStatus::NotCompiled,
            compilation_log: Vec::new(),
            #[cfg(not(target_arch = "wasm32"))]
            compilation_receiver: None,
            #[cfg(not(target_arch = "wasm32"))]
            compiled_sim: None,
            total_step_time: Duration::ZERO,
            step_count: 0,
            blocks: HashMap::new(),
            execution_order: Vec::new(),
        }
    }

    /// Check if simulation is running
    pub fn is_running(&self) -> bool {
        self.running
    }

    /// Start the simulation
    pub fn start(&mut self) {
        self.running = true;
        self.sim_time = 0.0;
        self.total_step_time = Duration::ZERO;
        self.step_count = 0;
    }

    /// Stop the simulation
    pub fn stop(&mut self) {
        self.running = false;
    }

    /// Reset simulation state
    pub fn reset(&mut self) {
        self.running = false;
        self.sim_time = 0.0;
        self.total_step_time = Duration::ZERO;
        self.step_count = 0;
        self.blocks.clear();
        self.execution_order.clear();
    }

    /// Invalidate compiled simulation when graph changes
    pub fn invalidate_compiled_sim(&mut self) {
        #[cfg(not(target_arch = "wasm32"))]
        {
            self.compiled_sim = None;
        }
        if self.compilation_status != CompilationStatus::NotCompiled {
            self.compilation_status = CompilationStatus::NotCompiled;
        }
    }

    /// Get average step execution time
    pub fn average_step_time(&self) -> Option<Duration> {
        if self.step_count > 0 {
            Some(self.total_step_time / self.step_count as u32)
        } else {
            None
        }
    }

    /// Get total step count
    pub fn get_step_count(&self) -> u64 {
        self.step_count
    }

    /// Create block instances from the current graph
    pub fn create_blocks_from_graph(&mut self, graph: &SimulationGraph) {
        self.blocks.clear();
        self.execution_order = Self::get_execution_order(graph);

        for node_id in &self.execution_order.clone() {
            if let Some(node) = graph.get_node(node_id) {
                let block = Self::create_block_from_node(node);
                if let Some(b) = block {
                    self.blocks.insert(node_id.clone(), b);
                }
            }
        }
    }

    /// Create a BlockKind instance from a node
    fn create_block_from_node(node: &NodeInstance) -> Option<BlockKind> {
        let block_type = node.block_type.as_str();

        Some(match block_type {
            "Constant" => {
                let value = node.get_param("value").and_then(|v| v.as_f64()).unwrap_or(1.0);
                BlockKind::Constant(Constant::new(value))
            }
            "Sinusoidal" => {
                let amp = node.get_param("amplitude").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let freq = node.get_param("frequency").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let phase = node.get_param("phase").and_then(|v| v.as_f64()).unwrap_or(0.0);
                BlockKind::Sinusoidal(Sinusoidal::new(amp, freq, phase))
            }
            "Step" => {
                let step_time = node.get_param("step_time").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let final_val = node.get_param("final_value").and_then(|v| v.as_f64()).unwrap_or(1.0);
                BlockKind::Step(Step::new(final_val, step_time))
            }
            "Ramp" => {
                let slope = node.get_param("slope").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let start = node.get_param("start_time").and_then(|v| v.as_f64()).unwrap_or(0.0);
                BlockKind::Ramp(Ramp::new(slope, start))
            }
            "Amplifier" | "Gain" => {
                let gain = node.get_param("gain")
                    .or_else(|| node.get_param("k"))
                    .and_then(|v| v.as_f64())
                    .unwrap_or(1.0);
                BlockKind::Amplifier(Amplifier::new(gain))
            }
            "Integrator" => {
                let initial = node.get_param("initial_value").and_then(|v| v.as_f64()).unwrap_or(0.0);
                BlockKind::Integrator(Integrator::new(initial))
            }
            "Differentiator" => {
                let tau = node.get_param("tau").and_then(|v| v.as_f64()).unwrap_or(0.01);
                BlockKind::Differentiator(Differentiator::new(tau))
            }
            "PID" => {
                let kp = node.get_param("kp").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let ki = node.get_param("ki").and_then(|v| v.as_f64()).unwrap_or(0.0);
                let kd = node.get_param("kd").and_then(|v| v.as_f64()).unwrap_or(0.0);
                BlockKind::PID(PID::new(kp, ki, kd))
            }
            "Adder" => {
                // Determine weights based on port names (+ or -)
                let n = node.inputs.len().max(2);
                let mut weights = vec![1.0; n];
                for (i, port) in node.inputs.iter().enumerate() {
                    if port.name.contains('-') {
                        weights[i] = -1.0;
                    }
                }
                match n {
                    2 => BlockKind::Adder2(Adder::<2>::with_weights([weights[0], weights[1]])),
                    3 => BlockKind::Adder3(Adder::<3>::with_weights([weights[0], weights[1], weights[2]])),
                    _ => BlockKind::Adder4(Adder::<4>::with_weights([weights[0], weights[1], weights[2], weights.get(3).copied().unwrap_or(1.0)])),
                }
            }
            "Multiplier" => {
                let n = node.inputs.len().max(2);
                match n {
                    2 => BlockKind::Multiplier2(Multiplier::<2>::new()),
                    3 => BlockKind::Multiplier3(Multiplier::<3>::new()),
                    _ => BlockKind::Multiplier4(Multiplier::<4>::new()),
                }
            }
            "Saturation" => {
                let min = node.get_param("min").and_then(|v| v.as_f64()).unwrap_or(-1.0);
                let max = node.get_param("max").and_then(|v| v.as_f64()).unwrap_or(1.0);
                BlockKind::Saturation(Saturation::new(min, max))
            }
            "Sin" => BlockKind::Sin(Sin::new()),
            "Cos" => BlockKind::Cos(Cos::new()),
            "Tan" => BlockKind::Tan(Tan::new()),
            "Abs" => BlockKind::Abs(Abs::new()),
            "Sqrt" => BlockKind::Sqrt(Sqrt::new()),
            "Exp" => BlockKind::Exp(Exp::new()),
            "Log" => BlockKind::Log(Log::new()),
            "Sign" => BlockKind::Sign(Sign::new()),
            "Pow" => {
                let exp = node.get_param("exponent").and_then(|v| v.as_f64()).unwrap_or(2.0);
                BlockKind::Pow(Pow::new(exp))
            }
            "Scope" => {
                let n = node.inputs.len().max(1);
                match n {
                    1 => BlockKind::Scope1(Scope::<1, 1000>::new()),
                    2 => BlockKind::Scope2(Scope::<2, 1000>::new()),
                    _ => BlockKind::Scope4(Scope::<4, 1000>::new()),
                }
            }
            "Pulse" => {
                let amp = node.get_param("amplitude").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let period = node.get_param("period").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let width = node.get_param("pulse_width").and_then(|v| v.as_f64()).unwrap_or(0.5);
                let rise = node.get_param("rise_time").and_then(|v| v.as_f64()).unwrap_or(0.0);
                let fall = node.get_param("fall_time").and_then(|v| v.as_f64()).unwrap_or(0.0);
                BlockKind::Pulse(Pulse::new(amp, period, width, rise, fall))
            }
            "SquareWave" => {
                let amp = node.get_param("amplitude").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let freq = node.get_param("frequency").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let duty = node.get_param("duty_cycle").and_then(|v| v.as_f64()).unwrap_or(0.5);
                BlockKind::SquareWave(SquareWave::new(amp, freq, duty))
            }
            "TriangleWave" => {
                let amp = node.get_param("amplitude").and_then(|v| v.as_f64()).unwrap_or(1.0);
                let freq = node.get_param("frequency").and_then(|v| v.as_f64()).unwrap_or(1.0);
                BlockKind::TriangleWave(TriangleWave::new(amp, freq))
            }
            "Comparator" => {
                let threshold = node.get_param("threshold").and_then(|v| v.as_f64()).unwrap_or(0.0);
                BlockKind::Comparator(Comparator::new(threshold))
            }
            "Relay" => {
                let off_thresh = node.get_param("off_threshold").and_then(|v| v.as_f64()).unwrap_or(0.0);
                let on_thresh = node.get_param("on_threshold").and_then(|v| v.as_f64()).unwrap_or(0.5);
                BlockKind::Relay(Relay::new(off_thresh, on_thresh))
            }
            "LowpassRC" => {
                let cutoff = node.get_param("cutoff").and_then(|v| v.as_f64()).unwrap_or(10.0);
                BlockKind::LowpassRC(LowpassRC::new(cutoff))
            }
            "HighpassRC" => {
                let cutoff = node.get_param("cutoff").and_then(|v| v.as_f64()).unwrap_or(10.0);
                BlockKind::HighpassRC(HighpassRC::new(cutoff))
            }
            "RateLimiter" => {
                let rate = node.get_param("rate").and_then(|v| v.as_f64()).unwrap_or(1.0);
                BlockKind::RateLimiter(RateLimiter::new(rate))
            }
            "Min" => BlockKind::Min2(Min::<2>::new()),
            "Max" => BlockKind::Max2(Max::<2>::new()),
            _ => return None, // Unknown block type
        })
    }

    /// Step the simulation forward by dt
    pub fn step(&mut self, graph: &SimulationGraph) {
        let step_start = Instant::now();
        let dt = self.settings.dt;
        let t = self.sim_time;

        // Use compiled simulation if available
        #[cfg(not(target_arch = "wasm32"))]
        if self.use_compiled_mode {
            if let Some(ref mut sim) = self.compiled_sim {
                sim.step(dt);
                self.sim_time += dt;
                self.total_step_time += step_start.elapsed();
                self.step_count += 1;
                return;
            } else {
                // No compiled sim available, fall back to interpreter
                self.use_compiled_mode = false;
            }
        }

        // Interpreter mode using real block instances
        let order = self.execution_order.clone();
        for node_id in &order {
            // Wire inputs from connected source blocks
            let connections: Vec<_> = graph.get_connections_to(node_id)
                .iter()
                .map(|c| (c.source_node_id.clone(), c.source_port_index, c.target_port_index))
                .collect();

            for (source_id, source_port, target_port) in connections {
                if let Some(source_block) = self.blocks.get(&source_id) {
                    let value = source_block.get_output(source_port);
                    if let Some(target_block) = self.blocks.get_mut(node_id) {
                        target_block.set_input(target_port, value);
                    }
                }
            }

            // Update and step the block
            if let Some(block) = self.blocks.get_mut(node_id) {
                block.update(t);
                block.step(t, dt);
            }
        }

        self.sim_time += dt;
        self.total_step_time += step_start.elapsed();
        self.step_count += 1;
    }

    /// Get block output for a specific node and port
    pub fn get_block_output(&self, node_id: &str, port: usize) -> Option<f64> {
        self.blocks.get(node_id).map(|b| b.get_output(port))
    }

    /// Get execution order
    pub fn execution_order(&self) -> &[String] {
        &self.execution_order
    }

    /// Get execution order using topological sort (stable ordering)
    fn get_execution_order(graph: &SimulationGraph) -> Vec<String> {
        let mut in_degree: HashMap<String, usize> = HashMap::new();
        let mut adjacency: HashMap<String, Vec<String>> = HashMap::new();

        // Initialize
        for id in graph.nodes.keys() {
            in_degree.insert(id.clone(), 0);
            adjacency.insert(id.clone(), Vec::new());
        }

        // Build graph
        for conn in &graph.connections {
            if let Some(adj) = adjacency.get_mut(&conn.source_node_id) {
                adj.push(conn.target_node_id.clone());
            }
            if let Some(deg) = in_degree.get_mut(&conn.target_node_id) {
                *deg += 1;
            }
        }

        // Kahn's algorithm with stable ordering (sort by ID for determinism)
        let mut ready: Vec<String> = in_degree
            .iter()
            .filter(|(_, &deg)| deg == 0)
            .map(|(id, _)| id.clone())
            .collect();
        ready.sort(); // Alphabetical order for stability

        let mut result = Vec::new();
        while let Some(node) = ready.pop() {
            result.push(node.clone());

            if let Some(neighbors) = adjacency.get(&node) {
                let mut new_ready = Vec::new();
                for neighbor in neighbors {
                    if let Some(deg) = in_degree.get_mut(neighbor) {
                        *deg -= 1;
                        if *deg == 0 {
                            new_ready.push(neighbor.clone());
                        }
                    }
                }
                // Sort new ready nodes and add them
                new_ready.sort();
                new_ready.reverse(); // Reverse so pop() gives smallest first
                ready.extend(new_ready);
            }
        }

        // Add any remaining nodes in sorted order (in case of cycles)
        let mut remaining: Vec<_> = graph
            .nodes
            .keys()
            .filter(|id| !result.contains(id))
            .cloned()
            .collect();
        remaining.sort();
        result.extend(remaining);

        result
    }

    // Compilation methods (only on native platforms)

    #[cfg(not(target_arch = "wasm32"))]
    pub fn compiled_sim(&self) -> Option<&CompiledSimulation> {
        self.compiled_sim.as_ref()
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub fn compiled_sim_mut(&mut self) -> Option<&mut CompiledSimulation> {
        self.compiled_sim.as_mut()
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub fn set_compiled_sim(&mut self, sim: CompiledSimulation) {
        self.compiled_sim = Some(sim);
        self.compilation_status = CompilationStatus::Ready;
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub fn take_compilation_receiver(&mut self) -> Option<Receiver<CompilationMessage>> {
        self.compilation_receiver.take()
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub fn set_compilation_receiver(&mut self, receiver: Receiver<CompilationMessage>) {
        self.compilation_receiver = Some(receiver);
    }
}

impl Default for SimulationState {
    fn default() -> Self {
        Self::new()
    }
}
