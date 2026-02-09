//! Rust code generator from simulation graphs.

use rustsim_types::{SimulationGraph, SimulationSettings};
use std::collections::HashMap;

/// Information about a runtime parameter
#[derive(Debug, Clone)]
pub struct ParamInfo {
    /// Parameter index in state.params
    pub index: usize,
    /// Node ID this parameter belongs to
    pub node_id: String,
    /// Parameter name (e.g., "value", "gain", "kp")
    pub name: String,
    /// Default value
    pub default_value: f64,
}

/// Parameter mapping for a generated simulation
#[derive(Debug, Clone, Default)]
pub struct ParamMapping {
    /// All parameters in order
    pub params: Vec<ParamInfo>,
    /// Map from "node_id.param_name" to parameter index
    pub name_to_index: HashMap<String, usize>,
}

impl ParamMapping {
    /// Get parameter index by node_id and param_name
    pub fn get_index(&self, node_id: &str, param_name: &str) -> Option<usize> {
        let key = format!("{}.{}", node_id, param_name);
        self.name_to_index.get(&key).copied()
    }

    /// Get all default values as a vector
    pub fn default_values(&self) -> Vec<f64> {
        self.params.iter().map(|p| p.default_value).collect()
    }
}

/// State allocation info for a block (legacy, kept for C ABI)
#[derive(Debug, Clone)]
#[allow(dead_code)]
struct StateAllocation {
    /// Starting index in the states vector
    start_idx: usize,
    /// Number of states used by this block
    count: usize,
}

/// Code generator for simulation graphs
pub struct CodeGenerator {
    /// Whether to generate C ABI exports
    pub generate_c_abi: bool,
    /// Whether to generate as a dynamic library (cdylib)
    pub generate_as_library: bool,
}

impl Default for CodeGenerator {
    fn default() -> Self {
        Self::new()
    }
}

impl CodeGenerator {
    pub fn new() -> Self {
        Self {
            generate_c_abi: false,
            generate_as_library: false,
        }
    }

    /// Enable or disable C ABI export generation
    pub fn with_c_abi(mut self, enabled: bool) -> Self {
        self.generate_c_abi = enabled;
        self
    }

    /// Enable or disable library generation (cdylib)
    pub fn with_library_mode(mut self, enabled: bool) -> Self {
        self.generate_as_library = enabled;
        self
    }

    /// Generate Rust code from a simulation graph
    pub fn generate(&self, graph: &SimulationGraph, settings: &SimulationSettings) -> String {
        self.generate_with_params(graph, settings).0
    }

    /// Generate Rust code and return parameter mapping
    pub fn generate_with_params(&self, graph: &SimulationGraph, settings: &SimulationSettings) -> (String, ParamMapping) {
        let mut code = String::new();

        // Header (includes library attributes if needed)
        code.push_str(&self.generate_header());

        // Compute execution order
        let execution_order = self.topological_sort(graph);

        // Collect all parameters
        let param_mapping = self.collect_parameters(graph, &execution_order);

        // Count I/O
        let (input_count, output_count, state_count) = self.count_io(graph);

        // Block structures
        code.push_str(&self.generate_block_structs(graph));

        // Initialization (now uses state.params)
        code.push_str(&self.generate_init_function_with_params(graph, &execution_order, &param_mapping));

        // Reinit function (recreates blocks from current params)
        code.push_str(&self.generate_reinit_function(graph, &execution_order, &param_mapping));

        // Step function
        code.push_str(&self.generate_step_function(graph, &execution_order, settings));

        // C ABI exports
        if self.generate_c_abi {
            code.push_str(&self.generate_c_exports(input_count, output_count, state_count, param_mapping.params.len()));
        }

        (code, param_mapping)
    }

    /// Collect all parameters from the graph
    fn collect_parameters(&self, graph: &SimulationGraph, order: &[String]) -> ParamMapping {
        let mut mapping = ParamMapping::default();
        let mut idx = 0;

        for node_id in order {
            if let Some(node) = graph.nodes.get(node_id) {
                let params_for_block: Vec<(&str, f64)> = match node.block_type.as_str() {
                    "Constant" => {
                        let value = node.get_param("value").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        vec![("value", value)]
                    }
                    "Sinusoidal" => {
                        let amp = node.get_param("amplitude").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        let freq = node.get_param("frequency").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        let phase = node.get_param("phase").and_then(|v| v.as_f64()).unwrap_or(0.0);
                        vec![("amplitude", amp), ("frequency", freq), ("phase", phase)]
                    }
                    "Step" => {
                        let step_time = node.get_param("step_time").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        let final_val = node.get_param("final_value").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        vec![("final_value", final_val), ("step_time", step_time)]
                    }
                    "Ramp" => {
                        let slope = node.get_param("slope").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        let start_time = node.get_param("start_time").and_then(|v| v.as_f64()).unwrap_or(0.0);
                        vec![("slope", slope), ("start_time", start_time)]
                    }
                    "Amplifier" | "Gain" => {
                        let gain = node.get_param("gain")
                            .or_else(|| node.get_param("k"))
                            .and_then(|v| v.as_f64())
                            .unwrap_or(1.0);
                        vec![("gain", gain)]
                    }
                    "Integrator" => {
                        let initial = node.get_param("initial_value").and_then(|v| v.as_f64()).unwrap_or(0.0);
                        vec![("initial_value", initial)]
                    }
                    "PID" => {
                        let kp = node.get_param("kp").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        let ki = node.get_param("ki").and_then(|v| v.as_f64()).unwrap_or(0.0);
                        let kd = node.get_param("kd").and_then(|v| v.as_f64()).unwrap_or(0.0);
                        vec![("kp", kp), ("ki", ki), ("kd", kd)]
                    }
                    "Differentiator" => {
                        let tau = node.get_param("tau").and_then(|v| v.as_f64()).unwrap_or(0.01);
                        vec![("tau", tau)]
                    }
                    "Saturation" => {
                        let min = node.get_param("min").and_then(|v| v.as_f64()).unwrap_or(-1.0);
                        let max = node.get_param("max").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        vec![("min", min), ("max", max)]
                    }
                    _ => vec![], // No parameters for this block type
                };

                for (name, default_value) in params_for_block {
                    let key = format!("{}.{}", node_id, name);
                    mapping.name_to_index.insert(key, idx);
                    mapping.params.push(ParamInfo {
                        index: idx,
                        node_id: node_id.clone(),
                        name: name.to_string(),
                        default_value,
                    });
                    idx += 1;
                }
            }
        }

        mapping
    }

    fn generate_header(&self) -> String {
        let mut header = String::from("//! Auto-generated simulation code\n");

        if self.generate_as_library {
            header.push_str("#![crate_type = \"cdylib\"]\n");
        }

        header.push_str("#![allow(unused)]\n");
        header.push_str("#![allow(non_snake_case)]\n\n");
        header.push_str("use rustsim_types::SimState;\n");
        header.push_str("use rustsim_core::blocks::*;\n");
        header.push_str("use rustsim_core::Block;\n\n");

        header
    }

    fn count_io(&self, graph: &SimulationGraph) -> (usize, usize, usize) {
        let input_count = 0;
        let mut output_count = 0;
        let mut state_count = 0;

        for node in graph.nodes.values() {
            // Count states for stateful blocks
            state_count += self.get_state_count(&node.block_type, node);

            // Count outputs
            match node.block_type.as_str() {
                "Scope" | "Spectrum" => {
                    // Recording blocks contribute to outputs
                    output_count += node.inputs.len().max(1);
                }
                "Constant" | "Sinusoidal" | "Step" | "Ramp" => {
                    output_count += 1;
                }
                _ => {
                    output_count += node.outputs.len().max(1);
                }
            }
        }

        (input_count.max(1), output_count.max(1), state_count.max(1))
    }

    /// Get the number of state variables needed for a block type
    fn get_state_count(&self, block_type: &str, node: &rustsim_types::NodeInstance) -> usize {
        match block_type {
            "Integrator" => 1,  // 1 state: integrated value
            "PID" => 2,          // 2 states: integral, prev_error
            "Differentiator" => 1, // 1 state: prev_input
            "Delay" => {
                // N states for delay buffer
                node.get_param("buffer_size")
                    .and_then(|v| v.as_u64())
                    .unwrap_or(10) as usize
            }
            _ => 0,  // Algebraic blocks have no state
        }
    }


    fn generate_block_structs(&self, graph: &SimulationGraph) -> String {
        let mut code = String::from("// Block index constants\n");

        // Generate block index constants
        let execution_order = self.topological_sort(graph);
        for (idx, node_id) in execution_order.iter().enumerate() {
            let safe_id = Self::safe_identifier(node_id);
            let node = graph.get_node(node_id).unwrap();
            code.push_str(&format!(
                "const BLK_{}: usize = {};  // {}\n",
                safe_id, idx, node.block_type
            ));
        }

        code.push('\n');
        code
    }

    /// Allocate state indices for all blocks
    #[allow(dead_code)]
    fn allocate_states(&self, graph: &SimulationGraph) -> HashMap<String, StateAllocation> {
        let mut allocations = HashMap::new();
        let mut next_idx = 0;

        // Process blocks in a stable order (sorted by ID for determinism)
        let mut node_ids: Vec<_> = graph.nodes.keys().cloned().collect();
        node_ids.sort();

        for id in node_ids {
            if let Some(node) = graph.nodes.get(&id) {
                let count = self.get_state_count(&node.block_type, node);
                allocations.insert(
                    id.clone(),
                    StateAllocation {
                        start_idx: next_idx,
                        count,
                    },
                );
                next_idx += count;
            }
        }

        allocations
    }


    fn generate_init_function(&self, graph: &SimulationGraph, order: &[String]) -> String {
        // Generate struct to hold all blocks
        let mut code = String::from("/// Simulation blocks\nstruct SimBlocks {\n");

        for node_id in order {
            if let Some(node) = graph.nodes.get(node_id) {
                let safe_id = Self::safe_identifier(node_id).to_lowercase();
                let block_type = self.get_block_type_name(&node.block_type, node);
                code.push_str(&format!("    {}: {},\n", safe_id, block_type));
            }
        }
        code.push_str("}\n\n");

        code.push_str("thread_local! {\n");
        code.push_str("    static BLOCKS: std::cell::RefCell<Option<SimBlocks>> = std::cell::RefCell::new(None);\n");
        code.push_str("}\n\n");

        code.push_str("/// Initialize simulation\n#[no_mangle]\npub fn sim_init(state: &mut SimState) {\n");
        code.push_str("    state.time = 0.0;\n\n");
        code.push_str("    // Create block instances\n");
        code.push_str("    let blocks = SimBlocks {\n");

        // Create each block field
        for node_id in order {
            if let Some(node) = graph.nodes.get(node_id) {
                let safe_id = Self::safe_identifier(node_id).to_lowercase();

                match node.block_type.as_str() {
                    "Constant" => {
                        let value = node
                            .get_param("value")
                            .and_then(|v| v.as_f64())
                            .unwrap_or(1.0);
                        code.push_str(&format!(
                            "        {}: Constant::new({:.6}),\n",
                            safe_id, value
                        ));
                    }
                    "Sinusoidal" => {
                        let amp = node
                            .get_param("amplitude")
                            .and_then(|v| v.as_f64())
                            .unwrap_or(1.0);
                        let freq = node
                            .get_param("frequency")
                            .and_then(|v| v.as_f64())
                            .unwrap_or(1.0);
                        let phase = node
                            .get_param("phase")
                            .and_then(|v| v.as_f64())
                            .unwrap_or(0.0);
                        code.push_str(&format!(
                            "        {}: Sinusoidal::new({:.6}, {:.6}, {:.6}),\n",
                            safe_id, amp, freq, phase
                        ));
                    }
                    "Step" => {
                        let step_time = node
                            .get_param("step_time")
                            .and_then(|v| v.as_f64())
                            .unwrap_or(1.0);
                        let final_val = node
                            .get_param("final_value")
                            .and_then(|v| v.as_f64())
                            .unwrap_or(1.0);
                        code.push_str(&format!(
                            "        {}: Step::new({:.6}, {:.6}),\n",
                            safe_id, final_val, step_time
                        ));
                    }
                    "Ramp" => {
                        let slope = node
                            .get_param("slope")
                            .and_then(|v| v.as_f64())
                            .unwrap_or(1.0);
                        let start_time = node
                            .get_param("start_time")
                            .and_then(|v| v.as_f64())
                            .unwrap_or(0.0);
                        code.push_str(&format!(
                            "        {}: Ramp::new({:.6}, {:.6}),\n",
                            safe_id, slope, start_time
                        ));
                    }
                    "Amplifier" | "Gain" => {
                        let gain = node
                            .get_param("gain")
                            .or_else(|| node.get_param("k"))
                            .and_then(|v| v.as_f64())
                            .unwrap_or(1.0);
                        code.push_str(&format!(
                            "        {}: Amplifier::new({:.6}),\n",
                            safe_id, gain
                        ));
                    }
                    "Integrator" => {
                        let initial = node
                            .get_param("initial_value")
                            .and_then(|v| v.as_f64())
                            .unwrap_or(0.0);
                        code.push_str(&format!(
                            "        {}: Integrator::new({:.6}),\n",
                            safe_id, initial
                        ));
                    }
                    "Adder" => {
                        let input_count = node.inputs.len().max(2);
                        // Determine weights based on port names (+ or -)
                        let weights: Vec<f64> = (0..input_count)
                            .map(|i| {
                                node.inputs
                                    .get(i)
                                    .map(|p| if p.name.contains('-') { -1.0 } else { 1.0 })
                                    .unwrap_or(1.0)
                            })
                            .collect();
                        let weights_str = weights
                            .iter()
                            .map(|w| format!("{:.1}", w))
                            .collect::<Vec<_>>()
                            .join(", ");
                        code.push_str(&format!(
                            "        {}: Adder::<{}>::with_weights([{}]),\n",
                            safe_id, input_count, weights_str
                        ));
                    }
                    "Multiplier" => {
                        let input_count = node.inputs.len().max(2);
                        code.push_str(&format!(
                            "        {}: Multiplier::<{}>::new(),\n",
                            safe_id, input_count
                        ));
                    }
                    "Scope" | "Spectrum" => {
                        let input_count = node.inputs.len().max(1);
                        code.push_str(&format!(
                            "        {}: Scope::<{}, 10000>::new(),\n",
                            safe_id, input_count
                        ));
                    }
                    "PID" => {
                        let kp = node.get_param("kp").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        let ki = node.get_param("ki").and_then(|v| v.as_f64()).unwrap_or(0.0);
                        let kd = node.get_param("kd").and_then(|v| v.as_f64()).unwrap_or(0.0);
                        code.push_str(&format!(
                            "        {}: PID::new({:.6}, {:.6}, {:.6}),\n",
                            safe_id, kp, ki, kd
                        ));
                    }
                    "Differentiator" => {
                        let tau = node.get_param("tau").and_then(|v| v.as_f64()).unwrap_or(0.01);
                        code.push_str(&format!(
                            "        {}: Differentiator::new({:.6}),\n",
                            safe_id, tau
                        ));
                    }
                    "Sin" => {
                        code.push_str(&format!("        {}: Sin::new(),\n", safe_id));
                    }
                    "Cos" => {
                        code.push_str(&format!("        {}: Cos::new(),\n", safe_id));
                    }
                    "Tan" => {
                        code.push_str(&format!("        {}: Tan::new(),\n", safe_id));
                    }
                    "Abs" => {
                        code.push_str(&format!("        {}: Abs::new(),\n", safe_id));
                    }
                    "Sqrt" => {
                        code.push_str(&format!("        {}: Sqrt::new(),\n", safe_id));
                    }
                    "Exp" => {
                        code.push_str(&format!("        {}: Exp::new(),\n", safe_id));
                    }
                    "Log" => {
                        code.push_str(&format!("        {}: Log::new(),\n", safe_id));
                    }
                    "Saturation" => {
                        let min = node.get_param("min").and_then(|v| v.as_f64()).unwrap_or(-1.0);
                        let max = node.get_param("max").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        code.push_str(&format!(
                            "        {}: Saturation::new({:.6}, {:.6}),\n",
                            safe_id, min, max
                        ));
                    }
                    _ => {
                        code.push_str(&format!(
                            "        {}: Constant::new(0.0), // TODO: '{}' not implemented\n",
                            safe_id, node.block_type
                        ));
                    }
                }
            }
        }

        code.push_str("    };\n\n");
        code.push_str("    // Store blocks in thread-local storage\n");
        code.push_str("    BLOCKS.with(|b| *b.borrow_mut() = Some(blocks));\n");
        code.push_str("}\n\n");
        code
    }

    /// Generate init function that reads from state.params
    fn generate_init_function_with_params(&self, graph: &SimulationGraph, order: &[String], param_mapping: &ParamMapping) -> String {
        // Generate struct to hold all blocks
        let mut code = String::from("/// Simulation blocks\nstruct SimBlocks {\n");

        for node_id in order {
            if let Some(node) = graph.nodes.get(node_id) {
                let safe_id = Self::safe_identifier(node_id).to_lowercase();
                let block_type = self.get_block_type_name(&node.block_type, node);
                code.push_str(&format!("    {}: {},\n", safe_id, block_type));
            }
        }
        code.push_str("}\n\n");

        code.push_str("thread_local! {\n");
        code.push_str("    static BLOCKS: std::cell::RefCell<Option<SimBlocks>> = std::cell::RefCell::new(None);\n");
        code.push_str("}\n\n");

        // Generate default params initialization
        code.push_str("/// Initialize simulation with default parameters\n");
        code.push_str("#[no_mangle]\n");
        code.push_str("pub fn sim_init(state: &mut SimState) {\n");
        code.push_str("    state.time = 0.0;\n\n");

        // Initialize params with default values
        code.push_str(&format!("    // Initialize {} parameters with defaults\n", param_mapping.params.len()));
        code.push_str(&format!("    state.params.resize({}, 0.0);\n", param_mapping.params.len()));
        for param in &param_mapping.params {
            code.push_str(&format!("    state.params[{}] = {:.6}; // {}.{}\n",
                param.index, param.default_value, param.node_id, param.name));
        }
        code.push_str("\n");

        // Create blocks
        code.push_str("    sim_reinit(state);\n");
        code.push_str("}\n\n");

        code
    }

    /// Generate reinit function that recreates blocks from current state.params
    fn generate_reinit_function(&self, graph: &SimulationGraph, order: &[String], param_mapping: &ParamMapping) -> String {
        let mut code = String::from("/// Reinitialize blocks with current parameters (preserves time and states)\n");
        code.push_str("#[no_mangle]\n");
        code.push_str("pub fn sim_reinit(state: &SimState) {\n");
        code.push_str("    let p = &state.params;\n\n");
        code.push_str("    // Create block instances from current parameters\n");
        code.push_str("    let blocks = SimBlocks {\n");

        // Create each block field reading from params
        for node_id in order {
            if let Some(node) = graph.nodes.get(node_id) {
                let safe_id = Self::safe_identifier(node_id).to_lowercase();

                match node.block_type.as_str() {
                    "Constant" => {
                        let idx = param_mapping.get_index(node_id, "value").unwrap_or(0);
                        code.push_str(&format!(
                            "        {}: Constant::new(p[{}]),\n",
                            safe_id, idx
                        ));
                    }
                    "Sinusoidal" => {
                        let amp_idx = param_mapping.get_index(node_id, "amplitude").unwrap_or(0);
                        let freq_idx = param_mapping.get_index(node_id, "frequency").unwrap_or(0);
                        let phase_idx = param_mapping.get_index(node_id, "phase").unwrap_or(0);
                        code.push_str(&format!(
                            "        {}: Sinusoidal::new(p[{}], p[{}], p[{}]),\n",
                            safe_id, amp_idx, freq_idx, phase_idx
                        ));
                    }
                    "Step" => {
                        let final_idx = param_mapping.get_index(node_id, "final_value").unwrap_or(0);
                        let time_idx = param_mapping.get_index(node_id, "step_time").unwrap_or(0);
                        code.push_str(&format!(
                            "        {}: Step::new(p[{}], p[{}]),\n",
                            safe_id, final_idx, time_idx
                        ));
                    }
                    "Ramp" => {
                        let slope_idx = param_mapping.get_index(node_id, "slope").unwrap_or(0);
                        let start_idx = param_mapping.get_index(node_id, "start_time").unwrap_or(0);
                        code.push_str(&format!(
                            "        {}: Ramp::new(p[{}], p[{}]),\n",
                            safe_id, slope_idx, start_idx
                        ));
                    }
                    "Amplifier" | "Gain" => {
                        let idx = param_mapping.get_index(node_id, "gain").unwrap_or(0);
                        code.push_str(&format!(
                            "        {}: Amplifier::new(p[{}]),\n",
                            safe_id, idx
                        ));
                    }
                    "Integrator" => {
                        let idx = param_mapping.get_index(node_id, "initial_value").unwrap_or(0);
                        code.push_str(&format!(
                            "        {}: Integrator::new(p[{}]),\n",
                            safe_id, idx
                        ));
                    }
                    "Adder" => {
                        let input_count = node.inputs.len().max(2);
                        let weights: Vec<f64> = (0..input_count)
                            .map(|i| {
                                node.inputs
                                    .get(i)
                                    .map(|p| if p.name.contains('-') { -1.0 } else { 1.0 })
                                    .unwrap_or(1.0)
                            })
                            .collect();
                        let weights_str = weights.iter().map(|w| format!("{:.1}", w)).collect::<Vec<_>>().join(", ");
                        code.push_str(&format!(
                            "        {}: Adder::<{}>::with_weights([{}]),\n",
                            safe_id, input_count, weights_str
                        ));
                    }
                    "Multiplier" => {
                        let input_count = node.inputs.len().max(2);
                        code.push_str(&format!(
                            "        {}: Multiplier::<{}>::new(),\n",
                            safe_id, input_count
                        ));
                    }
                    "Scope" | "Spectrum" => {
                        let input_count = node.inputs.len().max(1);
                        code.push_str(&format!(
                            "        {}: Scope::<{}, 10000>::new(),\n",
                            safe_id, input_count
                        ));
                    }
                    "PID" => {
                        let kp_idx = param_mapping.get_index(node_id, "kp").unwrap_or(0);
                        let ki_idx = param_mapping.get_index(node_id, "ki").unwrap_or(0);
                        let kd_idx = param_mapping.get_index(node_id, "kd").unwrap_or(0);
                        code.push_str(&format!(
                            "        {}: PID::new(p[{}], p[{}], p[{}]),\n",
                            safe_id, kp_idx, ki_idx, kd_idx
                        ));
                    }
                    "Differentiator" => {
                        let idx = param_mapping.get_index(node_id, "tau").unwrap_or(0);
                        code.push_str(&format!(
                            "        {}: Differentiator::new(p[{}]),\n",
                            safe_id, idx
                        ));
                    }
                    "Sin" => code.push_str(&format!("        {}: Sin::new(),\n", safe_id)),
                    "Cos" => code.push_str(&format!("        {}: Cos::new(),\n", safe_id)),
                    "Tan" => code.push_str(&format!("        {}: Tan::new(),\n", safe_id)),
                    "Abs" => code.push_str(&format!("        {}: Abs::new(),\n", safe_id)),
                    "Sqrt" => code.push_str(&format!("        {}: Sqrt::new(),\n", safe_id)),
                    "Exp" => code.push_str(&format!("        {}: Exp::new(),\n", safe_id)),
                    "Log" => code.push_str(&format!("        {}: Log::new(),\n", safe_id)),
                    "Saturation" => {
                        let min_idx = param_mapping.get_index(node_id, "min").unwrap_or(0);
                        let max_idx = param_mapping.get_index(node_id, "max").unwrap_or(0);
                        code.push_str(&format!(
                            "        {}: Saturation::new(p[{}], p[{}]),\n",
                            safe_id, min_idx, max_idx
                        ));
                    }
                    _ => {
                        code.push_str(&format!(
                            "        {}: Constant::new(0.0), // TODO: '{}' not implemented\n",
                            safe_id, node.block_type
                        ));
                    }
                }
            }
        }

        code.push_str("    };\n\n");
        code.push_str("    // Store blocks in thread-local storage\n");
        code.push_str("    BLOCKS.with(|b| *b.borrow_mut() = Some(blocks));\n");
        code.push_str("}\n\n");
        code
    }

    fn generate_step_function(
        &self,
        graph: &SimulationGraph,
        execution_order: &[String],
        _settings: &SimulationSettings,
    ) -> String {
        let mut code = String::from("/// Perform one simulation step\n#[no_mangle]\npub fn sim_step(state: &mut SimState, dt: f64) {\n    state.dt = dt;\n    let t = state.time;\n\n");
        code.push_str("    BLOCKS.with(|blocks_cell| {\n");
        code.push_str("        if let Some(ref mut blocks) = *blocks_cell.borrow_mut() {\n");

        let mut output_idx = 0;

        // Process blocks in topological order
        for node_id in execution_order {
            if let Some(node) = graph.get_node(node_id) {
                let safe_id = Self::safe_identifier(node_id).to_lowercase();

                code.push_str(&format!("            // {} - {}\n", node.name, node.block_type));

                // Wire inputs from connected blocks
                let inputs = graph.get_connections_to(node_id);
                for conn in &inputs {
                    let source_id = Self::safe_identifier(&conn.source_node_id).to_lowercase();
                    let target_port = conn.target_port_index;
                    let source_port = conn.source_port_index;
                    code.push_str(&format!(
                        "            blocks.{}.set_input({}, blocks.{}.get_output({}));\n",
                        safe_id, target_port, source_id, source_port
                    ));
                }

                // Update block (compute outputs from inputs)
                code.push_str(&format!("            blocks.{}.update(t);\n", safe_id));

                // Step dynamic blocks (integrate state)
                if matches!(
                    node.block_type.as_str(),
                    "Integrator" | "PID" | "Differentiator" | "Delay" | "ODE" | "StateSpace"
                ) {
                    code.push_str(&format!("            blocks.{}.step(t, dt);\n", safe_id));
                }

                // Copy outputs to state.outputs for recording
                match node.block_type.as_str() {
                    "Constant" | "Sinusoidal" | "Step" | "Ramp" | "Amplifier" | "Gain" | "Integrator" | "Adder" | "Multiplier" => {
                        code.push_str(&format!(
                            "            state.outputs[{}] = blocks.{}.get_output(0);\n",
                            output_idx, safe_id
                        ));
                        output_idx += 1;
                    }
                    "PID" | "Differentiator" | "Delay" => {
                        code.push_str(&format!(
                            "            state.outputs[{}] = blocks.{}.get_output(0);\n",
                            output_idx, safe_id
                        ));
                        output_idx += 1;
                    }
                    "Sin" | "Cos" | "Tan" | "Abs" | "Sqrt" | "Exp" | "Log" | "Saturation" => {
                        code.push_str(&format!(
                            "            state.outputs[{}] = blocks.{}.get_output(0);\n",
                            output_idx, safe_id
                        ));
                        output_idx += 1;
                    }
                    "Scope" | "Spectrum" => {
                        // Scope stores data internally, copy inputs to outputs for recording
                        let input_count = node.inputs.len().max(1);
                        for i in 0..input_count {
                            code.push_str(&format!(
                                "            state.outputs[{}] = blocks.{}.get_output({});\n",
                                output_idx, safe_id, i
                            ));
                            output_idx += 1;
                        }
                    }
                    _ => {
                        code.push_str(&format!("            // TODO: Block type '{}' not yet implemented\n", node.block_type));
                        code.push_str(&format!(
                            "            state.outputs[{}] = 0.0;\n",
                            output_idx
                        ));
                        output_idx += 1;
                    }
                }
                code.push('\n');
            }
        }

        code.push_str("        }\n");  // End if let Some(ref mut blocks)
        code.push_str("    });  // End BLOCKS.with\n");
        code.push_str("    state.time += dt;\n}\n\n");
        code
    }


    #[allow(dead_code)]
    fn get_input_expression(
        &self,
        graph: &SimulationGraph,
        target_node_id: &str,
        target_port: usize,
        output_map: &HashMap<String, String>,
    ) -> String {
        // Find connection to this input
        for conn in &graph.connections {
            if conn.target_node_id == target_node_id && conn.target_port_index == target_port {
                if let Some(expr) = output_map.get(&conn.source_node_id) {
                    return expr.clone();
                }
            }
        }
        // No connection found, return 0
        "0.0_f64".to_string()
    }

    fn generate_c_exports(&self, input_count: usize, output_count: usize, state_count: usize, param_count: usize) -> String {
        format!(
            r#"// C ABI exports for dynamic loading
// Uses opaque pointer to SimState for thread-safety

/// Create a new simulation state (caller owns the pointer)
#[no_mangle]
pub extern "C" fn simulation_create() -> *mut SimState {{
    Box::into_raw(Box::new(SimState::with_params({state_count}, {input_count}, {output_count}, {param_count})))
}}

/// Destroy a simulation state
#[no_mangle]
pub extern "C" fn simulation_destroy(state: *mut SimState) {{
    if !state.is_null() {{
        unsafe {{ drop(Box::from_raw(state)); }}
    }}
}}

/// Initialize simulation state
#[no_mangle]
pub extern "C" fn simulation_init(state: *mut SimState) {{
    if let Some(state) = unsafe {{ state.as_mut() }} {{
        sim_init(state);
    }}
}}

/// Reinitialize blocks with current parameters (call after changing params)
#[no_mangle]
pub extern "C" fn simulation_reinit(state: *const SimState) {{
    if let Some(state) = unsafe {{ state.as_ref() }} {{
        sim_reinit(state);
    }}
}}

/// Perform one simulation step
#[no_mangle]
pub extern "C" fn simulation_step(state: *mut SimState, dt: f64) {{
    if let Some(state) = unsafe {{ state.as_mut() }} {{
        sim_step(state, dt);
    }}
}}

/// Get an output value
#[no_mangle]
pub extern "C" fn simulation_get_output(state: *const SimState, index: usize) -> f64 {{
    unsafe {{ state.as_ref() }}
        .and_then(|s| s.outputs.get(index).copied())
        .unwrap_or(0.0)
}}

/// Set an input value
#[no_mangle]
pub extern "C" fn simulation_set_input(state: *mut SimState, index: usize, value: f64) {{
    if let Some(state) = unsafe {{ state.as_mut() }} {{
        if let Some(input) = state.inputs.get_mut(index) {{
            *input = value;
        }}
    }}
}}

/// Set a parameter value (call simulation_reinit after changing params)
#[no_mangle]
pub extern "C" fn simulation_set_param(state: *mut SimState, index: usize, value: f64) {{
    if let Some(state) = unsafe {{ state.as_mut() }} {{
        if let Some(param) = state.params.get_mut(index) {{
            *param = value;
        }}
    }}
}}

/// Get a parameter value
#[no_mangle]
pub extern "C" fn simulation_get_param(state: *const SimState, index: usize) -> f64 {{
    unsafe {{ state.as_ref() }}
        .and_then(|s| s.params.get(index).copied())
        .unwrap_or(0.0)
}}

/// Get current simulation time
#[no_mangle]
pub extern "C" fn simulation_get_time(state: *const SimState) -> f64 {{
    unsafe {{ state.as_ref() }}.map(|s| s.time).unwrap_or(0.0)
}}

/// Get output count
#[no_mangle]
pub extern "C" fn simulation_get_output_count() -> usize {{
    {output_count}
}}

/// Get input count
#[no_mangle]
pub extern "C" fn simulation_get_input_count() -> usize {{
    {input_count}
}}

/// Get parameter count
#[no_mangle]
pub extern "C" fn simulation_get_param_count() -> usize {{
    {param_count}
}}
"#
        )
    }

    /// Topological sort of nodes based on connections (stable ordering)
    fn topological_sort(&self, graph: &SimulationGraph) -> Vec<String> {
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

        // If some nodes weren't processed, there's a cycle - add them in sorted order
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

    /// Convert a node ID to a safe Rust identifier
    fn safe_identifier(id: &str) -> String {
        id.replace('-', "_").replace('.', "_").to_uppercase()
    }

    /// Get the Rust type name for a block
    fn get_block_type_name(&self, block_type: &str, node: &rustsim_types::NodeInstance) -> String {
        match block_type {
            "Constant" => "Constant".to_string(),
            "Sinusoidal" => "Sinusoidal".to_string(),
            "Step" => "Step".to_string(),
            "Ramp" => "Ramp".to_string(),
            "Amplifier" | "Gain" => "Amplifier".to_string(),
            "Integrator" => "Integrator".to_string(),
            "Adder" => {
                let input_count = node.inputs.len().max(2);
                format!("Adder<{}>", input_count)
            }
            "Multiplier" => {
                let input_count = node.inputs.len().max(2);
                format!("Multiplier<{}>", input_count)
            }
            "Scope" | "Spectrum" => {
                let input_count = node.inputs.len().max(1);
                format!("Scope<{}, 10000>", input_count)
            }
            "PID" => "PID".to_string(),
            "Differentiator" => "Differentiator".to_string(),
            "Sin" => "Sin".to_string(),
            "Cos" => "Cos".to_string(),
            "Tan" => "Tan".to_string(),
            "Abs" => "Abs".to_string(),
            "Sqrt" => "Sqrt".to_string(),
            "Exp" => "Exp".to_string(),
            "Log" => "Log".to_string(),
            "Saturation" => "Saturation".to_string(),
            _ => "Constant".to_string(), // Placeholder
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustsim_types::{Connection, NodeInstance, Position};

    #[test]
    fn test_generate_empty_graph() {
        let gen = CodeGenerator::new().with_c_abi(true);
        let graph = SimulationGraph::new();
        let settings = SimulationSettings::default();

        let code = gen.generate(&graph, &settings);

        assert!(code.contains("use rustsim_types::SimState"));
        assert!(code.contains("sim_init"));
        assert!(code.contains("sim_step"));
        assert!(code.contains("simulation_step"));
        assert!(code.contains("#[no_mangle]"));
    }

    #[test]
    fn test_generate_simple_graph() {
        let gen = CodeGenerator::new();
        let mut graph = SimulationGraph::new();
        let settings = SimulationSettings::default();

        // Add a constant block
        let mut node =
            NodeInstance::new("n1".to_string(), "Constant".to_string(), Position::zero());
        node.set_param("value", serde_json::json!(42.0));
        graph.add_node(node);

        let code = gen.generate(&graph, &settings);

        assert!(code.contains("42.0"));
        assert!(code.contains("Constant"));
    }

    #[test]
    fn test_topological_sort() {
        let gen = CodeGenerator::new();
        let mut graph = SimulationGraph::new();

        // Create: A -> B -> C
        graph.add_node(NodeInstance::new(
            "a".to_string(),
            "Constant".to_string(),
            Position::zero(),
        ));
        graph.add_node(NodeInstance::new(
            "b".to_string(),
            "Amplifier".to_string(),
            Position::zero(),
        ));
        graph.add_node(NodeInstance::new(
            "c".to_string(),
            "Scope".to_string(),
            Position::zero(),
        ));

        graph.add_connection(Connection::new(
            "c1".to_string(),
            "a".to_string(),
            0,
            "b".to_string(),
            0,
        ));
        graph.add_connection(Connection::new(
            "c2".to_string(),
            "b".to_string(),
            0,
            "c".to_string(),
            0,
        ));

        let order = gen.topological_sort(&graph);

        // A should come before B, B before C
        let pos_a = order.iter().position(|x| x == "a").unwrap();
        let pos_b = order.iter().position(|x| x == "b").unwrap();
        let pos_c = order.iter().position(|x| x == "c").unwrap();

        assert!(pos_a < pos_b);
        assert!(pos_b < pos_c);
    }

    #[test]
    fn test_safe_identifier() {
        assert_eq!(CodeGenerator::safe_identifier("node-1"), "NODE_1");
        assert_eq!(CodeGenerator::safe_identifier("my.block.id"), "MY_BLOCK_ID");
    }

    #[test]
    fn test_library_mode() {
        let gen = CodeGenerator::new()
            .with_library_mode(true)
            .with_c_abi(true);
        let graph = SimulationGraph::new();
        let settings = SimulationSettings::default();

        let code = gen.generate(&graph, &settings);

        assert!(code.contains("#![crate_type = \"cdylib\"]"));
        assert!(code.contains("use rustsim_types::SimState"));
        assert!(code.contains("#[no_mangle]"));
        assert!(code.contains("pub fn sim_init"));
        assert!(code.contains("pub fn sim_step"));
    }

    #[test]
    fn test_generated_functions_use_shared_simstate() {
        let gen = CodeGenerator::new();
        let mut graph = SimulationGraph::new();
        let settings = SimulationSettings::default();

        // Add a simple integrator block
        let mut node = NodeInstance::new(
            "integrator1".to_string(),
            "Integrator".to_string(),
            Position::zero(),
        );
        node.set_param("initial_value", serde_json::json!(1.0));
        graph.add_node(node);

        let code = gen.generate(&graph, &settings);

        // Verify no inline SimState definition
        assert!(!code.contains("pub struct SimState"));
        // Verify it imports from rustsim_types
        assert!(code.contains("use rustsim_types::SimState"));
        // Verify functions are exported with no_mangle
        assert!(code.contains("#[no_mangle]\npub fn sim_init"));
        assert!(code.contains("#[no_mangle]\npub fn sim_step"));
    }
}
