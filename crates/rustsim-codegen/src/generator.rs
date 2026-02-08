//! Rust code generator from simulation graphs.

use rustsim_types::{SimulationGraph, SimulationSettings};
use std::collections::HashMap;

/// Code generator for simulation graphs
pub struct CodeGenerator {
    /// Whether to generate C ABI exports
    pub generate_c_abi: bool,
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
        }
    }

    /// Enable or disable C ABI export generation
    pub fn with_c_abi(mut self, enabled: bool) -> Self {
        self.generate_c_abi = enabled;
        self
    }

    /// Generate Rust code from a simulation graph
    pub fn generate(&self, graph: &SimulationGraph, settings: &SimulationSettings) -> String {
        let mut code = String::new();

        // Header
        code.push_str(&self.generate_header());

        // Compute execution order
        let execution_order = self.topological_sort(graph);

        // Count I/O
        let (input_count, output_count, state_count) = self.count_io(graph);

        // State structure
        code.push_str(&self.generate_state_struct(input_count, output_count, state_count));

        // Block structures
        code.push_str(&self.generate_block_structs(graph));

        // Initialization
        code.push_str(&self.generate_init_function(graph, &execution_order));

        // Step function
        code.push_str(&self.generate_step_function(graph, &execution_order, settings));

        // C ABI exports
        if self.generate_c_abi {
            code.push_str(&self.generate_c_exports(input_count, output_count));
        }

        code
    }

    fn generate_header(&self) -> String {
        r#"//! Auto-generated simulation code
#![allow(unused)]
#![allow(non_snake_case)]

"#
        .to_string()
    }

    fn count_io(&self, graph: &SimulationGraph) -> (usize, usize, usize) {
        let input_count = 0;
        let mut output_count = 0;
        let mut state_count = 0;

        for node in graph.nodes.values() {
            match node.block_type.as_str() {
                "Scope" | "Spectrum" => {
                    // Recording blocks contribute to outputs
                    output_count += node.inputs.len().max(1);
                }
                "Integrator" => {
                    state_count += 1;
                    output_count += 1;
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

    fn generate_state_struct(&self, input_count: usize, output_count: usize, state_count: usize) -> String {
        format!(
            r#"/// Simulation state
#[repr(C)]
pub struct SimState {{
    pub time: f64,
    pub dt: f64,
    pub states: [f64; {state_count}],
    pub inputs: [f64; {input_count}],
    pub outputs: [f64; {output_count}],
}}

impl SimState {{
    pub const fn new() -> Self {{
        Self {{
            time: 0.0,
            dt: 0.01,
            states: [0.0; {state_count}],
            inputs: [0.0; {input_count}],
            outputs: [0.0; {output_count}],
        }}
    }}

    pub fn reset(&mut self) {{
        self.time = 0.0;
        self.states = [0.0; {state_count}];
        self.outputs = [0.0; {output_count}];
    }}
}}

impl Default for SimState {{
    fn default() -> Self {{
        Self::new()
    }}
}}

"#
        )
    }

    fn generate_block_structs(&self, graph: &SimulationGraph) -> String {
        let mut code = String::from("// Block state indices\n");

        // Assign state indices to dynamic blocks
        let mut state_idx = 0;
        for (id, node) in &graph.nodes {
            let safe_id = Self::safe_identifier(id);

            match node.block_type.as_str() {
                "Integrator" => {
                    code.push_str(&format!(
                        "const STATE_IDX_{}: usize = {}; // Integrator state\n",
                        safe_id, state_idx
                    ));
                    state_idx += 1;
                }
                "Differentiator" => {
                    code.push_str(&format!(
                        "const STATE_IDX_{}_PREV: usize = {}; // Previous value\n",
                        safe_id, state_idx
                    ));
                    state_idx += 1;
                }
                _ => {}
            }
        }

        code.push('\n');
        code
    }

    fn generate_init_function(&self, graph: &SimulationGraph, _order: &[String]) -> String {
        let mut code = String::from("/// Initialize simulation\npub fn init(state: &mut SimState) {\n    state.time = 0.0;\n");

        for (id, node) in &graph.nodes {
            let safe_id = Self::safe_identifier(id);

            match node.block_type.as_str() {
                "Integrator" => {
                    let initial = node
                        .get_param("initial_value")
                        .and_then(|v| v.as_f64())
                        .unwrap_or(0.0);
                    code.push_str(&format!(
                        "    state.states[STATE_IDX_{}] = {:.6};\n",
                        safe_id, initial
                    ));
                }
                "Differentiator" => {
                    code.push_str(&format!(
                        "    state.states[STATE_IDX_{}_PREV] = 0.0;\n",
                        safe_id
                    ));
                }
                _ => {}
            }
        }

        code.push_str("}\n\n");
        code
    }

    fn generate_step_function(
        &self,
        graph: &SimulationGraph,
        execution_order: &[String],
        _settings: &SimulationSettings,
    ) -> String {
        let mut code = String::from("/// Perform one simulation step\npub fn step(state: &mut SimState, dt: f64) {\n    state.dt = dt;\n\n");

        // Create a map of node outputs for wiring
        let mut output_map: HashMap<String, String> = HashMap::new();
        let mut output_idx = 0;

        // Process blocks in topological order
        for node_id in execution_order {
            if let Some(node) = graph.get_node(node_id) {
                let safe_id = Self::safe_identifier(node_id);

                code.push_str(&format!("    // {} ({})\n", node.name, node.block_type));

                match node.block_type.as_str() {
                    "Constant" => {
                        let value = node
                            .get_param("value")
                            .and_then(|v| v.as_f64())
                            .unwrap_or(1.0);
                        let var_name = format!("out_{}", safe_id);
                        code.push_str(&format!("    let {} = {:.6}_f64;\n", var_name, value));
                        output_map.insert(node_id.clone(), var_name.clone());
                        code.push_str(&format!("    state.outputs[{}] = {};\n", output_idx, var_name));
                        output_idx += 1;
                    }
                    "Sinusoidal" => {
                        let amp = node.get_param("amplitude").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        let freq = node.get_param("frequency").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        let phase = node.get_param("phase").and_then(|v| v.as_f64()).unwrap_or(0.0);
                        let var_name = format!("out_{}", safe_id);
                        code.push_str(&format!(
                            "let {} = {:.6}_f64 * ({:.6}_f64 * 2.0 * std::f64::consts::PI * state.time + {:.6}_f64).sin();\n",
                            var_name, amp, freq, phase
                        ));
                        output_map.insert(node_id.clone(), var_name.clone());
                        code.push_str(&format!("    state.outputs[{}] = {};\n", output_idx, var_name));
                        output_idx += 1;
                    }
                    "Step" => {
                        let step_time = node.get_param("step_time").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        let initial = node.get_param("initial").and_then(|v| v.as_f64()).unwrap_or(0.0);
                        let final_val = node.get_param("final_value").and_then(|v| v.as_f64()).unwrap_or(1.0);
                        let var_name = format!("out_{}", safe_id);
                        code.push_str(&format!(
                            "let {} = if state.time >= {:.6}_f64 {{ {:.6}_f64 }} else {{ {:.6}_f64 }};\n",
                            var_name, step_time, final_val, initial
                        ));
                        output_map.insert(node_id.clone(), var_name.clone());
                        code.push_str(&format!("    state.outputs[{}] = {};\n", output_idx, var_name));
                        output_idx += 1;
                    }
                    "Integrator" => {
                        // Get input from connected node
                        let input_expr = self.get_input_expression(graph, node_id, 0, &output_map);
                        let var_name = format!("out_{}", safe_id);

                        // Simple Euler integration
                        code.push_str(&format!("    state.states[STATE_IDX_{}] += {} * dt;\n", safe_id, input_expr));
                        code.push_str(&format!("    let {} = state.states[STATE_IDX_{}];\n", var_name, safe_id));
                        output_map.insert(node_id.clone(), var_name.clone());
                        code.push_str(&format!("    state.outputs[{}] = {};\n", output_idx, var_name));
                        output_idx += 1;
                    }
                    "Amplifier" | "Gain" => {
                        let gain = node
                            .get_param("gain")
                            .or_else(|| node.get_param("k"))
                            .and_then(|v| v.as_f64())
                            .unwrap_or(1.0);
                        let input_expr = self.get_input_expression(graph, node_id, 0, &output_map);
                        let var_name = format!("out_{}", safe_id);
                        code.push_str(&format!("    let {} = {} * {:.6}_f64;\n", var_name, input_expr, gain));
                        output_map.insert(node_id.clone(), var_name.clone());
                        code.push_str(&format!("    state.outputs[{}] = {};\n", output_idx, var_name));
                        output_idx += 1;
                    }
                    "Adder" => {
                        let input_count = node.inputs.len().max(2);
                        let mut terms = Vec::new();
                        for i in 0..input_count {
                            terms.push(self.get_input_expression(graph, node_id, i, &output_map));
                        }
                        let var_name = format!("out_{}", safe_id);
                        code.push_str(&format!("    let {} = {};\n", var_name, terms.join(" + ")));
                        output_map.insert(node_id.clone(), var_name.clone());
                        code.push_str(&format!("    state.outputs[{}] = {};\n", output_idx, var_name));
                        output_idx += 1;
                    }
                    "Multiplier" => {
                        let input_count = node.inputs.len().max(2);
                        let mut terms = Vec::new();
                        for i in 0..input_count {
                            terms.push(self.get_input_expression(graph, node_id, i, &output_map));
                        }
                        let var_name = format!("out_{}", safe_id);
                        code.push_str(&format!("    let {} = {};\n", var_name, terms.join(" * ")));
                        output_map.insert(node_id.clone(), var_name.clone());
                        code.push_str(&format!("    state.outputs[{}] = {};\n", output_idx, var_name));
                        output_idx += 1;
                    }
                    "Sin" => {
                        let input_expr = self.get_input_expression(graph, node_id, 0, &output_map);
                        let var_name = format!("out_{}", safe_id);
                        code.push_str(&format!("    let {} = ({}).sin();\n", var_name, input_expr));
                        output_map.insert(node_id.clone(), var_name.clone());
                        code.push_str(&format!("    state.outputs[{}] = {};\n", output_idx, var_name));
                        output_idx += 1;
                    }
                    "Cos" => {
                        let input_expr = self.get_input_expression(graph, node_id, 0, &output_map);
                        let var_name = format!("out_{}", safe_id);
                        code.push_str(&format!("    let {} = ({}).cos();\n", var_name, input_expr));
                        output_map.insert(node_id.clone(), var_name.clone());
                        code.push_str(&format!("    state.outputs[{}] = {};\n", output_idx, var_name));
                        output_idx += 1;
                    }
                    "Abs" => {
                        let input_expr = self.get_input_expression(graph, node_id, 0, &output_map);
                        let var_name = format!("out_{}", safe_id);
                        code.push_str(&format!("    let {} = ({}).abs();\n", var_name, input_expr));
                        output_map.insert(node_id.clone(), var_name.clone());
                        code.push_str(&format!("    state.outputs[{}] = {};\n", output_idx, var_name));
                        output_idx += 1;
                    }
                    "Scope" => {
                        // Scope just records its input
                        let input_expr = self.get_input_expression(graph, node_id, 0, &output_map);
                        code.push_str(&format!("    state.outputs[{}] = {};\n", output_idx, input_expr));
                        output_idx += 1;
                    }
                    _ => {
                        code.push_str(&format!("// TODO: Implement {} block\n", node.block_type));
                    }
                }
                code.push('\n');
            }
        }

        code.push_str("    state.time += dt;\n}\n\n");
        code
    }

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

    fn generate_c_exports(&self, input_count: usize, output_count: usize) -> String {
        format!(
            r#"// C ABI exports for dynamic loading
// Uses opaque pointer to SimState for thread-safety

/// Create a new simulation state (caller owns the pointer)
#[no_mangle]
pub extern "C" fn simulation_create() -> *mut SimState {{
    Box::into_raw(Box::new(SimState::new()))
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
        init(state);
    }}
}}

/// Perform one simulation step
#[no_mangle]
pub extern "C" fn simulation_step(state: *mut SimState, dt: f64) {{
    if let Some(state) = unsafe {{ state.as_mut() }} {{
        step(state, dt);
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
        let mut remaining: Vec<_> = graph.nodes.keys()
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustsim_types::{Connection, NodeInstance, Position};

    #[test]
    fn test_generate_empty_graph() {
        let gen = CodeGenerator::new();
        let graph = SimulationGraph::new();
        let settings = SimulationSettings::default();

        let code = gen.generate(&graph, &settings);

        assert!(code.contains("SimState"));
        assert!(code.contains("simulation_step"));
        assert!(code.contains("#[no_mangle]"));
    }

    #[test]
    fn test_generate_simple_graph() {
        let gen = CodeGenerator::new();
        let mut graph = SimulationGraph::new();
        let settings = SimulationSettings::default();

        // Add a constant block
        let mut node = NodeInstance::new("n1".to_string(), "Constant".to_string(), Position::zero());
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
        graph.add_node(NodeInstance::new("a".to_string(), "Constant".to_string(), Position::zero()));
        graph.add_node(NodeInstance::new("b".to_string(), "Amplifier".to_string(), Position::zero()));
        graph.add_node(NodeInstance::new("c".to_string(), "Scope".to_string(), Position::zero()));

        graph.add_connection(Connection::new("c1".to_string(), "a".to_string(), 0, "b".to_string(), 0));
        graph.add_connection(Connection::new("c2".to_string(), "b".to_string(), 0, "c".to_string(), 0));

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
}
