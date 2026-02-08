//! Rust code generator from simulation graphs.

use rustsim_types::{Connection, NodeInstance, SimulationGraph, SimulationSettings};
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
            generate_c_abi: true,
        }
    }

    /// Generate Rust code from a simulation graph
    pub fn generate(&self, graph: &SimulationGraph, settings: &SimulationSettings) -> String {
        let mut code = String::new();

        // Header
        code.push_str(&self.generate_header());

        // State structure
        code.push_str(&self.generate_state_struct(graph));

        // Block instantiation
        code.push_str(&self.generate_blocks(graph));

        // Connection wiring
        code.push_str(&self.generate_wiring(graph));

        // Step function
        code.push_str(&self.generate_step_function(graph, settings));

        // C ABI exports
        if self.generate_c_abi {
            code.push_str(&self.generate_c_exports());
        }

        code
    }

    fn generate_header(&self) -> String {
        r#"//! Auto-generated simulation code
#![allow(unused)]
#![allow(non_snake_case)]

use std::collections::HashMap;

"#
        .to_string()
    }

    fn generate_state_struct(&self, graph: &SimulationGraph) -> String {
        let num_states = graph.nodes.len() * 8; // Estimate max states
        let num_io = graph.nodes.len() * 4; // Estimate max I/O

        format!(
            r#"/// Simulation state
#[repr(C)]
pub struct SimState {{
    pub time: f64,
    pub states: [f64; {num_states}],
    pub inputs: [f64; {num_io}],
    pub outputs: [f64; {num_io}],
    state_count: usize,
    input_count: usize,
    output_count: usize,
}}

static mut STATE: SimState = SimState {{
    time: 0.0,
    states: [0.0; {num_states}],
    inputs: [0.0; {num_io}],
    outputs: [0.0; {num_io}],
    state_count: 0,
    input_count: 0,
    output_count: 0,
}};

"#
        )
    }

    fn generate_blocks(&self, graph: &SimulationGraph) -> String {
        let mut code = String::from("// Block definitions\n");

        for (id, node) in &graph.nodes {
            code.push_str(&format!(
                "// Block: {} ({})\n",
                node.name, node.block_type
            ));

            // Generate block-specific code based on type
            match node.block_type.as_str() {
                "Constant" => {
                    let value = node
                        .get_param("value")
                        .and_then(|v| v.as_f64())
                        .unwrap_or(1.0);
                    code.push_str(&format!(
                        "const BLOCK_{}_VALUE: f64 = {:.6};\n",
                        id.to_uppercase().replace("-", "_"),
                        value
                    ));
                }
                "Integrator" => {
                    let initial = node
                        .get_param("initial_value")
                        .and_then(|v| v.as_f64())
                        .unwrap_or(0.0);
                    code.push_str(&format!(
                        "const BLOCK_{}_INITIAL: f64 = {:.6};\n",
                        id.to_uppercase().replace("-", "_"),
                        initial
                    ));
                }
                "Amplifier" | "Gain" => {
                    let gain = node
                        .get_param("gain")
                        .or_else(|| node.get_param("k"))
                        .and_then(|v| v.as_f64())
                        .unwrap_or(1.0);
                    code.push_str(&format!(
                        "const BLOCK_{}_GAIN: f64 = {:.6};\n",
                        id.to_uppercase().replace("-", "_"),
                        gain
                    ));
                }
                _ => {
                    code.push_str(&format!(
                        "// TODO: Generate code for {} block\n",
                        node.block_type
                    ));
                }
            }
        }

        code.push('\n');
        code
    }

    fn generate_wiring(&self, graph: &SimulationGraph) -> String {
        let mut code = String::from("// Connection wiring\n");

        for conn in &graph.connections {
            code.push_str(&format!(
                "// {} [{}] -> {} [{}]\n",
                conn.source_node_id,
                conn.source_port_index,
                conn.target_node_id,
                conn.target_port_index
            ));
        }

        code.push('\n');
        code
    }

    fn generate_step_function(&self, graph: &SimulationGraph, settings: &SimulationSettings) -> String {
        let mut code = String::new();

        code.push_str(&format!(
            r#"/// Perform one simulation step
fn step_internal(dt: f64) -> i32 {{
    unsafe {{
        // TODO: Implement block updates in topological order

        STATE.time += dt;
        0 // Success
    }}
}}

"#
        ));

        code
    }

    fn generate_c_exports(&self) -> String {
        r#"// C ABI exports for dynamic loading

#[no_mangle]
pub extern "C" fn simulation_step(dt: f64) -> i32 {
    step_internal(dt)
}

#[no_mangle]
pub extern "C" fn get_output(index: usize) -> f64 {
    unsafe { STATE.outputs.get(index).copied().unwrap_or(0.0) }
}

#[no_mangle]
pub extern "C" fn set_input(index: usize, value: f64) {
    unsafe {
        if let Some(input) = STATE.inputs.get_mut(index) {
            *input = value;
        }
    }
}

#[no_mangle]
pub extern "C" fn reset() {
    unsafe {
        STATE.time = 0.0;
        STATE.states.fill(0.0);
        STATE.inputs.fill(0.0);
        STATE.outputs.fill(0.0);
    }
}

#[no_mangle]
pub extern "C" fn get_time() -> f64 {
    unsafe { STATE.time }
}

#[no_mangle]
pub extern "C" fn get_state_count() -> usize {
    unsafe { STATE.state_count }
}

#[no_mangle]
pub extern "C" fn get_input_count() -> usize {
    unsafe { STATE.input_count }
}

#[no_mangle]
pub extern "C" fn get_output_count() -> usize {
    unsafe { STATE.output_count }
}
"#
        .to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
