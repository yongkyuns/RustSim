//! Dynamic compilation and loading for RustSim compiled mode.
//!
//! This module provides functionality to compile generated Rust code to a shared library
//! and dynamically load it for high-performance simulation execution.

use rustsim_types::SimState;
use std::path::PathBuf;
use thiserror::Error;

#[cfg(not(target_arch = "wasm32"))]
use libloading::{Library, Symbol};
#[cfg(not(target_arch = "wasm32"))]
use std::process::Command;
#[cfg(not(target_arch = "wasm32"))]
use tempfile::TempDir;

/// Errors that can occur during compilation or loading
#[derive(Error, Debug)]
pub enum CompileError {
    #[error("Failed to create temporary directory: {0}")]
    TempDirCreation(#[from] std::io::Error),

    #[error("Cargo build failed: {0}")]
    CargoBuildFailed(String),

    #[error("Failed to find compiled library")]
    LibraryNotFound,

    #[cfg(not(target_arch = "wasm32"))]
    #[error("Failed to load dynamic library: {0}")]
    LibraryLoadFailed(#[from] libloading::Error),

    #[error("Function symbol not found: {0}")]
    SymbolNotFound(String),

    #[error("Compilation not supported on WASM target")]
    WasmNotSupported,
}

/// Compiles generated Rust code to a cdylib and returns the path to the compiled library.
///
/// # Arguments
///
/// * `code` - The generated Rust code to compile
/// * `workspace_root` - Path to the RustSim workspace root (for dependency resolution)
///
/// # Returns
///
/// Returns a tuple of (library_path, temp_dir). The temp_dir must be kept alive
/// as long as the library is in use.
#[cfg(not(target_arch = "wasm32"))]
pub fn compile_to_cdylib(code: &str, workspace_root: &str) -> Result<(PathBuf, TempDir), CompileError> {
    // Create temporary directory for compilation
    let temp_dir = TempDir::new()?;
    let project_dir = temp_dir.path();

    // Create src directory
    let src_dir = project_dir.join("src");
    std::fs::create_dir(&src_dir)?;

    // Write the generated code to src/lib.rs
    let lib_rs_path = src_dir.join("lib.rs");
    std::fs::write(&lib_rs_path, code)?;

    // Create Cargo.toml
    let cargo_toml = format!(
        r#"[package]
name = "rustsim-compiled"
version = "0.1.0"
edition = "2021"

[lib]
crate-type = ["cdylib"]

[dependencies]
rustsim-types = {{ path = "{}/crates/rustsim-types" }}
rustsim-core = {{ path = "{}/crates/rustsim-core" }}

[profile.release]
opt-level = 3
lto = true
codegen-units = 1
"#,
        workspace_root, workspace_root
    );

    let cargo_toml_path = project_dir.join("Cargo.toml");
    std::fs::write(&cargo_toml_path, cargo_toml)?;

    // Run cargo build --release
    let output = Command::new("cargo")
        .arg("build")
        .arg("--release")
        .arg("--lib")
        .current_dir(project_dir)
        .output()?;

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(CompileError::CargoBuildFailed(stderr.to_string()));
    }

    // Find the compiled library
    let target_dir = project_dir.join("target").join("release");
    let lib_name = if cfg!(target_os = "windows") {
        "rustsim_compiled.dll"
    } else if cfg!(target_os = "macos") {
        "librustsim_compiled.dylib"
    } else {
        "librustsim_compiled.so"
    };

    let lib_path = target_dir.join(lib_name);
    if !lib_path.exists() {
        return Err(CompileError::LibraryNotFound);
    }

    Ok((lib_path, temp_dir))
}

#[cfg(target_arch = "wasm32")]
pub fn compile_to_cdylib(_code: &str, _workspace_root: &str) -> Result<(PathBuf, ()), CompileError> {
    Err(CompileError::WasmNotSupported)
}

/// A compiled simulation that can be executed at high performance.
///
/// This structure loads a dynamically compiled simulation library and provides
/// a safe interface to initialize and step through the simulation.
#[cfg(not(target_arch = "wasm32"))]
pub struct CompiledSimulation {
    /// The loaded dynamic library (must be kept alive)
    _library: Library,

    /// Function pointer to sim_init
    sim_init_fn: Symbol<'static, unsafe fn(&mut SimState)>,

    /// Function pointer to sim_reinit (reinitialize blocks from params)
    sim_reinit_fn: Symbol<'static, unsafe fn(&SimState)>,

    /// Function pointer to sim_step
    sim_step_fn: Symbol<'static, unsafe fn(&mut SimState, f64)>,

    /// The simulation state
    state: SimState,

    /// Temporary directory (must be kept alive while library is in use)
    _temp_dir: TempDir,

    /// Parameter mapping for name-based access
    param_mapping: rustsim_codegen::ParamMapping,
}

#[cfg(not(target_arch = "wasm32"))]
impl std::fmt::Debug for CompiledSimulation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CompiledSimulation")
            .field("state", &self.state)
            .finish_non_exhaustive()
    }
}

#[cfg(not(target_arch = "wasm32"))]
impl CompiledSimulation {
    /// Create a new compiled simulation from generated Rust code.
    ///
    /// # Arguments
    ///
    /// * `code` - The generated Rust code containing sim_init and sim_step functions
    /// * `workspace_root` - Path to the RustSim workspace root
    /// * `num_states` - Number of internal state variables
    /// * `num_inputs` - Number of input values
    /// * `num_outputs` - Number of output values
    ///
    /// # Example
    ///
    /// ```no_run
    /// # use rustsim_app::compiler::CompiledSimulation;
    /// let code = r#"
    ///     use rustsim_types::SimState;
    ///
    ///     #[no_mangle]
    ///     pub fn sim_init(state: &mut SimState) {
    ///         state.time = 0.0;
    ///     }
    ///
    ///     #[no_mangle]
    ///     pub fn sim_step(state: &mut SimState, dt: f64) {
    ///         state.time += dt;
    ///     }
    /// "#;
    ///
    /// let mut sim = CompiledSimulation::new(code, ".", 0, 0, 0).unwrap();
    /// sim.init();
    /// sim.step(0.01);
    /// ```
    pub fn new(
        code: &str,
        workspace_root: &str,
        num_states: usize,
        num_inputs: usize,
        num_outputs: usize,
    ) -> Result<Self, CompileError> {
        Self::with_params(code, workspace_root, num_states, num_inputs, num_outputs, rustsim_codegen::ParamMapping::default())
    }

    /// Create a new compiled simulation with parameter mapping.
    pub fn with_params(
        code: &str,
        workspace_root: &str,
        num_states: usize,
        num_inputs: usize,
        num_outputs: usize,
        param_mapping: rustsim_codegen::ParamMapping,
    ) -> Result<Self, CompileError> {
        // Compile the code
        let (lib_path, temp_dir) = compile_to_cdylib(code, workspace_root)?;

        // Load the library
        // SAFETY: We just compiled this library ourselves, so we trust it.
        // The library will be unloaded when this struct is dropped.
        let library = unsafe { Library::new(&lib_path)? };

        // Load the function symbols
        // SAFETY: We're loading symbols that we know exist (sim_init, sim_reinit, and sim_step)
        // with the correct signatures. We extend the lifetime to 'static because
        // we own the Library and will keep it alive for the duration of this struct.
        let sim_init_fn: Symbol<unsafe fn(&mut SimState)> = unsafe {
            library.get(b"sim_init")
                .map_err(|_| CompileError::SymbolNotFound("sim_init".to_string()))?
        };
        let sim_init_fn: Symbol<'static, unsafe fn(&mut SimState)> = unsafe {
            std::mem::transmute(sim_init_fn)
        };

        let sim_reinit_fn: Symbol<unsafe fn(&SimState)> = unsafe {
            library.get(b"sim_reinit")
                .map_err(|_| CompileError::SymbolNotFound("sim_reinit".to_string()))?
        };
        let sim_reinit_fn: Symbol<'static, unsafe fn(&SimState)> = unsafe {
            std::mem::transmute(sim_reinit_fn)
        };

        let sim_step_fn: Symbol<unsafe fn(&mut SimState, f64)> = unsafe {
            library.get(b"sim_step")
                .map_err(|_| CompileError::SymbolNotFound("sim_step".to_string()))?
        };
        let sim_step_fn: Symbol<'static, unsafe fn(&mut SimState, f64)> = unsafe {
            std::mem::transmute(sim_step_fn)
        };

        // Create the simulation state with params
        let num_params = param_mapping.params.len();
        let state = SimState::with_params(num_states, num_inputs, num_outputs, num_params);

        Ok(Self {
            _library: library,
            sim_init_fn,
            sim_reinit_fn,
            sim_step_fn,
            state,
            _temp_dir: temp_dir,
            param_mapping,
        })
    }

    /// Initialize the simulation state.
    ///
    /// This calls the compiled `sim_init` function to set up initial conditions.
    pub fn init(&mut self) {
        // SAFETY: We know the function signature is correct and the state is valid
        unsafe {
            (self.sim_init_fn)(&mut self.state);
        }
    }

    /// Reinitialize blocks with current parameters.
    ///
    /// This recreates all block instances using the current values in state.params.
    /// Call this after changing parameters via set_param or set_param_by_name.
    pub fn reinit(&self) {
        // SAFETY: We know the function signature is correct and the state is valid
        unsafe {
            (self.sim_reinit_fn)(&self.state);
        }
    }

    /// Set a parameter by index.
    ///
    /// After setting parameters, call reinit() to apply changes.
    pub fn set_param_by_index(&mut self, index: usize, value: f64) {
        self.state.set_param(index, value);
    }

    /// Set a parameter by node ID and parameter name.
    ///
    /// After setting parameters, call reinit() to apply changes.
    pub fn set_param(&mut self, node_id: &str, param_name: &str, value: f64) -> bool {
        if let Some(idx) = self.param_mapping.get_index(node_id, param_name) {
            self.state.set_param(idx, value);
            true
        } else {
            false
        }
    }

    /// Get a parameter by index.
    pub fn get_param_by_index(&self, index: usize) -> f64 {
        self.state.get_param(index)
    }

    /// Get a parameter by node ID and parameter name.
    pub fn get_param(&self, node_id: &str, param_name: &str) -> Option<f64> {
        self.param_mapping.get_index(node_id, param_name)
            .map(|idx| self.state.get_param(idx))
    }

    /// Get the parameter mapping.
    pub fn param_mapping(&self) -> &rustsim_codegen::ParamMapping {
        &self.param_mapping
    }

    /// Step the simulation forward by one time step.
    ///
    /// # Arguments
    ///
    /// * `dt` - The time step in seconds
    pub fn step(&mut self, dt: f64) {
        // SAFETY: We know the function signature is correct and the state is valid
        unsafe {
            (self.sim_step_fn)(&mut self.state, dt);
        }
    }

    /// Get an output value by index.
    ///
    /// # Arguments
    ///
    /// * `index` - The output index
    ///
    /// # Returns
    ///
    /// Returns the output value, or 0.0 if the index is out of bounds.
    pub fn get_output(&self, index: usize) -> f64 {
        self.state.outputs.get(index).copied().unwrap_or(0.0)
    }

    /// Set an input value by index.
    ///
    /// # Arguments
    ///
    /// * `index` - The input index
    /// * `value` - The value to set
    pub fn set_input(&mut self, index: usize, value: f64) {
        if let Some(input) = self.state.inputs.get_mut(index) {
            *input = value;
        }
    }

    /// Get the current simulation time.
    pub fn time(&self) -> f64 {
        self.state.time
    }

    /// Get a reference to the simulation state.
    pub fn state(&self) -> &SimState {
        &self.state
    }

    /// Get a mutable reference to the simulation state.
    pub fn state_mut(&mut self) -> &mut SimState {
        &mut self.state
    }

    /// Reset the simulation to initial conditions and re-initialize.
    pub fn reset(&mut self) {
        self.state.reset();
        self.init();
    }
}

#[cfg(target_arch = "wasm32")]
pub struct CompiledSimulation;

#[cfg(target_arch = "wasm32")]
impl CompiledSimulation {
    pub fn new(
        _code: &str,
        _workspace_root: &str,
        _num_states: usize,
        _num_inputs: usize,
        _num_outputs: usize,
    ) -> Result<Self, CompileError> {
        Err(CompileError::WasmNotSupported)
    }
}

#[cfg(test)]
#[cfg(not(target_arch = "wasm32"))]
mod tests {
    use super::*;

    fn get_workspace_root() -> String {
        // Find workspace root by looking for Cargo.toml with [workspace]
        let mut current = std::env::current_dir().unwrap();
        loop {
            let cargo_toml = current.join("Cargo.toml");
            if cargo_toml.exists() {
                if let Ok(contents) = std::fs::read_to_string(&cargo_toml) {
                    if contents.contains("[workspace]") {
                        return current.to_string_lossy().to_string();
                    }
                }
            }
            if !current.pop() {
                panic!("Could not find workspace root");
            }
        }
    }

    #[test]
    fn test_compile_simple_code() {
        let code = r#"
            use rustsim_types::SimState;

            #[no_mangle]
            pub fn sim_init(state: &mut SimState) {
                state.time = 0.0;
            }

            #[no_mangle]
            pub fn sim_step(state: &mut SimState, dt: f64) {
                state.time += dt;
            }
        "#;

        let workspace_root = get_workspace_root();
        let result = compile_to_cdylib(code, &workspace_root);
        assert!(result.is_ok(), "Compilation failed: {:?}", result.err());

        let (lib_path, _temp_dir) = result.unwrap();
        assert!(lib_path.exists(), "Library file does not exist");
    }

    #[test]
    fn test_compiled_simulation_init_and_step() {
        let code = r#"
            use rustsim_types::SimState;

            #[no_mangle]
            pub fn sim_init(state: &mut SimState) {
                state.time = 0.0;
                if !state.states.is_empty() {
                    state.states[0] = 10.0;
                }
            }

            #[no_mangle]
            pub fn sim_step(state: &mut SimState, dt: f64) {
                state.time += dt;
                if !state.states.is_empty() {
                    state.states[0] += 1.0;
                }
                if !state.outputs.is_empty() {
                    state.outputs[0] = state.states[0];
                }
            }
        "#;

        let workspace_root = get_workspace_root();
        let mut sim = CompiledSimulation::new(code, &workspace_root, 1, 0, 1)
            .expect("Failed to create simulation");

        // Test initialization
        sim.init();
        assert_eq!(sim.time(), 0.0);
        assert_eq!(sim.state().states[0], 10.0);

        // Test stepping
        sim.step(0.01);
        assert!((sim.time() - 0.01).abs() < 1e-10);
        assert_eq!(sim.state().states[0], 11.0);
        assert_eq!(sim.get_output(0), 11.0);

        sim.step(0.01);
        assert!((sim.time() - 0.02).abs() < 1e-10);
        assert_eq!(sim.state().states[0], 12.0);
        assert_eq!(sim.get_output(0), 12.0);
    }

    #[test]
    fn test_compiled_simulation_inputs_outputs() {
        let code = r#"
            use rustsim_types::SimState;

            #[no_mangle]
            pub fn sim_init(state: &mut SimState) {
                state.time = 0.0;
            }

            #[no_mangle]
            pub fn sim_step(state: &mut SimState, dt: f64) {
                state.time += dt;
                // Output is 2x the input
                if !state.inputs.is_empty() && !state.outputs.is_empty() {
                    state.outputs[0] = state.inputs[0] * 2.0;
                }
            }
        "#;

        let workspace_root = get_workspace_root();
        let mut sim = CompiledSimulation::new(code, &workspace_root, 0, 1, 1)
            .expect("Failed to create simulation");

        sim.init();
        sim.set_input(0, 5.0);
        sim.step(0.01);

        assert_eq!(sim.get_output(0), 10.0);
    }

    #[test]
    fn test_compilation_error() {
        let bad_code = r#"
            this is not valid rust code
        "#;

        let workspace_root = get_workspace_root();
        let result = compile_to_cdylib(bad_code, &workspace_root);
        assert!(result.is_err());

        match result.unwrap_err() {
            CompileError::CargoBuildFailed(_) => {
                // Expected error
            }
            e => panic!("Expected CargoBuildFailed, got {:?}", e),
        }
    }

    #[test]
    fn test_missing_function_symbol() {
        let code = r#"
            use rustsim_types::SimState;

            // Missing sim_init function

            #[no_mangle]
            pub fn sim_step(state: &mut SimState, dt: f64) {
                state.time += dt;
            }
        "#;

        let workspace_root = get_workspace_root();
        let result = CompiledSimulation::new(code, &workspace_root, 0, 0, 0);
        assert!(result.is_err());

        match result.unwrap_err() {
            CompileError::SymbolNotFound(name) => {
                assert_eq!(name, "sim_init");
            }
            e => panic!("Expected SymbolNotFound, got {:?}", e),
        }
    }
}
