//! Native platform implementation using dynamic libraries.

use std::path::PathBuf;
use std::process::Command;

use libloading::{Library, Symbol};
use tempfile::TempDir;

use crate::{PluginError, PluginResult, SimulationPlugin};

/// Native simulation plugin loaded from a dynamic library
pub struct NativePlugin {
    _lib: Library,
    _temp_dir: Option<TempDir>,
    step_fn: Symbol<'static, unsafe extern "C" fn(f64) -> i32>,
    get_output_fn: Symbol<'static, unsafe extern "C" fn(usize) -> f64>,
    set_input_fn: Symbol<'static, unsafe extern "C" fn(usize, f64)>,
    reset_fn: Symbol<'static, unsafe extern "C" fn()>,
    get_time_fn: Symbol<'static, unsafe extern "C" fn() -> f64>,
    get_output_count_fn: Symbol<'static, unsafe extern "C" fn() -> usize>,
    get_input_count_fn: Symbol<'static, unsafe extern "C" fn() -> usize>,
}

impl NativePlugin {
    /// Compile source code and load as a plugin
    pub async fn compile_and_load(source: &str) -> PluginResult<Self> {
        // Create temp directory
        let temp_dir = TempDir::new()?;
        let src_path = temp_dir.path().join("simulation.rs");
        let lib_path = temp_dir.path().join(Self::lib_name());

        // Write source
        std::fs::write(&src_path, source)?;

        // Compile
        let status = Command::new("rustc")
            .args([
                "--crate-type=cdylib",
                "--edition=2021",
                "-C",
                "opt-level=3",
                "-C",
                "lto=thin",
                src_path.to_str().unwrap(),
                "-o",
                lib_path.to_str().unwrap(),
            ])
            .status()?;

        if !status.success() {
            return Err(PluginError::CompilationFailed(
                "rustc returned non-zero exit code".to_string(),
            ));
        }

        // Load the library
        Self::load_from_path(&lib_path, Some(temp_dir))
    }

    /// Load a pre-compiled library from a path
    pub fn load_from_path(path: &PathBuf, temp_dir: Option<TempDir>) -> PluginResult<Self> {
        unsafe {
            let lib = Library::new(path).map_err(|e| PluginError::LoadingFailed(e.to_string()))?;

            // Get function symbols
            // SAFETY: We keep _lib alive for the lifetime of the struct
            let step_fn: Symbol<unsafe extern "C" fn(f64) -> i32> = std::mem::transmute(
                lib.get::<unsafe extern "C" fn(f64) -> i32>(b"simulation_step")
                    .map_err(|e| PluginError::LoadingFailed(e.to_string()))?,
            );

            let get_output_fn: Symbol<unsafe extern "C" fn(usize) -> f64> = std::mem::transmute(
                lib.get::<unsafe extern "C" fn(usize) -> f64>(b"get_output")
                    .map_err(|e| PluginError::LoadingFailed(e.to_string()))?,
            );

            let set_input_fn: Symbol<unsafe extern "C" fn(usize, f64)> = std::mem::transmute(
                lib.get::<unsafe extern "C" fn(usize, f64)>(b"set_input")
                    .map_err(|e| PluginError::LoadingFailed(e.to_string()))?,
            );

            let reset_fn: Symbol<unsafe extern "C" fn()> = std::mem::transmute(
                lib.get::<unsafe extern "C" fn()>(b"reset")
                    .map_err(|e| PluginError::LoadingFailed(e.to_string()))?,
            );

            let get_time_fn: Symbol<unsafe extern "C" fn() -> f64> = std::mem::transmute(
                lib.get::<unsafe extern "C" fn() -> f64>(b"get_time")
                    .map_err(|e| PluginError::LoadingFailed(e.to_string()))?,
            );

            let get_output_count_fn: Symbol<unsafe extern "C" fn() -> usize> = std::mem::transmute(
                lib.get::<unsafe extern "C" fn() -> usize>(b"get_output_count")
                    .map_err(|e| PluginError::LoadingFailed(e.to_string()))?,
            );

            let get_input_count_fn: Symbol<unsafe extern "C" fn() -> usize> = std::mem::transmute(
                lib.get::<unsafe extern "C" fn() -> usize>(b"get_input_count")
                    .map_err(|e| PluginError::LoadingFailed(e.to_string()))?,
            );

            Ok(Self {
                _lib: lib,
                _temp_dir: temp_dir,
                step_fn,
                get_output_fn,
                set_input_fn,
                reset_fn,
                get_time_fn,
                get_output_count_fn,
                get_input_count_fn,
            })
        }
    }

    fn lib_name() -> &'static str {
        #[cfg(target_os = "linux")]
        {
            "libsimulation.so"
        }
        #[cfg(target_os = "macos")]
        {
            "libsimulation.dylib"
        }
        #[cfg(target_os = "windows")]
        {
            "simulation.dll"
        }
    }
}

impl SimulationPlugin for NativePlugin {
    fn step(&mut self, dt: f64) -> PluginResult<()> {
        let result = unsafe { (self.step_fn)(dt) };
        if result == 0 {
            Ok(())
        } else {
            Err(PluginError::StepFailed(format!(
                "Step returned error code: {}",
                result
            )))
        }
    }

    fn get_output(&self, index: usize) -> f64 {
        unsafe { (self.get_output_fn)(index) }
    }

    fn set_input(&mut self, index: usize, value: f64) {
        unsafe { (self.set_input_fn)(index, value) }
    }

    fn reset(&mut self) {
        unsafe { (self.reset_fn)() }
    }

    fn time(&self) -> f64 {
        unsafe { (self.get_time_fn)() }
    }

    fn output_count(&self) -> usize {
        unsafe { (self.get_output_count_fn)() }
    }

    fn input_count(&self) -> usize {
        unsafe { (self.get_input_count_fn)() }
    }
}

// SAFETY: The plugin can be sent between threads as long as only one thread uses it at a time
unsafe impl Send for NativePlugin {}
