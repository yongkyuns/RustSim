//! Platform-specific implementations.

#[cfg(not(target_arch = "wasm32"))]
mod native;

#[cfg(target_arch = "wasm32")]
mod wasm;

#[cfg(not(target_arch = "wasm32"))]
pub use native::*;

#[cfg(target_arch = "wasm32")]
pub use wasm::*;

use crate::{PluginResult, SimulationPlugin};

/// Platform-agnostic plugin loader
pub struct PluginLoader;

impl PluginLoader {
    /// Compile source code and load as a plugin
    #[cfg(not(target_arch = "wasm32"))]
    pub async fn compile_and_load(source: &str) -> PluginResult<Box<dyn SimulationPlugin>> {
        let plugin = NativePlugin::compile_and_load(source).await?;
        Ok(Box::new(plugin))
    }

    /// Compile source code and load as a plugin (WASM version)
    #[cfg(target_arch = "wasm32")]
    pub async fn compile_and_load(source: &str) -> PluginResult<Box<dyn SimulationPlugin>> {
        let plugin = WasmPlugin::compile_and_load(source).await?;
        Ok(Box::new(plugin))
    }
}
