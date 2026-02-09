//! WASM platform implementation.

use js_sys::{Function, Object, Reflect, Uint8Array};
use wasm_bindgen::prelude::*;
use wasm_bindgen_futures::JsFuture;
use web_sys::{Request, RequestInit, Response};

use crate::{PluginError, PluginResult, SimulationPlugin};

/// WASM simulation plugin loaded from a dynamically instantiated module
pub struct WasmPlugin {
    instance: web_sys::WebAssembly::Instance,
    output_count: usize,
    input_count: usize,
}

impl WasmPlugin {
    /// Compile source on server and load the resulting WASM
    pub async fn compile_and_load(source: &str) -> PluginResult<Self> {
        // Send compilation request to server
        let wasm_bytes = Self::request_compilation(source).await?;

        // Instantiate the received WASM
        Self::instantiate(&wasm_bytes).await
    }

    async fn request_compilation(source: &str) -> PluginResult<Vec<u8>> {
        let window = web_sys::window()
            .ok_or_else(|| PluginError::LoadingFailed("No window object".to_string()))?;

        let mut opts = RequestInit::new();
        opts.method("POST");

        let body = serde_json::json!({ "source": source }).to_string();
        opts.body(Some(&JsValue::from_str(&body)));

        let request = Request::new_with_str_and_init("/api/compile", &opts)
            .map_err(|e| PluginError::CompilationFailed(format!("{:?}", e)))?;

        request
            .headers()
            .set("Content-Type", "application/json")
            .map_err(|e| PluginError::CompilationFailed(format!("{:?}", e)))?;

        let resp_value = JsFuture::from(window.fetch_with_request(&request))
            .await
            .map_err(|e| PluginError::CompilationFailed(format!("{:?}", e)))?;

        let resp: Response = resp_value
            .dyn_into()
            .map_err(|e| PluginError::CompilationFailed(format!("{:?}", e)))?;

        let json = JsFuture::from(
            resp.json()
                .map_err(|e| PluginError::CompilationFailed(format!("{:?}", e)))?,
        )
        .await
        .map_err(|e| PluginError::CompilationFailed(format!("{:?}", e)))?;

        let success = Reflect::get(&json, &"success".into())
            .map_err(|e| PluginError::CompilationFailed(format!("{:?}", e)))?
            .as_bool()
            .unwrap_or(false);

        if !success {
            let error = Reflect::get(&json, &"error".into())
                .ok()
                .and_then(|v| v.as_string())
                .unwrap_or_else(|| "Unknown error".into());
            return Err(PluginError::CompilationFailed(error));
        }

        let wasm_array = Reflect::get(&json, &"wasm".into())
            .map_err(|e| PluginError::CompilationFailed(format!("{:?}", e)))?;

        let uint8_array: Uint8Array = wasm_array
            .dyn_into()
            .map_err(|e| PluginError::CompilationFailed(format!("{:?}", e)))?;

        Ok(uint8_array.to_vec())
    }

    async fn instantiate(wasm_bytes: &[u8]) -> PluginResult<Self> {
        let uint8_array = Uint8Array::from(wasm_bytes.as_ref());

        // Create import object
        let imports = Object::new();

        // Instantiate
        let result = JsFuture::from(web_sys::WebAssembly::instantiate_buffer(
            &uint8_array,
            &imports,
        ))
        .await
        .map_err(|e| PluginError::LoadingFailed(format!("{:?}", e)))?;

        let instance: web_sys::WebAssembly::Instance = Reflect::get(&result, &"instance".into())
            .map_err(|e| PluginError::LoadingFailed(format!("{:?}", e)))?
            .dyn_into()
            .map_err(|e| PluginError::LoadingFailed(format!("{:?}", e)))?;

        // Get counts
        let exports = instance.exports();
        let output_count = Self::call_count_fn(&exports, "get_output_count")?;
        let input_count = Self::call_count_fn(&exports, "get_input_count")?;

        Ok(Self {
            instance,
            output_count,
            input_count,
        })
    }

    fn call_count_fn(exports: &Object, name: &str) -> PluginResult<usize> {
        let func: Function = Reflect::get(exports, &name.into())
            .map_err(|e| PluginError::LoadingFailed(format!("{:?}", e)))?
            .dyn_into()
            .map_err(|e| PluginError::LoadingFailed(format!("{:?}", e)))?;

        let result = func
            .call0(&JsValue::NULL)
            .map_err(|e| PluginError::LoadingFailed(format!("{:?}", e)))?;

        Ok(result.as_f64().unwrap_or(0.0) as usize)
    }

    fn get_export_fn(&self, name: &str) -> PluginResult<Function> {
        let exports = self.instance.exports();
        Reflect::get(&exports, &name.into())
            .map_err(|e| PluginError::LoadingFailed(format!("{:?}", e)))?
            .dyn_into()
            .map_err(|e| PluginError::LoadingFailed(format!("{:?}", e)))
    }
}

impl SimulationPlugin for WasmPlugin {
    fn step(&mut self, dt: f64) -> PluginResult<()> {
        let step_fn = self.get_export_fn("simulation_step")?;
        let result = step_fn
            .call1(&JsValue::NULL, &JsValue::from(dt))
            .map_err(|e| PluginError::StepFailed(format!("{:?}", e)))?;

        let code = result.as_f64().unwrap_or(-1.0) as i32;
        if code == 0 {
            Ok(())
        } else {
            Err(PluginError::StepFailed(format!(
                "Step returned error code: {}",
                code
            )))
        }
    }

    fn get_output(&self, index: usize) -> f64 {
        self.get_export_fn("get_output")
            .and_then(|f| {
                f.call1(&JsValue::NULL, &JsValue::from(index as f64))
                    .map_err(|e| PluginError::LoadingFailed(format!("{:?}", e)))
            })
            .ok()
            .and_then(|v| v.as_f64())
            .unwrap_or(0.0)
    }

    fn set_input(&mut self, index: usize, value: f64) {
        if let Ok(f) = self.get_export_fn("set_input") {
            let _ = f.call2(
                &JsValue::NULL,
                &JsValue::from(index as f64),
                &JsValue::from(value),
            );
        }
    }

    fn reset(&mut self) {
        if let Ok(f) = self.get_export_fn("reset") {
            let _ = f.call0(&JsValue::NULL);
        }
    }

    fn time(&self) -> f64 {
        self.get_export_fn("get_time")
            .and_then(|f| {
                f.call0(&JsValue::NULL)
                    .map_err(|e| PluginError::LoadingFailed(format!("{:?}", e)))
            })
            .ok()
            .and_then(|v| v.as_f64())
            .unwrap_or(0.0)
    }

    fn output_count(&self) -> usize {
        self.output_count
    }

    fn input_count(&self) -> usize {
        self.input_count
    }
}

// SAFETY: WASM plugins are single-threaded by nature
unsafe impl Send for WasmPlugin {}
