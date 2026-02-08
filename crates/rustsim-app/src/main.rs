//! RustSim - Visual simulation environment

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod app;
mod state;
mod ui;

use app::RustSimApp;

fn main() -> eframe::Result<()> {
    // Log to stdout on native
    #[cfg(not(target_arch = "wasm32"))]
    {
        env_logger::init();
    }

    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1400.0, 900.0])
            .with_min_inner_size([800.0, 600.0])
            .with_title("RustSim"),
        ..Default::default()
    };

    eframe::run_native(
        "RustSim",
        native_options,
        Box::new(|cc| Ok(Box::new(RustSimApp::new(cc)))),
    )
}

// WASM entry point
#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen(start)]
pub async fn start() -> Result<(), JsValue> {
    // Redirect panic messages to console.error
    console_error_panic_hook::set_once();

    let web_options = eframe::WebOptions::default();

    wasm_bindgen_futures::spawn_local(async {
        eframe::WebRunner::new()
            .start(
                "rustsim-canvas",
                web_options,
                Box::new(|cc| Ok(Box::new(RustSimApp::new(cc)))),
            )
            .await
            .expect("Failed to start eframe");
    });

    Ok(())
}
