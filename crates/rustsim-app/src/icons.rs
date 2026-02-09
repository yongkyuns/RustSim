//! Block icons using embedded SVGs.
//!
//! This module provides SVG icons for different block types in the simulation.
//! Icons are embedded at compile time and rendered using egui_extras image support.

/// Embedded SVG icons for block types
mod embedded {
    pub const CONSTANT: &[u8] = include_bytes!("../assets/icons/constant.svg");
    pub const SINUSOIDAL: &[u8] = include_bytes!("../assets/icons/sinusoidal.svg");
    pub const STEP: &[u8] = include_bytes!("../assets/icons/step.svg");
    pub const GAIN: &[u8] = include_bytes!("../assets/icons/gain.svg");
    pub const INTEGRATOR: &[u8] = include_bytes!("../assets/icons/integrator.svg");
    pub const ADDER: &[u8] = include_bytes!("../assets/icons/adder.svg");
    pub const SCOPE: &[u8] = include_bytes!("../assets/icons/scope.svg");
    pub const PID: &[u8] = include_bytes!("../assets/icons/pid.svg");
    pub const DERIVATIVE: &[u8] = include_bytes!("../assets/icons/derivative.svg");
    pub const MULTIPLIER: &[u8] = include_bytes!("../assets/icons/multiplier.svg");
    pub const SUBSYSTEM: &[u8] = include_bytes!("../assets/icons/subsystem.svg");
    pub const TRANSFER_FUNCTION: &[u8] = include_bytes!("../assets/icons/transfer_function.svg");
    pub const DEFAULT: &[u8] = include_bytes!("../assets/icons/default.svg");
}

/// Get the URI for an SVG icon (for use with egui's Image widget)
pub fn get_icon_uri(block_type: &str) -> String {
    let normalized = block_type.to_lowercase();
    let icon_name = match normalized.as_str() {
        "sinusoidal" | "sine" => "sinusoidal",
        "gain" | "amplifier" => "gain",
        "adder" | "sum" | "add" => "adder",
        "scope" | "display" => "scope",
        "pid" | "pidcontroller" => "pid",
        "derivative" | "differentiator" => "derivative",
        "multiplier" | "multiply" | "product" => "multiplier",
        "transferfunction" | "transfer_function" | "tf" => "transfer_function",
        "constant" | "step" | "integrator" | "subsystem" => &normalized,
        _ => "default",
    };
    format!("bytes://icon_{}.svg", icon_name)
}

/// Register all icons with egui's image loader system
pub fn register_icons(ctx: &egui::Context) {
    let icons = [
        ("constant", embedded::CONSTANT),
        ("sinusoidal", embedded::SINUSOIDAL),
        ("step", embedded::STEP),
        ("gain", embedded::GAIN),
        ("integrator", embedded::INTEGRATOR),
        ("adder", embedded::ADDER),
        ("scope", embedded::SCOPE),
        ("pid", embedded::PID),
        ("derivative", embedded::DERIVATIVE),
        ("multiplier", embedded::MULTIPLIER),
        ("subsystem", embedded::SUBSYSTEM),
        ("transfer_function", embedded::TRANSFER_FUNCTION),
        ("default", embedded::DEFAULT),
    ];

    for (name, bytes) in icons {
        let uri = format!("bytes://icon_{}.svg", name);
        ctx.include_bytes(uri, bytes);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_icon_uri() {
        assert_eq!(get_icon_uri("Constant"), "bytes://icon_constant.svg");
        assert_eq!(get_icon_uri("Sinusoidal"), "bytes://icon_sinusoidal.svg");
        assert_eq!(get_icon_uri("sine"), "bytes://icon_sinusoidal.svg");
        assert_eq!(get_icon_uri("Amplifier"), "bytes://icon_gain.svg");
        assert_eq!(get_icon_uri("UnknownBlock"), "bytes://icon_default.svg");
    }
}
