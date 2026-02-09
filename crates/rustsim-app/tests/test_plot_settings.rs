//! Tests for plot settings and per-signal style configuration

use rustsim_app::ui::{LineStyle, MarkerStyle, PlotSettings, TraceStyle, get_default_color};
use std::collections::HashSet;

#[test]
fn test_default_plot_settings() {
    let settings = PlotSettings::default();
    assert!(!settings.y_axis_log);
    assert!(!settings.x_axis_log);
    assert!(settings.show_legend);
    assert!(settings.show_grid);
    assert_eq!(settings.ghost_trace_count, 0);
    assert!(settings.trace_styles.is_empty());
}

#[test]
fn test_default_trace_style() {
    let style = TraceStyle::default();
    assert_eq!(style.line_style, LineStyle::Solid);
    assert_eq!(style.marker_style, MarkerStyle::None);
    assert_eq!(style.color, None);
    assert!(style.visible);
}

#[test]
fn test_trace_style_with_color() {
    let style = TraceStyle::with_color([255, 0, 0]);
    assert_eq!(style.color, Some([255, 0, 0]));
    assert_eq!(style.line_style, LineStyle::Solid);
    assert!(style.visible);
}

#[test]
fn test_get_or_create_trace_style() {
    let mut settings = PlotSettings::default();

    // Create a new style with default color
    let style = settings.get_or_create_trace_style("node-1-0", Some([255, 0, 0]));
    assert_eq!(style.color, Some([255, 0, 0]));
    assert_eq!(settings.trace_styles.len(), 1);

    // Get the existing style (should not create a new one)
    let style2 = settings.get_or_create_trace_style("node-1-0", None);
    assert_eq!(style2.color, Some([255, 0, 0]));
    assert_eq!(settings.trace_styles.len(), 1);

    // Create another style
    settings.get_or_create_trace_style("node-2-0", Some([0, 255, 0]));
    assert_eq!(settings.trace_styles.len(), 2);
}

#[test]
fn test_cleanup_unused_traces() {
    let mut settings = PlotSettings::default();
    settings.trace_styles.insert("node-1-0".to_string(), TraceStyle::default());
    settings.trace_styles.insert("node-2-0".to_string(), TraceStyle::default());
    settings.trace_styles.insert("node-3-0".to_string(), TraceStyle::default());

    let mut active = HashSet::new();
    active.insert("node-1-0".to_string());
    active.insert("node-2-0".to_string());

    settings.cleanup_unused_traces(&active);
    assert_eq!(settings.trace_styles.len(), 2);
    assert!(settings.trace_styles.contains_key("node-1-0"));
    assert!(settings.trace_styles.contains_key("node-2-0"));
    assert!(!settings.trace_styles.contains_key("node-3-0"));
}

#[test]
fn test_get_trace_style() {
    let mut settings = PlotSettings::default();
    settings.trace_styles.insert("test-trace".to_string(), TraceStyle::default());

    assert!(settings.get_trace_style("test-trace").is_some());
    assert!(settings.get_trace_style("nonexistent").is_none());
}

#[test]
fn test_line_style_all() {
    let styles = LineStyle::all();
    assert_eq!(styles.len(), 4);
    assert!(styles.contains(&LineStyle::Solid));
    assert!(styles.contains(&LineStyle::Dash));
    assert!(styles.contains(&LineStyle::Dot));
    assert!(styles.contains(&LineStyle::None));
}

#[test]
fn test_line_style_names() {
    assert_eq!(LineStyle::Solid.name(), "Solid");
    assert_eq!(LineStyle::Dash.name(), "Dash");
    assert_eq!(LineStyle::Dot.name(), "Dot");
    assert_eq!(LineStyle::None.name(), "None");
}

#[test]
fn test_marker_style_all() {
    let styles = MarkerStyle::all();
    assert_eq!(styles.len(), 4);
    assert!(styles.contains(&MarkerStyle::Circle));
    assert!(styles.contains(&MarkerStyle::Square));
    assert!(styles.contains(&MarkerStyle::Triangle));
    assert!(styles.contains(&MarkerStyle::None));
}

#[test]
fn test_marker_style_names() {
    assert_eq!(MarkerStyle::Circle.name(), "Circle");
    assert_eq!(MarkerStyle::Square.name(), "Square");
    assert_eq!(MarkerStyle::Triangle.name(), "Triangle");
    assert_eq!(MarkerStyle::None.name(), "None");
}

#[test]
fn test_default_color_palette() {
    // Test that we get distinct colors for different indices
    let color0 = get_default_color(0);
    let color1 = get_default_color(1);
    let color9 = get_default_color(9);

    assert_ne!(color0, color1);
    assert_ne!(color0, color9);

    // Test wrapping (color 10 should be same as color 0)
    let color10 = get_default_color(10);
    assert_eq!(color0, color10);

    // Test wrapping at other indices
    let color15 = get_default_color(15);
    let color5 = get_default_color(5);
    assert_eq!(color15, color5);
}

#[test]
fn test_plot_settings_serialization() {
    let mut settings = PlotSettings::default();
    settings.y_axis_log = true;
    settings.x_axis_log = false;
    settings.show_legend = false;
    settings.show_grid = true;
    settings.ghost_trace_count = 3;

    let mut style = TraceStyle::default();
    style.color = Some([255, 0, 0]);
    style.visible = false;
    style.line_style = LineStyle::Dash;
    style.marker_style = MarkerStyle::Circle;
    settings.trace_styles.insert("test-trace".to_string(), style);

    // Serialize
    let json = serde_json::to_string(&settings).expect("Failed to serialize");

    // Deserialize
    let deserialized: PlotSettings = serde_json::from_str(&json).expect("Failed to deserialize");

    assert_eq!(deserialized.y_axis_log, settings.y_axis_log);
    assert_eq!(deserialized.x_axis_log, settings.x_axis_log);
    assert_eq!(deserialized.show_legend, settings.show_legend);
    assert_eq!(deserialized.show_grid, settings.show_grid);
    assert_eq!(deserialized.ghost_trace_count, settings.ghost_trace_count);

    let trace = deserialized.get_trace_style("test-trace").unwrap();
    assert_eq!(trace.color, Some([255, 0, 0]));
    assert!(!trace.visible);
    assert_eq!(trace.line_style, LineStyle::Dash);
    assert_eq!(trace.marker_style, MarkerStyle::Circle);
}

#[test]
fn test_trace_style_serialization() {
    let mut style = TraceStyle::default();
    style.line_style = LineStyle::Dot;
    style.marker_style = MarkerStyle::Square;
    style.color = Some([0, 128, 255]);
    style.visible = false;

    // Serialize
    let json = serde_json::to_string(&style).expect("Failed to serialize");

    // Deserialize
    let deserialized: TraceStyle = serde_json::from_str(&json).expect("Failed to deserialize");

    assert_eq!(deserialized.line_style, LineStyle::Dot);
    assert_eq!(deserialized.marker_style, MarkerStyle::Square);
    assert_eq!(deserialized.color, Some([0, 128, 255]));
    assert!(!deserialized.visible);
}

#[test]
fn test_plot_settings_log_scale_transformation() {
    let mut settings = PlotSettings::default();

    // Test X-axis log scale
    settings.x_axis_log = true;
    settings.y_axis_log = false;
    assert!(settings.x_axis_log);
    assert!(!settings.y_axis_log);

    // Test Y-axis log scale
    settings.x_axis_log = false;
    settings.y_axis_log = true;
    assert!(!settings.x_axis_log);
    assert!(settings.y_axis_log);

    // Test both axes log scale
    settings.x_axis_log = true;
    settings.y_axis_log = true;
    assert!(settings.x_axis_log);
    assert!(settings.y_axis_log);
}

#[test]
fn test_trace_style_visibility_toggle() {
    let mut style = TraceStyle::default();
    assert!(style.visible);

    style.visible = false;
    assert!(!style.visible);

    style.visible = true;
    assert!(style.visible);
}

#[test]
fn test_trace_style_color_reset() {
    let mut style = TraceStyle::with_color([100, 100, 100]);
    assert_eq!(style.color, Some([100, 100, 100]));

    // Reset color to None
    style.color = None;
    assert_eq!(style.color, None);
}

#[test]
fn test_multiple_trace_styles() {
    let mut settings = PlotSettings::default();

    // Create styles for multiple signals
    for i in 0..5 {
        let trace_key = format!("signal-{}", i);
        let color = get_default_color(i);
        settings.get_or_create_trace_style(&trace_key, Some(color));
    }

    assert_eq!(settings.trace_styles.len(), 5);

    // Verify each style has the correct color
    for i in 0..5 {
        let trace_key = format!("signal-{}", i);
        let style = settings.get_trace_style(&trace_key).unwrap();
        let expected_color = get_default_color(i);
        assert_eq!(style.color, Some(expected_color));
    }
}

#[test]
fn test_plot_settings_ghost_trace_count() {
    let mut settings = PlotSettings::default();
    assert_eq!(settings.ghost_trace_count, 0);

    settings.ghost_trace_count = 3;
    assert_eq!(settings.ghost_trace_count, 3);

    settings.ghost_trace_count = 6;
    assert_eq!(settings.ghost_trace_count, 6);
}
