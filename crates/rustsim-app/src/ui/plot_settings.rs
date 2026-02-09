//! Plot settings and style configuration.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Plot view mode
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum PlotViewMode {
    Time,
    Spectrum,
}

impl Default for PlotViewMode {
    fn default() -> Self {
        PlotViewMode::Time
    }
}

/// Line style for plot traces
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum LineStyle {
    Solid,
    Dash,
    Dot,
    None,
}

impl Default for LineStyle {
    fn default() -> Self {
        LineStyle::Solid
    }
}

impl LineStyle {
    /// Get all possible line styles
    pub fn all() -> Vec<LineStyle> {
        vec![
            LineStyle::Solid,
            LineStyle::Dash,
            LineStyle::Dot,
            LineStyle::None,
        ]
    }

    /// Get display name
    pub fn name(&self) -> &str {
        match self {
            LineStyle::Solid => "Solid",
            LineStyle::Dash => "Dash",
            LineStyle::Dot => "Dot",
            LineStyle::None => "None",
        }
    }
}

/// Marker style for plot points
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum MarkerStyle {
    Circle,
    Square,
    Triangle,
    None,
}

impl Default for MarkerStyle {
    fn default() -> Self {
        MarkerStyle::None
    }
}

impl MarkerStyle {
    /// Get all possible marker styles
    pub fn all() -> Vec<MarkerStyle> {
        vec![
            MarkerStyle::Circle,
            MarkerStyle::Square,
            MarkerStyle::Triangle,
            MarkerStyle::None,
        ]
    }

    /// Get display name
    pub fn name(&self) -> &str {
        match self {
            MarkerStyle::Circle => "Circle",
            MarkerStyle::Square => "Square",
            MarkerStyle::Triangle => "Triangle",
            MarkerStyle::None => "None",
        }
    }
}

/// Style configuration for a single trace
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TraceStyle {
    pub line_style: LineStyle,
    pub marker_style: MarkerStyle,
    /// RGB color override, None = use default color palette
    pub color: Option<[u8; 3]>,
    pub visible: bool,
}

impl Default for TraceStyle {
    fn default() -> Self {
        Self {
            line_style: LineStyle::Solid,
            marker_style: MarkerStyle::None,
            color: None,
            visible: true,
        }
    }
}

impl TraceStyle {
    /// Create a new trace style with a specific color
    pub fn with_color(color: [u8; 3]) -> Self {
        Self {
            color: Some(color),
            ..Default::default()
        }
    }
}

/// Global plot settings and per-trace style configuration
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PlotSettings {
    /// Per-trace style settings, keyed by "nodeId-signalIndex"
    #[serde(default)]
    pub trace_styles: HashMap<String, TraceStyle>,

    /// Current view mode (Time or Spectrum)
    #[serde(default)]
    pub view_mode: PlotViewMode,

    /// Use logarithmic scale for Y axis
    #[serde(default)]
    pub y_axis_log: bool,

    /// Use logarithmic scale for X axis
    #[serde(default)]
    pub x_axis_log: bool,

    /// Show legend
    #[serde(default = "default_true")]
    pub show_legend: bool,

    /// Show grid
    #[serde(default = "default_true")]
    pub show_grid: bool,

    /// Number of ghost traces to display (0-6)
    #[serde(default)]
    pub ghost_trace_count: usize,
}

fn default_true() -> bool {
    true
}

impl Default for PlotSettings {
    fn default() -> Self {
        Self {
            trace_styles: HashMap::new(),
            view_mode: PlotViewMode::Time,
            y_axis_log: false,
            x_axis_log: false,
            show_legend: true,
            show_grid: true,
            ghost_trace_count: 0,
        }
    }
}

impl PlotSettings {
    /// Get or create a trace style for a given trace key
    pub fn get_or_create_trace_style(
        &mut self,
        trace_key: &str,
        default_color: Option<[u8; 3]>,
    ) -> &mut TraceStyle {
        self.trace_styles
            .entry(trace_key.to_string())
            .or_insert_with(|| {
                if let Some(color) = default_color {
                    TraceStyle::with_color(color)
                } else {
                    TraceStyle::default()
                }
            })
    }

    /// Get a trace style if it exists
    pub fn get_trace_style(&self, trace_key: &str) -> Option<&TraceStyle> {
        self.trace_styles.get(trace_key)
    }

    /// Remove unused trace styles (those not in the provided set of active keys)
    pub fn cleanup_unused_traces(&mut self, active_keys: &std::collections::HashSet<String>) {
        self.trace_styles.retain(|key, _| active_keys.contains(key));
    }
}

/// Default color palette for traces (10 distinct colors)
pub const DEFAULT_COLORS: [[u8; 3]; 10] = [
    [31, 119, 180],  // Blue
    [255, 127, 14],  // Orange
    [44, 160, 44],   // Green
    [214, 39, 40],   // Red
    [148, 103, 189], // Purple
    [140, 86, 75],   // Brown
    [227, 119, 194], // Pink
    [127, 127, 127], // Gray
    [188, 189, 34],  // Yellow-green
    [23, 190, 207],  // Cyan
];

/// Get a default color from the palette based on an index
pub fn get_default_color(index: usize) -> [u8; 3] {
    DEFAULT_COLORS[index % DEFAULT_COLORS.len()]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_plot_settings() {
        let settings = PlotSettings::default();
        assert!(!settings.y_axis_log);
        assert!(!settings.x_axis_log);
        assert!(settings.show_legend);
        assert!(settings.show_grid);
        assert!(settings.trace_styles.is_empty());
    }

    #[test]
    fn test_trace_style_default() {
        let style = TraceStyle::default();
        assert_eq!(style.line_style, LineStyle::Solid);
        assert_eq!(style.marker_style, MarkerStyle::None);
        assert_eq!(style.color, None);
        assert!(style.visible);
    }

    #[test]
    fn test_get_or_create_trace_style() {
        let mut settings = PlotSettings::default();

        // Create a new style with default color
        let style = settings.get_or_create_trace_style("node-1-0", Some([255, 0, 0]));
        assert_eq!(style.color, Some([255, 0, 0]));

        // Get the existing style
        let style2 = settings.get_or_create_trace_style("node-1-0", None);
        assert_eq!(style2.color, Some([255, 0, 0]));
    }

    #[test]
    fn test_cleanup_unused_traces() {
        let mut settings = PlotSettings::default();
        settings
            .trace_styles
            .insert("node-1-0".to_string(), TraceStyle::default());
        settings
            .trace_styles
            .insert("node-2-0".to_string(), TraceStyle::default());
        settings
            .trace_styles
            .insert("node-3-0".to_string(), TraceStyle::default());

        let mut active = std::collections::HashSet::new();
        active.insert("node-1-0".to_string());
        active.insert("node-2-0".to_string());

        settings.cleanup_unused_traces(&active);
        assert_eq!(settings.trace_styles.len(), 2);
        assert!(settings.trace_styles.contains_key("node-1-0"));
        assert!(settings.trace_styles.contains_key("node-2-0"));
        assert!(!settings.trace_styles.contains_key("node-3-0"));
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
    fn test_marker_style_all() {
        let styles = MarkerStyle::all();
        assert_eq!(styles.len(), 4);
        assert!(styles.contains(&MarkerStyle::Circle));
        assert!(styles.contains(&MarkerStyle::Square));
        assert!(styles.contains(&MarkerStyle::Triangle));
        assert!(styles.contains(&MarkerStyle::None));
    }

    #[test]
    fn test_default_color_palette() {
        assert_eq!(DEFAULT_COLORS.len(), 10);

        // Test cycling through colors
        assert_eq!(get_default_color(0), DEFAULT_COLORS[0]);
        assert_eq!(get_default_color(9), DEFAULT_COLORS[9]);
        assert_eq!(get_default_color(10), DEFAULT_COLORS[0]); // Wraps around
        assert_eq!(get_default_color(15), DEFAULT_COLORS[5]);
    }

    #[test]
    fn test_serialization() {
        let mut settings = PlotSettings::default();
        settings.y_axis_log = true;
        settings.x_axis_log = false;
        settings.show_legend = false;
        settings.show_grid = true;

        let mut style = TraceStyle::default();
        style.color = Some([255, 0, 0]);
        style.visible = false;
        settings
            .trace_styles
            .insert("test-trace".to_string(), style);

        let json = serde_json::to_string(&settings).unwrap();
        let deserialized: PlotSettings = serde_json::from_str(&json).unwrap();

        assert_eq!(deserialized.y_axis_log, settings.y_axis_log);
        assert_eq!(deserialized.x_axis_log, settings.x_axis_log);
        assert_eq!(deserialized.show_legend, settings.show_legend);
        assert_eq!(deserialized.show_grid, settings.show_grid);

        let trace = deserialized.get_trace_style("test-trace").unwrap();
        assert_eq!(trace.color, Some([255, 0, 0]));
        assert!(!trace.visible);
    }
}
