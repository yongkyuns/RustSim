//! Tests for multi-signal plot functionality

use egui::Color32;

/// Color palette for signal traces (copied from plots.rs for testing)
const TRACE_COLORS: [Color32; 9] = [
    Color32::from_rgb(229, 115, 115), // Red
    Color32::from_rgb(129, 199, 132), // Green
    Color32::from_rgb(100, 181, 246), // Blue
    Color32::from_rgb(186, 104, 200), // Purple
    Color32::from_rgb(77, 208, 225),  // Cyan
    Color32::from_rgb(255, 183, 77),  // Orange
    Color32::from_rgb(240, 98, 146),  // Pink
    Color32::from_rgb(77, 182, 172),  // Teal
    Color32::from_rgb(144, 164, 174), // Grey
];

/// Get color for signal at given index
fn get_signal_color(index: usize) -> Color32 {
    TRACE_COLORS[index % TRACE_COLORS.len()]
}

#[test]
fn test_color_assignment() {
    // Test that signal indices map to correct colors
    assert_eq!(get_signal_color(0), TRACE_COLORS[0]); // Red
    assert_eq!(get_signal_color(1), TRACE_COLORS[1]); // Green
    assert_eq!(get_signal_color(2), TRACE_COLORS[2]); // Blue
    assert_eq!(get_signal_color(8), TRACE_COLORS[8]); // Grey

    // Test color wrapping for indices beyond palette size
    assert_eq!(get_signal_color(9), TRACE_COLORS[0]); // Wraps to Red
    assert_eq!(get_signal_color(10), TRACE_COLORS[1]); // Wraps to Green
    assert_eq!(get_signal_color(18), TRACE_COLORS[0]); // Wraps to Red
}

#[test]
fn test_color_palette_size() {
    // Verify we have exactly 9 colors
    assert_eq!(TRACE_COLORS.len(), 9);
}

#[test]
fn test_color_values() {
    // Test specific RGB values for each color
    let red = TRACE_COLORS[0];
    assert_eq!(red.r(), 229);
    assert_eq!(red.g(), 115);
    assert_eq!(red.b(), 115);

    let green = TRACE_COLORS[1];
    assert_eq!(green.r(), 129);
    assert_eq!(green.g(), 199);
    assert_eq!(green.b(), 132);

    let blue = TRACE_COLORS[2];
    assert_eq!(blue.r(), 100);
    assert_eq!(blue.g(), 181);
    assert_eq!(blue.b(), 246);
}

#[test]
fn test_multi_signal_data_extraction() {
    // Simulate multi-signal plot data
    let plot_data = vec![
        (0.0, vec![1.0, 2.0, 3.0]),
        (0.1, vec![1.1, 2.1, 3.1]),
        (0.2, vec![1.2, 2.2, 3.2]),
    ];

    // Extract signal 0
    let signal_0: Vec<f64> = plot_data
        .iter()
        .filter_map(|(_, outputs)| outputs.get(0).copied())
        .collect();
    assert_eq!(signal_0, vec![1.0, 1.1, 1.2]);

    // Extract signal 1
    let signal_1: Vec<f64> = plot_data
        .iter()
        .filter_map(|(_, outputs)| outputs.get(1).copied())
        .collect();
    assert_eq!(signal_1, vec![2.0, 2.1, 2.2]);

    // Extract signal 2
    let signal_2: Vec<f64> = plot_data
        .iter()
        .filter_map(|(_, outputs)| outputs.get(2).copied())
        .collect();
    assert_eq!(signal_2, vec![3.0, 3.1, 3.2]);
}

#[test]
fn test_label_generation() {
    // Test default label generation
    let labels: Vec<String> = (0..5).map(|i| format!("Signal {}", i)).collect();

    assert_eq!(labels[0], "Signal 0");
    assert_eq!(labels[1], "Signal 1");
    assert_eq!(labels[4], "Signal 4");
}

#[test]
fn test_label_from_port_names() {
    // Simulate labels from Scope input ports
    let port_names = vec!["in 0", "in 1", "voltage", "current"];
    let labels: Vec<String> = port_names.iter().map(|s| s.to_string()).collect();

    assert_eq!(labels[0], "in 0");
    assert_eq!(labels[1], "in 1");
    assert_eq!(labels[2], "voltage");
    assert_eq!(labels[3], "current");
}

#[test]
fn test_signal_count_detection() {
    // Test determining number of signals from first data point
    let plot_data = vec![
        (0.0, vec![1.0, 2.0, 3.0, 4.0]),
        (0.1, vec![1.1, 2.1, 3.1, 4.1]),
    ];

    let num_signals = plot_data
        .first()
        .map(|(_, outputs)| outputs.len())
        .unwrap_or(0);

    assert_eq!(num_signals, 4);
}

#[test]
fn test_empty_plot_data() {
    // Test with empty plot data
    let plot_data: Vec<(f64, Vec<f64>)> = Vec::new();

    let num_signals = plot_data
        .first()
        .map(|(_, outputs)| outputs.len())
        .unwrap_or(0);

    assert_eq!(num_signals, 0);
}

#[test]
fn test_single_signal() {
    // Test with single signal
    let plot_data = vec![(0.0, vec![1.0]), (0.1, vec![1.1]), (0.2, vec![1.2])];

    let num_signals = plot_data
        .first()
        .map(|(_, outputs)| outputs.len())
        .unwrap_or(0);

    assert_eq!(num_signals, 1);

    // Verify color assignment for single signal
    assert_eq!(get_signal_color(0), TRACE_COLORS[0]);
}

#[test]
fn test_varying_signal_lengths() {
    // Test handling of varying signal vector lengths (should use filter_map)
    let plot_data = vec![
        (0.0, vec![1.0, 2.0, 3.0]),
        (0.1, vec![1.1, 2.1]), // Missing third signal
        (0.2, vec![1.2, 2.2, 3.2]),
    ];

    // Extract signal 2 - should skip missing values
    let signal_2: Vec<f64> = plot_data
        .iter()
        .filter_map(|(_, outputs)| outputs.get(2).copied())
        .collect();

    assert_eq!(signal_2.len(), 2); // Only 2 values (skipped middle)
    assert_eq!(signal_2, vec![3.0, 3.2]);
}

#[test]
fn test_label_fallback() {
    // Test label fallback when labels vector is shorter than signals
    let labels = vec!["Signal A".to_string(), "Signal B".to_string()];
    let num_signals = 5;

    // Simulate getting label with fallback
    let get_label = |idx: usize| -> String {
        labels
            .get(idx)
            .cloned()
            .unwrap_or_else(|| format!("Signal {}", idx))
    };

    assert_eq!(get_label(0), "Signal A");
    assert_eq!(get_label(1), "Signal B");
    assert_eq!(get_label(2), "Signal 2"); // Fallback
    assert_eq!(get_label(3), "Signal 3"); // Fallback
    assert_eq!(get_label(4), "Signal 4"); // Fallback
}

#[test]
fn test_color_wrapping_with_many_signals() {
    // Test color assignment with more signals than colors
    let num_signals = 25;

    for i in 0..num_signals {
        let color = get_signal_color(i);
        let expected_color = TRACE_COLORS[i % TRACE_COLORS.len()];
        assert_eq!(color, expected_color);
    }
}
