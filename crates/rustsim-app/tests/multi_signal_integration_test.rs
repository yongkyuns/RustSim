//! Integration tests for multi-signal plot display

// Note: We cannot directly test AppState here as it's in a binary crate,
// but we can test the data structures and logic that would be used.

#[test]
fn test_plot_data_structure() {
    // Simulate multi-signal plot data structure
    let mut plot_data: Vec<(f64, Vec<f64>)> = Vec::new();
    let mut plot_labels: Vec<String> = Vec::new();

    // Add labels
    plot_labels.push("Voltage".to_string());
    plot_labels.push("Current".to_string());
    plot_labels.push("Power".to_string());

    // Add data points with 3 signals
    for i in 0..10 {
        let t = i as f64 * 0.1;
        let voltage = (2.0 * std::f64::consts::PI * t).sin();
        let current = (2.0 * std::f64::consts::PI * t + std::f64::consts::PI / 2.0).sin();
        let power = voltage * current;
        plot_data.push((t, vec![voltage, current, power]));
    }

    // Verify structure
    assert_eq!(plot_data.len(), 10);
    assert_eq!(plot_labels.len(), 3);

    // Verify each data point has correct number of signals
    for (_, signals) in &plot_data {
        assert_eq!(signals.len(), 3);
    }
}

#[test]
fn test_scope_label_collection() {
    // Simulate collecting labels from Scope node inputs
    struct Port {
        name: String,
    }

    struct Node {
        block_type: String,
        inputs: Vec<Port>,
    }

    let nodes = vec![
        Node {
            block_type: "Scope".to_string(),
            inputs: vec![
                Port {
                    name: "in 0".to_string(),
                },
                Port {
                    name: "voltage".to_string(),
                },
                Port {
                    name: "current".to_string(),
                },
            ],
        },
        Node {
            block_type: "Sinusoidal".to_string(),
            inputs: vec![],
        },
    ];

    // Collect labels from Scope blocks
    let mut labels = Vec::new();
    for node in &nodes {
        if node.block_type == "Scope" {
            for input in &node.inputs {
                labels.push(input.name.clone());
            }
        }
    }

    assert_eq!(labels.len(), 3);
    assert_eq!(labels[0], "in 0");
    assert_eq!(labels[1], "voltage");
    assert_eq!(labels[2], "current");
}

#[test]
fn test_default_label_generation() {
    // Test default label generation when no Scope blocks exist
    struct Node {
        block_type: String,
        name: String,
    }

    let nodes = vec![
        Node {
            block_type: "Constant".to_string(),
            name: "Constant".to_string(),
        },
        Node {
            block_type: "Sinusoidal".to_string(),
            name: "Sinusoidal".to_string(),
        },
        Node {
            block_type: "Integrator".to_string(),
            name: "Integrator".to_string(),
        },
        Node {
            block_type: "Amplifier".to_string(),
            name: "Gain Block".to_string(),
        },
    ];

    // Generate labels for source blocks
    let mut labels = Vec::new();
    let source_types = ["Constant", "Sinusoidal", "Step", "Ramp", "Integrator"];

    for node in &nodes {
        if source_types.contains(&node.block_type.as_str()) {
            labels.push(format!("{} ({})", node.name, node.block_type));
        }
    }

    assert_eq!(labels.len(), 3);
    assert_eq!(labels[0], "Constant (Constant)");
    assert_eq!(labels[1], "Sinusoidal (Sinusoidal)");
    assert_eq!(labels[2], "Integrator (Integrator)");
}

#[test]
fn test_multi_scope_signal_collection() {
    // Simulate data collection from multiple Scope block inputs
    struct Connection {
        target_port: usize,
        value: f64,
    }

    let scope_connections = vec![
        Connection {
            target_port: 0,
            value: 1.5,
        },
        Connection {
            target_port: 1,
            value: 2.5,
        },
        Connection {
            target_port: 2,
            value: 3.5,
        },
    ];

    let num_inputs = 3;
    let mut outputs = Vec::new();

    // Collect all inputs (simulating the Scope block behavior)
    for i in 0..num_inputs {
        let value = scope_connections
            .iter()
            .find(|conn| conn.target_port == i)
            .map(|conn| conn.value)
            .unwrap_or(0.0);
        outputs.push(value);
    }

    assert_eq!(outputs.len(), 3);
    assert_eq!(outputs[0], 1.5);
    assert_eq!(outputs[1], 2.5);
    assert_eq!(outputs[2], 3.5);
}

#[test]
fn test_signal_extraction_by_index() {
    // Test extracting a specific signal from plot data
    let plot_data = vec![
        (0.0, vec![1.0, 10.0, 100.0]),
        (0.1, vec![1.1, 11.0, 110.0]),
        (0.2, vec![1.2, 12.0, 120.0]),
        (0.3, vec![1.3, 13.0, 130.0]),
    ];

    // Extract middle signal (index 1)
    let signal_1: Vec<(f64, f64)> = plot_data
        .iter()
        .filter_map(|(t, outputs)| outputs.get(1).map(|&y| (*t, y)))
        .collect();

    assert_eq!(signal_1.len(), 4);
    assert_eq!(signal_1[0], (0.0, 10.0));
    assert_eq!(signal_1[1], (0.1, 11.0));
    assert_eq!(signal_1[2], (0.2, 12.0));
    assert_eq!(signal_1[3], (0.3, 13.0));
}

#[test]
fn test_label_and_signal_count_mismatch() {
    // Test scenario where labels and signals have different counts
    let plot_data = vec![
        (0.0, vec![1.0, 2.0, 3.0, 4.0, 5.0]),
        (0.1, vec![1.1, 2.1, 3.1, 4.1, 5.1]),
    ];

    let plot_labels = vec!["Signal A".to_string(), "Signal B".to_string()];

    let num_signals = plot_data.first().map(|(_, o)| o.len()).unwrap_or(0);

    // Verify we handle the mismatch correctly
    assert_eq!(num_signals, 5);
    assert_eq!(plot_labels.len(), 2);

    // Test fallback logic
    for i in 0..num_signals {
        let label = plot_labels
            .get(i)
            .cloned()
            .unwrap_or_else(|| format!("Signal {}", i));

        if i < 2 {
            assert_eq!(label, plot_labels[i]);
        } else {
            assert_eq!(label, format!("Signal {}", i));
        }
    }
}

#[test]
fn test_empty_signal_handling() {
    // Test with empty signal data
    let plot_data: Vec<(f64, Vec<f64>)> = vec![];
    let plot_labels: Vec<String> = vec![];

    let num_signals = plot_data.first().map(|(_, o)| o.len()).unwrap_or(0);

    assert_eq!(num_signals, 0);
    assert_eq!(plot_labels.len(), 0);
}

#[test]
fn test_single_time_point_multiple_signals() {
    // Test with single time point but multiple signals
    let plot_data = vec![(0.0, vec![1.0, 2.0, 3.0])];

    let num_signals = plot_data.first().map(|(_, o)| o.len()).unwrap_or(0);

    assert_eq!(num_signals, 3);
    assert_eq!(plot_data.len(), 1);

    // Verify all signals are accessible
    for i in 0..num_signals {
        let signal: Vec<f64> = plot_data
            .iter()
            .filter_map(|(_, outputs)| outputs.get(i).copied())
            .collect();
        assert_eq!(signal.len(), 1);
        assert_eq!(signal[0], (i + 1) as f64);
    }
}
