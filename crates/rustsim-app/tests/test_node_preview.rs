//! Tests for node preview rendering

use rustsim_app::ui::node_preview::*;

#[test]
fn test_decimate_minmax_empty() {
    let data: Vec<(f64, f64)> = vec![];
    let result = decimate_minmax(&data, 100);
    assert_eq!(result.len(), 0);
}

#[test]
fn test_decimate_minmax_small_data() {
    let data = vec![(0.0, 1.0), (1.0, 2.0), (2.0, 3.0)];
    let result = decimate_minmax(&data, 100);
    // Should return all data when less than target
    assert_eq!(result.len(), 3);
    assert_eq!(result, data);
}

#[test]
fn test_decimate_minmax_large_data() {
    let data: Vec<(f64, f64)> = (0..1000)
        .map(|i| (i as f64, (i as f64 * 0.1).sin()))
        .collect();

    let result = decimate_minmax(&data, 100);

    // Should reduce to approximately target points
    assert!(result.len() <= 100);
    assert!(result.len() > 0);
}

#[test]
fn test_decimate_minmax_preserves_range() {
    // Create data with known min/max
    let mut data = vec![];
    for i in 0..100 {
        data.push((i as f64, 5.0));
    }
    data.push((50.0, 10.0)); // Max
    data.push((51.0, 1.0)); // Min

    let result = decimate_minmax(&data, 20);

    // Find min and max in result
    let result_min = result
        .iter()
        .map(|(_, y)| y)
        .fold(f64::INFINITY, |a: f64, &b| a.min(b));
    let result_max = result
        .iter()
        .map(|(_, y)| y)
        .fold(f64::NEG_INFINITY, |a: f64, &b| a.max(b));

    // Should preserve extrema
    assert!((result_min - 1.0).abs() < 1e-6);
    assert!((result_max - 10.0).abs() < 1e-6);
}

#[test]
fn test_decimate_minmax_single_point() {
    let data = vec![(0.0, 1.0)];
    let result = decimate_minmax(&data, 100);
    assert_eq!(result.len(), 1);
    assert_eq!(result[0], (0.0, 1.0));
}

#[test]
fn test_decimate_minmax_two_points() {
    let data = vec![(0.0, 1.0), (1.0, 2.0)];
    let result = decimate_minmax(&data, 100);
    assert_eq!(result.len(), 2);
}

#[test]
fn test_decimate_minmax_monotonic_increasing() {
    let data: Vec<(f64, f64)> = (0..100).map(|i| (i as f64, i as f64)).collect();
    let result = decimate_minmax(&data, 20);

    // Result should still be monotonic increasing
    for i in 1..result.len() {
        assert!(result[i].0 >= result[i - 1].0);
    }
}

#[test]
fn test_decimate_minmax_preserves_spikes() {
    // Create a signal with a sharp spike
    let mut data = vec![];
    for i in 0..50 {
        data.push((i as f64, 1.0));
    }
    data.push((25.0, 100.0)); // Spike
    for i in 50..100 {
        data.push((i as f64, 1.0));
    }

    let result = decimate_minmax(&data, 20);

    // The spike should be preserved
    let max_val = result
        .iter()
        .map(|(_, y)| y)
        .fold(f64::NEG_INFINITY, |a: f64, &b| a.max(b));
    assert!((max_val - 100.0).abs() < 1e-6);
}
