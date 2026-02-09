//! Tests for ghost trace functionality

use rustsim_app::state::{AppState, SimulationResult};

#[test]
fn test_result_history_max_size() {
    let mut state = AppState::new();

    // Create 10 fake simulation results
    for i in 0..10 {
        let result = SimulationResult {
            plot_data: vec![(i as f64, vec![i as f64])],
            plot_labels: vec![format!("Signal {}", i)],
        };
        state.result_history_mut().push_front(result);

        // Trim to max 6
        while state.result_history().len() > 6 {
            state.result_history_mut().pop_back();
        }
    }

    // Should have max 6 results
    assert_eq!(state.result_history().len(), 6);

    // Check that we have the most recent ones (4-9)
    for (idx, result) in state.result_history().iter().enumerate() {
        let expected_value = 9 - idx;
        assert_eq!(result.plot_data[0].0, expected_value as f64);
    }
}

#[test]
fn test_result_history_fifo() {
    let mut state = AppState::new();

    // Add 3 results
    for i in 0..3 {
        let result = SimulationResult {
            plot_data: vec![(i as f64, vec![i as f64])],
            plot_labels: vec![format!("Signal {}", i)],
        };
        state.result_history_mut().push_front(result);
    }

    assert_eq!(state.result_history().len(), 3);

    // Newest should be first (2), oldest should be last (0)
    assert_eq!(state.result_history()[0].plot_data[0].0, 2.0);
    assert_eq!(state.result_history()[2].plot_data[0].0, 0.0);
}

#[test]
fn test_clear_result_history() {
    let mut state = AppState::new();

    // Add some results
    for i in 0..3 {
        let result = SimulationResult {
            plot_data: vec![(i as f64, vec![i as f64])],
            plot_labels: vec![format!("Signal {}", i)],
        };
        state.result_history_mut().push_front(result);
    }

    assert_eq!(state.result_history().len(), 3);

    // Clear history
    state.clear_result_history();

    assert_eq!(state.result_history().len(), 0);
}

#[test]
fn test_ghost_trace_count_bounds() {
    let state = AppState::new();

    // Ghost trace count should default to 0
    assert_eq!(state.plot_settings().ghost_trace_count, 0);

    // Check that it can be set from 0 to 6
    let mut state = state;
    for i in 0..=6 {
        state.plot_settings_mut().ghost_trace_count = i;
        assert_eq!(state.plot_settings().ghost_trace_count, i);
    }
}

#[test]
fn test_opacity_calculation() {
    // Test opacity calculation formula for different ghost counts
    fn calc_opacity(ghost_idx: usize, ghost_count: usize) -> f32 {
        if ghost_count > 1 {
            0.5 - (ghost_idx as f32 * (0.5 - 0.2) / (ghost_count as f32 - 1.0))
        } else {
            0.5
        }
    }

    // Single ghost: should be 0.5
    assert!((calc_opacity(0, 1) - 0.5).abs() < 0.01);

    // Two ghosts: newest=0.5, oldest=0.2
    assert!((calc_opacity(0, 2) - 0.5).abs() < 0.01);
    assert!((calc_opacity(1, 2) - 0.2).abs() < 0.01);

    // Six ghosts: should have linear interpolation
    assert!((calc_opacity(0, 6) - 0.5).abs() < 0.01);
    assert!((calc_opacity(5, 6) - 0.2).abs() < 0.01);
    // Middle values should be in between
    let mid_opacity = calc_opacity(2, 6);
    assert!(mid_opacity > 0.2 && mid_opacity < 0.5);
}

#[test]
fn test_simulation_result_cloning() {
    let result = SimulationResult {
        plot_data: vec![(1.0, vec![1.0, 2.0]), (2.0, vec![3.0, 4.0])],
        plot_labels: vec!["Signal 1".to_string(), "Signal 2".to_string()],
    };

    let cloned = result.clone();

    assert_eq!(cloned.plot_data.len(), 2);
    assert_eq!(cloned.plot_labels.len(), 2);
    assert_eq!(cloned.plot_data[0].0, 1.0);
    assert_eq!(cloned.plot_labels[0], "Signal 1");
}
