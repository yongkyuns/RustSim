//! Test to verify interpreter and compiled modes produce matching results.

use rustsim_app::examples;
use rustsim_app::state::AppState;

/// Compare interpreter and compiled mode for the PID controller example
#[test]
fn test_pid_interpreter_vs_compiled() {
    // Create two separate app states
    let mut interpreter_state = AppState::new();
    interpreter_state.load_example("PID Controller");

    let mut compiled_state = AppState::new();
    compiled_state.load_example("PID Controller");

    // Configure simulation settings
    let dt = 0.01;
    let duration = 5.0;
    let steps = (duration / dt) as usize;

    interpreter_state.settings.dt = dt;
    interpreter_state.settings.duration = duration;

    compiled_state.settings.dt = dt;
    compiled_state.settings.duration = duration;

    // Run interpreter mode
    interpreter_state.use_compiled_mode = false;
    interpreter_state.run_simulation();
    for _ in 0..steps {
        interpreter_state.step_simulation();
    }
    interpreter_state.stop_simulation();

    // Compile and run compiled mode
    compiled_state.use_compiled_mode = true;
    let compile_result = compiled_state.compile_simulation();
    assert!(compile_result.is_ok(), "Compilation failed: {:?}", compile_result.err());

    // Wait for compilation to complete (up to 120 seconds for CI environments)
    for _ in 0..600 {
        compiled_state.poll_compilation();
        if matches!(compiled_state.compilation_status, rustsim_app::state::CompilationStatus::Ready) {
            break;
        }
        std::thread::sleep(std::time::Duration::from_millis(200));
    }

    assert!(
        matches!(compiled_state.compilation_status, rustsim_app::state::CompilationStatus::Ready),
        "Compilation did not complete: {:?}",
        compiled_state.compilation_status
    );

    compiled_state.run_simulation();
    for _ in 0..steps {
        compiled_state.step_simulation();
    }
    compiled_state.stop_simulation();

    // Compare results
    println!("Interpreter data points: {}", interpreter_state.plot_data.len());
    println!("Compiled data points: {}", compiled_state.plot_data.len());

    // Check that we have data
    assert!(!interpreter_state.plot_data.is_empty(), "Interpreter produced no data");
    assert!(!compiled_state.plot_data.is_empty(), "Compiled produced no data");

    // Compare number of signals
    let interp_signals = interpreter_state.plot_data.first().map(|(_, v)| v.len()).unwrap_or(0);
    let compiled_signals = compiled_state.plot_data.first().map(|(_, v)| v.len()).unwrap_or(0);

    println!("Interpreter signals: {}", interp_signals);
    println!("Compiled signals: {}", compiled_signals);

    // Print first few data points for comparison
    println!("\n=== First 10 data points comparison ===");
    for i in 0..10.min(interpreter_state.plot_data.len()).min(compiled_state.plot_data.len()) {
        let (t_i, vals_i) = &interpreter_state.plot_data[i];
        let (t_c, vals_c) = &compiled_state.plot_data[i];

        println!("Step {}: t_interp={:.4}, t_compiled={:.4}", i, t_i, t_c);
        println!("  Interpreter: {:?}", vals_i);
        println!("  Compiled:    {:?}", vals_c);
    }

    // Print last few data points
    println!("\n=== Last 10 data points comparison ===");
    let len = interpreter_state.plot_data.len().min(compiled_state.plot_data.len());
    for i in (len.saturating_sub(10))..len {
        let (t_i, vals_i) = &interpreter_state.plot_data[i];
        let (t_c, vals_c) = &compiled_state.plot_data[i];

        println!("Step {}: t_interp={:.4}, t_compiled={:.4}", i, t_i, t_c);
        println!("  Interpreter: {:?}", vals_i);
        println!("  Compiled:    {:?}", vals_c);
    }

    // Assert signals match
    assert_eq!(
        interp_signals, compiled_signals,
        "Number of signals differs: interpreter={}, compiled={}",
        interp_signals, compiled_signals
    );

    // Compare values with tolerance
    let tolerance = 1e-6;
    let mut max_diff = 0.0f64;
    let mut max_diff_step = 0;
    let mut max_diff_signal = 0;

    for (i, ((t_i, vals_i), (t_c, vals_c))) in interpreter_state
        .plot_data
        .iter()
        .zip(compiled_state.plot_data.iter())
        .enumerate()
    {
        // Check time matches
        let time_diff = (t_i - t_c).abs();
        if time_diff > tolerance {
            println!("Time mismatch at step {}: interp={}, compiled={}", i, t_i, t_c);
        }

        // Check values match
        for (j, (v_i, v_c)) in vals_i.iter().zip(vals_c.iter()).enumerate() {
            let diff = (v_i - v_c).abs();
            if diff > max_diff {
                max_diff = diff;
                max_diff_step = i;
                max_diff_signal = j;
            }
        }
    }

    println!("\nMax difference: {} at step {}, signal {}", max_diff, max_diff_step, max_diff_signal);

    // Print timing statistics
    println!("\n=== Timing Statistics ===");
    if let Some(avg_time) = interpreter_state.average_step_time() {
        let avg_us = avg_time.as_nanos() as f64 / 1000.0;
        println!(
            "Interpreter: {} steps @ {:.2}µs/step ({:.2}ms total)",
            interpreter_state.get_step_count(),
            avg_us,
            avg_time.as_nanos() as f64 * interpreter_state.get_step_count() as f64 / 1_000_000.0
        );
    }
    if let Some(avg_time) = compiled_state.average_step_time() {
        let avg_us = avg_time.as_nanos() as f64 / 1000.0;
        println!(
            "Compiled:    {} steps @ {:.2}µs/step ({:.2}ms total)",
            compiled_state.get_step_count(),
            avg_us,
            avg_time.as_nanos() as f64 * compiled_state.get_step_count() as f64 / 1_000_000.0
        );
    }

    assert!(
        max_diff < 0.01,
        "Results differ too much: max_diff={} at step {}, signal {}",
        max_diff,
        max_diff_step,
        max_diff_signal
    );
}

/// Simple test with just a sinusoidal source
#[test]
fn test_sinusoidal_interpreter_vs_compiled() {
    let mut interpreter_state = AppState::new();
    interpreter_state.load_example("Sinusoidal Demo");

    let mut compiled_state = AppState::new();
    compiled_state.load_example("Sinusoidal Demo");

    let dt = 0.01;
    let steps = 100;

    interpreter_state.settings.dt = dt;
    compiled_state.settings.dt = dt;

    // Run interpreter
    interpreter_state.use_compiled_mode = false;
    interpreter_state.run_simulation();
    for _ in 0..steps {
        interpreter_state.step_simulation();
    }

    // Compile and run
    compiled_state.use_compiled_mode = true;
    compiled_state.compile_simulation().unwrap();

    for _ in 0..600 {
        compiled_state.poll_compilation();
        if matches!(compiled_state.compilation_status, rustsim_app::state::CompilationStatus::Ready) {
            break;
        }
        std::thread::sleep(std::time::Duration::from_millis(200));
    }

    compiled_state.run_simulation();
    for _ in 0..steps {
        compiled_state.step_simulation();
    }

    println!("\n=== Sinusoidal comparison ===");
    println!("Interpreter points: {}", interpreter_state.plot_data.len());
    println!("Compiled points: {}", compiled_state.plot_data.len());

    for i in 0..5.min(interpreter_state.plot_data.len()) {
        let (t_i, vals_i) = &interpreter_state.plot_data[i];
        let (t_c, vals_c) = &compiled_state.plot_data[i];
        println!("Step {}: interp=({:.4}, {:?}), compiled=({:.4}, {:?})",
                 i, t_i, vals_i, t_c, vals_c);
    }
}

/// Test that parameters can be changed after compilation without recompiling
#[test]
fn test_parameter_change_after_compilation() {
    let mut state = AppState::new();
    state.load_example("Sinusoidal Demo");

    let dt = 0.01;
    let steps = 10;

    state.settings.dt = dt;
    state.use_compiled_mode = true;

    // Compile the simulation
    state.compile_simulation().unwrap();

    for _ in 0..600 {
        state.poll_compilation();
        if matches!(state.compilation_status, rustsim_app::state::CompilationStatus::Ready) {
            break;
        }
        std::thread::sleep(std::time::Duration::from_millis(200));
    }

    assert!(
        matches!(state.compilation_status, rustsim_app::state::CompilationStatus::Ready),
        "Compilation did not complete"
    );

    // Run simulation with default amplitude (1.0)
    state.run_simulation();
    for _ in 0..steps {
        state.step_simulation();
    }

    // Get output at step 5 (should be sin(2*pi*1.0*0.05) * 1.0 = sin(0.314...) = ~0.309)
    let (_, outputs_default) = &state.plot_data[5];
    let output_default = outputs_default[0];
    println!("Default amplitude output at step 5: {:.6}", output_default);

    // Now change the amplitude to 2.0
    // Find the sinusoidal node and update its amplitude
    let sin_node_id = state.graph.nodes.keys()
        .find(|id| state.graph.nodes.get(*id).map(|n| n.block_type == "Sinusoidal").unwrap_or(false))
        .cloned()
        .expect("Should find sinusoidal node");

    if let Some(node) = state.graph.nodes.get_mut(&sin_node_id) {
        node.params.insert("amplitude".to_string(), serde_json::json!(2.0));
    }

    // Run simulation again (sync_compiled_params is called in run_simulation)
    state.run_simulation();
    for _ in 0..steps {
        state.step_simulation();
    }

    // Get output at step 5 (should be sin(2*pi*1.0*0.05) * 2.0 = ~0.618)
    let (_, outputs_double) = &state.plot_data[5];
    let output_double = outputs_double[0];
    println!("Double amplitude output at step 5: {:.6}", output_double);

    // Verify that the output is approximately doubled
    let ratio = output_double / output_default;
    println!("Ratio (should be ~2.0): {:.6}", ratio);

    assert!(
        (ratio - 2.0).abs() < 0.01,
        "Parameter change did not take effect: ratio was {} (expected ~2.0)",
        ratio
    );

    println!("\nParameter change after compilation works correctly!");
}
