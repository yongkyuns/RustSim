//! Comprehensive tests for Scope CSV export functionality
//!
//! Tests CSV export matching PathSim's Scope.save() behavior:
//! - CSV format with time column and channel data
//! - Custom channel labels
//! - Proper header formatting
//! - Edge cases (empty scope, full buffer, single channel)

use rustsim::prelude::*;
use std::fs;
use std::io::Read;

// Helper function to read CSV file contents
fn read_csv_file(filename: &str) -> String {
    let mut file = fs::File::open(filename).expect("Failed to open CSV file");
    let mut contents = String::new();
    file.read_to_string(&mut contents)
        .expect("Failed to read CSV file");
    contents
}

// Helper function to clean up test files
fn cleanup_file(filename: &str) {
    let _ = fs::remove_file(filename);
}

#[test]
fn test_save_single_channel() {
    // Test saving single-channel data with default labels
    let mut scope = Scope::<1, 100>::new();
    let dt = 0.01;

    // Record some data
    for i in 0..10 {
        let t = i as f64 * dt;
        scope.set_input(0, t);
        scope.update(t);
        scope.step(t, dt);
    }

    // Save to CSV
    let filename = "test_single_channel.csv";
    scope.save(filename).expect("Failed to save CSV");

    // Read and verify contents
    let contents = read_csv_file(filename);
    let lines: Vec<&str> = contents.lines().collect();

    // Check header
    assert_eq!(lines[0], "time [s],port 0");

    // Check first data row
    assert!(lines[1].starts_with("0,0"));

    // Check last data row (9th sample)
    let last_line = lines[10];
    assert!(last_line.starts_with("0.09"));

    // Total lines: 1 header + 10 data rows
    assert_eq!(lines.len(), 11);

    cleanup_file(filename);
}

#[test]
fn test_save_multi_channel() {
    // Test saving multi-channel data with default labels
    let mut scope = Scope::<3, 100>::new();
    let dt = 0.01;

    // Record some data with different values per channel
    for i in 0..5 {
        let t = i as f64 * dt;
        scope.set_input(0, t);
        scope.set_input(1, t * 2.0);
        scope.set_input(2, t * 3.0);
        scope.update(t);
        scope.step(t, dt);
    }

    // Save to CSV
    let filename = "test_multi_channel.csv";
    scope.save(filename).expect("Failed to save CSV");

    // Read and verify contents
    let contents = read_csv_file(filename);
    let lines: Vec<&str> = contents.lines().collect();

    // Check header
    assert_eq!(lines[0], "time [s],port 0,port 1,port 2");

    // Check first data row
    assert!(lines[1].starts_with("0,0,0,0"));

    // Check third data row (t=0.02)
    let third_row: Vec<&str> = lines[3].split(',').collect();
    assert!(third_row[0].starts_with("0.02"));
    assert!(third_row[1].starts_with("0.02"));
    assert!(third_row[2].starts_with("0.04"));
    assert!(third_row[3].starts_with("0.06"));

    cleanup_file(filename);
}

#[test]
fn test_save_with_custom_labels() {
    // Test saving with custom channel labels
    let mut scope = Scope::<2, 100>::new();
    let dt = 0.01;

    // Record some data
    for i in 0..5 {
        let t = i as f64 * dt;
        scope.set_input(0, t);
        scope.set_input(1, -t);
        scope.update(t);
        scope.step(t, dt);
    }

    // Save with custom labels
    let filename = "test_custom_labels.csv";
    scope
        .save_with_labels(filename, &["position", "velocity"])
        .expect("Failed to save CSV");

    // Read and verify contents
    let contents = read_csv_file(filename);
    let lines: Vec<&str> = contents.lines().collect();

    // Check header with custom labels
    assert_eq!(lines[0], "time [s],position,velocity");

    // Check data rows exist
    assert_eq!(lines.len(), 6); // 1 header + 5 data rows

    cleanup_file(filename);
}

#[test]
fn test_save_empty_scope() {
    // Test saving an empty scope (no data recorded)
    let scope = Scope::<2, 100>::new();

    // Save to CSV
    let filename = "test_empty_scope.csv";
    scope.save(filename).expect("Failed to save CSV");

    // Read and verify contents
    let contents = read_csv_file(filename);
    let lines: Vec<&str> = contents.lines().collect();

    // Should only have header, no data rows
    assert_eq!(lines.len(), 1);
    assert_eq!(lines[0], "time [s],port 0,port 1");

    cleanup_file(filename);
}

#[test]
fn test_save_full_buffer() {
    // Test saving when buffer is full (circular buffer wraps)
    let mut scope = Scope::<1, 10>::new();
    let dt = 0.01;

    // Record more than buffer size (20 samples)
    for i in 0..20 {
        let t = i as f64 * dt;
        scope.set_input(0, i as f64);
        scope.update(t);
        scope.step(t, dt);
    }

    // Save to CSV
    let filename = "test_full_buffer.csv";
    scope.save(filename).expect("Failed to save CSV");

    // Read and verify contents
    let contents = read_csv_file(filename);
    let lines: Vec<&str> = contents.lines().collect();

    // Should have header + 10 data rows (buffer size)
    assert_eq!(lines.len(), 11);

    // First data row should be oldest (sample 10)
    let first_row: Vec<&str> = lines[1].split(',').collect();
    assert!(first_row[1].starts_with("10"));

    // Last data row should be newest (sample 19)
    let last_row: Vec<&str> = lines[10].split(',').collect();
    assert!(last_row[1].starts_with("19"));

    cleanup_file(filename);
}

#[test]
fn test_save_auto_adds_csv_extension() {
    // Test that .csv extension is added automatically if missing
    let mut scope = Scope::<1, 100>::new();

    scope.set_input(0, 1.0);
    scope.update(0.0);
    scope.step(0.0, 0.01);

    // Save without .csv extension
    let filename = "test_auto_extension";
    scope.save(filename).expect("Failed to save CSV");

    // File should be created with .csv extension
    let full_filename = "test_auto_extension.csv";
    assert!(
        fs::metadata(full_filename).is_ok(),
        "File with .csv extension should exist"
    );

    cleanup_file(full_filename);
}

#[test]
fn test_save_with_labels_mismatch() {
    // Test that save_with_labels returns error when label count doesn't match channels
    let scope = Scope::<3, 100>::new();

    // Try to save with wrong number of labels
    let filename = "test_label_mismatch.csv";
    let result = scope.save_with_labels(filename, &["label1", "label2"]); // 2 labels for 3 channels

    // Should return an error
    assert!(result.is_err());

    // Verify error message
    let err = result.unwrap_err();
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidInput);
    assert!(err.to_string().contains("must match"));

    // File should not be created
    assert!(fs::metadata(filename).is_err());
}

#[test]
fn test_save_preserves_chronological_order() {
    // Test that data is saved in chronological order
    let mut scope = Scope::<1, 100>::new();
    let dt = 0.01;

    // Record data
    for i in 0..10 {
        let t = i as f64 * dt;
        scope.set_input(0, i as f64);
        scope.update(t);
        scope.step(t, dt);
    }

    // Save to CSV
    let filename = "test_chronological.csv";
    scope.save(filename).expect("Failed to save CSV");

    // Read and verify contents
    let contents = read_csv_file(filename);
    let lines: Vec<&str> = contents.lines().collect();

    // Verify data is in chronological order
    for i in 1..lines.len() {
        let row: Vec<&str> = lines[i].split(',').collect();
        let expected_value = (i - 1) as f64;
        let actual_value: f64 = row[1].parse().expect("Failed to parse value");
        assert!((actual_value - expected_value).abs() < 1e-10);
    }

    cleanup_file(filename);
}

#[test]
fn test_save_to_writer() {
    // Test saving to a writer (in-memory buffer) instead of a file
    let mut scope = Scope::<2, 100>::new();
    let dt = 0.01;

    // Record some data
    for i in 0..5 {
        let t = i as f64 * dt;
        scope.set_input(0, t);
        scope.set_input(1, t * 2.0);
        scope.update(t);
        scope.step(t, dt);
    }

    // Save to in-memory buffer
    let mut buffer = Vec::new();
    scope
        .save_to_writer(&mut buffer, &["x", "y"])
        .expect("Failed to write CSV");

    // Convert to string and verify
    let contents = String::from_utf8(buffer).expect("Invalid UTF-8");
    let lines: Vec<&str> = contents.lines().collect();

    // Check header
    assert_eq!(lines[0], "time [s],x,y");

    // Check data rows
    assert_eq!(lines.len(), 6); // 1 header + 5 data rows

    // Check first data row
    assert!(lines[1].starts_with("0,0,0"));
}

#[test]
fn test_save_numeric_precision() {
    // Test that numeric values are saved with appropriate precision
    let mut scope = Scope::<1, 100>::new();

    // Record a value with many decimal places
    scope.set_input(0, std::f64::consts::PI);
    scope.update(0.0);
    scope.step(0.0, 0.01);

    // Save to CSV
    let filename = "test_precision.csv";
    scope.save(filename).expect("Failed to save CSV");

    // Read and verify contents
    let contents = read_csv_file(filename);
    let lines: Vec<&str> = contents.lines().collect();

    // Parse the saved value
    let row: Vec<&str> = lines[1].split(',').collect();
    let saved_value: f64 = row[1].parse().expect("Failed to parse value");

    // Should match original value (within floating point precision)
    assert!((saved_value - std::f64::consts::PI).abs() < 1e-10);

    cleanup_file(filename);
}

#[test]
fn test_save_multiple_times() {
    // Test that save can be called multiple times (overwriting file)
    let mut scope = Scope::<1, 100>::new();

    // First save with some data
    scope.set_input(0, 1.0);
    scope.update(0.0);
    scope.step(0.0, 0.01);

    let filename = "test_multiple_saves.csv";
    scope.save(filename).expect("Failed to save CSV first time");

    // Read first save
    let contents1 = read_csv_file(filename);
    let lines1: Vec<&str> = contents1.lines().collect();
    assert_eq!(lines1.len(), 2); // header + 1 data row

    // Add more data and save again
    scope.set_input(0, 2.0);
    scope.update(0.01);
    scope.step(0.01, 0.01);

    scope.save(filename).expect("Failed to save CSV second time");

    // Read second save (should have both data points)
    let contents2 = read_csv_file(filename);
    let lines2: Vec<&str> = contents2.lines().collect();
    assert_eq!(lines2.len(), 3); // header + 2 data rows

    cleanup_file(filename);
}

#[test]
fn test_save_with_negative_values() {
    // Test saving data with negative values
    let mut scope = Scope::<2, 100>::new();
    let dt = 0.01;

    // Record data with negative values
    for i in 0..5 {
        let t = i as f64 * dt;
        scope.set_input(0, t);
        scope.set_input(1, -t * 10.0);
        scope.update(t);
        scope.step(t, dt);
    }

    // Save to CSV
    let filename = "test_negative_values.csv";
    scope.save(filename).expect("Failed to save CSV");

    // Read and verify contents
    let contents = read_csv_file(filename);
    let lines: Vec<&str> = contents.lines().collect();

    // Check a row with negative values
    let row: Vec<&str> = lines[3].split(',').collect();
    let value: f64 = row[2].parse().expect("Failed to parse value");
    assert!(value < 0.0, "Should have negative value");

    cleanup_file(filename);
}
