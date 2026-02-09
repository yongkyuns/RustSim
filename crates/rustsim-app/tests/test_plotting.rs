//! Integration tests for plotting functionality

use rustsim_app::plotting::{decimate_minmax, export_csv, PlotData};

#[test]
fn test_decimate_preserves_peaks() {
    // Create a signal with clear peaks and valleys
    let x: Vec<f64> = (0..100).map(|i| i as f64 * 0.1).collect();
    let y: Vec<f64> = x
        .iter()
        .map(|&t| (t * 2.0 * std::f64::consts::PI).sin())
        .collect();

    let (dec_x, dec_y) = decimate_minmax(&x, &y, 10);

    // Should have approximately 20 points (10 buckets * 2)
    assert!(dec_x.len() >= 10, "Too few points after decimation");
    assert!(dec_x.len() <= 25, "Too many points after decimation");
    assert_eq!(dec_x.len(), dec_y.len());

    // Find max and min in original
    let orig_max = y.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let orig_min = y.iter().copied().fold(f64::INFINITY, f64::min);

    // Find max and min in decimated
    let dec_max = dec_y.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let dec_min = dec_y.iter().copied().fold(f64::INFINITY, f64::min);

    // Decimated data should preserve the peaks and valleys
    assert!(
        (dec_max - orig_max).abs() < 0.1,
        "Max should be preserved: decimated = {}, original = {}",
        dec_max,
        orig_max
    );
    assert!(
        (dec_min - orig_min).abs() < 0.1,
        "Min should be preserved: decimated = {}, original = {}",
        dec_min,
        orig_min
    );
}

#[test]
fn test_decimate_reduces_points() {
    let n = 10000;
    let x: Vec<f64> = (0..n).map(|i| i as f64).collect();
    let y: Vec<f64> = (0..n).map(|i| (i as f64 * 0.01).sin()).collect();

    let target = 100;
    let (dec_x, dec_y) = decimate_minmax(&x, &y, target);

    // Should significantly reduce number of points
    assert!(
        dec_x.len() < n / 10,
        "Should significantly reduce points: {} < {}",
        dec_x.len(),
        n / 10
    );
    assert!(
        dec_x.len() >= target,
        "Should have at least target buckets: {} >= {}",
        dec_x.len(),
        target
    );
    assert!(
        dec_x.len() <= target * 2 + 10,
        "Should have at most 2 per bucket: {} <= {}",
        dec_x.len(),
        target * 2 + 10
    );
    assert_eq!(dec_x.len(), dec_y.len());
}

#[test]
fn test_decimate_chronological() {
    let x: Vec<f64> = (0..100).map(|i| i as f64).collect();
    let y: Vec<f64> = (0..100).map(|i| (i % 10) as f64).collect();

    let (dec_x, _) = decimate_minmax(&x, &y, 10);

    // X values should be monotonically increasing (chronological order)
    for i in 1..dec_x.len() {
        assert!(
            dec_x[i] >= dec_x[i - 1],
            "Points should be in chronological order: {} >= {}",
            dec_x[i],
            dec_x[i - 1]
        );
    }
}

#[test]
fn test_csv_format() {
    let data = PlotData::TimeDomain {
        time: vec![0.0, 0.1, 0.2],
        signals: vec![vec![1.0, 1.5, 2.0], vec![2.0, 2.5, 3.0]],
        labels: vec!["Signal 0".to_string(), "Signal 1".to_string()],
    };

    let csv = export_csv(&data);
    let lines: Vec<&str> = csv.lines().collect();

    assert_eq!(lines.len(), 4, "Should have header + 3 data rows");
    assert_eq!(lines[0], "time,Signal 0,Signal 1", "Header row incorrect");
    assert_eq!(lines[1], "0,1,2", "First data row incorrect");
    assert_eq!(lines[2], "0.1,1.5,2.5", "Second data row incorrect");
    assert_eq!(lines[3], "0.2,2,3", "Third data row incorrect");
}

#[test]
fn test_csv_escaping() {
    let data = PlotData::TimeDomain {
        time: vec![0.0],
        signals: vec![vec![1.0], vec![2.0]],
        labels: vec!["Signal, with comma".to_string(), "Signal \"quoted\"".to_string()],
    };

    let csv = export_csv(&data);

    // Comma should be escaped with quotes
    assert!(csv.contains("\"Signal, with comma\""), "Comma not escaped properly");

    // Quotes should be doubled
    assert!(csv.contains("\"Signal \"\"quoted\"\"\""), "Quotes not escaped properly");
}

#[test]
fn test_spectrum_csv_format() {
    let data = PlotData::Spectrum {
        frequency: vec![0.0, 10.0, 20.0],
        magnitudes: vec![vec![0.5, 0.6, 0.7], vec![0.3, 0.4, 0.5]],
        labels: vec!["Magnitude 0".to_string(), "Magnitude 1".to_string()],
    };

    let csv = export_csv(&data);
    let lines: Vec<&str> = csv.lines().collect();

    assert_eq!(lines.len(), 4, "Should have header + 3 data rows");
    assert_eq!(lines[0], "frequency,Magnitude 0,Magnitude 1", "Header row incorrect");
    assert_eq!(lines[1], "0,0.5,0.3", "First data row incorrect");
    assert_eq!(lines[2], "10,0.6,0.4", "Second data row incorrect");
    assert_eq!(lines[3], "20,0.7,0.5", "Third data row incorrect");
}
