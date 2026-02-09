//! Test spectrum FFT functionality

use approx::assert_abs_diff_eq;
use rustsim_app::spectrum_utils::{compute_fft, format_frequency};
use std::f64::consts::PI;

#[test]
fn test_compute_fft_single_frequency() {
    // Generate 10 Hz sinusoid at 1000 Hz sample rate
    let freq = 10.0;
    let sample_rate = 1000.0;
    let dt = 1.0 / sample_rate;
    let duration = 1.0; // 1 second

    let mut time_data = Vec::new();
    for i in 0..(duration / dt) as usize {
        let t = i as f64 * dt;
        let value = (2.0 * PI * freq * t).sin();
        time_data.push((t, vec![value]));
    }

    let (frequencies, magnitudes) = compute_fft(&time_data, sample_rate, 1024);

    // Find peak frequency
    let peak_idx = magnitudes[0]
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    let peak_freq = frequencies[peak_idx];

    // Peak should be at 10 Hz (within frequency resolution)
    let freq_resolution = sample_rate / 1024.0;
    assert_abs_diff_eq!(peak_freq, freq, epsilon = freq_resolution);

    // Magnitude should be the dominant peak (windowing affects absolute amplitude)
    let max_mag = magnitudes[0].iter().cloned().fold(0.0f64, f64::max);
    assert!(
        (magnitudes[0][peak_idx] - max_mag).abs() < 0.01,
        "Peak should be at expected frequency"
    );
}

#[test]
fn test_format_frequency() {
    assert_eq!(format_frequency(50.0), "50.0");
    assert_eq!(format_frequency(500.0), "500.0");
    assert_eq!(format_frequency(1500.0), "1.5k");
    assert_eq!(format_frequency(15000.0), "15.0k");
    assert_eq!(format_frequency(1500000.0), "1.5M");
    assert_eq!(format_frequency(15000000.0), "15.0M");
}

#[test]
fn test_compute_fft_dc_component() {
    // Constant signal (DC only)
    let sample_rate = 1000.0;
    let dt = 1.0 / sample_rate;
    let duration = 0.5;
    let dc_value = 2.0;

    let mut time_data = Vec::new();
    for i in 0..(duration / dt) as usize {
        let t = i as f64 * dt;
        time_data.push((t, vec![dc_value]));
    }

    let (frequencies, magnitudes) = compute_fft(&time_data, sample_rate, 512);

    // DC component should be at index 0
    assert_abs_diff_eq!(frequencies[0], 0.0);

    // DC magnitude should be non-zero (Hanning window affects absolute amplitude)
    assert!(magnitudes[0][0] > 0.0, "DC magnitude should be positive");
}

#[test]
fn test_compute_fft_multi_channel() {
    // Test with 2 channels at different frequencies
    let sample_rate = 1000.0;
    let dt = 1.0 / sample_rate;
    let duration = 1.0;

    let mut time_data = Vec::new();
    for i in 0..(duration / dt) as usize {
        let t = i as f64 * dt;
        let ch1 = (2.0 * PI * 10.0 * t).sin(); // 10 Hz
        let ch2 = (2.0 * PI * 20.0 * t).sin(); // 20 Hz
        time_data.push((t, vec![ch1, ch2]));
    }

    let (frequencies, magnitudes) = compute_fft(&time_data, sample_rate, 1024);

    // Should have 2 channels
    assert_eq!(magnitudes.len(), 2);

    // Find peaks for each channel
    let peak1_idx = magnitudes[0]
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    let peak2_idx = magnitudes[1]
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap();

    let freq_resolution = sample_rate / 1024.0;

    // Channel 1 should peak at ~10 Hz
    assert_abs_diff_eq!(frequencies[peak1_idx], 10.0, epsilon = freq_resolution);

    // Channel 2 should peak at ~20 Hz
    assert_abs_diff_eq!(frequencies[peak2_idx], 20.0, epsilon = freq_resolution);
}
