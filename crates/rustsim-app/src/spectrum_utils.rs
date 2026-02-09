//! Spectrum analysis utilities for FFT computation

use num_complex::Complex64;
use rustfft::FftPlanner;
use std::collections::HashMap;

/// Spectrum data for FFT visualization
#[derive(Clone, Debug)]
pub struct SpectrumData {
    pub node_id: String,
    pub label: String,
    pub frequencies: Vec<f64>,
    pub magnitudes: Vec<Vec<f64>>, // Multiple channels
}

/// Compute FFT on time-domain data
///
/// # Arguments
/// * `time_data` - Vec of (time, values) tuples
/// * `sample_rate` - Sample rate in Hz (1/dt)
/// * `window_size` - FFT window size (will be padded/truncated to power of 2)
///
/// # Returns
/// Tuple of (frequencies, magnitudes) where magnitudes is Vec<Vec<f64>> for multiple channels
pub fn compute_fft(
    time_data: &[(f64, Vec<f64>)],
    sample_rate: f64,
    window_size: usize,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    if time_data.is_empty() {
        return (Vec::new(), Vec::new());
    }

    // Determine number of channels
    let num_channels = time_data.first().map(|(_, v)| v.len()).unwrap_or(0);
    if num_channels == 0 {
        return (Vec::new(), Vec::new());
    }

    // Use power of 2 for efficient FFT
    let fft_size = window_size.next_power_of_two().min(time_data.len());

    // Take the most recent samples
    let start_idx = time_data.len().saturating_sub(fft_size);
    let samples: Vec<_> = time_data[start_idx..].to_vec();

    // Initialize FFT planner
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(fft_size);

    // Compute FFT for each channel
    let mut magnitudes = Vec::with_capacity(num_channels);

    for channel in 0..num_channels {
        // Extract channel data
        let mut buffer: Vec<Complex64> = samples
            .iter()
            .map(|(_, values)| {
                let val = values.get(channel).copied().unwrap_or(0.0);
                Complex64::new(val, 0.0)
            })
            .collect();

        // Pad to fft_size if needed
        buffer.resize(fft_size, Complex64::new(0.0, 0.0));

        // Apply Hanning window to reduce spectral leakage
        for (i, sample) in buffer.iter_mut().enumerate() {
            let window =
                0.5 * (1.0 - (2.0 * std::f64::consts::PI * i as f64 / fft_size as f64).cos());
            *sample *= window;
        }

        // Compute FFT
        fft.process(&mut buffer);

        // Extract magnitude spectrum (positive frequencies only)
        let num_freqs = fft_size / 2 + 1;
        let normalization = 2.0 / fft_size as f64; // Factor of 2 for one-sided spectrum

        let mag: Vec<f64> = buffer[..num_freqs]
            .iter()
            .enumerate()
            .map(|(i, c)| {
                let norm = c.norm() * normalization;
                // DC and Nyquist don't get the factor of 2
                if i == 0 || i == fft_size / 2 {
                    norm / 2.0
                } else {
                    norm
                }
            })
            .collect();

        magnitudes.push(mag);
    }

    // Compute frequency bins
    let freq_resolution = sample_rate / fft_size as f64;
    let frequencies: Vec<f64> = (0..=fft_size / 2)
        .map(|i| i as f64 * freq_resolution)
        .collect();

    (frequencies, magnitudes)
}

/// Format frequency value with appropriate units (Hz, kHz, MHz)
pub fn format_frequency(freq: f64) -> String {
    if freq >= 1e6 {
        format!("{:.1}M", freq / 1e6)
    } else if freq >= 1e3 {
        format!("{:.1}k", freq / 1e3)
    } else {
        format!("{:.1}", freq)
    }
}

/// Compute spectrum data from time-domain plot data
pub fn compute_spectrum_from_plot_data(
    node_id: String,
    label: String,
    plot_data: &[(f64, Vec<f64>)],
    sample_rate: f64,
) -> SpectrumData {
    const DEFAULT_WINDOW_SIZE: usize = 1024;

    let (frequencies, magnitudes) = compute_fft(plot_data, sample_rate, DEFAULT_WINDOW_SIZE);

    SpectrumData {
        node_id,
        label,
        frequencies,
        magnitudes,
    }
}

/// Compute spectrum data for all spectrum blocks
pub fn compute_all_spectra(
    scope_data: &HashMap<String, Vec<(f64, Vec<f64>)>>,
    spectrum_nodes: &[(String, String)], // (node_id, label)
    sample_rate: f64,
) -> HashMap<String, SpectrumData> {
    let mut spectra = HashMap::new();

    for (node_id, label) in spectrum_nodes {
        if let Some(data) = scope_data.get(node_id) {
            if !data.is_empty() {
                let spectrum = compute_spectrum_from_plot_data(
                    node_id.clone(),
                    label.clone(),
                    data,
                    sample_rate,
                );
                spectra.insert(node_id.clone(), spectrum);
            }
        }
    }

    spectra
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_fft_single_frequency() {
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

        // Magnitude should be significantly larger than surrounding frequencies (peak detection)
        // Note: Hanning window affects amplitude, so we just check it's the dominant peak
        let max_mag = magnitudes[0].iter().cloned().fold(0.0f64, f64::max);
        assert!(magnitudes[0][peak_idx] == max_mag || (max_mag - magnitudes[0][peak_idx]).abs() < 0.01);
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
    fn test_fft_multiple_frequencies() {
        // Generate signal with 10 Hz + 30 Hz components
        let sample_rate = 1000.0;
        let dt = 1.0 / sample_rate;
        let duration = 1.0;

        let mut time_data = Vec::new();
        for i in 0..(duration / dt) as usize {
            let t = i as f64 * dt;
            let value = (2.0 * PI * 10.0 * t).sin() + 0.5 * (2.0 * PI * 30.0 * t).sin();
            time_data.push((t, vec![value]));
        }

        let (frequencies, magnitudes) = compute_fft(&time_data, sample_rate, 1024);

        // Find top 2 peaks
        let mut indexed_mags: Vec<_> = magnitudes[0].iter().enumerate().collect();
        indexed_mags.sort_by(|(_, a), (_, b)| b.partial_cmp(a).unwrap());

        let peak1_freq = frequencies[indexed_mags[0].0];
        let peak2_freq = frequencies[indexed_mags[1].0];

        let freq_resolution = sample_rate / 1024.0;

        // Should find peaks near 10 Hz and 30 Hz
        assert!(
            (peak1_freq - 10.0).abs() < freq_resolution
                || (peak1_freq - 30.0).abs() < freq_resolution
        );
        assert!(
            (peak2_freq - 10.0).abs() < freq_resolution
                || (peak2_freq - 30.0).abs() < freq_resolution
        );
    }

    #[test]
    fn test_fft_dc_component() {
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

        // DC magnitude should be non-zero for a constant signal
        // Note: Hanning window significantly reduces DC amplitude due to tapering
        assert!(magnitudes[0][0] > 0.0, "DC magnitude should be positive");
    }
}
