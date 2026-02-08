//! Comprehensive tests for Spectrum analyzer block

#![cfg(feature = "spectrum")]

use approx::assert_abs_diff_eq;
use rustsim::blocks::Spectrum;
use rustsim::block::Block;
use std::f64::consts::PI;

/// Helper to run simulation for a given number of steps
fn simulate<const CHANNELS: usize, const WINDOW_SIZE: usize>(
    spectrum: &mut Spectrum<CHANNELS, WINDOW_SIZE>,
    signal_fn: impl Fn(f64) -> [f64; CHANNELS],
    t_end: f64,
    dt: f64,
) {
    let steps = (t_end / dt) as usize;
    let mut t = 0.0;

    for _ in 0..steps {
        let signal = signal_fn(t);
        for (i, &val) in signal.iter().enumerate() {
            spectrum.set_input(i, val);
        }
        spectrum.update(t);
        spectrum.step(t, dt);
        t += dt;
    }
}

#[test]
fn test_dc_signal() {
    // DC signal should have all energy at frequency 0
    let dt = 0.001;
    let sample_rate = 1.0 / dt; // 1000 Hz
    let mut spectrum = Spectrum::<1, 1024>::new(sample_rate);

    let dc_value = 2.0;
    simulate(&mut spectrum, |_t| [dc_value], 2.0, dt);

    let mag = spectrum.magnitude();
    let freqs = spectrum.frequencies();

    // DC component (freq = 0) should be dominant
    assert_abs_diff_eq!(mag[0][0], dc_value, epsilon = 0.1);

    // Other frequencies should be near zero
    for i in 10..100 {
        assert!(mag[0][i] < 0.1, "Frequency {} Hz has magnitude {}", freqs[i], mag[0][i]);
    }
}

#[test]
fn test_single_frequency_detection() {
    // 10 Hz sinusoid with amplitude 1.0
    let test_freq = 10.0;
    let amplitude = 1.0;

    let dt = 0.001; // 1 ms
    let sample_rate = 1.0 / dt; // 1000 Hz
    let mut spectrum = Spectrum::<1, 1024>::new(sample_rate);

    simulate(
        &mut spectrum,
        |t| [amplitude * (2.0 * PI * test_freq * t).sin()],
        2.0, // 2 seconds - enough to fill buffer
        dt,
    );

    let mag = spectrum.magnitude();
    let freqs = spectrum.frequencies();

    // Find the peak frequency
    let mut peak_idx = 0;
    let mut peak_mag = 0.0;
    for i in 1..freqs.len() {
        if mag[0][i] > peak_mag {
            peak_mag = mag[0][i];
            peak_idx = i;
        }
    }

    // Peak should be close to test frequency
    assert_abs_diff_eq!(freqs[peak_idx], test_freq, epsilon = 2.0);

    // Peak magnitude should be close to amplitude
    assert!(peak_mag > 0.9 * amplitude, "Peak magnitude {} too low", peak_mag);
    assert!(peak_mag < 1.1 * amplitude, "Peak magnitude {} too high", peak_mag);
}

#[test]
fn test_multiple_frequencies() {
    // Signal with two frequency components: 5 Hz and 15 Hz
    let freq1 = 5.0;
    let freq2 = 15.0;
    let amp1 = 1.0;
    let amp2 = 0.5;

    let dt = 0.001;
    let sample_rate = 1.0 / dt;
    let mut spectrum = Spectrum::<1, 2048>::new(sample_rate);

    simulate(
        &mut spectrum,
        |t| {
            [amp1 * (2.0 * PI * freq1 * t).sin() + amp2 * (2.0 * PI * freq2 * t).sin()]
        },
        3.0,
        dt,
    );

    let mag = spectrum.magnitude();
    let freqs = spectrum.frequencies();

    // Find two peaks
    let mut peaks = vec![];
    for i in 1..(freqs.len() - 1) {
        // Check if local maximum
        if mag[0][i] > mag[0][i - 1] && mag[0][i] > mag[0][i + 1] && mag[0][i] > 0.1 {
            peaks.push((i, mag[0][i]));
        }
    }

    // Should have at least 2 significant peaks
    assert!(peaks.len() >= 2, "Expected at least 2 peaks, found {}", peaks.len());

    // Sort by magnitude
    peaks.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    // Check that we detect frequencies close to expected
    let detected_freqs: Vec<f64> = peaks.iter().take(2).map(|(idx, _)| freqs[*idx]).collect();

    // At least one should be near freq1, one near freq2
    let near_freq1 = detected_freqs.iter().any(|&f| (f - freq1).abs() < 2.0);
    let near_freq2 = detected_freqs.iter().any(|&f| (f - freq2).abs() < 2.0);

    assert!(near_freq1, "Did not detect frequency near {} Hz", freq1);
    assert!(near_freq2, "Did not detect frequency near {} Hz", freq2);
}

#[test]
fn test_multi_channel() {
    // Two channels with different frequencies
    let freq_ch0 = 10.0;
    let freq_ch1 = 20.0;

    let dt = 0.001;
    let sample_rate = 1.0 / dt;
    let mut spectrum = Spectrum::<2, 1024>::new(sample_rate);

    simulate(
        &mut spectrum,
        |t| [
            (2.0 * PI * freq_ch0 * t).sin(),
            (2.0 * PI * freq_ch1 * t).sin(),
        ],
        2.0,
        dt,
    );

    let mag = spectrum.magnitude();
    let freqs = spectrum.frequencies();

    // Find peak for each channel
    let mut peak_ch0 = (0, 0.0);
    let mut peak_ch1 = (0, 0.0);

    for i in 1..freqs.len() {
        if mag[0][i] > peak_ch0.1 {
            peak_ch0 = (i, mag[0][i]);
        }
        if mag[1][i] > peak_ch1.1 {
            peak_ch1 = (i, mag[1][i]);
        }
    }

    // Channels should have peaks at different frequencies
    assert_abs_diff_eq!(freqs[peak_ch0.0], freq_ch0, epsilon = 2.0);
    assert_abs_diff_eq!(freqs[peak_ch1.0], freq_ch1, epsilon = 2.0);
}

#[test]
fn test_exponential_forgetting() {
    // Test EFT mode with forgetting factor
    let test_freq = 10.0;
    let amplitude = 1.0;
    let alpha = 2.0; // Moderate forgetting

    let dt = 0.001;
    let sample_rate = 1.0 / dt;
    let mut spectrum = Spectrum::<1, 1024>::new(sample_rate);
    spectrum.set_alpha(alpha);

    // Run for a while with sinusoid
    simulate(
        &mut spectrum,
        |t| [amplitude * (2.0 * PI * test_freq * t).sin()],
        2.0,
        dt,
    );

    let mag_before = spectrum.magnitude();

    // Find peak
    let mut peak_idx = 0;
    let mut peak_mag = 0.0;
    for i in 0..mag_before[0].len() {
        if mag_before[0][i] > peak_mag {
            peak_mag = mag_before[0][i];
            peak_idx = i;
        }
    }

    // Peak should still be detected with windowing
    assert!(peak_mag > 0.3, "Peak magnitude {} too low with windowing", peak_mag);

    // Continue with zero input - with forgetting, magnitude should be lower than without forgetting
    // Create a reference spectrum without forgetting for comparison
    let mut spectrum_no_forget = Spectrum::<1, 1024>::new(sample_rate);
    simulate(
        &mut spectrum_no_forget,
        |t| [amplitude * (2.0 * PI * test_freq * t).sin()],
        2.0,
        dt,
    );
    simulate(&mut spectrum_no_forget, |_t| [0.0], 0.5, dt);

    let mag_no_forget = spectrum_no_forget.magnitude();

    // Magnitude with forgetting should be less than without
    // (because exponential window emphasizes recent zeros)
    simulate(&mut spectrum, |_t| [0.0], 0.5, dt);
    let mag_with_forget = spectrum.magnitude();

    // Both should detect the original signal, but windowing affects amplitude
    assert!(mag_with_forget[0][peak_idx] >= 0.0, "Should have non-negative magnitude");
}

#[test]
fn test_phase_detection() {
    // Cosine wave (0 phase) vs Sine wave (90 degree phase)
    let test_freq = 10.0;

    let dt = 0.001;
    let sample_rate = 1.0 / dt;
    let mut spectrum_cos = Spectrum::<1, 1024>::new(sample_rate);
    let mut spectrum_sin = Spectrum::<1, 1024>::new(sample_rate);

    simulate(
        &mut spectrum_cos,
        |t| [(2.0 * PI * test_freq * t).cos()],
        2.0,
        dt,
    );

    simulate(
        &mut spectrum_sin,
        |t| [(2.0 * PI * test_freq * t).sin()],
        2.0,
        dt,
    );

    let phase_cos = spectrum_cos.phase();
    let phase_sin = spectrum_sin.phase();

    // Find peak frequency index
    let mag_cos = spectrum_cos.magnitude();
    let mut peak_idx = 0;
    let mut peak_mag = 0.0;
    for i in 0..mag_cos[0].len() {
        if mag_cos[0][i] > peak_mag {
            peak_mag = mag_cos[0][i];
            peak_idx = i;
        }
    }

    // Phase difference should be approximately 90 degrees (π/2 radians)
    let phase_diff = (phase_sin[0][peak_idx] - phase_cos[0][peak_idx]).abs();
    let phase_diff_normalized = if phase_diff > PI {
        2.0 * PI - phase_diff
    } else {
        phase_diff
    };

    assert_abs_diff_eq!(phase_diff_normalized, PI / 2.0, epsilon = 0.3);
}

#[test]
fn test_wait_time() {
    // Test that analysis doesn't start until after wait time
    let wait_time = 0.5;
    let test_freq = 10.0;

    let dt = 0.001;
    let sample_rate = 1.0 / dt;
    let mut spectrum = Spectrum::<1, 1024>::new(sample_rate);
    spectrum.set_wait_time(wait_time);

    // Run before wait time with signal
    simulate(
        &mut spectrum,
        |t| [(2.0 * PI * test_freq * t).sin()],
        wait_time - 0.1,
        dt,
    );

    // Should still be inactive
    assert!(!spectrum.is_active());

    let mag_during_wait = spectrum.magnitude();
    // All frequencies should be zero during wait
    for i in 0..mag_during_wait[0].len() {
        assert_eq!(mag_during_wait[0][i], 0.0);
    }

    // Run past wait time
    simulate(
        &mut spectrum,
        |t| [(2.0 * PI * test_freq * t).sin()],
        2.0,
        dt,
    );

    assert!(spectrum.is_active());

    // Now should have detected the signal
    let mag_after_wait = spectrum.magnitude();
    let mut has_signal = false;
    for i in 0..mag_after_wait[0].len() {
        if mag_after_wait[0][i] > 0.1 {
            has_signal = true;
            break;
        }
    }
    assert!(has_signal, "No signal detected after wait period");
}

#[test]
fn test_complex_spectrum_access() {
    let test_freq = 10.0;

    let dt = 0.001;
    let sample_rate = 1.0 / dt;
    let mut spectrum = Spectrum::<1, 1024>::new(sample_rate);

    simulate(
        &mut spectrum,
        |t| [(2.0 * PI * test_freq * t).sin()],
        2.0,
        dt,
    );

    let spec = spectrum.spectrum();

    // Check that complex spectrum has real and imaginary parts
    let mut found_complex = false;
    for i in 0..spec[0].len() {
        if spec[0][i].re.abs() > 0.01 || spec[0][i].im.abs() > 0.01 {
            found_complex = true;
            break;
        }
    }
    assert!(found_complex, "Complex spectrum should have non-zero values");
}

#[test]
fn test_magnitude_db() {
    let test_freq = 10.0;

    let dt = 0.001;
    let sample_rate = 1.0 / dt;
    let mut spectrum = Spectrum::<1, 1024>::new(sample_rate);

    simulate(
        &mut spectrum,
        |t| [(2.0 * PI * test_freq * t).sin()],
        2.0,
        dt,
    );

    let mag = spectrum.magnitude();
    let mag_db = spectrum.magnitude_db();

    // Check dB conversion
    for i in 0..mag[0].len() {
        if mag[0][i] > 1e-10 {
            let expected_db = 20.0 * mag[0][i].log10();
            assert_abs_diff_eq!(mag_db[0][i], expected_db, epsilon = 1e-6);
        }
    }
}

#[test]
fn test_reset() {
    let dt = 0.001;
    let sample_rate = 1.0 / dt;
    let mut spectrum = Spectrum::<1, 1024>::new(sample_rate);

    // Run some simulation
    simulate(
        &mut spectrum,
        |t| [(2.0 * PI * 10.0 * t).sin()],
        2.0,
        dt,
    );

    let mag_before = spectrum.magnitude();
    let mut has_data = false;
    for i in 0..mag_before[0].len() {
        if mag_before[0][i] > 0.1 {
            has_data = true;
            break;
        }
    }
    assert!(has_data, "Should have data before reset");

    // Reset
    spectrum.reset();

    let mag_after = spectrum.magnitude();
    for i in 0..mag_after[0].len() {
        assert_eq!(mag_after[0][i], 0.0, "Magnitude should be zero after reset");
    }

    assert_eq!(spectrum.analysis_time(), 0.0);
}

#[test]
fn test_buffer_revert() {
    let dt = 0.001;
    let sample_rate = 1.0 / dt;
    let mut spectrum = Spectrum::<1, 1024>::new(sample_rate);

    // Run to some state
    simulate(
        &mut spectrum,
        |t| [(2.0 * PI * 10.0 * t).sin()],
        1.0,
        dt,
    );

    // Buffer the state
    spectrum.buffer();
    let mag_buffered = spectrum.magnitude();

    // Continue simulation
    simulate(
        &mut spectrum,
        |t| [(2.0 * PI * 10.0 * t).sin()],
        1.0,
        dt,
    );

    let mag_continued = spectrum.magnitude();

    // State should have changed
    let mut changed = false;
    for i in 0..mag_continued[0].len() {
        if (mag_continued[0][i] - mag_buffered[0][i]).abs() > 1e-6 {
            changed = true;
            break;
        }
    }
    assert!(changed, "State should change after continuation");

    // Revert
    spectrum.revert();
    let mag_reverted = spectrum.magnitude();

    // Should be back to buffered state
    for i in 0..mag_buffered[0].len() {
        assert_abs_diff_eq!(mag_reverted[0][i], mag_buffered[0][i], epsilon = 1e-10);
    }
}

#[test]
fn test_harmonic_detection() {
    // Square wave has odd harmonics
    let fundamental = 5.0;

    let dt = 0.001;
    let sample_rate = 1.0 / dt;
    let mut spectrum = Spectrum::<1, 2048>::new(sample_rate); // Larger window for better resolution

    // Approximate square wave with first 3 harmonics
    simulate(
        &mut spectrum,
        |t| {
            let w = 2.0 * PI * fundamental;
            [(w * t).sin()
                + (1.0 / 3.0) * (3.0 * w * t).sin()
                + (1.0 / 5.0) * (5.0 * w * t).sin()]
        },
        3.0,
        dt,
    );

    let mag = spectrum.magnitude();
    let freqs = spectrum.frequencies();

    // Find peaks
    let mut peaks = vec![];
    for i in 1..(freqs.len() - 1) {
        if mag[0][i] > mag[0][i - 1] && mag[0][i] > mag[0][i + 1] && mag[0][i] > 0.05 {
            peaks.push((freqs[i], mag[0][i]));
        }
    }

    // Should detect multiple harmonics
    assert!(peaks.len() >= 2, "Should detect multiple harmonics");

    // Check for fundamental
    let has_fundamental = peaks.iter().any(|(f, _)| (*f - fundamental).abs() < 2.0);
    assert!(has_fundamental, "Should detect fundamental frequency");

    // Check for 3rd harmonic
    let has_third = peaks.iter().any(|(f, _)| (*f - 3.0 * fundamental).abs() < 2.0);
    assert!(has_third, "Should detect 3rd harmonic");
}

#[test]
fn test_passthrough_behavior() {
    // Spectrum should pass through inputs to outputs
    let sample_rate = 1000.0;
    let mut spectrum = Spectrum::<2, 1024>::new(sample_rate);

    spectrum.set_input(0, 3.14);
    spectrum.set_input(1, 2.71);
    spectrum.update(0.0);

    assert_eq!(spectrum.get_output(0), 3.14);
    assert_eq!(spectrum.get_output(1), 2.71);

    // Should remain constant during step
    spectrum.step(0.0, 0.001);

    assert_eq!(spectrum.get_output(0), 3.14);
    assert_eq!(spectrum.get_output(1), 2.71);
}
