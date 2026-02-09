//! Spectrum analyzer block for real-time Fourier analysis
//!
//! Uses FFT on a sliding window of samples for accurate spectral analysis.

use crate::block::{Block, DynamicBlock, StepResult};
use num_complex::Complex64;
use rustfft::{num_complex::Complex, FftPlanner};

/// Spectrum analyzer block using FFT on buffered samples
///
/// Stores input samples in a circular buffer and computes FFT when requested.
/// Supports exponential windowing for time-varying signals.
///
/// # Type Parameters
///
/// - `CHANNELS`: Number of input channels to analyze
/// - `WINDOW_SIZE`: Size of the sample buffer (must be power of 2 for efficient FFT)
///
/// # Example
///
/// ```rust,ignore
/// use rustsim::blocks::Spectrum;
///
/// // Analyze 2 channels with 1024-sample window
/// let mut spectrum = Spectrum::<2, 1024>::new(1000.0); // 1000 Hz sample rate
///
/// // Run simulation...
/// spectrum.set_input(0, signal1);
/// spectrum.set_input(1, signal2);
/// spectrum.update(t);
/// spectrum.step(t, dt);
///
/// // Get magnitude spectrum
/// let mag = spectrum.magnitude();
/// let freqs = spectrum.frequencies();
/// ```
#[derive(Clone)]
pub struct Spectrum<const CHANNELS: usize, const WINDOW_SIZE: usize> {
    /// Input signals
    inputs: [f64; CHANNELS],
    /// Pass-through outputs (same as inputs)
    outputs: [f64; CHANNELS],

    /// Sample rate in Hz
    sample_rate: f64,

    /// Circular buffer for each channel: [CHANNELS][WINDOW_SIZE]
    /// Stores recent samples for FFT computation
    buffers: Vec<Vec<f64>>,

    /// Current write position in circular buffer
    write_index: usize,

    /// Number of samples collected
    sample_count: usize,

    /// Exponential forgetting factor (0 = uniform window, >0 = exponential)
    alpha: f64,

    /// Wait time before starting analysis
    t_wait: f64,

    /// Time since wait period ended
    analysis_time: f64,

    /// Current time
    current_time: f64,

    /// Cached FFT results (computed on demand)
    cached_fft: Vec<Vec<Complex64>>,
    fft_valid: bool,

    // Buffered state for adaptive stepping
    buffered_write_index: usize,
    buffered_sample_count: usize,
    buffered_buffers: Vec<Vec<f64>>,
    buffered_analysis_time: f64,
}

impl<const CHANNELS: usize, const WINDOW_SIZE: usize> Spectrum<CHANNELS, WINDOW_SIZE> {
    /// Create spectrum analyzer with specified sample rate
    ///
    /// # Arguments
    ///
    /// - `sample_rate`: Sample rate in Hz (samples per second)
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// // 1024-sample window at 1000 Hz sample rate
    /// let spectrum = Spectrum::<1, 1024>::new(1000.0);
    /// ```
    pub fn new(sample_rate: f64) -> Self {
        assert!(CHANNELS > 0, "Must have at least one channel");
        assert!(WINDOW_SIZE > 0, "Window size must be positive");
        assert!(sample_rate > 0.0, "Sample rate must be positive");

        // Check if WINDOW_SIZE is power of 2 for efficient FFT
        if !WINDOW_SIZE.is_power_of_two() {
            eprintln!(
                "Warning: WINDOW_SIZE {} is not a power of 2. FFT will be slower.",
                WINDOW_SIZE
            );
        }

        Self {
            inputs: [0.0; CHANNELS],
            outputs: [0.0; CHANNELS],
            sample_rate,
            buffers: vec![vec![0.0; WINDOW_SIZE]; CHANNELS],
            write_index: 0,
            sample_count: 0,
            alpha: 0.0,
            t_wait: 0.0,
            analysis_time: 0.0,
            current_time: 0.0,
            cached_fft: vec![vec![Complex64::new(0.0, 0.0); WINDOW_SIZE / 2 + 1]; CHANNELS],
            fft_valid: false,
            buffered_write_index: 0,
            buffered_sample_count: 0,
            buffered_buffers: vec![vec![0.0; WINDOW_SIZE]; CHANNELS],
            buffered_analysis_time: 0.0,
        }
    }

    /// Set exponential forgetting factor
    ///
    /// - `alpha = 0.0`: Uniform window (all samples weighted equally)
    /// - `alpha > 0.0`: Exponential window (recent samples weighted more)
    pub fn set_alpha(&mut self, alpha: f64) {
        assert!(alpha >= 0.0, "Forgetting factor must be non-negative");
        self.alpha = alpha;
        self.fft_valid = false; // Invalidate cache
    }

    /// Set wait time before starting analysis
    pub fn set_wait_time(&mut self, t_wait: f64) {
        assert!(t_wait >= 0.0, "Wait time must be non-negative");
        self.t_wait = t_wait;
    }

    /// Get frequency bins in Hz
    ///
    /// Returns frequencies from DC (0 Hz) to Nyquist frequency (sample_rate / 2)
    pub fn frequencies(&self) -> Vec<f64> {
        let mut freqs = Vec::with_capacity(WINDOW_SIZE / 2 + 1);
        let freq_resolution = self.sample_rate / WINDOW_SIZE as f64;

        for i in 0..=(WINDOW_SIZE / 2) {
            freqs.push(i as f64 * freq_resolution);
        }

        freqs
    }

    /// Get frequency resolution (bin spacing) in Hz
    pub fn frequency_resolution(&self) -> f64 {
        self.sample_rate / WINDOW_SIZE as f64
    }

    /// Compute FFT for all channels
    fn compute_fft(&mut self) {
        if self.fft_valid {
            return; // Use cached result
        }

        let mut planner = FftPlanner::new();
        let fft = planner.plan_fft_forward(WINDOW_SIZE);

        for c in 0..CHANNELS {
            // Prepare buffer for FFT (with optional exponential windowing)
            let mut fft_buffer: Vec<Complex<f64>> = Vec::with_capacity(WINDOW_SIZE);

            for i in 0..WINDOW_SIZE {
                // Get sample in chronological order (oldest to newest)
                let buffer_idx = if self.sample_count < WINDOW_SIZE {
                    // Buffer not full yet, read from start
                    i
                } else {
                    // Buffer is full, oldest sample starts at write_index
                    (self.write_index + i) % WINDOW_SIZE
                };

                let sample = self.buffers[c][buffer_idx];

                // Apply exponential window if alpha > 0
                let weight = if self.alpha > 0.0 {
                    // Weight older samples less: exp(-alpha * (N - i))
                    // where i=0 is oldest, i=N-1 is newest
                    let age = (WINDOW_SIZE - 1 - i) as f64;
                    (-self.alpha * age / WINDOW_SIZE as f64).exp()
                } else {
                    1.0
                };

                fft_buffer.push(Complex::new(sample * weight, 0.0));
            }

            // Compute FFT
            fft.process(&mut fft_buffer);

            // Store positive frequencies (DC to Nyquist)
            // Normalize by window size to get correct amplitude
            let normalization = 1.0 / WINDOW_SIZE as f64;

            for i in 0..=(WINDOW_SIZE / 2) {
                let normalized = fft_buffer[i] * normalization;
                self.cached_fft[c][i] = Complex64::new(normalized.re, normalized.im);

                // For non-DC and non-Nyquist bins, multiply by 2 to account for negative frequencies
                if i > 0 && i < WINDOW_SIZE / 2 {
                    self.cached_fft[c][i] *= 2.0;
                }
            }
        }

        self.fft_valid = true;
    }

    /// Get complex spectrum for all channels
    ///
    /// Returns Vec of spectra, one per channel, each containing complex FFT values
    pub fn spectrum(&mut self) -> Vec<Vec<Complex64>> {
        if self.sample_count == 0 {
            return vec![vec![Complex64::new(0.0, 0.0); WINDOW_SIZE / 2 + 1]; CHANNELS];
        }

        self.compute_fft();
        self.cached_fft.clone()
    }

    /// Get magnitude spectrum for all channels
    pub fn magnitude(&mut self) -> Vec<Vec<f64>> {
        let spec = self.spectrum();
        spec.iter()
            .map(|channel_spec| channel_spec.iter().map(|c| c.norm()).collect())
            .collect()
    }

    /// Get phase spectrum in radians for all channels
    pub fn phase(&mut self) -> Vec<Vec<f64>> {
        let spec = self.spectrum();
        spec.iter()
            .map(|channel_spec| channel_spec.iter().map(|c| c.arg()).collect())
            .collect()
    }

    /// Get magnitude spectrum in decibels (20·log₁₀|X|)
    pub fn magnitude_db(&mut self) -> Vec<Vec<f64>> {
        let mag = self.magnitude();
        mag.iter()
            .map(|channel_mag| {
                channel_mag
                    .iter()
                    .map(|&m| {
                        if m > 0.0 {
                            20.0 * m.log10()
                        } else {
                            -std::f64::INFINITY
                        }
                    })
                    .collect()
            })
            .collect()
    }

    /// Get current analysis time (time since wait period ended)
    pub fn analysis_time(&self) -> f64 {
        self.analysis_time
    }

    /// Check if analysis has started (wait period has passed)
    pub fn is_active(&self) -> bool {
        self.current_time >= self.t_wait
    }

    /// Get number of samples collected
    pub fn sample_count(&self) -> usize {
        self.sample_count
    }

    /// Check if buffer is full
    pub fn is_full(&self) -> bool {
        self.sample_count >= WINDOW_SIZE
    }

    /// Clear all recorded data
    pub fn clear(&mut self) {
        for buffer in &mut self.buffers {
            buffer.fill(0.0);
        }
        self.write_index = 0;
        self.sample_count = 0;
        self.analysis_time = 0.0;
        self.fft_valid = false;
    }

    /// Record current inputs to buffer
    fn record_sample(&mut self) {
        // Only record if wait period has passed
        if !self.is_active() {
            return;
        }

        for c in 0..CHANNELS {
            self.buffers[c][self.write_index] = self.inputs[c];
        }

        self.write_index = (self.write_index + 1) % WINDOW_SIZE;

        if self.sample_count < WINDOW_SIZE {
            self.sample_count += 1;
        }

        // Invalidate FFT cache
        self.fft_valid = false;
    }
}

impl<const CHANNELS: usize, const WINDOW_SIZE: usize> Block for Spectrum<CHANNELS, WINDOW_SIZE> {
    const NUM_INPUTS: usize = CHANNELS;
    const NUM_OUTPUTS: usize = CHANNELS;
    const IS_DYNAMIC: bool = true;

    fn inputs(&self) -> &[f64] {
        &self.inputs
    }

    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut self.inputs
    }

    fn outputs(&self) -> &[f64] {
        &self.outputs
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        &mut self.outputs
    }

    fn update(&mut self, t: f64) {
        // Pass-through: outputs = inputs
        self.outputs.copy_from_slice(&self.inputs);

        // Update current time
        self.current_time = t;
    }

    fn step(&mut self, _t: f64, dt: f64) -> StepResult {
        // Record sample
        self.record_sample();

        // Update analysis time if active
        if self.is_active() {
            self.analysis_time += dt;
        }

        StepResult::default()
    }

    fn buffer(&mut self) {
        self.buffered_write_index = self.write_index;
        self.buffered_sample_count = self.sample_count;
        self.buffered_buffers.clone_from(&self.buffers);
        self.buffered_analysis_time = self.analysis_time;
    }

    fn revert(&mut self) {
        self.write_index = self.buffered_write_index;
        self.sample_count = self.buffered_sample_count;
        self.buffers.clone_from(&self.buffered_buffers);
        self.analysis_time = self.buffered_analysis_time;
        self.fft_valid = false;
    }

    fn reset(&mut self) {
        self.inputs.fill(0.0);
        self.outputs.fill(0.0);
        self.clear();
        self.buffered_write_index = 0;
        self.buffered_sample_count = 0;
        for buffer in &mut self.buffered_buffers {
            buffer.fill(0.0);
        }
        self.buffered_analysis_time = 0.0;
    }
}

impl<const CHANNELS: usize, const WINDOW_SIZE: usize> DynamicBlock
    for Spectrum<CHANNELS, WINDOW_SIZE>
{
    fn state(&self) -> &[f64] {
        // Return outputs for interface compliance
        &self.outputs
    }

    fn state_derivative(&self) -> &[f64] {
        // Spectrum has no continuous derivative
        &self.outputs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_spectrum_creation() {
        let spectrum = Spectrum::<1, 1024>::new(1000.0);

        let freqs = spectrum.frequencies();
        assert_eq!(freqs.len(), 513); // N/2 + 1
        assert_abs_diff_eq!(freqs[0], 0.0); // DC
        assert_abs_diff_eq!(freqs[512], 500.0, epsilon = 1e-10); // Nyquist at fs/2
    }

    #[test]
    fn test_frequency_resolution() {
        let spectrum = Spectrum::<1, 1024>::new(1000.0);
        let res = spectrum.frequency_resolution();
        assert_abs_diff_eq!(res, 1000.0 / 1024.0, epsilon = 1e-10);
    }

    #[test]
    fn test_passthrough() {
        let mut spectrum = Spectrum::<2, 1024>::new(1000.0);

        spectrum.set_input(0, 1.5);
        spectrum.set_input(1, 2.5);
        spectrum.update(0.0);

        assert_eq!(spectrum.get_output(0), 1.5);
        assert_eq!(spectrum.get_output(1), 2.5);
    }

    #[test]
    fn test_buffer_filling() {
        let mut spectrum = Spectrum::<1, 16>::new(1000.0);

        assert_eq!(spectrum.sample_count(), 0);
        assert!(!spectrum.is_full());

        // Add 10 samples
        for i in 0..10 {
            let t = i as f64 * 0.001;
            spectrum.set_input(0, i as f64);
            spectrum.update(t);
            spectrum.step(t, 0.001);
        }

        assert_eq!(spectrum.sample_count(), 10);
        assert!(!spectrum.is_full());

        // Fill buffer
        for i in 10..20 {
            let t = i as f64 * 0.001;
            spectrum.set_input(0, i as f64);
            spectrum.update(t);
            spectrum.step(t, 0.001);
        }

        assert_eq!(spectrum.sample_count(), 16);
        assert!(spectrum.is_full());
    }
}
