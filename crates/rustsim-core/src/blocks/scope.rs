//! Data recorder (scope) block
//!
//! Provides multi-channel data recording with CSV export functionality,
//! matching PathSim's Scope.save() behavior.

use crate::block::{Block, DynamicBlock, StepResult};
use std::io::{self, Write};

/// Scope: multi-channel data recorder with CSV export
///
/// Records input signals into a fixed-size ring buffer. Once the buffer
/// is full, oldest samples are overwritten.
///
/// # Type Parameters
///
/// - `CHANNELS`: Number of input channels (const generic)
/// - `BUFFER_SIZE`: Maximum number of samples per channel (const generic)
///
/// # Example
///
/// ```ignore
/// // 2 channels, 1000 samples each
/// let mut scope = Scope::<2, 1000>::new();
///
/// scope.set_input(0, x);
/// scope.set_input(1, y);
/// scope.update(t);
/// scope.step(t, dt);
///
/// // Get recorded data
/// let data = scope.data();
/// for (time, values) in data {
///     println!("t={}, x={}, y={}", time, values[0], values[1]);
/// }
///
/// // Export to CSV
/// scope.save("output.csv")?; // Default labels: "port 0", "port 1"
/// scope.save_with_labels("output.csv", &["position", "velocity"])?;
/// ```
#[derive(Debug, Clone)]
pub struct Scope<const CHANNELS: usize, const BUFFER_SIZE: usize> {
    inputs: [f64; CHANNELS],
    outputs: [f64; CHANNELS], // Pass-through
    /// Circular buffer storing (time, [channel values])
    buffer: [(f64, [f64; CHANNELS]); BUFFER_SIZE],
    /// Current write position
    write_index: usize,
    /// Number of valid samples
    count: usize,
    /// Current time
    current_time: f64,

    // For buffer/revert
    buffered_write_index: usize,
    buffered_count: usize,
    buffered_buffer: [(f64, [f64; CHANNELS]); BUFFER_SIZE],
}

impl<const CHANNELS: usize, const BUFFER_SIZE: usize> Scope<CHANNELS, BUFFER_SIZE> {
    /// Create new scope
    pub fn new() -> Self {
        assert!(CHANNELS > 0, "Must have at least one channel");
        assert!(BUFFER_SIZE > 0, "Buffer size must be positive");

        Self {
            inputs: [0.0; CHANNELS],
            outputs: [0.0; CHANNELS],
            buffer: [(0.0, [0.0; CHANNELS]); BUFFER_SIZE],
            write_index: 0,
            count: 0,
            current_time: 0.0,
            buffered_write_index: 0,
            buffered_count: 0,
            buffered_buffer: [(0.0, [0.0; CHANNELS]); BUFFER_SIZE],
        }
    }

    /// Get number of recorded samples
    pub fn len(&self) -> usize {
        self.count
    }

    /// Check if buffer is empty
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// Check if buffer is full
    pub fn is_full(&self) -> bool {
        self.count == BUFFER_SIZE
    }

    /// Clear all recorded data
    pub fn clear(&mut self) {
        self.write_index = 0;
        self.count = 0;
        self.current_time = 0.0;
    }

    /// Get recorded data as slice of (time, values) tuples
    ///
    /// Returns data in chronological order (oldest to newest)
    pub fn data(&self) -> Vec<(f64, [f64; CHANNELS])> {
        if self.count == 0 {
            return Vec::new();
        }

        let mut result = Vec::with_capacity(self.count);

        // If buffer not full, data is at start
        if self.count < BUFFER_SIZE {
            for i in 0..self.count {
                result.push(self.buffer[i]);
            }
        } else {
            // Buffer is full, oldest data starts at write_index
            for i in 0..BUFFER_SIZE {
                let idx = (self.write_index + i) % BUFFER_SIZE;
                result.push(self.buffer[idx]);
            }
        }

        result
    }

    /// Get data for a specific channel
    pub fn channel_data(&self, channel: usize) -> Vec<(f64, f64)> {
        assert!(channel < CHANNELS, "Channel index out of bounds");

        self.data()
            .iter()
            .map(|(t, values)| (*t, values[channel]))
            .collect()
    }

    /// Get most recent sample
    pub fn last(&self) -> Option<(f64, [f64; CHANNELS])> {
        if self.count == 0 {
            return None;
        }

        let idx = if self.write_index == 0 {
            BUFFER_SIZE - 1
        } else {
            self.write_index - 1
        };

        Some(self.buffer[idx])
    }

    /// Manually record current inputs (called automatically in step())
    pub fn record(&mut self, time: f64) {
        self.buffer[self.write_index] = (time, self.inputs);
        self.write_index = (self.write_index + 1) % BUFFER_SIZE;

        if self.count < BUFFER_SIZE {
            self.count += 1;
        }
    }

    /// Save recorded data to CSV file with default labels
    ///
    /// Creates a CSV file with time column and channel columns.
    /// Channel labels default to "port 0", "port 1", etc.
    ///
    /// # CSV Format
    ///
    /// ```csv
    /// time [s],port 0,port 1,...
    /// 0.0,1.0,2.0,...
    /// 0.01,1.1,2.1,...
    /// ```
    ///
    /// # Example
    ///
    /// ```ignore
    /// let mut scope = Scope::<2, 1000>::new();
    /// // ... run simulation ...
    /// scope.save("output.csv")?;
    /// ```
    pub fn save(&self, filename: &str) -> io::Result<()> {
        // Generate default labels
        let labels: Vec<String> = (0..CHANNELS)
            .map(|i| format!("port {}", i))
            .collect();
        let label_refs: Vec<&str> = labels.iter().map(|s| s.as_str()).collect();

        self.save_with_labels(filename, &label_refs)
    }

    /// Save recorded data to CSV file with custom channel labels
    ///
    /// Creates a CSV file with time column and labeled channel columns.
    ///
    /// # CSV Format
    ///
    /// ```csv
    /// time [s],position,velocity,...
    /// 0.0,1.0,0.0,...
    /// 0.01,0.99,0.1,...
    /// ```
    ///
    /// # Example
    ///
    /// ```ignore
    /// let mut scope = Scope::<2, 1000>::new();
    /// // ... run simulation ...
    /// scope.save_with_labels("output.csv", &["position", "velocity"])?;
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The file cannot be created or written
    /// - The number of labels doesn't match the number of channels
    pub fn save_with_labels(&self, filename: &str, labels: &[&str]) -> io::Result<()> {
        // Ensure labels match channel count
        if labels.len() != CHANNELS {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Number of labels ({}) must match number of channels ({})",
                    labels.len(),
                    CHANNELS
                ),
            ));
        }

        // Add .csv extension if not present
        let filename = if filename.to_lowercase().ends_with(".csv") {
            filename.to_string()
        } else {
            format!("{}.csv", filename)
        };

        // Create CSV writer
        let file = std::fs::File::create(&filename)?;
        let mut wtr = csv::Writer::from_writer(file);

        // Write header row
        let mut header = vec!["time [s]".to_string()];
        header.extend(labels.iter().map(|&s| s.to_string()));
        wtr.write_record(&header)?;

        // Get data in chronological order
        let data = self.data();

        // Write each sample
        for (time, values) in data {
            let mut record = vec![time.to_string()];
            record.extend(values.iter().map(|v| v.to_string()));
            wtr.write_record(&record)?;
        }

        wtr.flush()?;
        Ok(())
    }

    /// Save recorded data to a writer with custom labels
    ///
    /// This is useful for testing or when you want to write to something
    /// other than a file (e.g., a string buffer).
    ///
    /// # Example
    ///
    /// ```ignore
    /// let mut buffer = Vec::new();
    /// scope.save_to_writer(&mut buffer, &["x", "y"])?;
    /// let csv_string = String::from_utf8(buffer)?;
    /// ```
    pub fn save_to_writer<W: Write>(&self, writer: W, labels: &[&str]) -> io::Result<()> {
        // Ensure labels match channel count
        if labels.len() != CHANNELS {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Number of labels ({}) must match number of channels ({})",
                    labels.len(),
                    CHANNELS
                ),
            ));
        }

        // Create CSV writer
        let mut wtr = csv::Writer::from_writer(writer);

        // Write header row
        let mut header = vec!["time [s]".to_string()];
        header.extend(labels.iter().map(|&s| s.to_string()));
        wtr.write_record(&header)?;

        // Get data in chronological order
        let data = self.data();

        // Write each sample
        for (time, values) in data {
            let mut record = vec![time.to_string()];
            record.extend(values.iter().map(|v| v.to_string()));
            wtr.write_record(&record)?;
        }

        wtr.flush()?;
        Ok(())
    }
}

impl<const CHANNELS: usize, const BUFFER_SIZE: usize> Default for Scope<CHANNELS, BUFFER_SIZE> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const CHANNELS: usize, const BUFFER_SIZE: usize> Block for Scope<CHANNELS, BUFFER_SIZE> {
    const NUM_INPUTS: usize = CHANNELS;
    const NUM_OUTPUTS: usize = CHANNELS;
    const IS_DYNAMIC: bool = true; // Has state (recorded data)

    #[inline]
    fn inputs(&self) -> &[f64] {
        &self.inputs
    }

    #[inline]
    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut self.inputs
    }

    #[inline]
    fn outputs(&self) -> &[f64] {
        &self.outputs
    }

    #[inline]
    fn outputs_mut(&mut self) -> &mut [f64] {
        &mut self.outputs
    }

    fn update(&mut self, _t: f64) {
        // Pass-through: outputs = inputs
        self.outputs.copy_from_slice(&self.inputs);
    }

    fn step(&mut self, t: f64, dt: f64) -> StepResult {
        // Record current inputs with timestamp
        self.record(self.current_time);

        // Advance time
        self.current_time = t + dt;

        StepResult::default()
    }

    fn buffer(&mut self) {
        self.buffered_write_index = self.write_index;
        self.buffered_count = self.count;
        self.buffered_buffer.copy_from_slice(&self.buffer);
    }

    fn revert(&mut self) {
        self.write_index = self.buffered_write_index;
        self.count = self.buffered_count;
        self.buffer.copy_from_slice(&self.buffered_buffer);
    }

    fn reset(&mut self) {
        self.inputs.fill(0.0);
        self.outputs.fill(0.0);
        self.clear();
        self.buffered_write_index = 0;
        self.buffered_count = 0;
        self.buffered_buffer.fill((0.0, [0.0; CHANNELS]));
    }
}

impl<const CHANNELS: usize, const BUFFER_SIZE: usize> DynamicBlock
    for Scope<CHANNELS, BUFFER_SIZE>
{
    fn state(&self) -> &[f64] {
        // State is the recorded data, but we return outputs for interface compliance
        &self.outputs
    }

    fn state_derivative(&self) -> &[f64] {
        // Scope has no continuous derivative
        &self.outputs // Return zeros
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scope_single_channel() {
        let mut scope = Scope::<1, 100>::new();
        let dt = 0.01;

        // Record some data
        for i in 0..50 {
            let t = i as f64 * dt;
            scope.set_input(0, t);
            scope.update(t);
            scope.step(t, dt);
        }

        assert_eq!(scope.len(), 50);
        assert!(!scope.is_empty());
        assert!(!scope.is_full());

        let data = scope.data();
        assert_eq!(data.len(), 50);

        // Check first and last values
        assert_eq!(data[0].1[0], 0.0);
        assert!((data[49].1[0] - 0.49).abs() < 1e-10);
    }

    #[test]
    fn test_scope_multi_channel() {
        let mut scope = Scope::<3, 100>::new();
        let dt = 0.01;

        for i in 0..50 {
            let t = i as f64 * dt;
            scope.set_input(0, t);
            scope.set_input(1, t * 2.0);
            scope.set_input(2, t * 3.0);
            scope.update(t);
            scope.step(t, dt);
        }

        let data = scope.data();
        assert_eq!(data.len(), 50);

        // Check multi-channel values
        let last = data.last().unwrap();
        assert!((last.1[0] - 0.49).abs() < 1e-10);
        assert!((last.1[1] - 0.98).abs() < 1e-10);
        assert!((last.1[2] - 1.47).abs() < 1e-10);
    }

    #[test]
    fn test_scope_buffer_overflow() {
        let mut scope = Scope::<1, 10>::new();
        let dt = 0.01;

        // Record more than buffer size
        for i in 0..20 {
            let t = i as f64 * dt;
            scope.set_input(0, i as f64);
            scope.update(t);
            scope.step(t, dt);
        }

        // Should only have last 10 samples
        assert_eq!(scope.len(), 10);
        assert!(scope.is_full());

        let data = scope.data();
        assert_eq!(data.len(), 10);

        // Oldest sample should be #10, newest should be #19
        assert_eq!(data[0].1[0], 10.0);
        assert_eq!(data[9].1[0], 19.0);
    }

    #[test]
    fn test_scope_clear() {
        let mut scope = Scope::<1, 100>::new();

        scope.set_input(0, 5.0);
        scope.update(0.0);
        scope.step(0.0, 0.01);

        assert_eq!(scope.len(), 1);

        scope.clear();
        assert_eq!(scope.len(), 0);
        assert!(scope.is_empty());
    }

    #[test]
    fn test_scope_channel_data() {
        let mut scope = Scope::<2, 100>::new();
        let dt = 0.01;

        for i in 0..10 {
            let t = i as f64 * dt;
            scope.set_input(0, i as f64);
            scope.set_input(1, i as f64 * 10.0);
            scope.update(t);
            scope.step(t, dt);
        }

        let ch0 = scope.channel_data(0);
        let ch1 = scope.channel_data(1);

        assert_eq!(ch0.len(), 10);
        assert_eq!(ch1.len(), 10);

        assert_eq!(ch0[5].1, 5.0);
        assert_eq!(ch1[5].1, 50.0);
    }

    #[test]
    fn test_scope_last() {
        let mut scope = Scope::<1, 100>::new();

        assert!(scope.last().is_none());

        scope.set_input(0, 1.0);
        scope.update(0.0);
        scope.step(0.0, 0.01);

        let last = scope.last().unwrap();
        assert_eq!(last.1[0], 1.0);

        scope.set_input(0, 2.0);
        scope.update(0.01);
        scope.step(0.01, 0.01);

        let last = scope.last().unwrap();
        assert_eq!(last.1[0], 2.0);
    }

    #[test]
    fn test_scope_buffer_revert() {
        let mut scope = Scope::<1, 100>::new();

        scope.set_input(0, 1.0);
        scope.update(0.0);
        scope.step(0.0, 0.01);

        scope.buffer();

        scope.set_input(0, 2.0);
        scope.update(0.01);
        scope.step(0.01, 0.01);

        assert_eq!(scope.len(), 2);

        scope.revert();
        assert_eq!(scope.len(), 1);
    }

    #[test]
    fn test_scope_reset() {
        let mut scope = Scope::<2, 100>::new();

        scope.set_input(0, 1.0);
        scope.set_input(1, 2.0);
        scope.update(0.0);
        scope.step(0.0, 0.01);

        scope.reset();

        assert_eq!(scope.len(), 0);
        assert_eq!(scope.get_output(0), 0.0);
        assert_eq!(scope.get_output(1), 0.0);
    }

    #[test]
    fn test_scope_passthrough() {
        let mut scope = Scope::<2, 100>::new();

        scope.set_input(0, 1.5);
        scope.set_input(1, 2.5);
        scope.update(0.0);

        // Outputs should equal inputs (pass-through)
        assert_eq!(scope.get_output(0), 1.5);
        assert_eq!(scope.get_output(1), 2.5);
    }
}
