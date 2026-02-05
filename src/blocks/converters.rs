//! ADC and DAC converter blocks

use crate::block::{AlgebraicBlock, Block};

/// Analog-to-Digital Converter (ADC)
///
/// Models an ideal ADC that samples an analog input signal periodically,
/// quantizes it according to the specified number of bits and input span,
/// and outputs the resulting digital code on multiple output ports.
///
/// # Type Parameters
///
/// - `N_BITS`: Number of bits for the digital output code (const generic)
///
/// # Functionality
///
/// 1. Samples the analog input at intervals of `period`, starting after delay `tau`.
/// 2. Clips the input voltage to the defined `span` [min_voltage, max_voltage].
/// 3. Scales the clipped voltage to the range [0, 1].
/// 4. Quantizes the scaled value to an integer code between 0 and 2^N_BITS - 1 using flooring.
/// 5. Converts the integer code to an N_BITS binary representation.
/// 6. Outputs the binary code on ports 0 (LSB) to N_BITS-1 (MSB).
///
/// # Ideal Characteristics
///
/// - Instantaneous sampling at scheduled times.
/// - Perfect, noise-free quantization.
/// - No aperture jitter or other dynamic errors.
///
/// # Parameters
///
/// - `span`: The valid analog input value range [min_voltage, max_voltage].
///   Inputs outside this range will be clipped. Default is [-1, 1].
/// - `period`: Sampling period (time between samples). Default is 1 time unit.
/// - `tau`: Initial delay before the first sample is taken. Default is 0.
///
/// # Example
///
/// ```ignore
/// let mut adc = ADC::<4>::new([-1.0, 1.0], 1.0, 0.0);
/// adc.set_input(0, 0.5);
/// adc.update(0.0); // First sample at t=0
/// // Outputs on ports 0-3 represent the 4-bit digital code (LSB to MSB)
/// ```
#[derive(Debug, Clone)]
pub struct ADC<const N_BITS: usize> {
    input: f64,
    outputs: [f64; N_BITS],
    span: [f64; 2],
    period: f64,
    tau: f64,
    last_sample_time: Option<f64>,
}

impl<const N_BITS: usize> ADC<N_BITS> {
    /// Create a new ADC with span, sampling period, and delay
    pub fn new(span: [f64; 2], period: f64, tau: f64) -> Self {
        assert!(N_BITS > 0, "Number of bits must be positive");
        assert!(period > 0.0, "Sampling period must be positive");
        assert!(tau >= 0.0, "Delay must be non-negative");
        assert!(span[0] < span[1], "Invalid span: min must be less than max");

        Self {
            input: 0.0,
            outputs: [0.0; N_BITS],
            span,
            period,
            tau,
            last_sample_time: None,
        }
    }

    /// Create an ADC with default parameters
    pub fn default_config() -> Self {
        Self::new([-1.0, 1.0], 1.0, 0.0)
    }

    /// Get the span
    pub fn span(&self) -> [f64; 2] {
        self.span
    }

    /// Get the sampling period
    pub fn period(&self) -> f64 {
        self.period
    }

    /// Get the delay
    pub fn tau(&self) -> f64 {
        self.tau
    }

    /// Get the number of bits
    pub fn n_bits(&self) -> usize {
        N_BITS
    }

    /// Check if we should sample at this time
    fn should_sample(&self, t: f64) -> bool {
        if t < self.tau {
            return false;
        }

        match self.last_sample_time {
            None => t >= self.tau,
            Some(last_t) => {
                let elapsed = t - last_t;
                elapsed >= self.period - 1e-10
            }
        }
    }

    /// Perform the ADC conversion
    fn convert(&mut self) {
        // Clip input to span
        let lower = self.span[0];
        let upper = self.span[1];
        let clipped = self.input.clamp(lower, upper);

        // Scale to [0, 1]
        let scaled = (clipped - lower) / (upper - lower);

        // Quantize to integer code
        let levels = 1 << N_BITS; // 2^N_BITS
        let int_val = (scaled * levels as f64).floor() as usize;
        let int_val = int_val.min(levels - 1);

        // Convert to binary and set outputs (LSB to MSB)
        for i in 0..N_BITS {
            self.outputs[i] = ((int_val >> i) & 1) as f64;
        }
    }
}

impl<const N_BITS: usize> Default for ADC<N_BITS> {
    fn default() -> Self {
        Self::new([-1.0, 1.0], 1.0, 0.0)
    }
}

impl<const N_BITS: usize> Block for ADC<N_BITS> {
    const NUM_INPUTS: usize = 1;
    const NUM_OUTPUTS: usize = N_BITS;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        std::slice::from_ref(&self.input)
    }

    fn inputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.input)
    }

    fn outputs(&self) -> &[f64] {
        &self.outputs
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        &mut self.outputs
    }

    fn update(&mut self, t: f64) {
        if self.should_sample(t) {
            self.convert();
            self.last_sample_time = Some(t);
        }
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.outputs = [0.0; N_BITS];
        self.last_sample_time = None;
    }
}

impl<const N_BITS: usize> AlgebraicBlock for ADC<N_BITS> {}

/// Digital-to-Analog Converter (DAC)
///
/// Models an ideal DAC that reads a digital input code periodically from its input ports,
/// reconstructs the corresponding analog value based on the number of bits and output span,
/// and holds the output constant between updates.
///
/// # Type Parameters
///
/// - `N_BITS`: Number of digital input bits expected (const generic)
///
/// # Functionality
///
/// 1. Reads the digital code from input ports 0 (LSB) to N_BITS-1 (MSB) at intervals of `period`, starting after delay `tau`.
/// 2. Interprets the inputs as an unsigned binary integer code.
/// 3. Converts the integer code to a fractional value between 0 and (2^N_BITS - 1) / (2^N_BITS - 1).
/// 4. Scales this fractional value to the specified analog output `span`.
/// 5. Outputs the resulting analog value on `outputs[0]`.
/// 6. Holds the output value constant until the next scheduled update.
///
/// # Ideal Characteristics
///
/// - Instantaneous update at scheduled times.
/// - Perfect, noise-free reconstruction.
/// - No glitches or settling time.
///
/// # Parameters
///
/// - `span`: The analog output value range [min_voltage, max_voltage] corresponding
///   to the digital codes 0 and 2^N_BITS - 1, respectively.
///   Default is [-1, 1].
/// - `period`: Update period (time between output updates). Default is 1 time unit.
/// - `tau`: Initial delay before the first output update. Default is 0.
///
/// # Example
///
/// ```ignore
/// let mut dac = DAC::<4>::new([-1.0, 1.0], 1.0, 0.0);
/// // Set digital code (LSB to MSB)
/// dac.set_input(0, 1.0); // LSB = 1
/// dac.set_input(1, 0.0);
/// dac.set_input(2, 0.0);
/// dac.set_input(3, 1.0); // MSB = 1 (code = 9)
/// dac.update(0.0);
/// let analog_out = dac.get_output(0);
/// ```
#[derive(Debug, Clone)]
pub struct DAC<const N_BITS: usize> {
    inputs: [f64; N_BITS],
    output: f64,
    span: [f64; 2],
    period: f64,
    tau: f64,
    last_sample_time: Option<f64>,
}

impl<const N_BITS: usize> DAC<N_BITS> {
    /// Create a new DAC with span, update period, and delay
    pub fn new(span: [f64; 2], period: f64, tau: f64) -> Self {
        assert!(N_BITS > 0, "Number of bits must be positive");
        assert!(period > 0.0, "Update period must be positive");
        assert!(tau >= 0.0, "Delay must be non-negative");
        assert!(span[0] < span[1], "Invalid span: min must be less than max");

        Self {
            inputs: [0.0; N_BITS],
            output: 0.0,
            span,
            period,
            tau,
            last_sample_time: None,
        }
    }

    /// Create a DAC with default parameters
    pub fn default_config() -> Self {
        Self::new([-1.0, 1.0], 1.0, 0.0)
    }

    /// Get the span
    pub fn span(&self) -> [f64; 2] {
        self.span
    }

    /// Get the update period
    pub fn period(&self) -> f64 {
        self.period
    }

    /// Get the delay
    pub fn tau(&self) -> f64 {
        self.tau
    }

    /// Get the number of bits
    pub fn n_bits(&self) -> usize {
        N_BITS
    }

    /// Check if we should update at this time
    fn should_update(&self, t: f64) -> bool {
        if t < self.tau {
            return false;
        }

        match self.last_sample_time {
            None => t >= self.tau,
            Some(last_t) => {
                let elapsed = t - last_t;
                elapsed >= self.period - 1e-10
            }
        }
    }

    /// Perform the DAC conversion
    fn convert(&mut self) {
        // Convert binary inputs to integer (LSB to MSB)
        let mut val = 0;
        for i in 0..N_BITS {
            if self.inputs[i] != 0.0 {
                val |= 1 << i;
            }
        }

        // Scale to output span
        let lower = self.span[0];
        let upper = self.span[1];
        let levels = (1 << N_BITS) - 1; // 2^N_BITS - 1

        let scaled_val = if levels > 0 {
            val as f64 / levels as f64
        } else {
            0.0
        };

        self.output = lower + (upper - lower) * scaled_val;
    }
}

impl<const N_BITS: usize> Default for DAC<N_BITS> {
    fn default() -> Self {
        Self::new([-1.0, 1.0], 1.0, 0.0)
    }
}

impl<const N_BITS: usize> Block for DAC<N_BITS> {
    const NUM_INPUTS: usize = N_BITS;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        &self.inputs
    }

    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut self.inputs
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, t: f64) {
        if self.should_update(t) {
            self.convert();
            self.last_sample_time = Some(t);
        }
    }

    fn reset(&mut self) {
        self.inputs = [0.0; N_BITS];
        self.output = 0.0;
        self.last_sample_time = None;
    }
}

impl<const N_BITS: usize> AlgebraicBlock for DAC<N_BITS> {}

#[cfg(test)]
mod tests {
    use super::*;

    // ADC Tests

    #[test]
    fn test_adc_init() {
        // Default initialization
        let adc = ADC::<4>::default();
        assert_eq!(adc.n_bits(), 4);
        assert_eq!(adc.span(), [-1.0, 1.0]);
        assert_eq!(adc.period(), 1.0);
        assert_eq!(adc.tau(), 0.0);

        // Custom initialization
        let adc = ADC::<8>::new([0.0, 5.0], 0.1, 0.05);
        assert_eq!(adc.n_bits(), 8);
        assert_eq!(adc.span(), [0.0, 5.0]);
        assert_eq!(adc.period(), 0.1);
        assert_eq!(adc.tau(), 0.05);
    }

    #[test]
    fn test_adc_length() {
        // ADC has no direct passthrough (sampling only happens at specific times)
        let adc = ADC::<4>::default();
        assert_eq!(adc.inputs().len(), 1);
        assert_eq!(adc.outputs().len(), 4);
    }

    #[test]
    fn test_adc_output_ports() {
        // Test 4-bit ADC
        let mut adc = ADC::<4>::default();
        adc.update(0.0); // Trigger sampling
        assert_eq!(adc.outputs().len(), 4);

        // Test 8-bit ADC
        let mut adc = ADC::<8>::default();
        adc.update(0.0); // Trigger sampling
        assert_eq!(adc.outputs().len(), 8);
    }

    #[test]
    fn test_adc_sampling_midrange() {
        let mut adc = ADC::<4>::new([-1.0, 1.0], 1.0, 0.0);
        adc.set_input(0, 0.0); // Midrange

        // Trigger sampling at t=0
        adc.update(0.0);

        // For 4-bit ADC with span [-1,1]:
        // 0.0 maps to middle of range
        // scaled = (0 - (-1)) / 2 = 0.5
        // int_val = floor(0.5 * 16) = 8
        // binary: 1000 (MSB to LSB)
        // outputs: [0,0,0,1] (LSB to MSB)
        assert_eq!(adc.get_output(0), 0.0); // LSB
        assert_eq!(adc.get_output(1), 0.0);
        assert_eq!(adc.get_output(2), 0.0);
        assert_eq!(adc.get_output(3), 1.0); // MSB
    }

    #[test]
    fn test_adc_sampling_minimum() {
        let mut adc = ADC::<4>::new([-1.0, 1.0], 1.0, 0.0);
        adc.set_input(0, -1.0); // Minimum

        // Trigger sampling
        adc.update(0.0);

        // -1.0 maps to 0
        // binary: 0000
        assert_eq!(adc.get_output(0), 0.0);
        assert_eq!(adc.get_output(1), 0.0);
        assert_eq!(adc.get_output(2), 0.0);
        assert_eq!(adc.get_output(3), 0.0);
    }

    #[test]
    fn test_adc_sampling_maximum() {
        let mut adc = ADC::<4>::new([-1.0, 1.0], 1.0, 0.0);
        adc.set_input(0, 1.0); // Maximum

        // Trigger sampling
        adc.update(0.0);

        // 1.0 maps to max code (15)
        // binary: 1111
        assert_eq!(adc.get_output(0), 1.0); // LSB
        assert_eq!(adc.get_output(1), 1.0);
        assert_eq!(adc.get_output(2), 1.0);
        assert_eq!(adc.get_output(3), 1.0); // MSB
    }

    #[test]
    fn test_adc_clipping_above() {
        let mut adc = ADC::<4>::new([-1.0, 1.0], 1.0, 0.0);
        adc.set_input(0, 2.0); // Above maximum

        // Trigger sampling
        adc.update(0.0);

        // Should clip to 1.0, which maps to 1111
        assert_eq!(adc.get_output(0), 1.0);
        assert_eq!(adc.get_output(1), 1.0);
        assert_eq!(adc.get_output(2), 1.0);
        assert_eq!(adc.get_output(3), 1.0);
    }

    #[test]
    fn test_adc_clipping_below() {
        let mut adc = ADC::<4>::new([-1.0, 1.0], 1.0, 0.0);
        adc.set_input(0, -2.0); // Below minimum

        // Trigger sampling
        adc.update(0.0);

        // Should clip to -1.0, which maps to 0000
        assert_eq!(adc.get_output(0), 0.0);
        assert_eq!(adc.get_output(1), 0.0);
        assert_eq!(adc.get_output(2), 0.0);
        assert_eq!(adc.get_output(3), 0.0);
    }

    #[test]
    fn test_adc_different_spans() {
        let mut adc = ADC::<2>::new([0.0, 10.0], 1.0, 0.0);
        adc.set_input(0, 5.0); // Midrange

        // Trigger sampling
        adc.update(0.0);

        // scaled = (5 - 0) / 10 = 0.5
        // int_val = floor(0.5 * 4) = 2
        // binary: 10
        assert_eq!(adc.get_output(0), 0.0); // LSB
        assert_eq!(adc.get_output(1), 1.0); // MSB
    }

    // DAC Tests

    #[test]
    fn test_dac_init() {
        // Default initialization
        let dac = DAC::<4>::default();
        assert_eq!(dac.n_bits(), 4);
        assert_eq!(dac.span(), [-1.0, 1.0]);
        assert_eq!(dac.period(), 1.0);
        assert_eq!(dac.tau(), 0.0);

        // Custom initialization
        let dac = DAC::<8>::new([0.0, 5.0], 0.1, 0.05);
        assert_eq!(dac.n_bits(), 8);
        assert_eq!(dac.span(), [0.0, 5.0]);
        assert_eq!(dac.period(), 0.1);
        assert_eq!(dac.tau(), 0.05);
    }

    #[test]
    fn test_dac_length() {
        let dac = DAC::<4>::default();
        assert_eq!(dac.inputs().len(), 4);
        assert_eq!(dac.outputs().len(), 1);
    }

    #[test]
    fn test_dac_input_ports() {
        // Test 4-bit DAC
        let mut dac = DAC::<4>::default();
        for i in 0..4 {
            dac.set_input(i, 0.0);
        }
        assert_eq!(dac.inputs().len(), 4);

        // Test 8-bit DAC
        let mut dac = DAC::<8>::default();
        for i in 0..8 {
            dac.set_input(i, 0.0);
        }
        assert_eq!(dac.inputs().len(), 8);
    }

    #[test]
    fn test_dac_zero_code() {
        let mut dac = DAC::<4>::new([-1.0, 1.0], 1.0, 0.0);

        // Set all bits to 0
        for i in 0..4 {
            dac.set_input(i, 0.0);
        }

        // Trigger update
        dac.update(0.0);

        // Code 0 maps to minimum of span
        assert_eq!(dac.get_output(0), -1.0);
    }

    #[test]
    fn test_dac_max_code() {
        let mut dac = DAC::<4>::new([-1.0, 1.0], 1.0, 0.0);

        // Set all bits to 1 (code = 15)
        for i in 0..4 {
            dac.set_input(i, 1.0);
        }

        // Trigger update
        dac.update(0.0);

        // Code 15 maps to maximum of span
        assert_eq!(dac.get_output(0), 1.0);
    }

    #[test]
    fn test_dac_midrange_code() {
        let mut dac = DAC::<4>::new([-1.0, 1.0], 1.0, 0.0);

        // Set code to 8 (binary: 1000)
        // LSB=0, bit1=0, bit2=0, MSB=1
        dac.set_input(0, 0.0); // LSB
        dac.set_input(1, 0.0);
        dac.set_input(2, 0.0);
        dac.set_input(3, 1.0); // MSB

        // Trigger update
        dac.update(0.0);

        // Code 8 out of 15 levels
        // scaled_val = 8 / 15 = 0.533...
        // output = -1 + 2 * 0.533... = 0.0666...
        let expected = -1.0 + 2.0 * (8.0 / 15.0);
        assert!((dac.get_output(0) - expected).abs() < 1e-10);
    }

    #[test]
    fn test_dac_different_spans() {
        let mut dac = DAC::<2>::new([0.0, 10.0], 1.0, 0.0);

        // Set code to 2 (binary: 10)
        dac.set_input(0, 0.0); // LSB
        dac.set_input(1, 1.0); // MSB

        // Trigger update
        dac.update(0.0);

        // Code 2 out of 3 levels (2^2 - 1 = 3)
        // scaled_val = 2 / 3 = 0.666...
        // output = 0 + 10 * 0.666... = 6.666...
        let expected = 0.0 + 10.0 * (2.0 / 3.0);
        assert!((dac.get_output(0) - expected).abs() < 1e-10);
    }

    #[test]
    fn test_dac_1bit() {
        let mut dac = DAC::<1>::new([0.0, 1.0], 1.0, 0.0);

        // Code 0
        dac.set_input(0, 0.0);
        dac.update(0.0);
        assert_eq!(dac.get_output(0), 0.0);

        // Code 1
        dac.set_input(0, 1.0);
        dac.update(1.0);
        assert_eq!(dac.get_output(0), 1.0);
    }
}
