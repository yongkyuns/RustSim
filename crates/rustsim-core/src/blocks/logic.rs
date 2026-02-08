//! Logic and comparison blocks

use crate::block::{AlgebraicBlock, Block};

/// Comparator block with optional hysteresis
///
/// Outputs 1.0 when input exceeds threshold (+ hysteresis/2), 0.0 when input
/// falls below threshold (- hysteresis/2).
///
/// # Example
///
/// ```ignore
/// // Simple threshold at 0.5
/// let mut comp = Comparator::new(0.5);
/// comp.set_input(0, 0.6);
/// comp.update(0.0);
/// assert_eq!(comp.get_output(0), 1.0);
///
/// // With hysteresis
/// let mut comp = Comparator::with_hysteresis(0.5, 0.2);
/// // Turns on at 0.6, turns off at 0.4
/// ```
#[derive(Debug, Clone, Copy)]
pub struct Comparator {
    input: f64,
    output: f64,
    threshold: f64,
    hysteresis: f64,
}

impl Comparator {
    /// Create comparator with threshold, no hysteresis
    pub fn new(threshold: f64) -> Self {
        Self {
            input: 0.0,
            output: 0.0,
            threshold,
            hysteresis: 0.0,
        }
    }

    /// Create comparator with threshold and hysteresis
    ///
    /// - Turns on when input > threshold + hysteresis/2
    /// - Turns off when input < threshold - hysteresis/2
    pub fn with_hysteresis(threshold: f64, hysteresis: f64) -> Self {
        assert!(hysteresis >= 0.0, "hysteresis must be non-negative");
        Self {
            input: 0.0,
            output: 0.0,
            threshold,
            hysteresis,
        }
    }

    /// Get threshold
    pub fn threshold(&self) -> f64 {
        self.threshold
    }

    /// Set threshold
    pub fn set_threshold(&mut self, threshold: f64) {
        self.threshold = threshold;
    }

    /// Get hysteresis
    pub fn hysteresis(&self) -> f64 {
        self.hysteresis
    }

    /// Set hysteresis
    pub fn set_hysteresis(&mut self, hysteresis: f64) {
        assert!(hysteresis >= 0.0, "hysteresis must be non-negative");
        self.hysteresis = hysteresis;
    }
}

impl Block for Comparator {
    const NUM_INPUTS: usize = 1;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    #[inline]
    fn inputs(&self) -> &[f64] {
        std::slice::from_ref(&self.input)
    }

    #[inline]
    fn inputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.input)
    }

    #[inline]
    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    #[inline]
    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    #[inline]
    fn update(&mut self, _t: f64) {
        let half_hyst = self.hysteresis / 2.0;
        if self.input > self.threshold + half_hyst {
            self.output = 1.0;
        } else if self.input < self.threshold - half_hyst {
            self.output = 0.0;
        }
        // Else maintain current state (hysteresis region)
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = 0.0;
    }
}

impl AlgebraicBlock for Comparator {}

/// Switch block - select between inputs based on control signal
///
/// # Ports
/// - Input 0: Control signal
/// - Input 1: Value when control < threshold
/// - Input 2: Value when control >= threshold
/// - Output: Selected value
///
/// # Example
///
/// ```ignore
/// let mut sw = Switch::new(0.5);
/// sw.inputs_mut()[0] = 0.3; // control < 0.5
/// sw.inputs_mut()[1] = 10.0;
/// sw.inputs_mut()[2] = 20.0;
/// sw.update(0.0);
/// assert_eq!(sw.get_output(0), 10.0);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct Switch {
    inputs: [f64; 3],
    output: f64,
    threshold: f64,
}

impl Switch {
    /// Create switch with given threshold
    pub fn new(threshold: f64) -> Self {
        Self {
            inputs: [0.0; 3],
            output: 0.0,
            threshold,
        }
    }

    /// Get threshold
    pub fn threshold(&self) -> f64 {
        self.threshold
    }

    /// Set threshold
    pub fn set_threshold(&mut self, threshold: f64) {
        self.threshold = threshold;
    }
}

impl Default for Switch {
    fn default() -> Self {
        Self::new(0.0)
    }
}

impl Block for Switch {
    const NUM_INPUTS: usize = 3;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

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
        std::slice::from_ref(&self.output)
    }

    #[inline]
    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    #[inline]
    fn update(&mut self, _t: f64) {
        let control = self.inputs[0];
        self.output = if control < self.threshold {
            self.inputs[1]
        } else {
            self.inputs[2]
        };
    }

    fn reset(&mut self) {
        self.inputs.fill(0.0);
        self.output = 0.0;
    }
}

impl AlgebraicBlock for Switch {}

/// Relay (Schmitt trigger) - on/off output with hysteresis
///
/// A relay that switches on at `on_threshold` and off at `off_threshold`.
/// Requires `off_threshold < on_threshold` for proper hysteresis.
///
/// # Example
///
/// ```ignore
/// // Relay that turns on at 1.0, off at 0.5
/// let mut relay = Relay::new(0.5, 1.0);
/// relay.set_input(0, 0.7); // Below on threshold
/// relay.update(0.0);
/// assert_eq!(relay.get_output(0), 0.0);
///
/// relay.set_input(0, 1.2); // Above on threshold
/// relay.update(0.0);
/// assert_eq!(relay.get_output(0), 1.0);
///
/// relay.set_input(0, 0.7); // Above off threshold but below on
/// relay.update(0.0);
/// assert_eq!(relay.get_output(0), 1.0); // Still on
///
/// relay.set_input(0, 0.3); // Below off threshold
/// relay.update(0.0);
/// assert_eq!(relay.get_output(0), 0.0); // Now off
/// ```
#[derive(Debug, Clone, Copy)]
pub struct Relay {
    input: f64,
    output: f64,
    off_threshold: f64,
    on_threshold: f64,
}

impl Relay {
    /// Create relay with off and on thresholds
    pub fn new(off_threshold: f64, on_threshold: f64) -> Self {
        assert!(
            off_threshold < on_threshold,
            "off_threshold must be < on_threshold"
        );
        Self {
            input: 0.0,
            output: 0.0,
            off_threshold,
            on_threshold,
        }
    }

    /// Get off threshold
    pub fn off_threshold(&self) -> f64 {
        self.off_threshold
    }

    /// Get on threshold
    pub fn on_threshold(&self) -> f64 {
        self.on_threshold
    }

    /// Set thresholds
    pub fn set_thresholds(&mut self, off_threshold: f64, on_threshold: f64) {
        assert!(
            off_threshold < on_threshold,
            "off_threshold must be < on_threshold"
        );
        self.off_threshold = off_threshold;
        self.on_threshold = on_threshold;
    }
}

impl Block for Relay {
    const NUM_INPUTS: usize = 1;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    #[inline]
    fn inputs(&self) -> &[f64] {
        std::slice::from_ref(&self.input)
    }

    #[inline]
    fn inputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.input)
    }

    #[inline]
    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    #[inline]
    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    #[inline]
    fn update(&mut self, _t: f64) {
        if self.input >= self.on_threshold {
            self.output = 1.0;
        } else if self.input <= self.off_threshold {
            self.output = 0.0;
        }
        // Else maintain current state (between thresholds)
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = 0.0;
    }
}

impl AlgebraicBlock for Relay {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_comparator_basic() {
        let mut comp = Comparator::new(0.5);

        comp.set_input(0, 0.3);
        comp.update(0.0);
        assert_eq!(comp.get_output(0), 0.0);

        comp.set_input(0, 0.7);
        comp.update(0.0);
        assert_eq!(comp.get_output(0), 1.0);
    }

    #[test]
    fn test_comparator_hysteresis() {
        let mut comp = Comparator::with_hysteresis(0.5, 0.2);
        // Turns on at 0.6, turns off at 0.4

        comp.set_input(0, 0.3);
        comp.update(0.0);
        assert_eq!(comp.get_output(0), 0.0);

        comp.set_input(0, 0.7);
        comp.update(0.0);
        assert_eq!(comp.get_output(0), 1.0);

        // In hysteresis region, should stay on
        comp.set_input(0, 0.5);
        comp.update(0.0);
        assert_eq!(comp.get_output(0), 1.0);

        // Below off threshold
        comp.set_input(0, 0.3);
        comp.update(0.0);
        assert_eq!(comp.get_output(0), 0.0);

        // In hysteresis region, should stay off
        comp.set_input(0, 0.5);
        comp.update(0.0);
        assert_eq!(comp.get_output(0), 0.0);
    }

    #[test]
    fn test_switch() {
        let mut sw = Switch::new(0.5);

        // Control below threshold, select input 1
        sw.inputs_mut()[0] = 0.3;
        sw.inputs_mut()[1] = 10.0;
        sw.inputs_mut()[2] = 20.0;
        sw.update(0.0);
        assert_eq!(sw.get_output(0), 10.0);

        // Control above threshold, select input 2
        sw.inputs_mut()[0] = 0.7;
        sw.update(0.0);
        assert_eq!(sw.get_output(0), 20.0);

        // Control at threshold, select input 2
        sw.inputs_mut()[0] = 0.5;
        sw.update(0.0);
        assert_eq!(sw.get_output(0), 20.0);
    }

    #[test]
    fn test_relay() {
        let mut relay = Relay::new(0.5, 1.0);

        // Below on threshold - stays off
        relay.set_input(0, 0.7);
        relay.update(0.0);
        assert_eq!(relay.get_output(0), 0.0);

        // Above on threshold - turns on
        relay.set_input(0, 1.2);
        relay.update(0.0);
        assert_eq!(relay.get_output(0), 1.0);

        // Between thresholds - stays on
        relay.set_input(0, 0.7);
        relay.update(0.0);
        assert_eq!(relay.get_output(0), 1.0);

        // Below off threshold - turns off
        relay.set_input(0, 0.3);
        relay.update(0.0);
        assert_eq!(relay.get_output(0), 0.0);

        // Between thresholds - stays off
        relay.set_input(0, 0.7);
        relay.update(0.0);
        assert_eq!(relay.get_output(0), 0.0);
    }

    #[test]
    fn test_comparator_reset() {
        let mut comp = Comparator::new(0.5);
        comp.set_input(0, 1.0);
        comp.update(0.0);
        comp.reset();
        assert_eq!(comp.get_input(0), 0.0);
        assert_eq!(comp.get_output(0), 0.0);
    }

    #[test]
    fn test_switch_reset() {
        let mut sw = Switch::new(0.5);
        sw.inputs_mut()[0] = 0.7;
        sw.inputs_mut()[1] = 10.0;
        sw.inputs_mut()[2] = 20.0;
        sw.update(0.0);
        sw.reset();
        assert_eq!(sw.get_input(0), 0.0);
        assert_eq!(sw.get_input(1), 0.0);
        assert_eq!(sw.get_input(2), 0.0);
        assert_eq!(sw.get_output(0), 0.0);
    }

    #[test]
    fn test_relay_reset() {
        let mut relay = Relay::new(0.5, 1.0);
        relay.set_input(0, 1.5);
        relay.update(0.0);
        relay.reset();
        assert_eq!(relay.get_input(0), 0.0);
        assert_eq!(relay.get_output(0), 0.0);
    }

    #[test]
    #[should_panic(expected = "off_threshold must be < on_threshold")]
    fn test_relay_invalid_thresholds() {
        let _relay = Relay::new(1.0, 0.5);
    }

    #[test]
    #[should_panic(expected = "hysteresis must be non-negative")]
    fn test_comparator_negative_hysteresis() {
        let _comp = Comparator::with_hysteresis(0.5, -0.1);
    }
}
