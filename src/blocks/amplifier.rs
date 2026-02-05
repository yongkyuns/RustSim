//! Scalar amplifier block: y = gain * u

use crate::block::{AlgebraicBlock, Block, StepResult};

/// Scalar amplifier: y = gain * u
///
/// # Example
///
/// ```ignore
/// let mut amp = Amplifier::new(2.0);
/// amp.set_input(0, 3.0);
/// amp.update(0.0);
/// assert_eq!(amp.get_output(0), 6.0);
/// ```
#[derive(Debug, Clone)]
pub struct Amplifier {
    input: f64,
    output: f64,
    gain: f64,
}

impl Amplifier {
    /// Create amplifier with given gain
    pub fn new(gain: f64) -> Self {
        Self {
            input: 0.0,
            output: 0.0,
            gain,
        }
    }

    /// Get current gain
    pub fn gain(&self) -> f64 {
        self.gain
    }

    /// Set gain
    pub fn set_gain(&mut self, gain: f64) {
        self.gain = gain;
    }
}

impl Block for Amplifier {
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
        self.output = self.gain * self.input;
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = 0.0;
    }
}

impl AlgebraicBlock for Amplifier {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_amplifier_basic() {
        let mut amp = Amplifier::new(2.0);
        amp.set_input(0, 3.0);
        amp.update(0.0);
        assert_eq!(amp.get_output(0), 6.0);
    }

    #[test]
    fn test_amplifier_negative_gain() {
        let mut amp = Amplifier::new(-1.5);
        amp.set_input(0, 4.0);
        amp.update(0.0);
        assert_eq!(amp.get_output(0), -6.0);
    }

    #[test]
    fn test_amplifier_reset() {
        let mut amp = Amplifier::new(2.0);
        amp.set_input(0, 5.0);
        amp.update(0.0);
        amp.reset();
        assert_eq!(amp.get_input(0), 0.0);
        assert_eq!(amp.get_output(0), 0.0);
    }
}
