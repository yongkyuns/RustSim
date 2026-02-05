//! N-input multiplier block with const generics

use crate::block::{AlgebraicBlock, Block};

/// N-input multiplier: y = product(inputs[i])
///
/// The number of inputs N is fixed at compile time via const generics.
///
/// # Example
///
/// ```ignore
/// // 3-input multiplier
/// let mut mult = Multiplier::<3>::new();
/// mult.inputs_mut()[0] = 2.0;
/// mult.inputs_mut()[1] = 3.0;
/// mult.inputs_mut()[2] = 4.0;
/// mult.update(0.0);
/// assert_eq!(mult.get_output(0), 24.0);
/// ```
#[derive(Debug, Clone)]
pub struct Multiplier<const N: usize> {
    inputs: [f64; N],
    output: f64,
}

impl<const N: usize> Multiplier<N> {
    /// Create multiplier
    pub fn new() -> Self {
        Self {
            inputs: [0.0; N],
            output: 0.0,
        }
    }
}

impl<const N: usize> Default for Multiplier<N> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const N: usize> Block for Multiplier<N> {
    const NUM_INPUTS: usize = N;
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
        self.output = self.inputs.iter().product();
    }

    fn reset(&mut self) {
        self.inputs.fill(0.0);
        self.output = 0.0;
    }
}

impl<const N: usize> AlgebraicBlock for Multiplier<N> {}

// Common type aliases
pub type Multiplier2 = Multiplier<2>;
pub type Multiplier3 = Multiplier<3>;
pub type Multiplier4 = Multiplier<4>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multiplier_basic() {
        let mut mult = Multiplier::<3>::new();
        mult.inputs_mut()[0] = 2.0;
        mult.inputs_mut()[1] = 3.0;
        mult.inputs_mut()[2] = 4.0;
        mult.update(0.0);
        assert_eq!(mult.get_output(0), 24.0);
    }

    #[test]
    fn test_multiplier_two_inputs() {
        let mut mult = Multiplier::<2>::new();
        mult.inputs_mut()[0] = 5.0;
        mult.inputs_mut()[1] = 7.0;
        mult.update(0.0);
        assert_eq!(mult.get_output(0), 35.0);
    }

    #[test]
    fn test_multiplier_with_negative() {
        let mut mult = Multiplier::<2>::new();
        mult.inputs_mut()[0] = -3.0;
        mult.inputs_mut()[1] = 4.0;
        mult.update(0.0);
        assert_eq!(mult.get_output(0), -12.0);
    }

    #[test]
    fn test_multiplier_with_zero() {
        let mut mult = Multiplier::<3>::new();
        mult.inputs_mut()[0] = 5.0;
        mult.inputs_mut()[1] = 0.0;
        mult.inputs_mut()[2] = 7.0;
        mult.update(0.0);
        assert_eq!(mult.get_output(0), 0.0);
    }

    #[test]
    fn test_multiplier_single() {
        let mut mult = Multiplier::<1>::new();
        mult.set_input(0, 42.0);
        mult.update(0.0);
        assert_eq!(mult.get_output(0), 42.0);
    }

    #[test]
    fn test_multiplier_reset() {
        let mut mult = Multiplier::<2>::new();
        mult.inputs_mut()[0] = 3.0;
        mult.inputs_mut()[1] = 4.0;
        mult.update(0.0);
        mult.reset();
        assert_eq!(mult.get_input(0), 0.0);
        assert_eq!(mult.get_input(1), 0.0);
        assert_eq!(mult.get_output(0), 0.0);
    }

    #[test]
    fn test_multiplier_fractional() {
        let mut mult = Multiplier::<2>::new();
        mult.inputs_mut()[0] = 0.5;
        mult.inputs_mut()[1] = 0.25;
        mult.update(0.0);
        assert_eq!(mult.get_output(0), 0.125);
    }
}
