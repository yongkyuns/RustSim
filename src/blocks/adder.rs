//! N-input adder block with const generics

use crate::block::{AlgebraicBlock, Block};

/// N-input weighted adder: y = sum(weights[i] * inputs[i])
///
/// The number of inputs N is fixed at compile time via const generics.
///
/// # Example
///
/// ```ignore
/// // 3-input adder with default weights (all 1.0)
/// let mut adder = Adder::<3>::new();
/// adder.inputs_mut()[0] = 1.0;
/// adder.inputs_mut()[1] = 2.0;
/// adder.inputs_mut()[2] = 3.0;
/// adder.update(0.0);
/// assert_eq!(adder.get_output(0), 6.0);
///
/// // 2-input subtractor
/// let mut sub = Adder::<2>::with_weights([1.0, -1.0]);
/// ```
#[derive(Debug, Clone)]
pub struct Adder<const N: usize> {
    inputs: [f64; N],
    output: f64,
    weights: [f64; N],
}

impl<const N: usize> Adder<N> {
    /// Create adder with all weights = 1.0
    pub fn new() -> Self {
        Self {
            inputs: [0.0; N],
            output: 0.0,
            weights: [1.0; N],
        }
    }

    /// Create adder with specified weights
    pub fn with_weights(weights: [f64; N]) -> Self {
        Self {
            inputs: [0.0; N],
            output: 0.0,
            weights,
        }
    }

    /// Create subtractor: first input positive, rest negative
    pub fn subtractor() -> Self {
        let mut weights = [-1.0; N];
        if N > 0 {
            weights[0] = 1.0;
        }
        Self::with_weights(weights)
    }

    /// Get weights
    pub fn weights(&self) -> &[f64; N] {
        &self.weights
    }

    /// Set weights
    pub fn set_weights(&mut self, weights: [f64; N]) {
        self.weights = weights;
    }
}

impl<const N: usize> Default for Adder<N> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const N: usize> Block for Adder<N> {
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
        self.output = self
            .inputs
            .iter()
            .zip(self.weights.iter())
            .map(|(&x, &w)| x * w)
            .sum();
    }

    fn reset(&mut self) {
        self.inputs.fill(0.0);
        self.output = 0.0;
    }
}

impl<const N: usize> AlgebraicBlock for Adder<N> {}

// Common type aliases
#[allow(dead_code)]
pub type Adder2 = Adder<2>;
#[allow(dead_code)]
pub type Adder3 = Adder<3>;
#[allow(dead_code)]
pub type Adder4 = Adder<4>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_adder_basic() {
        let mut adder = Adder::<3>::new();
        adder.inputs_mut()[0] = 1.0;
        adder.inputs_mut()[1] = 2.0;
        adder.inputs_mut()[2] = 3.0;
        adder.update(0.0);
        assert_eq!(adder.get_output(0), 6.0);
    }

    #[test]
    fn test_adder_weighted() {
        let mut adder = Adder::<2>::with_weights([2.0, 3.0]);
        adder.inputs_mut()[0] = 1.0;
        adder.inputs_mut()[1] = 2.0;
        adder.update(0.0);
        assert_eq!(adder.get_output(0), 8.0); // 2*1 + 3*2
    }

    #[test]
    fn test_subtractor() {
        let mut sub = Adder::<2>::subtractor();
        sub.inputs_mut()[0] = 5.0;
        sub.inputs_mut()[1] = 3.0;
        sub.update(0.0);
        assert_eq!(sub.get_output(0), 2.0); // 5 - 3
    }

    #[test]
    fn test_adder_single() {
        let mut adder = Adder::<1>::new();
        adder.set_input(0, 42.0);
        adder.update(0.0);
        assert_eq!(adder.get_output(0), 42.0);
    }
}
