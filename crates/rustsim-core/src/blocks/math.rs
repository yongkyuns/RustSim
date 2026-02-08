//! Mathematical function blocks (1→1 operations)
//!
//! This module provides unary mathematical operations like sin, cos, exp, sqrt, etc.
//! All blocks are 1-input, 1-output algebraic operations.

use crate::block::{AlgebraicBlock, Block};

/// Macro to generate simple 1→1 mathematical function blocks
///
/// Reduces boilerplate for blocks that apply a single mathematical operation.
macro_rules! math_block {
    ($name:ident, $doc:expr, $op:expr) => {
        #[doc = $doc]
        #[derive(Debug, Clone, Copy, Default)]
        pub struct $name {
            input: f64,
            output: f64,
        }

        impl $name {
            /// Create a new instance
            pub fn new() -> Self {
                Self {
                    input: 0.0,
                    output: 0.0,
                }
            }
        }

        impl Block for $name {
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
                self.output = $op(self.input);
            }

            fn reset(&mut self) {
                self.input = 0.0;
                self.output = 0.0;
            }
        }

        impl AlgebraicBlock for $name {}
    };
}

// Trigonometric functions
math_block!(Sin, "Sine: y = sin(u)", f64::sin);
math_block!(Cos, "Cosine: y = cos(u)", f64::cos);
math_block!(Tan, "Tangent: y = tan(u)", f64::tan);

// Inverse trigonometric functions
math_block!(Asin, "Arcsine: y = asin(u)", f64::asin);
math_block!(Acos, "Arccosine: y = acos(u)", f64::acos);
math_block!(Atan, "Arctangent: y = atan(u)", f64::atan);

// Hyperbolic functions
math_block!(Sinh, "Hyperbolic sine: y = sinh(u)", f64::sinh);
math_block!(Cosh, "Hyperbolic cosine: y = cosh(u)", f64::cosh);
math_block!(Tanh, "Hyperbolic tangent: y = tanh(u)", f64::tanh);

// Exponential and logarithmic functions
math_block!(Exp, "Exponential: y = e^u", f64::exp);
math_block!(Log, "Natural logarithm: y = ln(u)", f64::ln);
math_block!(Log10, "Base-10 logarithm: y = log10(u)", f64::log10);

// Root and power functions
math_block!(Sqrt, "Square root: y = sqrt(u)", f64::sqrt);
math_block!(Cbrt, "Cube root: y = cbrt(u)", f64::cbrt);

// Sign and rounding functions
math_block!(Abs, "Absolute value: y = |u|", f64::abs);
math_block!(Floor, "Floor: y = floor(u)", f64::floor);
math_block!(Ceil, "Ceiling: y = ceil(u)", f64::ceil);
math_block!(Round, "Round: y = round(u)", f64::round);

/// Sign function: y = sign(u)
///
/// Returns:
/// - 1.0 if u > 0
/// - 0.0 if u = 0
/// - -1.0 if u < 0
#[derive(Debug, Clone, Copy, Default)]
pub struct Sign {
    input: f64,
    output: f64,
}

impl Sign {
    pub fn new() -> Self {
        Self {
            input: 0.0,
            output: 0.0,
        }
    }
}

impl Block for Sign {
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
        self.output = if self.input > 0.0 {
            1.0
        } else if self.input < 0.0 {
            -1.0
        } else {
            0.0
        };
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = 0.0;
    }
}

impl AlgebraicBlock for Sign {}

/// Power function: y = u^n
///
/// Raises input to a configurable exponent.
#[derive(Debug, Clone, Copy)]
pub struct Pow {
    input: f64,
    output: f64,
    exponent: f64,
}

impl Pow {
    /// Create power block with given exponent
    pub fn new(exponent: f64) -> Self {
        Self {
            input: 0.0,
            output: 0.0,
            exponent,
        }
    }

    /// Get current exponent
    pub fn exponent(&self) -> f64 {
        self.exponent
    }

    /// Set exponent
    pub fn set_exponent(&mut self, exponent: f64) {
        self.exponent = exponent;
    }
}

impl Block for Pow {
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
        self.output = self.input.powf(self.exponent);
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = 0.0;
    }
}

impl AlgebraicBlock for Pow {}

/// Saturation/clipping block: y = clamp(u, min, max)
///
/// Limits the input to the range [min, max].
#[derive(Debug, Clone, Copy)]
pub struct Clip {
    input: f64,
    output: f64,
    min: f64,
    max: f64,
}

impl Clip {
    /// Create saturation block with given limits
    pub fn new(min: f64, max: f64) -> Self {
        assert!(min <= max, "min must be <= max");
        Self {
            input: 0.0,
            output: 0.0,
            min,
            max,
        }
    }

    /// Get min limit
    pub fn min(&self) -> f64 {
        self.min
    }

    /// Get max limit
    pub fn max(&self) -> f64 {
        self.max
    }

    /// Set limits
    pub fn set_limits(&mut self, min: f64, max: f64) {
        assert!(min <= max, "min must be <= max");
        self.min = min;
        self.max = max;
    }
}

impl Block for Clip {
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
        self.output = self.input.clamp(self.min, self.max);
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = 0.0;
    }
}

impl AlgebraicBlock for Clip {}

/// Modulo operation: y = u mod m
///
/// Computes the modulo of the input with respect to a configurable modulus.
#[derive(Debug, Clone, Copy)]
pub struct Mod {
    input: f64,
    output: f64,
    modulus: f64,
}

impl Mod {
    /// Create modulo block with given modulus
    pub fn new(modulus: f64) -> Self {
        Self {
            input: 0.0,
            output: 0.0,
            modulus,
        }
    }

    /// Get current modulus
    pub fn modulus(&self) -> f64 {
        self.modulus
    }

    /// Set modulus
    pub fn set_modulus(&mut self, modulus: f64) {
        self.modulus = modulus;
    }
}

impl Default for Mod {
    fn default() -> Self {
        Self::new(1.0)
    }
}

impl Block for Mod {
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
        self.output = self.input.rem_euclid(self.modulus);
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = 0.0;
    }
}

impl AlgebraicBlock for Mod {}

/// Vector L2 norm: y = sqrt(sum(u_i^2))
///
/// Computes the Euclidean norm of N inputs.
#[derive(Debug, Clone)]
pub struct Norm<const N: usize> {
    inputs: [f64; N],
    output: f64,
}

impl<const N: usize> Norm<N> {
    /// Create a new norm block
    pub fn new() -> Self {
        Self {
            inputs: [0.0; N],
            output: 0.0,
        }
    }
}

impl<const N: usize> Default for Norm<N> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const N: usize> Block for Norm<N> {
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
        let sum_squares: f64 = self.inputs.iter().map(|&x| x * x).sum();
        self.output = sum_squares.sqrt();
    }

    fn reset(&mut self) {
        self.inputs.fill(0.0);
        self.output = 0.0;
    }
}

impl<const N: usize> AlgebraicBlock for Norm<N> {}

// Common type aliases for Norm
pub type Norm2 = Norm<2>;
pub type Norm3 = Norm<3>;

/// Minimum of N inputs
///
/// Computes the minimum value among all inputs.
#[derive(Debug, Clone)]
pub struct Min<const N: usize> {
    inputs: [f64; N],
    output: f64,
}

impl<const N: usize> Min<N> {
    /// Create a new min block
    pub fn new() -> Self {
        Self {
            inputs: [0.0; N],
            output: 0.0,
        }
    }
}

impl<const N: usize> Default for Min<N> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const N: usize> Block for Min<N> {
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
        self.output = self.inputs.iter().copied().fold(f64::INFINITY, f64::min);
    }

    fn reset(&mut self) {
        self.inputs.fill(0.0);
        self.output = 0.0;
    }
}

impl<const N: usize> AlgebraicBlock for Min<N> {}

// Common type aliases for Min
pub type Min2 = Min<2>;
pub type Min3 = Min<3>;

/// Maximum of N inputs
///
/// Computes the maximum value among all inputs.
#[derive(Debug, Clone)]
pub struct Max<const N: usize> {
    inputs: [f64; N],
    output: f64,
}

impl<const N: usize> Max<N> {
    /// Create a new max block
    pub fn new() -> Self {
        Self {
            inputs: [0.0; N],
            output: 0.0,
        }
    }
}

impl<const N: usize> Default for Max<N> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const N: usize> Block for Max<N> {
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
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
    }

    fn reset(&mut self) {
        self.inputs.fill(0.0);
        self.output = 0.0;
    }
}

impl<const N: usize> AlgebraicBlock for Max<N> {}

// Common type aliases for Max
pub type Max2 = Max<2>;
pub type Max3 = Max<3>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sin() {
        let mut block = Sin::new();
        block.set_input(0, std::f64::consts::PI / 2.0);
        block.update(0.0);
        assert!((block.get_output(0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cos() {
        let mut block = Cos::new();
        block.set_input(0, std::f64::consts::PI);
        block.update(0.0);
        assert!((block.get_output(0) - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_tan() {
        let mut block = Tan::new();
        block.set_input(0, std::f64::consts::PI / 4.0);
        block.update(0.0);
        assert!((block.get_output(0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_exp() {
        let mut block = Exp::new();
        block.set_input(0, 1.0);
        block.update(0.0);
        assert!((block.get_output(0) - std::f64::consts::E).abs() < 1e-10);
    }

    #[test]
    fn test_log() {
        let mut block = Log::new();
        block.set_input(0, std::f64::consts::E);
        block.update(0.0);
        assert!((block.get_output(0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_sqrt() {
        let mut block = Sqrt::new();
        block.set_input(0, 16.0);
        block.update(0.0);
        assert_eq!(block.get_output(0), 4.0);
    }

    #[test]
    fn test_abs() {
        let mut block = Abs::new();
        block.set_input(0, -5.5);
        block.update(0.0);
        assert_eq!(block.get_output(0), 5.5);
    }

    #[test]
    fn test_sign() {
        let mut block = Sign::new();

        block.set_input(0, 3.5);
        block.update(0.0);
        assert_eq!(block.get_output(0), 1.0);

        block.set_input(0, -2.1);
        block.update(0.0);
        assert_eq!(block.get_output(0), -1.0);

        block.set_input(0, 0.0);
        block.update(0.0);
        assert_eq!(block.get_output(0), 0.0);
    }

    #[test]
    fn test_pow() {
        let mut block = Pow::new(3.0);
        block.set_input(0, 2.0);
        block.update(0.0);
        assert_eq!(block.get_output(0), 8.0);

        block.set_exponent(0.5);
        block.set_input(0, 16.0);
        block.update(0.0);
        assert_eq!(block.get_output(0), 4.0);
    }

    #[test]
    fn test_clip() {
        let mut block = Clip::new(-1.0, 1.0);

        block.set_input(0, 0.5);
        block.update(0.0);
        assert_eq!(block.get_output(0), 0.5);

        block.set_input(0, 2.0);
        block.update(0.0);
        assert_eq!(block.get_output(0), 1.0);

        block.set_input(0, -3.0);
        block.update(0.0);
        assert_eq!(block.get_output(0), -1.0);
    }

    #[test]
    fn test_floor_ceil_round() {
        let mut floor = Floor::new();
        floor.set_input(0, 2.7);
        floor.update(0.0);
        assert_eq!(floor.get_output(0), 2.0);

        let mut ceil = Ceil::new();
        ceil.set_input(0, 2.3);
        ceil.update(0.0);
        assert_eq!(ceil.get_output(0), 3.0);

        let mut round = Round::new();
        round.set_input(0, 2.5);
        round.update(0.0);
        assert_eq!(round.get_output(0), 3.0);
    }

    #[test]
    fn test_hyperbolic() {
        let mut sinh = Sinh::new();
        sinh.set_input(0, 0.0);
        sinh.update(0.0);
        assert_eq!(sinh.get_output(0), 0.0);

        let mut cosh = Cosh::new();
        cosh.set_input(0, 0.0);
        cosh.update(0.0);
        assert_eq!(cosh.get_output(0), 1.0);

        let mut tanh = Tanh::new();
        tanh.set_input(0, 0.0);
        tanh.update(0.0);
        assert_eq!(tanh.get_output(0), 0.0);
    }

    #[test]
    fn test_inverse_trig() {
        let mut asin = Asin::new();
        asin.set_input(0, 0.5);
        asin.update(0.0);
        assert!((asin.get_output(0) - std::f64::consts::PI / 6.0).abs() < 1e-10);

        let mut acos = Acos::new();
        acos.set_input(0, 0.5);
        acos.update(0.0);
        assert!((acos.get_output(0) - std::f64::consts::PI / 3.0).abs() < 1e-10);

        let mut atan = Atan::new();
        atan.set_input(0, 1.0);
        atan.update(0.0);
        assert!((atan.get_output(0) - std::f64::consts::PI / 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_reset() {
        let mut block = Sin::new();
        block.set_input(0, 1.0);
        block.update(0.0);
        block.reset();
        assert_eq!(block.get_input(0), 0.0);
        assert_eq!(block.get_output(0), 0.0);
    }

    // Log10 tests
    #[test]
    fn test_log10_basic() {
        let mut block = Log10::new();
        block.set_input(0, 100.0);
        block.update(0.0);
        assert!((block.get_output(0) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_log10_powers_of_10() {
        let mut block = Log10::new();

        // Test various powers of 10
        let test_cases = vec![
            (1.0, 0.0),
            (10.0, 1.0),
            (100.0, 2.0),
            (1000.0, 3.0),
            (0.1, -1.0),
            (0.01, -2.0),
        ];

        for (input, expected) in test_cases {
            block.set_input(0, input);
            block.update(0.0);
            assert!(
                (block.get_output(0) - expected).abs() < 1e-10,
                "log10({}) = {}, expected {}",
                input,
                block.get_output(0),
                expected
            );
        }
    }

    // Norm tests
    #[test]
    fn test_norm_single() {
        let mut block = Norm::<1>::new();
        block.set_input(0, 5.0);
        block.update(0.0);
        assert!((block.get_output(0) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_norm_vector() {
        let mut block = Norm::<2>::new();
        block.set_input(0, 3.0);
        block.set_input(1, 4.0);
        block.update(0.0);
        assert!((block.get_output(0) - 5.0).abs() < 1e-10); // sqrt(3^2 + 4^2) = 5
    }

    #[test]
    fn test_norm_3d() {
        let mut block = Norm::<3>::new();
        block.set_input(0, 1.0);
        block.set_input(1, 2.0);
        block.set_input(2, 2.0);
        block.update(0.0);
        let expected = (1.0_f64 + 4.0 + 4.0).sqrt(); // sqrt(1 + 4 + 4) = 3
        assert!((block.get_output(0) - expected).abs() < 1e-10);
    }

    // Mod tests
    #[test]
    fn test_mod_basic() {
        let mut block = Mod::new(1.0);
        block.set_input(0, 2.5);
        block.update(0.0);
        assert!((block.get_output(0) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_mod_configurable_modulus() {
        let mut block = Mod::new(3.0);
        block.set_input(0, 7.0);
        block.update(0.0);
        assert!((block.get_output(0) - 1.0).abs() < 1e-10);

        block.set_modulus(5.0);
        block.set_input(0, 12.0);
        block.update(0.0);
        assert!((block.get_output(0) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_mod_negative() {
        let mut block = Mod::new(2.0);
        block.set_input(0, -3.0);
        block.update(0.0);
        assert!((block.get_output(0) - 1.0).abs() < 1e-10); // -3 mod 2 = 1 (Euclidean)
    }

    // Min tests
    #[test]
    fn test_min_two_inputs() {
        let mut block = Min::<2>::new();
        block.set_input(0, 5.0);
        block.set_input(1, 3.0);
        block.update(0.0);
        assert_eq!(block.get_output(0), 3.0);

        block.set_input(0, 2.0);
        block.set_input(1, 8.0);
        block.update(0.0);
        assert_eq!(block.get_output(0), 2.0);
    }

    #[test]
    fn test_min_three_inputs() {
        let mut block = Min::<3>::new();
        block.set_input(0, 5.0);
        block.set_input(1, 3.0);
        block.set_input(2, 7.0);
        block.update(0.0);
        assert_eq!(block.get_output(0), 3.0);

        block.set_input(0, -1.0);
        block.set_input(1, 3.0);
        block.set_input(2, 0.0);
        block.update(0.0);
        assert_eq!(block.get_output(0), -1.0);
    }

    // Max tests
    #[test]
    fn test_max_two_inputs() {
        let mut block = Max::<2>::new();
        block.set_input(0, 5.0);
        block.set_input(1, 3.0);
        block.update(0.0);
        assert_eq!(block.get_output(0), 5.0);

        block.set_input(0, 2.0);
        block.set_input(1, 8.0);
        block.update(0.0);
        assert_eq!(block.get_output(0), 8.0);
    }

    #[test]
    fn test_max_three_inputs() {
        let mut block = Max::<3>::new();
        block.set_input(0, 5.0);
        block.set_input(1, 3.0);
        block.set_input(2, 7.0);
        block.update(0.0);
        assert_eq!(block.get_output(0), 7.0);

        block.set_input(0, -1.0);
        block.set_input(1, 3.0);
        block.set_input(2, 0.0);
        block.update(0.0);
        assert_eq!(block.get_output(0), 3.0);
    }
}
