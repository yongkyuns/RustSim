//! Noise source blocks for stochastic simulations
//!
//! This module requires the `rand-support` feature.

use crate::block::{AlgebraicBlock, Block};
use rand::distributions::Distribution;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_distr::{Normal, Uniform};

/// Gaussian white noise source
///
/// Generates normally distributed random values with configurable
/// standard deviation. Supports reproducible sequences via optional seed.
///
/// # Example
///
/// ```ignore
/// let mut noise = WhiteNoise::new(1.0, Some(42));
/// noise.update(0.0);
/// let sample = noise.get_output(0);
/// ```
#[derive(Debug, Clone)]
pub struct WhiteNoise {
    output: f64,
    std_dev: f64,
    rng: StdRng,
    distribution: Normal<f64>,
}

impl WhiteNoise {
    /// Create white noise with standard deviation and optional seed
    pub fn new(std_dev: f64, seed: Option<u64>) -> Self {
        assert!(std_dev >= 0.0, "Standard deviation must be non-negative");

        let rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        let distribution = Normal::new(0.0, std_dev).unwrap();

        Self {
            output: 0.0,
            std_dev,
            rng,
            distribution,
        }
    }

    /// Create unit variance white noise (std_dev = 1.0)
    pub fn unit(seed: Option<u64>) -> Self {
        Self::new(1.0, seed)
    }

    /// Get current standard deviation
    pub fn std_dev(&self) -> f64 {
        self.std_dev
    }

    /// Set new standard deviation
    pub fn set_std_dev(&mut self, std_dev: f64) {
        assert!(std_dev >= 0.0, "Standard deviation must be non-negative");
        self.std_dev = std_dev;
        self.distribution = Normal::new(0.0, std_dev).unwrap();
    }

    /// Reset with new seed
    pub fn reseed(&mut self, seed: u64) {
        self.rng = StdRng::seed_from_u64(seed);
    }
}

impl Block for WhiteNoise {
    const NUM_INPUTS: usize = 0;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        &[]
    }
    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut []
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, _t: f64) {
        self.output = self.distribution.sample(&mut self.rng);
    }

    fn reset(&mut self) {
        // Reset does not reseed - preserves current RNG state
        self.output = 0.0;
    }
}

impl AlgebraicBlock for WhiteNoise {}

/// Uniform random noise source
///
/// Generates uniformly distributed random values in [min, max].
/// Supports reproducible sequences via optional seed.
///
/// # Example
///
/// ```ignore
/// let mut noise = UniformNoise::new(-1.0, 1.0, Some(42));
/// noise.update(0.0);
/// let sample = noise.get_output(0);
/// ```
#[derive(Debug, Clone)]
pub struct UniformNoise {
    output: f64,
    min: f64,
    max: f64,
    rng: StdRng,
    distribution: Uniform<f64>,
}

impl UniformNoise {
    /// Create uniform noise in range [min, max] with optional seed
    pub fn new(min: f64, max: f64, seed: Option<u64>) -> Self {
        assert!(min < max, "min must be less than max");

        let rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        let distribution = Uniform::new(min, max);

        Self {
            output: min,
            min,
            max,
            rng,
            distribution,
        }
    }

    /// Create uniform noise in [0, 1]
    pub fn unit(seed: Option<u64>) -> Self {
        Self::new(0.0, 1.0, seed)
    }

    /// Get range bounds
    pub fn range(&self) -> (f64, f64) {
        (self.min, self.max)
    }

    /// Set new range
    pub fn set_range(&mut self, min: f64, max: f64) {
        assert!(min < max, "min must be less than max");
        self.min = min;
        self.max = max;
        self.distribution = Uniform::new(min, max);
    }

    /// Reset with new seed
    pub fn reseed(&mut self, seed: u64) {
        self.rng = StdRng::seed_from_u64(seed);
    }
}

impl Block for UniformNoise {
    const NUM_INPUTS: usize = 0;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        &[]
    }
    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut []
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, _t: f64) {
        self.output = self.distribution.sample(&mut self.rng);
    }

    fn reset(&mut self) {
        // Reset does not reseed - preserves current RNG state
        self.output = self.min;
    }
}

impl AlgebraicBlock for UniformNoise {}

/// Pink noise (1/f noise) source
///
/// Generates pink noise using the Voss-McCartney algorithm with
/// a fixed-size filter bank. Pink noise has equal power per octave.
///
/// # Example
///
/// ```ignore
/// let mut noise = PinkNoise::new(1.0, Some(42));
/// noise.update(0.0);
/// let sample = noise.get_output(0);
/// ```
#[derive(Debug, Clone)]
pub struct PinkNoise {
    output: f64,
    amplitude: f64,
    rng: StdRng,
    white_distribution: Uniform<f64>,
    // Voss-McCartney algorithm state
    rows: [f64; 16],
    running_sum: f64,
    max_key: usize,
}

impl PinkNoise {
    /// Create pink noise with amplitude scaling and optional seed
    pub fn new(amplitude: f64, seed: Option<u64>) -> Self {
        let rng = match seed {
            Some(s) => StdRng::seed_from_u64(s),
            None => StdRng::from_entropy(),
        };

        let white_distribution = Uniform::new(-1.0, 1.0);

        let mut noise = Self {
            output: 0.0,
            amplitude,
            rng,
            white_distribution,
            rows: [0.0; 16],
            running_sum: 0.0,
            max_key: 0,
        };

        // Initialize all rows with random values
        for i in 0..16 {
            noise.rows[i] = noise.white_distribution.sample(&mut noise.rng);
            noise.running_sum += noise.rows[i];
        }

        noise
    }

    /// Create unit amplitude pink noise
    pub fn unit(seed: Option<u64>) -> Self {
        Self::new(1.0, seed)
    }

    /// Get current amplitude
    pub fn amplitude(&self) -> f64 {
        self.amplitude
    }

    /// Set amplitude
    pub fn set_amplitude(&mut self, amplitude: f64) {
        self.amplitude = amplitude;
    }

    /// Reset with new seed
    pub fn reseed(&mut self, seed: u64) {
        self.rng = StdRng::seed_from_u64(seed);
        // Re-initialize state
        self.running_sum = 0.0;
        for i in 0..16 {
            self.rows[i] = self.white_distribution.sample(&mut self.rng);
            self.running_sum += self.rows[i];
        }
        self.max_key = 0;
    }
}

impl Block for PinkNoise {
    const NUM_INPUTS: usize = 0;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        &[]
    }
    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut []
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, _t: f64) {
        // Voss-McCartney algorithm
        // Count trailing zeros to determine which rows to update
        let last_key = self.max_key;
        self.max_key = self.max_key.wrapping_add(1);

        // XOR to find changed bits
        let diff = last_key ^ self.max_key;

        // Update rows corresponding to changed bits
        for i in 0..16 {
            if (diff & (1 << i)) != 0 {
                self.running_sum -= self.rows[i];
                self.rows[i] = self.white_distribution.sample(&mut self.rng);
                self.running_sum += self.rows[i];
            }
        }

        // Scale output (divide by number of rows for normalization)
        self.output = self.amplitude * self.running_sum / 16.0;
    }

    fn reset(&mut self) {
        // Reset does not reseed - preserves current RNG state
        self.output = 0.0;
    }
}

impl AlgebraicBlock for PinkNoise {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_white_noise_deterministic() {
        let mut n1 = WhiteNoise::new(1.0, Some(42));
        let mut n2 = WhiteNoise::new(1.0, Some(42));

        // With same seed, should produce identical sequences
        for _ in 0..10 {
            n1.update(0.0);
            n2.update(0.0);
            assert_eq!(n1.get_output(0), n2.get_output(0));
        }
    }

    #[test]
    fn test_white_noise_statistics() {
        let mut noise = WhiteNoise::new(1.0, Some(12345));
        let n = 10000;
        let mut sum = 0.0;
        let mut sum_sq = 0.0;

        for _ in 0..n {
            noise.update(0.0);
            let x = noise.get_output(0);
            sum += x;
            sum_sq += x * x;
        }

        let mean = sum / n as f64;
        let variance = sum_sq / n as f64 - mean * mean;

        // Mean should be close to 0
        assert!(mean.abs() < 0.05, "Mean = {}, expected ~0", mean);

        // Variance should be close to 1
        assert!(
            (variance - 1.0).abs() < 0.1,
            "Variance = {}, expected ~1",
            variance
        );
    }

    #[test]
    fn test_white_noise_set_std_dev() {
        let mut noise = WhiteNoise::new(1.0, Some(42));
        assert_eq!(noise.std_dev(), 1.0);

        noise.set_std_dev(2.0);
        assert_eq!(noise.std_dev(), 2.0);
    }

    #[test]
    fn test_uniform_noise_deterministic() {
        let mut n1 = UniformNoise::new(-1.0, 1.0, Some(42));
        let mut n2 = UniformNoise::new(-1.0, 1.0, Some(42));

        // With same seed, should produce identical sequences
        for _ in 0..10 {
            n1.update(0.0);
            n2.update(0.0);
            assert_eq!(n1.get_output(0), n2.get_output(0));
        }
    }

    #[test]
    fn test_uniform_noise_range() {
        let mut noise = UniformNoise::new(-2.0, 3.0, Some(12345));

        for _ in 0..100 {
            noise.update(0.0);
            let x = noise.get_output(0);
            assert!(x >= -2.0 && x < 3.0, "Value {} out of range [-2, 3)", x);
        }
    }

    #[test]
    fn test_uniform_noise_statistics() {
        let mut noise = UniformNoise::new(0.0, 10.0, Some(12345));
        let n = 10000;
        let mut sum = 0.0;

        for _ in 0..n {
            noise.update(0.0);
            sum += noise.get_output(0);
        }

        let mean = sum / n as f64;

        // Mean of uniform [0, 10] should be 5
        assert!((mean - 5.0).abs() < 0.1, "Mean = {}, expected ~5", mean);
    }

    #[test]
    fn test_uniform_noise_set_range() {
        let mut noise = UniformNoise::new(0.0, 1.0, Some(42));
        assert_eq!(noise.range(), (0.0, 1.0));

        noise.set_range(-5.0, 5.0);
        assert_eq!(noise.range(), (-5.0, 5.0));

        // Verify new range is enforced
        for _ in 0..100 {
            noise.update(0.0);
            let x = noise.get_output(0);
            assert!(x >= -5.0 && x < 5.0);
        }
    }

    #[test]
    fn test_pink_noise_deterministic() {
        let mut n1 = PinkNoise::new(1.0, Some(42));
        let mut n2 = PinkNoise::new(1.0, Some(42));

        // With same seed, should produce identical sequences
        for _ in 0..10 {
            n1.update(0.0);
            n2.update(0.0);
            assert_eq!(n1.get_output(0), n2.get_output(0));
        }
    }

    #[test]
    fn test_pink_noise_bounded() {
        let mut noise = PinkNoise::new(1.0, Some(12345));

        // Pink noise should generally stay within reasonable bounds
        // (approximately -amplitude to +amplitude)
        for _ in 0..1000 {
            noise.update(0.0);
            let x = noise.get_output(0);
            // With 16 averaging terms, should be well bounded
            assert!(
                x.abs() < 5.0,
                "Pink noise value {} exceeds expected bounds",
                x
            );
        }
    }

    #[test]
    fn test_pink_noise_set_amplitude() {
        let mut noise = PinkNoise::new(1.0, Some(42));
        assert_eq!(noise.amplitude(), 1.0);

        noise.set_amplitude(2.0);
        assert_eq!(noise.amplitude(), 2.0);
    }

    #[test]
    fn test_pink_noise_reseed() {
        let mut n1 = PinkNoise::new(1.0, Some(42));
        let mut n2 = PinkNoise::new(1.0, Some(99));

        // Generate some values
        for _ in 0..5 {
            n1.update(0.0);
            n2.update(0.0);
        }

        // Reseed n2 with same seed as n1's original
        n2.reseed(42);
        n1 = PinkNoise::new(1.0, Some(42));

        // Should now produce same sequence
        for _ in 0..10 {
            n1.update(0.0);
            n2.update(0.0);
            assert_eq!(n1.get_output(0), n2.get_output(0));
        }
    }

    #[test]
    fn test_noise_blocks_have_zero_inputs() {
        assert_eq!(WhiteNoise::NUM_INPUTS, 0);
        assert_eq!(UniformNoise::NUM_INPUTS, 0);
        assert_eq!(PinkNoise::NUM_INPUTS, 0);
    }

    #[test]
    fn test_noise_blocks_are_algebraic() {
        assert!(!WhiteNoise::IS_DYNAMIC);
        assert!(!UniformNoise::IS_DYNAMIC);
        assert!(!PinkNoise::IS_DYNAMIC);
    }
}
