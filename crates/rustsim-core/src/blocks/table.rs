//! Lookup table blocks with linear interpolation
//!
//! This module implements:
//! - LUT1D: One-dimensional lookup table with linear interpolation
//! - LUT: Multi-dimensional lookup table with linear interpolation

use crate::block::{AlgebraicBlock, Block};

/// One-dimensional lookup table with linear interpolation
///
/// Performs piecewise linear interpolation on 1D data. Supports both
/// single and multiple outputs. Extrapolates linearly beyond data range.
///
/// # Example
///
/// ```ignore
/// // Single output: y = 2*x
/// let points = vec![0.0, 1.0, 2.0];
/// let values = vec![0.0, 2.0, 4.0];
/// let mut lut = LUT1D::new(points, values);
/// ```
#[derive(Clone, Debug)]
pub struct LUT1D {
    inputs: [f64; 1],
    outputs: Vec<f64>,
    points: Vec<f64>,
    values: Vec<Vec<f64>>, // values[output_idx][point_idx]
    #[allow(dead_code)]
    num_outputs: usize,
}

impl LUT1D {
    /// Create a new 1D lookup table
    ///
    /// # Arguments
    ///
    /// * `points` - Sorted array of x values (must be monotonically increasing)
    /// * `values` - Either a 1D array (single output) or 2D array (multiple outputs).
    ///   For multiple outputs, shape is (num_points, num_outputs)
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - Points and values have different lengths
    /// - Points array is empty
    pub fn new(points: Vec<f64>, values: Vec<f64>) -> Self {
        if points.len() != values.len() {
            panic!("Points and values must have the same length");
        }
        if points.is_empty() {
            panic!("Points array cannot be empty");
        }

        Self {
            inputs: [0.0],
            outputs: vec![0.0],
            points,
            values: vec![values],
            num_outputs: 1,
        }
    }

    /// Create a new 1D lookup table with multiple outputs
    ///
    /// # Arguments
    ///
    /// * `points` - Sorted array of x values (must be monotonically increasing)
    /// * `values` - 2D array where values[i] is the output vector at points[i].
    ///   Shape: (num_points, num_outputs)
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - Points and values have different lengths
    /// - Points array is empty
    /// - Values array is empty or jagged
    pub fn new_multi(points: Vec<f64>, values: Vec<Vec<f64>>) -> Self {
        if points.len() != values.len() {
            panic!("Points and values must have the same length");
        }
        if points.is_empty() || values.is_empty() {
            panic!("Points and values arrays cannot be empty");
        }

        let num_outputs = values[0].len();
        if num_outputs == 0 {
            panic!("Values must have at least one output");
        }

        // Check that all value vectors have the same length
        for v in &values {
            if v.len() != num_outputs {
                panic!("All value vectors must have the same length");
            }
        }

        // Transpose values to values[output_idx][point_idx]
        let mut transposed = vec![vec![0.0; points.len()]; num_outputs];
        for (point_idx, value_vec) in values.iter().enumerate() {
            for (output_idx, &val) in value_vec.iter().enumerate() {
                transposed[output_idx][point_idx] = val;
            }
        }

        Self {
            inputs: [0.0],
            outputs: vec![0.0; num_outputs],
            points,
            values: transposed,
            num_outputs,
        }
    }

    /// Perform linear interpolation/extrapolation
    fn interpolate(&self, x: f64) -> Vec<f64> {
        let n = self.points.len();

        // Handle edge cases
        if n == 1 {
            return self.values.iter().map(|v| v[0]).collect();
        }

        // Find the interval containing x
        let idx = self.points.iter().position(|&p| p >= x);

        let (i0, i1) = match idx {
            None => {
                // x is beyond the last point - extrapolate
                (n - 2, n - 1)
            }
            Some(0) => {
                // x is before the first point - extrapolate
                (0, 1)
            }
            Some(i) => {
                // x is between points
                (i - 1, i)
            }
        };

        let x0 = self.points[i0];
        let x1 = self.points[i1];
        let t = (x - x0) / (x1 - x0);

        // Interpolate each output
        self.values
            .iter()
            .map(|values| {
                let y0 = values[i0];
                let y1 = values[i1];
                y0 + t * (y1 - y0)
            })
            .collect()
    }
}

impl Block for LUT1D {
    const NUM_INPUTS: usize = 1;
    const NUM_OUTPUTS: usize = 0; // Dynamic
    const IS_DYNAMIC: bool = false;

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

    fn update(&mut self, _t: f64) {
        let x = self.inputs[0];
        let interpolated = self.interpolate(x);
        self.outputs.copy_from_slice(&interpolated);
    }

    fn reset(&mut self) {
        self.inputs.fill(0.0);
        self.outputs.fill(0.0);
    }
}

impl AlgebraicBlock for LUT1D {}

/// Multi-dimensional lookup table with linear interpolation
///
/// Uses linear interpolation in N-dimensional space based on a grid of points.
/// For simplicity, this implementation supports 2D lookup tables.
///
/// # Example
///
/// ```ignore
/// // 2D lookup table: z = x + y
/// let points = vec![
///     vec![0.0, 0.0],
///     vec![1.0, 0.0],
///     vec![0.0, 1.0],
///     vec![1.0, 1.0],
/// ];
/// let values = vec![0.0, 1.0, 1.0, 2.0];
/// let mut lut = LUT::new(points, values);
/// ```
#[derive(Clone, Debug)]
pub struct LUT {
    inputs: Vec<f64>,
    outputs: Vec<f64>,
    points: Vec<Vec<f64>>, // points[i] is the i-th point in N-D space
    values: Vec<Vec<f64>>, // values[i] is the output vector at points[i]
    num_inputs: usize,
    num_outputs: usize,
}

impl LUT {
    /// Create a new multi-dimensional lookup table
    ///
    /// # Arguments
    ///
    /// * `points` - Array of N-dimensional points, shape: (num_points, ndim)
    /// * `values` - Either a 1D array (single output) or 2D array (multiple outputs).
    ///   For single output: shape is (num_points,).
    ///   For multiple outputs: shape is (num_points, num_outputs)
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - Points and values have different lengths
    /// - Points array is empty
    pub fn new(points: Vec<Vec<f64>>, values: Vec<f64>) -> Self {
        if points.len() != values.len() {
            panic!("Points and values must have the same length");
        }
        if points.is_empty() {
            panic!("Points array cannot be empty");
        }

        let num_inputs = points[0].len();

        // Check that all points have the same dimension
        for p in &points {
            if p.len() != num_inputs {
                panic!("All points must have the same dimension");
            }
        }

        // Convert single-output values to multi-output format
        let values_2d: Vec<Vec<f64>> = values.iter().map(|&v| vec![v]).collect();

        Self {
            inputs: vec![0.0; num_inputs],
            outputs: vec![0.0],
            points,
            values: values_2d,
            num_inputs,
            num_outputs: 1,
        }
    }

    /// Create a new multi-dimensional lookup table with multiple outputs
    ///
    /// # Arguments
    ///
    /// * `points` - Array of N-dimensional points, shape: (num_points, ndim)
    /// * `values` - 2D array where values[i] is the output vector at points[i].
    ///   Shape: (num_points, num_outputs)
    pub fn new_multi(points: Vec<Vec<f64>>, values: Vec<Vec<f64>>) -> Self {
        if points.len() != values.len() {
            panic!("Points and values must have the same length");
        }
        if points.is_empty() || values.is_empty() {
            panic!("Points and values arrays cannot be empty");
        }

        let num_inputs = points[0].len();
        let num_outputs = values[0].len();

        // Check that all points have the same dimension
        for p in &points {
            if p.len() != num_inputs {
                panic!("All points must have the same dimension");
            }
        }

        // Check that all value vectors have the same length
        for v in &values {
            if v.len() != num_outputs {
                panic!("All value vectors must have the same length");
            }
        }

        Self {
            inputs: vec![0.0; num_inputs],
            outputs: vec![0.0; num_outputs],
            points,
            values,
            num_inputs,
            num_outputs,
        }
    }

    /// Perform linear interpolation in N-D space
    ///
    /// Uses inverse distance weighting for simplicity
    fn interpolate(&self, x: &[f64]) -> Vec<f64> {
        // Find nearest neighbors and compute interpolation weights
        let mut distances: Vec<(usize, f64)> = self
            .points
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let dist_sq: f64 = p
                    .iter()
                    .zip(x.iter())
                    .map(|(pi, xi)| (pi - xi).powi(2))
                    .sum();
                (i, dist_sq.sqrt())
            })
            .collect();

        // Sort by distance
        distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        // If we have an exact match, return it
        if distances[0].1 < 1e-10 {
            return self.values[distances[0].0].clone();
        }

        // Use inverse distance weighting with the nearest points
        // For 2D, use 4 nearest neighbors; for higher dimensions, use more
        let k = (self.num_inputs + 2).max(4).min(self.points.len());
        let neighbors = &distances[..k];

        let mut result = vec![0.0; self.num_outputs];
        let mut weight_sum = 0.0;

        for &(idx, dist) in neighbors {
            let weight = 1.0 / (dist + 1e-10); // Add small epsilon to avoid division by zero
            weight_sum += weight;

            for (out_idx, &val) in self.values[idx].iter().enumerate() {
                result[out_idx] += weight * val;
            }
        }

        // Normalize by total weight
        for val in &mut result {
            *val /= weight_sum;
        }

        result
    }
}

impl Block for LUT {
    const NUM_INPUTS: usize = 0; // Dynamic
    const NUM_OUTPUTS: usize = 0; // Dynamic
    const IS_DYNAMIC: bool = false;

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

    fn update(&mut self, _t: f64) {
        let interpolated = self.interpolate(&self.inputs);
        self.outputs.copy_from_slice(&interpolated);
    }

    fn reset(&mut self) {
        self.inputs.fill(0.0);
        self.outputs.fill(0.0);
    }
}

impl AlgebraicBlock for LUT {}

#[cfg(test)]
mod tests {
    use super::*;

    // ============================================================================
    // LUT1D Tests
    // ============================================================================

    #[test]
    fn test_lut1d_init() {
        // Test initialization with 1D data and single output
        let points = vec![0.0, 1.0, 2.0, 3.0];
        let values = vec![0.0, 1.0, 4.0, 9.0]; // y = x^2

        let lut = LUT1D::new(points, values);

        // Verify initialization
        assert_eq!(lut.outputs().len(), 1);
        assert_eq!(lut.inputs.len(), 1);
        assert_eq!(lut.outputs.len(), 1);
    }

    #[test]
    #[should_panic(expected = "Points and values must have the same length")]
    fn test_lut1d_validation() {
        // Test that validation catches length mismatch
        let points = vec![0.0, 1.0, 2.0];
        let values = vec![0.0, 1.0]; // Wrong length

        let _lut = LUT1D::new(points, values);
    }

    #[test]
    fn test_lut1d_siso_embedding() {
        // Test 1D SISO case (y = 2*x)
        let points = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let values = vec![0.0, 2.0, 4.0, 6.0, 8.0];

        let mut lut = LUT1D::new(points, values);

        // Test exact points
        for t in [0.0, 1.0, 2.0, 3.0, 4.0] {
            lut.set_input(0, t);
            lut.update(0.0);
            let expected = 2.0 * t;
            assert!((lut.get_output(0) - expected).abs() < 1e-6);
        }
    }

    #[test]
    fn test_lut1d_simo_embedding() {
        // Test 1D SIMO case (multiple outputs)
        let points = vec![0.0, 1.0, 2.0];
        let values = vec![vec![0.0, 0.0], vec![1.0, 2.0], vec![4.0, 4.0]];

        let mut lut = LUT1D::new_multi(points, values);

        // Test at t=0
        lut.set_input(0, 0.0);
        lut.update(0.0);
        assert!((lut.outputs()[0] - 0.0).abs() < 1e-6);
        assert!((lut.outputs()[1] - 0.0).abs() < 1e-6);

        // Test at t=1
        lut.set_input(0, 1.0);
        lut.update(0.0);
        assert!((lut.outputs()[0] - 1.0).abs() < 1e-6);
        assert!((lut.outputs()[1] - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_lut1d_interpolation() {
        // Test interpolation between points
        let points = vec![0.0, 1.0, 2.0];
        let values = vec![vec![0.0, 0.0], vec![1.0, 10.0], vec![4.0, 40.0]];

        let mut lut = LUT1D::new_multi(points, values);

        // Test exact points
        lut.set_input(0, 1.0);
        lut.update(0.0);
        assert!((lut.outputs()[0] - 1.0).abs() < 1e-6);
        assert!((lut.outputs()[1] - 10.0).abs() < 1e-6);

        // Test interpolation at midpoint
        lut.set_input(0, 0.5);
        lut.update(0.0);
        assert!((lut.outputs()[0] - 0.5).abs() < 1e-6); // (0+1)/2
        assert!((lut.outputs()[1] - 5.0).abs() < 1e-6); // (0+10)/2
    }

    #[test]
    fn test_lut1d_extrapolation() {
        // Test extrapolation behavior
        let points = vec![1.0, 2.0, 3.0];
        let values = vec![1.0, 4.0, 9.0];

        let mut lut = LUT1D::new(points, values);

        // Test extrapolation below domain
        lut.set_input(0, 0.0);
        lut.update(0.0);
        assert!(!lut.get_output(0).is_nan());

        // Test extrapolation above domain
        lut.set_input(0, 4.0);
        lut.update(0.0);
        assert!(!lut.get_output(0).is_nan());
    }

    #[test]
    fn test_lut1d_multiple_outputs() {
        // Test initialization with 1D data and multiple outputs
        let points = vec![0.0, 1.0, 2.0];
        let values = vec![vec![0.0, 0.0], vec![1.0, 10.0], vec![4.0, 40.0]];

        let lut = LUT1D::new_multi(points, values);

        // Test that it was properly initialized
        assert_eq!(lut.outputs().len(), 2);
        assert_eq!(lut.inputs.len(), 1);
        assert_eq!(lut.outputs.len(), 2);
    }

    #[test]
    fn test_lut1d_update_siso() {
        // Test update method for SISO case
        let points = vec![0.0, 1.0, 2.0];
        let values = vec![0.0, 10.0, 20.0];

        let mut lut = LUT1D::new(points, values);

        // Set block inputs
        lut.inputs_mut()[0] = 1.0;

        // Update block
        lut.update(0.0);

        // Test if update was correct
        assert_eq!(lut.outputs()[0], 10.0);
    }

    #[test]
    fn test_lut1d_update_simo() {
        // Test update method for SIMO case
        let points = vec![0.0, 1.0, 2.0];
        let values = vec![vec![0.0, 0.0], vec![1.0, 10.0], vec![4.0, 40.0]];

        let mut lut = LUT1D::new_multi(points, values);

        // Set block inputs
        lut.inputs_mut()[0] = 1.0;

        // Update block
        lut.update(0.0);

        // Test if update was correct
        assert!((lut.outputs()[0] - 1.0).abs() < 1e-6);
        assert!((lut.outputs()[1] - 10.0).abs() < 1e-6);
    }

    // ============================================================================
    // LUT Tests (Multi-dimensional)
    // ============================================================================

    #[test]
    fn test_lut_init() {
        // Test initialization of 2D lookup table
        let points = vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 1.0],
        ];
        let values = vec![0.0, 1.0, 1.0, 2.0]; // z = x + y

        let lut = LUT::new(points, values);

        // Verify initialization
        assert_eq!(lut.inputs().len(), 2);
        assert_eq!(lut.outputs().len(), 1);
    }

    #[test]
    fn test_lut_2d_update() {
        // Test update method for 2D MIMO case
        let points = vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 1.0],
        ];
        let values = vec![
            vec![0.0, 0.0],
            vec![1.0, 2.0],
            vec![1.0, 2.0],
            vec![2.0, 4.0],
        ]; // [x+y, 2*(x+y)]

        let mut lut = LUT::new_multi(points, values);

        // Set block inputs
        lut.inputs_mut()[0] = 1.0;
        lut.inputs_mut()[1] = 1.0;

        // Update block
        lut.update(0.0);

        // Test if update was correct (exact match at grid point)
        assert!((lut.outputs()[0] - 2.0).abs() < 1e-6);
        assert!((lut.outputs()[1] - 4.0).abs() < 1e-6);
    }
}
