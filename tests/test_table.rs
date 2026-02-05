//! Comprehensive tests for lookup table blocks with complete PathSim test parity
//!
//! This test suite matches all tests from:
//! /home/yongkyunshin/personal/pathsim/tests/pathsim/blocks/test_table.py

use rustsim::prelude::*;

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
    assert_eq!(lut.inputs().len(), 1);
    assert_eq!(lut.outputs().len(), 1);
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

    // Test for several time points
    for t in [0.0, 1.0, 2.0, 3.0] {
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
    let values = vec![
        vec![0.0, 0.0],
        vec![1.0, 2.0],
        vec![4.0, 4.0],
    ]; // [x^2, 2*x]

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
    let values = vec![
        vec![0.0, 0.0],
        vec![1.0, 10.0],
        vec![4.0, 40.0],
    ]; // [x^2 approx, x*10]

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
    assert!(!lut.get_output(0).is_nan()); // Should extrapolate, not NaN

    // Test extrapolation above domain
    lut.set_input(0, 4.0);
    lut.update(0.0);
    assert!(!lut.get_output(0).is_nan()); // Should extrapolate, not NaN
}

#[test]
fn test_lut1d_multiple_outputs() {
    // Test initialization with 1D data and multiple outputs
    let points = vec![0.0, 1.0, 2.0];
    let values = vec![
        vec![0.0, 0.0],
        vec![1.0, 10.0],
        vec![4.0, 40.0],
    ]; // [x^2 approx, x*10]

    let lut = LUT1D::new_multi(points, values);

    // Test that it was properly initialized
    assert_eq!(lut.outputs().len(), 2);
    assert_eq!(lut.inputs().len(), 1);
    assert_eq!(lut.outputs().len(), 2);
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
    let values = vec![
        vec![0.0, 0.0],
        vec![1.0, 10.0],
        vec![4.0, 40.0],
    ]; // [x^2 approx, x*10]

    let mut lut = LUT1D::new_multi(points, values);

    // Set block inputs
    lut.inputs_mut()[0] = 1.0;

    // Update block
    lut.update(0.0);

    // Test if update was correct
    assert!((lut.outputs()[0] - 1.0).abs() < 1e-6);
    assert!((lut.outputs()[1] - 10.0).abs() < 1e-6);
}

#[test]
fn test_lut1d_init_single_output() {
    // Test initialization with 1D data and single output
    let points = vec![0.0, 1.0, 2.0, 3.0];
    let values = vec![0.0, 1.0, 4.0, 9.0]; // y = x^2

    let lut = LUT1D::new(points, values);

    // Test that it was properly initialized
    assert_eq!(lut.outputs().len(), 1);
    assert_eq!(lut.inputs().len(), 1);
}

#[test]
fn test_lut1d_init_multiple_outputs() {
    // Test initialization with 1D data and multiple outputs
    let points = vec![0.0, 1.0, 2.0];
    let values = vec![
        vec![0.0, 0.0],
        vec![1.0, 10.0],
        vec![4.0, 40.0],
    ]; // [x^2 approx, x*10]

    let lut = LUT1D::new_multi(points, values);

    // Test that it was properly initialized
    assert_eq!(lut.outputs().len(), 2);
    assert_eq!(lut.inputs().len(), 1);
}

#[test]
fn test_lut1d_interpolation_multiple_outputs() {
    // Test 1D interpolation functionality with multiple outputs
    let points = vec![0.0, 1.0, 2.0];
    let values = vec![
        vec![0.0, 0.0],
        vec![1.0, 10.0],
        vec![4.0, 40.0],
    ]; // [x^2 approx, x*10]

    let mut lut = LUT1D::new_multi(points, values);

    // Test exact points
    lut.set_input(0, 1.0);
    lut.update(0.0);
    assert!((lut.outputs()[0] - 1.0).abs() < 1e-6);
    assert!((lut.outputs()[1] - 10.0).abs() < 1e-6);

    // Test interpolation
    lut.set_input(0, 0.5); // Halfway between 0 and 1
    lut.update(0.0);
    assert!((lut.outputs()[0] - 0.5).abs() < 1e-6); // (0+1)/2
    assert!((lut.outputs()[1] - 5.0).abs() < 1e-6); // (0+10)/2
}

#[test]
fn test_lut1d_update_functionality_multiple_outputs() {
    // Test the update method directly with multiple outputs
    let points = vec![0.0, 1.0, 2.0];
    let values = vec![
        vec![0.0, 0.0],
        vec![1.0, 10.0],
        vec![4.0, 40.0],
    ]; // [x^2 approx, x*10]

    let mut lut = LUT1D::new_multi(points, values);

    // Test update with exact points
    lut.set_input(0, 1.0);
    lut.update(0.0);
    assert!((lut.outputs()[0] - 1.0).abs() < 1e-6);
    assert!((lut.outputs()[1] - 10.0).abs() < 1e-6);

    // Test interpolation
    lut.set_input(0, 0.5);
    lut.update(0.0);
    assert!((lut.outputs()[0] - 0.5).abs() < 1e-6);
    assert!((lut.outputs()[1] - 5.0).abs() < 1e-6);
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

    // Set block inputs to exact grid point
    lut.inputs_mut()[0] = 1.0;
    lut.inputs_mut()[1] = 1.0;

    // Update block
    lut.update(0.0);

    // Test if update was correct (exact match at grid point)
    assert!((lut.outputs()[0] - 2.0).abs() < 1e-6);
    assert!((lut.outputs()[1] - 4.0).abs() < 1e-6);
}

// ============================================================================
// Additional edge case tests
// ============================================================================

#[test]
fn test_lut1d_single_point_extrapolation() {
    // Edge case: single point should return constant value
    let points = vec![1.0];
    let values = vec![5.0];

    let mut lut = LUT1D::new(points, values);

    // Below the point
    lut.set_input(0, 0.0);
    lut.update(0.0);
    assert_eq!(lut.get_output(0), 5.0);

    // At the point
    lut.set_input(0, 1.0);
    lut.update(0.0);
    assert_eq!(lut.get_output(0), 5.0);

    // Above the point
    lut.set_input(0, 2.0);
    lut.update(0.0);
    assert_eq!(lut.get_output(0), 5.0);
}

#[test]
fn test_lut1d_reset() {
    // Test reset functionality
    let points = vec![0.0, 1.0, 2.0];
    let values = vec![0.0, 10.0, 20.0];

    let mut lut = LUT1D::new(points, values);

    // Set inputs and update
    lut.set_input(0, 1.0);
    lut.update(0.0);
    assert_eq!(lut.get_output(0), 10.0);

    // Reset
    lut.reset();
    assert_eq!(lut.inputs()[0], 0.0);
    assert_eq!(lut.outputs()[0], 0.0);
}

#[test]
fn test_lut_reset() {
    // Test reset functionality for multi-dimensional LUT
    let points = vec![
        vec![0.0, 0.0],
        vec![1.0, 0.0],
        vec![0.0, 1.0],
        vec![1.0, 1.0],
    ];
    let values = vec![0.0, 1.0, 1.0, 2.0];

    let mut lut = LUT::new(points, values);

    // Set inputs and update
    lut.inputs_mut()[0] = 1.0;
    lut.inputs_mut()[1] = 1.0;
    lut.update(0.0);

    // Reset
    lut.reset();
    assert_eq!(lut.inputs()[0], 0.0);
    assert_eq!(lut.inputs()[1], 0.0);
    assert_eq!(lut.outputs()[0], 0.0);
}

#[test]
fn test_lut1d_linear_extrapolation_slope() {
    // Verify that extrapolation follows linear trend
    let points = vec![0.0, 1.0, 2.0];
    let values = vec![0.0, 2.0, 4.0]; // y = 2*x

    let mut lut = LUT1D::new(points, values);

    // Extrapolate below
    lut.set_input(0, -1.0);
    lut.update(0.0);
    assert!((lut.get_output(0) - (-2.0)).abs() < 1e-6);

    // Extrapolate above
    lut.set_input(0, 3.0);
    lut.update(0.0);
    assert!((lut.get_output(0) - 6.0).abs() < 1e-6);
}

#[test]
fn test_lut1d_multiple_outputs_extrapolation() {
    // Test that extrapolation works correctly with multiple outputs
    let points = vec![0.0, 1.0];
    let values = vec![
        vec![0.0, 0.0],
        vec![2.0, 4.0],
    ]; // [2*x, 4*x]

    let mut lut = LUT1D::new_multi(points, values);

    // Extrapolate above
    lut.set_input(0, 2.0);
    lut.update(0.0);
    assert!((lut.outputs()[0] - 4.0).abs() < 1e-6);
    assert!((lut.outputs()[1] - 8.0).abs() < 1e-6);
}

#[test]
fn test_lut_interpolation_interior() {
    // Test interpolation at interior points
    let points = vec![
        vec![0.0, 0.0],
        vec![2.0, 0.0],
        vec![0.0, 2.0],
        vec![2.0, 2.0],
    ];
    let values = vec![0.0, 2.0, 2.0, 4.0]; // z = x + y

    let mut lut = LUT::new(points, values);

    // Test at center point (should interpolate to 2.0)
    lut.inputs_mut()[0] = 1.0;
    lut.inputs_mut()[1] = 1.0;
    lut.update(0.0);

    // With inverse distance weighting, center should be close to 2.0
    assert!((lut.outputs()[0] - 2.0).abs() < 0.5);
}
