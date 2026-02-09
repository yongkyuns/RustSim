//! Min-max decimation algorithm for efficient plot rendering.

/// Decimates data points while preserving visual peaks and valleys.
///
/// This algorithm divides the data into buckets and outputs the minimum and maximum
/// y values for each bucket, preserving the visual appearance of the signal while
/// reducing the number of points that need to be rendered.
///
/// # Arguments
///
/// * `x` - X-axis data (typically time)
/// * `y` - Y-axis data (signal values)
/// * `target_points` - Target number of buckets (output will be ~2x this)
///
/// # Returns
///
/// A tuple of (decimated_x, decimated_y) vectors with approximately `target_points * 2` points.
/// Points are returned in chronological order.
///
/// # Example
///
/// ```
/// use rustsim_app::plotting::decimate_minmax;
///
/// let x = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
/// let y = vec![0.0, 1.0, -1.0, 2.0, -2.0, 0.0];
/// let (dec_x, dec_y) = decimate_minmax(&x, &y, 2);
///
/// // Should preserve peaks and valleys
/// assert_eq!(dec_x.len(), dec_y.len());
/// ```
pub fn decimate_minmax(x: &[f64], y: &[f64], target_points: usize) -> (Vec<f64>, Vec<f64>) {
    // Handle edge cases
    if x.is_empty() || y.is_empty() {
        return (Vec::new(), Vec::new());
    }

    if x.len() != y.len() {
        // Mismatched lengths - return empty
        return (Vec::new(), Vec::new());
    }

    if x.len() <= target_points * 2 {
        // Already small enough, return as-is
        return (x.to_vec(), y.to_vec());
    }

    if target_points == 0 {
        return (Vec::new(), Vec::new());
    }

    // Calculate bucket size
    let bucket_size = x.len() / target_points;
    if bucket_size == 0 {
        return (x.to_vec(), y.to_vec());
    }

    let mut result_x = Vec::with_capacity(target_points * 2);
    let mut result_y = Vec::with_capacity(target_points * 2);

    // Process each bucket
    for bucket_idx in 0..target_points {
        let start = bucket_idx * bucket_size;
        let end = if bucket_idx == target_points - 1 {
            // Last bucket includes any remainder
            x.len()
        } else {
            (bucket_idx + 1) * bucket_size
        };

        if start >= x.len() {
            break;
        }

        // Find min and max in this bucket
        let mut min_idx = start;
        let mut max_idx = start;
        let mut min_val = y[start];
        let mut max_val = y[start];

        for i in start..end.min(x.len()) {
            if y[i] < min_val {
                min_val = y[i];
                min_idx = i;
            }
            if y[i] > max_val {
                max_val = y[i];
                max_idx = i;
            }
        }

        // Add points in chronological order
        if min_idx < max_idx {
            // Min comes first
            result_x.push(x[min_idx]);
            result_y.push(y[min_idx]);
            result_x.push(x[max_idx]);
            result_y.push(y[max_idx]);
        } else if max_idx < min_idx {
            // Max comes first
            result_x.push(x[max_idx]);
            result_y.push(y[max_idx]);
            result_x.push(x[min_idx]);
            result_y.push(y[min_idx]);
        } else {
            // Same point is both min and max (flat region)
            result_x.push(x[min_idx]);
            result_y.push(y[min_idx]);
        }
    }

    (result_x, result_y)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decimate_empty() {
        let x: Vec<f64> = vec![];
        let y: Vec<f64> = vec![];
        let (dec_x, dec_y) = decimate_minmax(&x, &y, 10);
        assert!(dec_x.is_empty());
        assert!(dec_y.is_empty());
    }

    #[test]
    fn test_decimate_mismatched_lengths() {
        let x = vec![0.0, 1.0, 2.0];
        let y = vec![0.0, 1.0];
        let (dec_x, dec_y) = decimate_minmax(&x, &y, 10);
        assert!(dec_x.is_empty());
        assert!(dec_y.is_empty());
    }

    #[test]
    fn test_decimate_already_small() {
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let y = vec![0.0, 1.0, -1.0, 2.0];
        let (dec_x, dec_y) = decimate_minmax(&x, &y, 10);
        assert_eq!(dec_x, x);
        assert_eq!(dec_y, y);
    }

    #[test]
    fn test_decimate_preserves_peaks() {
        // Create a signal with clear peaks and valleys
        let x: Vec<f64> = (0..100).map(|i| i as f64 * 0.1).collect();
        let y: Vec<f64> = x
            .iter()
            .map(|&t| (t * 2.0 * std::f64::consts::PI).sin())
            .collect();

        let (dec_x, dec_y) = decimate_minmax(&x, &y, 10);

        // Should have approximately 20 points (10 buckets * 2)
        assert!(dec_x.len() >= 10);
        assert!(dec_x.len() <= 25); // Some wiggle room
        assert_eq!(dec_x.len(), dec_y.len());

        // Find max and min in original
        let orig_max = y.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let orig_min = y.iter().copied().fold(f64::INFINITY, f64::min);

        // Find max and min in decimated
        let dec_max = dec_y.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let dec_min = dec_y.iter().copied().fold(f64::INFINITY, f64::min);

        // Decimated data should preserve the peaks and valleys
        assert!((dec_max - orig_max).abs() < 0.1, "Max should be preserved");
        assert!((dec_min - orig_min).abs() < 0.1, "Min should be preserved");
    }

    #[test]
    fn test_decimate_reduces_points() {
        let n = 10000;
        let x: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let y: Vec<f64> = (0..n).map(|i| (i as f64 * 0.01).sin()).collect();

        let target = 100;
        let (dec_x, dec_y) = decimate_minmax(&x, &y, target);

        // Should significantly reduce number of points
        assert!(dec_x.len() < n / 10);
        assert!(dec_x.len() >= target); // At least target buckets
        assert!(dec_x.len() <= target * 2 + 10); // At most 2 per bucket plus some wiggle room
        assert_eq!(dec_x.len(), dec_y.len());
    }

    #[test]
    fn test_decimate_chronological() {
        let x: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let y: Vec<f64> = (0..100).map(|i| (i % 10) as f64).collect();

        let (dec_x, _) = decimate_minmax(&x, &y, 10);

        // X values should be monotonically increasing (chronological order)
        for i in 1..dec_x.len() {
            assert!(
                dec_x[i] >= dec_x[i - 1],
                "Points should be in chronological order: {} >= {}",
                dec_x[i],
                dec_x[i - 1]
            );
        }
    }

    #[test]
    fn test_decimate_zero_target() {
        let x = vec![0.0, 1.0, 2.0];
        let y = vec![0.0, 1.0, -1.0];
        let (dec_x, dec_y) = decimate_minmax(&x, &y, 0);
        assert!(dec_x.is_empty());
        assert!(dec_y.is_empty());
    }

    #[test]
    fn test_decimate_flat_signal() {
        // All values the same
        let x: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let y: Vec<f64> = vec![5.0; 100];

        let (dec_x, dec_y) = decimate_minmax(&x, &y, 10);

        // Should still produce output
        assert!(!dec_x.is_empty());
        assert_eq!(dec_x.len(), dec_y.len());

        // All y values should be 5.0
        for &val in &dec_y {
            assert_eq!(val, 5.0);
        }
    }
}
