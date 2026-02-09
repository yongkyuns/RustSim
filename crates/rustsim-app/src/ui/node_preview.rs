//! Plot preview rendering for Scope and Spectrum nodes

use egui::{Color32, Painter, Pos2, Rect, Stroke};

/// Decimate data using min-max algorithm to limit points for rendering
///
/// Reduces data to approximately `target_points` by picking min/max pairs
/// within each bucket, preserving the visual envelope of the signal.
pub fn decimate_minmax(data: &[(f64, f64)], target_points: usize) -> Vec<(f64, f64)> {
    if data.len() <= target_points {
        return data.to_vec();
    }

    let bucket_size = data.len() / (target_points / 2).max(1);
    let mut result = Vec::with_capacity(target_points);

    for chunk in data.chunks(bucket_size) {
        if chunk.is_empty() {
            continue;
        }

        // Find min and max in this bucket
        let mut min_point = chunk[0];
        let mut max_point = chunk[0];

        for &point in chunk.iter() {
            if point.1 < min_point.1 {
                min_point = point;
            }
            if point.1 > max_point.1 {
                max_point = point;
            }
        }

        // Add points in time order
        if min_point.0 < max_point.0 {
            result.push(min_point);
            result.push(max_point);
        } else {
            result.push(max_point);
            result.push(min_point);
        }
    }

    result
}

/// Draw a small plot preview inside a rectangle
///
/// # Arguments
///
/// * `painter` - egui Painter for drawing
/// * `rect` - Rectangle bounds for the preview
/// * `data` - Time-series data as (x, y) points
/// * `color` - Line color
pub fn draw_plot_preview(painter: &Painter, rect: Rect, data: &[(f64, f64)], color: Color32) {
    if data.is_empty() {
        // Draw placeholder for empty data
        painter.rect_stroke(rect, 2.0, Stroke::new(1.0, Color32::from_gray(60)));
        return;
    }

    // Draw background
    painter.rect_filled(rect, 2.0, Color32::from_rgba_unmultiplied(30, 30, 35, 200));

    // Decimate data to ~100 points for performance
    let decimated = decimate_minmax(data, 100);

    if decimated.len() < 2 {
        return;
    }

    // Find data bounds
    let mut min_x = decimated[0].0;
    let mut max_x = decimated[0].0;
    let mut min_y = decimated[0].1;
    let mut max_y = decimated[0].1;

    for &(x, y) in &decimated {
        min_x = min_x.min(x);
        max_x = max_x.max(x);
        min_y = min_y.min(y);
        max_y = max_y.max(y);
    }

    // Add small margin to y-axis
    let y_range = (max_y - min_y).max(1e-6);
    min_y -= y_range * 0.05;
    let max_y_adjusted = max_y + y_range * 0.05;

    let x_range = (max_x - min_x).max(1e-6);

    // Transform data points to screen coordinates
    let to_screen = |x: f64, y: f64| -> Pos2 {
        let norm_x = ((x - min_x) / x_range) as f32;
        let norm_y = 1.0 - ((y - min_y) / (max_y_adjusted - min_y)) as f32; // Flip y-axis

        let screen_x = rect.min.x + norm_x * rect.width();
        let screen_y = rect.min.y + norm_y * rect.height();

        Pos2::new(screen_x, screen_y)
    };

    // Draw polyline
    let points: Vec<Pos2> = decimated.iter().map(|&(x, y)| to_screen(x, y)).collect();

    for window in points.windows(2) {
        painter.line_segment([window[0], window[1]], Stroke::new(1.0, color));
    }

    // Draw border
    painter.rect_stroke(rect, 2.0, Stroke::new(1.0, Color32::from_gray(80)));
}

/// Draw spectrum (frequency magnitude) preview
///
/// Similar to draw_plot_preview but optimized for frequency-domain data
pub fn draw_spectrum_preview(
    painter: &Painter,
    rect: Rect,
    frequencies: &[f64],
    magnitudes: &[f64],
    color: Color32,
) {
    if frequencies.is_empty() || magnitudes.is_empty() {
        painter.rect_stroke(rect, 2.0, Stroke::new(1.0, Color32::from_gray(60)));
        return;
    }

    // Combine into (freq, mag) pairs
    let data: Vec<(f64, f64)> = frequencies
        .iter()
        .zip(magnitudes.iter())
        .map(|(&f, &m)| (f, m))
        .collect();

    // Draw using standard plot preview
    draw_plot_preview(painter, rect, &data, color);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decimate_empty() {
        let data: Vec<(f64, f64)> = vec![];
        let result = decimate_minmax(&data, 100);
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_decimate_small() {
        let data = vec![(0.0, 1.0), (1.0, 2.0), (2.0, 3.0)];
        let result = decimate_minmax(&data, 100);
        assert_eq!(result.len(), 3);
    }

    #[test]
    fn test_decimate_large() {
        let data: Vec<(f64, f64)> = (0..1000)
            .map(|i| (i as f64, (i as f64 * 0.1).sin()))
            .collect();
        let result = decimate_minmax(&data, 100);
        assert!(result.len() <= 100);
        assert!(result.len() > 0);
    }

    #[test]
    fn test_decimate_preserves_extrema() {
        // Create data with a spike
        let mut data = vec![];
        for i in 0..100 {
            data.push((i as f64, 1.0));
        }
        data.push((50.0, 10.0)); // Spike

        let result = decimate_minmax(&data, 10);

        // Check that spike is preserved
        let max_val = result
            .iter()
            .map(|(_, y)| y)
            .fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        assert!((max_val - 10.0).abs() < 1e-6);
    }
}
