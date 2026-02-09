//! CSV export functionality for plot data.

use std::path::Path;

/// Represents plot data that can be exported to CSV
#[derive(Debug, Clone)]
pub enum PlotData {
    /// Time-domain data: (time, [signal values])
    TimeDomain {
        time: Vec<f64>,
        signals: Vec<Vec<f64>>,
        labels: Vec<String>,
    },
    /// Frequency-domain data: (frequency, [magnitude values])
    Spectrum {
        frequency: Vec<f64>,
        magnitudes: Vec<Vec<f64>>,
        labels: Vec<String>,
    },
}

/// Exports plot data to CSV format.
///
/// # Arguments
///
/// * `data` - The plot data to export
///
/// # Returns
///
/// A CSV-formatted string
///
/// # Example
///
/// ```
/// use rustsim_app::plotting::{PlotData, export_csv};
///
/// let data = PlotData::TimeDomain {
///     time: vec![0.0, 0.1, 0.2],
///     signals: vec![
///         vec![1.0, 1.1, 1.2],
///         vec![2.0, 2.1, 2.2],
///     ],
///     labels: vec!["Signal 0".to_string(), "Signal 1".to_string()],
/// };
///
/// let csv = export_csv(&data);
/// assert!(csv.contains("time,Signal 0,Signal 1"));
/// ```
pub fn export_csv(data: &PlotData) -> String {
    match data {
        PlotData::TimeDomain {
            time,
            signals,
            labels,
        } => export_time_domain_csv(time, signals, labels),
        PlotData::Spectrum {
            frequency,
            magnitudes,
            labels,
        } => export_spectrum_csv(frequency, magnitudes, labels),
    }
}

/// Export time-domain data to CSV
fn export_time_domain_csv(time: &[f64], signals: &[Vec<f64>], labels: &[String]) -> String {
    let mut csv = String::new();

    // Header row
    csv.push_str("time");
    for label in labels {
        csv.push(',');
        csv.push_str(&escape_csv_field(label));
    }
    csv.push('\n');

    // Data rows
    for (i, &t) in time.iter().enumerate() {
        csv.push_str(&format_number(t));

        for signal in signals {
            csv.push(',');
            if i < signal.len() {
                csv.push_str(&format_number(signal[i]));
            }
        }
        csv.push('\n');
    }

    csv
}

/// Export spectrum data to CSV
fn export_spectrum_csv(frequency: &[f64], magnitudes: &[Vec<f64>], labels: &[String]) -> String {
    let mut csv = String::new();

    // Header row
    csv.push_str("frequency");
    for label in labels {
        csv.push(',');
        csv.push_str(&escape_csv_field(label));
    }
    csv.push('\n');

    // Data rows
    for (i, &f) in frequency.iter().enumerate() {
        csv.push_str(&format_number(f));

        for magnitude in magnitudes {
            csv.push(',');
            if i < magnitude.len() {
                csv.push_str(&format_number(magnitude[i]));
            }
        }
        csv.push('\n');
    }

    csv
}

/// Escape a CSV field if it contains special characters
fn escape_csv_field(field: &str) -> String {
    if field.contains(',') || field.contains('"') || field.contains('\n') {
        // Quote the field and escape internal quotes
        format!("\"{}\"", field.replace('"', "\"\""))
    } else {
        field.to_string()
    }
}

/// Format a number for CSV output
fn format_number(num: f64) -> String {
    // Use scientific notation for very large or very small numbers
    if num.abs() < 1e-10 && num != 0.0 {
        format!("{:e}", num)
    } else if num.abs() > 1e10 {
        format!("{:e}", num)
    } else {
        num.to_string()
    }
}

/// Sanitize a filename by replacing invalid characters
pub fn sanitize_filename(name: &str) -> String {
    name.chars()
        .map(|c| match c {
            '/' | '\\' | ':' | '*' | '?' | '"' | '<' | '>' | '|' => '_',
            c => c,
        })
        .collect()
}

/// Save CSV data to a file (native only)
#[cfg(not(target_arch = "wasm32"))]
pub fn save_csv_to_file(csv_data: &str, path: &Path) -> std::io::Result<()> {
    std::fs::write(path, csv_data)
}

/// Trigger browser download of CSV data (WASM only)
#[cfg(target_arch = "wasm32")]
pub fn download_csv_browser(csv_data: &str, filename: &str) {
    use wasm_bindgen::JsCast;

    if let Some(window) = web_sys::window() {
        if let Some(document) = window.document() {
            // Create a blob
            let array = js_sys::Array::new();
            array.push(&wasm_bindgen::JsValue::from_str(csv_data));

            let mut blob_options = web_sys::BlobPropertyBag::new();
            blob_options.type_("text/csv;charset=utf-8;");

            if let Ok(blob) =
                web_sys::Blob::new_with_str_sequence_and_options(&array, &blob_options)
            {
                // Create object URL
                if let Ok(url) = web_sys::Url::create_object_url_with_blob(&blob) {
                    // Create download link
                    if let Ok(Some(element)) = document.create_element("a") {
                        if let Ok(link) = element.dyn_into::<web_sys::HtmlAnchorElement>() {
                            link.set_href(&url);
                            link.set_download(filename);
                            link.style().set_property("display", "none").ok();

                            // Append, click, and remove
                            if let Some(body) = document.body() {
                                body.append_child(&link).ok();
                                link.click();
                                body.remove_child(&link).ok();
                            }

                            // Revoke object URL
                            web_sys::Url::revoke_object_url(&url).ok();
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_csv_format_time_domain() {
        let data = PlotData::TimeDomain {
            time: vec![0.0, 0.1, 0.2],
            signals: vec![vec![1.0, 1.5, 2.0], vec![2.0, 2.5, 3.0]],
            labels: vec!["Signal 0".to_string(), "Signal 1".to_string()],
        };

        let csv = export_csv(&data);
        let lines: Vec<&str> = csv.lines().collect();

        assert_eq!(lines.len(), 4); // Header + 3 data rows
        assert_eq!(lines[0], "time,Signal 0,Signal 1");
        assert_eq!(lines[1], "0,1,2");
        assert_eq!(lines[2], "0.1,1.5,2.5");
        assert_eq!(lines[3], "0.2,2,3");
    }

    #[test]
    fn test_csv_format_spectrum() {
        let data = PlotData::Spectrum {
            frequency: vec![0.0, 10.0, 20.0],
            magnitudes: vec![vec![0.5, 0.6, 0.7], vec![0.3, 0.4, 0.5]],
            labels: vec!["Magnitude 0".to_string(), "Magnitude 1".to_string()],
        };

        let csv = export_csv(&data);
        let lines: Vec<&str> = csv.lines().collect();

        assert_eq!(lines.len(), 4); // Header + 3 data rows
        assert_eq!(lines[0], "frequency,Magnitude 0,Magnitude 1");
        assert_eq!(lines[1], "0,0.5,0.3");
        assert_eq!(lines[2], "10,0.6,0.4");
        assert_eq!(lines[3], "20,0.7,0.5");
    }

    #[test]
    fn test_csv_escaping_commas() {
        let data = PlotData::TimeDomain {
            time: vec![0.0],
            signals: vec![vec![1.0]],
            labels: vec!["Signal, with comma".to_string()],
        };

        let csv = export_csv(&data);
        assert!(csv.contains("\"Signal, with comma\""));
    }

    #[test]
    fn test_csv_escaping_quotes() {
        let data = PlotData::TimeDomain {
            time: vec![0.0],
            signals: vec![vec![1.0]],
            labels: vec!["Signal \"quoted\"".to_string()],
        };

        let csv = export_csv(&data);
        assert!(csv.contains("\"Signal \"\"quoted\"\"\""));
    }

    #[test]
    fn test_csv_escaping_newlines() {
        let data = PlotData::TimeDomain {
            time: vec![0.0],
            signals: vec![vec![1.0]],
            labels: vec!["Signal\nwith newline".to_string()],
        };

        let csv = export_csv(&data);
        assert!(csv.contains("\"Signal\nwith newline\""));
    }

    #[test]
    fn test_sanitize_filename() {
        assert_eq!(sanitize_filename("normal_name.csv"), "normal_name.csv");
        assert_eq!(sanitize_filename("name/with\\slashes"), "name_with_slashes");
        assert_eq!(
            sanitize_filename("name:with*special?chars"),
            "name_with_special_chars"
        );
        assert_eq!(sanitize_filename("name\"<>|chars"), "name____chars");
    }

    #[test]
    fn test_format_number_normal() {
        assert_eq!(format_number(1.234), "1.234");
        assert_eq!(format_number(0.0), "0");
        assert_eq!(format_number(-5.678), "-5.678");
    }

    #[test]
    fn test_format_number_scientific() {
        let large = format_number(1e15);
        assert!(large.contains('e') || large.contains('E'));

        let small = format_number(1e-15);
        assert!(small.contains('e') || small.contains('E'));
    }

    #[test]
    fn test_empty_signals() {
        let data = PlotData::TimeDomain {
            time: vec![0.0, 0.1, 0.2],
            signals: vec![],
            labels: vec![],
        };

        let csv = export_csv(&data);
        let lines: Vec<&str> = csv.lines().collect();

        assert_eq!(lines.len(), 4); // Header + 3 data rows
        assert_eq!(lines[0], "time");
    }

    #[test]
    fn test_mismatched_lengths() {
        // Signal shorter than time
        let data = PlotData::TimeDomain {
            time: vec![0.0, 0.1, 0.2, 0.3],
            signals: vec![vec![1.0, 2.0]],
            labels: vec!["Signal 0".to_string()],
        };

        let csv = export_csv(&data);
        let lines: Vec<&str> = csv.lines().collect();

        assert_eq!(lines.len(), 5); // Header + 4 data rows
                                    // Should handle gracefully with empty values for missing data
        assert!(lines[3].ends_with(',') || lines[3].split(',').count() == 2);
    }
}
