//! Base event trait and types

/// Result of event detection
#[derive(Debug, Clone, Copy)]
pub struct EventDetection {
    /// Was an event detected?
    pub detected: bool,
    /// Are we close enough to resolve the event?
    pub close: bool,
    /// Interpolation ratio for event time (0.0 to 1.0)
    pub ratio: f64,
}

impl Default for EventDetection {
    fn default() -> Self {
        Self {
            detected: false,
            close: false,
            ratio: 1.0,
        }
    }
}

/// Core event trait for hybrid system events
pub trait Event: Send + Sync {
    /// Buffer state before timestep for interpolation
    fn buffer(&mut self, t: f64);

    /// Detect if an event occurred during the timestep
    fn detect(&mut self, t: f64) -> EventDetection;

    /// Resolve the event (execute action function)
    fn resolve(&mut self, t: f64);

    /// Estimate time until next event (if predictable)
    fn estimate(&self, t: f64) -> Option<f64> {
        let _ = t;
        None
    }

    /// Check if event is active
    fn is_active(&self) -> bool {
        true
    }

    /// Reset event state
    fn reset(&mut self);
}
