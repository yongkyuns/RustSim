//! Condition-based event detection
//!
//! Triggers when a boolean condition becomes true, using bisection
//! to locate the event time accurately.

use super::base::{Event, EventDetection};

/// Condition-based event detector
///
/// Monitors a boolean condition and triggers when it transitions from
/// false to true. Uses bisection method to narrow down the event time.
/// Automatically deactivates after the first event resolution.
pub struct Condition<F, A>
where
    F: Fn(f64) -> bool + Send + Sync,
    A: FnMut(f64) + Send + Sync,
{
    func_evt: F,
    func_act: Option<A>,
    tolerance: f64,
    history: Option<(bool, f64)>, // (result, time)
    times: Vec<f64>,
    active: bool,
}

impl<F, A> Condition<F, A>
where
    F: Fn(f64) -> bool + Send + Sync,
    A: FnMut(f64) + Send + Sync,
{
    /// Create a new condition event detector
    ///
    /// # Arguments
    /// * `func_evt` - Condition function that returns true when event should trigger
    /// * `func_act` - Optional action function to execute on event resolution
    /// * `tolerance` - Tolerance for event detection (time interval)
    pub fn new(func_evt: F, func_act: Option<A>, tolerance: f64) -> Self {
        Self {
            func_evt,
            func_act,
            tolerance,
            history: None,
            times: Vec::new(),
            active: true,
        }
    }

    /// Get recorded event times
    pub fn event_times(&self) -> &[f64] {
        &self.times
    }

    /// Get number of detected events
    pub fn len(&self) -> usize {
        self.times.len()
    }

    /// Check if no events have been detected
    pub fn is_empty(&self) -> bool {
        self.times.is_empty()
    }

    /// Activate event detection
    pub fn on(&mut self) {
        self.active = true;
    }

    /// Deactivate event detection
    pub fn off(&mut self) {
        self.active = false;
    }
}

impl<F, A> Event for Condition<F, A>
where
    F: Fn(f64) -> bool + Send + Sync,
    A: FnMut(f64) + Send + Sync,
{
    fn buffer(&mut self, t: f64) {
        let result = (self.func_evt)(t);
        self.history = Some((result, t));
    }

    fn detect(&mut self, t: f64) -> EventDetection {
        let Some((prev_result, prev_t)) = self.history else {
            return EventDetection::default();
        };

        let result = (self.func_evt)(t);

        // Check if interval narrowed down sufficiently
        let close = result && (t - prev_t) < self.tolerance;

        // Close enough to event
        if close {
            return EventDetection {
                detected: true,
                close: true,
                ratio: 1.0,
            };
        }

        // Half the stepsize to creep closer to event (bisection)
        EventDetection {
            detected: result,
            close: false,
            ratio: 0.5,
        }
    }

    fn resolve(&mut self, t: f64) {
        self.times.push(t);
        if let Some(ref mut act) = self.func_act {
            act(t);
        }
        // Deactivate condition tracking after first resolution
        self.off();
    }

    fn is_active(&self) -> bool {
        self.active
    }

    fn reset(&mut self) {
        self.history = None;
        self.times.clear();
        self.active = true;
    }
}
