//! Zero-crossing event detectors for hybrid systems
//!
//! Provides bidirectional and unidirectional zero-crossing detection
//! with secant-based interpolation for accurate event localization.

use super::base::{Event, EventDetection};
use crate::utils::constants::TOLERANCE;

/// Bidirectional zero-crossing event detector
///
/// Triggers when the event function crosses zero in either direction
/// (positive to negative or negative to positive).
pub struct ZeroCrossing<F, A>
where
    F: Fn(f64) -> f64 + Send + Sync,
    A: FnMut(f64) + Send + Sync,
{
    func_evt: F,
    func_act: Option<A>,
    tolerance: f64,
    history: Option<(f64, f64)>, // (result, time)
    times: Vec<f64>,
    active: bool,
}

impl<F, A> ZeroCrossing<F, A>
where
    F: Fn(f64) -> f64 + Send + Sync,
    A: FnMut(f64) + Send + Sync,
{
    /// Create a new bidirectional zero-crossing event detector
    ///
    /// # Arguments
    /// * `func_evt` - Event function where zeros represent events
    /// * `func_act` - Optional action function to execute on event resolution
    /// * `tolerance` - Tolerance for event detection (default: EVT_TOLERANCE = 1e-4)
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
}

impl<F, A> Event for ZeroCrossing<F, A>
where
    F: Fn(f64) -> f64 + Send + Sync,
    A: FnMut(f64) + Send + Sync,
{
    fn buffer(&mut self, t: f64) {
        let result = (self.func_evt)(t);
        self.history = Some((result, t));
    }

    fn detect(&mut self, t: f64) -> EventDetection {
        let Some((prev_result, _prev_t)) = self.history else {
            return EventDetection::default();
        };

        let result = (self.func_evt)(t);

        // Exactly hit zero
        if result == 0.0 && prev_result != 0.0 {
            return EventDetection {
                detected: true,
                close: true,
                ratio: 1.0,
            };
        }

        // Check if we're close to the event
        let close = result.abs() <= self.tolerance;

        // Check for sign change (bidirectional)
        let is_event = (result * prev_result).signum() < 0.0;

        if !is_event {
            return EventDetection {
                detected: false,
                close: false,
                ratio: 1.0,
            };
        }

        // Secant method interpolation
        let ratio = prev_result.abs() / (prev_result - result).abs().max(TOLERANCE);

        EventDetection {
            detected: true,
            close,
            ratio,
        }
    }

    fn resolve(&mut self, t: f64) {
        self.times.push(t);
        if let Some(ref mut act) = self.func_act {
            act(t);
        }
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

/// Upward zero-crossing event detector
///
/// Triggers only when the event function crosses zero from negative to positive.
pub struct ZeroCrossingUp<F, A>
where
    F: Fn(f64) -> f64 + Send + Sync,
    A: FnMut(f64) + Send + Sync,
{
    func_evt: F,
    func_act: Option<A>,
    tolerance: f64,
    history: Option<(f64, f64)>, // (result, time)
    times: Vec<f64>,
    active: bool,
}

impl<F, A> ZeroCrossingUp<F, A>
where
    F: Fn(f64) -> f64 + Send + Sync,
    A: FnMut(f64) + Send + Sync,
{
    /// Create a new upward zero-crossing event detector
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

    pub fn event_times(&self) -> &[f64] {
        &self.times
    }

    pub fn len(&self) -> usize {
        self.times.len()
    }

    pub fn is_empty(&self) -> bool {
        self.times.is_empty()
    }
}

impl<F, A> Event for ZeroCrossingUp<F, A>
where
    F: Fn(f64) -> f64 + Send + Sync,
    A: FnMut(f64) + Send + Sync,
{
    fn buffer(&mut self, t: f64) {
        let result = (self.func_evt)(t);
        self.history = Some((result, t));
    }

    fn detect(&mut self, t: f64) -> EventDetection {
        let result = (self.func_evt)(t);

        let Some((prev_result, _prev_t)) = self.history else {
            return EventDetection::default();
        };

        // Exactly hit zero (from negative side)
        if result == 0.0 && prev_result < 0.0 {
            return EventDetection {
                detected: true,
                close: true,
                ratio: 1.0,
            };
        }

        // Check if we're close to the event
        let close = result.abs() <= self.tolerance;

        // Check for sign change AND upward direction
        let is_event = (result * prev_result).signum() < 0.0 && result > prev_result;

        // No event or wrong direction
        if !is_event || prev_result >= 0.0 {
            return EventDetection {
                detected: false,
                close: false,
                ratio: 1.0,
            };
        }

        // Secant method interpolation
        let ratio = prev_result.abs() / (prev_result - result).abs().max(TOLERANCE);

        EventDetection {
            detected: true,
            close,
            ratio,
        }
    }

    fn resolve(&mut self, t: f64) {
        self.times.push(t);
        if let Some(ref mut act) = self.func_act {
            act(t);
        }
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

/// Downward zero-crossing event detector
///
/// Triggers only when the event function crosses zero from positive to negative.
pub struct ZeroCrossingDown<F, A>
where
    F: Fn(f64) -> f64 + Send + Sync,
    A: FnMut(f64) + Send + Sync,
{
    func_evt: F,
    func_act: Option<A>,
    tolerance: f64,
    history: Option<(f64, f64)>, // (result, time)
    times: Vec<f64>,
    active: bool,
}

impl<F, A> ZeroCrossingDown<F, A>
where
    F: Fn(f64) -> f64 + Send + Sync,
    A: FnMut(f64) + Send + Sync,
{
    /// Create a new downward zero-crossing event detector
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

    pub fn event_times(&self) -> &[f64] {
        &self.times
    }

    pub fn len(&self) -> usize {
        self.times.len()
    }

    pub fn is_empty(&self) -> bool {
        self.times.is_empty()
    }
}

impl<F, A> Event for ZeroCrossingDown<F, A>
where
    F: Fn(f64) -> f64 + Send + Sync,
    A: FnMut(f64) + Send + Sync,
{
    fn buffer(&mut self, t: f64) {
        let result = (self.func_evt)(t);
        self.history = Some((result, t));
    }

    fn detect(&mut self, t: f64) -> EventDetection {
        let result = (self.func_evt)(t);

        let Some((prev_result, _prev_t)) = self.history else {
            return EventDetection::default();
        };

        // Exactly hit zero (from positive side)
        if result == 0.0 && prev_result > 0.0 {
            return EventDetection {
                detected: true,
                close: true,
                ratio: 1.0,
            };
        }

        // Check if we're close to the event
        let close = result.abs() <= self.tolerance;

        // Check for sign change AND downward direction
        let is_event = (result * prev_result).signum() < 0.0 && result < prev_result;

        // No event or wrong direction
        if !is_event || prev_result <= 0.0 {
            return EventDetection {
                detected: false,
                close: false,
                ratio: 1.0,
            };
        }

        // Secant method interpolation
        let ratio = prev_result.abs() / (prev_result - result).abs().max(TOLERANCE);

        EventDetection {
            detected: true,
            close,
            ratio,
        }
    }

    fn resolve(&mut self, t: f64) {
        self.times.push(t);
        if let Some(ref mut act) = self.func_act {
            act(t);
        }
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
