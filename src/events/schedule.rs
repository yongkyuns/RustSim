//! Scheduled event detection
//!
//! Provides periodic and list-based scheduled events for time-dependent
//! event triggering.

use super::base::{Event, EventDetection};
use crate::utils::constants::TOLERANCE;

/// Periodic scheduled event detector
///
/// Triggers at regular intervals starting from a specified time.
/// Optionally supports an end time after which events stop.
pub struct Schedule<A>
where
    A: FnMut(f64) + Send + Sync,
{
    t_start: f64,
    t_end: Option<f64>,
    t_period: f64,
    func_act: Option<A>,
    tolerance: f64,
    history: Option<((), f64)>, // (_, time)
    times: Vec<f64>,
    active: bool,
}

impl<A> Schedule<A>
where
    A: FnMut(f64) + Send + Sync,
{
    /// Create a new periodic scheduled event detector
    ///
    /// # Arguments
    /// * `t_start` - Starting time for schedule
    /// * `t_end` - Optional termination time for schedule
    /// * `t_period` - Time period between events
    /// * `func_act` - Optional action function to execute on event resolution
    /// * `tolerance` - Tolerance for event detection
    pub fn new(
        t_start: f64,
        t_end: Option<f64>,
        t_period: f64,
        func_act: Option<A>,
        tolerance: f64,
    ) -> Self {
        Self {
            t_start,
            t_end,
            t_period,
            func_act,
            tolerance,
            history: None,
            times: Vec::new(),
            active: true,
        }
    }

    /// Get the next scheduled event time
    pub fn next(&self) -> f64 {
        self.t_start + (self.times.len() as f64) * self.t_period
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

impl<A> Event for Schedule<A>
where
    A: FnMut(f64) + Send + Sync,
{
    fn buffer(&mut self, t: f64) {
        self.history = Some(((), t));
    }

    fn detect(&mut self, t: f64) -> EventDetection {
        let t_next = self.next();

        // End time reached? -> deactivate event
        if let Some(t_end) = self.t_end {
            if t_next > t_end {
                self.active = false;
                return EventDetection::default();
            }
        }

        // No event yet
        if t_next > t {
            return EventDetection::default();
        }

        // Are we close enough to the scheduled event?
        if (t_next - t).abs() <= self.tolerance {
            return EventDetection {
                detected: true,
                close: true,
                ratio: 0.0,
            };
        }

        let Some((_, prev_t)) = self.history else {
            return EventDetection {
                detected: true,
                close: true,
                ratio: 0.0,
            };
        };

        // Have we already passed the event? (first timestep)
        if prev_t >= t_next {
            return EventDetection {
                detected: true,
                close: true,
                ratio: 0.0,
            };
        }

        // Calculate timestep ratio
        let ratio = (t_next - prev_t) / (t - prev_t).max(TOLERANCE);

        EventDetection {
            detected: true,
            close: false,
            ratio,
        }
    }

    fn resolve(&mut self, t: f64) {
        self.times.push(t);
        if let Some(ref mut act) = self.func_act {
            act(t);
        }
    }

    fn estimate(&self, t: f64) -> Option<f64> {
        Some(self.next() - t)
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

/// List-based scheduled event detector
///
/// Triggers at specific times from a provided list.
/// Automatically deactivates after all events have been processed.
pub struct ScheduleList<A>
where
    A: FnMut(f64) + Send + Sync,
{
    pub times_evt: Vec<f64>,
    pub func_act: Option<A>,
    tolerance: f64,
    history: Option<((), f64)>, // (_, time)
    times: Vec<f64>,
    active: bool,
}

impl<A> ScheduleList<A>
where
    A: FnMut(f64) + Send + Sync,
{
    /// Create a new list-based scheduled event detector
    ///
    /// # Arguments
    /// * `times_evt` - List of event times in ascending order
    /// * `func_act` - Optional action function to execute on event resolution
    /// * `tolerance` - Tolerance for event detection
    ///
    /// # Errors
    /// Returns an error if `times_evt` is not in ascending order
    pub fn new(times_evt: Vec<f64>, func_act: Option<A>, tolerance: f64) -> Result<Self, String> {
        // Validate times are in ascending order
        if times_evt.len() > 1 {
            for i in 1..times_evt.len() {
                if times_evt[i] <= times_evt[i - 1] {
                    return Err("'times_evt' need to be in ascending order!".to_string());
                }
            }
        }

        Ok(Self {
            times_evt,
            func_act,
            tolerance,
            history: None,
            times: Vec::new(),
            active: true,
        })
    }

    /// Get the next scheduled event time from the list
    pub fn next(&self) -> f64 {
        let n = self.times.len();
        if n < self.times_evt.len() {
            self.times_evt[n]
        } else {
            self.times_evt[self.times_evt.len() - 1]
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

impl<A> Event for ScheduleList<A>
where
    A: FnMut(f64) + Send + Sync,
{
    fn buffer(&mut self, t: f64) {
        self.history = Some(((), t));
    }

    fn detect(&mut self, t: f64) -> EventDetection {
        // Check if out of bounds
        let n = self.times.len();
        if n >= self.times_evt.len() {
            self.active = false;
            return EventDetection::default();
        }

        let t_next = self.next();

        // No event yet
        if t_next > t {
            return EventDetection::default();
        }

        // Are we close enough to the scheduled event?
        if (t_next - t).abs() <= self.tolerance {
            return EventDetection {
                detected: true,
                close: true,
                ratio: 0.0,
            };
        }

        let Some((_, prev_t)) = self.history else {
            return EventDetection {
                detected: true,
                close: true,
                ratio: 0.0,
            };
        };

        // Have we already passed the event? (first timestep)
        if prev_t >= t_next {
            return EventDetection {
                detected: true,
                close: true,
                ratio: 0.0,
            };
        }

        // Calculate timestep ratio
        let ratio = (t_next - prev_t) / (t - prev_t).max(TOLERANCE);

        EventDetection {
            detected: true,
            close: false,
            ratio,
        }
    }

    fn resolve(&mut self, t: f64) {
        self.times.push(t);
        if let Some(ref mut act) = self.func_act {
            act(t);
        }
    }

    fn estimate(&self, t: f64) -> Option<f64> {
        Some(self.next() - t)
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
