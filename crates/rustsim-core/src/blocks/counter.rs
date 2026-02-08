//! Counter blocks for event counting

use crate::block::{AlgebraicBlock, Block};

/// Bidirectional counter block
///
/// Counts the number of times the input crosses a threshold (in either direction).
/// The output equals start + number_of_crossings.
///
/// # Parameters
///
/// - `start`: Initial counter value
/// - `threshold`: Threshold for zero-crossing detection
///
/// # Behavior
///
/// Detects when the input crosses the threshold in either direction (up or down)
/// and increments or decrements the counter accordingly.
///
/// # Example
///
/// ```ignore
/// let mut counter = Counter::new(0, 5.0);
/// counter.set_input(0, 0.0);
/// counter.update(0.0);
/// assert_eq!(counter.get_output(0), 0); // start value
///
/// counter.set_input(0, 10.0); // Crosses threshold upward
/// counter.update(0.1);
/// assert_eq!(counter.get_output(0), 1); // start + 1
///
/// counter.set_input(0, 0.0); // Crosses threshold downward
/// counter.update(0.2);
/// assert_eq!(counter.get_output(0), 0); // start + 1 - 1
/// ```
#[derive(Debug, Clone)]
pub struct Counter {
    input: f64,
    output: f64,
    start: i32,
    threshold: f64,
    count: i32,
    last_state: Option<bool>, // None = not initialized, Some(true) = above threshold, Some(false) = below
}

impl Counter {
    /// Create a new bidirectional counter with start value and threshold
    pub fn new(start: i32, threshold: f64) -> Self {
        Self {
            input: 0.0,
            output: start as f64,
            start,
            threshold,
            count: 0,
            last_state: None,
        }
    }

    /// Get the start value
    pub fn start(&self) -> i32 {
        self.start
    }

    /// Get the threshold
    pub fn threshold(&self) -> f64 {
        self.threshold
    }

    /// Get the current count (number of events)
    pub fn count(&self) -> i32 {
        self.count
    }

    /// Get the event function value (input - threshold)
    pub fn event_value(&self) -> f64 {
        self.input - self.threshold
    }
}

impl Default for Counter {
    fn default() -> Self {
        Self::new(0, 0.0)
    }
}

impl Block for Counter {
    const NUM_INPUTS: usize = 1;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        std::slice::from_ref(&self.input)
    }

    fn inputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.input)
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, _t: f64) {
        let current_state = self.input >= self.threshold;

        // Check for zero crossing
        if let Some(last) = self.last_state {
            if last != current_state {
                // Zero crossing detected
                if current_state {
                    // Upward crossing
                    self.count += 1;
                } else {
                    // Downward crossing
                    self.count -= 1;
                }
            }
        }

        self.last_state = Some(current_state);

        // Update output: start + count
        self.output = (self.start + self.count) as f64;
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = self.start as f64;
        self.count = 0;
        self.last_state = None;
    }
}

impl AlgebraicBlock for Counter {}

/// Upward-only counter block
///
/// Counts only upward (low-to-high) threshold crossings.
///
/// # Example
///
/// ```ignore
/// let mut counter = CounterUp::new(0, 5.0);
/// counter.set_input(0, 10.0); // Crosses upward
/// counter.update(0.1);
/// assert_eq!(counter.get_output(0), 1);
///
/// counter.set_input(0, 0.0); // Crosses downward (not counted)
/// counter.update(0.2);
/// assert_eq!(counter.get_output(0), 1); // Still 1
/// ```
#[derive(Debug, Clone)]
pub struct CounterUp {
    input: f64,
    output: f64,
    start: i32,
    threshold: f64,
    count: i32,
    last_state: Option<bool>,
}

impl CounterUp {
    /// Create a new upward counter with start value and threshold
    pub fn new(start: i32, threshold: f64) -> Self {
        Self {
            input: 0.0,
            output: start as f64,
            start,
            threshold,
            count: 0,
            last_state: None,
        }
    }

    /// Get the start value
    pub fn start(&self) -> i32 {
        self.start
    }

    /// Get the threshold
    pub fn threshold(&self) -> f64 {
        self.threshold
    }

    /// Get the current count (number of events)
    pub fn count(&self) -> i32 {
        self.count
    }

    /// Get the event function value (input - threshold)
    pub fn event_value(&self) -> f64 {
        self.input - self.threshold
    }
}

impl Default for CounterUp {
    fn default() -> Self {
        Self::new(0, 0.0)
    }
}

impl Block for CounterUp {
    const NUM_INPUTS: usize = 1;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        std::slice::from_ref(&self.input)
    }

    fn inputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.input)
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, _t: f64) {
        let current_state = self.input >= self.threshold;

        // Check for upward zero crossing only
        if let Some(last) = self.last_state {
            if !last && current_state {
                // Upward crossing (low to high)
                self.count += 1;
            }
        }

        self.last_state = Some(current_state);

        // Update output: start + count
        self.output = (self.start + self.count) as f64;
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = self.start as f64;
        self.count = 0;
        self.last_state = None;
    }
}

impl AlgebraicBlock for CounterUp {}

/// Downward-only counter block
///
/// Counts only downward (high-to-low) threshold crossings.
///
/// # Example
///
/// ```ignore
/// let mut counter = CounterDown::new(0, 5.0);
/// counter.set_input(0, 10.0); // Start above threshold
/// counter.update(0.0);
///
/// counter.set_input(0, 0.0); // Crosses downward
/// counter.update(0.1);
/// assert_eq!(counter.get_output(0), 1);
///
/// counter.set_input(0, 10.0); // Crosses upward (not counted)
/// counter.update(0.2);
/// assert_eq!(counter.get_output(0), 1); // Still 1
/// ```
#[derive(Debug, Clone)]
pub struct CounterDown {
    input: f64,
    output: f64,
    start: i32,
    threshold: f64,
    count: i32,
    last_state: Option<bool>,
}

impl CounterDown {
    /// Create a new downward counter with start value and threshold
    pub fn new(start: i32, threshold: f64) -> Self {
        Self {
            input: 0.0,
            output: start as f64,
            start,
            threshold,
            count: 0,
            last_state: None,
        }
    }

    /// Get the start value
    pub fn start(&self) -> i32 {
        self.start
    }

    /// Get the threshold
    pub fn threshold(&self) -> f64 {
        self.threshold
    }

    /// Get the current count (number of events)
    pub fn count(&self) -> i32 {
        self.count
    }

    /// Get the event function value (input - threshold)
    pub fn event_value(&self) -> f64 {
        self.input - self.threshold
    }
}

impl Default for CounterDown {
    fn default() -> Self {
        Self::new(0, 0.0)
    }
}

impl Block for CounterDown {
    const NUM_INPUTS: usize = 1;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        std::slice::from_ref(&self.input)
    }

    fn inputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.input)
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, _t: f64) {
        let current_state = self.input >= self.threshold;

        // Check for downward zero crossing only
        if let Some(last) = self.last_state {
            if last && !current_state {
                // Downward crossing (high to low)
                self.count += 1;
            }
        }

        self.last_state = Some(current_state);

        // Update output: start + count
        self.output = (self.start + self.count) as f64;
    }

    fn reset(&mut self) {
        self.input = 0.0;
        self.output = self.start as f64;
        self.count = 0;
        self.last_state = None;
    }
}

impl AlgebraicBlock for CounterDown {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_counter_init() {
        // Default initialization
        let c = Counter::default();
        assert_eq!(c.start(), 0);
        assert_eq!(c.threshold(), 0.0);

        // Custom initialization
        let c = Counter::new(10, 5.0);
        assert_eq!(c.start(), 10);
        assert_eq!(c.threshold(), 5.0);
    }

    #[test]
    fn test_counter_initial_output() {
        let mut c = Counter::new(5, 0.0);
        c.set_input(0, 0.0);
        c.update(0.0);

        // Output should be start value
        assert_eq!(c.get_output(0), 5.0);
    }

    #[test]
    fn test_counter_output_formula() {
        let mut c = Counter::new(5, 0.0);

        // Initially no events, output = start
        c.update(0.0);
        assert_eq!(c.get_output(0), (5 + c.count()) as f64);

        // After initialization, count should be 0
        assert_eq!(c.get_output(0), 5.0);
    }

    #[test]
    fn test_counter_custom_start() {
        let start_values = [0, 10, -5, 100];

        for start_val in start_values {
            let mut c = Counter::new(start_val, 0.0);
            c.update(0.0);
            // Output should be start + number of events (initially 0)
            assert_eq!(c.get_output(0), start_val as f64);
        }
    }

    #[test]
    fn test_counter_threshold() {
        let mut c = Counter::new(0, 10.0);

        c.set_input(0, 15.0);
        let event_val = c.event_value();
        assert_eq!(event_val, 5.0); // 15 - 10

        c.set_input(0, 5.0);
        let event_val = c.event_value();
        assert_eq!(event_val, -5.0); // 5 - 10
    }

    #[test]
    fn test_counter_bidirectional() {
        let mut c = Counter::new(0, 5.0);

        // Start below threshold
        c.set_input(0, 0.0);
        c.update(0.0);
        assert_eq!(c.get_output(0), 0.0);

        // Cross upward
        c.set_input(0, 10.0);
        c.update(0.1);
        assert_eq!(c.get_output(0), 1.0);

        // Cross downward
        c.set_input(0, 0.0);
        c.update(0.2);
        assert_eq!(c.get_output(0), 0.0);

        // Cross upward again
        c.set_input(0, 10.0);
        c.update(0.3);
        assert_eq!(c.get_output(0), 1.0);
    }

    #[test]
    fn test_counterup_init() {
        let cu = CounterUp::new(5, 2.0);
        assert_eq!(cu.start(), 5);
        assert_eq!(cu.threshold(), 2.0);
    }

    #[test]
    fn test_counterup_update() {
        let mut cu = CounterUp::new(10, 5.0);

        cu.update(0.0);
        // Output should be start + number of events (initially 0)
        assert_eq!(cu.get_output(0), 10.0);

        // Start below threshold
        cu.set_input(0, 0.0);
        cu.update(0.0);
        assert_eq!(cu.get_output(0), 10.0);

        // Cross upward
        cu.set_input(0, 10.0);
        cu.update(0.1);
        assert_eq!(cu.get_output(0), 11.0);

        // Cross downward (should not count)
        cu.set_input(0, 0.0);
        cu.update(0.2);
        assert_eq!(cu.get_output(0), 11.0);

        // Cross upward again
        cu.set_input(0, 10.0);
        cu.update(0.3);
        assert_eq!(cu.get_output(0), 12.0);
    }

    #[test]
    fn test_counterdown_init() {
        let cd = CounterDown::new(5, 2.0);
        assert_eq!(cd.start(), 5);
        assert_eq!(cd.threshold(), 2.0);
    }

    #[test]
    fn test_counterdown_update() {
        let mut cd = CounterDown::new(20, 5.0);

        cd.update(0.0);
        // Output should be start + number of events (initially 0)
        assert_eq!(cd.get_output(0), 20.0);

        // Start above threshold
        cd.set_input(0, 10.0);
        cd.update(0.0);
        assert_eq!(cd.get_output(0), 20.0);

        // Cross downward
        cd.set_input(0, 0.0);
        cd.update(0.1);
        assert_eq!(cd.get_output(0), 21.0);

        // Cross upward (should not count)
        cd.set_input(0, 10.0);
        cd.update(0.2);
        assert_eq!(cd.get_output(0), 21.0);

        // Cross downward again
        cd.set_input(0, 0.0);
        cd.update(0.3);
        assert_eq!(cd.get_output(0), 22.0);
    }

    #[test]
    fn test_counter_reset() {
        let mut c = Counter::new(5, 2.0);

        c.set_input(0, 10.0);
        c.update(0.0);
        c.set_input(0, 0.0);
        c.update(0.1);

        c.reset();
        assert_eq!(c.get_input(0), 0.0);
        assert_eq!(c.get_output(0), 5.0);
        assert_eq!(c.count(), 0);
    }

    #[test]
    fn test_counterup_reset() {
        let mut cu = CounterUp::new(10, 5.0);

        cu.set_input(0, 10.0);
        cu.update(0.0);

        cu.reset();
        assert_eq!(cu.get_input(0), 0.0);
        assert_eq!(cu.get_output(0), 10.0);
        assert_eq!(cu.count(), 0);
    }

    #[test]
    fn test_counterdown_reset() {
        let mut cd = CounterDown::new(20, 5.0);

        cd.set_input(0, 10.0);
        cd.update(0.0);

        cd.reset();
        assert_eq!(cd.get_input(0), 0.0);
        assert_eq!(cd.get_output(0), 20.0);
        assert_eq!(cd.count(), 0);
    }
}
