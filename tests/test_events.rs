//! Comprehensive tests for event detection system
//!
//! Test parity with PathSim event tests

use rustsim::events::*;
use rustsim::utils::constants::EVT_TOLERANCE;

// ==================================================================================
// ZERO CROSSING TESTS
// ==================================================================================

#[test]
fn test_zerocrossing_detect_up() {
    // Upward crossing: event function crosses from negative to positive
    let mut e = ZeroCrossing::new(|t| t - 2.0, None::<fn(f64)>, EVT_TOLERANCE);

    // Before crossing
    e.buffer(0.0);
    let det = e.detect(1.0);
    assert!(!det.detected);
    assert!(!det.close);
    assert_eq!(det.ratio, 1.0);

    // Crossing 1: from t=1 to t=3, crosses at t=2
    // ratio = |prev_result| / |prev_result - result| = |-1| / |-1 - 1| = 1/2
    // But history is at t=1, result=-1
    // At t=3, result=1
    // ratio = 1 / 2 = 0.5, but PathSim expects 2/3
    // Actually the history should be at 0, not 1
    // Let me recalculate: buffer(0) -> result = -2
    // detect(3) -> result = 1
    // ratio = |-2| / |-2 - 1| = 2/3
    let det = e.detect(3.0);
    assert!(det.detected);
    assert!(!det.close);
    assert!((det.ratio - 2.0 / 3.0).abs() < 1e-10);

    // Crossing 2: buffer at t=1, detect at t=3
    // buffer(1) -> result = -1
    // detect(3) -> result = 1
    // ratio = |-1| / |-1 - 1| = 1/2
    e.buffer(1.0);
    let det = e.detect(3.0);
    assert!(det.detected);
    assert!(!det.close);
    assert!((det.ratio - 1.0 / 2.0).abs() < 1e-10);

    // After crossing
    e.buffer(3.0);
    let det = e.detect(4.0);
    assert!(!det.detected);
    assert!(!det.close);
    assert_eq!(det.ratio, 1.0);
}

#[test]
fn test_zerocrossing_detect_down() {
    // Downward crossing is also detected by bidirectional ZeroCrossing
    // Using a function that crosses from positive to negative
    // But the test uses the same function t-2, which crosses upward
    // PathSim test seems to be testing the same upward crossing
    // Let me check the Python test again - it's actually the same test!
    // So test_detect_down for ZeroCrossing is identical to test_detect_up

    let mut e = ZeroCrossing::new(|t| t - 2.0, None::<fn(f64)>, EVT_TOLERANCE);

    // Before crossing
    e.buffer(0.0);
    let det = e.detect(1.0);
    assert!(!det.detected);
    assert!(!det.close);
    assert_eq!(det.ratio, 1.0);

    // Crossing 1
    let det = e.detect(3.0);
    assert!(det.detected);
    assert!(!det.close);
    assert!((det.ratio - 2.0 / 3.0).abs() < 1e-10);

    // Crossing 2
    e.buffer(1.0);
    let det = e.detect(3.0);
    assert!(det.detected);
    assert!(!det.close);
    assert!((det.ratio - 1.0 / 2.0).abs() < 1e-10);

    // After crossing
    e.buffer(3.0);
    let det = e.detect(4.0);
    assert!(!det.detected);
    assert!(!det.close);
    assert_eq!(det.ratio, 1.0);
}

#[test]
fn test_zerocrossingup_detect_up() {
    // Upward only: should trigger on upward crossing
    let mut e = ZeroCrossingUp::new(|t| t - 2.0, None::<fn(f64)>, EVT_TOLERANCE);

    // Before crossing
    e.buffer(0.0);
    let det = e.detect(1.0);
    assert!(!det.detected);
    assert!(!det.close);
    assert_eq!(det.ratio, 1.0);

    // Crossing 1
    let det = e.detect(3.0);
    assert!(det.detected);
    assert!(!det.close);
    assert!((det.ratio - 2.0 / 3.0).abs() < 1e-10);

    // Crossing 2
    e.buffer(1.0);
    let det = e.detect(3.0);
    assert!(det.detected);
    assert!(!det.close);
    assert!((det.ratio - 1.0 / 2.0).abs() < 1e-10);

    // After crossing
    e.buffer(3.0);
    let det = e.detect(4.0);
    assert!(!det.detected);
    assert!(!det.close);
    assert_eq!(det.ratio, 1.0);
}

#[test]
fn test_zerocrossingup_detect_down() {
    // Upward only: should NOT trigger on downward crossing
    // However, the test function t-2 crosses upward, not downward
    // The test just checks that it doesn't trigger after the crossing
    let mut e = ZeroCrossingUp::new(|t| t - 2.0, None::<fn(f64)>, EVT_TOLERANCE);

    // Before crossing
    e.buffer(0.0);
    let det = e.detect(1.0);
    assert!(!det.detected);
    assert!(!det.close);
    assert_eq!(det.ratio, 1.0);

    // After crossing
    e.buffer(3.0);
    let det = e.detect(4.0);
    assert!(!det.detected);
    assert!(!det.close);
    assert_eq!(det.ratio, 1.0);
}

#[test]
fn test_zerocrossingdown_detect_up() {
    // Downward only: should NOT trigger on upward crossing
    let mut e = ZeroCrossingDown::new(|t| t - 2.0, None::<fn(f64)>, EVT_TOLERANCE);

    // Before crossing
    e.buffer(0.0);
    let det = e.detect(1.0);
    assert!(!det.detected);
    assert!(!det.close);
    assert_eq!(det.ratio, 1.0);

    // Crossing 1 -> not triggering
    let det = e.detect(3.0);
    assert!(!det.detected);
    assert!(!det.close);
    assert_eq!(det.ratio, 1.0);

    // Crossing 2 -> not triggering
    e.buffer(1.0);
    let det = e.detect(3.0);
    assert!(!det.detected);
    assert!(!det.close);
    assert_eq!(det.ratio, 1.0);

    // After crossing
    e.buffer(3.0);
    let det = e.detect(4.0);
    assert!(!det.detected);
    assert!(!det.close);
    assert_eq!(det.ratio, 1.0);
}

#[test]
fn test_zerocrossingdown_detect_down() {
    // Downward only: should trigger on downward crossing
    // But the function t-2 crosses upward, so this won't trigger
    // The PathSim test just checks it doesn't trigger
    let mut e = ZeroCrossingDown::new(|t| t - 2.0, None::<fn(f64)>, EVT_TOLERANCE);

    // Before crossing
    e.buffer(0.0);
    let det = e.detect(1.0);
    assert!(!det.detected);
    assert!(!det.close);
    assert_eq!(det.ratio, 1.0);

    // After crossing
    e.buffer(3.0);
    let det = e.detect(4.0);
    assert!(!det.detected);
    assert!(!det.close);
    assert_eq!(det.ratio, 1.0);
}

// ==================================================================================
// CONDITION TESTS
// ==================================================================================

#[test]
fn test_condition_detect_false() {
    // Test detect when condition is false
    let mut e = Condition::new(|t| t > 5.0, None::<fn(f64)>, 0.1);

    // Before condition is met
    e.buffer(1.0);
    let det = e.detect(2.0);

    assert!(!det.detected);
    assert!(!det.close);
    assert_eq!(det.ratio, 0.5); // Bisection halves the step
}

#[test]
fn test_condition_detect_true_not_close() {
    // Test detect when condition becomes true but not close enough
    let mut e = Condition::new(|t| t > 5.0, None::<fn(f64)>, 0.1);

    // Condition just became true, but time gap is large
    e.buffer(4.0);
    let det = e.detect(6.0);

    assert!(det.detected);
    assert!(!det.close); // 6.0 - 4.0 = 2.0 > tolerance
    assert_eq!(det.ratio, 0.5);
}

#[test]
fn test_condition_detect_true_close() {
    // Test detect when condition is true and we're close enough
    let mut e = Condition::new(|t| t > 5.0, None::<fn(f64)>, 0.1);

    // Condition true and time gap is small
    e.buffer(5.05);
    let det = e.detect(5.08);

    assert!(det.detected);
    assert!(det.close); // 5.08 - 5.05 = 0.03 < tolerance
    assert_eq!(det.ratio, 1.0);
}

#[test]
fn test_condition_resolve_without_action() {
    // Test resolve without action function
    let mut e = Condition::new(|t| t > 5.0, None::<fn(f64)>, 0.1);

    // Initially active
    assert!(e.is_active());

    // Resolve event
    e.resolve(5.5);

    // Should be deactivated after resolution
    assert!(!e.is_active());

    // Event time should be recorded
    assert_eq!(e.len(), 1);
    assert_eq!(e.event_times()[0], 5.5);
}

#[test]
fn test_condition_resolve_with_action() {
    // Test resolve with action function
    use std::sync::{Arc, Mutex};

    let action_called = Arc::new(Mutex::new(Vec::new()));
    let action_called_clone = action_called.clone();

    let mut e = Condition::new(
        |t| t > 5.0,
        Some(move |t| {
            action_called_clone.lock().unwrap().push(t);
        }),
        0.1,
    );

    // Initially active
    assert!(e.is_active());

    // Resolve event
    e.resolve(5.5);

    // Action should have been called
    let called = action_called.lock().unwrap();
    assert_eq!(called.len(), 1);
    assert_eq!(called[0], 5.5);

    // Should be deactivated
    assert!(!e.is_active());
}

#[test]
fn test_condition_bisection() {
    // Test that bisection narrows down to tolerance
    let mut e = Condition::new(|t| t > 10.0, None::<fn(f64)>, 0.1);

    // Start before condition (large gap)
    e.buffer(9.0);
    let det = e.detect(11.0);
    assert!(det.detected);
    assert!(!det.close); // 11.0 - 9.0 = 2.0 > tolerance
    assert_eq!(det.ratio, 0.5);

    // Narrow down (small gap)
    e.buffer(10.05);
    let det = e.detect(10.08);
    assert!(det.detected);
    assert!(det.close); // 10.08 - 10.05 = 0.03 < tolerance
    assert_eq!(det.ratio, 1.0);
}

// ==================================================================================
// SCHEDULE TESTS
// ==================================================================================

#[test]
fn test_schedule_next() {
    let s = Schedule::new(0.0, None, 20.0, None::<fn(f64)>, EVT_TOLERANCE);

    assert_eq!(s.next(), 0.0);

    let mut s = s;
    s.resolve(0.0);

    assert_eq!(s.next(), 20.0);
}

#[test]
fn test_schedule_estimate() {
    let s = Schedule::new(2.0, None, 20.0, None::<fn(f64)>, EVT_TOLERANCE);

    assert_eq!(s.estimate(0.0), Some(2.0));
    assert_eq!(s.estimate(1.0), Some(1.0));

    let mut s = s;
    s.resolve(2.0);

    assert_eq!(s.estimate(2.0), Some(20.0));
    assert_eq!(s.estimate(13.0), Some(9.0));
}

#[test]
fn test_schedule_detect() {
    let mut s = Schedule::new(2.0, None, 20.0, None::<fn(f64)>, EVT_TOLERANCE);

    s.buffer(0.0);

    let det = s.detect(0.0);
    assert!(!det.detected);
    assert!(!det.close);

    let det = s.detect(4.0);
    assert!(det.detected);
    assert!(!det.close);
    assert!((det.ratio - 0.5).abs() < 1e-10);
}

#[test]
fn test_schedulelist_init() {
    // Test that creation fails with non-ascending times
    let result = ScheduleList::new(
        vec![1.0, 3.0, 5.0, 2.0, 7.0],
        None::<fn(f64)>,
        EVT_TOLERANCE,
    );
    assert!(result.is_err());

    // Test successful creation with ascending times
    let result = ScheduleList::new(vec![1.0, 3.0, 5.0, 7.0], None::<fn(f64)>, EVT_TOLERANCE);
    assert!(result.is_ok());

    let s = result.unwrap();
    assert_eq!(s.times_evt, vec![1.0, 3.0, 5.0, 7.0]);
}

#[test]
fn test_schedulelist_next() {
    let s = ScheduleList::new(vec![1.0, 3.0, 5.0, 7.0], None::<fn(f64)>, EVT_TOLERANCE).unwrap();

    assert_eq!(s.next(), 1.0);

    let mut s = s;
    s.resolve(1.0);

    assert_eq!(s.next(), 3.0);

    s.resolve(3.0);

    assert_eq!(s.next(), 5.0);
}

#[test]
fn test_schedulelist_detect() {
    let mut s =
        ScheduleList::new(vec![1.0, 3.0, 5.0, 7.0], None::<fn(f64)>, EVT_TOLERANCE).unwrap();

    s.buffer(0.0);

    let det = s.detect(0.0);
    assert!(!det.detected);
    assert!(!det.close);

    let det = s.detect(2.0);
    assert!(det.detected);
    assert!(!det.close);
    assert!((det.ratio - 0.5).abs() < 1e-10);
}

#[test]
fn test_schedulelist_estimate() {
    let s = ScheduleList::new(vec![1.0, 3.0, 5.0, 7.0], None::<fn(f64)>, EVT_TOLERANCE).unwrap();

    assert_eq!(s.estimate(0.0), Some(1.0));
    assert_eq!(s.estimate(0.5), Some(0.5));

    let mut s = s;
    s.resolve(1.0);

    assert_eq!(s.estimate(1.0), Some(2.0));
    assert_eq!(s.estimate(2.0), Some(1.0));
}

#[test]
fn test_schedulelist_func_act_is_not_none() {
    use std::sync::{Arc, Mutex};

    let called = Arc::new(Mutex::new(false));
    let called_clone = called.clone();

    let event = ScheduleList::new(
        vec![1.0, 2.0, 3.0],
        Some(move |_| {
            *called_clone.lock().unwrap() = true;
        }),
        EVT_TOLERANCE,
    )
    .unwrap();

    // func_act should not be None
    assert!(event.func_act.is_some());
}

#[test]
fn test_condition_len() {
    // Test that len() returns number of times event occurred
    let mut e = Condition::new(|t| t > 5.0, None::<fn(f64)>, 0.1);

    assert_eq!(e.len(), 0);

    e.resolve(5.5);
    assert_eq!(e.len(), 1);

    e.resolve(10.5);
    assert_eq!(e.len(), 2);
}

// Additional tests for coverage and edge cases

#[test]
fn test_zerocrossing_exactly_zero() {
    let mut e = ZeroCrossing::new(|t| t - 2.0, None::<fn(f64)>, EVT_TOLERANCE);

    e.buffer(1.0);
    let det = e.detect(2.0); // Exactly hits zero

    assert!(det.detected);
    assert!(det.close);
    assert_eq!(det.ratio, 1.0);
}

#[test]
fn test_zerocrossingup_exactly_zero() {
    let mut e = ZeroCrossingUp::new(|t| t - 2.0, None::<fn(f64)>, EVT_TOLERANCE);

    e.buffer(1.0);
    let det = e.detect(2.0); // Exactly hits zero from negative

    assert!(det.detected);
    assert!(det.close);
    assert_eq!(det.ratio, 1.0);
}

#[test]
fn test_zerocrossingdown_exactly_zero() {
    let mut e = ZeroCrossingDown::new(|t| 2.0 - t, None::<fn(f64)>, EVT_TOLERANCE);

    e.buffer(1.0);
    let det = e.detect(2.0); // Exactly hits zero from positive

    assert!(det.detected);
    assert!(det.close);
    assert_eq!(det.ratio, 1.0);
}

#[test]
fn test_event_reset() {
    let mut e = ZeroCrossing::new(|t| t - 2.0, None::<fn(f64)>, EVT_TOLERANCE);

    e.buffer(0.0);
    e.detect(3.0);
    e.resolve(2.0);

    assert_eq!(e.len(), 1);

    e.reset();

    assert_eq!(e.len(), 0);
    assert!(e.is_active());
}

#[test]
fn test_schedule_with_end_time() {
    let mut s = Schedule::new(0.0, Some(50.0), 20.0, None::<fn(f64)>, EVT_TOLERANCE);

    // First two events should trigger
    s.buffer(0.0);
    s.detect(1.0);
    s.resolve(0.0);

    s.buffer(10.0);
    s.detect(25.0);
    s.resolve(20.0);

    // Third event (at t=40) should trigger
    s.buffer(30.0);
    let det = s.detect(45.0);
    assert!(det.detected);
    s.resolve(40.0);

    // Fourth event would be at t=60, which is > t_end=50
    s.buffer(50.0);
    let det = s.detect(65.0);
    assert!(!det.detected);
    assert!(!s.is_active());
}

#[test]
fn test_schedulelist_deactivates_after_all_events() {
    let mut s = ScheduleList::new(vec![1.0, 2.0], None::<fn(f64)>, EVT_TOLERANCE).unwrap();

    s.buffer(0.0);
    s.detect(1.5);
    s.resolve(1.0);

    s.buffer(1.5);
    s.detect(2.5);
    s.resolve(2.0);

    // All events processed, should deactivate
    s.buffer(2.5);
    let det = s.detect(3.0);
    assert!(!det.detected);
    assert!(!s.is_active());
}

#[test]
fn test_condition_on_off() {
    let mut e = Condition::new(|t| t > 5.0, None::<fn(f64)>, 0.1);

    assert!(e.is_active());

    e.off();
    assert!(!e.is_active());

    e.on();
    assert!(e.is_active());
}

#[test]
fn test_zerocrossing_with_action() {
    use std::sync::{Arc, Mutex};

    let called = Arc::new(Mutex::new(Vec::new()));
    let called_clone = called.clone();

    let mut e = ZeroCrossing::new(
        |t| t - 2.0,
        Some(move |t| {
            called_clone.lock().unwrap().push(t);
        }),
        EVT_TOLERANCE,
    );

    e.buffer(0.0);
    e.detect(3.0);
    e.resolve(2.0);

    let times = called.lock().unwrap();
    assert_eq!(times.len(), 1);
    assert_eq!(times[0], 2.0);
}
