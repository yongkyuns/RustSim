//! Signal source blocks (zero inputs, one or more outputs)

use crate::block::{AlgebraicBlock, Block};
use std::f64::consts::PI;

/// Constant signal source
///
/// # Example
///
/// ```ignore
/// let c = Constant::new(5.0);
/// assert_eq!(c.get_output(0), 5.0);
/// ```
#[derive(Debug, Clone)]
pub struct Constant {
    output: f64,
    value: f64,
}

impl Constant {
    pub fn new(value: f64) -> Self {
        Self {
            output: value,
            value,
        }
    }

    pub fn value(&self) -> f64 {
        self.value
    }

    pub fn set_value(&mut self, value: f64) {
        self.value = value;
        self.output = value;
    }
}

impl Block for Constant {
    const NUM_INPUTS: usize = 0;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        &[]
    }
    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut []
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, _t: f64) {
        self.output = self.value;
    }

    fn reset(&mut self) {
        self.output = self.value;
    }
}

impl AlgebraicBlock for Constant {}

/// Sinusoidal signal source: y = amplitude * sin(2*pi*frequency*t + phase)
#[derive(Debug, Clone)]
pub struct Sinusoidal {
    output: f64,
    amplitude: f64,
    frequency: f64,
    phase: f64,
}

impl Sinusoidal {
    /// Create sinusoid with amplitude, frequency (Hz), and phase (radians)
    pub fn new(amplitude: f64, frequency: f64, phase: f64) -> Self {
        Self {
            output: amplitude * phase.sin(),
            amplitude,
            frequency,
            phase,
        }
    }

    /// Create unit sinusoid at given frequency
    pub fn unit(frequency: f64) -> Self {
        Self::new(1.0, frequency, 0.0)
    }

    pub fn amplitude(&self) -> f64 {
        self.amplitude
    }
    pub fn frequency(&self) -> f64 {
        self.frequency
    }
    pub fn phase(&self) -> f64 {
        self.phase
    }
}

impl Block for Sinusoidal {
    const NUM_INPUTS: usize = 0;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        &[]
    }
    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut []
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, t: f64) {
        self.output = self.amplitude * (2.0 * PI * self.frequency * t + self.phase).sin();
    }

    fn reset(&mut self) {
        self.output = self.amplitude * self.phase.sin();
    }
}

impl AlgebraicBlock for Sinusoidal {}

/// Step signal source: y = 0 for t < t0, y = amplitude for t >= t0
#[derive(Debug, Clone)]
pub struct Step {
    output: f64,
    amplitude: f64,
    t0: f64,
}

impl Step {
    /// Create step at time t0 with given amplitude
    pub fn new(amplitude: f64, t0: f64) -> Self {
        Self {
            output: if t0 <= 0.0 { amplitude } else { 0.0 },
            amplitude,
            t0,
        }
    }

    /// Create unit step at t=0
    pub fn unit() -> Self {
        Self::new(1.0, 0.0)
    }

    pub fn amplitude(&self) -> f64 {
        self.amplitude
    }
    pub fn t0(&self) -> f64 {
        self.t0
    }
}

impl Block for Step {
    const NUM_INPUTS: usize = 0;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        &[]
    }
    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut []
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, t: f64) {
        self.output = if t >= self.t0 { self.amplitude } else { 0.0 };
    }

    fn reset(&mut self) {
        self.output = if self.t0 <= 0.0 { self.amplitude } else { 0.0 };
    }
}

impl AlgebraicBlock for Step {}

/// Ramp signal source: y = slope * (t - t0) for t >= t0, else 0
#[derive(Debug, Clone)]
pub struct Ramp {
    output: f64,
    slope: f64,
    t0: f64,
}

impl Ramp {
    /// Create ramp starting at time t0 with given slope
    pub fn new(slope: f64, t0: f64) -> Self {
        Self {
            output: 0.0,
            slope,
            t0,
        }
    }

    /// Create unit ramp starting at t=0
    pub fn unit() -> Self {
        Self::new(1.0, 0.0)
    }

    pub fn slope(&self) -> f64 {
        self.slope
    }
    pub fn t0(&self) -> f64 {
        self.t0
    }
}

impl Block for Ramp {
    const NUM_INPUTS: usize = 0;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        &[]
    }
    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut []
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, t: f64) {
        self.output = if t >= self.t0 {
            self.slope * (t - self.t0)
        } else {
            0.0
        };
    }

    fn reset(&mut self) {
        self.output = 0.0;
    }
}

impl AlgebraicBlock for Ramp {}

/// Square wave signal source with configurable duty cycle
#[derive(Debug, Clone)]
pub struct SquareWave {
    output: f64,
    amplitude: f64,
    frequency: f64,
    duty_cycle: f64, // 0.0 to 1.0
}

impl SquareWave {
    /// Create square wave with amplitude, frequency (Hz), and duty cycle (0-1)
    pub fn new(amplitude: f64, frequency: f64, duty_cycle: f64) -> Self {
        assert!(
            (0.0..=1.0).contains(&duty_cycle),
            "Duty cycle must be in [0, 1]"
        );
        Self {
            output: amplitude, // Start high
            amplitude,
            frequency,
            duty_cycle,
        }
    }

    /// Create square wave with 50% duty cycle
    pub fn symmetric(amplitude: f64, frequency: f64) -> Self {
        Self::new(amplitude, frequency, 0.5)
    }

    pub fn amplitude(&self) -> f64 {
        self.amplitude
    }
    pub fn frequency(&self) -> f64 {
        self.frequency
    }
    pub fn duty_cycle(&self) -> f64 {
        self.duty_cycle
    }
}

impl Block for SquareWave {
    const NUM_INPUTS: usize = 0;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        &[]
    }
    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut []
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, t: f64) {
        let period = 1.0 / self.frequency;
        let phase = (t % period) / period;
        self.output = if phase < self.duty_cycle {
            self.amplitude
        } else {
            -self.amplitude
        };
    }

    fn reset(&mut self) {
        self.output = self.amplitude;
    }
}

impl AlgebraicBlock for SquareWave {}

/// Triangle wave signal source
#[derive(Debug, Clone)]
pub struct TriangleWave {
    output: f64,
    amplitude: f64,
    frequency: f64,
}

impl TriangleWave {
    /// Create triangle wave with amplitude and frequency (Hz)
    pub fn new(amplitude: f64, frequency: f64) -> Self {
        Self {
            output: 0.0,
            amplitude,
            frequency,
        }
    }

    pub fn amplitude(&self) -> f64 {
        self.amplitude
    }
    pub fn frequency(&self) -> f64 {
        self.frequency
    }
}

impl Block for TriangleWave {
    const NUM_INPUTS: usize = 0;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        &[]
    }
    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut []
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, t: f64) {
        let period = 1.0 / self.frequency;
        let phase = (t % period) / period;
        // Triangle: 0 -> 1 -> 0 -> -1 -> 0
        self.output = if phase < 0.25 {
            4.0 * phase * self.amplitude
        } else if phase < 0.75 {
            (2.0 - 4.0 * phase) * self.amplitude
        } else {
            (4.0 * phase - 4.0) * self.amplitude
        };
    }

    fn reset(&mut self) {
        self.output = 0.0;
    }
}

impl AlgebraicBlock for TriangleWave {}

/// Periodic pulse signal with rise/fall times
#[derive(Debug, Clone)]
pub struct Pulse {
    output: f64,
    amplitude: f64,
    period: f64,
    pulse_width: f64,
    rise_time: f64,
    fall_time: f64,
}

impl Pulse {
    /// Create pulse with amplitude, period, pulse_width, rise_time, fall_time
    pub fn new(
        amplitude: f64,
        period: f64,
        pulse_width: f64,
        rise_time: f64,
        fall_time: f64,
    ) -> Self {
        assert!(
            pulse_width > 0.0 && pulse_width < period,
            "Pulse width must be in (0, period)"
        );
        assert!(
            rise_time >= 0.0 && fall_time >= 0.0,
            "Rise and fall times must be non-negative"
        );
        assert!(
            rise_time + fall_time < pulse_width,
            "Rise + fall time must be less than pulse width"
        );

        Self {
            output: 0.0,
            amplitude,
            period,
            pulse_width,
            rise_time,
            fall_time,
        }
    }

    /// Create ideal pulse (zero rise/fall time)
    pub fn ideal(amplitude: f64, period: f64, pulse_width: f64) -> Self {
        Self::new(amplitude, period, pulse_width, 0.0, 0.0)
    }

    pub fn amplitude(&self) -> f64 {
        self.amplitude
    }
    pub fn period(&self) -> f64 {
        self.period
    }
    pub fn pulse_width(&self) -> f64 {
        self.pulse_width
    }
}

impl Block for Pulse {
    const NUM_INPUTS: usize = 0;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        &[]
    }
    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut []
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, t: f64) {
        let t_in_period = t % self.period;

        if t_in_period < self.rise_time {
            // Rising edge
            self.output = self.amplitude * (t_in_period / self.rise_time);
        } else if t_in_period < self.pulse_width - self.fall_time {
            // High plateau
            self.output = self.amplitude;
        } else if t_in_period < self.pulse_width {
            // Falling edge
            let fall_progress =
                (t_in_period - (self.pulse_width - self.fall_time)) / self.fall_time;
            self.output = self.amplitude * (1.0 - fall_progress);
        } else {
            // Low
            self.output = 0.0;
        }
    }

    fn reset(&mut self) {
        self.output = 0.0;
    }
}

impl AlgebraicBlock for Pulse {}

/// Digital clock signal with configurable period and duty cycle
#[derive(Debug, Clone)]
pub struct Clock {
    output: f64,
    period: f64,
    duty_cycle: f64, // 0.0 to 1.0
}

impl Clock {
    /// Create clock with period and duty cycle (0-1)
    pub fn new(period: f64, duty_cycle: f64) -> Self {
        assert!(period > 0.0, "Period must be positive");
        assert!(
            (0.0..=1.0).contains(&duty_cycle),
            "Duty cycle must be in [0, 1]"
        );
        Self {
            output: 1.0, // Start high
            period,
            duty_cycle,
        }
    }

    /// Create symmetric clock (50% duty cycle)
    pub fn symmetric(period: f64) -> Self {
        Self::new(period, 0.5)
    }

    pub fn period(&self) -> f64 {
        self.period
    }
    pub fn duty_cycle(&self) -> f64 {
        self.duty_cycle
    }
}

impl Block for Clock {
    const NUM_INPUTS: usize = 0;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        &[]
    }
    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut []
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, t: f64) {
        let phase = (t % self.period) / self.period;
        self.output = if phase < self.duty_cycle { 1.0 } else { 0.0 };
    }

    fn reset(&mut self) {
        self.output = 1.0;
    }
}

impl AlgebraicBlock for Clock {}

/// Gaussian pulse signal source
///
/// Generates a Gaussian-shaped pulse: y = amplitude * exp(-((t-tau)/sigma)^2)
/// where sigma = 0.5 / f_max
#[derive(Debug, Clone)]
pub struct GaussianPulse {
    output: f64,
    amplitude: f64,
    f_max: f64,
    tau: f64,
}

impl GaussianPulse {
    /// Create Gaussian pulse with amplitude, bandwidth (f_max), and center time (tau)
    ///
    /// # Parameters
    /// - `amplitude`: Peak amplitude of the pulse
    /// - `f_max`: Maximum frequency component (controls pulse width, sigma = 0.5 / f_max)
    /// - `tau`: Center time of the pulse (peak location)
    pub fn new(amplitude: f64, f_max: f64, tau: f64) -> Self {
        Self {
            output: 0.0,
            amplitude,
            f_max,
            tau,
        }
    }

    /// Create unit Gaussian pulse at t=0 with f_max=1e3
    pub fn unit() -> Self {
        Self::new(1.0, 1e3, 0.0)
    }

    pub fn amplitude(&self) -> f64 {
        self.amplitude
    }
    pub fn f_max(&self) -> f64 {
        self.f_max
    }
    pub fn tau(&self) -> f64 {
        self.tau
    }
}

impl Block for GaussianPulse {
    const NUM_INPUTS: usize = 0;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        &[]
    }
    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut []
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, t: f64) {
        let sigma = 0.5 / self.f_max;
        let t_shifted = t - self.tau;
        self.output = self.amplitude * (-(t_shifted / sigma).powi(2)).exp();
    }

    fn reset(&mut self) {
        let sigma = 0.5 / self.f_max;
        let t_shifted = -self.tau;
        self.output = self.amplitude * (-(t_shifted / sigma).powi(2)).exp();
    }
}

impl AlgebraicBlock for GaussianPulse {}

/// Chirp signal source (linear frequency sweep)
///
/// Generates a chirp signal with linearly varying frequency from f0 to f0+BW over period T.
/// The instantaneous frequency varies as a triangle wave, sweeping up and down.
#[derive(Debug, Clone)]
pub struct Chirp {
    output: f64,
    amplitude: f64,
    f0: f64,
    bandwidth: f64,
    sweep_time: f64,
    phase: f64,
}

impl Chirp {
    /// Create chirp with start frequency, bandwidth, and sweep period
    ///
    /// # Parameters
    /// - `amplitude`: Signal amplitude
    /// - `f0`: Start frequency (Hz)
    /// - `bandwidth`: Frequency sweep range (Hz)
    /// - `sweep_time`: Period of the frequency sweep (s)
    /// - `phase`: Initial phase offset (radians)
    pub fn new(amplitude: f64, f0: f64, bandwidth: f64, sweep_time: f64, phase: f64) -> Self {
        Self {
            output: amplitude * phase.sin(),
            amplitude,
            f0,
            bandwidth,
            sweep_time,
            phase,
        }
    }

    /// Create unit chirp with default parameters
    pub fn unit() -> Self {
        Self::new(1.0, 1.0, 1.0, 1.0, 0.0)
    }

    pub fn amplitude(&self) -> f64 {
        self.amplitude
    }
    pub fn f0(&self) -> f64 {
        self.f0
    }
    pub fn bandwidth(&self) -> f64 {
        self.bandwidth
    }
    pub fn sweep_time(&self) -> f64 {
        self.sweep_time
    }
    pub fn phase(&self) -> f64 {
        self.phase
    }

    /// Triangle wave helper: -1 to +1 with frequency f
    fn triangle_wave(t: f64, f: f64) -> f64 {
        2.0 * (2.0 * ((t * f) % 1.0) - 1.0).abs() - 1.0
    }

    /// Compute instantaneous frequency at time t
    pub fn instantaneous_frequency(&self, t: f64) -> f64 {
        let tri = Self::triangle_wave(t, 1.0 / self.sweep_time);
        self.f0 + self.bandwidth * (1.0 + tri) / 2.0
    }
}

impl Block for Chirp {
    const NUM_INPUTS: usize = 0;
    const NUM_OUTPUTS: usize = 1;
    const IS_DYNAMIC: bool = false;

    fn inputs(&self) -> &[f64] {
        &[]
    }
    fn inputs_mut(&mut self) -> &mut [f64] {
        &mut []
    }

    fn outputs(&self) -> &[f64] {
        std::slice::from_ref(&self.output)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        std::slice::from_mut(&mut self.output)
    }

    fn update(&mut self, t: f64) {
        // Simplified linear chirp: instantaneous frequency sweeps linearly
        // f(t) = f0 + (BW / T) * t (modulo sweep period)
        // Phase is integral of frequency

        let t_mod = t % self.sweep_time;
        let chirp_rate = self.bandwidth / self.sweep_time;

        // Phase from linear chirp formula: phi = 2*pi*(f0*t + 0.5*rate*t^2)
        let phase_integral = self.f0 * t_mod + 0.5 * chirp_rate * t_mod * t_mod;

        self.output = self.amplitude * (2.0 * PI * phase_integral + self.phase).sin();
    }

    fn reset(&mut self) {
        self.output = self.amplitude * self.phase.sin();
    }
}

impl AlgebraicBlock for Chirp {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constant() {
        let mut c = Constant::new(42.0);
        c.update(100.0); // Time doesn't matter
        assert_eq!(c.get_output(0), 42.0);
    }

    #[test]
    fn test_sinusoidal() {
        let mut s = Sinusoidal::new(2.0, 1.0, 0.0);

        // At t=0, sin(0) = 0
        s.update(0.0);
        assert!((s.get_output(0) - 0.0).abs() < 1e-10);

        // At t=0.25, sin(pi/2) = 1
        s.update(0.25);
        assert!((s.get_output(0) - 2.0).abs() < 1e-10);

        // At t=0.5, sin(pi) = 0
        s.update(0.5);
        assert!((s.get_output(0) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_step() {
        let mut s = Step::new(5.0, 1.0);

        s.update(0.5);
        assert_eq!(s.get_output(0), 0.0);

        s.update(1.0);
        assert_eq!(s.get_output(0), 5.0);

        s.update(2.0);
        assert_eq!(s.get_output(0), 5.0);
    }

    #[test]
    fn test_ramp() {
        let mut r = Ramp::new(2.0, 1.0);

        r.update(0.5);
        assert_eq!(r.get_output(0), 0.0);

        r.update(1.0);
        assert_eq!(r.get_output(0), 0.0);

        r.update(2.0);
        assert_eq!(r.get_output(0), 2.0);

        r.update(3.0);
        assert_eq!(r.get_output(0), 4.0);
    }

    #[test]
    fn test_square_wave() {
        let mut sq = SquareWave::new(1.0, 1.0, 0.5);

        // At t=0, should be high
        sq.update(0.0);
        assert_eq!(sq.get_output(0), 1.0);

        // At t=0.25 (first quarter), still high
        sq.update(0.25);
        assert_eq!(sq.get_output(0), 1.0);

        // At t=0.75 (third quarter), should be low
        sq.update(0.75);
        assert_eq!(sq.get_output(0), -1.0);

        // At t=1.0 (new period), should be high again
        sq.update(1.0);
        assert_eq!(sq.get_output(0), 1.0);
    }

    #[test]
    fn test_triangle_wave() {
        let mut tri = TriangleWave::new(2.0, 1.0);

        // At t=0, should be 0
        tri.update(0.0);
        assert!((tri.get_output(0) - 0.0).abs() < 1e-10);

        // At t=0.25 (peak), should be amplitude
        tri.update(0.25);
        assert!((tri.get_output(0) - 2.0).abs() < 1e-10);

        // At t=0.5, should be 0
        tri.update(0.5);
        assert!((tri.get_output(0) - 0.0).abs() < 1e-10);

        // At t=0.75 (trough), should be -amplitude
        tri.update(0.75);
        assert!((tri.get_output(0) + 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_pulse_ideal() {
        let mut p = Pulse::ideal(5.0, 10.0, 3.0);

        // Before pulse
        p.update(0.5);
        assert_eq!(p.get_output(0), 5.0);

        // During pulse
        p.update(2.0);
        assert_eq!(p.get_output(0), 5.0);

        // After pulse
        p.update(5.0);
        assert_eq!(p.get_output(0), 0.0);

        // Second period
        p.update(11.0);
        assert_eq!(p.get_output(0), 5.0);
    }

    #[test]
    fn test_pulse_with_edges() {
        let mut p = Pulse::new(10.0, 10.0, 5.0, 1.0, 1.0);

        // During rise
        p.update(0.5);
        assert_eq!(p.get_output(0), 5.0);

        // At plateau
        p.update(2.0);
        assert_eq!(p.get_output(0), 10.0);

        // During fall (at t=4.0, fall starts; at t=4.5, half-way down)
        p.update(4.5);
        assert_eq!(p.get_output(0), 5.0);

        // After pulse
        p.update(6.0);
        assert_eq!(p.get_output(0), 0.0);
    }

    #[test]
    fn test_clock() {
        let mut clk = Clock::symmetric(2.0);

        // High phase
        clk.update(0.5);
        assert_eq!(clk.get_output(0), 1.0);

        // Low phase
        clk.update(1.5);
        assert_eq!(clk.get_output(0), 0.0);

        // Second period, high again
        clk.update(2.5);
        assert_eq!(clk.get_output(0), 1.0);
    }

    #[test]
    fn test_clock_asymmetric() {
        let mut clk = Clock::new(10.0, 0.3);

        // High phase (first 30%)
        clk.update(2.0);
        assert_eq!(clk.get_output(0), 1.0);

        // Low phase (last 70%)
        clk.update(5.0);
        assert_eq!(clk.get_output(0), 0.0);
    }

    // GaussianPulse tests - matching PathSim test_sources.py

    #[test]
    fn test_gaussian_pulse_init() {
        // Default
        let gp = GaussianPulse::new(1.0, 1e3, 0.0);
        assert_eq!(gp.amplitude(), 1.0);
        assert_eq!(gp.f_max(), 1e3);
        assert_eq!(gp.tau(), 0.0);

        // Specific
        let gp2 = GaussianPulse::new(5.0, 2e3, 1.0);
        assert_eq!(gp2.amplitude(), 5.0);
        assert_eq!(gp2.f_max(), 2e3);
        assert_eq!(gp2.tau(), 1.0);
    }

    #[test]
    fn test_gaussian_pulse_peak_timing() {
        let mut gp = GaussianPulse::new(2.0, 100.0, 0.0);

        // Test at peak (t=0, which is tau)
        gp.update(0.0);
        assert!((gp.get_output(0) - 2.0).abs() < 1e-10);

        // Test with tau offset
        let mut gp2 = GaussianPulse::new(1.0, 1000.0, 5.0);
        gp2.update(5.0); // Peak should be at t=5
        assert!((gp2.get_output(0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_gaussian_pulse_amplitude() {
        let mut gp = GaussianPulse::new(10.0, 1000.0, 2.0);

        // At the center (tau=2), output should equal amplitude
        gp.update(2.0);
        assert!((gp.get_output(0) - 10.0).abs() < 1e-10);

        // Away from center, output should be less than amplitude
        gp.update(0.0);
        assert!(gp.get_output(0) < 10.0);

        gp.update(4.0);
        assert!(gp.get_output(0) < 10.0);
    }

    #[test]
    fn test_gaussian_pulse_bandwidth() {
        // Higher f_max means narrower pulse (smaller sigma)
        let mut gp_narrow = GaussianPulse::new(1.0, 1000.0, 0.0);
        let mut gp_wide = GaussianPulse::new(1.0, 100.0, 0.0);

        // At same time offset from center, narrow pulse decays faster
        let t_offset = 0.001;

        gp_narrow.update(t_offset);
        let val_narrow = gp_narrow.get_output(0);

        gp_wide.update(t_offset);
        let val_wide = gp_wide.get_output(0);

        // Narrow pulse (high f_max) should decay faster
        assert!(val_narrow < val_wide);
    }

    // Chirp tests - matching PathSim expectations

    #[test]
    fn test_chirp_init() {
        // Default-like
        let c = Chirp::new(1.0, 1.0, 1.0, 1.0, 0.0);
        assert_eq!(c.amplitude(), 1.0);
        assert_eq!(c.f0(), 1.0);
        assert_eq!(c.bandwidth(), 1.0);
        assert_eq!(c.sweep_time(), 1.0);
        assert_eq!(c.phase(), 0.0);

        // Specific
        let c2 = Chirp::new(2.0, 100.0, 200.0, 1.0, PI / 2.0);
        assert_eq!(c2.amplitude(), 2.0);
        assert_eq!(c2.f0(), 100.0);
        assert_eq!(c2.bandwidth(), 200.0);
        assert_eq!(c2.sweep_time(), 1.0);
        assert!((c2.phase() - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_chirp_frequency_sweep() {
        let c = Chirp::new(1.0, 10.0, 20.0, 2.0, 0.0);

        // Triangle wave starts at +1, so at t=0, frequency is f0 + BW
        let f_start = c.instantaneous_frequency(0.0);
        assert!((f_start - 30.0).abs() < 1.0); // f0 + BW = 10 + 20 = 30

        // At t=T/2 (t=1.0), triangle is at -1, so frequency is f0
        let f_mid = c.instantaneous_frequency(1.0);
        assert!((f_mid - 10.0).abs() < 1.0); // Should be f0 = 10

        // At end of sweep (t=2.0), triangle is back at +1
        let f_end = c.instantaneous_frequency(2.0);
        assert!((f_end - 30.0).abs() < 1.0); // Back to f0 + BW
    }

    #[test]
    fn test_chirp_start_frequency() {
        let mut c = Chirp::new(1.0, 50.0, 100.0, 1.0, 0.0);

        // At t=0, with phase=0, chirp output depends on integrated phase
        c.update(0.0);
        assert!((c.get_output(0) - 0.0).abs() < 1e-10);

        // At t=0, triangle=+1, so instantaneous frequency is f0 + BW
        let f_inst = c.instantaneous_frequency(0.0);
        assert!((f_inst - 150.0).abs() < 5.0); // Should be f0 + BW = 50 + 100 = 150
    }

    #[test]
    fn test_chirp_end_frequency() {
        let c = Chirp::new(1.0, 100.0, 200.0, 1.0, 0.0);

        // At t=0, triangle=+1, so frequency is f0 + BW = 300
        let f_start = c.instantaneous_frequency(0.0);
        assert!((f_start - 300.0).abs() < 5.0);

        // At t=0.5 (T/2), triangle=-1, so frequency is f0 = 100
        let f_mid = c.instantaneous_frequency(0.5);
        assert!((f_mid - 100.0).abs() < 5.0);

        // Verify frequency sweeps from high to low in first half
        assert!(f_start > f_mid); // Frequency decreases in first half period
    }
}
