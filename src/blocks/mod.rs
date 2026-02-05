//! Block implementations for v2 architecture

mod adder;
mod amplifier;
mod butterworth;
mod delay;
mod differentiator;
mod filters;
mod fir;
mod function;
mod integrator;
mod logic;
mod math;
mod multiplier;
#[cfg(feature = "rand-support")]
mod noise;
mod ode;
mod pid;
mod scope;
mod sources;
// mod statespace;  // TODO: Fix type ambiguity in tests
mod converters;
mod counter;
mod kalman;
mod samplehold;
mod table;

pub use adder::Adder;
pub use amplifier::Amplifier;
pub use butterworth::{Allpass, ButterworthBandpass, ButterworthHighpass, ButterworthLowpass};
pub use delay::Delay;
pub use differentiator::Differentiator;
pub use filters::{HighpassRC, LowpassRC};
pub use fir::FIR;
pub use function::Function;
pub use integrator::Integrator;
pub use logic::{Comparator, Relay, Switch};
pub use math::*;
pub use multiplier::Multiplier;
#[cfg(feature = "rand-support")]
pub use noise::{PinkNoise, UniformNoise, WhiteNoise};
pub use ode::ODE;
pub use pid::{AntiWindupPID, RateLimiter, Saturation, PID};
pub use scope::Scope;
pub use sources::{
    Chirp, Clock, Constant, GaussianPulse, Pulse, Ramp, Sinusoidal, SquareWave, Step, TriangleWave,
};
// pub use statespace::StateSpace;  // TODO: Fix type ambiguity in tests
pub use converters::{ADC, DAC};
pub use counter::{Counter, CounterDown, CounterUp};
pub use kalman::KalmanFilter;
pub use samplehold::SampleHold;
pub use table::{LUT, LUT1D};
