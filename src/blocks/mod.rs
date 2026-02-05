//! Block implementations for v2 architecture

mod amplifier;
mod adder;
mod integrator;
mod sources;
mod function;
mod math;
#[cfg(feature = "rand-support")]
mod noise;
mod multiplier;
mod logic;
mod pid;
mod filters;
mod butterworth;
mod ode;
mod delay;
mod differentiator;
mod scope;
mod fir;
// mod statespace;  // TODO: Fix type ambiguity in tests
mod kalman;
mod samplehold;
mod counter;
mod converters;
mod table;

pub use amplifier::Amplifier;
pub use adder::Adder;
pub use integrator::Integrator;
pub use sources::{Constant, Sinusoidal, Step, Ramp, SquareWave, TriangleWave, Pulse, Clock, GaussianPulse, Chirp};
pub use function::Function;
pub use math::*;
#[cfg(feature = "rand-support")]
pub use noise::{WhiteNoise, UniformNoise, PinkNoise};
pub use multiplier::Multiplier;
pub use logic::{Comparator, Switch, Relay};
pub use pid::{PID, AntiWindupPID, RateLimiter, Saturation};
pub use filters::{LowpassRC, HighpassRC};
pub use butterworth::{ButterworthLowpass, ButterworthHighpass, ButterworthBandpass, Allpass};
pub use ode::ODE;
pub use delay::Delay;
pub use differentiator::Differentiator;
pub use scope::Scope;
pub use fir::FIR;
// pub use statespace::StateSpace;  // TODO: Fix type ambiguity in tests
pub use kalman::KalmanFilter;
pub use samplehold::SampleHold;
pub use counter::{Counter, CounterUp, CounterDown};
pub use converters::{ADC, DAC};
pub use table::{LUT1D, LUT};
