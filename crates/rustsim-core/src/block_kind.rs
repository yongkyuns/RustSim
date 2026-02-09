//! BlockKind enum for efficient dispatch of heterogeneous block types
//!
//! This module provides a type-erased enum wrapper around all concrete block types,
//! enabling storage of different block types in a single collection while maintaining
//! compile-time type safety and efficient dispatch via match expressions.

use crate::block::{Block, StepResult};
use crate::blocks::*;

/// Macro to reduce boilerplate for Block trait method delegation
macro_rules! dispatch_method {
    ($self:ident, $method:ident, $($args:expr),*) => {
        match $self {
            BlockKind::Constant(b) => b.$method($($args),*),
            BlockKind::Amplifier(b) => b.$method($($args),*),
            BlockKind::Integrator(b) => b.$method($($args),*),
            BlockKind::Differentiator(b) => b.$method($($args),*),
            BlockKind::Adder2(b) => b.$method($($args),*),
            BlockKind::Adder3(b) => b.$method($($args),*),
            BlockKind::Adder4(b) => b.$method($($args),*),
            BlockKind::Multiplier2(b) => b.$method($($args),*),
            BlockKind::Multiplier3(b) => b.$method($($args),*),
            BlockKind::Multiplier4(b) => b.$method($($args),*),
            BlockKind::Scope1(b) => b.$method($($args),*),
            BlockKind::Scope2(b) => b.$method($($args),*),
            BlockKind::Scope4(b) => b.$method($($args),*),
            BlockKind::Sinusoidal(b) => b.$method($($args),*),
            BlockKind::Step(b) => b.$method($($args),*),
            BlockKind::Ramp(b) => b.$method($($args),*),
            BlockKind::Pulse(b) => b.$method($($args),*),
            BlockKind::SquareWave(b) => b.$method($($args),*),
            BlockKind::TriangleWave(b) => b.$method($($args),*),
            BlockKind::Chirp(b) => b.$method($($args),*),
            BlockKind::GaussianPulse(b) => b.$method($($args),*),
            BlockKind::Clock(b) => b.$method($($args),*),
            BlockKind::PID(b) => b.$method($($args),*),
            BlockKind::AntiWindupPID(b) => b.$method($($args),*),
            BlockKind::Saturation(b) => b.$method($($args),*),
            BlockKind::RateLimiter(b) => b.$method($($args),*),
            BlockKind::LowpassRC(b) => b.$method($($args),*),
            BlockKind::HighpassRC(b) => b.$method($($args),*),
            BlockKind::Comparator(b) => b.$method($($args),*),
            BlockKind::Relay(b) => b.$method($($args),*),
            BlockKind::Switch(b) => b.$method($($args),*),
            BlockKind::Abs(b) => b.$method($($args),*),
            BlockKind::Sign(b) => b.$method($($args),*),
            BlockKind::Sqrt(b) => b.$method($($args),*),
            BlockKind::Pow(b) => b.$method($($args),*),
            BlockKind::Exp(b) => b.$method($($args),*),
            BlockKind::Log(b) => b.$method($($args),*),
            BlockKind::Sin(b) => b.$method($($args),*),
            BlockKind::Cos(b) => b.$method($($args),*),
            BlockKind::Tan(b) => b.$method($($args),*),
            BlockKind::Min2(b) => b.$method($($args),*),
            BlockKind::Max2(b) => b.$method($($args),*),
            #[cfg(feature = "rand-support")]
            BlockKind::WhiteNoise(b) => b.$method($($args),*),
            #[cfg(feature = "rand-support")]
            BlockKind::UniformNoise(b) => b.$method($($args),*),
            #[cfg(feature = "rand-support")]
            BlockKind::PinkNoise(b) => b.$method($($args),*),
            BlockKind::Counter(b) => b.$method($($args),*),
            BlockKind::CounterUp(b) => b.$method($($args),*),
            BlockKind::CounterDown(b) => b.$method($($args),*),
            BlockKind::SampleHold(b) => b.$method($($args),*),
            BlockKind::LUT1D(b) => b.$method($($args),*),
        }
    };
}

/// Type-erased enum wrapping all block types
///
/// This enum enables heterogeneous collections of blocks while maintaining
/// efficient dispatch via match expressions rather than dynamic dispatch.
///
/// Note: Some blocks with const generics (like Delay, FIR, etc.) are included
/// with common default sizes. For custom sizes, use the specific block types directly.
///
/// # Example
///
/// ```ignore
/// use rustsim_core::block_kind::BlockKind;
/// use rustsim_core::blocks::*;
///
/// let blocks: Vec<BlockKind> = vec![
///     BlockKind::from(Constant::new(1.0)),
///     BlockKind::from(Amplifier::new(2.0)),
///     BlockKind::from(Integrator::new(0.0)),
/// ];
///
/// for block in &mut blocks {
///     block.update(0.0);
/// }
/// ```
#[derive(Clone)]
pub enum BlockKind {
    // Sources
    Constant(Constant),
    Sinusoidal(Sinusoidal),
    Step(Step),
    Ramp(Ramp),
    Pulse(Pulse),
    SquareWave(SquareWave),
    TriangleWave(TriangleWave),
    Chirp(Chirp),
    GaussianPulse(GaussianPulse),
    Clock(Clock),

    // Basic operations
    Amplifier(Amplifier),
    Adder2(Adder<2>),
    Adder3(Adder<3>),
    Adder4(Adder<4>),
    Multiplier2(Multiplier<2>),
    Multiplier3(Multiplier<3>),
    Multiplier4(Multiplier<4>),

    // Dynamic blocks
    Integrator(Integrator),
    Differentiator(Differentiator),

    // Control systems
    PID(PID),
    AntiWindupPID(AntiWindupPID),

    // Nonlinear
    Saturation(Saturation),
    RateLimiter(RateLimiter),

    // Filters
    LowpassRC(LowpassRC),
    HighpassRC(HighpassRC),

    // Logic
    Comparator(Comparator),
    Relay(Relay),
    Switch(Switch),

    // Math operations
    Abs(Abs),
    Sign(Sign),
    Sqrt(Sqrt),
    Pow(Pow),
    Exp(Exp),
    Log(Log),
    Sin(Sin),
    Cos(Cos),
    Tan(Tan),
    Min2(Min<2>),
    Max2(Max<2>),

    // Noise sources (feature-gated)
    #[cfg(feature = "rand-support")]
    WhiteNoise(WhiteNoise),
    #[cfg(feature = "rand-support")]
    UniformNoise(UniformNoise),
    #[cfg(feature = "rand-support")]
    PinkNoise(PinkNoise),

    // Data acquisition
    Scope1(Scope<1, 1000>),
    Scope2(Scope<2, 1000>),
    Scope4(Scope<4, 1000>),
    SampleHold(SampleHold),

    // Counters
    Counter(Counter),
    CounterUp(CounterUp),
    CounterDown(CounterDown),

    // Lookup tables
    LUT1D(LUT1D),
}

impl Block for BlockKind {
    const NUM_INPUTS: usize = 0; // Runtime-determined
    const NUM_OUTPUTS: usize = 0; // Runtime-determined
    const IS_DYNAMIC: bool = false; // Runtime-determined

    fn num_inputs(&self) -> usize {
        dispatch_method!(self, num_inputs,)
    }

    fn num_outputs(&self) -> usize {
        dispatch_method!(self, num_outputs,)
    }

    fn inputs(&self) -> &[f64] {
        dispatch_method!(self, inputs,)
    }

    fn inputs_mut(&mut self) -> &mut [f64] {
        dispatch_method!(self, inputs_mut,)
    }

    fn outputs(&self) -> &[f64] {
        dispatch_method!(self, outputs,)
    }

    fn outputs_mut(&mut self) -> &mut [f64] {
        dispatch_method!(self, outputs_mut,)
    }

    fn update(&mut self, t: f64) {
        dispatch_method!(self, update, t)
    }

    fn step(&mut self, t: f64, dt: f64) -> StepResult {
        dispatch_method!(self, step, t, dt)
    }

    fn buffer(&mut self) {
        dispatch_method!(self, buffer,)
    }

    fn revert(&mut self) {
        dispatch_method!(self, revert,)
    }

    fn reset(&mut self) {
        dispatch_method!(self, reset,)
    }

    fn is_dynamic(&self) -> bool {
        dispatch_method!(self, is_dynamic,)
    }
}

// From implementations for automatic conversion

impl From<Constant> for BlockKind {
    fn from(block: Constant) -> Self {
        BlockKind::Constant(block)
    }
}

impl From<Amplifier> for BlockKind {
    fn from(block: Amplifier) -> Self {
        BlockKind::Amplifier(block)
    }
}

impl From<Integrator> for BlockKind {
    fn from(block: Integrator) -> Self {
        BlockKind::Integrator(block)
    }
}

impl From<Differentiator> for BlockKind {
    fn from(block: Differentiator) -> Self {
        BlockKind::Differentiator(block)
    }
}

impl From<Adder<2>> for BlockKind {
    fn from(block: Adder<2>) -> Self {
        BlockKind::Adder2(block)
    }
}

impl From<Adder<3>> for BlockKind {
    fn from(block: Adder<3>) -> Self {
        BlockKind::Adder3(block)
    }
}

impl From<Adder<4>> for BlockKind {
    fn from(block: Adder<4>) -> Self {
        BlockKind::Adder4(block)
    }
}

impl From<Multiplier<2>> for BlockKind {
    fn from(block: Multiplier<2>) -> Self {
        BlockKind::Multiplier2(block)
    }
}

impl From<Multiplier<3>> for BlockKind {
    fn from(block: Multiplier<3>) -> Self {
        BlockKind::Multiplier3(block)
    }
}

impl From<Multiplier<4>> for BlockKind {
    fn from(block: Multiplier<4>) -> Self {
        BlockKind::Multiplier4(block)
    }
}

impl From<Scope<1, 1000>> for BlockKind {
    fn from(block: Scope<1, 1000>) -> Self {
        BlockKind::Scope1(block)
    }
}

impl From<Scope<2, 1000>> for BlockKind {
    fn from(block: Scope<2, 1000>) -> Self {
        BlockKind::Scope2(block)
    }
}

impl From<Scope<4, 1000>> for BlockKind {
    fn from(block: Scope<4, 1000>) -> Self {
        BlockKind::Scope4(block)
    }
}

impl From<Sinusoidal> for BlockKind {
    fn from(block: Sinusoidal) -> Self {
        BlockKind::Sinusoidal(block)
    }
}

impl From<Step> for BlockKind {
    fn from(block: Step) -> Self {
        BlockKind::Step(block)
    }
}

impl From<Ramp> for BlockKind {
    fn from(block: Ramp) -> Self {
        BlockKind::Ramp(block)
    }
}

impl From<Pulse> for BlockKind {
    fn from(block: Pulse) -> Self {
        BlockKind::Pulse(block)
    }
}

impl From<SquareWave> for BlockKind {
    fn from(block: SquareWave) -> Self {
        BlockKind::SquareWave(block)
    }
}

impl From<TriangleWave> for BlockKind {
    fn from(block: TriangleWave) -> Self {
        BlockKind::TriangleWave(block)
    }
}

impl From<Chirp> for BlockKind {
    fn from(block: Chirp) -> Self {
        BlockKind::Chirp(block)
    }
}

impl From<GaussianPulse> for BlockKind {
    fn from(block: GaussianPulse) -> Self {
        BlockKind::GaussianPulse(block)
    }
}

impl From<Clock> for BlockKind {
    fn from(block: Clock) -> Self {
        BlockKind::Clock(block)
    }
}

impl From<PID> for BlockKind {
    fn from(block: PID) -> Self {
        BlockKind::PID(block)
    }
}

impl From<AntiWindupPID> for BlockKind {
    fn from(block: AntiWindupPID) -> Self {
        BlockKind::AntiWindupPID(block)
    }
}

impl From<Saturation> for BlockKind {
    fn from(block: Saturation) -> Self {
        BlockKind::Saturation(block)
    }
}

impl From<RateLimiter> for BlockKind {
    fn from(block: RateLimiter) -> Self {
        BlockKind::RateLimiter(block)
    }
}

impl From<LowpassRC> for BlockKind {
    fn from(block: LowpassRC) -> Self {
        BlockKind::LowpassRC(block)
    }
}

impl From<HighpassRC> for BlockKind {
    fn from(block: HighpassRC) -> Self {
        BlockKind::HighpassRC(block)
    }
}

impl From<Comparator> for BlockKind {
    fn from(block: Comparator) -> Self {
        BlockKind::Comparator(block)
    }
}

impl From<Relay> for BlockKind {
    fn from(block: Relay) -> Self {
        BlockKind::Relay(block)
    }
}

impl From<Switch> for BlockKind {
    fn from(block: Switch) -> Self {
        BlockKind::Switch(block)
    }
}

impl From<Abs> for BlockKind {
    fn from(block: Abs) -> Self {
        BlockKind::Abs(block)
    }
}

impl From<Sign> for BlockKind {
    fn from(block: Sign) -> Self {
        BlockKind::Sign(block)
    }
}

impl From<Sqrt> for BlockKind {
    fn from(block: Sqrt) -> Self {
        BlockKind::Sqrt(block)
    }
}

impl From<Pow> for BlockKind {
    fn from(block: Pow) -> Self {
        BlockKind::Pow(block)
    }
}

impl From<Exp> for BlockKind {
    fn from(block: Exp) -> Self {
        BlockKind::Exp(block)
    }
}

impl From<Log> for BlockKind {
    fn from(block: Log) -> Self {
        BlockKind::Log(block)
    }
}

impl From<Sin> for BlockKind {
    fn from(block: Sin) -> Self {
        BlockKind::Sin(block)
    }
}

impl From<Cos> for BlockKind {
    fn from(block: Cos) -> Self {
        BlockKind::Cos(block)
    }
}

impl From<Tan> for BlockKind {
    fn from(block: Tan) -> Self {
        BlockKind::Tan(block)
    }
}

impl From<Min<2>> for BlockKind {
    fn from(block: Min<2>) -> Self {
        BlockKind::Min2(block)
    }
}

impl From<Max<2>> for BlockKind {
    fn from(block: Max<2>) -> Self {
        BlockKind::Max2(block)
    }
}

#[cfg(feature = "rand-support")]
impl From<WhiteNoise> for BlockKind {
    fn from(block: WhiteNoise) -> Self {
        BlockKind::WhiteNoise(block)
    }
}

#[cfg(feature = "rand-support")]
impl From<UniformNoise> for BlockKind {
    fn from(block: UniformNoise) -> Self {
        BlockKind::UniformNoise(block)
    }
}

#[cfg(feature = "rand-support")]
impl From<PinkNoise> for BlockKind {
    fn from(block: PinkNoise) -> Self {
        BlockKind::PinkNoise(block)
    }
}

impl From<Counter> for BlockKind {
    fn from(block: Counter) -> Self {
        BlockKind::Counter(block)
    }
}

impl From<CounterUp> for BlockKind {
    fn from(block: CounterUp) -> Self {
        BlockKind::CounterUp(block)
    }
}

impl From<CounterDown> for BlockKind {
    fn from(block: CounterDown) -> Self {
        BlockKind::CounterDown(block)
    }
}

impl From<SampleHold> for BlockKind {
    fn from(block: SampleHold) -> Self {
        BlockKind::SampleHold(block)
    }
}

impl From<LUT1D> for BlockKind {
    fn from(block: LUT1D) -> Self {
        BlockKind::LUT1D(block)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_block_kind_from_conversion() {
        let constant = Constant::new(5.0);
        let block: BlockKind = constant.into();
        assert_eq!(block.num_outputs(), 1);
    }

    #[test]
    fn test_block_kind_heterogeneous_collection() {
        let blocks: Vec<BlockKind> = vec![
            Constant::new(1.0).into(),
            Amplifier::new(2.0).into(),
            Integrator::new(0.0).into(),
        ];

        assert_eq!(blocks.len(), 3);
        assert!(!blocks[0].is_dynamic());
        assert!(!blocks[1].is_dynamic());
        assert!(blocks[2].is_dynamic());
    }

    #[test]
    fn test_block_kind_update() {
        let mut block: BlockKind = Amplifier::new(3.0).into();
        block.set_input(0, 2.0);
        block.update(0.0);
        assert_eq!(block.get_output(0), 6.0);
    }

    #[test]
    fn test_block_kind_reset() {
        let mut block: BlockKind = Integrator::new(5.0).into();
        block.set_input(0, 1.0);
        block.update(0.0);
        block.step(0.0, 0.1);
        block.reset();
        assert_eq!(block.get_output(0), 5.0);
    }

    #[test]
    fn test_block_kind_adder_variants() {
        let adder2: BlockKind = Adder::<2>::new().into();
        let adder3: BlockKind = Adder::<3>::new().into();
        let adder4: BlockKind = Adder::<4>::new().into();

        assert_eq!(adder2.num_inputs(), 2);
        assert_eq!(adder3.num_inputs(), 3);
        assert_eq!(adder4.num_inputs(), 4);
    }

    #[test]
    fn test_block_kind_scope_variants() {
        let scope1: BlockKind = Scope::<1, 1000>::new().into();
        let scope2: BlockKind = Scope::<2, 1000>::new().into();
        let scope4: BlockKind = Scope::<4, 1000>::new().into();

        assert_eq!(scope1.num_inputs(), 1);
        assert_eq!(scope2.num_inputs(), 2);
        assert_eq!(scope4.num_inputs(), 4);
    }
}
