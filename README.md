# RustSim

A high-performance, block-based time-domain simulation framework for dynamical systems, written in Rust.

**RustSim is a Rust port of [PathSim](https://github.com/pathsim/pathsim)**, bringing the power of Python's PathSim simulation framework to Rust with compile-time safety, zero-cost abstractions, and exceptional performance.

[![CI](https://github.com/your-username/rustsim/workflows/CI/badge.svg)](https://github.com/your-username/rustsim/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Features

- **Compile-time block sizing**: Fixed input/output sizes via const generics eliminate runtime overhead
- **Zero dynamic dispatch**: All block types are concrete at compile time
- **Comprehensive solver library**: 24+ ODE solvers from simple Euler to high-order adaptive methods
- **Event detection**: Zero-crossing detection with bisection refinement
- **WebAssembly support**: Build for browsers with `--target wasm32-unknown-unknown`
- **Cross-platform**: Runs on Linux, macOS, Windows, and WebAssembly

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
rustsim = "0.1"
```

## Quick Example

```rust
use rustsim::prelude::*;

fn main() {
    // Create blocks
    let mut source = Sinusoidal::new(1.0, 1.0, 0.0);  // amplitude, frequency, phase
    let mut integrator = Integrator::new(0.0);        // initial condition
    let mut scope = Scope::<1, 1000>::new();          // 1 channel, 1000 samples

    // Simulation loop
    let dt = 0.01;
    for i in 0..1000 {
        let t = i as f64 * dt;

        source.update(t);
        integrator.set_input(0, source.get_output(0));
        integrator.step(t, dt);
        scope.set_input(0, integrator.get_output(0));
        scope.update(t);
    }

    // Export results
    scope.to_csv("output.csv").unwrap();
}
```

## Implementation Status

### Blocks (50+ implementations)

| Category | Blocks | Status |
|----------|--------|--------|
| **Sources** | `Constant`, `Sinusoidal`, `Step`, `Ramp`, `SquareWave`, `TriangleWave`, `Pulse`, `Clock`, `GaussianPulse`, `Chirp` | Complete |
| **Math** | `Adder<N>`, `Multiplier<N>`, `Amplifier`, `Function`, `Sign`, `Pow`, `Clip`, `Mod`, `Norm<N>`, `Min<N>`, `Max<N>` | Complete |
| **Trigonometric** | `Sin`, `Cos`, `Tan`, `Asin`, `Acos`, `Atan`, `Atan2` | Complete |
| **Transcendental** | `Exp`, `Log`, `Log10`, `Sqrt`, `Abs` | Complete |
| **Dynamic** | `Integrator`, `Differentiator`, `ODE<N>`, `Delay<N>`, `StateSpace<N,M,P>` | Complete |
| **Filters** | `LowpassRC`, `HighpassRC`, `ButterworthLowpass`, `ButterworthHighpass`, `ButterworthBandpass`, `Allpass`, `FIR<N>` | Complete |
| **Control** | `PID`, `AntiWindupPID`, `RateLimiter`, `Saturation`, `KalmanFilter` | Complete |
| **Logic** | `Comparator`, `Switch`, `Relay` | Complete |
| **Noise** | `WhiteNoise`, `UniformNoise`, `PinkNoise` | Complete |
| **Utility** | `Scope<C,N>`, `SampleHold`, `Counter`, `CounterUp`, `CounterDown` | Complete |
| **Converters** | `ADC<N>`, `DAC<N>` | Complete |
| **Lookup Tables** | `LUT1D`, `LUT` | Complete |

### ODE Solvers (24 implementations)

| Family | Solvers | Order | Type |
|--------|---------|-------|------|
| **Basic** | `Euler` | 1 | Explicit |
| **Runge-Kutta** | `RK4` | 4 | Explicit, Fixed-step |
| **Embedded RK** | `RKF21`, `RKBS32`, `RKF45`, `RKCK54`, `RKDP54`, `RKV65`, `RKF78`, `RKDP87` | 2-8 | Explicit, Adaptive |
| **SSP RK** | `SSPRK22`, `SSPRK33`, `SSPRK34` | 2-3 | Explicit, SSP |
| **ESDIRK** | `ESDIRK32`, `ESDIRK4`, `ESDIRK43`, `ESDIRK54`, `ESDIRK85` | 3-8 | Implicit, L-stable |
| **DIRK** | `DIRK2`, `DIRK3`, `EulerBackward` | 1-3 | Implicit |
| **BDF** | `BDF2`, `BDF3`, `BDF4`, `BDF5`, `BDF6` | 2-6 | Implicit, Multistep |
| **GEAR** | `GEAR21`, `GEAR32`, `GEAR43`, `GEAR54`, `GEAR52A` | 2-5 | Implicit, Variable-order |

### Event System

| Component | Description | Status |
|-----------|-------------|--------|
| `ZeroCrossing` | Detect signal zero crossings | Complete |
| `ZeroCrossingUp` | Detect rising zero crossings | Complete |
| `ZeroCrossingDown` | Detect falling zero crossings | Complete |
| `Condition` | General condition-based events | Complete |
| `Schedule` | Time-scheduled events | Complete |
| `ScheduleList` | Multiple scheduled events | Complete |

### Optimization

| Component | Description | Status |
|-----------|-------------|--------|
| `Anderson` | Anderson acceleration for fixed-point iteration | Complete |
| `NewtonAnderson` | Newton-Anderson hybrid solver | Complete |

## Feature Flags

```toml
[features]
default = ["rand-support"]     # Random number generation for noise blocks
full = ["serialization", "spectrum", "parallel", "rand-support"]
parallel = ["rayon"]           # Parallel iteration support
serialization = ["serde"]      # Serialize/deserialize blocks
spectrum = ["rustfft"]         # FFT-based spectrum analysis
wasm = ["getrandom/js"]        # WebAssembly with JS random source
```

## WebAssembly

Build for web browsers:

```bash
# Core features (no random)
cargo build --lib --target wasm32-unknown-unknown --no-default-features

# Full features with JS random
cargo build --lib --target wasm32-unknown-unknown --features wasm
```

## Test Coverage

```
253 tests passing
6 tests ignored (solver edge cases under investigation)
```

## Comparison with PathSim

| Aspect | PathSim (Python) | RustSim |
|--------|------------------|---------|
| Language | Python 3.8+ | Rust 1.70+ |
| Block sizing | Runtime | Compile-time (const generics) |
| Type safety | Duck typing | Static, zero-cost |
| Performance | Interpreted + NumPy | Native, SIMD-ready |
| Memory | GC-managed | Stack-allocated where possible |
| WASM support | Limited (Pyodide) | Native |

## Architecture

RustSim uses a compile-time static architecture where block I/O sizes are fixed at compile time:

```rust
pub trait Block {
    const NUM_INPUTS: usize;
    const NUM_OUTPUTS: usize;
    const IS_DYNAMIC: bool;

    fn update(&mut self, t: f64);
    fn get_output(&self, index: usize) -> f64;
    fn set_input(&mut self, index: usize, value: f64);
}

pub trait DynamicBlock: Block {
    fn step(&mut self, t: f64, dt: f64) -> StepResult;
    fn buffer(&mut self);
    fn revert(&mut self);
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome! Please see the [PathSim repository](https://github.com/pathsim/pathsim) for reference implementations and test cases.

## Acknowledgments

RustSim is a port of [PathSim](https://github.com/pathsim/pathsim), a Python-based simulation framework. Thanks to the PathSim authors for the excellent reference implementation.
