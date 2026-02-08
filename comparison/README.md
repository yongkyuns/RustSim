# RustSim Benchmarks

This directory contains benchmark utilities for comparing RustSim (Rust) performance against PathSim-style (Python) implementations.

## Quick Start

Run the complete benchmark suite:

```bash
./run_benchmarks.sh
```

This will:
1. Compile RustSim benchmarks in release mode with full optimizations
2. Run both Python and Rust implementations of the harmonic oscillator
3. Compare execution times and calculate speedup ratios
4. Verify numerical accuracy between implementations
5. Generate a detailed summary report

## Prerequisites

- **Rust**: Version 1.70 or later (install from [rustup.rs](https://rustup.rs))
- **Python**: Version 3.8 or later
- **jq** (optional): For enhanced JSON parsing in the comparison output

### Checking Dependencies

The benchmark script will automatically check for required dependencies and provide helpful error messages if anything is missing.

## Available Benchmarks

### Harmonic Oscillator

**System**: d²x/dt² = -x

A classic test problem with known analytical solution: x(t) = cos(t) for initial conditions x(0) = 1, v(0) = 0.

- **Python implementation**: `harmonic_oscillator_python.py`
- **Rust implementation**: `harmonic_oscillator_rust.rs`
- **Duration**: 2π seconds (one complete period)
- **Time step**: 0.001 seconds
- **Total steps**: ~6,283

This benchmark tests:
- Raw computational performance
- Numerical integration accuracy
- Memory efficiency
- Energy conservation

### Coupled Oscillators (NEW)

**System**: m*x_i'' = -k*(x_i - x_{i-1}) - k*(x_i - x_{i+1}) - c*x_i'

A chain of N=20 coupled harmonic oscillators connected by springs with fixed boundary conditions.

- **PathSim implementation**: `coupled_oscillators_pathsim.py`
- **RustSim implementation**: `coupled_oscillators_rustsim.rs` (also in `examples/`)
- **State dimension**: 40 (20 positions + 20 velocities)
- **Duration**: 100 seconds
- **Time step**: 0.001 seconds
- **Total steps**: 100,000
- **Solver**: RK4 (4th order Runge-Kutta)

**Parameters**:
- Number of oscillators: N = 20
- Mass: m = 1.0
- Spring constant: k = 100.0
- Damping coefficient: c = 0.1

**Running the benchmark**:

RustSim (Rust):
```bash
cargo build --release --example coupled_oscillators_rustsim
./target/release/examples/coupled_oscillators_rustsim
```

PathSim (Python):
```bash
python3 comparison/coupled_oscillators_pathsim.py
```

**Typical Results**:

| Implementation | Wall-clock time | Steps/second | Time per step | Speedup |
|----------------|----------------|--------------|---------------|---------|
| **RustSim**    | ~0.03 seconds  | ~3,000,000   | ~0.3 μs       | 260x    |
| **PathSim**    | ~8.7 seconds   | ~11,500      | ~87 μs        | 1x      |

Both implementations produce identical numerical results:
- Initial energy: 21.233905
- Final energy: 0.000534
- Energy decreases due to damping (c = 0.1)

This benchmark tests:
- Large-scale ODE system performance (40 state variables)
- RK4 solver efficiency
- Vector operations and memory access patterns
- Comparison against PathSim (Python block diagram framework)

## Usage

### Run Complete Benchmark Suite

```bash
./run_benchmarks.sh
```

### Command-Line Options

```bash
./run_benchmarks.sh [OPTIONS]

Options:
  --skip-compile    Skip Rust compilation step (use existing binary)
  --python-only     Run only Python benchmark
  --rust-only       Run only Rust benchmark
  --keep-temp       Keep temporary result files
  --help, -h        Show help message
```

### Examples

Run only the Python benchmark:
```bash
./run_benchmarks.sh --python-only
```

Run only the Rust benchmark (assumes already compiled):
```bash
./run_benchmarks.sh --rust-only --skip-compile
```

Keep temporary JSON result files for further analysis:
```bash
./run_benchmarks.sh --keep-temp
```

### Running Individual Benchmarks

#### Python Benchmark

```bash
python3 harmonic_oscillator_python.py
```

With JSON output:
```bash
python3 harmonic_oscillator_python.py --json
```

#### Rust Benchmark

First, compile in release mode:
```bash
cd ..
cargo build --release --bin harmonic_oscillator_rust
```

Then run:
```bash
./target/release/harmonic_oscillator_rust
```

With JSON output:
```bash
./target/release/harmonic_oscillator_rust --json
```

## Understanding the Results

### Performance Metrics

- **Execution Time**: Total wall-clock time to complete the simulation
- **Steps/Second**: Number of integration steps performed per second (higher is better)
- **Speedup Ratio**: Python time ÷ Rust time (e.g., 50x means Rust is 50 times faster)

### Accuracy Metrics

Both implementations use the same Euler integration method, so accuracy should be comparable:

- **Position Error**: |computed - exact| for final position
- **Velocity Error**: |computed - exact| for final velocity
- **Energy Error**: Deviation from theoretical energy (0.5 for this system)

Small differences may occur due to:
- Floating-point arithmetic differences between languages
- Compiler optimizations
- Order of operations

### Expected Results

Typical performance characteristics:

- **Rust**: 0.001-0.003 seconds for ~6,283 steps
- **Python**: 0.05-0.15 seconds for ~6,283 steps
- **Speedup**: 30-100x (depends on system and Python implementation)

## Output Files

- `benchmark_report.txt`: Complete report with all results
- `/tmp/python_benchmark_results.json`: Python results in JSON format
- `/tmp/rust_benchmark_results.json`: Rust results in JSON format
- `/tmp/python_output.txt`: Python console output
- `/tmp/rust_output.txt`: Rust console output

Note: Temporary files are automatically cleaned up unless `--keep-temp` is specified.

## Benchmark Implementation Details

### Block-Based Architecture

Both implementations use the same block-based simulation architecture:

```
  ┌──────────────────────────────────┐
  │                                  │
  ▼                                  │
[gain: -1] ──► [velocity] ──► [position]
```

### Integration Method

Both use the Euler method for time integration:

```
x(t + dt) = x(t) + dx/dt * dt
```

This simple method was chosen for:
- Clear, comparable implementation
- Minimal complexity
- Focus on computational overhead rather than solver sophistication

### Optimization Settings

**Rust** (release mode):
- Link-Time Optimization (LTO): enabled
- Codegen units: 1 (maximum optimization)
- Optimization level: 3
- Target: native CPU features

**Python**:
- Standard CPython interpreter
- No NumPy or external acceleration libraries
- Pure Python implementation for fair comparison

## Adding New Benchmarks

To add a new benchmark:

1. **Create Python implementation**: `comparison/your_benchmark_python.py`
   - Follow the structure of `harmonic_oscillator_python.py`
   - Support `--json` flag for machine-readable output

2. **Create Rust implementation**: `comparison/your_benchmark_rust.rs`
   - Follow the structure of `harmonic_oscillator_rust.rs`
   - Support `--json` flag for machine-readable output

3. **Add binary to Cargo.toml**:
   ```toml
   [[bin]]
   name = "your_benchmark_rust"
   path = "comparison/your_benchmark_rust.rs"
   ```

4. **Update run_benchmarks.sh** to include your new benchmark

## Troubleshooting

### Rust Compilation Fails

If you see compilation errors:

1. Make sure you're using Rust 1.70 or later:
   ```bash
   rustc --version
   ```

2. Update Rust toolchain:
   ```bash
   rustup update
   ```

3. Clean and rebuild:
   ```bash
   cargo clean
   cargo build --release
   ```

### Python Script Fails

1. Check Python version (requires 3.8+):
   ```bash
   python3 --version
   ```

2. Make sure the script is executable:
   ```bash
   chmod +x harmonic_oscillator_python.py
   ```

### Large Performance Variations

If you see inconsistent results:

1. Close unnecessary applications
2. Ensure system is not under heavy load
3. Run benchmark multiple times and average results
4. Check CPU governor settings (use "performance" mode if available)

## Performance Analysis Tips

### Profiling Rust Code

Use `cargo-flamegraph` for detailed profiling:

```bash
cargo install flamegraph
cargo flamegraph --bin harmonic_oscillator_rust
```

### Profiling Python Code

Use `cProfile` for Python profiling:

```bash
python3 -m cProfile -s cumtime harmonic_oscillator_python.py
```

### Benchmarking with Different Problem Sizes

Modify the simulation parameters in the source files:

- **Duration**: Change from `2π` to larger values
- **Time step**: Smaller dt = more steps = longer runtime
- **System complexity**: Add more state variables

## Contributing

Contributions of new benchmarks are welcome! Please ensure:

1. Both Python and Rust implementations produce identical results
2. Benchmarks have clear documentation
3. Results are reproducible
4. Code follows project style guidelines

## License

MIT License - see main project LICENSE file for details.
