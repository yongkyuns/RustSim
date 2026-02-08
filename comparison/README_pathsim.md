# PathSim Benchmark Suite

This comprehensive benchmark suite tests PathSim's performance on computationally challenging dynamical systems.

## Benchmark Script

**File**: `pathsim_benchmark.py`

## Overview

The benchmark suite evaluates PathSim's solvers on four different types of challenging systems:

### 1. Robertson's Stiff Chemical Kinetics
- **Type**: Very stiff ODE system
- **Characteristics**: Time scales vary over 10+ orders of magnitude
- **Duration**: 100 time units
- **State Dimension**: 3
- **Description**: Classic test case for stiff ODE solvers, modeling a chemical reaction with species concentrations:
  - `dy1/dt = -0.04*y1 + 1e4*y2*y3`
  - `dy2/dt = 0.04*y1 - 1e4*y2*y3 - 3e7*y2^2`
  - `dy3/dt = 3e7*y2^2`

### 2. Van der Pol Oscillator (Large mu)
- **Type**: Very stiff nonlinear oscillator
- **Characteristics**: Sharp transitions, limit cycle behavior
- **Parameters**: mu = 1000 (extreme stiffness)
- **Duration**: 3000 time units
- **State Dimension**: 2
- **Description**: Nonlinear oscillator with:
  - `dx1/dt = x2`
  - `dx2/dt = mu*(1 - x1^2)*x2 - x1`
- **Note**: Only implicit solvers tested (too stiff for explicit methods)

### 3. Lorenz Attractor
- **Type**: Chaotic dynamical system
- **Characteristics**: Sensitive to initial conditions, strange attractor
- **Parameters**: sigma=10, rho=28, beta=8/3
- **Duration**: 200 time units
- **State Dimension**: 3
- **Description**: Famous chaotic system demonstrating butterfly effect
- **Requirements**: High precision needed for accurate long-term evolution

### 4. Coupled Oscillator Chain
- **Type**: Large coupled linear system
- **Characteristics**: 20 coupled harmonic oscillators
- **Duration**: 50 time units
- **State Dimension**: 40 (20 positions + 20 velocities)
- **Description**: Chain of spring-mass-damper systems coupled to neighbors
- **Tests**: Performance on larger dimensional systems

## Solvers Tested

### Explicit Runge-Kutta Methods
- **RK4**: Classic 4th-order Runge-Kutta
- **RKBS32**: Bogacki-Shampine 3(2) adaptive
- **RKCK54**: Cash-Karp 5(4) adaptive
- **RKDP54**: Dormand-Prince 5(4) adaptive

### Implicit ESDIRK Methods
- **ESDIRK32**: 3rd-order L-stable ESDIRK
- **ESDIRK43**: 4th-order L-stable ESDIRK
- **ESDIRK54**: 5th-order L-stable ESDIRK

### Implicit Multistep Methods
- **GEAR21**: 2nd-order GEAR (Adams)
- **GEAR32**: 3rd-order GEAR
- **GEAR43**: 4th-order GEAR
- **GEAR52A**: 5th-order GEAR with adaptive stepping
- **BDF3**: 3rd-order Backward Differentiation Formula

## Running the Benchmark

### Prerequisites
```bash
# Ensure PathSim is installed or available at:
# /home/yongkyunshin/personal/pathsim/src

# The script automatically adds PathSim to the Python path
```

### Execution
```bash
# Run from the benchmarks directory
cd /home/yongkyunshin/personal/RustSim/comparison
python3 pathsim_benchmark.py

# Or run directly if executable
./pathsim_benchmark.py
```

### Expected Runtime
- **Total duration**: 10-30 minutes (depending on system)
- Individual benchmark times vary by solver and problem stiffness

## Output

### Console Output
- Real-time progress for each benchmark and solver
- Elapsed time and number of steps for each run
- Summary table ranking solvers by performance for each benchmark

### JSON Results File
**Location**: `pathsim_benchmark_results.json`

**Structure**:
```json
{
  "timestamp": "ISO timestamp",
  "benchmarks": {
    "robertson": [
      {
        "solver": "GEAR52A",
        "elapsed_time": 1.234,
        "num_steps": 5678,
        "final_state": {
          "x": 0.123,
          "y": 0.456,
          "z": 0.789,
          "time": 100.0
        }
      },
      ...
    ],
    "vanderpol": [...],
    "lorenz": [...],
    "coupled_oscillators": [...]
  }
}
```

## Performance Metrics

For each benchmark run, the following metrics are recorded:

1. **elapsed_time**: Wall-clock time in seconds
2. **num_steps**: Total number of simulation steps completed
3. **final_state**: Final values of state variables (for accuracy verification)
4. **error**: Error message if the solver failed

## Expected Results

### Typical Performance Characteristics

**Robertson (Stiff)**:
- Implicit solvers (GEAR, ESDIRK) should significantly outperform explicit methods
- GEAR52A typically fastest for this extreme stiffness

**Van der Pol (Stiff)**:
- Only implicit solvers can handle mu=1000 efficiently
- Higher-order GEAR methods usually best

**Lorenz (Chaotic)**:
- Both explicit and implicit methods viable
- Adaptive methods (RKDP54, GEAR52A) often most efficient
- Higher-order methods needed for accuracy

**Coupled Oscillators (Large)**:
- Performance depends on system size and coupling strength
- Implicit methods may have advantage due to coupling stiffness
- Fixed-step methods may be competitive for moderate coupling

## Use Cases

This benchmark suite is useful for:

1. **Solver Selection**: Identifying the best solver for different problem types
2. **Performance Comparison**: Comparing PathSim against other ODE solver libraries
3. **Regression Testing**: Ensuring performance doesn't degrade with code changes
4. **Hardware Benchmarking**: Comparing performance across different machines
5. **Optimization Validation**: Verifying that code optimizations improve performance

## Customization

You can modify the benchmark parameters by editing the script:

- **duration**: Simulation duration for each benchmark
- **n_oscillators**: Number of oscillators in coupled system (default: 20)
- **mu**: Stiffness parameter for Van der Pol (default: 1000)
- **tolerance_lte_abs/rel**: Error tolerances for adaptive solvers
- **solvers**: List of solvers to test

## Accuracy Verification

The final states can be compared against reference solutions:

**Robertson at t=100**:
- y1 ≈ 0.617 (approximate, varies slightly with solver)
- y2 ≈ 1.6e-5
- y3 ≈ 0.383

**Lorenz** (chaotic, no exact solution):
- Verify state remains bounded
- Check for numerical divergence

## Notes

- Some solvers may fail on certain problems (e.g., explicit methods on very stiff systems)
- Failed runs are logged with error messages in the results file
- The benchmark uses Jacobian functions where available to improve implicit solver performance
- Memory usage is not tracked but can be monitored separately with tools like `memory_profiler`

## Future Enhancements

Potential additions to the benchmark suite:

1. Memory usage profiling
2. Accuracy metrics (compare against reference solutions)
3. Adaptive vs. fixed timestep comparison
4. Problem scaling studies (vary system size)
5. Additional problem types (DAEs, PDEs discretized as ODEs)
6. Parallel solver benchmarking (if available)
