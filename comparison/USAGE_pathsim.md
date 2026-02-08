# PathSim Benchmark - Usage Guide

## Quick Start

### Run the Full Benchmark Suite

```bash
cd /home/yongkyunshin/personal/RustSim/comparison
python3 pathsim_benchmark.py
```

This will:
- Run 4 different benchmark problems
- Test each problem with multiple solvers (12 solvers total)
- Save results to `pathsim_benchmark_results.json`
- Display timing comparisons and rankings
- Take approximately 10-30 minutes to complete

### Run Quick Tests

To verify the benchmark works without running the full suite:

```bash
cd /home/yongkyunshin/personal/RustSim/comparison
python3 test_pathsim_benchmark.py
```

This runs a smaller version with reduced duration for quick validation.

## Understanding the Output

### Console Output Example

```
================================================================================
BENCHMARK 1: Robertson's Stiff Chemical Kinetics (duration=100)
================================================================================
  Running Robertson with RK4... 45.234s (123456 steps)
  Running Robertson with GEAR52A... 2.345s (567 steps)
  ...

BENCHMARK SUMMARY
================================================================================

ROBERTSON:
--------------------------------------------------------------------------------
Solver          Time (s)     Steps      Speed Rank
--------------------------------------------------------------------------------
GEAR52A            2.345       567  #1
ESDIRK43           3.456       678  #2
RK4               45.234    123456  #12
```

### JSON Results File

Results are saved to `pathsim_benchmark_results.json`:

```json
{
  "timestamp": "2024-02-05T12:00:00",
  "benchmarks": {
    "robertson": [
      {
        "solver": "GEAR52A",
        "elapsed_time": 2.345,
        "num_steps": 567,
        "final_state": {
          "x": 0.617,
          "y": 1.6e-05,
          "z": 0.383,
          "time": 100.0
        }
      },
      ...
    ]
  }
}
```

## Customizing the Benchmark

### Modify Duration

Edit `pathsim_benchmark.py` and change the duration parameters in `run_all_benchmarks()`:

```python
# Original
result = benchmark_robertson(solver_class, solver_name, duration=100)

# Shorter for testing
result = benchmark_robertson(solver_class, solver_name, duration=10)
```

### Test Specific Solvers

Modify the `solvers` list in `run_all_benchmarks()`:

```python
# Original - all solvers
solvers = [
    (RK4, "RK4"),
    (RKBS32, "RKBS32"),
    # ... many more
]

# Custom - only test a few
solvers = [
    (RK4, "RK4"),
    (GEAR52A, "GEAR52A"),
    (RKDP54, "RKDP54"),
]
```

### Run Individual Benchmarks

You can import and run specific benchmarks programmatically:

```python
import sys
sys.path.insert(0, "/home/yongkyunshin/personal/pathsim/src")

from pathsim_benchmark import benchmark_robertson
from pathsim.solvers import GEAR52A

result = benchmark_robertson(GEAR52A, "GEAR52A", duration=100)
print(f"Elapsed: {result['elapsed_time']:.3f}s")
print(f"Steps: {result['num_steps']}")
print(f"Final state: {result['final_state']}")
```

### Adjust System Size

For the coupled oscillator benchmark, you can change the number of oscillators:

```python
# Original - 20 oscillators (40 state variables)
result = benchmark_coupled_oscillators(solver, name, n_oscillators=20, duration=50)

# Larger system - 50 oscillators (100 state variables)
result = benchmark_coupled_oscillators(solver, name, n_oscillators=50, duration=50)

# Smaller system - 10 oscillators (20 state variables)
result = benchmark_coupled_oscillators(solver, name, n_oscillators=10, duration=50)
```

### Adjust Stiffness

For Van der Pol, you can change the stiffness parameter mu:

```python
# Original - very stiff
result = benchmark_vanderpol(solver, name, mu=1000, duration=3000)

# Moderately stiff
result = benchmark_vanderpol(solver, name, mu=100, duration=300)

# Mildly stiff
result = benchmark_vanderpol(solver, name, mu=10, duration=30)
```

## Performance Analysis

### Interpreting Results

**For Stiff Problems (Robertson, Van der Pol)**:
- Implicit methods (GEAR, ESDIRK, BDF) should be much faster
- Explicit methods (RK4, RKDP54) will take many more steps
- Higher-order implicit methods usually win for extreme stiffness

**For Non-Stiff Problems (Lorenz, Coupled Oscillators)**:
- Explicit methods (RKDP54, RKCK54) often competitive
- Adaptive methods typically beat fixed-step methods
- Implicit methods may be slower due to Newton iteration overhead

**Number of Steps**:
- Fewer steps generally means better efficiency
- Adaptive methods adjust step size for efficiency
- Stiff problems require smaller steps for explicit methods

### Comparing Against Other Libraries

To compare PathSim against scipy.integrate.solve_ivp or other libraries:

1. Use the same problem definitions (func, x0, duration)
2. Use similar tolerances (rtol=1e-4, atol=1e-6)
3. Measure wall-clock time with `time.perf_counter()`
4. Compare final state values to verify accuracy
5. Record number of function evaluations if available

Example comparison script structure:

```python
import time
import numpy as np
from scipy.integrate import solve_ivp

# Robertson problem
a, b, c = 0.04, 1e4, 3e7
x0 = np.array([1.0, 0.0, 0.0])

def func(t, x):
    return np.array([
        -a*x[0] + b*x[1]*x[2],
         a*x[0] - b*x[1]*x[2] - c*x[1]**2,
                                c*x[1]**2
    ])

# Scipy benchmark
start = time.perf_counter()
sol = solve_ivp(func, [0, 100], x0, method='BDF', rtol=1e-4, atol=1e-6)
scipy_time = time.perf_counter() - start

print(f"Scipy BDF: {scipy_time:.3f}s, {sol.nfev} evaluations")
print(f"Final state: {sol.y[:, -1]}")

# PathSim benchmark
from pathsim_benchmark import benchmark_robertson
from pathsim.solvers import BDF3

result = benchmark_robertson(BDF3, "BDF3", duration=100)
print(f"PathSim BDF3: {result['elapsed_time']:.3f}s, {result['num_steps']} steps")
print(f"Final state: {result['final_state']}")
```

## Troubleshooting

### Import Errors

If you get import errors, verify PathSim path:

```python
import sys
import os

pathsim_path = "/home/yongkyunshin/personal/pathsim/src"
if os.path.exists(pathsim_path):
    print(f"PathSim found at: {pathsim_path}")
    sys.path.insert(0, pathsim_path)
else:
    print(f"PathSim NOT found at: {pathsim_path}")
    print("Update the path in pathsim_benchmark.py")
```

### Solver Failures

Some solvers may fail on certain problems:
- Explicit methods often fail on very stiff problems (mu=1000)
- This is expected and logged in the results with an "error" field
- Failed runs are excluded from the summary rankings

### Memory Issues

For very large systems or long durations:
- Scope stores all data points in memory
- Consider reducing duration or sampling_period
- Or modify Scope to save data incrementally to disk

### Performance Issues

If benchmarks run too slowly:
- Reduce duration parameters
- Reduce n_oscillators for coupled system
- Reduce mu for Van der Pol
- Test fewer solvers
- Use fixed-step solvers instead of adaptive

## Advanced Usage

### Adding New Benchmarks

To add a new benchmark problem:

1. Create a new benchmark function following the pattern:

```python
def benchmark_my_problem(solver_class, solver_name, duration=100):
    """Your problem description."""
    print(f"  Running MyProblem with {solver_name}...", end="", flush=True)

    # Define problem
    x0 = np.array([...])
    def func(x, u, t):
        return ...

    # Build system
    ODE_block = ODE(func, x0)
    Sco = Scope(labels=[...])
    blocks = [ODE_block, Sco]
    connections = [...]

    # Create simulation
    Sim = Simulation(blocks, connections, Solver=solver_class, ...)

    # Time it
    start_time = time.perf_counter()
    Sim.run(duration)
    elapsed = time.perf_counter() - start_time

    # Extract results
    rec_time, rec_data = Sco.read()
    final_state = {...}

    print(f" {elapsed:.3f}s ({len(rec_time)} steps)")

    return {
        "solver": solver_name,
        "elapsed_time": elapsed,
        "num_steps": len(rec_time),
        "final_state": final_state
    }
```

2. Add it to `run_all_benchmarks()`:

```python
# Benchmark 5: My Problem
print("\n" + "="*80)
print("BENCHMARK 5: My Problem Description")
print("="*80)
results["benchmarks"]["my_problem"] = []

for solver_class, solver_name in solvers:
    try:
        result = benchmark_my_problem(solver_class, solver_name, duration=100)
        results["benchmarks"]["my_problem"].append(result)
    except Exception as e:
        print(f"  FAILED: {solver_name} - {str(e)}")
        results["benchmarks"]["my_problem"].append({
            "solver": solver_name,
            "error": str(e)
        })
```

### Profiling

To profile specific benchmarks with cProfile:

```python
import cProfile
import pstats
from pathsim_benchmark import benchmark_robertson
from pathsim.solvers import GEAR52A

profiler = cProfile.Profile()
profiler.enable()

result = benchmark_robertson(GEAR52A, "GEAR52A", duration=100)

profiler.disable()
stats = pstats.Stats(profiler)
stats.sort_stats('cumulative')
stats.print_stats(20)  # Top 20 functions
```

### Memory Profiling

To profile memory usage with memory_profiler:

```bash
pip install memory_profiler
python -m memory_profiler pathsim_benchmark.py
```

Or programmatically:

```python
from memory_profiler import profile

@profile
def run_benchmark():
    from pathsim_benchmark import benchmark_robertson
    from pathsim.solvers import GEAR52A
    return benchmark_robertson(GEAR52A, "GEAR52A", duration=100)

result = run_benchmark()
```
