# PathSim Benchmark - Quick Reference Card

## Files Overview

| File | Purpose | Size |
|------|---------|------|
| `pathsim_benchmark.py` | Main benchmark script | 18 KB |
| `test_pathsim_benchmark.py` | Quick validation test | 2.2 KB |
| `README_pathsim.md` | Technical documentation | 6.6 KB |
| `USAGE_pathsim.md` | Usage guide | 9.4 KB |
| `PATHSIM_BENCHMARK_SUMMARY.md` | Overview & findings | 8.8 KB |
| `QUICK_REFERENCE_pathsim.md` | This file | - |

## One-Line Commands

```bash
# Quick test (1-2 min)
python3 /home/yongkyunshin/personal/RustSim/comparison/test_pathsim_benchmark.py

# Full benchmark (10-30 min)
python3 /home/yongkyunshin/personal/RustSim/comparison/pathsim_benchmark.py

# View results
cat /home/yongkyunshin/personal/RustSim/comparison/pathsim_benchmark_results.json
```

## Benchmark Problems

| Problem | Type | Dimension | Best Solver Type | Challenge |
|---------|------|-----------|------------------|-----------|
| Robertson | Stiff ODE | 3 | Implicit (GEAR) | Extreme stiffness |
| Van der Pol | Stiff ODE | 2 | Implicit (GEAR) | Sharp transitions |
| Lorenz | Chaotic | 3 | Explicit (RKDP54) | Long-term accuracy |
| Coupled Oscillators | Large System | 40 | Adaptive | Scalability |

## Solvers Tested

### Explicit (Good for non-stiff)
- RK4, RKBS32, RKCK54, RKDP54

### Implicit (Good for stiff)
- ESDIRK32, ESDIRK43, ESDIRK54
- GEAR21, GEAR32, GEAR43, GEAR52A
- BDF3

## Expected Performance

### Robertson (Stiff)
```
Winner: GEAR52A, ESDIRK43
Loser:  RK4, RKDP54
Ratio:  ~20x speedup for implicit
```

### Van der Pol (Very Stiff)
```
Winner: GEAR52A, GEAR43
Loser:  Explicit methods (fail)
Note:   Only implicit methods tested
```

### Lorenz (Non-Stiff)
```
Winner: RKDP54, RK4
Loser:  GEAR methods
Ratio:  ~10x speedup for explicit
```

### Coupled Oscillators
```
Winner: Adaptive methods
Varies: Depends on coupling strength
```

## Quick Customization

### Change Duration
```python
# In pathsim_benchmark.py, find:
result = benchmark_robertson(solver, name, duration=100)

# Change to:
result = benchmark_robertson(solver, name, duration=10)  # Faster
```

### Test Fewer Solvers
```python
# In run_all_benchmarks(), change solvers list:
solvers = [
    (RK4, "RK4"),           # Simple explicit
    (GEAR52A, "GEAR52A"),   # Best implicit
    (RKDP54, "RKDP54"),     # Best explicit
]
```

### Change System Size
```python
# For coupled oscillators:
result = benchmark_coupled_oscillators(solver, name,
    n_oscillators=10,   # Smaller (was 20)
    duration=50)
```

## Interpreting Results

### Timing
- Lower is better
- Compare similar problem types
- Hardware dependent

### Steps
- Fewer steps usually means more efficient
- Adaptive methods adjust automatically
- Very different step counts indicate solver mismatch

### Final State
- Use for accuracy verification
- Should be similar across solvers
- Large differences indicate error

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Import error | Check PathSim path in script |
| Solver fails | Expected for wrong solver type (explicit on stiff) |
| Too slow | Reduce duration or test fewer solvers |
| Out of memory | Reduce duration or system size |

## Common Use Cases

### Compare Solvers for Your Problem
1. Add your problem as new benchmark function
2. Run with all solvers
3. Check rankings for your problem type

### Performance Regression Testing
1. Run benchmark before code changes
2. Run benchmark after code changes
3. Compare timing results

### Hardware Benchmarking
1. Run on different machines
2. Compare wall-clock times
3. Normalize by fastest solver

### Research/Publication
1. Run full benchmark
2. Extract timing tables from JSON
3. Create comparison plots

## Performance Tips

### For Stiff Problems
- Use GEAR52A or ESDIRK43
- Provide Jacobian function
- Expect 10-50x speedup over explicit

### For Non-Stiff Problems
- Use RKDP54 or RKCK54
- Explicit methods faster
- Adaptive better than fixed-step

### For Large Systems
- Use adaptive methods
- Consider sparse Jacobians
- Monitor memory usage

### For High Accuracy
- Use higher-order methods (ESDIRK54, RKDP54)
- Tighten tolerances
- Verify with reference solution

## Key Metrics

```
elapsed_time:  Wall-clock seconds (lower = faster)
num_steps:     Simulation steps (fewer = more efficient)
final_state:   End values (for accuracy check)
error:         Failure message (if solver failed)
```

## Example Output

```
BENCHMARK 1: Robertson's Stiff Chemical Kinetics
  Running Robertson with GEAR52A... 0.022s (26 steps)
  Running Robertson with RK4... 0.375s (10002 steps)

Ranking:
  #1: GEAR52A     0.022s  (   26 steps)
  #2: RK4         0.375s  (10002 steps)
```

## Integration with Other Tools

### Profile with cProfile
```python
import cProfile
from pathsim_benchmark import benchmark_robertson
from pathsim.solvers import GEAR52A

cProfile.run('benchmark_robertson(GEAR52A, "GEAR52A", 100)')
```

### Memory Profile
```bash
pip install memory_profiler
python -m memory_profiler pathsim_benchmark.py
```

### Compare with SciPy
```python
from scipy.integrate import solve_ivp
import time

start = time.perf_counter()
sol = solve_ivp(func, [0, 100], x0, method='BDF')
scipy_time = time.perf_counter() - start

# Compare with PathSim result
```

## JSON Results Structure

```json
{
  "timestamp": "2024-02-05T12:00:00",
  "benchmarks": {
    "robertson": [
      {
        "solver": "GEAR52A",
        "elapsed_time": 2.345,
        "num_steps": 567,
        "final_state": {"x": 0.617, "y": 1.6e-05, "z": 0.383, "time": 100.0}
      }
    ]
  }
}
```

## Documentation Locations

- **Technical Details**: README_pathsim.md
- **How-To Guide**: USAGE_pathsim.md
- **Overview**: PATHSIM_BENCHMARK_SUMMARY.md
- **This Card**: QUICK_REFERENCE_pathsim.md

## Support

1. Check documentation files
2. Verify PathSim installation
3. Review error messages
4. Test with quick validation script

## Version Info

- **Created**: 2026-02-05
- **PathSim**: /home/yongkyunshin/personal/pathsim/src
- **Python**: 3.x required
- **Dependencies**: NumPy, PathSim

---

**TIP**: Start with `test_pathsim_benchmark.py` before running the full benchmark!
