# Solver Quick Reference Card

## Which Solver Should I Use?

```
┌─────────────────────────────────────────────────────────────┐
│                    SOLVER DECISION TREE                      │
└─────────────────────────────────────────────────────────────┘

Do you know your timestep in advance?
│
├─ YES (Fixed timestep)
│  │
│  ├─ Need highest accuracy?        → RK4
│  ├─ PDE with TVD/SSP property?    → SSPRK33 or SSPRK22
│  └─ Just testing/prototyping?     → Euler (low accuracy)
│
└─ NO (Adaptive timestep)
   │
   └─ Non-stiff system?              → RKDP54
      Stiff system?                  → (Future: ESDIRK/BDF)
```

## Solver Comparison Table

| Solver   | Order | Stages | Adaptive | Cost  | Accuracy | Use Case              |
|----------|-------|--------|----------|-------|----------|-----------------------|
| Euler    | 1     | 1      | ✗        | ★☆☆☆☆ | ★☆☆☆☆    | Testing only          |
| SSPRK22  | 2     | 2      | ✗        | ★★☆☆☆ | ★★☆☆☆    | TVD schemes           |
| SSPRK33  | 3     | 3      | ✗        | ★★★☆☆ | ★★★☆☆    | TVD schemes (better)  |
| RK4      | 4     | 4      | ✗        | ★★★★☆ | ★★★★☆    | **Default fixed**     |
| RKDP54   | 5/4   | 7      | ✓        | ★★★★★ | ★★★★★    | **Default adaptive**  |

## Code Templates

### Template 1: Fixed-Step Integration

```rust
use nalgebra::DVector;
use rustsim::solvers::{RK4, ExplicitSolver, Solver};

// Initialize
let x0 = DVector::from_vec(vec![/* initial state */]);
let mut solver = RK4::new(x0);

// Time parameters
let dt = 0.01;
let n_steps = 1000;

// Integration loop
for _ in 0..n_steps {
    solver.buffer(dt);

    for _ in 0..solver.stages() {
        solver.step(|x, t| {
            // Your ODE: dx/dt = f(x, t)
            let mut dxdt = DVector::zeros(x.len());
            // ... compute dxdt ...
            dxdt
        }, dt);
    }
}

// Extract result
let final_state = solver.state();
```

### Template 2: Adaptive Integration

```rust
use nalgebra::DVector;
use rustsim::solvers::{RKDP54, ExplicitSolver, Solver};

// Initialize
let x0 = DVector::from_vec(vec![/* initial state */]);
let mut solver = RKDP54::with_tolerances(x0, 1e-8, 1e-4);

// Time parameters
let mut dt = 0.1;
let dt_min = 1e-6;
let dt_max = 1.0;
let mut t = 0.0;
let t_final = 10.0;

// Adaptive loop
while t < t_final {
    solver.buffer(dt);

    let mut result = None;
    for _ in 0..solver.stages() {
        result = Some(solver.step(|x, t| {
            // Your ODE: dx/dt = f(x, t)
            let mut dxdt = DVector::zeros(x.len());
            // ... compute dxdt ...
            dxdt
        }, dt));
    }

    let result = result.unwrap();

    if result.success {
        // Accept step
        t += dt;

        // Optional: store solution at t
        // solution.push((t, solver.state().clone()));

        // Adjust timestep
        if let Some(scale) = result.scale {
            dt = (dt * scale).clamp(dt_min, dt_max);
        }
    } else {
        // Reject step and retry with smaller dt
        solver.revert().unwrap();
        if let Some(scale) = result.scale {
            dt = (dt * scale).max(dt_min);
        }
    }
}
```

### Template 3: Multiple Coupled ODEs

```rust
use nalgebra::DVector;
use rustsim::solvers::{RK4, ExplicitSolver, Solver};

// Example: Lorenz system
// dx/dt = σ(y - x)
// dy/dt = x(ρ - z) - y
// dz/dt = xy - βz

let x0 = DVector::from_vec(vec![1.0, 1.0, 1.0]);
let mut solver = RK4::new(x0);

let sigma = 10.0;
let rho = 28.0;
let beta = 8.0 / 3.0;

for _ in 0..10000 {
    solver.buffer(0.01);

    for _ in 0..4 {
        solver.step(|state, _t| {
            let x = state[0];
            let y = state[1];
            let z = state[2];

            let mut dxdt = DVector::zeros(3);
            dxdt[0] = sigma * (y - x);
            dxdt[1] = x * (rho - z) - y;
            dxdt[2] = x * y - beta * z;
            dxdt
        }, 0.01);
    }
}
```

## Common Patterns

### Pattern: Save Solution at Regular Intervals

```rust
let save_interval = 0.1;
let mut next_save = save_interval;
let mut solution = vec![];

while t < t_final {
    // ... integration step ...

    if t >= next_save {
        solution.push((t, solver.state().clone()));
        next_save += save_interval;
    }
}
```

### Pattern: Event Detection (Zero-Crossing)

```rust
let mut previous_value = f(solver.state());

while t < t_final {
    // ... integration step ...

    let current_value = f(solver.state());

    if previous_value * current_value < 0.0 {
        // Zero crossing detected!
        println!("Event at t ≈ {}", t);
    }

    previous_value = current_value;
}
```

### Pattern: Stop on Condition

```rust
while t < t_final {
    // ... integration step ...

    // Stop if state exceeds threshold
    if solver.state().norm() > 1000.0 {
        println!("Solution diverged at t = {}", t);
        break;
    }
}
```

## Troubleshooting

### Problem: Solution is Inaccurate

**For fixed-step methods:**
- ✓ Decrease timestep `dt`
- ✓ Use higher-order method (RK4 instead of Euler)

**For adaptive methods:**
- ✓ Tighten tolerances: `RKDP54::with_tolerances(x0, 1e-10, 1e-6)`
- ✓ Decrease `dt_min`

### Problem: Simulation is Too Slow

**For fixed-step methods:**
- ✓ Increase timestep `dt` (if stable)
- ✓ Use lower-order method if accuracy permits

**For adaptive methods:**
- ✓ Loosen tolerances
- ✓ Increase `dt_max`
- ✓ Check if system is stiff (needs implicit solver)

### Problem: Solution is Unstable

**Signs of stiffness:**
- Very small timesteps required
- Frequent step rejections
- Better at large `dt` but unstable

**Solutions:**
- ✓ For now: use very small fixed timestep
- ✓ Future: use implicit solver (ESDIRK, BDF)

### Problem: "EmptyHistory" Error

```
Error: EmptyHistory
```

**Cause:** Called `revert()` or accessed buffered state without calling `buffer()` first

**Fix:** Always call `solver.buffer(dt)` before stepping:

```rust
// ✗ Wrong
solver.step(|x, _t| -x, dt);

// ✓ Correct
solver.buffer(dt);
solver.step(|x, _t| -x, dt);
```

## Performance Tips

### Tip 1: Pre-allocate Workspace

```rust
// ✓ Good: reuse allocation
let mut workspace = DVector::zeros(n);

for _ in 0..n_steps {
    solver.buffer(dt);
    for _ in 0..solver.stages() {
        solver.step(|x, _t| {
            // Compute into workspace
            workspace.copy_from(x);
            // ... modify workspace ...
            workspace.clone()  // Only clone when returning
        }, dt);
    }
}
```

### Tip 2: Avoid Unnecessary Clones

```rust
// ✗ Bad: clones state every access
let x = solver.state().clone();
let y = x.clone();

// ✓ Good: use references
let x = solver.state();
let y = x;  // Just a reference
```

### Tip 3: Use Fixed-Step When Possible

Adaptive methods have overhead from error estimation. If you know your timestep in advance, use RK4.

## Testing Your ODE Implementation

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_my_ode() {
        let x0 = DVector::from_vec(vec![1.0]);
        let mut solver = RK4::new(x0);

        // Integrate dx/dt = -x from 0 to 1
        for _ in 0..100 {
            solver.buffer(0.01);
            for _ in 0..4 {
                solver.step(|x, _t| -x, 0.01);
            }
        }

        // Exact solution: x(1) = exp(-1)
        let exact = (-1.0_f64).exp();
        assert_relative_eq!(
            solver.state()[0],
            exact,
            epsilon = 1e-6
        );
    }
}
```

## References

- **Full Documentation**: `README.md`
- **Mathematical Details**: `BUTCHER_TABLEAUX.md`
- **Implementation Notes**: `/home/yongkyunshin/personal/RustSim/SOLVER_IMPLEMENTATION.md`
- **Examples**: `examples/test_solvers.rs`
- **Tests**: `tests/test_solvers.rs`

## API Cheat Sheet

```rust
// Create solver
let solver = RK4::new(initial_state);
let solver = RKDP54::new(initial_state);
let solver = RKDP54::with_tolerances(initial_state, abs_tol, rel_tol);

// Solver properties
solver.order()         // Integration order
solver.stages()        // Number of stages
solver.is_adaptive()   // Has error control?
solver.is_explicit()   // Is explicit method?

// State management
solver.state()         // Get current state
solver.state_mut()     // Get mutable state
solver.set_state(x)    // Set state directly
solver.buffer(dt)      // Save state for revert
solver.revert()        // Restore buffered state
solver.reset()         // Reset to initial state

// Integration
let result = solver.step(|x, t| {
    // Return dx/dt
}, dt);

// Result fields
result.success         // Step accepted?
result.error_norm      // Error estimate
result.scale           // Timestep multiplier
```
