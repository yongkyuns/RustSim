# Numerical Integration Solvers

This module provides various numerical integration methods for solving ordinary differential equations (ODEs).

## Available Solvers

### Explicit Methods

#### Euler (Forward Euler)
- **Order**: 1
- **Stages**: 1
- **Type**: Fixed timestep, explicit
- **Use Case**: Testing only; generally too inaccurate for production use

The simplest integration method. While cheap per step, it requires very small timesteps for stability and accuracy, making higher-order methods more efficient in practice.

```rust
use nalgebra::DVector;
use rustsim::solvers::{Euler, ExplicitSolver, Solver};

let x0 = DVector::from_vec(vec![1.0]);
let mut solver = Euler::new(x0);

solver.buffer(dt);
solver.step(|x, _t| -x, dt);
```

#### RK4 (Classic Runge-Kutta)
- **Order**: 4
- **Stages**: 4
- **Type**: Fixed timestep, explicit
- **Use Case**: Standard choice for non-stiff ODEs with fixed timestep

The workhorse of fixed-step integration. Provides excellent accuracy-per-cost for smooth, non-stiff systems. Recommended default for deterministic simulations with fixed timestep.

```rust
use nalgebra::DVector;
use rustsim::solvers::{RK4, ExplicitSolver, Solver};

let x0 = DVector::from_vec(vec![1.0]);
let mut solver = RK4::new(x0);

solver.buffer(dt);
for _ in 0..4 {
    solver.step(|x, _t| -x, dt);
}
```

#### RKDP54 (Dormand-Prince 5(4))
- **Order**: 5 (propagating) / 4 (embedded)
- **Stages**: 7
- **Type**: Adaptive timestep, explicit
- **Use Case**: Default for non-stiff ODEs with adaptive timestepping

The industry standard adaptive solver, used in MATLAB's `ode45`. Automatically adjusts timestep based on error estimates. Recommended default for general-purpose integration when you don't know the system dynamics in advance.

```rust
use nalgebra::DVector;
use rustsim::solvers::{RKDP54, ExplicitSolver, Solver};

let x0 = DVector::from_vec(vec![1.0]);
let mut solver = RKDP54::new(x0);
// Or with custom tolerances:
// let mut solver = RKDP54::with_tolerances(x0, 1e-8, 1e-4);

let mut dt = 0.1;
loop {
    solver.buffer(dt);

    let mut result = None;
    for _ in 0..7 {
        result = Some(solver.step(|x, _t| -x, dt));
    }

    let result = result.unwrap();

    if result.success {
        // Accept step and adjust timestep
        if let Some(scale) = result.scale {
            dt *= scale;
        }
    } else {
        // Reject step and retry with smaller timestep
        solver.revert().unwrap();
        if let Some(scale) = result.scale {
            dt *= scale;
        }
    }
}
```

#### SSPRK22 (SSP Runge-Kutta 2-2)
- **Order**: 2
- **Stages**: 2
- **Type**: Fixed timestep, explicit
- **SSP Coefficient**: 1
- **Use Case**: TVD schemes for hyperbolic PDEs

Strong Stability Preserving method. Preserves total variation diminishing (TVD) properties when used with appropriate spatial discretizations. Use for method-of-lines discretizations of conservation laws.

```rust
use nalgebra::DVector;
use rustsim::solvers::{SSPRK22, ExplicitSolver, Solver};

let x0 = DVector::from_vec(vec![1.0]);
let mut solver = SSPRK22::new(x0);

solver.buffer(dt);
for _ in 0..2 {
    solver.step(|x, _t| spatial_operator(x), dt);
}
```

#### SSPRK33 (SSP Runge-Kutta 3-3)
- **Order**: 3
- **Stages**: 3
- **Type**: Fixed timestep, explicit
- **SSP Coefficient**: 1
- **Use Case**: TVD schemes for hyperbolic PDEs

The optimal three-stage SSP method. Standard choice for shock-capturing schemes (WENO, ENO) in computational fluid dynamics.

```rust
use nalgebra::DVector;
use rustsim::solvers::{SSPRK33, ExplicitSolver, Solver};

let x0 = DVector::from_vec(vec![1.0]);
let mut solver = SSPRK33::new(x0);

solver.buffer(dt);
for _ in 0..3 {
    solver.step(|x, _t| spatial_operator(x), dt);
}
```

## Solver Selection Guide

### For Non-Stiff Systems

1. **Need adaptive timestep?** → Use **RKDP54**
2. **Fixed timestep, smooth dynamics?** → Use **RK4**
3. **PDE method-of-lines with TVD property?** → Use **SSPRK33** or **SSPRK22**

### For Stiff Systems

Currently, only explicit methods are available. For stiff systems (fast transients, high-gain feedback, chemical kinetics), implicit methods like ESDIRK or BDF will be needed (future work).

### Signs of Stiffness

If you observe:
- Very small timesteps required for stability
- Frequent step rejections in adaptive methods
- Better accuracy at larger timesteps but instability

Then your system is likely stiff and needs an implicit solver.

## Solver Traits

All solvers implement the `Solver` trait:

```rust
pub trait Solver {
    fn state(&self) -> &DVector<f64>;
    fn state_mut(&mut self) -> &mut DVector<f64>;
    fn set_state(&mut self, state: DVector<f64>);
    fn buffer(&mut self, dt: f64);
    fn revert(&mut self) -> Result<(), SolverError>;
    fn reset(&mut self);
    fn order(&self) -> usize;
    fn stages(&self) -> usize;
    fn is_adaptive(&self) -> bool;
    fn is_explicit(&self) -> bool;
}
```

Explicit methods additionally implement `ExplicitSolver`:

```rust
pub trait ExplicitSolver: Solver {
    fn step<F>(&mut self, f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>;
}
```

## Error Handling

Solvers use the `SolverError` enum for error reporting:

- `ConvergenceFailure(usize)`: Implicit solver failed to converge
- `TimestepTooSmall { dt, dt_min }`: Timestep below minimum threshold
- `EmptyHistory`: Attempted to revert with no buffered state
- `InvalidStage { stage, max_stages }`: Invalid stage index

## Step Results

The `SolverStepResult` struct contains:

```rust
pub struct SolverStepResult {
    pub success: bool,      // True if step should be accepted
    pub error_norm: f64,    // Estimated error (for adaptive methods)
    pub scale: Option<f64>, // Timestep scale factor (for adaptive methods)
}
```

For fixed-step methods, `success` is always `true` and `scale` is `None`.

For adaptive methods like RKDP54, the error controller returns:
- `success = true` if `error_norm <= 1.0`
- `scale` = suggested timestep multiplier (clamped to [0.1, 10.0])

## Implementation Details

### Butcher Tableaux

All Runge-Kutta methods are defined by their Butcher tableau:

```
c₁ | a₁₁  a₁₂  ...
c₂ | a₂₁  a₂₂  ...
...
---+---------------
   | b₁   b₂   ...
```

Where:
- `c` = evaluation time ratios
- `a` = stage coefficients
- `b` = output weights

### Error Control (RKDP54)

The adaptive timestep controller uses:

```
error_norm = max(|dt * error_slope| / scale)
scale = tolerance_abs + tolerance_rel * |x|
timestep_scale = beta / error_norm^(1/(p+1))
```

Where:
- `p` = min(propagating_order, embedded_order) = 4
- `beta` = safety factor = 0.9
- Default tolerances: `abs = 1e-8`, `rel = 1e-4`

## Testing

Run solver tests:

```bash
cargo test --lib solvers
```

All solvers are tested against:
- Exponential decay: `dx/dt = -x`
- Harmonic oscillator: `d²x/dt² = -x`

Expected convergence rates are verified for each method.

## References

1. Hairer, E., Nørsett, S. P., & Wanner, G. (1993). "Solving Ordinary Differential Equations I: Nonstiff Problems". Springer.

2. Dormand, J. R., & Prince, P. J. (1980). "A family of embedded Runge-Kutta formulae". Journal of Computational and Applied Mathematics, 6(1), 19-26.

3. Shu, C.-W., & Osher, S. (1988). "Efficient implementation of essentially non-oscillatory shock-capturing schemes". Journal of Computational Physics, 77(2), 439-471.

4. Gottlieb, S., Shu, C.-W., & Tadmor, E. (2001). "Strong stability-preserving high-order time discretization methods". SIAM Review, 43(1), 89-112.

5. Shampine, L. F., & Reichelt, M. W. (1997). "The MATLAB ODE Suite". SIAM Journal on Scientific Computing, 18(1), 1-22.
