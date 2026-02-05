# Solver Architecture Diagram

## Trait Hierarchy

```
┌─────────────────────────────────────────────────────────────┐
│                        Solver Trait                         │
│                     (base interface)                        │
├─────────────────────────────────────────────────────────────┤
│  • state() / state_mut() / set_state()                      │
│  • buffer(dt) / revert() / reset()                          │
│  • order() / stages() / is_adaptive() / is_explicit()       │
└─────────────────────────────────────────────────────────────┘
                            ▲
                            │
            ┌───────────────┴───────────────┐
            │                               │
┌───────────────────────┐       ┌───────────────────────┐
│   ExplicitSolver      │       │   ImplicitSolver      │
│                       │       │   (future)            │
├───────────────────────┤       ├───────────────────────┤
│  • step(f, dt)        │       │  • solve(f, J, dt)    │
└───────────────────────┘       └───────────────────────┘
            ▲
            │
            │
    ┌───────┴────────────────────────────────┐
    │                                        │
    │                                        │
┌───────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐
│   Euler   │  │ SSPRK22  │  │ SSPRK33  │  │   RK4    │
│           │  │          │  │          │  │          │
│ Order: 1  │  │ Order: 2 │  │ Order: 3 │  │ Order: 4 │
│ Stages: 1 │  │ Stages: 2│  │ Stages: 3│  │ Stages: 4│
│ Fixed     │  │ Fixed    │  │ Fixed    │  │ Fixed    │
│           │  │ SSP C=1  │  │ SSP C=1  │  │          │
└───────────┘  └──────────┘  └──────────┘  └──────────┘

    ┌──────────┐
    │  RKDP54  │
    │          │
    │ Order:5/4│
    │ Stages: 7│
    │ Adaptive │
    │ Error    │
    │ Control  │
    └──────────┘
```

## Data Flow: Fixed-Step Integration

```
┌─────────────┐
│ User Code   │
│ (ODE Block) │
└─────┬───────┘
      │ 1. Create solver
      │    solver = RK4::new(x0)
      ▼
┌─────────────┐
│   Solver    │
│  Instance   │
│             │
│  state: x   │
│  history: []│
│  stage: 0   │
└─────┬───────┘
      │ 2. Start step
      │    solver.buffer(dt)
      ▼
┌─────────────┐
│   Solver    │
│             │
│ state: x    │
│ history:[x] │◄─── Buffered for potential revert
│ stage: 0    │
└─────┬───────┘
      │ 3. Evaluate stages
      │    for stage in 0..N
      │       solver.step(f, dt)
      ▼
┌─────────────┐
│  RK Stage   │
│  Evaluation │
│             │
│  k[i] = f(x)│  ◄─── Compute slope
│  x += ...   │  ◄─── Update state
│  stage++    │
└─────┬───────┘
      │ 4. After all stages
      ▼
┌─────────────┐
│   Solver    │
│             │
│ state: x_n+1│  ◄─── Final state
│ stage: 0    │  ◄─── Reset for next step
└─────────────┘
```

## Data Flow: Adaptive Integration (RKDP54)

```
┌─────────────┐
│ User Code   │
└─────┬───────┘
      │ 1. Create solver
      │    solver = RKDP54::new(x0)
      ▼
┌─────────────┐
│   Solver    │
│             │
│ state: x    │
│ history: [] │
│ stage: 0    │
│ dt: dt0     │
└─────┬───────┘
      │ 2. Start step
      │    solver.buffer(dt)
      ▼
┌─────────────┐
│   Solver    │
│ history:[x] │◄─── Save for revert
└─────┬───────┘
      │ 3. Evaluate 7 stages
      │    for stage in 0..7
      ▼
┌─────────────┐
│  RK Stages  │
│  0,1,2,...6 │
│             │
│  k[0..6]    │  ◄─── Store all slopes
└─────┬───────┘
      │ 4. After stage 6
      │    (last stage)
      ▼
┌─────────────┐
│   Error     │
│ Controller  │
│             │
│ error_slope │  ◄─── Σ TR[i] * k[i]
│ error_norm  │  ◄─── max(|error| / scale)
│             │
│ success?    │  ◄─── error_norm ≤ 1.0?
│ scale       │  ◄─── β / err^(1/(p+1))
└─────┬───────┘
      │
      ├─── success = true ───┐
      │                      ▼
      │              ┌─────────────┐
      │              │   Accept    │
      │              │   Step      │
      │              │             │
      │              │ t += dt     │
      │              │ dt *= scale │
      │              └─────────────┘
      │
      └─── success = false ──┐
                             ▼
                     ┌─────────────┐
                     │   Reject    │
                     │   Step      │
                     │             │
                     │ revert()    │◄─── Restore from history
                     │ dt *= scale │◄─── Reduce timestep
                     └─────────────┘
```

## Module Organization

```
rustsim/
│
└── src/
    └── solvers/
        │
        ├── base.rs ─────────────┐
        │   • Solver trait       │
        │   • ExplicitSolver     │  ◄─── Core Traits
        │   • ImplicitSolver     │
        │   • SolverError        │
        │   • SolverStepResult   │
        │                        │
        ├── euler.rs ────────────┤
        │   • Euler struct       │
        │   • impl Solver        │  ◄─── Fixed-Step
        │   • impl ExplicitSolver│       Methods
        │   • tests              │
        │                        │
        ├── rk4.rs ──────────────┤
        │   • RK4 struct         │
        │   • impl Solver        │
        │   • impl ExplicitSolver│
        │   • tests              │
        │                        │
        ├── ssprk.rs ────────────┤
        │   • SSPRK22 struct     │
        │   • SSPRK33 struct     │
        │   • impl Solver (×2)   │
        │   • impl ExplicitSolver│
        │   • tests              │
        │                        │
        ├── rkdp54.rs ───────────┤
        │   • RKDP54 struct      │
        │   • impl Solver        │  ◄─── Adaptive
        │   • impl ExplicitSolver│       Method
        │   • error_controller() │
        │   • tests              │
        │                        │
        ├── mod.rs ──────────────┘
        │   • pub use ...
        │
        ├── README.md
        ├── BUTCHER_TABLEAUX.md
        ├── QUICK_REFERENCE.md
        └── ARCHITECTURE.md (this file)
```

## State Vector Flow

```
                    ┌──────────────┐
                    │ Initial State│
                    │   x0: [x, y] │
                    └───────┬──────┘
                            │
                    ┌───────▼──────┐
                    │  Solver Init │
                    │ state = x0   │
                    └───────┬──────┘
                            │
            ┌───────────────┼───────────────┐
            │               │               │
    ┌───────▼──────┐ ┌─────▼─────┐ ┌───────▼──────┐
    │ Buffer State │ │   Stage   │ │  Get State   │
    │ history ← x  │ │  Updates  │ │  x = state() │
    └──────────────┘ └─────┬─────┘ └──────────────┘
                           │
                    ┌──────▼──────┐
                    │ Evaluate f  │
                    │ f(x, t)     │
                    └──────┬──────┘
                           │
                    ┌──────▼──────┐
                    │Update State │
                    │ x += dt*Σk  │
                    └──────┬──────┘
                           │
            ┌──────────────┼──────────────┐
            │              │              │
    ┌───────▼──────┐ ┌─────▼─────┐ ┌─────▼─────┐
    │   Success    │ │  Failure  │ │   Reset   │
    │ Continue     │ │  Revert   │ │ x = x0    │
    │              │ │ x←history │ │           │
    └──────────────┘ └───────────┘ └───────────┘
```

## Error Estimation (RKDP54)

```
Stage Evaluations:
k0 = f(x0, t0)
k1 = f(x0 + a10*k0, t0 + c1*dt)
k2 = f(x0 + a20*k0 + a21*k1, t0 + c2*dt)
...
k6 = f(x0 + Σ a6i*ki, t0 + c6*dt)

          ↓

5th Order Solution (propagating):
x_{n+1} = x_n + dt * Σ b_i * k_i

          ↓

4th Order Solution (embedded):
x*_{n+1} = x_n + dt * Σ b*_i * k_i

          ↓

Local Truncation Error:
LTE = x_{n+1} - x*_{n+1}
    = dt * Σ (b_i - b*_i) * k_i
    = dt * Σ TR_i * k_i

          ↓

Error Norm:
scale = tol_abs + tol_rel * |x|
error_norm = max(|LTE| / scale)

          ↓

Timestep Control:
success = (error_norm ≤ 1)
scale = β / error_norm^(1/(p+1))
scale = clamp(scale, 0.1, 10.0)
dt_new = dt * scale
```

## Type Relationships

```rust
// Core trait that all solvers implement
pub trait Solver: Send + Sync {
    fn state(&self) -> &DVector<f64>;
    // ... other methods
}

// Explicit methods extend Solver
pub trait ExplicitSolver: Solver {
    fn step<F>(&mut self, f: F, dt: f64) -> SolverStepResult
    where
        F: FnMut(&DVector<f64>, f64) -> DVector<f64>;
}

// Concrete types implement both traits
impl Solver for RK4 { /* ... */ }
impl ExplicitSolver for RK4 { /* ... */ }

impl Solver for RKDP54 { /* ... */ }
impl ExplicitSolver for RKDP54 { /* ... */ }

// etc.
```

## Memory Layout

```
Euler (smallest):
┌─────────────────┐
│ state: DVector  │  ← Current state
│ initial: DVector│  ← Initial condition
│ history: Deque  │  ← VecDeque<DVector>
└─────────────────┘

RK4 (multi-stage):
┌─────────────────┐
│ state: DVector  │
│ initial: DVector│
│ history: Deque  │
│ slopes: Vec[4]  │  ← Vec<DVector> (k0, k1, k2, k3)
│ stage: usize    │  ← Current stage index
└─────────────────┘

RKDP54 (adaptive):
┌─────────────────┐
│ state: DVector  │
│ initial: DVector│
│ history: Deque  │
│ slopes: Vec[7]  │  ← Vec<DVector> (k0..k6)
│ stage: usize    │
│ tol_abs: f64    │  ← Error tolerances
│ tol_rel: f64    │
│ beta: f64       │  ← Safety factor (0.9)
└─────────────────┘
```

## Comparison with PathSim

```
PathSim (Python)          RustSim (Rust)
───────────────────────   ───────────────────────
Solver base class    ───▶ Solver trait
ExplicitSolver       ───▶ ExplicitSolver trait
ImplicitSolver       ───▶ ImplicitSolver trait

euler.py             ───▶ euler.rs
rkdp54.py            ───▶ rkdp54.rs
ssprk22.py           ───▶ ssprk.rs::SSPRK22
ssprk33.py           ───▶ ssprk.rs::SSPRK33
_rungekutta.py (RK4) ───▶ rk4.rs

numpy.array          ───▶ nalgebra::DVector
deque                ───▶ std::collections::VecDeque
```

## Thread Safety

All solvers are `Send + Sync`:

```rust
impl Solver for RK4 // automatically Send + Sync
    because:
    - DVector<f64>: Send + Sync
    - VecDeque<T>: Send + Sync where T: Send + Sync
    - Vec<T>: Send + Sync where T: Send + Sync
    - usize: Send + Sync
```

This allows:
- Parallel evaluation of independent ODEs
- Safe sharing across threads (with proper synchronization)
- Integration in multi-threaded simulations

## Summary

The solver architecture provides:

✓ **Clean separation of concerns** (traits vs implementations)
✓ **Type safety** (Rust's type system prevents misuse)
✓ **Extensibility** (easy to add new methods)
✓ **Performance** (zero-cost abstractions)
✓ **Testability** (each component independently testable)
✓ **Documentation** (extensive inline and external docs)
