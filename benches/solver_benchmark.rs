//! Solver benchmarks
//!
//! Benchmarks numerical integration performance for various solvers.

use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use nalgebra::DVector;
use rustsim::solvers::{ExplicitSolver, RK4, Solver};

/// Simple exponential decay: dx/dt = -k*x
fn exponential_decay(x: &DVector<f64>, _t: f64, k: f64) -> DVector<f64> {
    -k * x
}

/// Benchmark RK4 solver with different state vector sizes
fn bench_rk4_step(c: &mut Criterion) {
    let mut group = c.benchmark_group("RK4 Step");

    for size in [1, 5, 10, 20, 50, 100].iter() {
        let initial = DVector::from_element(*size, 1.0);
        let dt = 0.001;
        let k = 0.5;

        group.bench_with_input(BenchmarkId::new("state_size", size), size, |b, _| {
            b.iter(|| {
                let mut solver = RK4::new(initial.clone());
                solver.buffer(dt);

                // Perform one complete RK4 step (4 stages)
                for _ in 0..4 {
                    solver.step(|x, t| exponential_decay(x, t, k), black_box(dt));
                }
            });
        });
    }

    group.finish();
}

/// Benchmark RK4 solver for a fixed problem over many steps
fn bench_rk4_integration(c: &mut Criterion) {
    let initial = DVector::from_element(10, 1.0);
    let dt = 0.001;
    let k = 0.5;
    let steps = 1000;

    c.bench_function("RK4 1000 steps (n=10)", |b| {
        b.iter(|| {
            let mut solver = RK4::new(initial.clone());

            for _ in 0..steps {
                solver.buffer(dt);
                for _ in 0..4 {
                    solver.step(|x, t| exponential_decay(x, t, k), black_box(dt));
                }
            }

            black_box(solver.state());
        });
    });
}

/// Benchmark harmonic oscillator (2D system)
fn bench_rk4_harmonic_oscillator(c: &mut Criterion) {
    let initial = DVector::from_vec(vec![1.0, 0.0]); // [x, dx/dt]
    let dt = 0.001;
    let omega = 2.0 * std::f64::consts::PI;

    c.bench_function("RK4 Harmonic Oscillator", |b| {
        b.iter(|| {
            let mut solver = RK4::new(initial.clone());

            for _ in 0..1000 {
                solver.buffer(dt);
                for _ in 0..4 {
                    solver.step(
                        |state, _t| {
                            DVector::from_vec(vec![state[1], -omega * omega * state[0]])
                        },
                        black_box(dt),
                    );
                }
            }

            black_box(solver.state());
        });
    });
}

/// Benchmark solver buffer and revert operations
fn bench_solver_buffer_revert(c: &mut Criterion) {
    let initial = DVector::from_element(10, 1.0);
    let dt = 0.001;

    c.bench_function("RK4 Buffer & Revert", |b| {
        let mut solver = RK4::new(initial.clone());

        b.iter(|| {
            solver.buffer(black_box(dt));
            let _ = solver.revert();
        });
    });
}

criterion_group!(
    benches,
    bench_rk4_step,
    bench_rk4_integration,
    bench_rk4_harmonic_oscillator,
    bench_solver_buffer_revert
);
criterion_main!(benches);
