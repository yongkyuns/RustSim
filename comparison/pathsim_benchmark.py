#!/usr/bin/env python3
"""
PathSim Comprehensive Benchmark Suite

This benchmark suite tests PathSim's performance on computationally challenging systems:
1. Robertson's Stiff Chemical Kinetics (very stiff ODE system)
2. Van der Pol Oscillator with large mu (stiff system)
3. Lorenz Attractor (chaotic system, long duration)
4. Coupled Oscillator Chain (large coupled system, 20 oscillators)

Each benchmark tests multiple solvers and records:
- Wall-clock execution time
- Final state values (for accuracy comparison)
- Number of simulation steps completed

Results are saved to a JSON file and printed to stdout.
"""

import sys
import os
import time
import json
import numpy as np
from datetime import datetime

# Add PathSim to path if needed
pathsim_path = "/home/yongkyunshin/personal/pathsim/src"
if pathsim_path not in sys.path:
    sys.path.insert(0, pathsim_path)

from pathsim import Simulation, Connection
from pathsim.blocks import Scope, ODE, Integrator, Adder, Amplifier, Multiplier, Constant
from pathsim.solvers import (
    RK4, RKBS32, RKCK54, RKDP54,  # Explicit RK methods
    ESDIRK32, ESDIRK43, ESDIRK54,  # Implicit ESDIRK methods
    GEAR21, GEAR32, GEAR43, GEAR52A,  # Implicit GEAR methods
    BDF3  # BDF method
)


# ======================================================================================
# BENCHMARK 1: Robertson's Stiff Chemical Kinetics
# ======================================================================================

def benchmark_robertson(solver_class, solver_name, duration=100):
    """
    Robertson's stiff chemical kinetics problem.
    Classic stiff ODE test case with time scales varying over 10+ orders of magnitude.
    """
    print(f"  Running Robertson with {solver_name}...", end="", flush=True)

    # Parameters
    a, b, c = 0.04, 1e4, 3e7

    # Initial condition
    x0 = np.array([1.0, 0.0, 0.0])

    def func(x, u, t):
        return np.array([
            -a*x[0] + b*x[1]*x[2],
             a*x[0] - b*x[1]*x[2] - c*x[1]**2,
                                    c*x[1]**2
        ])

    # Jacobian for implicit solvers
    def jac(x, u, t):
        return np.array([
            [-a,              b*x[2],         b*x[1]],
            [a,               -b*x[2]-2*c*x[1], -b*x[1]],
            [0,               2*c*x[1],       0]
        ])

    # Build system
    Rob = ODE(func, x0, jac)
    Sco = Scope(labels=["x", "y", "z"])

    blocks = [Rob, Sco]
    connections = [
        Connection(Rob[0], Sco[0]),
        Connection(Rob[1], Sco[1]),
        Connection(Rob[2], Sco[2])
    ]

    # Create simulation
    Sim = Simulation(
        blocks,
        connections,
        dt=0.001,
        Solver=solver_class,
        tolerance_lte_abs=1e-6,
        tolerance_lte_rel=1e-4,
        tolerance_fpi=1e-9
    )

    # Time the simulation
    start_time = time.perf_counter()
    Sim.run(duration)
    end_time = time.perf_counter()

    elapsed = end_time - start_time

    # Get final state
    rec_time, rec_data = Sco.read()
    if rec_time is not None and rec_data is not None:
        final_state = {
            "x": float(rec_data[0][-1]),
            "y": float(rec_data[1][-1]),
            "z": float(rec_data[2][-1]),
            "time": float(rec_time[-1])
        }
        num_steps = len(rec_time)
    else:
        final_state = None
        num_steps = 0

    print(f" {elapsed:.3f}s ({num_steps} steps)")

    return {
        "solver": solver_name,
        "elapsed_time": elapsed,
        "num_steps": num_steps,
        "final_state": final_state
    }


# ======================================================================================
# BENCHMARK 2: Van der Pol Oscillator (Stiff)
# ======================================================================================

def benchmark_vanderpol(solver_class, solver_name, mu=1000, duration=3000):
    """
    Van der Pol oscillator with large mu parameter.
    Creates a very stiff system with sharp transitions.
    """
    print(f"  Running Van der Pol (mu={mu}) with {solver_name}...", end="", flush=True)

    # Initial condition
    x0 = np.array([2.0, 0.0])

    def func(x, u, t):
        return np.array([x[1], mu*(1 - x[0]**2)*x[1] - x[0]])

    # Analytical Jacobian
    def jac(x, u, t):
        return np.array([
            [0, 1],
            [-mu*2*x[0]*x[1]-1, mu*(1 - x[0]**2)]
        ])

    # Build system
    VDP = ODE(func, x0, jac)
    Sco = Scope(labels=["x1", "x2"])

    blocks = [VDP, Sco]
    connections = [
        Connection(VDP[0], Sco[0]),
        Connection(VDP[1], Sco[1])
    ]

    # Create simulation
    Sim = Simulation(
        blocks,
        connections,
        Solver=solver_class,
        tolerance_lte_abs=1e-5,
        tolerance_lte_rel=1e-3,
        tolerance_fpi=1e-8
    )

    # Time the simulation
    start_time = time.perf_counter()
    Sim.run(duration)
    end_time = time.perf_counter()

    elapsed = end_time - start_time

    # Get final state
    rec_time, rec_data = Sco.read()
    if rec_time is not None and rec_data is not None:
        final_state = {
            "x1": float(rec_data[0][-1]),
            "x2": float(rec_data[1][-1]),
            "time": float(rec_time[-1])
        }
        num_steps = len(rec_time)
    else:
        final_state = None
        num_steps = 0

    print(f" {elapsed:.3f}s ({num_steps} steps)")

    return {
        "solver": solver_name,
        "elapsed_time": elapsed,
        "num_steps": num_steps,
        "final_state": final_state
    }


# ======================================================================================
# BENCHMARK 3: Lorenz Attractor (Chaotic System)
# ======================================================================================

def benchmark_lorenz(solver_class, solver_name, duration=200):
    """
    Lorenz attractor - chaotic system requiring high precision for long simulations.
    """
    print(f"  Running Lorenz with {solver_name}...", end="", flush=True)

    # Parameters
    sigma, rho, beta = 10, 28, 8/3

    # Initial conditions
    x0, y0, z0 = 1.0, 1.0, 1.0

    # Integrators
    itg_x = Integrator(x0)
    itg_y = Integrator(y0)
    itg_z = Integrator(z0)

    # Components for dx/dt = sigma * (y - x)
    amp_sigma = Amplifier(sigma)
    add_x = Adder("+-")

    # Components for dy/dt = x * (rho - z) - y
    cns_rho = Constant(rho)
    add_rho_z = Adder("+-")
    mul_x_rho_z = Multiplier()
    add_y = Adder("-+")

    # Components for dz/dt = x * y - beta * z
    mul_xy = Multiplier()
    amp_beta = Amplifier(beta)
    add_z = Adder("+-")

    # Scope
    sco = Scope(labels=["x", "y", "z"])

    # All blocks
    blocks = [
        itg_x, itg_y, itg_z,
        amp_sigma, add_x,
        cns_rho, add_rho_z, mul_x_rho_z, add_y,
        mul_xy, amp_beta, add_z,
        sco
    ]

    # Connections
    connections = [
        # Output signals
        Connection(itg_x, add_x[1], mul_x_rho_z[0], mul_xy[0], sco[0]),
        Connection(itg_y, add_x[0], add_y[0], mul_xy[1], sco[1]),
        Connection(itg_z, add_rho_z[1], amp_beta, sco[2]),

        # dx/dt path
        Connection(add_x, amp_sigma),
        Connection(amp_sigma, itg_x),

        # dy/dt path
        Connection(cns_rho, add_rho_z[0]),
        Connection(add_rho_z, mul_x_rho_z[1]),
        Connection(mul_x_rho_z, add_y[1]),
        Connection(add_y, itg_y),

        # dz/dt path
        Connection(mul_xy, add_z[0]),
        Connection(amp_beta, add_z[1]),
        Connection(add_z, itg_z)
    ]

    # Create simulation
    Sim = Simulation(
        blocks,
        connections,
        Solver=solver_class,
        tolerance_lte_rel=1e-6,
        tolerance_lte_abs=1e-8,
        tolerance_fpi=1e-8
    )

    # Time the simulation
    start_time = time.perf_counter()
    Sim.run(duration)
    end_time = time.perf_counter()

    elapsed = end_time - start_time

    # Get final state
    rec_time, rec_data = sco.read()
    if rec_time is not None and rec_data is not None:
        final_state = {
            "x": float(rec_data[0][-1]),
            "y": float(rec_data[1][-1]),
            "z": float(rec_data[2][-1]),
            "time": float(rec_time[-1])
        }
        num_steps = len(rec_time)
    else:
        final_state = None
        num_steps = 0

    print(f" {elapsed:.3f}s ({num_steps} steps)")

    return {
        "solver": solver_name,
        "elapsed_time": elapsed,
        "num_steps": num_steps,
        "final_state": final_state
    }


# ======================================================================================
# BENCHMARK 4: Large Coupled Oscillator Chain
# ======================================================================================

def benchmark_coupled_oscillators(solver_class, solver_name, n_oscillators=20, duration=50):
    """
    Chain of N coupled harmonic oscillators.
    Tests performance on large coupled systems.

    Each oscillator is coupled to its neighbors with spring constants.
    """
    print(f"  Running {n_oscillators} Coupled Oscillators with {solver_name}...", end="", flush=True)

    # Parameters
    m = 1.0  # mass
    k_internal = 2.0  # spring constant within each oscillator
    k_coupling = 0.5  # spring constant between oscillators
    c = 0.1  # damping

    # Initial conditions (alternating displacements)
    x0 = np.zeros(2 * n_oscillators)
    for i in range(n_oscillators):
        x0[2*i] = np.sin(i * np.pi / n_oscillators)  # position
        x0[2*i + 1] = 0  # velocity

    def func(x, u, t):
        """
        State vector: [x0, v0, x1, v1, x2, v2, ...]
        For each oscillator i:
        - x[2*i] = position
        - x[2*i+1] = velocity
        """
        dx = np.zeros_like(x)

        for i in range(n_oscillators):
            pos_idx = 2 * i
            vel_idx = 2 * i + 1

            # Position derivative is velocity
            dx[pos_idx] = x[vel_idx]

            # Acceleration from internal spring and damping
            force = -k_internal * x[pos_idx] - c * x[vel_idx]

            # Coupling with left neighbor
            if i > 0:
                force += k_coupling * (x[2*(i-1)] - x[pos_idx])

            # Coupling with right neighbor
            if i < n_oscillators - 1:
                force += k_coupling * (x[2*(i+1)] - x[pos_idx])

            # Velocity derivative is acceleration
            dx[vel_idx] = force / m

        return dx

    # Jacobian for implicit solvers
    def jac(x, u, t):
        J = np.zeros((2*n_oscillators, 2*n_oscillators))

        for i in range(n_oscillators):
            pos_idx = 2 * i
            vel_idx = 2 * i + 1

            # d(dx_i)/dx_i
            J[pos_idx, vel_idx] = 1.0

            # d(dv_i)/dx_i
            total_k = k_internal
            if i > 0:
                total_k += k_coupling
            if i < n_oscillators - 1:
                total_k += k_coupling

            J[vel_idx, pos_idx] = -total_k / m
            J[vel_idx, vel_idx] = -c / m

            # Coupling terms
            if i > 0:
                J[vel_idx, 2*(i-1)] = k_coupling / m
            if i < n_oscillators - 1:
                J[vel_idx, 2*(i+1)] = k_coupling / m

        return J

    # Build system
    Osc = ODE(func, x0, jac)

    # Create scope to monitor first, middle, and last oscillators
    Sco = Scope(labels=[f"osc_0_pos", f"osc_{n_oscillators//2}_pos", f"osc_{n_oscillators-1}_pos"])

    blocks = [Osc, Sco]
    connections = [
        Connection(Osc[0], Sco[0]),  # First oscillator position
        Connection(Osc[2*(n_oscillators//2)], Sco[1]),  # Middle oscillator position
        Connection(Osc[2*(n_oscillators-1)], Sco[2])  # Last oscillator position
    ]

    # Create simulation
    Sim = Simulation(
        blocks,
        connections,
        Solver=solver_class,
        tolerance_lte_abs=1e-6,
        tolerance_lte_rel=1e-4,
        tolerance_fpi=1e-8
    )

    # Time the simulation
    start_time = time.perf_counter()
    Sim.run(duration)
    end_time = time.perf_counter()

    elapsed = end_time - start_time

    # Get final state
    rec_time, rec_data = Sco.read()
    if rec_time is not None and rec_data is not None:
        final_state = {
            "osc_0_pos": float(rec_data[0][-1]),
            f"osc_{n_oscillators//2}_pos": float(rec_data[1][-1]),
            f"osc_{n_oscillators-1}_pos": float(rec_data[2][-1]),
            "time": float(rec_time[-1])
        }
        num_steps = len(rec_time)
    else:
        final_state = None
        num_steps = 0

    print(f" {elapsed:.3f}s ({num_steps} steps)")

    return {
        "solver": solver_name,
        "elapsed_time": elapsed,
        "num_steps": num_steps,
        "final_state": final_state
    }


# ======================================================================================
# MAIN BENCHMARK RUNNER
# ======================================================================================

def run_all_benchmarks():
    """
    Run all benchmarks with multiple solvers and save results.
    """

    print("="*80)
    print("PathSim Comprehensive Benchmark Suite")
    print("="*80)
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    # Define solvers to test
    # Note: Some solvers may be better suited for different problem types
    solvers = [
        # Explicit RK methods (good for non-stiff problems)
        (RK4, "RK4"),
        (RKBS32, "RKBS32"),
        (RKCK54, "RKCK54"),
        (RKDP54, "RKDP54"),

        # Implicit ESDIRK methods (good for stiff problems)
        (ESDIRK32, "ESDIRK32"),
        (ESDIRK43, "ESDIRK43"),
        (ESDIRK54, "ESDIRK54"),

        # Implicit GEAR methods (excellent for stiff problems)
        (GEAR21, "GEAR21"),
        (GEAR32, "GEAR32"),
        (GEAR43, "GEAR43"),
        (GEAR52A, "GEAR52A"),

        # BDF method (also good for stiff problems)
        (BDF3, "BDF3"),
    ]

    results = {
        "timestamp": datetime.now().isoformat(),
        "benchmarks": {}
    }

    # Benchmark 1: Robertson (stiff)
    print("\n" + "="*80)
    print("BENCHMARK 1: Robertson's Stiff Chemical Kinetics (duration=100)")
    print("="*80)
    results["benchmarks"]["robertson"] = []

    for solver_class, solver_name in solvers:
        try:
            result = benchmark_robertson(solver_class, solver_name, duration=100)
            results["benchmarks"]["robertson"].append(result)
        except Exception as e:
            print(f"  FAILED: {solver_name} - {str(e)}")
            results["benchmarks"]["robertson"].append({
                "solver": solver_name,
                "error": str(e)
            })

    # Benchmark 2: Van der Pol (stiff)
    print("\n" + "="*80)
    print("BENCHMARK 2: Van der Pol Oscillator (mu=1000, duration=3000)")
    print("="*80)
    results["benchmarks"]["vanderpol"] = []

    # Only test implicit solvers for Van der Pol (too stiff for explicit methods)
    stiff_solvers = [
        (ESDIRK32, "ESDIRK32"),
        (ESDIRK43, "ESDIRK43"),
        (ESDIRK54, "ESDIRK54"),
        (GEAR21, "GEAR21"),
        (GEAR32, "GEAR32"),
        (GEAR43, "GEAR43"),
        (GEAR52A, "GEAR52A"),
        (BDF3, "BDF3"),
    ]

    for solver_class, solver_name in stiff_solvers:
        try:
            result = benchmark_vanderpol(solver_class, solver_name, mu=1000, duration=3000)
            results["benchmarks"]["vanderpol"].append(result)
        except Exception as e:
            print(f"  FAILED: {solver_name} - {str(e)}")
            results["benchmarks"]["vanderpol"].append({
                "solver": solver_name,
                "error": str(e)
            })

    # Benchmark 3: Lorenz (chaotic)
    print("\n" + "="*80)
    print("BENCHMARK 3: Lorenz Attractor (duration=200)")
    print("="*80)
    results["benchmarks"]["lorenz"] = []

    for solver_class, solver_name in solvers:
        try:
            result = benchmark_lorenz(solver_class, solver_name, duration=200)
            results["benchmarks"]["lorenz"].append(result)
        except Exception as e:
            print(f"  FAILED: {solver_name} - {str(e)}")
            results["benchmarks"]["lorenz"].append({
                "solver": solver_name,
                "error": str(e)
            })

    # Benchmark 4: Coupled Oscillators (large system)
    print("\n" + "="*80)
    print("BENCHMARK 4: Coupled Oscillator Chain (20 oscillators, duration=50)")
    print("="*80)
    results["benchmarks"]["coupled_oscillators"] = []

    for solver_class, solver_name in solvers:
        try:
            result = benchmark_coupled_oscillators(solver_class, solver_name, n_oscillators=20, duration=50)
            results["benchmarks"]["coupled_oscillators"].append(result)
        except Exception as e:
            print(f"  FAILED: {solver_name} - {str(e)}")
            results["benchmarks"]["coupled_oscillators"].append({
                "solver": solver_name,
                "error": str(e)
            })

    # Save results to JSON file
    output_file = "/home/yongkyunshin/personal/RustSim/benchmarks/pathsim_benchmark_results.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)

    print("\n" + "="*80)
    print("BENCHMARK SUMMARY")
    print("="*80)

    for benchmark_name, benchmark_results in results["benchmarks"].items():
        print(f"\n{benchmark_name.upper()}:")
        print("-" * 80)

        # Find successful runs
        successful = [r for r in benchmark_results if "error" not in r]

        if successful:
            # Sort by elapsed time
            successful.sort(key=lambda x: x["elapsed_time"])

            print(f"{'Solver':<15} {'Time (s)':<12} {'Steps':<10} {'Speed Rank'}")
            print("-" * 80)
            for i, result in enumerate(successful):
                print(f"{result['solver']:<15} {result['elapsed_time']:>10.3f}  {result['num_steps']:>8}  #{i+1}")
        else:
            print("No successful runs!")

    print("\n" + "="*80)
    print(f"Results saved to: {output_file}")
    print(f"Completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*80)

    return results


if __name__ == "__main__":
    results = run_all_benchmarks()
