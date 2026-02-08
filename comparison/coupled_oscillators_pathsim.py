#!/usr/bin/env python3
"""
Coupled Oscillator Benchmark - PathSim Implementation

Simulates N=20 coupled harmonic oscillators (chain of masses connected by springs).

System equations:
    m*x_i'' = -k*(x_i - x_{i-1}) - k*(x_i - x_{i+1}) - c*x_i'

Where:
    - x_i is the position of mass i
    - m is the mass
    - k is the spring constant
    - c is the damping coefficient
    - Fixed boundary conditions: x_0 = x_{N+1} = 0

State variables: 40 total (20 positions + 20 velocities)
"""

import time
import numpy as np
import pathsim
from pathsim.blocks import ODE
from pathsim.solvers import RK4


class CoupledOscillators:
    """Coupled oscillator system for PathSim"""

    def __init__(self, n=20, m=1.0, k=1.0, c=0.1):
        self.n = n
        self.m = m
        self.k = k
        self.c = c

    def initial_state(self):
        """Initialize state with sinusoidal displacement pattern"""
        state = np.zeros(2 * self.n)

        # Set initial positions (first N elements)
        for i in range(self.n):
            state[i] = np.sin(np.pi * (i + 1) / (self.n + 1))

        # Velocities (second N elements) are zero
        # Already initialized to zero

        return state

    def derivatives(self, state, _u, _t):
        """
        Compute derivatives: [dx/dt, dv/dt]

        State vector layout:
            state[0..N] = positions x_i
            state[N..2N] = velocities v_i

        Derivatives:
            dx_i/dt = v_i
            dv_i/dt = (1/m) * Force_i

        Force_i = -k*(x_i - x_{i-1}) - k*(x_i - x_{i+1}) - c*v_i

        With boundary conditions:
            x_{-1} = 0 (left boundary)
            x_N = 0 (right boundary)
        """
        n = self.n
        dstate = np.zeros(2 * n)

        # Extract positions and velocities
        positions = state[:n]
        velocities = state[n:]

        # Compute derivatives
        for i in range(n):
            x_i = positions[i]
            v_i = velocities[i]

            # dx_i/dt = v_i
            dstate[i] = v_i

            # Compute force on oscillator i
            force = 0.0

            # Left spring force: -k*(x_i - x_{i-1})
            x_left = 0.0 if i == 0 else positions[i - 1]
            force -= self.k * (x_i - x_left)

            # Right spring force: -k*(x_i - x_{i+1})
            x_right = 0.0 if i == n - 1 else positions[i + 1]
            force -= self.k * (x_i - x_right)

            # Damping force: -c*v_i
            force -= self.c * v_i

            # dv_i/dt = force / m
            dstate[n + i] = force / self.m

        return dstate

    def energy(self, state):
        """Compute total energy (kinetic + potential)"""
        n = self.n
        positions = state[:n]
        velocities = state[n:]

        # Kinetic energy: (1/2) * m * sum(v_i^2)
        kinetic = 0.5 * self.m * np.sum(velocities**2)

        # Potential energy: (1/2) * k * sum((x_i - x_{i-1})^2 + (x_i - x_{i+1})^2)
        # With boundary conditions x_{-1} = 0, x_N = 0
        potential = 0.0
        for i in range(n):
            x_i = positions[i]

            # Left spring
            x_left = 0.0 if i == 0 else positions[i - 1]
            potential += 0.5 * self.k * (x_i - x_left)**2

            # Right spring
            x_right = 0.0 if i == n - 1 else positions[i + 1]
            potential += 0.5 * self.k * (x_i - x_right)**2

        return kinetic + potential


def main():
    print("=" * 70)
    print("Coupled Oscillators Benchmark - PathSim")
    print("=" * 70)
    print()

    # Parameters
    N = 20
    m = 1.0
    k = 100.0
    c = 0.1

    print(f"System parameters:")
    print(f"  Number of oscillators: N = {N}")
    print(f"  Mass: m = {m}")
    print(f"  Spring constant: k = {k}")
    print(f"  Damping coefficient: c = {c}")
    print(f"  State dimension: {2*N}")
    print()

    # Simulation parameters
    t_start_sim = 0.0
    t_end_sim = 100.0
    dt = 0.001

    print(f"Simulation parameters:")
    print(f"  Time span: [{t_start_sim}, {t_end_sim}] seconds")
    print(f"  Timestep: dt = {dt} seconds")
    print(f"  Total steps: {int((t_end_sim - t_start_sim) / dt)}")
    print(f"  Solver: RK4")
    print()

    # Create system
    print("Creating coupled oscillator system...")
    system = CoupledOscillators(n=N, m=m, k=k, c=c)
    initial_state = system.initial_state()

    # Create ODE block
    ode_block = ODE(
        func=system.derivatives,
        initial_value=initial_state
    )

    print(f"  State dimension: {len(initial_state)}")
    print()

    # Compute initial energy
    initial_energy = system.energy(initial_state)
    print(f"Initial energy: {initial_energy:.6f}")
    print()

    # Create simulation
    sim = pathsim.Simulation(
        blocks=[ode_block],
        connections=[],
        dt=dt,
        Solver=RK4,
        log=False
    )

    # Run simulation
    print("Running simulation...")
    wall_start = time.perf_counter()

    sim.run(t_end_sim)

    wall_end = time.perf_counter()
    wall_time = wall_end - wall_start

    print(f"Wall-clock time: {wall_time:.4f} seconds")
    print()

    # Compute final energy
    final_state = ode_block.outputs
    final_energy = system.energy(final_state)
    energy_change = abs(final_energy - initial_energy)
    energy_relative_change = energy_change / initial_energy if initial_energy > 0 else 0.0

    print(f"Results:")
    print(f"  Final energy: {final_energy:.6f}")
    print(f"  Energy change: {energy_change:.6e}")
    print(f"  Relative energy change: {energy_relative_change:.6e}")
    print()

    # Performance metrics
    steps = int((t_end_sim - t_start_sim) / dt)
    steps_per_second = steps / wall_time
    time_per_step = wall_time / steps * 1e6  # microseconds

    print(f"Performance:")
    print(f"  Steps per second: {steps_per_second:.1f}")
    print(f"  Time per step: {time_per_step:.3f} microseconds")
    print()

    print("=" * 70)


if __name__ == "__main__":
    main()
