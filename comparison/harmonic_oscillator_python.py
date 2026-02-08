#!/usr/bin/env python3
"""
Harmonic oscillator benchmark in Python (PathSim-style)

This implements the same harmonic oscillator simulation as the RustSim example
for performance comparison.

System: d²x/dt² = -x
Solution: x(t) = x0*cos(t) + v0*sin(t)
With x0=1, v0=0: x(t) = cos(t)
"""

import time
import math
import json
import sys


class Integrator:
    """Simple integrator using Euler method"""

    def __init__(self, initial_value):
        self.state = initial_value
        self.input = 0.0

    def set_input(self, value):
        self.input = value

    def get_output(self):
        return self.state

    def step(self, dt):
        self.state += self.input * dt


class Amplifier:
    """Simple amplifier/gain block"""

    def __init__(self, gain):
        self.gain = gain
        self.input = 0.0
        self.output = 0.0

    def set_input(self, value):
        self.input = value

    def get_output(self):
        return self.output

    def update(self):
        self.output = self.gain * self.input


class HarmonicOscillator:
    """
    Harmonic oscillator simulation

    Block diagram:
      ┌──────────────────────────────────┐
      │                                  │
      ▼                                  │
     [gain: -1] ──► [velocity] ──► [position]
    """

    def __init__(self, x0, v0):
        self.gain = Amplifier(-1.0)
        self.velocity = Integrator(v0)
        self.position = Integrator(x0)
        self.time = 0.0

    def propagate(self):
        """Wire outputs to inputs"""
        self.gain.set_input(self.position.get_output())
        self.velocity.set_input(self.gain.get_output())
        self.position.set_input(self.velocity.get_output())

    def update(self):
        """Evaluate algebraic blocks"""
        self.propagate()
        self.gain.update()
        self.propagate()

    def step(self, dt):
        """Advance simulation by dt"""
        self.update()
        self.velocity.step(dt)
        self.position.step(dt)
        self.time += dt

    def run(self, duration, dt):
        """Run for duration"""
        end = self.time + duration
        while self.time < end:
            self.step(dt)

    def get_position(self):
        return self.position.get_output()

    def get_velocity(self):
        return self.velocity.get_output()

    def get_time(self):
        return self.time


def main():
    print("Harmonic Oscillator - Python Benchmark")
    print("=" * 50)
    print()

    # Simulation parameters
    x0 = 1.0
    v0 = 0.0
    dt = 0.001
    duration = 2.0 * math.pi  # One period

    print(f"System: d²x/dt² = -x")
    print(f"Initial: x(0) = {x0}, v(0) = {v0}")
    print(f"Duration: {duration:.6f} seconds")
    print(f"Time step: {dt}")
    print(f"Total steps: {int(duration / dt)}")
    print()

    # Create simulator
    sim = HarmonicOscillator(x0, v0)

    # Run benchmark
    start_time = time.perf_counter()
    sim.run(duration, dt)
    end_time = time.perf_counter()

    elapsed = end_time - start_time

    # Calculate final values and errors
    final_position = sim.get_position()
    final_velocity = sim.get_velocity()
    exact_position = math.cos(sim.get_time())
    exact_velocity = -math.sin(sim.get_time())

    position_error = abs(final_position - exact_position)
    velocity_error = abs(final_velocity - exact_velocity)

    # Energy (should be conserved at 0.5)
    energy = 0.5 * (final_position**2 + final_velocity**2)
    energy_error = abs(energy - 0.5)

    # Print results
    print("Results:")
    print("-" * 50)
    print(f"  Execution time: {elapsed:.6f} seconds")
    print(f"  Steps per second: {int(duration / dt) / elapsed:,.0f}")
    print()
    print(f"Final state (t = {sim.get_time():.6f}):")
    print(f"  Position: {final_position:.6f} (exact: {exact_position:.6f})")
    print(f"  Velocity: {final_velocity:.6f} (exact: {exact_velocity:.6f})")
    print(f"  Position error: {position_error:.2e}")
    print(f"  Velocity error: {velocity_error:.2e}")
    print(f"  Energy: {energy:.6f} (should be 0.5)")
    print(f"  Energy error: {energy_error:.2e}")
    print()

    # Output JSON for benchmark script
    results = {
        "language": "Python",
        "elapsed_time": elapsed,
        "steps": int(duration / dt),
        "steps_per_second": int(duration / dt) / elapsed,
        "final_state": {
            "time": sim.get_time(),
            "position": final_position,
            "velocity": final_velocity,
        },
        "exact_state": {
            "position": exact_position,
            "velocity": exact_velocity,
        },
        "errors": {
            "position": position_error,
            "velocity": velocity_error,
            "energy": energy_error,
        },
        "energy": energy,
    }

    # Write results to file if requested
    if len(sys.argv) > 1 and sys.argv[1] == "--json":
        with open("/tmp/python_benchmark_results.json", "w") as f:
            json.dump(results, f, indent=2)
        print("Results written to /tmp/python_benchmark_results.json")


if __name__ == "__main__":
    main()
