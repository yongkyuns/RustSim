#!/usr/bin/env python3
"""
Quick test of the PathSim benchmark suite.
Runs a single benchmark with one solver to verify everything works.
"""

import sys
sys.path.insert(0, "/home/yongkyunshin/personal/pathsim/src")

from pathsim_benchmark import benchmark_robertson, benchmark_lorenz, benchmark_vanderpol, benchmark_coupled_oscillators
from pathsim.solvers import GEAR52A, RKDP54

print("Testing PathSim benchmarks...")
print("="*80)

# Test 1: Robertson with GEAR52A (implicit, good for stiff)
print("\n1. Testing Robertson with GEAR52A...")
try:
    result = benchmark_robertson(GEAR52A, "GEAR52A", duration=10)
    print(f"   SUCCESS: {result['elapsed_time']:.3f}s, {result['num_steps']} steps")
    print(f"   Final state: {result['final_state']}")
except Exception as e:
    print(f"   FAILED: {e}")

# Test 2: Lorenz with RKDP54 (explicit, good for non-stiff)
print("\n2. Testing Lorenz with RKDP54...")
try:
    result = benchmark_lorenz(RKDP54, "RKDP54", duration=20)
    print(f"   SUCCESS: {result['elapsed_time']:.3f}s, {result['num_steps']} steps")
    print(f"   Final state: {result['final_state']}")
except Exception as e:
    print(f"   FAILED: {e}")

# Test 3: Coupled oscillators with GEAR52A
print("\n3. Testing 10 Coupled Oscillators with GEAR52A...")
try:
    result = benchmark_coupled_oscillators(GEAR52A, "GEAR52A", n_oscillators=10, duration=10)
    print(f"   SUCCESS: {result['elapsed_time']:.3f}s, {result['num_steps']} steps")
    print(f"   Final state: {result['final_state']}")
except Exception as e:
    print(f"   FAILED: {e}")

# Test 4: Van der Pol with GEAR52A (only test with smaller mu for quick test)
print("\n4. Testing Van der Pol (mu=100) with GEAR52A...")
try:
    result = benchmark_vanderpol(GEAR52A, "GEAR52A", mu=100, duration=300)
    print(f"   SUCCESS: {result['elapsed_time']:.3f}s, {result['num_steps']} steps")
    print(f"   Final state: {result['final_state']}")
except Exception as e:
    print(f"   FAILED: {e}")

print("\n" + "="*80)
print("All tests completed!")
print("If all tests succeeded, the full benchmark should work correctly.")
print("Run 'python3 pathsim_benchmark.py' to execute the complete benchmark suite.")
