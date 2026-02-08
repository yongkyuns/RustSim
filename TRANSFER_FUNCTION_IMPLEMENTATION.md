# Transfer Function Implementation for RustSim

## Overview

Implemented a comprehensive `TransferFunction` block for RustSim that converts numerator-denominator polynomial representations to state-space form and simulates them using RK4 integration.

## Implementation Details

### File Structure

- **`src/blocks/transfer_function.rs`**: Main implementation
  - `TransferFunction<N>` struct with const generic for system order
  - Observable canonical form state-space realization (matching scipy.signal)
  - RK4 integration for state evolution
  - Full Block and DynamicBlock trait implementations

- **`tests/test_transfer_function.rs`**: Comprehensive test suite
  - 15 test cases covering various transfer function types
  - Validation against analytical solutions
  - Testing of first to fourth-order systems

- **`examples/transfer_function.rs`**: Demonstration examples
  - First-order low-pass filter
  - Second-order critically damped system
  - Lead compensator with direct feedthrough
  - Pure integrator

- **`src/blocks/mod.rs`**: Module export

### State-Space Realization

The implementation uses **Observable Canonical Form** (OCF), matching scipy.signal.TransferFunction.to_ss():

For H(s) = (b_n s^n + ... + b_0) / (s^n + a_{n-1} s^{n-1} + ... + a_0):

```
A = [-a_{n-1}  -a_{n-2}  ...  -a_1  -a_0 ]
    [   1         0      ...   0     0   ]
    [   0         1      ...   0     0   ]
    [   ⋮         ⋮      ⋱     ⋮     ⋮   ]
    [   0         0      ...   1     0   ]

B = [1]
    [0]
    [⋮]
    [0]

C = [b_{n-1}  b_{n-2}  ...  b_1  b_0]  (after removing D component)

D = b_n  (if numerator degree equals denominator degree)
```

### Key Features

1. **Polynomial Convention**: Descending powers of s
   - `num = [b_n, b_{n-1}, ..., b_0]` represents `b_n*s^n + ... + b_0`
   - `den = [a_m, a_{m-1}, ..., a_0]` represents `a_m*s^m + ... + a_0`

2. **Automatic Normalization**: Leading denominator coefficient normalized to 1

3. **Direct Feedthrough**: Properly handles D matrix when num and den have equal degree

4. **Proper vs Improper**: Only supports strictly proper and proper transfer functions (numerator degree ≤ denominator degree)

5. **Integration**: Uses RK4 (4th-order Runge-Kutta) for numerical integration

## Test Coverage

### Test Suite (15 comprehensive tests)

1. **test_tf_01_first_order_lowpass**: H(s) = 1/(s+1)
   - Validates step response against analytical solution
   - Tests at multiple time constants (1τ, 2τ, 3τ, 5τ)

2. **test_tf_02_first_order_highpass**: H(s) = s/(s+1)
   - Tests system with direct feedthrough
   - Validates exponential decay response

3. **test_tf_03_second_order_underdamped**: ζ = 0.2, ω_n = 1
   - Tests oscillatory response
   - Validates damped oscillation formula

4. **test_tf_04_second_order_critically_damped**: H(s) = 1/(s+1)^2
   - Tests critical damping case
   - Validates no overshoot behavior

5. **test_tf_05_second_order_overdamped**: H(s) = 1/((s+1)(s+2))
   - Tests heavily damped response
   - Validates partial fraction decomposition

6. **test_tf_06_lead_compensator**: H(s) = (s+2)/(s+10)
   - Tests system with zero
   - Validates lead compensation behavior

7. **test_tf_07_lag_compensator**: H(s) = (s+1)/(s+0.1)
   - Tests large DC gain
   - Validates steady-state value

8. **test_tf_08_pure_integrator**: H(s) = 1/s
   - Tests marginally stable system
   - Validates linear ramp response

9. **test_tf_09_double_integrator**: H(s) = 1/s^2
   - Tests unstable system under step input
   - Validates quadratic response

10. **test_tf_10_notch_filter**: H(s) = (s^2+1)/(s^2+0.1s+1)
    - Tests complex numerator dynamics
    - Validates notch frequency behavior

11. **test_tf_11_third_order**: H(s) = 1/((s+1)(s+2)(s+3))
    - Tests higher-order system
    - Validates DC gain

12. **test_tf_12_fourth_order_butterworth**: 4th-order Butterworth
    - Tests maximum order system
    - Validates filter response

13. **test_tf_13_resonant_system**: ζ = 0.05 (very lightly damped)
    - Tests near-oscillatory behavior
    - Validates overshoot prediction

14. **test_tf_14_impulse_response**: Impulse approximation test
    - Tests transient response
    - Validates impulse-to-step relationship

15. **test_tf_15_dc_gain**: DC gain verification
    - Tests multiple system orders
    - Validates steady-state accuracy

## Comparison with PathSim

The implementation follows PathSim's approach:

```python
# PathSim (from pathsim/src/pathsim/blocks/lti.py)
class TransferFunctionNumDen(StateSpace):
    def __init__(self, Num=[1], Den=[1, 1]):
        self.Num, self.Den = Num, Den
        sp_SS = _TransferFunction(Num, Den).to_ss()
        super().__init__(sp_SS.A, sp_SS.B, sp_SS.C, sp_SS.D)
```

```rust
// RustSim equivalent
let tf = TransferFunction::<1>::new(&[1.0], &[1.0, 1.0]);
// Internally converts to state-space using same algorithm as scipy
```

### Verified Compatibility

The implementation produces identical state-space matrices and simulation results to scipy.signal.TransferFunction.to_ss(), ensuring compatibility with PathSim's behavior.

## Usage Examples

### Basic First-Order System

```rust
use rustsim::blocks::TransferFunction;
use rustsim::Block;

// H(s) = 1/(s+1)
let mut tf = TransferFunction::<1>::new(&[1.0], &[1.0, 1.0]);

// Apply step input
tf.set_input(0, 1.0);
tf.update(0.0);

// Simulate
let dt = 0.01;
for _ in 0..500 {
    tf.step(0.0, dt);
}

// Check output (should be near 0.993 after 5 time constants)
let output = tf.get_output(0);
```

### Second-Order System

```rust
// H(s) = 1/(s^2 + 2s + 1)
let mut tf = TransferFunction::<2>::new(&[1.0], &[1.0, 2.0, 1.0]);

// Rest of simulation same as above
```

### Lead Compensator with Direct Feedthrough

```rust
// H(s) = (s+2)/(s+10)
let mut tf = TransferFunction::<1>::new(&[1.0, 2.0], &[1.0, 10.0]);

// Check for direct feedthrough
assert!(tf.has_passthrough());
```

## Performance

- **Compilation**: Const generics allow compile-time size determination
- **Runtime**: Zero-overhead abstraction with inline functions
- **Integration**: RK4 provides good accuracy for typical control systems
- **Memory**: Fixed-size arrays (no heap allocations during simulation)

## Error Handling

The implementation panics on:
- Empty denominator
- Zero leading coefficient in denominator
- Improper transfer function (num degree > den degree)
- Order mismatch between denominator and const generic parameter

## Future Enhancements

Potential additions:
1. Support for MIMO transfer functions (multiple inputs/outputs)
2. Zeros-Poles-Gain (ZPG) representation
3. Frequency response analysis
4. Bode plot generation
5. Pole-Residue-Constant (PRC) form (like PathSim's TransferFunctionPRC)

## References

- Ogata, K. (2010). Modern Control Engineering (5th ed.). Section 5.6
- Chen, C.T. (1999). Linear System Theory and Design (3rd ed.). Section 5.5
- scipy.signal.TransferFunction documentation
- PathSim implementation: pathsim/src/pathsim/blocks/lti.py
