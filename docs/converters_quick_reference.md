# ADC/DAC Quick Reference

## Import

```rust
use rustsim::blocks::{ADC, DAC};
use rustsim::block::Block;
```

## ADC - Analog to Digital Converter

### Creation

```rust
// 4-bit ADC, span [-1, 1], period 1.0, delay 0.0
let adc = ADC::<4>::new([-1.0, 1.0], 1.0, 0.0);

// 8-bit ADC, span [0, 5], period 0.1, delay 0.05
let adc = ADC::<8>::new([0.0, 5.0], 0.1, 0.05);

// Default: 4-bit, [-1, 1], period 1.0, delay 0.0
let adc = ADC::<4>::default();
```

### Usage

```rust
let mut adc = ADC::<8>::new([0.0, 5.0], 1.0, 0.0);

// Set analog input
adc.set_input(0, 3.7);

// Sample (at t=0, first sample time)
adc.update(0.0);

// Read digital bits (LSB to MSB)
let lsb = adc.get_output(0);
let msb = adc.get_output(7);

// Read all bits
for i in 0..8 {
    println!("Bit {}: {}", i, adc.get_output(i));
}
```

### Ports

- **Inputs**: 1 (analog input)
- **Outputs**: N_BITS (digital bits, LSB at port 0)

### Behavior

```
Input Range: [span[0], span[1]]
Quantization: floor(scaled * 2^N_BITS)
Output: Binary code on N_BITS ports
Sampling: Periodic at times: tau, tau+T, tau+2T, ...
```

## DAC - Digital to Analog Converter

### Creation

```rust
// 4-bit DAC, span [-1, 1], period 1.0, delay 0.0
let dac = DAC::<4>::new([-1.0, 1.0], 1.0, 0.0);

// 8-bit DAC, span [0, 5], period 0.1, delay 0.05
let dac = DAC::<8>::new([0.0, 5.0], 0.1, 0.05);

// Default: 4-bit, [-1, 1], period 1.0, delay 0.0
let dac = DAC::<4>::default();
```

### Usage

```rust
let mut dac = DAC::<8>::new([0.0, 5.0], 1.0, 0.0);

// Set digital bits (LSB to MSB)
dac.set_input(0, 1.0); // LSB
dac.set_input(1, 0.0);
dac.set_input(2, 1.0);
// ... set remaining bits
dac.set_input(7, 1.0); // MSB

// Convert (at t=0, first conversion time)
dac.update(0.0);

// Read analog output
let analog_out = dac.get_output(0);
```

### Ports

- **Inputs**: N_BITS (digital bits, LSB at port 0)
- **Outputs**: 1 (analog output)

### Behavior

```
Digital Code: Σ(bit[i] * 2^i) for i in 0..N_BITS
Scaling: code / (2^N_BITS - 1)
Output: span[0] + (span[1] - span[0]) * scaling
Update: Periodic at times: tau, tau+T, tau+2T, ...
```

## Examples

### Example 1: Basic ADC

```rust
let mut adc = ADC::<4>::new([-1.0, 1.0], 1.0, 0.0);

adc.set_input(0, 0.0);  // Midrange
adc.update(0.0);

// Expected: code 8 = binary 1000 = [0,0,0,1] LSB to MSB
assert_eq!(adc.get_output(0), 0.0); // LSB
assert_eq!(adc.get_output(3), 1.0); // MSB
```

### Example 2: Basic DAC

```rust
let mut dac = DAC::<4>::new([-1.0, 1.0], 1.0, 0.0);

// Set code to 15 (all 1s)
for i in 0..4 {
    dac.set_input(i, 1.0);
}

dac.update(0.0);

// Expected: maximum of span
assert_eq!(dac.get_output(0), 1.0);
```

### Example 3: Round-Trip Conversion

```rust
let mut adc = ADC::<8>::new([0.0, 5.0], 1.0, 0.0);
let mut dac = DAC::<8>::new([0.0, 5.0], 1.0, 0.0);

// Original signal
let original = 3.14;

// ADC: Analog → Digital
adc.set_input(0, original);
adc.update(0.0);

// Transfer bits
for i in 0..8 {
    dac.set_input(i, adc.get_output(i));
}

// DAC: Digital → Analog
dac.update(0.0);
let reconstructed = dac.get_output(0);

// Quantization error should be small
let error = (original - reconstructed).abs();
let max_error = 5.0 / 255.0; // LSB for 8-bit
assert!(error <= max_error);
```

### Example 4: Periodic Sampling

```rust
let mut adc = ADC::<4>::new([-1.0, 1.0], 1.0, 0.0);

adc.set_input(0, 0.5);

// Sample at t=0
adc.update(0.0);
let code1: Vec<f64> = (0..4).map(|i| adc.get_output(i)).collect();

// Change input, but t=0.5 is not a sample time
adc.set_input(0, -0.5);
adc.update(0.5);
let code2: Vec<f64> = (0..4).map(|i| adc.get_output(i)).collect();

// Output should be unchanged (holding previous sample)
assert_eq!(code1, code2);

// Sample at t=1.0 (next sample time)
adc.update(1.0);
let code3: Vec<f64> = (0..4).map(|i| adc.get_output(i)).collect();

// Now it should be different
assert_ne!(code1, code3);
```

### Example 5: Different Bit Widths

```rust
// 1-bit (binary)
let adc1 = ADC::<1>::new([0.0, 1.0], 1.0, 0.0);
let dac1 = DAC::<1>::new([0.0, 1.0], 1.0, 0.0);

// 4-bit (common)
let adc4 = ADC::<4>::new([0.0, 1.0], 1.0, 0.0);
let dac4 = DAC::<4>::new([0.0, 1.0], 1.0, 0.0);

// 8-bit (typical)
let adc8 = ADC::<8>::new([0.0, 1.0], 1.0, 0.0);
let dac8 = DAC::<8>::new([0.0, 1.0], 1.0, 0.0);

// 12-bit (high precision)
let adc12 = ADC::<12>::new([0.0, 1.0], 1.0, 0.0);
let dac12 = DAC::<12>::new([0.0, 1.0], 1.0, 0.0);

// 16-bit (very high precision)
let adc16 = ADC::<16>::new([0.0, 1.0], 1.0, 0.0);
let dac16 = DAC::<16>::new([0.0, 1.0], 1.0, 0.0);
```

## Quantization Levels

| Bits | Levels | Step Size (0-1V) | ENOB (dB) |
|------|--------|------------------|-----------|
| 1    | 2      | 1.000            | 6.02      |
| 2    | 4      | 0.333            | 12.04     |
| 4    | 16     | 0.067            | 24.08     |
| 8    | 256    | 0.004            | 48.16     |
| 10   | 1024   | 0.001            | 60.21     |
| 12   | 4096   | 0.0002           | 72.25     |
| 16   | 65536  | 0.00002          | 96.33     |

## Common Use Cases

### Audio ADC/DAC (16-bit, ±1V)
```rust
let adc = ADC::<16>::new([-1.0, 1.0], 1.0/44100.0, 0.0); // 44.1kHz
let dac = DAC::<16>::new([-1.0, 1.0], 1.0/44100.0, 0.0);
```

### Control System (12-bit, 0-10V)
```rust
let adc = ADC::<12>::new([0.0, 10.0], 0.001, 0.0); // 1kHz
let dac = DAC::<12>::new([0.0, 10.0], 0.001, 0.0);
```

### Digital Communication (1-bit)
```rust
let adc = ADC::<1>::new([0.0, 1.0], 1e-6, 0.0); // 1MHz
let dac = DAC::<1>::new([0.0, 1.0], 1e-6, 0.0);
```

### Sensor Interface (10-bit, 0-5V)
```rust
let adc = ADC::<10>::new([0.0, 5.0], 0.01, 0.0); // 100Hz
let dac = DAC::<10>::new([0.0, 5.0], 0.01, 0.0);
```
