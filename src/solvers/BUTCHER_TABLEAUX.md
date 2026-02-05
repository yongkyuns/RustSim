# Butcher Tableaux Reference

This document provides the Butcher tableaux for all implemented Runge-Kutta methods.

## Notation

A Butcher tableau has the form:

```
c₁ | a₁₁  a₁₂  a₁₃  ...
c₂ | a₂₁  a₂₂  a₂₃  ...
c₃ | a₃₁  a₃₂  a₃₃  ...
...
---+--------------------
   | b₁   b₂   b₃   ...
```

For explicit methods, the matrix A is strictly lower triangular (aᵢⱼ = 0 for j ≥ i).

The update formulas are:

```
kᵢ = f(xₙ + dt * Σⱼ aᵢⱼ kⱼ, tₙ + cᵢ * dt)
xₙ₊₁ = xₙ + dt * Σᵢ bᵢ kᵢ
```

## Forward Euler

**Order**: 1, **Stages**: 1

```
0 |
---+---
  | 1
```

## SSPRK22 (Heun's Method)

**Order**: 2, **Stages**: 2, **SSP Coefficient**: C = 1

```
  0 |
  1 | 1
----+-------
    | 1/2  1/2
```

Equivalent form (Shu-Osher):
```
k₁ = f(xₙ)
k₂ = f(xₙ + dt * k₁)
xₙ₊₁ = xₙ + dt/2 * (k₁ + k₂)
```

## SSPRK33

**Order**: 3, **Stages**: 3, **SSP Coefficient**: C = 1

```
  0 |
  1 | 1
1/2 | 1/4  1/4
----+-------------
    | 1/6  1/6  2/3
```

Equivalent form (Shu-Osher):
```
x⁽¹⁾ = xₙ + dt * f(xₙ)
x⁽²⁾ = 3/4 * xₙ + 1/4 * x⁽¹⁾ + 1/4 * dt * f(x⁽¹⁾)
xₙ₊₁ = 1/3 * xₙ + 2/3 * x⁽²⁾ + 2/3 * dt * f(x⁽²⁾)
```

## RK4 (Classical Runge-Kutta)

**Order**: 4, **Stages**: 4

```
  0 |
1/2 | 1/2
1/2 | 0    1/2
  1 | 0    0    1
----+-------------------
    | 1/6  1/3  1/3  1/6
```

Note: The weights for stages 2 and 3 are 2/6 = 1/3.

## RKDP54 (Dormand-Prince 5(4))

**Order**: 5 (propagating) / 4 (embedded), **Stages**: 7

```
   0 |
 1/5 | 1/5
3/10 | 3/40       9/40
 4/5 | 44/45     -56/15      32/9
 8/9 | 19372/6561 -25360/2187 64448/6561 -212/729
   1 | 9017/3168  -355/33     46732/5247  49/176  -5103/18656
   1 | 35/384     0           500/1113    125/192 -2187/6784  11/84
-----+------------------------------------------------------------------------
  b  | 35/384     0           500/1113    125/192 -2187/6784  11/84    0
  b* | 5179/57600 0           7571/16695  393/640 -92097/339200 187/2100 1/40
```

Where:
- `b` = 5th order solution (propagating)
- `b*` = 4th order solution (embedded)

**Error estimate coefficients** (b - b*):

```
TR = [71/57600, 0, -71/16695, 71/1920, -17253/339200, 22/525, -1/40]
```

The local truncation error is estimated as:

```
LTE ≈ dt * Σᵢ TRᵢ * kᵢ
```

### Timestep Control

The error controller uses a PI-type controller:

```
error_norm = max(|LTE| / scale)
where scale = atol + rtol * |x|

success = (error_norm ≤ 1)

timestep_scale = β / error_norm^(1/(p+1))
where β = 0.9 (safety factor)
      p = min(5, 4) = 4 (min of propagating and embedded orders)

timestep_scale is clamped to [0.1, 10.0]
```

Default tolerances:
- Absolute: `1e-8`
- Relative: `1e-4`

## Properties Summary

| Method   | Order | Stages | Adaptive | FSAL | SSP | Use Case                    |
|----------|-------|--------|----------|------|-----|----------------------------|
| Euler    | 1     | 1      | No       | Yes  | Yes | Testing only               |
| SSPRK22  | 2     | 2      | No       | Yes  | Yes | TVD schemes, C=1           |
| SSPRK33  | 3     | 3      | No       | No   | Yes | TVD schemes, C=1           |
| RK4      | 4     | 4      | No       | No   | No  | General fixed-step         |
| RKDP54   | 5/4   | 7      | Yes      | Yes  | No  | General adaptive (default) |

**FSAL**: First Same As Last (last stage of step n can be reused as first stage of step n+1)
**SSP**: Strong Stability Preserving

## References

### RK4
- Kutta, W. (1901). "Beitrag zur näherungsweisen Integration totaler Differentialgleichungen".
  Zeitschrift für Mathematik und Physik, 46, 435-453.

### RKDP54
- Dormand, J. R., & Prince, P. J. (1980). "A family of embedded Runge-Kutta formulae".
  Journal of Computational and Applied Mathematics, 6(1), 19-26.
  DOI: 10.1016/0771-050X(80)90013-3

- Shampine, L. F., & Reichelt, M. W. (1997). "The MATLAB ODE Suite".
  SIAM Journal on Scientific Computing, 18(1), 1-22.
  DOI: 10.1137/S1064827594276424

### SSP Methods
- Shu, C.-W., & Osher, S. (1988). "Efficient implementation of essentially non-oscillatory
  shock-capturing schemes". Journal of Computational Physics, 77(2), 439-471.
  DOI: 10.1016/0021-9991(88)90177-5

- Gottlieb, S., Shu, C.-W., & Tadmor, E. (2001). "Strong stability-preserving high-order
  time discretization methods". SIAM Review, 43(1), 89-112.
  DOI: 10.1137/S003614450036757X

- Gottlieb, S., Ketcheson, D. I., & Shu, C.-W. (2011). "Strong Stability Preserving
  Runge-Kutta and Multistep Time Discretizations". World Scientific.
  DOI: 10.1142/7498
