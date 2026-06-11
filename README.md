# ukf-motor-cpp-STM32H755

Adaptive Unscented Kalman Filter for DC motor state estimation — 4-state augmented filter implemented in C++ targeting the **STM32H755ZI** at 480 MHz.

This is the embedded C++ port of the MATLAB prototype. For full filter design details, tuning methodology, experimental data, and validation plots refer to the companion repo: **[ukf-motor-sysid](https://github.com/ramdanebrahim9/ukf-motor-sysid)** — also available in my repos.

---

## Target

| | |
|---|---|
| **MCU** | STM32H755ZI — Cortex-M7 @ 480 MHz |
| **Sampling rate** | 20 kHz (Ts = 50 µs) |
| **Math library** | Eigen — fixed-size types only |
| **Memory** | No heap, no dynamic allocation, no STL containers |
| **FPU** | Full double-precision, all operations in registers |

This is written for a microcontroller. That means: **no `new`, no `malloc`, no `std::vector`, no exceptions**. Every matrix is a fixed-size Eigen type resolved entirely at compile time. The filter runs in a hard real-time ISR context with deterministic timing enforced by TIM1 + DWT cycle counter on the CM7 core.

---

## What this estimates

| State | Symbol | Unit | Description |
|-------|--------|------|-------------|
| Position | `θ` | rad | Output shaft angle (gear ratio N = 13.7335) |
| Velocity | `ω` | rad/s | Angular velocity from 64 PPR encoder |
| Current | `i` | A | Armature current from INA240 sense amplifier |
| Disturbance | `d` | A | Lumped model error — friction residual + unmodeled dynamics |

---

## Repo structure

```
ukf-motor-cpp-STM32H755/
│
├── UKF.hpp / UKF.cpp              # top-level filter class — single step() call per loop
│
├── Models/
│   ├── UKF_motor_ode.cpp/.hpp     # 2-state motor ODE (ω, i) — Stribeck friction
│   ├── UKF_motor_ode_Aug.cpp/.hpp # 4-state augmented ODE (θ, ω, i, d)
│   ├── UKF_rk4_step.cpp/.hpp      # RK4 integrator — 2-state open-loop
│   └── UKF_rk4_step_Aug.cpp/.hpp  # RK4 integrator — 4-state augmented
│
└── utilities/
    ├── UKF_predict.cpp/.hpp       # UKF predict step — sigma point propagation
    ├── UKF_correct.cpp/.hpp       # UKF correct step — Kalman gain + update
    ├── UKF_sigma_points.cpp/.hpp  # 2n+1 sigma points via Cholesky
    ├── UKF_weights.cpp/.hpp       # Wm / Wc weights (alpha, beta, kappa)
    ├── UKF_unscented_transform.hpp# weighted mean + covariance
    ├── UKF_get_Q.cpp/.hpp         # adaptive process noise — voltage-scheduled
    ├── UKF_get_R.cpp/.hpp         # adaptive measurement noise — ω + tick scheduled
    └── UKF_utils.hpp              # interp1, clamp, helpers
```

---

## Filter design

### Process model — augmented ODE (`UKF_motor_ode_Aug.cpp`)

```
dθ/dt = ω
dω/dt = (Ke/J)·i_rk4 − friction(ω)        ← i injected from open-loop RK4
di/dt = (1/L)·(V − R·i − Ke·ω_rk4 + d)    ← ω injected from open-loop RK4
dd/dt = −d / τ_d                            τ_d = 0.1 s
```

Friction model: Stribeck + viscous
```
f(ω) = tanh(ω/ε)·[tc/J + (ts/J − tc/J)·exp(−(ω/ωs)²)] + (B/J)·ω
```

RK4 cross-injection decouples the electrical and mechanical states — prevents the sigma point instability that arises from a fully coupled 4-state propagation at 20 kHz.

### Motor parameters

| Parameter | Value | Unit |
|-----------|-------|------|
| R | 1.438907 | Ω |
| L | 0.415 | mH |
| Ke | 0.137961 | V·s/rad |
| J | 1.263×10⁻⁴ | kg·m² |
| Gear ratio N | 13.7335 | — |
| Encoder | 64 | PPR |

### Adaptive Q

Process noise on `i` is voltage-scheduled via a compile-time lookup table. `θ` and `ω` are fixed. `d` uses constant `q_d = 0.01`.

### Adaptive R

Dual-scheduled measurement noise:
- **`r_i`** — ω-indexed lookup table (ADC noise characterized vs operating point)
- **`r_θ`, `r_ω`** — tick-event scheduled:
  - On tick fired: `r_θ = 1e-6`, `r_ω = k1(V)` — trust the encoder
  - Between ticks: `r_θ = tick_rad`, `r_ω = k2(V)` — coast with high uncertainty

### Disturbance clipping

```cpp
drift   = mean_i_iir − i_rk4;
d_clip  = (|drift| > thresh(V)) ? |drift| : 0.0;
d       = clamp(d, −d_clip, +d_clip);
```

IIR mean current filtered at 200 Hz (`α = 0.9391`, computed from `exp(−2π·200·50e-6)`).

---

## Usage

One call per control loop tick:

```cpp
UKF_InitParams ic;
ic.theta_0 = ticks_to_rad(encoder_ticks);
ic.omega_0 = 0.0f;
ic.i_0     = i_measured;
ic.d_0     = 0.0f;

UKF::init(ic);

// inside ISR or control loop at 20 kHz:
const UKF_State& est = UKF::step(V_cmd, i_measured, encoder_ticks, omega_dirty);

double theta = est.x_est(0);
double omega = est.x_est(1);
double i     = est.x_est(2);
double d     = est.x_est(3);
```

---

## Dependencies

- **[Eigen](https://eigen.tuxfamily.org/)** — header-only, fixed-size types (`Vector2d`, `Vector4d`, `Matrix4d`, `Matrix3d`). No dynamic Eigen usage anywhere.
- ARM GCC toolchain with `-O2` or `-O3` and FPU flags: `-mfpu=fpv5-d16 -mfloat-abi=hard`

---

## Author

**Ibrahim Ramdane** — M1 VISTA (Vision, Signal, Trajectographie et Automatique), Université de Toulon  
Control systems · State estimation · Embedded real-time implementation
