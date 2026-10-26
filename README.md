# Maglev: analytical PID design (Routh-Hurwitz and pole placement)

An analytical design of a PID controller for a simplified magnetic-levitation model, checked
numerically. This project uses its own simplified model (constant coil inductance, no
friction) with its own parameters. It is separate from my maglev MPC project
(`maglev-mpc`), which uses a different plant, so the numbers here and there should not be
compared directly.

## Contents

| File | Purpose |
|---|---|
| `derivation.md` | Typed step-by-step derivation (model, linearisation, G(s), Routh-Hurwitz, pole placement) |
| `plant.py` | Model as a state-space system (`ct.ss`), transfer function (`ct.ss2tf`), PID, PD and I-PD closed loops (`ct.feedback`), Routh-Hurwitz checks, pole placement |
| `stability_analysis.py` | Prints G(s) and the stability conditions; compares Routh-Hurwitz with the closed-loop poles on thousands of random gain sets (5000 by default, see `N_TRIALS`); draws the PD region |
| `design_and_simulate.py` | Pole-placement design, step responses (standard PID vs I-PD), time-domain simulation |
| `results/` | Figures and text summaries produced by the scripts |
| `docs/` | Scanned hand-written derivation pages |
| `tests/` | Automatic checks (model numbers, transfer function, Routh-Hurwitz against roots, designed poles) |

## Model and results

Parameters: m = 0.1 kg, z0 = 0.05 m, i0 = 0.7 A, R = 2 ohm, L = 0.05 H, g = 9.81 m/s^2.
Linearised model: A = [[0, 1, 0], [392.4, 0, -28.03], [0, 0, -40]], B = [0, 0, 20]^T.
Open-loop eigenvalues +19.81, -19.81, -40 (unstable); controllable (rank 3).

Transfer function: G(s) = -560.6 / (s^3 + 40 s^2 - 392.4 s - 15696). The plant gain is
negative, so the stabilising gains are negative.

With p = -Kp, d = -Kd, i = -Ki:

- PD: d > 0.70, p > 28.0, p < 40 d (40 = R/L)
- PID: the same, and i < (560.6 p - 15696)(40 d - p) / 1600

The conditions agree with the closed-loop poles computed by the library on 5000 random gain sets (PD and PID); the number is set by `N_TRIALS` in `stability_analysis.py`.

Design by pole placement (poles -8 +/- 12j, -10, -14): **Kp = -40.90, Kd = -2.006, Ki = -51.95**.

| Structure | Overshoot | Settling (2 %) |
|---|---|---|
| Standard PID | 268 % | 0.80 s |
| I-PD (same poles, no zeros) | 0 % | 0.47 s |

For a 1 mm step with the I-PD controller the peak voltage deviation is 0.03 V (equilibrium
voltage 1.4 V).

Figures: `results/pd_region.png`, `results/pid_i_limit.png`, `results/step_response.png`.

## Requirements

The Python Control Systems Library (`control`), NumPy, SciPy and Matplotlib (`pip install -r requirements.txt`).
The model, transfer function, closed loops, poles and simulations all use the control library.

## Run

```
pip install -r requirements.txt
python stability_analysis.py
python design_and_simulate.py
python -m pytest tests        # optional
```

## Limitations

- Linear model valid near the equilibrium; no voltage saturation, noise or disturbances.
- Constant inductance and no friction; the parameters are my own choice.
- The pure derivative on a measured signal needs filtering in practice.
- Stability conditions say nothing about the recoverable region under a voltage limit.

## Author

Muhammad Ibrahim Zubair, Kano, Nigeria.
