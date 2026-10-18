"""Simplified maglev model for the analytical PID design (constant coil inductance, no friction).

States: x1 = gap z (m, measured downward from the magnet), x2 = velocity (m/s), x3 = coil current (A).
Input:  u = coil voltage deviation (V).  Everything is a deviation from the equilibrium.

Parameters (own choice):
    m = 0.1 kg, z0 = 0.05 m, i0 = 0.7 A, R = 2 ohm, L = 0.05 H, g = 9.81 m/s^2
    k = m g z0^2 / i0^2        (force constant, F = k i^2 / z^2)
    u0 = R i0 = 1.4 V          (equilibrium voltage)
"""
import numpy as np
from scipy import signal

M, G_ACC, Z0, I0, R, L = 0.1, 9.81, 0.05, 0.7, 2.0, 0.05
K = M * G_ACC * Z0**2 / I0**2
U0 = R * I0

A = np.array([[0.0, 1.0, 0.0],
              [2 * K * I0**2 / (M * Z0**3), 0.0, -2 * K * I0 / (M * Z0**2)],
              [0.0, 0.0, -R / L]])
B = np.array([[0.0], [0.0], [1.0 / L]])
C = np.array([[1.0, 0.0, 0.0]])
D = np.zeros((1, 1))


def transfer_function():
    """G(s) = C (sI - A)^-1 B = N / DEN with a constant numerator N."""
    num, den = signal.ss2tf(A, B, C, D)
    return float(num[0][-1]), den


N_GAIN, DEN = transfer_function()      # N_GAIN about -560.6, DEN = s^3 + 40 s^2 - 392.4 s - 15696


def char_poly_pd(kp, kd):
    """den(s) + N (kd s + kp): cubic."""
    p = DEN.copy()
    p[2] += N_GAIN * kd
    p[3] += N_GAIN * kp
    return p


def char_poly_pid(kp, kd, ki):
    """s den(s) + N (kd s^2 + kp s + ki): quartic."""
    p = np.concatenate([DEN, [0.0]])
    p[2] += N_GAIN * kd
    p[3] += N_GAIN * kp
    p[4] += N_GAIN * ki
    return p


def hurwitz_pd(kp, kd):
    """Routh-Hurwitz for the cubic s^3 + b2 s^2 + b1 s + b0."""
    _, b2, b1, b0 = char_poly_pd(kp, kd)
    return b2 > 0 and b1 > 0 and b0 > 0 and b2 * b1 > b0


def hurwitz_pid(kp, kd, ki):
    """Routh-Hurwitz for the quartic s^4 + a1 s^3 + a2 s^2 + a3 s + a4 (a0 = 1)."""
    _, a1, a2, a3, a4 = char_poly_pid(kp, kd, ki)
    return (a1 > 0 and a2 > 0 and a3 > 0 and a4 > 0
            and a1 * a2 - a3 > 0 and a3 * (a1 * a2 - a3) - a1**2 * a4 > 0)


def is_stable_roots(poly):
    return bool(np.all(np.roots(poly).real < 0))


def i_max(p, d):
    """Largest i = -Ki allowed by the last Routh condition, for p = -Kp, d = -Kd (p, d > 0)."""
    a1, a2, a3 = 40.0, -N_GAIN * d + DEN[2], -N_GAIN * p + DEN[3]
    return a3 * (a1 * a2 - a3) / (a1**2 * -N_GAIN)
