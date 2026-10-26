"""Simplified maglev model for the analytical PID design, built with the Python Control Systems Library.

States: x1 = gap z (m, measured downward from the magnet), x2 = velocity (m/s), x3 = coil current (A).
Input:  u = coil voltage deviation (V).  Everything is a deviation from the equilibrium.

Parameters (own choice):
    m = 0.1 kg, z0 = 0.05 m, i0 = 0.7 A, R = 2 ohm, L = 0.05 H, g = 9.81 m/s^2
    k = m g z0^2 / i0^2        (force constant, F = k i^2 / z^2)
    u0 = R i0 = 1.4 V          (equilibrium voltage)
"""
import numpy as np
import control as ct


def poles_of(system):
    """Poles of a system. The function name differs between versions of the library."""
    try:
        return np.asarray(ct.poles(system))
    except AttributeError:
        return np.asarray(ct.pole(system))


M, G_ACC, Z0, I0, R, L = 0.1, 9.81, 0.05, 0.7, 2.0, 0.05
K = M * G_ACC * Z0**2 / I0**2
U0 = R * I0

A = np.array([[0.0, 1.0, 0.0],
              [2 * K * I0**2 / (M * Z0**3), 0.0, -2 * K * I0 / (M * Z0**2)],
              [0.0, 0.0, -R / L]])
B = np.array([[0.0], [0.0], [1.0 / L]])
C = np.array([[1.0, 0.0, 0.0]])
D = np.zeros((1, 1))

plant = ct.ss(A, B, C, D)                      # state-space model
G = ct.ss2tf(plant)                            # G(s) = C (sI - A)^-1 B
_den = np.asarray(G.den[0][0], float)
DEN = _den / _den[0]                           # s^3 + 40 s^2 - 392.4 s - 15696
N_GAIN = float(np.asarray(G.num[0][0], float)[-1]) / _den[0]   # constant numerator, about -560.6


# ---------------------------------------------------------------------------------------------
# Controllers and closed loops
# ---------------------------------------------------------------------------------------------
def pid(kp, kd, ki):
    """K(s) = Kp + Kd s + Ki / s = (Kd s^2 + Kp s + Ki) / s"""
    return ct.tf([kd, kp, ki], [1, 0])


def pd(kp, kd):
    """K(s) = Kp + Kd s"""
    return ct.tf([kd, kp], [1])


def closed_loop_pid(kp, kd, ki):
    """Unity feedback with the PID in the forward path: T = G K / (1 + G K)."""
    return ct.feedback(G * pid(kp, kd, ki), 1)


def closed_loop_pd(kp, kd):
    return ct.feedback(G * pd(kp, kd), 1)


def closed_loop_ipd(kp, kd, ki):
    """I-PD: Kp and Kd act on the measurement, only Ki acts on the error. Same poles as the PID, no zeros.

    ct.feedback(G, K) = G / (1 + G K); multiplying by Ki / s gives the reference-to-output transfer function.
    """
    return ct.feedback(G, pid(kp, kd, ki)) * ct.tf([ki], [1, 0])


# ---------------------------------------------------------------------------------------------
# Hand-derived Routh-Hurwitz conditions (p = -Kp, d = -Kd, i = -Ki are positive for stabilising gains)
# ---------------------------------------------------------------------------------------------
def routh_stable_pd(kp, kd):
    """Cubic s^3 + b2 s^2 + b1 s + b0."""
    b2, b1, b0 = DEN[1], DEN[2] + N_GAIN * kd, DEN[3] + N_GAIN * kp
    return b2 > 0 and b1 > 0 and b0 > 0 and b2 * b1 > b0


def routh_stable_pid(kp, kd, ki):
    """Quartic s^4 + a1 s^3 + a2 s^2 + a3 s + a4 (a0 = 1)."""
    a1, a2, a3, a4 = DEN[1], DEN[2] + N_GAIN * kd, DEN[3] + N_GAIN * kp, N_GAIN * ki
    return (a1 > 0 and a2 > 0 and a3 > 0 and a4 > 0
            and a1 * a2 - a3 > 0 and a3 * (a1 * a2 - a3) - a1**2 * a4 > 0)


def poles_stable(system):
    return bool(np.all(poles_of(system).real < 0))


def i_max(p, d):
    """Largest i = -Ki allowed by the last Routh condition, for p = -Kp and d = -Kd."""
    a1, a2, a3 = DEN[1], -N_GAIN * d + DEN[2], -N_GAIN * p + DEN[3]
    return a3 * (a1 * a2 - a3) / (a1**2 * -N_GAIN)


# ---------------------------------------------------------------------------------------------
# Pole placement
# ---------------------------------------------------------------------------------------------
def design_gains(desired_poles):
    """Match s^4 + 40 s^3 + a2 s^2 + a3 s + a4 to the polynomial of the desired poles.

    The s^3 coefficient (40 = R/L) is fixed by the plant, so the four poles must sum to -40.
    Returns (kp, kd, ki) and the desired polynomial.
    """
    des = np.real(np.poly(desired_poles))
    if abs(des[1] - DEN[1]) > 1e-9:
        raise ValueError("the poles must sum to %.1f" % -DEN[1])
    kd = (des[2] - DEN[2]) / N_GAIN
    kp = (des[3] - DEN[3]) / N_GAIN
    ki = des[4] / N_GAIN
    return kp, kd, ki, des
