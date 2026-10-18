"""Analytical PID design by pole placement, checked in simulation.

The closed-loop polynomial is s^4 + 40 s^3 + a2 s^2 + a3 s + a4 with
a2 = -392.4 + N kd,  a3 = -15696 + N kp,  a4 = N ki   (N = plant gain, about -560.6).
Choosing the four closed-loop poles fixes the polynomial; the three gains follow by matching
coefficients. The s^3 coefficient (40 = R/L) cannot be chosen: the poles must sum to -40.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import signal
from plant import A, B, DEN, N_GAIN, char_poly_pid, hurwitz_pid

lines = []
def log(s=""):
    print(s); lines.append(s)

# Desired poles: -8 +/- 12j, -10, -14 (sum = -40)
poles = [-8 + 12j, -8 - 12j, -10, -14]
des = np.real(np.poly(poles))
assert abs(des[1] - 40) < 1e-9
kd = (des[2] - DEN[2]) / N_GAIN
kp = (des[3] - DEN[3]) / N_GAIN
ki = des[4] / N_GAIN
log("Desired polynomial: " + str(np.round(des, 1)))
log("Gains: Kp = %.2f, Kd = %.3f, Ki = %.2f" % (kp, kd, ki))
log("Closed-loop roots: " + str(np.round(np.roots(char_poly_pid(kp, kd, ki)), 3)))
log("Inside the Routh-Hurwitz region: %s" % hurwitz_pid(kp, kd, ki))

char = char_poly_pid(kp, kd, ki)
t = np.linspace(0, 3, 3001)

def metrics(y):
    os_ = max(0.0, (y.max() - 1) * 100)
    ts = t[np.where(np.abs(y - 1) > 0.02)[0][-1]] if (np.abs(y - 1) > 0.02).any() else 0.0
    return os_, ts

# Standard PID: u = (Kp + Kd s + Ki/s)(r - y)  -> closed loop has the controller zeros
_, y_pid = signal.step(signal.TransferFunction(N_GAIN * np.array([kd, kp, ki]), char), T=t)
# I-PD: integral on the error, Kp and Kd on the measurement -> same poles, no zeros
_, y_ipd = signal.step(signal.TransferFunction([N_GAIN * ki], char), T=t)
for name, y in (("Standard PID", y_pid), ("I-PD", y_ipd)):
    os_, ts = metrics(y)
    log("%-13s step response: overshoot %.1f %%, settling (2 %%) %.2f s" % (name, os_, ts))

# Time-domain simulation of the I-PD loop on the state-space plant, 1 mm step
Acl = np.zeros((4, 4))
Acl[:3, :3] = A + B @ np.array([[-kp, -kd, 0.0]])
Acl[:3, 3] = (B * ki).ravel()
Acl[3, 0] = -1.0
Bcl = np.array([[0], [0], [0], [1.0]])
sys = signal.StateSpace(Acl, Bcl, np.array([[1.0, 0, 0, 0]]), np.zeros((1, 1)))
r = 0.001 * np.ones_like(t)
_, y, x = signal.lsim(sys, r, t)
u = ki * x[:, 3] - kp * x[:, 0] - kd * x[:, 1]
log("1 mm step, I-PD: final %.6f mm, peak |voltage deviation| %.2f V (equilibrium voltage is 1.4 V; no saturation modelled)"
    % (y[-1] * 1000, np.abs(u).max()))

# A gain set outside the region
bad = (-20.0, -2.0, -10.0)
log("Outside the region %s: roots %s" % (bad, np.round(np.roots(char_poly_pid(*bad)), 2)))

plt.figure()
plt.plot(t, y_pid, label="standard PID"); plt.plot(t, y_ipd, label="I-PD (same poles)")
plt.axhline(1, color="gray", ls="--"); plt.xlabel("time (s)"); plt.ylabel("normalised position")
plt.title("Step response, designed poles -8 +/- 12j, -10, -14"); plt.legend(); plt.grid(True)
plt.savefig("results/step_response.png", dpi=150); plt.close()
open("results/design_summary.txt", "w").write("\n".join(lines) + "\n")
