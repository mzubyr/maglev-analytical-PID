"""Analytical PID design by pole placement, checked in simulation, with the control library.

The closed-loop polynomial is s^4 + 40 s^3 + a2 s^2 + a3 s + a4 with
a2 = -392.4 + N kd,  a3 = -15696 + N kp,  a4 = N ki   (N = plant gain, about -560.6).
Choosing the four closed-loop poles fixes the polynomial; the three gains follow by matching coefficients.

Run: python design_and_simulate.py      (a figure and a summary are written to results/)
"""
import numpy as np
import control as ct
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plant import (A, B, N_GAIN, poles_of, closed_loop_pid, closed_loop_ipd, routh_stable_pid, design_gains)

lines = []
def log(s=""):
    print(s); lines.append(s)

# Desired poles: -8 +/- 12j, -10, -14 (sum = -40)
kp, kd, ki, des = design_gains([-8 + 12j, -8 - 12j, -10, -14])
log("Desired polynomial: " + str(np.round(des, 1)))
log("Gains: Kp = %.2f, Kd = %.3f, Ki = %.2f" % (kp, kd, ki))
log("Closed-loop poles: " + str(np.round(poles_of(closed_loop_pid(kp, kd, ki)), 3)))
log("Inside the Routh-Hurwitz region: %s" % routh_stable_pid(kp, kd, ki))

t = np.linspace(0, 3, 3001)

def metrics(y):
    overshoot = max(0.0, (y.max() - 1) * 100)
    outside = np.where(np.abs(y - 1) > 0.02)[0]
    return overshoot, (t[outside[-1]] if len(outside) else 0.0)

plt.figure()
for name, system in (("Standard PID", closed_loop_pid(kp, kd, ki)), ("I-PD", closed_loop_ipd(kp, kd, ki))):
    resp = ct.step_response(system, T=t)
    y = np.asarray(resp.outputs).ravel()
    os_, ts = metrics(y)
    log("%-13s step response: overshoot %.1f %%, settling (2 %%) %.2f s" % (name, os_, ts))
    plt.plot(resp.time, y, label=name)

# Time-domain simulation of the I-PD loop on the state-space plant, 1 mm step.
# States [x1, x2, x3, xi]; u = ki*xi - kp*x1 - kd*x2; xi' = r - x1.
Acl = np.zeros((4, 4))
Acl[:3, :3] = A + B @ np.array([[-kp, -kd, 0.0]])
Acl[:3, 3] = (B * ki).ravel()
Acl[3, 0] = -1.0
loop = ct.ss(Acl, np.array([[0.0], [0.0], [0.0], [1.0]]), np.array([[1.0, 0, 0, 0]]), np.zeros((1, 1)))
sim = ct.forced_response(loop, T=t, U=0.001 * np.ones_like(t))
x = np.asarray(sim.states)
u = ki * x[3] - kp * x[0] - kd * x[1]
log("1 mm step, I-PD: final %.6f mm, peak |voltage deviation| %.2f V (equilibrium voltage is 1.4 V; no saturation modelled)"
    % (np.asarray(sim.outputs).ravel()[-1] * 1000, np.abs(u).max()))

bad = (-20.0, -2.0, -10.0)
log("Outside the region %s: poles %s" % (bad, np.round(poles_of(closed_loop_pid(*bad)), 2)))

plt.axhline(1, color="gray", ls="--"); plt.xlabel("time (s)"); plt.ylabel("normalised position")
plt.title("Step response, designed poles -8 +/- 12j, -10, -14"); plt.legend(); plt.grid(True)
plt.savefig("results/step_response.png", dpi=150); plt.close()
open("results/design_summary.txt", "w").write("\n".join(lines) + "\n")
