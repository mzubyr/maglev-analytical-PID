"""Model, transfer function and Routh-Hurwitz stability region (PD and PID), with the control library.

Run: python stability_analysis.py      (figures and a summary are written to results/)
"""
import numpy as np
import control as ct
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plant import (A, B, plant, N_GAIN, DEN, poles_of, closed_loop_pid, closed_loop_pd,
                   routh_stable_pd, routh_stable_pid, poles_stable, i_max)

N_TRIALS = 5000          # random gain sets for the check; you can raise this, for example to 20000

lines = []
def log(s=""):
    print(s); lines.append(s)

log("Plant: G(s) = %.1f / (s^3 + %.0f s^2 %.1f s %.0f)" % (N_GAIN, DEN[1], DEN[2], DEN[3]))
log("Open-loop poles (eigenvalues of A): " + str(np.round(poles_of(plant), 3)))
ctrb = ct.ctrb(A, B)
log("Controllability: rank %d of 3, det = %.3e" % (np.linalg.matrix_rank(ctrb), np.linalg.det(ctrb)))

# Hand-derived stability conditions (p = -Kp, d = -Kd, i = -Ki; all positive)
log("\nPD  (cubic):  d > %.2f,  p > %.1f,  p < 40 d" % (-DEN[2] / -N_GAIN, -DEN[3] / -N_GAIN))
log("PID (quartic): the same, and  i < (40 d - p)(%.1f p - %.0f) / 1600" % (-N_GAIN, -DEN[3]))

# Compare the conditions with the closed-loop poles computed by the library
rng = np.random.default_rng(0)
agree_pd = agree_pid = 0
for _ in range(N_TRIALS):
    kp, kd, ki = -rng.uniform(0, 150), -rng.uniform(0, 6), -rng.uniform(0, 300)
    agree_pd += routh_stable_pd(kp, kd) == poles_stable(closed_loop_pd(kp, kd))
    agree_pid += routh_stable_pid(kp, kd, ki) == poles_stable(closed_loop_pid(kp, kd, ki))
log("\nRouth-Hurwitz vs closed-loop poles on %d random gain sets: PD agree %d, PID agree %d"
    % (N_TRIALS, agree_pd, agree_pid))

# Figure 1: PD stable region in the (d, p) plane, from the closed-loop poles
d_vals = np.linspace(0.05, 6, 80); p_vals = np.linspace(1, 250, 80)
region = np.array([[poles_stable(closed_loop_pd(-p, -d)) for d in d_vals] for p in p_vals], dtype=float)
plt.figure()
plt.contourf(d_vals, p_vals, region, levels=[-0.5, 0.5, 1.5], colors=["white", "lightgreen"])
plt.plot(d_vals, 40 * d_vals, "b--", label="p = 40 d")
plt.axhline(-DEN[3] / -N_GAIN, color="r", ls="--", label="p = 28.0")
plt.axvline(-DEN[2] / -N_GAIN, color="m", ls="--", label="d = 0.70")
plt.xlabel("d = -Kd"); plt.ylabel("p = -Kp"); plt.title("PD stability region (green)")
plt.legend(); plt.grid(True); plt.savefig("results/pd_region.png", dpi=150); plt.close()

# Figure 2: largest allowed i for d = 2
p_line = np.linspace(28, 79, 200)
plt.figure()
plt.plot(p_line, [i_max(p, 2.0) for p in p_line])
plt.xlabel("p = -Kp"); plt.ylabel("largest i = -Ki"); plt.title("PID: integral gain limit for d = 2")
plt.grid(True); plt.savefig("results/pid_i_limit.png", dpi=150); plt.close()

log("\nExamples (Kp, Kd, Ki):")
for g in [(-60, -2, -50), (-60, -2, -300), (-20, -2, -10), (20, 1, 5)]:
    log("  %s  Routh: %-5s  poles: %s" % (g, routh_stable_pid(*g), poles_stable(closed_loop_pid(*g))))
open("results/stability_summary.txt", "w").write("\n".join(lines) + "\n")
