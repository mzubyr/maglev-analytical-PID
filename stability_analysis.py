"""Transfer function, closed-loop polynomial and Routh-Hurwitz stability region (PD and PID)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plant import (A, B, N_GAIN, DEN, char_poly_pd, char_poly_pid, hurwitz_pd, hurwitz_pid,
                   is_stable_roots, i_max)

lines = []
def log(s=""):
    print(s); lines.append(s)

log("Plant: G(s) = %.1f / (s^3 + %.0f s^2 %.1f s %.0f)" % (N_GAIN, DEN[1], DEN[2], DEN[3]))
log("Open-loop eigenvalues of A: " + str(np.round(np.linalg.eigvals(A), 3)))
ctrb = np.hstack([np.linalg.matrix_power(A, k) @ B for k in range(3)])
log("Controllability: rank %d of 3, det = %.3e" % (np.linalg.matrix_rank(ctrb), np.linalg.det(ctrb)))

# Hand-derived stability conditions (p = -Kp, d = -Kd, i = -Ki; all positive)
log("\nPD  (cubic):  d > %.2f,  p > %.1f,  p < 40 d" % (-DEN[2] / -N_GAIN, -DEN[3] / -N_GAIN))
log("PID (quartic): the same, and  i < (40 d - p)(%.1f p - %.0f) / 1600" % (-N_GAIN, -DEN[3]))

# Check Routh-Hurwitz against the roots on random gains
rng = np.random.default_rng(0)
agree_pd = agree_pid = 0
n = 20000
for _ in range(n):
    kp, kd, ki = -rng.uniform(0, 150), -rng.uniform(0, 6), -rng.uniform(0, 300)
    agree_pd += hurwitz_pd(kp, kd) == is_stable_roots(char_poly_pd(kp, kd))
    agree_pid += hurwitz_pid(kp, kd, ki) == is_stable_roots(char_poly_pid(kp, kd, ki))
log("\nRouth-Hurwitz vs numerical roots on %d random gain sets: PD agree %d, PID agree %d" % (n, agree_pd, agree_pid))

# Figure 1: PD stable region in the (d, p) plane
d_vals = np.linspace(0, 6, 241); p_vals = np.linspace(0, 250, 251)
Z = np.array([[hurwitz_pd(-p, -d) for d in d_vals] for p in p_vals])
plt.figure()
plt.contourf(d_vals, p_vals, Z.astype(float), levels=[-0.5, 0.5, 1.5], colors=["white", "lightgreen"])
dd = np.linspace(0, 6, 100)
plt.plot(dd, 40 * dd, "b--", label="p = 40 d")
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

# A few explicit checks
log("\nExamples (Kp, Kd, Ki):")
for g in [(-60, -2, -50), (-60, -2, -300), (-20, -2, -10), (20, 1, 5)]:
    log("  %s  Routh: %-5s  roots: %s" % (g, hurwitz_pid(*g), is_stable_roots(char_poly_pid(*g))))
open("results/stability_summary.txt", "w").write("\n".join(lines) + "\n")
