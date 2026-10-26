# Maglev: analytical PID design, worked derivation

Typed version of the hand-written derivation (see `docs/` for the scanned pages). Notation follows the
hand-written pages: closed-loop polynomial a0 s^4 + a1 s^3 + a2 s^2 + a3 s + a4 = 0.

## 1. Physical model

A steel ball of mass m is held below an electromagnet. z is the gap, measured downward.

- Coil: u = L di/dt + R i
- Mechanics: m z'' = m g - F, with F = k i^2 / z^2

Parameters: m = 0.1 kg, g = 9.81 m/s^2, z0 = 0.05 m, i0 = 0.7 A, R = 2 ohm, L = 0.05 H.
Assumptions: no disturbance (d = 0), constant inductance L, no friction.

## 2. State space

x = [x1, x2, x3] = [z, z', i], input u (voltage):

```
x1' = x2
x2' = g - k x3^2 / (m x1^2)
x3' = -(R/L) x3 + u / L
```

## 3. Equilibrium

With x2 = 0 and all derivatives zero:

- i0 = z0 sqrt(m g / k), so k = m g z0^2 / i0^2 = 0.1 x 9.81 x 0.05^2 / 0.7^2 = 0.005005
- u0 = R i0 = 2 x 0.7 = 1.4 V

## 4. Linearisation

A = df/dx and B = df/du at (x0, u0):

```
A = [ 0                    1    0
      2 k i0^2/(m z0^3)    0    -2 k i0/(m z0^2)
      0                    0    -R/L ]          B = [0, 0, 1/L]^T          C = [1 0 0]
```

Since k i0^2 / (m z0^2) = g: A(2,1) = 2g/z0 = 392.4 and A(2,3) = -2g/i0 = -28.03. With R/L = 40 and 1/L = 20:

```
A = [ 0     1    0          B = [ 0
      392.4 0   -28.03            0
      0     0   -40 ]            20 ]
```

## 5. Stability and controllability

det(sI - A) = (s + 40)(s^2 - 392.4). Eigenvalues +19.81, -19.81, -40: one is positive, so the open loop is unstable.

AB = [0, -560.6, -800]^T, A^2 B = [-560.6, 22423, 32000]^T and det [B AB A^2B] = -6.28e6, non-zero: fully controllable.

## 6. Transfer function

G(s) = C (sI - A)^-1 B. Because B has only a third entry and C picks the first row, only the cofactor
C31 = (-1)^(3+1) det[[-1, 0], [s, 28.03]] = -28.03 is needed:

```
G(s) = (-28.03 x 20) / (s^3 + 40 s^2 - 392.4 s - 15696) = -560.6 / (s^3 + 40 s^2 - 392.4 s - 15696)
```

The plant gain is negative: more voltage means more current, which pulls the ball toward the magnet and
reduces the gap. The stabilising controller gains are therefore negative.

## 7. Closed loop

E = z_ref - z, U = K E, z = G U, so z/z_ref = G K / (1 + G K), and 1 + G K = 0 gives the poles.
With K = Kp + Kd s + Ki/s:

```
s (s^3 + 40 s^2 - 392.4 s - 15696) + N (Kd s^2 + Kp s + Ki) = 0,      N = -560.6
s^4 + 40 s^3 + (-392.4 - 560.6 Kd) s^2 + (-15696 - 560.6 Kp) s - 560.6 Ki = 0
```

```
a0 = 1     a1 = 40     a2 = -392.4 - 560.6 Kd     a3 = -15696 - 560.6 Kp     a4 = -560.6 Ki
```

## 8. Routh array and stability conditions

```
s^4 | 1                          a2     a4
s^3 | 40                         a3     0
s^2 | b1 = (40 a2 - a3)/40       a4
s^1 | c1 = (b1 a3 - 40 a4)/b1    0
s^0 | a4
```

Stable if every number in the first column is positive. Let p = -Kp, d = -Kd, i = -Ki (all positive):

- a2 > 0: 560.6 d - 392.4 > 0, so d > 0.70
- a3 > 0: 560.6 p - 15696 > 0, so p > 27.998
- a4 > 0: i > 0
- b1 > 0: 40 a2 > a3, i.e. 22424 d - 15696 > 560.6 p - 15696, so p < 40 d (40 = R/L)
- c1 > 0: b1 a3 > 40 a4, i.e. i < (40 d - p)(560.6 p - 15696) / 1600

For a PD controller the polynomial is the cubic s^3 + 40 s^2 + a2 s + a3 and the conditions are the same
except that the last one is dropped.

## 9. Design by pole placement

The poles must sum to -a1 = -40. Chosen poles: -8 +/- 12j, -10, -14:

```
(s^2 + 16 s + 208)(s^2 + 24 s + 140) = s^4 + 40 s^3 + 732 s^2 + 7232 s + 29120
```

Matching coefficients:

```
-392.4 - 560.6 Kd = 732     ->  Kd = -(732 + 392.4) / 560.6   = -2.006
-15696 - 560.6 Kp = 7232    ->  Kp = -(7232 + 15696) / 560.6  = -40.90
-560.6 Ki = 29120           ->  Ki = -29120 / 560.6           = -51.9
```

## 10. Check

d = 2.006 > 0.70, p = 40.90 > 27.998, i = 51.9 > 0, p < 40 d = 80.2, and
i < (80.2 - 40.9)(560.6 x 40.9 - 15696)/1600 = 177.6. Also b1 = 551.2 > 0 and c1 = 7232 - 40 x 29120/551.2 = 5119 > 0.

`stability_analysis.py` compares these conditions with the closed-loop poles computed by the control library on thousands of random gain sets (5000 by default; all agree), and
`design_and_simulate.py` simulates the design.

## 11. Standard PID vs I-PD

With the standard PID on the error, the reference-to-output transfer function contains the controller zeros, and
the step response overshoots by 268 % (0.80 s settling). Putting Kp and Kd on the measurement and only the integral
on the error (I-PD) gives the same characteristic polynomial, so the same poles and the same stability region,
but no zeros: overshoot 0 %, settling 0.47 s.
