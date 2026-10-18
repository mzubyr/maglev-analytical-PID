import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from plant import (A, B, N_GAIN, DEN, char_poly_pd, char_poly_pid, hurwitz_pd, hurwitz_pid,
                   is_stable_roots)


def test_model_numbers():
    assert abs(A[1, 0] - 392.4) < 1e-6
    assert abs(A[1, 2] + 28.03) < 0.01
    assert abs(A[2, 2] + 40) < 1e-9 and abs(B[2, 0] - 20) < 1e-9


def test_transfer_function():
    assert abs(N_GAIN + 560.6) < 0.1
    assert np.allclose(DEN, [1, 40, -392.4, -15696], atol=1e-6)


def test_open_loop_unstable_and_controllable():
    assert np.linalg.eigvals(A).real.max() > 19
    ctrb = np.hstack([np.linalg.matrix_power(A, k) @ B for k in range(3)])
    assert np.linalg.matrix_rank(ctrb) == 3 and np.linalg.det(ctrb) < 0


def test_routh_hurwitz_matches_roots():
    rng = np.random.default_rng(1)
    for _ in range(3000):
        kp, kd, ki = -rng.uniform(0, 150), -rng.uniform(0, 6), -rng.uniform(0, 300)
        assert hurwitz_pd(kp, kd) == is_stable_roots(char_poly_pd(kp, kd))
        assert hurwitz_pid(kp, kd, ki) == is_stable_roots(char_poly_pid(kp, kd, ki))


def test_designed_gains_place_the_poles():
    kd, kp, ki = -2.006, -40.90, -51.95
    roots = np.sort_complex(np.roots(char_poly_pid(kp, kd, ki)))
    assert np.all(roots.real < -3.9)
    assert hurwitz_pid(kp, kd, ki)


def test_wrong_sign_gains_are_unstable():
    assert not is_stable_roots(char_poly_pid(20, 1, 5))
