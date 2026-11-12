import numpy as np
from vpstm.q import score


def test_q_zero_for_identical_source_and_ansatz():
    rng = np.random.default_rng(3)
    f = rng.normal(size=(20, 4)) + 1j * rng.normal(size=(20, 4))
    fz = rng.normal(size=(20, 4)) + 1j * rng.normal(size=(20, 4))
    fzz = rng.normal(size=(20, 4)) + 1j * rng.normal(size=(20, 4))
    q, _ = score(f, fz, fzz, f, fz, fzz,
                 np.ones(4), np.ones(20), np.array([-1., -.5, 0., .5]),
                 vacuum_reference=5.0, distance=6.0)
    assert q < 1e-12
