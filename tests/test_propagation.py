import numpy as np
from vpstm.propagation import propagate, height_slices


def test_complex_boundary_and_height_sampling():
    potential = np.full((8, 6, 6), 5.0)
    boundary = np.ones((6, 6), dtype=complex) * (1.0 + 0.2j)
    values, meta = propagate(boundary, potential, -1.0,
                              [0.5, 0.5, 0.5], tolerance=1e-7)
    assert values.dtype.kind == 'c'
    assert np.allclose(values[0], boundary)
    assert meta['relative_residual'] < 1e-5
    sampled = height_slices(values, np.arange(8) * 0.5, [0.5, 1.0])
    assert sampled.shape == (2, 6, 6)
