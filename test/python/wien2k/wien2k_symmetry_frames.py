################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

import numpy as np

from triqs_dftkit.wien2k._symmetry import rotate_to_local_frames


# Two angular momenta on each atom make an atom/sort-only lookup ambiguous.
orbits = [
    dict(atom=1, sort=1, l=0, dim=1),
    dict(atom=1, sort=1, l=1, dim=2),
    dict(atom=2, sort=1, l=0, dim=1),
    dict(atom=2, sort=1, l=1, dim=2),
]
r0 = np.array([[1j]])
r1 = np.array([[0, 1], [1j, 0]], dtype=complex)
r2 = np.array([[1.]])
r3 = np.array([[1, 1j], [1j, 1]], dtype=complex) / np.sqrt(2)
rotations = [r0, r1, r2, r3]
permutations = [[2, 1], [1, 2]]
rng = np.random.default_rng(1729)


def random_matrix(dim):
    return rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))


def action(matrix, density, antiunitary):
    source = density.conj() if antiunitary else density
    return matrix @ source @ matrix.conj().T

for antiunitary in (False, True):
    matrices = [[np.linalg.qr(random_matrix(orbit['dim']))[0] for orbit in orbits]
                for _ in permutations]
    global_matrices = [[matrix.copy() for matrix in row] for row in matrices]
    rotate_to_local_frames(matrices, permutations,
                           [antiunitary, antiunitary], orbits, rotations)
    # Independent covariance oracle: acting in the global frame and then
    # changing frame must equal changing frame first and acting locally.
    for isym, targets in enumerate(((2, 3, 0, 1), (0, 1, 2, 3))):
        for iorb, jorb in enumerate(targets):
            seed = random_matrix(orbits[iorb]['dim'])
            density = seed + seed.conj().T
            local_density = rotations[iorb].conj().T @ density @ rotations[iorb]
            global_result = action(global_matrices[isym][iorb], density, antiunitary)
            expected = rotations[jorb].conj().T @ global_result @ rotations[jorb]
            actual = action(matrices[isym][iorb], local_density, antiunitary)
            np.testing.assert_allclose(actual, expected, atol=1e-13, rtol=0)

print('wien2k_symmetry_frames: ok')
