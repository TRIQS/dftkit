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
permutations = [[2, 1], [2, 1]]

for antiunitary in (False, True):
    matrices = [[np.eye(orbit['dim'], dtype=complex) for orbit in orbits]
                for _ in permutations]
    rotate_to_local_frames(matrices, permutations,
                           [antiunitary, antiunitary], orbits, rotations)
    for iorb, jorb in enumerate((2, 3, 0, 1)):
        source = rotations[iorb].conj() if antiunitary else rotations[iorb]
        expected = rotations[jorb].conj().T @ source
        np.testing.assert_allclose(matrices[0][iorb], expected,
                                   atol=1e-14, rtol=0)

print('wien2k_symmetry_frames: ok')
