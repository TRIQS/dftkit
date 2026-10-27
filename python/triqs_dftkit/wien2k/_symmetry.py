################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

"""Small, NumPy-only helpers for Wien2k symmetry conversion."""

import numpy


def rotate_to_local_frames(mat, perm, time_inv, orbits, rot_mat):
    """Rotate global-frame symmetry matrices into each shell's solver frame.

    A symmetry maps source shell ``i`` to ``j``.  In local frames its matrix is
    ``R_j^dag M_ji R_i``.  For an antiunitary operation Symmetry.symmetrize
    conjugates the source quantity, so the source rotation is conjugated too.
    """
    lookup = {}
    for iorb, orbit in enumerate(orbits):
        key = (orbit['atom'], orbit['sort'], orbit['l'], orbit['dim'])
        if key in lookup:
            raise ValueError('duplicate shell identity while mapping symmetry frames')
        lookup[key] = iorb

    for isym, atom_perm in enumerate(perm):
        targets = []
        for orbit in orbits:
            key = (atom_perm[orbit['atom'] - 1], orbit['sort'],
                   orbit['l'], orbit['dim'])
            if key not in lookup:
                raise ValueError('symmetry operation maps to an unknown shell')
            targets.append(lookup[key])
        for iorb, jorb in enumerate(targets):
            source = rot_mat[iorb].conjugate() if time_inv[isym] else rot_mat[iorb]
            mat[isym][iorb] = (rot_mat[jorb].conjugate().transpose()
                               @ mat[isym][iorb] @ source)
    return mat
