################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.symqmc generator on a non-centrosymmetric
# spin-orbit cell (tellurium, P3_1 21, point group 32, a correlated p shell),
# the regime dft_tools #148 is about and the centrosymmetric d-shell fixtures
# (CaOs2, SrVO3) never exercised. Two convention bugs hid there:
#
#   - the cubic l = 1 (p) transform must be unitary and px,py,pz ordered;
#   - the magnetic (time-inversion) operations need the antiunitary basis
#     change P D P^T, not P D P^dag, which differs only for a complex transform
#     (the p_y row), so the real cubic-d path never showed it.
#
# Reference: the dmftproj Fortran case.symqmc for the same symmetry input. The
# Fortran reads the cubic template at single precision, so its floor is ~1e-7;
# the generator is double precision. We compare the matrix data at 1e-6 and
# separately assert every spinor matrix is unitary.

import os
import shutil
import tempfile

import numpy as np

from triqs_dftkit.wien2k import symqmc
from triqs_dftkit.wien2k._dmftproj import reptrans

HERE = os.path.dirname(os.path.abspath(__file__))

# The cubic l = 1 transform itself must be unitary (P P^dag = I). The d (l = 2)
# table was correct; only this one carried the un-normalized, mis-ordered rows.
P1 = reptrans('cubic', 1)
assert np.max(np.abs(P1 @ np.conj(P1.T) - np.eye(3))) < 1e-12, 'cubic l=1 not unitary'


def _matrix_data(path, nsym, natom):
    """Flat array of the matrix blocks, skipping header/perm/timeflag lines."""
    toks = iter(open(path).read().split())
    assert int(next(toks)) == nsym and int(next(toks)) == natom
    for _ in range(nsym * natom + nsym):     # perm rows + the time-reversal flags
        next(toks)
    return np.array([float(x) for x in toks])


tmp = tempfile.mkdtemp()
try:
    for src, dst in (('Te.dmftsym', 'case.dmftsym'),
                     ('Te.struct', 'case.struct'),
                     ('Te.indmftpr', 'case.indmftpr')):
        shutil.copy(os.path.join(HERE, src), os.path.join(tmp, dst))
    case = os.path.join(tmp, 'case')
    symqmc.write_symqmc(case)
    nsym, ops = symqmc._read_dmftsym(case + '.dmftsym')
    shells, so = symqmc._read_correlated_shells(case + '.indmftpr', case + '.struct')
    natom = len(ops[0]['perm'])
    data = _matrix_data(case + '.symqmc', nsym, natom)
finally:
    shutil.rmtree(tmp)

ref = np.load(os.path.join(HERE, 'wien2k_symqmc_te.ref.npy'))
assert data.shape == ref.shape, (data.shape, ref.shape)
assert np.max(np.abs(data - ref)) < 1e-6, np.max(np.abs(data - ref))


def _timeinv(op, so):
    k = op['krotm']
    return 1 if (so and k[0, 0] * k[1, 1] - k[0, 1] * k[1, 0] < 0.0) else 0


for op in ops:
    for sh in shells:
        mat = symqmc._shell_matrix(op, sh, _timeinv(op, so), bool(so))
        dev = np.max(np.abs(mat @ np.conj(mat.T) - np.eye(mat.shape[0])))
        assert dev < 1e-7, ('not unitary', dev)

print('wien2k_symqmc_noncentrosym: ok')
