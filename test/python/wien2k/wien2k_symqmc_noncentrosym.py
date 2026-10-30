################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Check unitarity of the cubic p transform and spinor symmetry matrices for
# non-centrosymmetric tellurium (P3_1 21, point group 32). The comparison with
# double-precision native dmftproj output is covered by wien2k_te_python.py.

import os

import numpy as np

from triqs_dftkit.wien2k import symqmc
from triqs_dftkit.wien2k._dmftproj import reptrans

HERE = os.path.dirname(os.path.abspath(__file__))

P1 = reptrans('cubic', 1)
assert np.max(np.abs(P1 @ np.conj(P1.T) - np.eye(3))) < 1e-12, 'cubic l=1 not unitary'

_, ops = symqmc._read_dmftsym(os.path.join(HERE, 'Te.dmftsym'))
shells, so = symqmc._read_correlated_shells(os.path.join(HERE, 'Te.indmftpr'),
                                        os.path.join(HERE, 'Te.struct'))


def _timeinv(op, so):
    k = op['krotm']
    return 1 if (so and k[0, 0] * k[1, 1] - k[0, 1] * k[1, 0] < 0.0) else 0


for op in ops:
    for sh in shells:
        mat = symqmc._shell_matrix(op, sh, _timeinv(op, so), bool(so))
        dev = np.max(np.abs(mat @ np.conj(mat.T) - np.eye(mat.shape[0])))
        assert dev < 1e-7, ('not unitary', dev)

print('wien2k_symqmc_noncentrosym: ok')
