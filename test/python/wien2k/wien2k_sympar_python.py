################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.sympar generator (triqs_dftkit.wien2k.sympar)
# reproduces the dmftproj Fortran output. sympar holds the symmetry matrices for
# all INCLUDED orbitals (correlated + partial), so the CaOs2_partial fixture adds
# an uncorrelated Ca d-shell to the correlated Os d-shells. Inputs are the shared
# CaOs2 almblm/struct/dmftsym renamed to the CaOs2_partial case.
#
# Floats compared at 1e-11 against native dmftproj output (full-precision
# cubic-d template).

import gzip
import os
import shutil
import tempfile
import numpy as np

from triqs_dftkit.wien2k import sympar

HERE = os.path.dirname(os.path.abspath(__file__))


def _gunzip(src, dst):
    with gzip.open(src, 'rb') as fi, open(dst, 'wb') as fo:
        shutil.copyfileobj(fi, fo)


def _compare(got_path, ref_lines):
    got = open(got_path).read().split()
    ref = ''.join(ref_lines).split()
    assert len(got) == len(ref), (len(got), len(ref))
    # 16 operations, three atoms, time flags, then three 10x10 shell matrices.
    header = 2 + 16 * 3 + 16
    np.testing.assert_array_equal(np.array(got[:header], int),
                                  np.array(ref[:header], int))
    data = np.array(got[header:], float).reshape(16, 3, 2, 10, 10)
    reference = np.array(ref[header:], float).reshape(data.shape)
    matrices = data[:, :, 0] + 1j * data[:, :, 1]
    expected = reference[:, :, 0] + 1j * reference[:, :, 1]
    np.testing.assert_allclose(matrices, expected, atol=1e-11, rtol=0,
                               equal_nan=False)
    for operation in matrices:
        for matrix in operation:
            np.testing.assert_allclose(matrix @ matrix.conj().T, np.eye(10),
                                       atol=1e-11, rtol=0)


def _check(case):
    tmp = tempfile.mkdtemp()
    try:
        shutil.copy(os.path.join(HERE, f'{case}.indmftpr'),
                    os.path.join(tmp, f'{case}.indmftpr'))
        for ext in ('struct', 'dmftsym'):
            shutil.copy(os.path.join(HERE, f'CaOs2.{ext}'),
                        os.path.join(tmp, f'{case}.{ext}'))
        for ext in ('almblmup', 'almblmdn'):
            _gunzip(os.path.join(HERE, f'CaOs2.{ext}.gz'),
                    os.path.join(tmp, f'{case}.{ext}'))
        sympar.write_sympar(os.path.join(tmp, case))
        with gzip.open(os.path.join(HERE, f'{case}.sympar.gz'), 'rt') as fh:
            _compare(os.path.join(tmp, f'{case}.sympar'), fh.readlines())
    finally:
        shutil.rmtree(tmp)


# Inversion has proper part I but parity -1: odd-l matrices must be -I.
inversion = dict(a=0., b=0., c=0., iprop=-1, krotm=np.eye(3))
for l in (1, 3):
    shell = dict(l=l, mixing=False, P=np.eye(2 * l + 1, dtype=complex))
    np.testing.assert_allclose(sympar._shell_matrix(inversion, shell, False),
                               -np.eye(2 * (2 * l + 1)), atol=1e-14, rtol=0)
p = sympar.reptrans('cubic', 1)
np.testing.assert_allclose(p @ p.conj().T, np.eye(3), atol=1e-11, rtol=0)

_check('CaOs2_partial')

print('wien2k_sympar_python: ok')
