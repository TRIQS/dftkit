################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.parproj generator (triqs_dftkit.wien2k.parproj)
# reproduces the dmftproj Fortran output. parproj holds the partial projectors
# and density matrices for the INCLUDED orbitals; the CaOs2_partial fixture adds
# an uncorrelated Ca d-shell. Inputs are the shared CaOs2 almblm/struct/dmftsym
# renamed to the CaOs2_partial case.
#
# Floats compared at 1e-7 (dmftproj single-precision basis-transform floor).

import gzip
import os
import shutil
import tempfile
import numpy as np

from triqs_dftkit.wien2k import parproj

HERE = os.path.dirname(os.path.abspath(__file__))


def _gunzip(src, dst):
    with gzip.open(src, 'rb') as fi, open(dst, 'wb') as fo:
        shutil.copyfileobj(fi, fo)


def _compare(got_path, ref_lines):
    got = open(got_path).read().split()
    ref = ''.join(ref_lines).split()
    assert len(got) == len(ref), (len(got), len(ref))
    floats = []
    for a, b in zip(got, ref):
        if b.lstrip('-').isdigit():
            assert int(a) == int(b), (a, b)
        else:
            floats.append((float(a), float(b)))
    a, b = np.array(floats).T
    np.testing.assert_allclose(a, b, atol=1e-11, rtol=0, equal_nan=False)


def _check(case, mode=0):
    tmp = tempfile.mkdtemp()
    try:
        shutil.copy(os.path.join(HERE, f'{case}.indmftpr'),
                    os.path.join(tmp, f'{case}.indmftpr'))
        if mode:
            path = os.path.join(tmp, f'{case}.indmftpr')
            with open(path) as fh:
                lines = fh.read().splitlines()
            lines[-1] = lines[-1] + ' 1' if mode == 1 else '60 76 2'
            with open(path, 'w') as fh:
                fh.write('\n'.join(lines) + '\n')
        shutil.copy(os.path.join(HERE, 'jbasis_d.dat'), os.path.join(tmp, 'jbasis_d.dat'))
        for ext in ('struct', 'dmftsym'):
            shutil.copy(os.path.join(HERE, f'CaOs2.{ext}'),
                        os.path.join(tmp, f'{case}.{ext}'))
        for ext in ('almblmup', 'almblmdn'):
            _gunzip(os.path.join(HERE, f'CaOs2.{ext}.gz'),
                    os.path.join(tmp, f'{case}.{ext}'))
        parproj.write_parproj(os.path.join(tmp, case))
        with gzip.open(os.path.join(HERE, f'{case}.parproj.gz'), 'rt') as fh:
            _compare(os.path.join(tmp, f'{case}.parproj'), fh.readlines())
    finally:
        shutil.rmtree(tmp)


# Identity symmetry must preserve distinct s and p density blocks on one atom.
orbs = [dict(atom=1, sort=1, l=l) for l in (0, 1)]
raw = {(1, l): [np.eye(2 * l + 1, dtype=complex) * value
                for value in (1., 2., 0.25, 0.25)] for l in (0, 1)}
representations = {(l, 1, 0): (np.eye(2 * l + 1), False, 0.) for l in (0, 1)}
result = parproj._symmetrize_densmat(orbs, [dict(perm=[1])], 1,
                                    representations, raw)
for shell in raw:
    for actual, expected in zip(result[shell], raw[shell]):
        np.testing.assert_allclose(actual, expected, atol=1e-14, rtol=0)

_check('CaOs2_partial')
_check('CaOs2_jbasis_full')
# All definitions select bands 60..76 here; the partial projectors and
# below-window density must match the native energy-window reference.
_check('CaOs2_partial', mode=1)
_check('CaOs2_partial', mode=2)

print('wien2k_parproj_python: ok')

# Equivalent atoms must share the average after identity and swap operations.
raw_mixing = {1: np.diag([1., 3.]), 2: np.diag([5., 7.])}
actual = parproj._symmetrize_densmat_mixing(
    [dict(atom=1), dict(atom=2)], [dict(perm=[1, 2]), dict(perm=[2, 1])], 2,
    [(np.eye(2), False)] * 2, raw_mixing)
for matrix in actual.values():
    np.testing.assert_allclose(matrix, np.diag([3., 5.]), atol=1e-14, rtol=0)
