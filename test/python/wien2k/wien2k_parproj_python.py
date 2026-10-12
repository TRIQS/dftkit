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
    np.testing.assert_allclose(a, b, atol=1e-7, rtol=0, equal_nan=False)


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
result = parproj._symmetrize_densmat(orbs, {}, [dict(perm=[1])], 1,
                                    representations, raw, 2, True, True)
for shell in raw:
    for actual, expected in zip(result[shell], raw[shell]):
        np.testing.assert_allclose(actual, expected, atol=1e-14, rtol=0)

_check('CaOs2_partial')

with tempfile.TemporaryDirectory() as tmp:
    path = os.path.join(tmp, 'P1.dmftsym')
    with open(path, 'w') as fh:
        fh.write('1\n1 2\nSym. op. 1\n0 0 0 1\n'
                 '1 0 0\n0 1 0\n0 0 1\nGlobal->local\n'
                 '1\n1 0 0\n0 1 0\n0 0 1\n0 0 0 1\n'
                 '2\n1 0 0\n0 1 0\n0 0 1\n90 0 0 1\n')
    _, _, refs = parproj._read_dmftsym(path)
    assert len(refs) == 2
    np.testing.assert_allclose(refs[1]['a'], np.pi / 2, atol=1e-14, rtol=0)

print('wien2k_parproj_python: ok')
