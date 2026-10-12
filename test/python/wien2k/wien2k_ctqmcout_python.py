################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.ctqmcout generator (triqs_dftkit.wien2k.ctqmcout)
# reproduces the dmftproj Fortran correlated-shell projectors: regenerate
# case.ctqmcout from case.almblm{up,dn} + case.indmftpr + case.struct +
# case.dmftsym and compare to the committed dmftproj reference.
#
# CaOs2 is spin-orbit + spin-polarized, cubic, two equivalent correlated d atoms.
# Integer fields must match exactly; floats are compared at 1e-7: dmftproj stores
# the cubic basis transform through a single-precision CMPLX cast, so its
# projectors carry ~1e-7 noise that the double-precision generator reproduces to
# the same floor (cf. the symqmc mixing test).

import gzip
import os
import shutil
import tempfile
import numpy as np

from triqs_dftkit.wien2k import ctqmcout

HERE = os.path.dirname(os.path.abspath(__file__))


def _gunzip(src, dst):
    with gzip.open(src, 'rb') as fi, open(dst, 'wb') as fo:
        shutil.copyfileobj(fi, fo)


def _compare(got_path, ref_path):
    got = open(got_path).read().split()
    ref = open(ref_path).read().split()
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
        for ext in ('indmftpr', 'struct', 'dmftsym'):
            shutil.copy(os.path.join(HERE, f'{case}.{ext}'),
                        os.path.join(tmp, f'{case}.{ext}'))
        for ext in ('almblmup', 'almblmdn'):
            _gunzip(os.path.join(HERE, f'{case}.{ext}.gz'),
                    os.path.join(tmp, f'{case}.{ext}'))
        ctqmcout.write_ctqmcout(os.path.join(tmp, case))
        _compare(os.path.join(tmp, f'{case}.ctqmcout'),
                 os.path.join(HERE, f'{case}.ctqmcout'))
        # A mismatched SO spin window must fail before replacing the output.
        dn = os.path.join(tmp, f'{case}.almblmdn')
        with open(dn) as f:
            lines = f.readlines()
        lines[3] = '100.0\n'  # Shift all down-spin bands outside the window.
        with open(dn, 'w') as f:
            f.writelines(lines)
        try:
            ctqmcout.write_ctqmcout(os.path.join(tmp, case))
        except ValueError as error:
            assert 'identical up/dn band windows' in str(error), error
        else:
            raise AssertionError('mismatched SO spin windows were accepted')
    finally:
        shutil.rmtree(tmp)


assert ctqmcout._select_window(10, 13, [-3., -2., -1., 0.5], 0., 1.) == (False, 0, 0)
assert ctqmcout._select_window(5, 5, [0.5], 0., 1.) == (False, 0, 0)
_check('CaOs2')

for scale in (1., 1e-20):
    overlap = np.diag([scale, 2 * scale]).astype(complex)
    inverse = ctqmcout._sqrt_inv(overlap)
    np.testing.assert_allclose(inverse @ overlap @ inverse.conj().T,
                               np.eye(2), atol=1e-14, rtol=0)
try:
    ctqmcout._sqrt_inv(np.diag([0., 1.]).astype(complex))
except ValueError as error:
    assert 'rank deficient' in str(error)
else:
    raise AssertionError('singular overlap was accepted')

with tempfile.TemporaryDirectory() as tmp:
    path = os.path.join(tmp, 'multi.indmftpr')
    with open(path, 'w') as fh:
        fh.write('1\n1\n2\ncomplex\n0 2 2\n0 1 1\n1\n1\n1\n-2 3\n')
    parsed = ctqmcout._read_indmftpr(path)
    assert parsed['so'] == 1 and parsed['sorts'][0]['correlated_ls'] == [1, 2]

print('wien2k_ctqmcout_python: ok')
