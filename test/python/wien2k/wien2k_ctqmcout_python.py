################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.ctqmcout generator (triqs_dftkit.wien2k.ctqmcout)
# reproduces the dmftproj Fortran correlated-shell projectors to machine
# precision: regenerate case.ctqmcout from the almblm/struct/dmftsym + a wide
# energy window and compare to the precision-fixed dmftproj reference.
#
# The window (-2.0, 3.0) keeps 50 bands, more than the 20 correlated spin-orbitals
# of the two SOC Os d shells, so the Loewdin overlap O is full rank. (The physical
# narrow-window CaOs2 of the SOC converter test makes O rank-deficient: O^{-1/2}
# then amplifies the last-ULP libm difference between gfortran and numpy, which is
# a numerical property of that degenerate case, not of the port.) Reference built
# with the precision-fixed dmftproj (PR #14) and full-precision templates.

import gzip
import os
import shutil
import tempfile

from triqs_dftkit.wien2k import ctqmcout

HERE = os.path.dirname(os.path.abspath(__file__))


def _gunzip(src, dst):
    with gzip.open(src, 'rb') as fi, open(dst, 'wb') as fo:
        shutil.copyfileobj(fi, fo)


def _compare(got_path, ref_lines):
    got = open(got_path).read().split()
    ref = ''.join(ref_lines).split()
    assert len(got) == len(ref), (len(got), len(ref))
    max_int, max_float = 0, 0.0
    for a, b in zip(got, ref):
        if a.lstrip('-').isdigit() and b.lstrip('-').isdigit():
            max_int = max(max_int, abs(int(a) - int(b)))
        else:
            max_float = max(max_float, abs(float(a) - float(b)))
    assert max_int == 0, f'integer field differs (max {max_int})'
    assert max_float < 1e-11, f'float field differs (max {max_float:.2e})'


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
        ctqmcout.write_ctqmcout(os.path.join(tmp, case))
        with gzip.open(os.path.join(HERE, f'{case}.ctqmcout.gz'), 'rt') as fh:
            _compare(os.path.join(tmp, f'{case}.ctqmcout'), fh.readlines())
    finally:
        shutil.rmtree(tmp)


_check('CaOs2_full')

print('wien2k_ctqmcout_python: ok')
