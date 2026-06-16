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
# Floats compared at 1e-6: dmftproj stores the cubic basis transform through a
# single-precision CMPLX cast, so its matrices carry ~1e-7 noise the double-
# precision generator reproduces to the same floor.

import gzip
import os
import shutil
import tempfile

from triqs_dftkit.wien2k import sympar

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
        sympar.write_sympar(os.path.join(tmp, case))
        with gzip.open(os.path.join(HERE, f'{case}.sympar.gz'), 'rt') as fh:
            _compare(os.path.join(tmp, f'{case}.sympar'), fh.readlines())
    finally:
        shutil.rmtree(tmp)


_check('CaOs2_partial')

print('wien2k_sympar_python: ok')
