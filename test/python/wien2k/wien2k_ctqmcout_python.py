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
# Integer fields must match exactly; floats are compared at 1e-6: dmftproj stores
# the cubic basis transform through a single-precision CMPLX cast, so its
# projectors carry ~1e-7 noise that the double-precision generator reproduces to
# the same floor (cf. the symqmc mixing test).

import gzip
import os
import shutil
import tempfile

from triqs_dftkit.wien2k import ctqmcout

HERE = os.path.dirname(os.path.abspath(__file__))


def _gunzip(src, dst):
    with gzip.open(src, 'rb') as fi, open(dst, 'wb') as fo:
        shutil.copyfileobj(fi, fo)


def _compare(got_path, ref_path):
    got = open(got_path).read().split()
    ref = open(ref_path).read().split()
    assert len(got) == len(ref), (len(got), len(ref))
    max_int, max_float = 0, 0.0
    for a, b in zip(got, ref):
        ia, ib = a.lstrip('-').isdigit(), b.lstrip('-').isdigit()
        if ia and ib:
            max_int = max(max_int, abs(int(a) - int(b)))
        else:
            max_float = max(max_float, abs(float(a) - float(b)))
    assert max_int == 0, f'integer field differs (max {max_int})'
    assert max_float < 1e-6, f'float field differs (max {max_float:.2e})'


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
    finally:
        shutil.rmtree(tmp)


_check('CaOs2')

print('wien2k_ctqmcout_python: ok')
