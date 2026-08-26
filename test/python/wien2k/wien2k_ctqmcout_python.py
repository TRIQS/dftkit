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
#
# CaOs2_jbasis_full exercises the spin-mixing fromfile (|j, m_j>) path: the Os d
# shell is a 2(2l+1)=10 spinor basis (jbasis_d.dat), so the projector, Rloc
# rotrep and complex-harmonics transform blocks are the full mixing matrices the
# dft_tools #148 path singles out. Same wide window for the full-rank Loewdin.
#
# CaOs2_mode1_full and CaOs2_mode2_full exercise the band-index projection modes
# (set_projections.f:70-88). Mode 1 ("-2.0 3.0 1") scans all spins/k for the
# global band-index window; mode 2 ("43 92 2") takes the indices straight from
# the window line. Both resolve to bands 43..92 (50 bands > 20 correlated
# spin-orbitals), so the Loewdin overlap stays full rank and the projectors
# match the precision-fixed dmftproj to machine precision.

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
        if os.path.exists(os.path.join(HERE, 'jbasis_d.dat')):
            shutil.copy(os.path.join(HERE, 'jbasis_d.dat'),
                        os.path.join(tmp, 'jbasis_d.dat'))
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
_check('CaOs2_jbasis_full')
_check('CaOs2_mode1_full')            # proj_mode 1: band-index window 43..92
_check('CaOs2_mode2_full')            # proj_mode 2: explicit band indices 43 92

print('wien2k_ctqmcout_python: ok')
