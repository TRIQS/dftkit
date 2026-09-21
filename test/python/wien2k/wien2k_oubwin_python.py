################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.oubwin generator (triqs_dftkit.wien2k.oubwin)
# reproduces the dmftproj Fortran output: regenerate case.oubwinup / case.oubwindn
# from case.almblm{up,dn} + case.indmftpr and compare byte for byte to the
# committed dmftproj references (the same files the SOC converter test consumes).
#
# CaOs2 is spin-orbit + spin-polarized. The almblm fixtures are gzipped to keep
# the tree small (they are shared with the ctqmcout test, which needs the full
# projector payload).
#
# CaOs2_mode1 repeats the same physical band selection through the band-index
# projection path (proj_mode 1): the window line "-0.15 0.30 1" is scanned over
# all spins/k for the global min/max band index (bands 60..76), the same range
# the energy window picks, so the reference is identical and the mode-1 code in
# set_projections.f:70-88 is exercised against the dmftproj Fortran output.

import gzip
import os
import shutil
import tempfile

from triqs_dftkit.wien2k import oubwin

HERE = os.path.dirname(os.path.abspath(__file__))


def _gunzip(src, dst):
    with gzip.open(src, 'rb') as fi, open(dst, 'wb') as fo:
        shutil.copyfileobj(fi, fo)


def _check(case, almblm='CaOs2'):
    """Regenerate <case>.oubwin{up,dn} from the <almblm> almblm fixtures and
    compare byte for byte. `almblm` selects which gzipped almblm payload to feed
    (the band-index fixtures reuse the CaOs2 coefficients with a different
    indmftpr window line)."""
    tmp = tempfile.mkdtemp()
    try:
        shutil.copy(os.path.join(HERE, f'{case}.indmftpr'),
                    os.path.join(tmp, f'{case}.indmftpr'))
        for ext in ('almblmup', 'almblmdn'):
            _gunzip(os.path.join(HERE, f'{almblm}.{ext}.gz'),
                    os.path.join(tmp, f'{case}.{ext}'))
        oubwin.write_oubwin(os.path.join(tmp, case))
        for ext in ('oubwinup', 'oubwindn'):
            got = open(os.path.join(tmp, f'{case}.{ext}')).read()
            ref = open(os.path.join(HERE, f'{case}.{ext}')).read()
            assert got == ref, f'{case}.{ext} differs from dmftproj reference'
    finally:
        shutil.rmtree(tmp)


_check('CaOs2')
_check('CaOs2_mode1')                 # proj_mode 1: band-index window 60..76

print('wien2k_oubwin_python: ok')
