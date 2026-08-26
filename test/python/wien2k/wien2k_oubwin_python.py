################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.oubwin generator (triqs_dftkit.wien2k.oubwin)
# reproduces the dmftproj Fortran output: regenerate case.oubwinup / case.oubwindn
# from case.almblm{up,dn} + case.indmftpr and compare byte for byte to the
# committed dmftproj references (the same files the SOC converter test consumes).
#
# CaOs2 is spin-orbit + spin-polarized; the almblm fixtures are truncated to the
# header and first-sort band block the window selection reads (the projector
# payload that follows is not needed for oubwin).

import os
import shutil
import tempfile

from triqs_dftkit.wien2k import oubwin

HERE = os.path.dirname(os.path.abspath(__file__))


def _check(case):
    tmp = tempfile.mkdtemp()
    try:
        for ext in ('almblmup', 'almblmdn', 'indmftpr'):
            shutil.copy(os.path.join(HERE, f'{case}.{ext}'),
                        os.path.join(tmp, f'{case}.{ext}'))
        oubwin.write_oubwin(os.path.join(tmp, case))
        for ext in ('oubwinup', 'oubwindn'):
            got = open(os.path.join(tmp, f'{case}.{ext}')).read()
            ref = open(os.path.join(HERE, f'{case}.{ext}')).read()
            assert got == ref, f'{case}.{ext} differs from dmftproj reference'
    finally:
        shutil.rmtree(tmp)


_check('CaOs2')

print('wien2k_oubwin_python: ok')
