################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.oubwin generator (triqs_dftkit.wien2k.oubwin)
# reproduces the dmftproj Fortran output: regenerate case.oubwinup / case.oubwindn
# from case.almblm{up,dn} + case.indmftpr and compare to the native reference,
# ignoring only trailing whitespace removed from the shipped fixture.
#
# CaOs2 is spin-orbit + spin-polarized, so both outputs match CaOs2.oubwinup.
# Reuse the compressed almblm inputs shipped with the SOC fixture.

import gzip
import os
import shutil
import tempfile
import numpy as np

from triqs_dftkit.wien2k import oubwin
from triqs_dftkit.wien2k import _dmftproj as dp

HERE = os.path.dirname(os.path.abspath(__file__))


def _check(case):
    tmp = tempfile.mkdtemp()
    try:
        for ext in ('almblmup', 'almblmdn'):
            with gzip.open(os.path.join(HERE, f'{case}.{ext}.gz'), 'rb') as src, \
                    open(os.path.join(tmp, f'{case}.{ext}'), 'wb') as dst:
                shutil.copyfileobj(src, dst)
        shutil.copy(os.path.join(HERE, f'{case}.indmftpr'), tmp)
        oubwin.write_oubwin(os.path.join(tmp, case))
        with open(os.path.join(HERE, f'{case}.oubwinup')) as f:
            ref = [line.rstrip() for line in f]
        for ext in ('oubwinup', 'oubwindn'):
            with open(os.path.join(tmp, f'{case}.{ext}')) as f:
                got = [line.rstrip() for line in f]
            assert got == ref, f'{case}.{ext} differs from dmftproj reference'
    finally:
        shutil.rmtree(tmp)


# Independent expected band ranges for set_projections.f's branch ordering.
for nbmin, energies, expected in (
        (10, [-3., -2., -1., 0.5], (False, 0, 0)),  # opens only at nbmax
        (5, [0.5], (False, 0, 0)),                  # nbmin == nbmax
        (10, [0., 0.5, 1., 2.], (True, 11, 12)),   # strict lower, inclusive upper
        (10, [0.5, 1.], (True, 10, 11)),           # open window reaches nbmax
        (10, [-1., 0., 2.], (False, 0, 0))):       # no in-window band
    got = dp.select_window(nbmin, nbmin + len(energies) - 1,
                                np.array(energies), 0., 1.)
    assert got == expected, (energies, got, expected)

# Three k-points with distinct weights, a Fermi shift of 10 Ry, and one
# excluded point. The lowest-band weight is used even if that band is outside
# the window; an excluded point emits no band range or weight record.
with tempfile.TemporaryDirectory() as tmp:
    almblm = os.path.join(tmp, 'case.almblm')
    with open(almblm, 'w') as f:
        f.write('0\n3\n0\n10.0\n1.0\n0\n')
        for ik, nbmin, energies, weights in (
                (1, 10, [9., 10.25, 10.75, 12.], [.25, .80, .60, .25]),
                (2, 20, [8., 9.], [.10, .10]),
                (3, 30, [10.25, 10.75], [.75, .50])):
            f.write(f'IK = {ik}\n0 0 0\n{ik} {nbmin} {nbmin + len(energies) - 1}\n')
            for energy, weight in zip(energies, weights):
                f.write(f'{weight} {energy}\n')
            f.write('ATOM = 1\n1\n')
            f.write('(0.,0.) (0.,0.)\n' * len(energies))
    windows = oubwin._windows_from_spin(dp.read_almblm(almblm, dict(nsort=1, lmax=0, mult=[1])), dict(proj_mode=0, e_bot=0., e_top=1.), None)
    assert windows == [(True, 11, 12, 0.25), (False, 0, 0, 0.10),
                       (True, 30, 31, 0.75)], windows
    output = os.path.join(tmp, 'case.oubwin')
    oubwin._write_oubwin_file(output, 0, windows)
    with open(output) as f:
        assert [line.split() for line in f] == [
            ['3'], ['0'], ['1'], ['11', '12'], ['0.2500000000000000'],
            ['0'], ['1'], ['30', '31'], ['0.7500000000000000']]

_check('CaOs2')

print('wien2k_oubwin_python: ok')
