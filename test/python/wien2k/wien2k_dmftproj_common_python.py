################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Unit tests for the shared dmftproj machinery (triqs_dftkit.wien2k._dmftproj),
# the module the five case.* generators are built on. The generators' golden
# tests cover the assembled output; these exercise the isolated primitives
# directly so a regression in a shared unit is caught at its source.

import os
import tempfile

import numpy as np

from triqs_dftkit.wien2k import _dmftproj as dp

HERE = os.path.dirname(os.path.abspath(__file__))


# --- list-directed numeric parsing -------------------------------------------

assert dp.to_float('1.0D-3') == 1.0e-3
assert dp.to_float('4.57E-002') == 4.57e-2
assert dp.to_float('-2.5d0') == -2.5

assert dp.to_complex('(1.0,-2.0)') == complex(1.0, -2.0)


def _reader(text):
    fh = tempfile.NamedTemporaryFile('w', suffix='.tmp', delete=False)
    fh.write(text)
    fh.close()
    return dp.Reader(fh.name)


r = _reader('(1.0,2.0)\n(3.0,-4.0) (5.0,6.0)\n')
assert dp.read_complex(r) == complex(1.0, 2.0)
a, b = dp.read_two_complex(r)
assert a == complex(3.0, -4.0) and b == complex(5.0, 6.0)


# --- angular basis: cubic transform, with and without the float32 cast -------

exact = dp.reptrans('cubic', 2, cast=False)
cast = dp.reptrans('cubic', 2, cast=True)
assert exact.shape == (5, 5)
# unitary
assert np.max(np.abs(exact @ np.conj(exact.T) - np.eye(5))) < 1e-12
# the cast truncates 1/sqrt2 to single precision (the dmftproj #148 noise)
assert abs(abs(exact[1, 0]) - 2 ** -0.5) < 1e-15
assert 0 < abs(abs(cast[1, 0]) - 2 ** -0.5) < 1e-6
assert dp.reptrans('complex', 2).shape == (5, 5)
assert np.allclose(dp.reptrans('complex', 2), np.eye(5))


# --- Wigner D and orbital time reversal --------------------------------------

for l in (1, 2):
    D0 = dp.dmat(l, 0.0, 0.0, 0.0, 1.0)
    assert np.max(np.abs(D0 - np.eye(2 * l + 1))) < 1e-12        # zero rotation
    D = dp.dmat(l, 0.3, 0.7, 1.1, 1.0)
    assert np.max(np.abs(D @ np.conj(D.T) - np.eye(2 * l + 1))) < 1e-12  # unitary
    T = dp.tmat(l)
    for m in range(-l, l + 1):
        assert T[-m + l, m + l] == (-1) ** m


# --- select_window: contiguous band range in (e1, e2] ------------------------

eband = np.array([-1.0, -0.1, 0.05, 0.2, 0.5])     # bands 10..14, window (-0.15, 0.30]
incl, lo, hi = dp.select_window(10, 14, eband, -0.15, 0.30)
assert incl and lo == 11 and hi == 13


# --- band-index window (proj_mode 1/2) ---------------------------------------

# set_projections.f:70-88: every k included, indices clamped to nbmin/nbmax.
incl, lo, hi = dp.select_band_window(1, 92, 60, 76)
assert incl and lo == 60 and hi == 76
incl, lo, hi = dp.select_band_window(50, 80, 43, 92)     # clamp both ends
assert incl and lo == 50 and hi == 80

# proj_mode 2 takes b_bot/b_top straight from the (rounded) window line.
m2 = {'proj_mode': 2, 'e_bot': 43.0, 'e_top': 92.0}
assert dp.band_index_window(m2, []) == (43, 92)

# proj_mode 1 scans the (e_bot, e_top] energies for the global band-index range.
m1 = {'proj_mode': 1, 'e_bot': -0.15, 'e_top': 0.30}
spins = [{'kp': [{'nbmin': 1, 'nbmax': 5,
                  'eband': np.array([-1.0, -0.1, 0.05, 0.2, 0.5])}]}]
assert dp.band_index_window(m1, spins) == (2, 4)


# --- case.indmftpr / case.dmftsym structured parse ---------------------------

info = dp.read_indmftpr(os.path.join(HERE, 'CaOs2.indmftpr'))
assert info['nsort'] == 2 and info['mult'] == [1, 2] and info['lmax'] == 3
assert info['so'] == 1
assert info['sorts'][1]['correlated_ls'] == [2]      # Os d correlated
assert info['sorts'][0]['correlated_ls'] == []       # Ca nothing

part = dp.read_indmftpr(os.path.join(HERE, 'CaOs2_partial.indmftpr'))
assert part['sorts'][0]['included_ls'] == [2]        # Ca d now an included shell
assert part['sorts'][0]['correlated_ls'] == []
assert part['sorts'][1]['correlated_ls'] == [2]

nsym, ops = dp.read_dmftsym(os.path.join(HERE, 'CaOs2.dmftsym'))
assert nsym == 16 and len(ops) == 16
assert all(o['krotm'].shape == (3, 3) and len(o['perm']) == 3 for o in ops)

print('wien2k_dmftproj_common: ok')
