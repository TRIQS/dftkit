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


# --- angular basis: cubic transform -------------------------------------------

exact = dp.reptrans('cubic', 2)
assert exact.shape == (5, 5)
# unitary
assert np.max(np.abs(exact @ np.conj(exact.T) - np.eye(5))) < 1e-12
assert abs(abs(exact[1, 0]) - 2 ** -0.5) < 1e-15
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

# set_projections.f:70-88: indices clamped to nbmin/nbmax.
incl, lo, hi = dp.select_band_window(1, 92, 60, 76)
assert incl and lo == 60 and hi == 76
incl, lo, hi = dp.select_band_window(50, 80, 43, 92)     # clamp both ends
assert incl and lo == 50 and hi == 80
# the below-window (1, b_bot - 1) is empty once b_bot <= nbmin
assert not dp.select_band_window(5, 92, 1, 4)[0]
assert not dp.select_band_window(1, 92, 1, -1)[0]

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

# Atomic sorts and symmetry operations are independent counts (P1 can have
# many inequivalent atoms but only the identity operation).
with tempfile.TemporaryDirectory() as tmp:
    path = os.path.join(tmp, 'P1.dmftsym')
    with open(path, 'w') as fh:
        fh.write('1\n1 2\nSym. op. 1\n0 0 0 1\n'
                 '1 0 0\n0 1 0\n0 0 1\nGlobal->local\n'
                 '1\n1 0 0\n0 1 0\n0 0 1\n0 0 0 1\n'
                 '2\n1 0 0\n0 1 0\n0 0 1\n90 0 0 1\n')
    _, _, refs = dp.read_dmftsym(path, rotloc=True)
    assert len(refs) == 2
    np.testing.assert_allclose(refs[1]['a'], np.pi / 2, atol=1e-14, rtol=0)

    # A tiny but nonzero spin coupling must retain the full spinor basis.
    theta = 1e-9
    basis = np.eye(6)
    basis[0, 0] = basis[3, 3] = np.cos(theta)
    basis[0, 3], basis[3, 0] = np.sin(theta), -np.sin(theta)
    path = os.path.join(tmp, 'basis.dat')
    with open(path, 'w') as fh:
        for row in basis:
            fh.write(' ' + ' '.join(f'{value:.17E} 0' for value in row).replace('E', 'D') + '\n')
    transform, mixing = dp.read_fromfile(path, 1)
    assert mixing and transform.shape == (6, 6)
    np.testing.assert_array_equal(transform, basis.conjugate())

    # Each correlated l with an explicit irep owns its own correps record.
    path = os.path.join(tmp, 'multi.indmftpr')
    with open(path, 'w') as fh:
        fh.write('1\n1\n2\ncomplex\n0 2 2\n0 1 1\n1\n1\n1\n-2 3\n')
    parsed = dp.read_indmftpr(path)
    assert parsed['so'] == 1 and parsed['sorts'][0]['correlated_ls'] == [1, 2]

    # Selecting individual irreps is not implemented: do not silently export
    # the entire shell in place of the requested correlated subspace.
    with open(path, 'w') as fh:
        fh.write('1\n1\n2\ncubic\n0 0 2\n0 0 2\n10\n1\n-2 3\n')
    try:
        dp.read_indmftpr(path)
    except NotImplementedError as error:
        assert 'irreducible' in str(error)
    else:
        raise AssertionError('split correlated shell was silently accepted')

    # The input lacks the cross overlaps needed for two local orbitals;
    # the native implementation rejects this case as well.
    path = os.path.join(tmp, 'multi.almblm')
    with open(path, 'w') as fh:
        fh.write('1\n1\n2\n0\n1\n2\n')
    try:
        dp.read_almblm(path, dict(nsort=1, lmax=0, mult=[1]))
    except NotImplementedError as error:
        assert 'local orbital' in str(error)
    else:
        raise AssertionError('multiple local orbitals were silently accepted')

print('wien2k_dmftproj_common: ok')
