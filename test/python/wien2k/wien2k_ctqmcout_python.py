################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.ctqmcout generator (triqs_dftkit.wien2k.ctqmcout)
# reproduces the dmftproj Fortran correlated-shell projectors to machine
# precision: regenerate case.ctqmcout from the almblm/struct/dmftsym + a wide
# energy window and compare to the native dmftproj reference.
#
# The window (-2.0, 3.0) keeps 50 bands, more than the 20 correlated spin-orbitals
# of the two SOC Os d shells, so the Loewdin overlap O is full rank. The
# references are native dmftproj output with double-precision casts and a
# full-precision cubic-d template.
#
# CaOs2_jbasis_full exercises the spin-mixing fromfile (|j, m_j>) path: the Os d
# shell is a 2(2l+1)=10 spinor basis (jbasis_d.dat), so the projector, Rloc
# rotrep and complex-harmonics transform blocks are the full mixing matrices.
# Same wide window for the full-rank Loewdin.
#
# CaOs2_mode1_full and CaOs2_mode2_full exercise the band-index projection modes
# (set_projections.f:70-88). Mode 1 ("-2.0 3.0 1") scans all spins/k for the
# global band-index window; mode 2 ("43 92 2") takes the indices straight from
# the window line. Both resolve to bands 43..92 (50 bands > 20 correlated
# spin-orbitals), so the Loewdin overlap stays full rank and the projectors
# match the native dmftproj to machine precision. CaOs2_mode2_below ("0 92 2")
# puts b_bot below nbmin, so no bands lie below the window and qbbot is zero.

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


def _compare(got_path, ref_lines):
    got = open(got_path).read().split()
    ref = ''.join(ref_lines).split()
    assert len(got) == len(ref), (len(got), len(ref))
    floats = []
    for a, b in zip(got, ref):
        if b.lstrip('-').isdigit():
            assert int(a) == int(b), (a, b)
        else:
            floats.append((float(a), float(b)))
    a, b = np.array(floats).T
    np.testing.assert_allclose(a, b, atol=1e-11, rtol=0, equal_nan=False)


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
_check('CaOs2_mode2_below')           # proj_mode 2 with b_bot below nbmin

print('wien2k_ctqmcout_python: ok')

# A singular overlap must fail clearly rather than writing NaN projectors.
try:
    ctqmcout._sqrt_inv(np.ones((2, 2), dtype=complex))
except ValueError as error:
    assert 'rank deficient' in str(error)
else:
    raise AssertionError('singular overlap was accepted')

# Rank is invariant under a uniform projector scaling. Small, well-conditioned
# overlaps must remain valid and normalize to the identity.
for scale in (1., 1e-20):
    overlap = np.diag([scale, 2 * scale]).astype(complex)
    inverse_root = ctqmcout._sqrt_inv(overlap)
    np.testing.assert_allclose(inverse_root @ overlap @ inverse_root.conj().T,
                               np.eye(2), atol=1e-14, rtol=0)
