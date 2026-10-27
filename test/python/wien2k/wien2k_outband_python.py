################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.outband generator (triqs_dftkit.wien2k.outband)
# reproduces the dmftproj Fortran band-structure projectors to machine
# precision: regenerate case.outband from the band-mode almblm/struct/dmftsym,
# the k-path (klist_band) and the band indmftpr, and compare to the native
# dmftproj -band reference (full-precision cubic-d template).
#
# case.outband feeds convert_bands_input. It is the band analog of
# case.ctqmcout: the same correlated-shell projector construction (raw Alm/Clm
# projector -> rot_projectmat -> fromfile/cubic transmat -> Loewdin
# orthonormalization) plus the raw Theta projectors of the parproj path, written
# in the outband layout. The CaOs2_band fixture is SP+SO with the spin-mixing
# fromfile (|j, m_j>) Os d shell (jbasis_d.dat); the wide window (-2.0, 3.0)
# keeps 50 bands > 20 correlated spin-orbitals, so the Loewdin overlap is full
# rank and the projectors match to machine precision.
#
# Band-mode specifics (dmftproj.f:557-581): nkband, the number of k-points along
# the plotted path, is read from CaOs2_band.klist_band; the Fermi energy is the
# last line of CaOs2_band.indmftpr (the band almblm keeps a placeholder eferm
# record that the reader skips).
#
# CaOs2_band_cubic runs the same band inputs with a cubic (spin-diagonal) Os d
# shell and an uncorrelated included Ca d shell, covering the non-mixing
# projector path; its reference is native dmftproj -band output.

import gzip
import os
import shutil
import tempfile
import numpy as np

from triqs_dftkit.wien2k import outband

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
        if a.lstrip('-').isdigit() and b.lstrip('-').isdigit():
            assert int(a) == int(b), (a, b)
        elif _isfloat(a) and _isfloat(b):
            floats.append((float(a), float(b)))
        else:
            assert a == b, (a, b)        # k-label tokens
    a, b = np.array(floats).T
    np.testing.assert_allclose(a, b, atol=1e-11, rtol=0, equal_nan=False)


def _isfloat(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def _check(case, inputs=None):
    inputs = inputs or case
    tmp = tempfile.mkdtemp()
    try:
        shutil.copy(os.path.join(HERE, f'{case}.indmftpr'),
                    os.path.join(tmp, f'{case}.indmftpr'))
        shutil.copy(os.path.join(HERE, f'{inputs}.klist_band'),
                    os.path.join(tmp, f'{case}.klist_band'))
        shutil.copy(os.path.join(HERE, 'jbasis_d.dat'),
                    os.path.join(tmp, 'jbasis_d.dat'))
        for ext in ('struct', 'dmftsym'):
            shutil.copy(os.path.join(HERE, f'CaOs2.{ext}'),
                        os.path.join(tmp, f'{case}.{ext}'))
        for ext in ('almblmup', 'almblmdn'):
            _gunzip(os.path.join(HERE, f'{inputs}.{ext}.gz'),
                    os.path.join(tmp, f'{case}.{ext}'))
        outband.write_outband(os.path.join(tmp, case))
        with gzip.open(os.path.join(HERE, f'{case}.outband.gz'), 'rt') as fh:
            _compare(os.path.join(tmp, f'{case}.outband'), fh.readlines())
        with open(os.path.join(tmp, f'{case}.klist_band'), 'w') as fh:
            fh.write('GAMMA            0    0    0    8  2.0\nX                4    0    0    8  2.0\nEND\n')
        assert outband.read_k_list(os.path.join(tmp, f'{case}.klist_band')) == (2, [(1, 'GAMMA'), (2, 'X')])
        try:
            outband.write_outband(os.path.join(tmp, case))
        except ValueError as exc:
            assert 'k-point counts differ' in str(exc)
        else:
            raise AssertionError('inconsistent band-path dimensions must be rejected')
    finally:
        shutil.rmtree(tmp)


_check('CaOs2_band')
_check('CaOs2_band_cubic', inputs='CaOs2_band')

print('wien2k_outband_python: ok')
