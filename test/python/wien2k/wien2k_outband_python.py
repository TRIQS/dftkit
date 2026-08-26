################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.outband generator (triqs_dftkit.wien2k.outband)
# reproduces the dmftproj Fortran band-structure projectors to machine
# precision: regenerate case.outband from the band-mode almblm/struct/dmftsym,
# the k-path (klist_band) and the band indmftpr, and compare to the
# precision-fixed dmftproj -band reference.
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

import gzip
import os
import shutil
import tempfile

from triqs_dftkit.wien2k import outband

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
        elif _isfloat(a) and _isfloat(b):
            max_float = max(max_float, abs(float(a) - float(b)))
        else:
            assert a == b, (a, b)        # k-label tokens
    assert max_int == 0, f'integer field differs (max {max_int})'
    assert max_float < 1e-11, f'float field differs (max {max_float:.2e})'


def _isfloat(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def _check(case):
    tmp = tempfile.mkdtemp()
    try:
        for ext in ('indmftpr', 'klist_band'):
            shutil.copy(os.path.join(HERE, f'{case}.{ext}'),
                        os.path.join(tmp, f'{case}.{ext}'))
        shutil.copy(os.path.join(HERE, 'jbasis_d.dat'),
                    os.path.join(tmp, 'jbasis_d.dat'))
        for ext in ('struct', 'dmftsym'):
            shutil.copy(os.path.join(HERE, f'CaOs2.{ext}'),
                        os.path.join(tmp, f'{case}.{ext}'))
        for ext in ('almblmup', 'almblmdn'):
            _gunzip(os.path.join(HERE, f'{case}.{ext}.gz'),
                    os.path.join(tmp, f'{case}.{ext}'))
        outband.write_outband(os.path.join(tmp, case))
        with gzip.open(os.path.join(HERE, f'{case}.outband.gz'), 'rt') as fh:
            _compare(os.path.join(tmp, f'{case}.outband'), fh.readlines())
    finally:
        shutil.rmtree(tmp)


_check('CaOs2_band')

print('wien2k_outband_python: ok')
