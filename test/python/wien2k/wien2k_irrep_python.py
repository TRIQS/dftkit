################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.* generators when case.indmftpr selects
# individual irreducible representations of a correlated shell, against native
# dmftproj output in double precision (see CaOs2.README.md, "Irrep selection
# references"). The Os d shell of CaOs2 splits into eg (2) and t2g (3) in the
# cubic basis and into j=3/2 (4) and j=5/2 (6) in jbasis_d.dat.
#
#   CaOs2_irrep_t2g          t2g, no SO, spin-polarized
#   CaOs2_irrep_eg           eg, no SO, spin-polarized
#   CaOs2_irrep_t2g_paramag  t2g, no SO, paramagnetic
#   CaOs2_irrep_so           t2g selected under SO in the cubic basis: the
#                            whole shell is kept, also for a blank selection
#                            line, which the Fortran '(14i1)' read takes as 00
#   CaOs2_irrep_j52          j=5/2 of the spin-mixing basis, SO
#   CaOs2_band_j52           as CaOs2_irrep_j52 along the band path (outband)

import gzip
import os
import shutil
import tempfile

import numpy as np

from triqs_dftkit.wien2k import ctqmcout, outband, symqmc
from triqs_dftkit.wien2k._dmftproj import read_indmftpr

HERE = os.path.dirname(os.path.abspath(__file__))

WRITERS = dict(ctqmcout=ctqmcout.write_ctqmcout, symqmc=symqmc.write_symqmc,
               outband=outband.write_outband)


def _gunzip(src, dst):
    with gzip.open(src, 'rb') as fi, open(dst, 'wb') as fo:
        shutil.copyfileobj(fi, fo)


def _compare(got_path, ref_path):
    got = open(got_path).read().split()
    with gzip.open(ref_path, 'rt') as fh:
        ref = fh.read().split()
    assert len(got) == len(ref), (got_path, len(got), len(ref))
    floats = []
    for a, b in zip(got, ref):
        if b.lstrip('-').isdigit():
            assert int(a) == int(b), (got_path, a, b)
        else:
            try:
                floats.append((float(a), float(b)))
            except ValueError:
                assert a == b, (got_path, a, b)    # k-path labels
    actual, expected = np.array(floats).T
    np.testing.assert_allclose(actual, expected, atol=1e-11, rtol=0,
                               equal_nan=False, err_msg=got_path)


def _setup(tmp, indmftpr, almblm, paramagnetic=False):
    case = os.path.join(tmp, 'case')
    shutil.copy(os.path.join(HERE, f'{indmftpr}.indmftpr'), case + '.indmftpr')
    shutil.copy(os.path.join(HERE, 'CaOs2.struct'), case + '.struct')
    shutil.copy(os.path.join(HERE, 'CaOs2.dmftsym'), case + '.dmftsym')
    shutil.copy(os.path.join(HERE, 'jbasis_d.dat'), os.path.join(tmp, 'jbasis_d.dat'))
    shutil.copy(os.path.join(HERE, 'CaOs2_band.klist_band'), case + '.klist_band')
    if paramagnetic:
        _gunzip(os.path.join(HERE, f'{almblm}.almblmup.gz'), case + '.almblm')
    else:
        for ext in ('almblmup', 'almblmdn'):
            _gunzip(os.path.join(HERE, f'{almblm}.{ext}.gz'), f'{case}.{ext}')
    return case


def _check(reference, indmftpr, almblm, exts, paramagnetic=False, correps_line=None):
    tmp = tempfile.mkdtemp()
    try:
        case = _setup(tmp, indmftpr, almblm, paramagnetic)
        if correps_line is not None:
            lines = open(case + '.indmftpr').read().splitlines()
            lines[9] = correps_line
            open(case + '.indmftpr', 'w').write('\n'.join(lines) + '\n')
        for ext in exts:
            WRITERS[ext](case)
            _compare(f'{case}.{ext}', os.path.join(HERE, f'{reference}.{ext}.gz'))
    finally:
        shutil.rmtree(tmp)


def _rejects(correps_lines, message):
    """An indmftpr whose irrep selection native dmftproj refuses must raise."""
    tmp = tempfile.mkdtemp()
    try:
        case = _setup(tmp, 'CaOs2_irrep_t2g', 'CaOs2_noso')
        lines = open(case + '.indmftpr').read().splitlines()
        lines[8:10] = correps_lines
        open(case + '.indmftpr', 'w').write('\n'.join(lines) + '\n')
        try:
            ctqmcout.write_ctqmcout(case)
        except ValueError as error:
            assert message in str(error), error
        else:
            raise AssertionError(f'{correps_lines} was accepted')
    finally:
        shutil.rmtree(tmp)


both = ('ctqmcout', 'symqmc')
_check('CaOs2_irrep_t2g', 'CaOs2_irrep_t2g', 'CaOs2_noso', both)
_check('CaOs2_irrep_eg', 'CaOs2_irrep_eg', 'CaOs2_noso', both)
_check('CaOs2_irrep_t2g_paramag', 'CaOs2_irrep_t2g', 'CaOs2_noso', both,
       paramagnetic=True)
_check('CaOs2_irrep_so', 'CaOs2_irrep_so', 'CaOs2', both)
_check('CaOs2_irrep_so', 'CaOs2_irrep_so', 'CaOs2', both, correps_line='')
_check('CaOs2_irrep_j52', 'CaOs2_irrep_j52', 'CaOs2', both)
_check('CaOs2_band_j52', 'CaOs2_band_j52', 'CaOs2_band', ('outband',))

_rejects(['0 0 3 0', '011'], 'the basis has 2')
_rejects(['0 0 2 0', '00'], 'no irrep')
_rejects(['0 0 2 0', ''], 'no irrep')
_rejects(['0 0 2 0', '11'], 'one irrep per shell')
_rejects(['0 0 6 0', '000001'], 'not possible')

# Selection lines are read one physical line each, also after a blank one.
_lines = ['1', '1', '3', 'cubic', '0 2 2 0', '0 1 2 0', '', '01', '0', '-2.0 3.0']
assert read_indmftpr('case.indmftpr', lines=_lines)['sorts'][0]['correps'] == {1: [0], 2: [0, 1]}
_lines[6:8] = ['1', '']
assert read_indmftpr('case.indmftpr', lines=_lines)['sorts'][0]['correps'] == {1: [1], 2: [0, 0]}

print('wien2k_irrep_python: ok')
