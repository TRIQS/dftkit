################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.* generators on tellurium (P3_121), a
# non-centrosymmetric cell with three equivalent atoms. The local frame of each
# atom is the representative rotloc composed with the symmetry operation that
# maps the representative onto it, which the centrosymmetric fixtures never
# exercise. Inputs: Te.almblm{up,dn}.gz (one k-point, see Te.README.md) with
# Te.struct and Te.dmftsym. References are native dmftproj output, double
# precision, with the full-precision cubic templates.
#
#   Te              p correlated, SO (dmftproj -sp -so), window -0.5 0.5 Ry
#   Te_sp           s and p correlated, SO, window -0.8 0.5 Ry
#   Te_sp_mag       as Te_sp, with Te_mag.dmftsym: the operations of Te.dmftsym
#                   reordered so that a two-fold axis, a magnetic operation under
#                   SO, maps the representative onto each atom
#   Te_noso         s included, p correlated, no SO (dmftproj -sp)
#   Te_noso_paramag as Te_noso from Te.almblmup alone (dmftproj without flags)

import gzip
import os
import shutil
import tempfile

import numpy as np

from triqs_dftkit.wien2k import ctqmcout, symqmc, sympar, parproj

HERE = os.path.dirname(os.path.abspath(__file__))

WRITERS = dict(ctqmcout=ctqmcout.write_ctqmcout, symqmc=symqmc.write_symqmc,
               sympar=sympar.write_sympar, parproj=parproj.write_parproj)


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
            floats.append((float(a), float(b)))
    actual, expected = np.array(floats).T
    np.testing.assert_allclose(actual, expected, atol=1e-11, rtol=0,
                               equal_nan=False, err_msg=got_path)


def _check(reference, indmftpr, exts, dmftsym='Te', paramagnetic=False):
    tmp = tempfile.mkdtemp()
    case = os.path.join(tmp, 'case')
    try:
        shutil.copy(os.path.join(HERE, f'{indmftpr}.indmftpr'), case + '.indmftpr')
        shutil.copy(os.path.join(HERE, 'Te.struct'), case + '.struct')
        shutil.copy(os.path.join(HERE, f'{dmftsym}.dmftsym'), case + '.dmftsym')
        if paramagnetic:
            _gunzip(os.path.join(HERE, 'Te.almblmup.gz'), case + '.almblm')
        else:
            for ext in ('almblmup', 'almblmdn'):
                _gunzip(os.path.join(HERE, f'Te.{ext}.gz'), f'{case}.{ext}')
        for ext in exts:
            WRITERS[ext](case)
            _compare(f'{case}.{ext}', os.path.join(HERE, f'{reference}.{ext}.gz'))
    finally:
        shutil.rmtree(tmp)


_check('Te', 'Te', ('ctqmcout', 'symqmc', 'parproj'))
_check('Te_sp', 'Te_sp', ('ctqmcout', 'symqmc', 'parproj'))
_check('Te_sp_mag', 'Te_sp', ('ctqmcout', 'parproj'), dmftsym='Te_mag')
_check('Te_noso', 'Te_noso', ('ctqmcout',))
_check('Te_noso_paramag', 'Te_noso', ('ctqmcout', 'symqmc', 'sympar', 'parproj'),
       paramagnetic=True)

print('wien2k_te_python: ok')
