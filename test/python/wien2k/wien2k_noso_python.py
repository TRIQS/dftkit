################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.* generators on the non-spin-orbit path (ifSO=0),
# the common DMFT case (e.g. SrVO3). dmftproj writes (2l+1) matrices there, not
# the 2*(2l+1) spinor matrices of the SO path: the orthonormalization is
# per-spin (orthogonal_wannier), the symmetry and Rloc matrices are the bare
# (2l+1) representation, the density matrix splits into two independent spin
# blocks, and srot%timeinv is always false (setsym.f sets it only under SP+SO).
#
# CaOs2_noso is the spin-polarized, non-SO CaOs2 (indmftpr SO flag 0, wide
# window -2.0 3.0): two equivalent Os d shells, cubic harmonics. The wide window
# keeps the Loewdin overlap full rank. ctqmcout and symqmc are checked against
# the precision-fixed dmftproj reference for this case.
#
# CaOs2_noso_partial adds an uncorrelated Ca d shell, so sympar and parproj have
# an included-but-not-correlated shell to symmetrize. Same struct/dmftsym/almblm
# as the correlated-only case. References built with the precision-fixed
# dmftproj run without -so.

import gzip
import os
import shutil
import tempfile

from triqs_dftkit.wien2k import ctqmcout, symqmc, sympar, parproj

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
        else:
            max_float = max(max_float, abs(float(a) - float(b)))
    assert max_int == 0, f'integer field differs (max {max_int})'
    assert max_float < 1e-11, f'float field differs (max {max_float:.2e})'


def _stage(case):
    tmp = tempfile.mkdtemp()
    shutil.copy(os.path.join(HERE, f'{case}.indmftpr'),
                os.path.join(tmp, f'{case}.indmftpr'))
    for ext in ('struct', 'dmftsym'):
        shutil.copy(os.path.join(HERE, f'CaOs2.{ext}'),
                    os.path.join(tmp, f'{case}.{ext}'))
    for ext in ('almblmup', 'almblmdn'):
        _gunzip(os.path.join(HERE, f'CaOs2_noso.{ext}.gz'),
                os.path.join(tmp, f'{case}.{ext}'))
    return tmp


def _check_output(case, ext, writer):
    tmp = _stage(case)
    try:
        writer(os.path.join(tmp, case))
        with gzip.open(os.path.join(HERE, f'{case}.{ext}.gz'), 'rt') as fh:
            _compare(os.path.join(tmp, f'{case}.{ext}'), fh.readlines())
    finally:
        shutil.rmtree(tmp)


# Correlated-only: the projector and correlated-shell symmetry files.
_check_output('CaOs2_noso', 'ctqmcout', ctqmcout.write_ctqmcout)
_check_output('CaOs2_noso', 'symqmc', symqmc.write_symqmc)

# Partial: the included-shell symmetry and partial-projector files.
_check_output('CaOs2_noso_partial', 'sympar', sympar.write_sympar)
_check_output('CaOs2_noso_partial', 'parproj', parproj.write_parproj)

print('wien2k_noso_python: ok')
