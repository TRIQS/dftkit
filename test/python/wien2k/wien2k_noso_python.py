################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.* generators on the non-spin-orbit path (ifSO=0),
# the common DMFT case (e.g. SrVO3). dmftproj writes (2l+1) matrices there, not
# the 2*(2l+1) spinor matrices of the SO path: the orthonormalization is
# per-spin (orthogonal_wannier), the symmetry and Rloc matrices are the bare
# (2l+1) representation, there is one density block per input spin channel,
# and srot%timeinv is always false (setsym.f sets it only under SP+SO).
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

import numpy as np

from triqs_dftkit.wien2k import ctqmcout, symqmc, sympar, parproj

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
    actual, expected = np.array(floats).T
    np.testing.assert_allclose(actual, expected, atol=1e-11, rtol=0,
                               equal_nan=False)


def _stage(case, paramagnetic=False):
    tmp = tempfile.mkdtemp()
    shutil.copy(os.path.join(HERE, f'{case}.indmftpr'),
                os.path.join(tmp, f'{case}.indmftpr'))
    for ext in ('struct', 'dmftsym'):
        shutil.copy(os.path.join(HERE, f'CaOs2.{ext}'),
                    os.path.join(tmp, f'{case}.{ext}'))
    if paramagnetic:
        _gunzip(os.path.join(HERE, 'CaOs2_noso.almblmup.gz'),
                os.path.join(tmp, f'{case}.almblm'))
    else:
        for ext in ('almblmup', 'almblmdn'):
            _gunzip(os.path.join(HERE, f'CaOs2_noso.{ext}.gz'),
                    os.path.join(tmp, f'{case}.{ext}'))
    return tmp


def _check_output(case, ext, writer, paramagnetic=False):
    tmp = _stage(case, paramagnetic)
    try:
        writer(os.path.join(tmp, case))
        reference = case + '_paramag' if paramagnetic else case
        with gzip.open(os.path.join(HERE, f'{reference}.{ext}.gz'), 'rt') as fh:
            _compare(os.path.join(tmp, f'{case}.{ext}'), fh.readlines())
    finally:
        shutil.rmtree(tmp)


# Correlated-only: the projector and correlated-shell symmetry files.
_check_output('CaOs2_noso', 'ctqmcout', ctqmcout.write_ctqmcout)
_check_output('CaOs2_noso', 'symqmc', symqmc.write_symqmc)

# Partial: the included-shell symmetry and partial-projector files.
_check_output('CaOs2_noso_partial', 'sympar', sympar.write_sympar)
_check_output('CaOs2_noso_partial', 'parproj', parproj.write_parproj)

# A single non-magnetic input carries one density block and no spin flag.
_check_output('CaOs2_noso_partial', 'parproj', parproj.write_parproj,
              paramagnetic=True)

# Shells with different l on one atom remain independent in the non-SO path.
orbs = [dict(atom=1, sort=1, l=l) for l in (0, 1)]
raw = {(1, l): [np.eye(2 * l + 1) * value for value in (1., 2.)]
       for l in (0, 1)}
representations = {(l, 1, 0): (np.eye(2 * l + 1), False, 0.) for l in (0, 1)}
result = parproj._symmetrize_densmat_nonso(
    orbs, {}, [dict(perm=[1])], 1, representations, raw)
for shell in raw:
    for actual, expected in zip(result[shell], raw[shell]):
        np.testing.assert_allclose(actual, expected, atol=1e-14, rtol=0)

# SO metadata requires both spin channels, including when a single-channel
# almblm happens to be present. Never export a half-sized SO projector.
tmp = _stage('CaOs2_noso')
try:
    case = os.path.join(tmp, 'CaOs2_noso')
    os.rename(case + '.almblmup', case + '.almblm')
    os.unlink(case + '.almblmdn')
    with open(case + '.indmftpr') as fh:
        lines = fh.read().splitlines()
    lines[-2] = '1'
    with open(case + '.indmftpr', 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    for writer in (ctqmcout.write_ctqmcout, parproj.write_parproj):
        try:
            writer(case)
        except ValueError as error:
            assert 'both up and down' in str(error)
        else:
            raise AssertionError('single-channel SO projectors were accepted')
finally:
    shutil.rmtree(tmp)

for module in (symqmc, sympar):
    try:
        module._shell_matrix({}, dict(l=1, mixing=True, P=np.eye(6)), 0, False)
    except NotImplementedError as error:
        assert 'spin-mixing' in str(error)
    else:
        raise AssertionError('non-SO symmetry accepted a spinor basis')

print('wien2k_noso_python: ok')
