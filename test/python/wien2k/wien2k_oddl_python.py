################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the generators on the odd-l shells (p, l=1; f, l=3). The d-only fixtures
# (l=2) cannot catch the improper-rotation parity: the Wigner-D parity factor is
# (-1)^l, which is +1 for even l and only flips sign for odd l. case.dmftsym
# stores the proper part of each operation in krotm (det(krotm)=+1 even for
# improper ops) and the true parity in iprop; dmat must use iprop. The p- and
# f-shell symqmc here pin that sign.
#
# Both are the spin-orbit CaOs2 cell (two equivalent Os atoms) with the Os shell
# promoted to a complex p or f shell.
#   p (l=1): the standard almblm, wide window -2..3 (6 spin-orbitals/atom), so
#     the projector path (ctqmcout) is checked as well.
#   f (l=3): complex-basis symmetry, plus all generators in Wien2k's mm2
#     cubic basis over the full retained band range, with and without SO.

import gzip
import os
import shutil
import tempfile

import numpy as np

from triqs_dftkit.wien2k import ctqmcout, parproj, symqmc, sympar
from triqs_dftkit.wien2k._dmftproj import basis_dreps, reptrans

HERE = os.path.dirname(os.path.abspath(__file__))


def _gunzip(src, dst):
    with gzip.open(src, 'rb') as fi, open(dst, 'wb') as fo:
        shutil.copyfileobj(fi, fo)


def _compare(got_path, ref_gz):
    got = open(got_path).read().split()
    with gzip.open(ref_gz, 'rt') as fh:
        ref = fh.read().split()
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


def _check(case, almblm_case, exts, paramagnetic=False):
    tmp = tempfile.mkdtemp()
    try:
        shutil.copy(os.path.join(HERE, f'{case}.indmftpr'),
                    os.path.join(tmp, f'{case}.indmftpr'))
        for ext in ('struct', 'dmftsym'):
            shutil.copy(os.path.join(HERE, f'CaOs2.{ext}'),
                        os.path.join(tmp, f'{case}.{ext}'))
        for spin in (('up',) if paramagnetic else ('up', 'dn')) if almblm_case else ():
            _gunzip(os.path.join(HERE, f'{almblm_case}.almblm{spin}.gz'),
                    os.path.join(tmp, f'{case}.almblm' + ('' if paramagnetic else spin)))
        case_path = os.path.join(tmp, case)
        for ext, writer in exts:
            writer(case_path)
            _compare(case_path + '.' + ext,
                     os.path.join(HERE, f'{case}.{"symqmc" if ext == "sympar" else ext}.gz'))
    finally:
        shutil.rmtree(tmp)


# p shell: standard SO almblm, full rank on the wide window.
_check('CaOs2_pshell', 'CaOs2',
       [('ctqmcout', ctqmcout.write_ctqmcout), ('symqmc', symqmc.write_symqmc)])

# f shell: symmetry output only.
_check('CaOs2_fshell', None, [('symqmc', symqmc.write_symqmc)])

all_outputs = [('ctqmcout', ctqmcout.write_ctqmcout), ('parproj', parproj.write_parproj),
               ('symqmc', symqmc.write_symqmc), ('sympar', sympar.write_sympar)]
_check('CaOs2_fcubic_so', 'CaOs2', all_outputs)
_check('CaOs2_fcubic_sp', 'CaOs2_noso', all_outputs)
_check('CaOs2_fcubic_pm', 'CaOs2_noso', all_outputs, paramagnetic=True)
_check('CaOs2_fcubic_a2', 'CaOs2_noso',
       [('ctqmcout', ctqmcout.write_ctqmcout), ('symqmc', symqmc.write_symqmc)])


def _rejects_irrep(selection):
    """An irrep selection the Os site symmetry mixes with the rest of the f
    shell must be refused by every generator that writes irrep blocks."""
    tmp = tempfile.mkdtemp()
    try:
        case = os.path.join(tmp, 'case')
        lines = open(os.path.join(HERE, 'CaOs2_fcubic_a2.indmftpr')).read().splitlines()
        lines[9] = selection
        open(case + '.indmftpr', 'w').write('\n'.join(lines) + '\n')
        for ext in ('struct', 'dmftsym'):
            shutil.copy(os.path.join(HERE, f'CaOs2.{ext}'), f'{case}.{ext}')
        for spin in ('up', 'dn'):
            _gunzip(os.path.join(HERE, f'CaOs2_noso.almblm{spin}.gz'), f'{case}.almblm{spin}')
        for writer in (ctqmcout.write_ctqmcout, symqmc.write_symqmc):
            try:
                writer(case)
            except ValueError as error:
                assert 'not invariant' in str(error), error
            else:
                raise AssertionError(f'{writer.__name__} accepted f irreps {selection}')
    finally:
        shutil.rmtree(tmp)


_rejects_irrep('0010')
_rejects_irrep('0001')

np.testing.assert_allclose(reptrans('cubic', 3) @ reptrans('cubic', 3).conj().T,
                           np.eye(7), atol=1e-14, rtol=0)
assert basis_dreps('cubic', 3) == [2, 1, 2, 2]

print('wien2k_oddl_python: ok')
