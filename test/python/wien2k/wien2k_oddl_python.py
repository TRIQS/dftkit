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
# promoted to a complex p or f shell. The Loewdin overlap must be full rank:
#   p (l=1): the standard almblm, wide window -2..3 (6 spin-orbitals/atom).
#   f (l=3): Os carries no f weight near E_F, so the projector is rank-deficient
#     on the standard bands. CaOs2_fshell.almblm is recomputed with a high band
#     cutoff (lapw1/lapwso EMAX raised on the converged density); the window
#     3..10 then sits on high-energy bands with real l=3 character, full rank
#     (overlap condition number ~1.3).

import gzip
import os
import shutil
import tempfile

from triqs_dftkit.wien2k import ctqmcout, symqmc

HERE = os.path.dirname(os.path.abspath(__file__))


def _gunzip(src, dst):
    with gzip.open(src, 'rb') as fi, open(dst, 'wb') as fo:
        shutil.copyfileobj(fi, fo)


def _compare(got_path, ref_gz):
    got = open(got_path).read().split()
    with gzip.open(ref_gz, 'rt') as fh:
        ref = fh.read().split()
    assert len(got) == len(ref), (len(got), len(ref))
    max_int, max_float = 0, 0.0
    for a, b in zip(got, ref):
        if a.lstrip('-').isdigit() and b.lstrip('-').isdigit():
            max_int = max(max_int, abs(int(a) - int(b)))
        else:
            max_float = max(max_float, abs(float(a) - float(b)))
    assert max_int == 0, f'integer field differs (max {max_int})'
    assert max_float < 1e-11, f'float field differs (max {max_float:.2e})'


def _check(case, almblm_case, exts):
    tmp = tempfile.mkdtemp()
    try:
        shutil.copy(os.path.join(HERE, f'{case}.indmftpr'),
                    os.path.join(tmp, f'{case}.indmftpr'))
        for ext in ('struct', 'dmftsym'):
            shutil.copy(os.path.join(HERE, f'CaOs2.{ext}'),
                        os.path.join(tmp, f'{case}.{ext}'))
        for spin in ('up', 'dn'):
            _gunzip(os.path.join(HERE, f'{almblm_case}.almblm{spin}.gz'),
                    os.path.join(tmp, f'{case}.almblm{spin}'))
        case_path = os.path.join(tmp, case)
        for ext, writer in exts:
            writer(case_path)
            _compare(case_path + '.' + ext,
                     os.path.join(HERE, f'{case}.{ext}.gz'))
    finally:
        shutil.rmtree(tmp)


# p shell: standard SO almblm, full rank on the wide window.
_check('CaOs2_pshell', 'CaOs2',
       [('ctqmcout', ctqmcout.write_ctqmcout), ('symqmc', symqmc.write_symqmc)])

# f shell: high-band-cutoff almblm, full rank on the high-energy window.
_check('CaOs2_fshell', 'CaOs2_fshell',
       [('ctqmcout', ctqmcout.write_ctqmcout), ('symqmc', symqmc.write_symqmc)])

print('wien2k_oddl_python: ok')
