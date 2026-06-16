################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.parproj generator (triqs_dftkit.wien2k.parproj)
# reproduces the dmftproj Fortran output. parproj holds the partial projectors
# and density matrices for the INCLUDED orbitals; the CaOs2_partial fixture adds
# an uncorrelated Ca d-shell. Inputs are the shared CaOs2 almblm/struct/dmftsym
# renamed to the CaOs2_partial case.
#
# Floats compared at 1e-6 (dmftproj single-precision basis-transform floor).

import gzip
import os
import shutil
import tempfile

from triqs_dftkit.wien2k import parproj

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


def _check(case):
    tmp = tempfile.mkdtemp()
    try:
        shutil.copy(os.path.join(HERE, f'{case}.indmftpr'),
                    os.path.join(tmp, f'{case}.indmftpr'))
        for ext in ('struct', 'dmftsym'):
            shutil.copy(os.path.join(HERE, f'CaOs2.{ext}'),
                        os.path.join(tmp, f'{case}.{ext}'))
        for ext in ('almblmup', 'almblmdn'):
            _gunzip(os.path.join(HERE, f'CaOs2.{ext}.gz'),
                    os.path.join(tmp, f'{case}.{ext}'))
        parproj.write_parproj(os.path.join(tmp, case))
        with gzip.open(os.path.join(HERE, f'{case}.parproj.gz'), 'rt') as fh:
            _compare(os.path.join(tmp, f'{case}.parproj'), fh.readlines())
    finally:
        shutil.rmtree(tmp)


_check('CaOs2_partial')

print('wien2k_parproj_python: ok')
