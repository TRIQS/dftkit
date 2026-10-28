################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.symqmc generator on the paths the spin-diagonal
# cubic test (wien2k_symqmc_python) does not exercise:
#
#   - a spin-mixing fromfile basis (the |j, m_j> basis), which uses the full
#     2(2l+1) spinor representation and the spinor time-reversal operator;
#   - an l = 0 (s) correlated shell, the 2x2 spin phase block.
#   - cubic p harmonics, including the complex basis change of magnetic operations.
#
# These reuse the CaOs2 symmetry input (16 operations, 8 of them time-reversal).
# The reference is the matrix data of the dmftproj case.symqmc for the same
# input, stored in double precision and compared at 1e-11; each Python spinor
# matrix is separately checked for unitarity.

import gzip
import os
import shutil
import tempfile
import numpy as np

from triqs_dftkit.wien2k import symqmc

HERE = os.path.dirname(os.path.abspath(__file__))


def _matrix_data(symqmc_path, nsym, natom):
    """Flat array of the matrix blocks, skipping header/perm/timeflag lines."""
    toks = iter(open(symqmc_path).read().split())
    assert int(next(toks)) == nsym and int(next(toks)) == natom
    for _ in range(nsym * natom + nsym):     # perm rows + the time-reversal flags
        next(toks)
    return np.array([float(x) for x in toks])


def _run(indmftpr, extra=(), ref_gz=None):
    """Generate case.symqmc in a scratch dir from the shared CaOs2 symmetry
    input and the given indmftpr, returning (matrix data, shells, ops, so)."""
    tmp = tempfile.mkdtemp()
    try:
        for src, dst in (('CaOs2.dmftsym', 'case.dmftsym'),
                         ('CaOs2.struct', 'case.struct'),
                         (indmftpr, 'case.indmftpr')) + extra:
            shutil.copy(os.path.join(HERE, src), os.path.join(tmp, dst))
        case = os.path.join(tmp, 'case')
        symqmc.write_symqmc(case)
        nsym, ops = symqmc._read_dmftsym(case + '.dmftsym')
        shells, so = symqmc._read_correlated_shells(
            case + '.indmftpr', case + '.struct')
        natom = len(ops[0]['perm'])
        if ref_gz is not None:
            got = open(case + '.symqmc').read().split()
            with gzip.open(os.path.join(HERE, ref_gz), 'rt') as fh:
                ref = fh.read().split()
            header = 2 + nsym * natom + nsym
            assert len(got) == len(ref)
            assert [int(x) for x in got[:header]] == [int(x) for x in ref[:header]]
            np.testing.assert_allclose(np.array(got[header:], dtype=float),
                                       np.array(ref[header:], dtype=float),
                                       atol=1e-11, rtol=0, equal_nan=False)
        return _matrix_data(case + '.symqmc', nsym, natom), shells, ops, so
    finally:
        shutil.rmtree(tmp)


def _timeinv(op, so):
    k = op['krotm']
    return 1 if (so and k[0, 0] * k[1, 1] - k[0, 1] * k[1, 0] < 0.0) else 0


def _check(indmftpr, ref_npy, extra=(), unit_tol=1e-7):
    data, shells, ops, so = _run(indmftpr, extra)
    ref = np.load(os.path.join(HERE, ref_npy))
    assert data.shape == ref.shape, (data.shape, ref.shape)
    assert np.max(np.abs(data - ref)) < 1e-11, np.max(np.abs(data - ref))
    _check_unitary(shells, ops, so, unit_tol)


def _check_unitary(shells, ops, so, unit_tol):
    for op in ops:
        for sh in shells:
            mat = symqmc._shell_matrix(op, sh, _timeinv(op, so), bool(so))
            dev = np.max(np.abs(mat @ np.conj(mat.T) - np.eye(mat.shape[0])))
            assert dev < unit_tol, ('not unitary', dev)


# spin-mixing |j, m_j> basis
data, shells, _, _ = _run('CaOs2_jbasis.indmftpr',
                          (('jbasis_d.dat', 'jbasis_d.dat'),))
assert shells[0]['mixing'] and shells[0]['P'].shape == (10, 10)
_check('CaOs2_jbasis.indmftpr', 'wien2k_symqmc_jbasis.ref.npy',
       (('jbasis_d.dat', 'jbasis_d.dat'),))

# l = 0 (s) correlated shell
_, shells, _, _ = _run('CaOs2_l0.indmftpr')
assert shells[0]['l'] == 0
_check('CaOs2_l0.indmftpr', 'wien2k_symqmc_l0.ref.npy', unit_tol=1e-12)

# Cubic p basis: full native header, all operations and both correlated atoms.
_, shells, ops, so = _run('CaOs2_pcubic.indmftpr', ref_gz='CaOs2_pcubic.symqmc.gz')
assert all(sh['l'] == 1 and not sh['mixing'] for sh in shells)
_check_unitary(shells, ops, so, unit_tol=1e-12)

print('wien2k_symqmc_mixing: ok')
