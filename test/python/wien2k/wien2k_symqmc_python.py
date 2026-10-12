################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.symqmc generator (triqs_dftkit.wien2k.symqmc)
# reproduces the shipped dmftproj Fortran output within its precision limits.
# Generate in a scratch directory to preserve the SOC converter's fixture.

from pathlib import Path
import shutil
import tempfile
import numpy as np

from triqs_dftkit.wien2k.symqmc import write_symqmc

HERE = Path(__file__).resolve().parent


def read_symqmc(path):
    with open(path) as f:
        tokens = f.read().split()
    # This fixture has 16 operations, 3 atoms and two 10-dimensional shells.
    header_size = 2 + 16 * 3 + 16
    header = np.array(tokens[:header_size], dtype=int)
    blocks = np.array(tokens[header_size:], dtype=float).reshape(16, 2, 2, 10, 10)
    return header, blocks[:, :, 0] + 1j * blocks[:, :, 1]


with tempfile.TemporaryDirectory() as tmp:
    case = Path(tmp) / 'CaOs2'
    for suffix in ('dmftsym', 'indmftpr', 'struct'):
        shutil.copyfile(HERE / f'CaOs2.{suffix}', f'{case}.{suffix}')
    write_symqmc(str(case))
    header, matrices = read_symqmc(f'{case}.symqmc')

ref_header, ref_matrices = read_symqmc(HERE / 'CaOs2.symqmc')
np.testing.assert_array_equal(header, ref_header)
# set_ang_trans.f reads 0.70710678 through single-precision CMPLX; magnetic
# operations use the rounded basis twice. outputqmc.f also evaluates the spin
# phase in single precision. The complex matrix residual is about 8.6e-8.
np.testing.assert_allclose(matrices, ref_matrices, atol=1e-7, rtol=0)
for operation in matrices:
    for matrix in operation:
        np.testing.assert_allclose(matrix @ matrix.conj().T, np.eye(10),
                                   atol=1e-12, rtol=0)

print('wien2k_symqmc_python: ok')
