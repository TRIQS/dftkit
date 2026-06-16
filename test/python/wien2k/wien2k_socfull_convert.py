################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# End-to-end spin-orbit converter test at machine precision: the Python driver
# writes every case.* file, the Wien2k converter ingests them, and the HDF5 is
# h5diffed against the Fortran-derived reference at 1e-12. The Loewdin overlap
# must be full rank; the physical narrow window is rank-deficient (the overlap
# has near-null eigenvalues whose eigenvectors are not unique across LAPACK
# builds), so this test uses the wide window -2..3. The narrow physical case is
# locked exactly by wien2k_soc_convert (converter on the Fortran output).
#
# CaOs2_socfull is spin-orbit + spin-polarized CaOs2 with the Os d shell in the
# complex basis (no cubic template), wide window. convert_dft_input only; the
# cell has no uncorrelated shell.

import gzip
import os
import shutil
import tempfile

from triqs.utility.h5diff import h5diff
import triqs.utility.mpi as mpi
from triqs_dftkit.wien2k import Converter
from triqs_dftkit.wien2k.dmftproj import run_dmftproj

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = 'CaOs2_socfull'

tmp = tempfile.mkdtemp()
shutil.copy(os.path.join(HERE, f'{CASE}.indmftpr'),
            os.path.join(tmp, f'{CASE}.indmftpr'))
for ext in ('struct', 'dmftsym', 'outputs'):
    shutil.copy(os.path.join(HERE, f'CaOs2.{ext}'),
                os.path.join(tmp, f'{CASE}.{ext}'))
for spin in ('up', 'dn'):
    with gzip.open(os.path.join(HERE, f'CaOs2.almblm{spin}.gz'), 'rb') as fi, \
            open(os.path.join(tmp, f'{CASE}.almblm{spin}'), 'wb') as fo:
        shutil.copyfileobj(fi, fo)

run_dmftproj(os.path.join(tmp, CASE))

conv = Converter(filename=os.path.join(tmp, CASE))
conv.hdf_file = os.path.join(tmp, 'wien2k_socfull_convert.out.h5')
conv.convert_dft_input()

if mpi.is_master_node():
    h5diff(os.path.join(tmp, 'wien2k_socfull_convert.out.h5'),
           os.path.join(HERE, 'wien2k_socfull_convert.ref.h5'), precision=1e-12)
    print('wien2k_socfull_convert: ok')
