################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# End-to-end non-spin-orbit converter test. The Python dmftproj driver writes
# every case.* file, the Wien2k converter ingests them, and the HDF5 is compared
# at machine precision (1e-12) against the reference built from the Fortran
# dmftproj output. Beyond the SO converter test, this exercises convert_dft_input
# on the non-SO path (SP=1, SO=0: the .oubwinup/.oubwindn band windows) and
# convert_parproj_input (case.parproj + case.sympar) end to end.
#
# CaOs2_noso_partial is the spin-polarized non-SO CaOs2 with an uncorrelated Ca
# d shell added to the correlated Os d shells; the wide window (-2..3) keeps the
# Loewdin overlap full rank. The reference is precision-fixed dmftproj run
# without -so, so the agreement is double precision, not the 1e-6 SO floor.

import gzip
import os
import shutil
import tempfile

from triqs.utility.h5diff import h5diff
import triqs.utility.mpi as mpi
from triqs_dftkit.wien2k import Converter
from triqs_dftkit.wien2k.dmftproj import run_dmftproj

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = 'CaOs2_noso_partial'

tmp = tempfile.mkdtemp()
shutil.copy(os.path.join(HERE, f'{CASE}.indmftpr'),
            os.path.join(tmp, f'{CASE}.indmftpr'))
for ext in ('struct', 'dmftsym', 'outputs'):
    shutil.copy(os.path.join(HERE, f'CaOs2.{ext}'),
                os.path.join(tmp, f'{CASE}.{ext}'))
for spin in ('up', 'dn'):
    with gzip.open(os.path.join(HERE, f'CaOs2_noso.almblm{spin}.gz'), 'rb') as fi, \
            open(os.path.join(tmp, f'{CASE}.almblm{spin}'), 'wb') as fo:
        shutil.copyfileobj(fi, fo)

run_dmftproj(os.path.join(tmp, CASE))

conv = Converter(filename=os.path.join(tmp, CASE))
conv.hdf_file = os.path.join(tmp, 'wien2k_noso_convert.out.h5')
conv.convert_dft_input()
conv.convert_parproj_input()

if mpi.is_master_node():
    h5diff(os.path.join(tmp, 'wien2k_noso_convert.out.h5'),
           os.path.join(HERE, 'wien2k_noso_convert.ref.h5'), precision=1e-12)
    print('wien2k_noso_convert: ok')
