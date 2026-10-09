################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.symqmc generator (triqs_dftkit.wien2k.symqmc)
# reproduces the dmftproj Fortran output: regenerate case.symqmc in Python from
# case.dmftsym + case.indmftpr + case.struct, convert, and compare the resulting
# HDF5 to the reference produced from the Fortran dmftproj symqmc.

from triqs_dftkit.wien2k.symqmc import write_symqmc
from triqs_dftkit.wien2k import Converter
from triqs.utility.h5diff import h5diff
import triqs.utility.mpi as mpi

write_symqmc('CaOs2')

Converter = Converter(filename='CaOs2')
Converter.hdf_file = 'wien2k_symqmc_python.out.h5'
Converter.convert_dft_input()

if mpi.is_master_node():
    h5diff('wien2k_symqmc_python.out.h5', 'wien2k_soc_convert.ref.h5')
