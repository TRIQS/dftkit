
################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################


"""
Collinear SrVO3 tests: two correlated sites (V t2g and V eg) and uncorrelated
orbitals (O p), two impurities
"""

import os
from triqs_dftkit.wannier90 import Converter
from triqs.utility import h5diff
from triqs.utility import mpi

subfolder = 'w90_convert/'

# Test 1: Wannier basis (bloch_basis=False)
seedname = subfolder+'SrVO3_col'
converter = Converter(seedname=seedname, hdf_filename=seedname+'_wannierbasis.out.h5',
                               rot_mat_type='hloc_diag', bloch_basis=False)
converter.convert_dft_input()

if mpi.is_master_node():
    h5diff.h5diff(seedname+'_wannierbasis.out.h5', seedname+'_wannierbasis.ref.h5')

# Test 2: Bloch basis (bloch_basis=True)
# Need to temporarily move OUTCAR and LOCPROJ files for converter to find them
if mpi.is_master_node():
    os.rename(seedname + '.OUTCAR', subfolder + 'OUTCAR')
    os.rename(seedname + '.LOCPROJ', subfolder + 'LOCPROJ')
mpi.barrier()

try:
    converter = Converter(seedname=seedname, hdf_filename=seedname+'_blochbasis.out.h5',
                                   rot_mat_type='hloc_diag', bloch_basis=True)
    converter.convert_dft_input()
finally:
    if mpi.is_master_node():
        os.rename(subfolder + 'OUTCAR', seedname + '.OUTCAR')
        os.rename(subfolder + 'LOCPROJ', seedname + '.LOCPROJ')

if mpi.is_master_node():
    h5diff.h5diff(seedname+'_blochbasis.out.h5', seedname+'_blochbasis.ref.h5')

# Test 3: the two modes must describe the same operator
#
# h5diff only compares each archive with its own frozen reference, so an error
# present in both modes would pass. The two are related by an identity that can
# be checked directly: the Bloch basis stores the Kohn-Sham eigenvalues plus the
# projectors P(k) = (U_dis(k) U(k))^dag, the Wannier basis stores H_W(k), and
#
#     P(k) H_bloch(k) P(k)^dag = U^dag eps(k) U = H_W(k)
#
# for every k, also with band disentanglement (25 bands, 14 Wannier functions
# here). This is what lets postprocessing that works in the Wannier basis, such
# as solid_dmft's plot_correlated_bands, use either archive.
#
# Two things have to be taken care of: the two modes subtract a different Fermi
# energy from the diagonal, and the Bloch mode takes its k mesh from
# seedname_u.mat, so the k-points come in a different order.
if mpi.is_master_node():
    import numpy as np
    from h5 import HDFArchive

    def read_archive(filename):
        with HDFArchive(filename, 'r') as ar:
            dft_input = ar['dft_input']
            fermi_energy = 0.0
            if 'dft_misc_input' in ar and 'dft_fermi_energy' in ar['dft_misc_input']:
                fermi_energy = float(ar['dft_misc_input']['dft_fermi_energy'])
            return (dft_input['hopping'], dft_input['proj_mat'], dft_input['n_orbitals'],
                    dft_input['kpts'], [sh['dim'] for sh in dft_input['corr_shells']], fermi_energy)

    h_wan, _, _, kpts_wan, dims, ef_wan = read_archive(seedname+'_wannierbasis.out.h5')
    h_blo, proj_blo, n_orbitals_blo, kpts_blo, dims_blo, ef_blo = read_archive(seedname+'_blochbasis.out.h5')

    assert dims == dims_blo, 'correlated shells differ between the two modes'

    def kpt_key(kpt):
        return tuple(np.round(np.mod(kpt + 1e-9, 1.0), 6))

    wan_of_key = {kpt_key(kpt): ik for ik, kpt in enumerate(kpts_wan)}
    kpt_pairs = [(ik_blo, wan_of_key[kpt_key(kpt)]) for ik_blo, kpt in enumerate(kpts_blo)
                 if kpt_key(kpt) in wan_of_key]
    assert len(kpt_pairs) == len(kpts_blo), \
        f'only {len(kpt_pairs)} of {len(kpts_blo)} k-points of the Bloch archive found in the Wannier one'

    max_deviation = 0.0
    for ik_blo, ik_wan in kpt_pairs:
        n_bands = n_orbitals_blo[ik_blo, 0]
        for ish, dim in enumerate(dims):
            proj = proj_blo[ik_blo, 0, ish, :dim, :n_bands]
            downfolded = np.dot(proj, np.dot(h_blo[ik_blo, 0, :n_bands, :n_bands],
                                             proj.conjugate().transpose()))
            # the two modes subtract a different Fermi energy from the diagonal
            downfolded += np.eye(dim) * (ef_blo - ef_wan)
            offset = sum(dims[:ish])
            max_deviation = max(max_deviation, np.max(np.abs(
                downfolded - h_wan[ik_wan, 0, offset:offset+dim, offset:offset+dim])))

    mpi.report('Bloch vs Wannier basis: max_k |P(k) H(k) P(k)^dag - H_W(k)| = '
               '{:.2e} eV over {} k-points'.format(max_deviation, len(kpt_pairs)))
    # limited by the precision seedname_hr.dat is written with
    assert max_deviation < 1e-3, \
        f'Bloch and Wannier basis disagree by {max_deviation:.3e} eV'
