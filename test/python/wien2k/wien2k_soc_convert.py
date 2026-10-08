# Spin-orbit + spin-polarized Wien2k converter test.
# See CaOs2.README.md for the native producer inputs and band window.

from h5 import HDFArchive
import numpy as np
from triqs.utility.h5diff import h5diff
import triqs.utility.mpi as mpi
from triqs_dftkit.wien2k import Converter

converter = Converter(filename='CaOs2')
converter.hdf_file = 'wien2k_soc_convert.out.h5'
converter.convert_dft_input()

if mpi.is_master_node():
    with HDFArchive(converter.hdf_file, 'r') as archive:
        data = archive['dft_input']
        assert data['SP'] == data['SO'] == 1
        assert data['proj_mat'].shape[1] == 1
        assert data['n_k'] == 1
        assert list(data['rot_mat_time_inv']) == [0, 0]
        assert data['n_corr_shells'] == 2
        assert [shell['dim'] for shell in data['corr_shells']] == [10, 10]
        np.testing.assert_allclose(data['bz_weights'].sum(), 1.0, atol=1e-14)
        for ik in range(data['n_k']):
            n_band = data['n_orbitals'][ik, 0]
            blocks = [data['proj_mat'][ik, 0, ish, :shell['dim'], :n_band]
                      for ish, shell in enumerate(data['corr_shells'])]
            # Test the full correlated subspace, including cross-shell overlaps.
            projector = np.vstack(blocks)
            assert n_band >= projector.shape[0]
            np.testing.assert_allclose(projector @ projector.conj().T,
                                       np.eye(projector.shape[0]),
                                       atol=1e-6, rtol=0)
        symmetry = archive[converter.symmcorr_subgrp]
        assert symmetry['n_symm'] == 16
        assert list(symmetry['time_inv']) == [0] * 8 + [1] * 8
        assert sum(perm[1:] == [3, 2] for perm in symmetry['perm']) == 8
        for operation in symmetry['mat']:
            for matrix in operation:
                np.testing.assert_allclose(matrix @ matrix.conj().T,
                                           np.eye(matrix.shape[0]),
                                           atol=1e-6, rtol=0)
    h5diff(converter.hdf_file, 'wien2k_soc_convert.ref.h5', precision=1e-12)
