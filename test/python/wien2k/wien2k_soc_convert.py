# Spin-orbit + spin-polarized Wien2k converter test.
# See CaOs2.README.md for the native producer inputs and band window.

from h5 import HDFArchive
import numpy as np
import triqs.utility.mpi as mpi
from triqs_dftkit.wien2k import Converter


def assert_archives_close(a, b, path=''):
    """Recursive comparison at absolute tolerance 1e-12.

    triqs.utility.h5diff compares arrays at 1e-6 whatever precision is
    passed, so it cannot check this tolerance.
    """
    if hasattr(a, 'keys'):
        assert set(a.keys()) == set(b.keys()), path
        for key in a.keys():
            assert_archives_close(a[key], b[key], f'{path}/{key}')
    elif isinstance(a, (list, tuple)):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            assert_archives_close(x, y, f'{path}[{i}]')
    elif isinstance(a, (np.ndarray, np.number, int, float, complex)):
        np.testing.assert_allclose(a, b, atol=1e-12, rtol=0, err_msg=path)
    else:
        assert a == b, path


converter = Converter(filename='CaOs2')
converter.hdf_file = 'wien2k_soc_convert.out.h5'
converter.convert_dft_input()
converter.convert_parproj_input()

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

        parproj = archive[converter.parproj_subgrp]
        assert list(parproj['n_parproj']) == [2, 2]
        assert list(parproj['rot_mat_all_time_inv']) == [0, 0]
        n_band = data['n_orbitals'][0, 0]
        for ish in range(data['n_shells']):
            np.testing.assert_array_equal(parproj['rot_mat_all'][ish],
                                          data['rot_mat'][ish])
            dens_below = parproj['dens_mat_below'][0][ish]
            np.testing.assert_allclose(dens_below, dens_below.conj().T,
                                       atol=1e-12, rtol=0)
        # The two Os shells are symmetry-equivalent.
        for ir in range(2):
            norms = [np.linalg.norm(parproj['proj_mat_all'][0, 0, ish, ir, :, :n_band])
                     for ish in range(data['n_shells'])]
            np.testing.assert_allclose(norms[0], norms[1], rtol=1e-8)

        # Partial shells coincide with the correlated shells, so sympar
        # must reproduce the symqmc symmetry data.
        symmpar = archive[converter.symmpar_subgrp]
        assert symmpar['n_symm'] == symmetry['n_symm']
        assert list(symmpar['time_inv']) == list(symmetry['time_inv'])
        for key in ('perm', 'mat', 'mat_tinv'):
            for a, b in zip(symmpar[key], symmetry[key]):
                np.testing.assert_array_equal(np.array(a), np.array(b))

    with HDFArchive(converter.hdf_file, 'r') as out, \
            HDFArchive('wien2k_soc_convert.ref.h5', 'r') as ref:
        assert_archives_close(out, ref)
