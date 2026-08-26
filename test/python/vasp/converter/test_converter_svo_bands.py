import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

from h5 import HDFArchive
import numpy as np

from triqs_dftkit.vasp.plovasp.converter import generate_and_output_as_text
from triqs_dftkit.vasp import Converter
import mytest


class TestConverterSVOBands(mytest.MyTestCase):
    """
    Test conversion of KPOINTS_OPT + LOCPROJ_OPT into dft_bands_input.
    """

    def _check_bands_payload(self, test_file, vasp_dir):
        with HDFArchive(test_file, 'r') as ar:
            assert 'dft_bands_input' in ar, "Missing dft_bands_input group"
            bands = ar['dft_bands_input']

            things = ['n_k', 'n_orbitals', 'proj_mat', 'hopping', 'n_parproj', 'proj_mat_all']
            for it in things:
                assert it in bands, "Missing key in dft_bands_input: %s" % it

            n_k = int(bands['n_k'])
            n_orbitals = bands['n_orbitals']
            proj_mat = bands['proj_mat']
            hopping = bands['hopping']
            n_orb_min = int(np.min(n_orbitals))
            n_orb_max = int(np.max(n_orbitals))

            assert n_k == 200, "Unexpected number of k-points in bands data"
            self.assertEqual(n_orbitals.shape, (200, 1))
            self.assertEqual(n_orb_min, 3)
            self.assertEqual(n_orb_max, 5)
            self.assertEqual(proj_mat.shape[0:4], (200, 1, 1, 3))
            self.assertEqual(proj_mat.shape[4], n_orb_max)
            self.assertEqual(hopping.shape[0:2], (200, 1))
            self.assertEqual(hopping.shape[2], n_orb_max)
            self.assertEqual(hopping.shape[3], n_orb_max)

            # High-symmetry k-path labels read from vaspout.h5 (/input/kpoints_opt).
            # The path is GAMMA-X-M-GAMMA-R with 50 points per segment; the shared
            # endpoints of adjacent segments are collapsed into a single tick.
            assert 'kpts_labels' in bands, "Missing kpts_labels in dft_bands_input"
            assert 'kpts_labels_idx' in bands, "Missing kpts_labels_idx in dft_bands_input"
            self.assertEqual(list(bands['kpts_labels']), ['GAMMA', 'X', 'M', 'GAMMA', 'R'])
            np.testing.assert_array_equal(bands['kpts_labels_idx'], [0, 49, 99, 149, 199])

        with HDFArchive(vasp_dir + 'vaspout.h5', 'r') as src:
            eig = src['results/electron_eigenvalues_kpoints_opt/eigenvalues']
            efermi = float(src['results/electron_dos/efermi'])
        with HDFArchive(test_file, 'r') as ar:
            ib1 = int(ar['dft_misc_input']['band_window'][0][0, 0]) - 1
            expected_h = eig[0, 0, ib1] - efermi

        with HDFArchive(test_file, 'r') as ar:
            h00 = ar['dft_bands_input']['hopping'][0, 0, 0, 0]
            self.assertAlmostEqual(h00.real, expected_h)
            self.assertAlmostEqual(h00.imag, 0.0)

        self._check_locproj_opt_scale(vasp_dir)
        self._check_downfolded_t2g(test_file)

    def _check_locproj_opt_scale(self, vasp_dir):
        """
        Compare the KPOINTS_OPT projectors with the regular-mesh ones at GAMMA,
        the one k-point the line-mode path and the 15x15x15 mesh have in common.

        Both groups must express the same localized orbitals, so at a shared
        k-point sum_orb sum_band |<phi|psi>|^2 has to agree. That sum is
        invariant under the arbitrary band phase and under any unitary mixing
        inside a degenerate multiplet, hence independent of the MPI
        decomposition (verified to 1e-8 relative between 4 and 8 ranks).

        PLOVasp orthonormalizes the projectors, so an overall factor on the raw
        KPOINTS_OPT amplitudes divides out again and is invisible in
        dft_bands_input. Two kinds of error are therefore only detectable here:
        an amplitude that scales with the MPI decomposition, and a LORBIT=14
        "optimal" PAW channel rebuilt from the interpolated KPOINTS_OPT
        k-points instead of kept from the ground-state mesh. Both have occurred
        during development of the VASP side of this interface.
        """
        def to_complex(raw):
            arr = np.asarray(raw)
            # (proj, spin, k, band, 2) -> (proj, k, band), single spin channel
            return arr[:, 0, ..., 0] + 1j * arr[:, 0, ..., 1]

        with HDFArchive(vasp_dir + 'vaspout.h5', 'r') as src:
            proj_opt = to_complex(src['results/locproj_opt/data'])
            proj_mesh = to_complex(src['results/locproj/data'])
            kpts_opt = np.asarray(src['results/electron_eigenvalues_kpoints_opt/kpoint_coords'])
            kpts_mesh = np.asarray(src['results/projectors/kpoints'])

        self.assertEqual(proj_opt.shape, (5, 200, 24))

        ik_opt = int(np.argmin(np.abs(kpts_opt).sum(axis=1)))
        ik_mesh = int(np.argmin(np.abs(kpts_mesh).sum(axis=1)))
        for label, kpt in (('path', kpts_opt[ik_opt]), ('mesh', kpts_mesh[ik_mesh])):
            np.testing.assert_allclose(kpt, 0.0, atol=1e-8,
                                       err_msg='GAMMA not found in the %s k-points' % label)

        weight_opt = (np.abs(proj_opt[:, ik_opt, :]) ** 2).sum()
        weight_mesh = (np.abs(proj_mesh[:, ik_mesh, :]) ** 2).sum()
        np.testing.assert_allclose(weight_opt, weight_mesh, rtol=1e-6)

    def _check_downfolded_t2g(self, test_file):
        """
        Downfold hopping onto the t2g shell with proj_mat and check the SrVO3
        t2g dispersion along GAMMA-X-M-GAMMA-R. This validates the projectors
        themselves (not just shapes): a wrong orbital character, a k/band
        misalignment or a broken label mapping all show up here, while the
        overall projector normalization does not.
        """
        with HDFArchive(test_file, 'r') as ar:
            bands = ar['dft_bands_input']
            proj_mat = bands['proj_mat']
            hopping = bands['hopping']
            n_orbitals = bands['n_orbitals']
            label_idx = np.asarray(bands['kpts_labels_idx'])

        eigs = []
        for ik in range(proj_mat.shape[0]):
            nb = int(n_orbitals[ik, 0])
            p = proj_mat[ik, 0, 0, :, :nb]
            h = hopping[ik, 0, :nb, :nb]
            eigs.append(np.linalg.eigvalsh(p @ h @ p.conj().T))
        eigs = np.array(eigs)

        # t2g bandwidth of SrVO3.
        np.testing.assert_allclose(eigs.max() - eigs.min(), 2.4866, atol=1e-3)

        # Cubic symmetry: threefold degenerate at GAMMA and at R.
        for tick in (0, 3, 4):   # GAMMA, GAMMA (second pass), R
            ik = int(label_idx[tick])
            self.assertAlmostEqual(float(np.ptp(eigs[ik])), 0.0, places=5)

        # Both GAMMA ticks are the same k-point and must give the same bands.
        np.testing.assert_allclose(eigs[int(label_idx[0])], eigs[int(label_idx[3])], atol=1e-6)

        # Band bottom at GAMMA, band top at R.
        np.testing.assert_allclose(eigs[int(label_idx[0])].mean(), -0.9745, atol=1e-3)
        np.testing.assert_allclose(eigs[int(label_idx[4])].mean(), 1.5121, atol=1e-3)
        self.assertAlmostEqual(float(eigs.min()), float(eigs[int(label_idx[0])].min()), places=5)
        self.assertAlmostEqual(float(eigs.max()), float(eigs[int(label_idx[4])].max()), places=5)

    def test_convert_svo_bands_auto(self):
        vasp_dir = _rpath + 'svo/kpoints_opt/'

        generate_and_output_as_text(vasp_dir + 'plo.cfg', vasp_dir)

        test_file = _rpath + 'svo_bands_auto.test.h5'
        converter = Converter(filename=vasp_dir + 'plo_full', hdf_filename=test_file)

        converter.convert_dft_input()

        self._check_bands_payload(test_file, vasp_dir)

    def test_convert_svo_bands_explicit(self):
        vasp_dir = _rpath + 'svo/kpoints_opt/'

        generate_and_output_as_text(vasp_dir + 'plo.cfg', vasp_dir)

        test_file = _rpath + 'svo_bands_explicit.test.h5'
        converter = Converter(filename=vasp_dir + 'plo_full', hdf_filename=test_file)

        converter.convert_dft_input()
        converter.convert_bands_input(cfg_filename=vasp_dir + 'plo.cfg')

        self._check_bands_payload(test_file, vasp_dir)


if __name__ == '__main__':
    import unittest
    unittest.main(verbosity=2, buffer=False)
