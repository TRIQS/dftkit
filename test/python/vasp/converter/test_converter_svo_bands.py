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
