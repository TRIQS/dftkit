
import os
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

from triqs_dftkit.vasp.plovasp.converter import generate_and_output_as_text
from triqs_dftkit.vasp import Converter
import mytest

################################################################################
#
# TestConverterIBZ
#
# Converts VASP output to the irreducible-BZ representation with symmetrization
# (use_ibz=True): the data are reduced to the irreducible k-points and the
# correlated-shell symmetry operations are written to dft_symmcorr_input so that
# DMFT can run on the IBZ instead of the unfolded full grid.
#
################################################################################
class TestConverterIBZ(mytest.MyTestCase):

    def test_convert_svo_ibz(self):
        """SrVO3: single V, cubic O_h, t2g (l=2) projection."""
        generate_and_output_as_text(_rpath + 'svo.cfg', _rpath + 'svo/')

        test_file = _rpath + 'svo_ibz.test.h5'
        converter = Converter(filename=_rpath + 'svo', hdf_filename=test_file)
        converter.convert_dft_input(use_ibz=True,
                                    vasp_h5=_rpath + 'svo/vaspout.h5',
                                    plo_cfg=_rpath + 'svo.cfg')

        self.assertH5FileEqual(test_file, _rpath + 'svo_ibz.ref.h5')

    def test_convert_nio_afm_ibz(self):
        """AFM NiO: two equivalent Ni, non-cubic cell, full d (l=2) projection."""
        generate_and_output_as_text(_rpath + 'nio_afm.cfg', _rpath + 'nio_afm/')

        test_file = _rpath + 'nio_afm_ibz.test.h5'
        converter = Converter(filename=_rpath + 'nio_afm', hdf_filename=test_file)
        converter.convert_dft_input(use_ibz=True,
                                    vasp_h5=_rpath + 'nio_afm/vaspout.h5',
                                    plo_cfg=_rpath + 'nio_afm.cfg')

        self.assertH5FileEqual(test_file, _rpath + 'nio_afm_ibz.ref.h5')


if __name__ == '__main__':
    import unittest
    unittest.main(verbosity=2, buffer=False)
