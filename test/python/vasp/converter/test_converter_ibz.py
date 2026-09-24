
import os
import unittest.mock
import numpy as np
import rpath
_rpath = os.path.dirname(rpath.__file__) + '/'

from h5 import HDFArchive
from triqs_dftkit.vasp.plovasp.converter import generate_and_output_as_text
from triqs_dftkit.vasp import Converter
from triqs_dftkit.vasp import symmetry
import mytest


def local_occupations(h5_file):
    """
    DFT occupation matrix of every correlated shell, sum_k w_k P_k f_k P_k^dag.
    """
    with HDFArchive(h5_file, 'r') as ar:
        dft = ar['dft_input']
        proj_mat, bz_weights, n_orbitals = dft['proj_mat'], dft['bz_weights'], dft['n_orbitals']
        corr_shells = dft['corr_shells']
        f_weights = ar['dft_misc_input']['dft_fermi_weights']

    n_k, n_spin = proj_mat.shape[:2]
    occ = []
    for icrsh, csh in enumerate(corr_shells):
        o = np.zeros((n_spin, csh['dim'], csh['dim']), dtype=complex)
        for ik in range(n_k):
            for isp in range(n_spin):
                nb = n_orbitals[ik, isp]
                P = proj_mat[ik, isp, icrsh, :csh['dim'], :nb]
                o[isp] += bz_weights[ik] * (P * f_weights[ik, isp, :nb]) @ P.conj().T
        occ.append(o)
    return occ


def symmetrize(occ, h5_file):
    """
    Symmetrize per-shell matrices with the dft_symmcorr_input of ``h5_file``,
    following triqs_dft_tools.symmetry.Symmetry, including how the image of a
    shell is found (comparison of the whole orbit dict, 'sort' included), so the
    archive is checked in the form a DMFT code consumes it.
    """
    with HDFArchive(h5_file, 'r') as ar:
        symm = ar['dft_symmcorr_input']
    orbits = symm['orbits']
    symm_occ = [np.zeros_like(o) for o in occ]
    for isym in range(symm['n_symm']):
        for iorb, orb in enumerate(orbits):
            jorb = orbits.index(dict(orb, atom=symm['perm'][isym][orb['atom'] - 1]))
            mat = symm['mat'][isym][iorb]
            symm_occ[jorb] += mat @ occ[iorb] @ mat.conj().T / symm['n_symm']
    return symm_occ


################################################################################
#
# TestConverterIBZ
#
# Converts VASP output to the irreducible-BZ representation with symmetrization
# (use_ibz=True): the data are reduced to the irreducible k-points and the
# correlated-shell symmetry operations are written to dft_symmcorr_input so that
# DMFT can run on the IBZ instead of the unfolded full grid.
#
# Besides the comparison with the reference archive, the symmetrized IBZ
# occupation matrices must reproduce the symmetrized full-grid ones, which only
# holds if the weights, the orbital rotations and the atom permutations are all
# correct. The unsymmetrized full grid agrees only to the symmetry of the VASP
# projectors themselves, hence the looser tolerance for that comparison.
#
################################################################################
class TestConverterIBZ(mytest.MyTestCase):

    def convert(self, name, use_ibz):
        generate_and_output_as_text(_rpath + name + '.cfg', _rpath + name + '/')
        test_file = _rpath + name + ('_ibz' if use_ibz else '_full') + '.test.h5'
        if os.path.exists(test_file): os.remove(test_file)
        converter = Converter(filename=_rpath + name, hdf_filename=test_file)
        converter.convert_dft_input(use_ibz=use_ibz, vasp_h5=_rpath + name + '/vaspout.h5')
        return test_file

    def check_ibz(self, name, tol):
        full_file = self.convert(name, use_ibz=False)
        ibz_file = self.convert(name, use_ibz=True)

        self.assertH5FileEqual(ibz_file, _rpath + name + '_ibz.ref.h5')

        occ_full = local_occupations(full_file)
        occ_ibz = symmetrize(local_occupations(ibz_file), ibz_file)
        occ_full_symm = symmetrize(occ_full, ibz_file)
        self.assertEqual(len(occ_full), len(occ_ibz))
        for icrsh in range(len(occ_full)):
            np.testing.assert_allclose(occ_ibz[icrsh], occ_full_symm[icrsh], rtol=0, atol=1e-8,
                                       err_msg=f"correlated shell {icrsh}")
            np.testing.assert_allclose(occ_ibz[icrsh], occ_full[icrsh], rtol=0, atol=tol,
                                       err_msg=f"correlated shell {icrsh}")

    def test_convert_svo_ibz(self):
        """SrVO3: single V, cubic O_h, t2g (l=2) projection."""
        self.check_ibz('svo', tol=1e-10)

    def test_convert_nio_afm_ibz(self):
        """
        NiO in the AFM cell (non-spin-polarized): full Ni d (l=2) and O p (l=1),
        non-cubic cell, operations that exchange the two O. The full-grid VASP
        projectors are only symmetric to ~1e-5 here.
        """
        self.check_ibz('nio_afm', tol=1e-4)


    def test_rotation_representation(self):
        """D^l(R) is orthogonal and a representation, D(R1 R2) = D(R1) D(R2)."""
        rng = np.random.default_rng(1)
        def random_op():
            q, r = np.linalg.qr(rng.normal(size=(3, 3)))
            return q * np.sign(np.diag(r))       # random O(3) element, proper or improper
        for l in symmetry.SUPPORTED_L:
            for _ in range(5):
                R1, R2 = random_op(), random_op()
                D1, D2 = (symmetry.real_harmonic_rotation(R, l) for R in (R1, R2))
                np.testing.assert_allclose(D1 @ D1.T, np.eye(2 * l + 1), atol=1e-12)
                np.testing.assert_allclose(symmetry.real_harmonic_rotation(R1 @ R2, l), D1 @ D2, atol=1e-12)

    def test_unsupported_shells(self):
        """f shells and non-orthonormal transforms are refused."""
        vasp_h5 = _rpath + 'svo/vaspout.h5'
        f_shell = [{'atom': 2, 'sort': 2, 'l': 3, 'dim': 7, 'SO': 0, 'irep': 0}]
        with self.assertRaises(NotImplementedError):
            symmetry.build_symmcorr(vasp_h5, f_shell, [np.eye(7)], 0, 0)
        d_shell = [{'atom': 2, 'sort': 2, 'l': 2, 'dim': 5, 'SO': 0, 'irep': 0}]
        with self.assertRaises(RuntimeError):
            symmetry.build_symmcorr(vasp_h5, d_shell, [2 * np.eye(5)], 0, 0)
        with self.assertRaises(NotImplementedError):
            symmetry.build_symmcorr(vasp_h5, d_shell, [np.eye(5)], 0, 1)

    def check_fallback(self, name, cfg, **kwargs):
        """use_ibz=None falls back to the full grid, use_ibz=True raises."""
        generate_and_output_as_text(_rpath + cfg + '.cfg', _rpath + name + '/')
        test_file = _rpath + cfg + '_fallback.test.h5'
        if os.path.exists(test_file): os.remove(test_file)
        converter = Converter(filename=_rpath + cfg, hdf_filename=test_file)
        converter.convert_dft_input(vasp_h5=_rpath + name + '/vaspout.h5', **kwargs)
        with HDFArchive(test_file, 'r') as ar:
            self.assertEqual(ar['dft_input']['symm_op'], 0)
            self.assertGreater(ar['dft_input']['n_k'], ar['dft_misc_input']['n_k_ibz'])
        with self.assertRaises(RuntimeError):
            Converter(filename=_rpath + cfg, hdf_filename=test_file).convert_dft_input(
                use_ibz=True, vasp_h5=_rpath + name + '/vaspout.h5', **kwargs)

    def test_fallback_missing_equivalent_atom(self):
        """Only one of the two equivalent O is projected: no IBZ."""
        self.check_fallback('nio_afm', 'nio_afm_one_o')

    def test_fallback_wrong_operations(self):
        """
        Wrong d rotation matrices (z^2 and x^2-y^2 exchanged) pass all structural
        checks but not the comparison with the full grid.
        """
        rotation = symmetry.real_harmonic_rotation
        def wrong_rotation(R, l):
            D = rotation(R, l)
            if l == 2:
                idx = [0, 1, 4, 3, 2]
                D = D[np.ix_(idx, idx)]
            return D
        with unittest.mock.patch.object(symmetry, 'real_harmonic_rotation', wrong_rotation):
            self.check_fallback('nio_afm', 'nio_afm')


if __name__ == '__main__':
    import unittest
    unittest.main(verbosity=2, buffer=False)
