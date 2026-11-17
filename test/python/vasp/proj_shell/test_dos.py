r"""
Tests of ProjectorShell.density_of_states, density_matrix and local_hamiltonian
on small synthetic inputs.
"""
from types import SimpleNamespace

import numpy as np
from triqs_dftkit.vasp.plovasp.proj_shell import ProjectorShell
import mytest

# one tetrahedron spanned by the four (IBZ) k-points 1..4, multiplicity 1
ITET = np.array([[1, 1, 2, 3, 4]])
VOLT = 1.0 / 6.0
EMESH = np.linspace(-2.0, 2.0, 801)
# two dispersive bands, one below and one above the Fermi level (0)
EIGS = np.array([[-1.2, 0.7], [-1.0, 0.9], [-0.9, 1.0], [-0.7, 1.2]])


def make_shell(proj_win, ib_win, ib_min, ib_max):
    shell = ProjectorShell.__new__(ProjectorShell)
    shell.proj_win = proj_win
    shell.ib_win = ib_win
    shell.ib_min = ib_min
    shell.ib_max = ib_max
    return shell


def make_el_struct(eigs, nkibz, ksymmap=None):
    kmesh = {'itet': ITET, 'volt': VOLT, 'nkibz': nkibz, 'ksymmap': ksymmap}
    return SimpleNamespace(eigvals=eigs[:, :, None], efermi=0.0, kmesh=kmesh)


def integrate(dos, mask):
    return np.trapezoid(dos[mask], EMESH[mask], axis=0)


class TestDensityOfStates(mytest.MyTestCase):
    """
    Scenarios:
    - the band index of the projectors counts from the first band of the window at each
      k-point, so a band that leaves the window must not shift the weights to another band
    - with symmetry the orbital weights are averaged over the star of each IBZ k-point
    """

    def test_k_dependent_window(self):
        # band 0 is in the window at k-points 0 and 1 only; the projector is on band 1 everywhere
        nk, nlm = 4, 1
        ib_win = np.array([[[0, 1]], [[0, 1]], [[1, 1]], [[1, 1]]])
        proj_win = np.zeros((1, 1, nk, nlm, 2), dtype=complex)
        proj_win[0, 0, :2, 0, 1] = 1.0   # band 1 is the second band of the window
        proj_win[0, 0, 2:, 0, 0] = 1.0   # band 1 is the first band of the window
        shell = make_shell(proj_win, ib_win, ib_min=0, ib_max=1)

        dos = shell.density_of_states(make_el_struct(EIGS, nkibz=nk), EMESH)[:, 0, 0, 0]

        # reference: the same projectors with both bands in the window at all k-points
        proj_ref = np.zeros_like(proj_win)
        proj_ref[0, 0, :, 0, 1] = 1.0
        ib_ref = np.zeros_like(ib_win)
        ib_ref[:, 0, 1] = 1
        dos_ref = make_shell(proj_ref, ib_ref, 0, 1).density_of_states(make_el_struct(EIGS, nkibz=nk), EMESH)[:, 0, 0, 0]

        self.assertTrue(np.all(dos[EMESH < 0] == 0.0))
        self.assertTrue(integrate(dos, EMESH > 0) > 0.0)
        self.assertTrue(np.allclose(dos, dos_ref))

    def test_star_average(self):
        # 8 full-grid k-points, two per IBZ k-point; the two members of a star carry
        # the weight in different orbitals, as related by a symmetry operation
        nkibz, nlm = 4, 2
        ksymmap = np.array([0, 1, 2, 3, 0, 1, 2, 3])
        eigs_full = EIGS[ksymmap]
        ib_win = np.zeros((8, 1, 2), dtype=int)
        ib_win[:, 0, 1] = 1
        proj_win = np.zeros((1, 1, 8, nlm, 2), dtype=complex)
        proj_win[0, 0, :4, 0, 0] = 1.0
        proj_win[0, 0, 4:, 1, 0] = 1.0
        shell = make_shell(proj_win, ib_win, ib_min=0, ib_max=1)

        dos = shell.density_of_states(make_el_struct(eigs_full, nkibz, ksymmap), EMESH)[:, 0, 0, :]

        # reference: no symmetry, the weight in a single orbital at the 4 k-points
        ib_ref = np.zeros((4, 1, 2), dtype=int)
        ib_ref[:, 0, 1] = 1
        proj_ref = np.zeros((1, 1, 4, 1, 2), dtype=complex)
        proj_ref[0, 0, :, 0, 0] = 1.0
        dos_ref = make_shell(proj_ref, ib_ref, 0, 1).density_of_states(make_el_struct(EIGS, nkibz=4), EMESH)[:, 0, 0, 0]

        self.assertTrue(np.allclose(dos[:, 0], dos[:, 1]))
        self.assertTrue(np.allclose(dos[:, 0] + dos[:, 1], dos_ref))
        self.assertTrue(integrate(dos_ref, EMESH < 0) > 0.0)


class TestWindowBandIndex(mytest.MyTestCase):
    """
    Scenario:
    - the band window starts at different bands at different k-points: density matrix and
      local Hamiltonian pair each projector with its own band (proj_win counts from ib1)
    """

    def test_density_matrix_and_local_hamiltonian(self):
        # band 0 (occupied) is in the window at k-points 0 and 1 only; the projector
        # is on band 1 (empty) everywhere
        nk = 4
        ib_win = np.array([[[0, 1]], [[0, 1]], [[1, 1]], [[1, 1]]])
        proj_win = np.zeros((1, 1, nk, 1, 2), dtype=complex)
        proj_win[0, 0, :2, 0, 1] = 1.0
        proj_win[0, 0, 2:, 0, 0] = 1.0
        shell = make_shell(proj_win, ib_win, ib_min=0, ib_max=1)
        shell.nc_flag = False

        el_struct = make_el_struct(EIGS, nkibz=nk)
        el_struct.kmesh['kweights'] = np.full(nk, 1.0 / nk)
        el_struct.ferw = np.zeros((1, nk, 2))
        el_struct.ferw[0, :, 0] = 1.0

        occ, overlap = shell.density_matrix(el_struct)
        loc_ham = shell.local_hamiltonian(el_struct)

        self.assertTrue(np.allclose(occ, 0.0))
        self.assertTrue(np.allclose(overlap, 1.0))
        self.assertTrue(np.allclose(loc_ham, EIGS[:, 1].mean()))

