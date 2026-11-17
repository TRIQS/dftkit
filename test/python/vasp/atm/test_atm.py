
import os

import numpy as np
from triqs_dftkit.vasp.plovasp.lintetra import dos_tetra_weights_3d
import mytest

################################################################################
#
# TestProjectorShell
#
################################################################################
class TestProjectorShell(mytest.MyTestCase):
    """
    Class:

    ProjectorShell(sh_pars, proj_raw)

    Scenarios:
    - **if** a correct input is given **compare** output arrays
    """
# Scenario 1
    def test_example(self):
        eigs = np.array([-1.5, -1.309017, -1.0, -0.5])
        en = -0.55
        itt = np.array([[1, 0, 1, 2, 3]]).T

        res = dos_tetra_weights_3d(eigs, en, itt)[:, 0]

        r_should = np.zeros(4)
        r_should[0] = 0.000309016992226;
        r_should[1] = 0.000381966005939;
        r_should[2] = 0.000618033984453;
        r_should[3] = 0.017232002550965;

        self.assertEqual(res, r_should)

# Scenario 2
    def test_energy_array(self):
        rng = np.random.default_rng(1)
        n_k, n_tetra = 30, 200
        eigs = rng.normal(size=n_k)
        eigs[:6] = 0.3   # degenerate corners: uniform and regularized branches
        itt = np.vstack([rng.integers(1, 7, n_tetra), rng.integers(0, n_k, (4, n_tetra))])
        itt[1:, :10] = rng.integers(0, 6, (4, 10))
        # include energies exactly at corner values, where the cases share boundaries
        energies = np.concatenate([np.linspace(-3, 3, 101), eigs[:15]])

        res = dos_tetra_weights_3d(eigs, energies, itt)
        res_single = np.array([dos_tetra_weights_3d(eigs, en, itt) for en in energies])

        self.assertEqual(res.shape, (len(energies), 4, n_tetra))
        self.assertEqual(dos_tetra_weights_3d(eigs, 0.1, itt).shape, (4, n_tetra))
        self.assertTrue(np.array_equal(res, res_single))
