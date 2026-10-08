# ============================================================================
# Charge self-consistent (CSC) DFT+DMFT for SrVO3 with Wien2k + TRIQS/modest
# ============================================================================

import numpy as np

import triqs.utility.mpi as mpi

from triqs.gfs import MeshImFreq

from triqs_cthyb import solve_generic, TailFitParams

import triqs_modest as tm
from triqs_modest.utils import Checkpointer, IterationData
from triqs_modest.dft_driver import DftDriver

from triqs_dftkit.wien2k import Driver as Wien2kDriver


# --- physical and run parameters ---------------------------------------------
seedname = "SrVO3"
U, J     = 4.5, 0.65        # eV
Up       = U - 2 * J        # rotationally invariant Kanamori
beta     = 20.0             # 1/eV
n_iw     = 1000

n_csc_loops        = 50     # outer DFT charge updates
n_dmft_loops       = 1      # DMFT iterations per charge update
n_dmft_loops_first = 5      # ... except the first, where Sigma is converged
                            # against the unmodified DFT charge density
dc_method = "cHeld"         # Held's DC, the usual choice for a t2g-only shell

mesh = MeshImFreq(beta=beta, statistic="Fermion", n_iw=n_iw)

# --- solver setup ------------------------------------------------------------
solver_params = dict(
    n_tau=10 * n_iw,
    length_cycle=500,
    n_cycles=int(7e6 / mpi.size),
    n_warmup_cycles=int(1e4),
)
tail_fit_params = TailFitParams(fit_min_w=20, fit_max_w=24, fit_max_moment=4)


def dmft_cycle(obe, target_density, Sigma_imp_dyn, Sigma_imp_hf, Sigma_dc, n_loops):
    """Run `n_loops` DMFT iterations at fixed H(k), then prepare the DFT feedback.

    The run-level constants (embedding `E`, `h_int`, `dc`, `deg_blocks`, mesh and
    solver parameters) are read from module scope; everything that changes from
    one CSC iteration to the next is passed in explicitly.

    After the last impurity solve the double counting is recomputed from the new
    Gimp and mu is re-found, so that the returned N_k is consistent with the
    self-energy and DC that go with it.

    Returns
    -------
    (IterationData, float, numpy.ndarray)
        The iteration data (mu, self-energies, DC, Gimp/Gloc/Delta), the
        impurity interaction energy minus the DC energy in eV, and the
        band-basis charge density correction N_k[k, sigma, nu, nu'].
    """
    if n_loops < 1:
        raise ValueError("n_loops must be at least 1: the double counting, Eint-Edc "
                         "and N_k below are all built from the impurity solution")

    epsilon_d = E.extract(tm.impurity_levels(obe))[0]

    for n in range(n_loops):
        Sigma_imp_hf_m_dc = [hf - sig_dc for (hf, sig_dc) in zip(Sigma_imp_hf, Sigma_dc)]

        Sigma_C = E.embed([Sigma_imp_dyn], [Sigma_imp_hf_m_dc])                                    # Embed self-energy

        mu = tm.find_chemical_potential(target_density, obe, *Sigma_C, verbosity=False)            # Find chemical potential

        Gloc = E.extract(tm.gloc(obe, mu, *Sigma_C))[0]                                            # Compute Gloc

        ed = [(eps - mu * np.eye(eps.shape[0]) - sig_dc).real for (eps, sig_dc) in zip(epsilon_d, Sigma_dc)]

        Delta = tm.symmetrize(tm.hybridization(ed, Gloc, Sigma_imp_dyn, Sigma_imp_hf), deg_blocks)  # Compute Delta

        res = solve_generic(Delta, ed, h_int, postprocess=tail_fit_params, **solver_params)         # Solve the impurity problem

        Sigma_imp_dyn = tm.symmetrize(res.Sigma_dynamic, deg_blocks)                                # Update Sigma_imp
        Sigma_imp_hf  = tm.symmetrize(res.Sigma_HartreeFock, deg_blocks)
        Sigma_imp_iw  = tm.symmetrize(res.Sigma_iw, deg_blocks)
        Gimp          = tm.symmetrize(res.G_iw, deg_blocks)

        mpi.report(f"  [dmft {n + 1}/{n_loops}] mu = {mu:.6f}  "
                   f"n_loc = {Gloc.total_density().real:.6f}  "
                   f"n_imp = {Gimp.total_density().real:.6f}")

    # --- double counting and the charge density correction ---
    Sigma_dc  = dc.dc_self_energy(Gimp)
    Eint      = 0.5 * np.real((Sigma_imp_iw * Gimp).total_density())
    Eint_m_dc = Eint - dc.dc_energy(Gimp)

    Sigma_C = E.embed([Sigma_imp_dyn], [[hf - sig_dc for (hf, sig_dc) in zip(Sigma_imp_hf, Sigma_dc)]])
    mu   = tm.find_chemical_potential(target_density, obe, *Sigma_C, verbosity=False)
    N_k  = tm.charge_density_correction(obe, mu, *Sigma_C)

    it_data = IterationData(mu=mu,
                            Sigma_imp_list=[Sigma_imp_dyn],
                            Sigma_hartree_list=[Sigma_imp_hf],
                            Sigma_dc_list=[Sigma_dc],
                            Gimp_list=[Gimp], Gloc_list=[Gloc], Delta_list=[Delta])
    return it_data, Eint_m_dc, N_k


# --- DFT driver --------------------------------------------------------------
# The Wien2k driver owns the lapw0 -> lapw1 -> lapw2 -> lcore -> mixer chain in
# python, replacing run_lapw, so that the charge correction can be injected
# between cycles. run_initial_stage() converges the SCF (reusing a converged
# case.scf if there is one), runs dmftproj and converts to SrVO3.h5.
driver = DftDriver(Wien2kDriver(seedname=seedname, ecut=1e-3, ccut=1e-3))

target_density, obe = driver.one_body_elements_from_dft()
mpi.report(f"target_density= {target_density}")
mpi.report(obe)

# --- embedding: one impurity, three degenerate t2g orbitals ------------------
E = tm.make_embedding(obe.C_space)
mpi.report(E.description(True))

# --- interaction and double counting ----------------------------------------
h_int = tm.make_kanamori(E.sigma_names, E.imp_decomposition(0), U, Up, J)
dc    = tm.DcSolver("NonPolarized", dc_method, U, J)

# --- DFT-only pass: block structure and the initial DC -----------------------
mu_dft = tm.find_chemical_potential(target_density, obe, beta, verbosity=False)
Gdft   = E.extract(tm.gloc(mesh, obe, mu_dft))[0]
mpi.report(f"mu_dft= {mu_dft:.6f}  n_dft= {Gdft.total_density().real:.6f}")

deg_blocks = tm.analyze_degenerate_blocks(Gdft)
mpi.report(f"degenerate blocks= {deg_blocks}")

# --- checkpoint: restart if there is something to restart from ---------------
ckpt = Checkpointer(f"svo_csc_beta{beta}_U{U}_J{J}.ckpt")

if (prev := ckpt.restart()):
    Sigma_imp_dyn, Sigma_imp_hf = prev.Sigma_imp_list[0], prev.Sigma_hartree_list[0]
    Sigma_dc = prev.Sigma_dc_list[0]
    mpi.report(f"restarting from checkpoint at iteration {len(ckpt)}")
else:
    Sigma_dc = dc.dc_self_energy(Gdft)
    Sigma_imp_dyn, Sigma_imp_hf = E.make_zero_imp_self_energies(mesh)[0]
    for ibl in range(len(Sigma_imp_hf)):
        Sigma_imp_hf[ibl] += Sigma_dc[ibl]

# --- DFT + DMFT loop ---------------------------------------------------------
for it in range(len(ckpt), n_csc_loops):
    mpi.report(f"\n=== DFT+DMFT iteration {it + 1}/{n_csc_loops} ===")

    n_loops = n_dmft_loops_first if it == 0 else n_dmft_loops
    it_data, Eint_m_dc, N_k = dmft_cycle(obe, target_density, Sigma_imp_dyn, Sigma_imp_hf, Sigma_dc, n_loops)

    Sigma_imp_dyn = it_data.Sigma_imp_list[0]
    Sigma_imp_hf  = it_data.Sigma_hartree_list[0]
    Sigma_dc      = it_data.Sigma_dc_list[0]

    ckpt.append(it_data, Eint_m_dc=Eint_m_dc)
    mpi.report(f"[csc {it + 1}] mu = {it_data.mu:.6f}  Eint-Edc = {Eint_m_dc:.6f} eV")

    # Wien2k charge update. Skipped on the last iteration so the calculation
    # ends on a DMFT step. lapw2 -qdmft folds Eint-Edc into the total energy
    # itself, so nothing further has to be added on this side.
    if it < n_csc_loops - 1:
        target_density, obe = driver.update_one_body_elements_with_charge_correction(N_k, Eint_m_dc, mu=it_data.mu, beta=beta)

ckpt.summarize()
