import numpy as np

import triqs.utility.mpi as mpi

from triqs.gfs import BlockGf, MeshImFreq

from triqs_ctseg.solve_generic import solve_generic

import triqs_modest as M

from triqs_dftkit.vasp import Driver as VaspDriver, MPIHandler
from triqs_modest.dft_driver import DftDriver

from h5 import HDFArchive


seedname = "vasp"
beta = 10.0
U    = 4.50
J    = 0.65
Up   = U -2*J
n_iw = 1000
n_total_loops = 10
n_dmft_loops  =  1

driver = DftDriver(VaspDriver(seedname=seedname, plo_cfg="plo.cfg",
                               mpi_handler=MPIHandler(mpi_exec="mpirun -np 8"),
                               vasp_command="vasp_std"))

target_density, obe = driver.one_body_elements_from_dft()
mpi.report(obe)

E = M.make_embedding(obe.C_space); mpi.report(E.description(True))

h_int = M.make_kanamori(E.sigma_names, E.imp_decomposition(0), U, Up, J, False, False)

DcTerm = M.DcSolver("NonPolarized", "cHeld", U, J)

mesh = MeshImFreq(beta=beta, S = "Fermion", n_iw=n_iw)

mu_dft = M.find_chemical_potential(target_density, obe, beta, verbosity=False)
Gdft = E.extract(M.local_gf.gloc(mesh, obe, mu_dft))[0]
mpi.report(f"Gdft density= {Gdft.total_density().real}")

deg_blocks = M.analyze_degenerate_blocks(Gdft)
mpi.report(f"deg_blocks= {deg_blocks}")

Sigma_imp_dc = DcTerm.dc_self_energy(Gdft)

Sigma_imp_dynamic, Sigma_imp_static = E.make_zero_imp_self_energies(mesh)[0]
for i in range(len(Sigma_imp_static)): Sigma_imp_static[i] += Sigma_imp_dc[i]

# DFT + DMFT loop
try:
    for n_iter in range(n_total_loops):
        mpi.report(f"\n{'='*60}\n=== Global DFT+DMFT iteration {n_iter+1}/{n_total_loops} ===\n{'='*60}")

        Hloc0 = E.extract(M.atomic_levels_and_delta.impurity_levels(obe))[0]
        mpi.report(f"Hloc0= {[h[0,0].real for h in Hloc0]}")

        # for first iteration converge Sigma imp first
        if n_iter == 0:
            n_dmft_loops_loc = 5
        else:
            n_dmft_loops_loc = n_dmft_loops
        # Begin DMFT loop
        for n_dmft_iter in range(n_dmft_loops_loc):
            mpi.report(f"\n--- DMFT iteration {n_dmft_iter+1}/{n_dmft_loops_loc} "
                       f"(global iter {n_iter+1}/{n_total_loops}) ---")

            Sigma_imp_static_minus_dc = [block-dc for (block, dc) in zip(Sigma_imp_static, Sigma_imp_dc)]

            Sigma_C_dynamic, Sigma_C_static = E.embed([Sigma_imp_dynamic], [Sigma_imp_static_minus_dc])

            mu = M.find_chemical_potential(target_density, obe, Sigma_C_dynamic, Sigma_C_static)

            Gloc = E.extract(M.local_gf.gloc(obe, mu, Sigma_C_dynamic, Sigma_C_static))[0]

            Eimp = [h-mu*np.eye(h.shape[0])-dc for (h,dc) in zip(Hloc0, Sigma_imp_dc)]

            Delta = M.symmetrize(M.hybridization(Eimp, Gloc, Sigma_imp_dynamic, Sigma_imp_static), deg_blocks)

            solver_params = dict(n_iw=n_iw, n_tau=10*n_iw, length_cycle=50,
                                 n_cycles = int(1e+6/mpi.size),
                                 n_warmup_cycles = int(1e+4),
                                 perform_tail_fit=True,
                                 fit_min_w=10, fit_max_w=14,
                                 )
            # CT-SEG needs a real h_loc0: op_from_block_matrix flags coefficients complex
            # whenever the input matrices are complex dtype (work_data.cpp extracts mu via a
            # real_or_complex->double cast that throws on complex-flagged values). The impurity
            # levels are Hermitian/real, so drop the ~1e-18 numerical imaginary noise here.
            Eimp_real = [h.real.copy() for h in Eimp]
            solver_results = solve_generic(Delta, Eimp_real, h_int, **solver_params)

            Sigma_imp_dynamic, Sigma_imp_static = solver_results.Sigma_dynamic, solver_results.Sigma_HartreeFock

            solver_results.G_iw         <<  M.symmetrize(solver_results.G_iw,          deg_blocks)
            solver_results.Sigma_iw     <<  M.symmetrize(solver_results.Sigma_iw,      deg_blocks)
            solver_results.Sigma_dynamic << M.symmetrize(solver_results.Sigma_dynamic, deg_blocks)
        # End of DMFT loop
        # Update Double-counting Term
        Sigma_imp_dc = DcTerm.dc_self_energy(solver_results.G_iw)
        Eimp_dc      = DcTerm.dc_energy(solver_results.G_iw)

        Sigma_imp_static_minus_dc = [block-dc for (block, dc) in zip(Sigma_imp_static, Sigma_imp_dc)]

        # Re-embed the self-energy
        Sigma_C_dynamic, Sigma_C_static = E.embed([Sigma_imp_dynamic], [Sigma_imp_static_minus_dc])

        # Re-compute the chemical potential
        mu = M.find_chemical_potential(target_density, obe, Sigma_C_dynamic, Sigma_C_static)

        # Compute the charge density correction
        N_k = M.charge_density_correction(obe, mu, Sigma_C_dynamic, Sigma_C_static)

        # Impurity interaction energy
        Eint = 0.5 * np.real((solver_results.Sigma_iw*solver_results.G_iw).total_density())
        mpi.report(f"Eint= {Eint}")
        Eint_m_dc = (Eint - Eimp_dc)
        mpi.report(f"Eint-Edc= {Eint_m_dc}")

        mpi.report("Saving DFT + DMFT iteration...")
        if mpi.is_master_node():
            with HDFArchive("checkpoint.h5", "a") as ar:
                path = f"it={len(ar) + 1}"
                ar.create_group(path)
                ar[path]["mu"]              = mu
                ar[path]["nloc"]            = Gloc.total_density().real
                ar[path]["nimp"]            = solver_results.G_iw.total_density().real
                ar[path]["Delta_iw"]        = Delta
                ar[path]["Eimp"]            = Eimp
                ar[path]["Gloc_iw"]         = Gloc
                ar[path]["Sigma_iw_static"] = solver_results.Sigma_HartreeFock
                ar[path]["Sigma_dc"]        = Sigma_imp_dc

        # Update the one-body Hamiltonian with the charge density correction
        mpi.report(f"Calling VASP charge update / DFT driver "
                   f"(global iter {n_iter+1}/{n_total_loops})...")
        obe = driver.update_one_body_elements_with_charge_correction(N_k, Eint_m_dc)[1]

finally:
    # Ensure VASP is killed even if the script crashes
    driver.driver.kill()
