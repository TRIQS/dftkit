"""
VASP driver for TRIQS+DFT workflow automation.

This module provides a driver class for automating VASP calculations using
PLO projectors in the context of DMFT calculations with TRIQS/modest.

VASP runs as a persistent forked process communicating via a lock file mechanism,
unlike QE which is called as separate subprocesses for each step.
"""
import os
import shlex
import signal
import time
from dataclasses import dataclass
from datetime import datetime

import numpy as np
from h5 import HDFArchive
import triqs.utility.mpi as mpi

from .converter import Converter
from .plovasp.converter import generate_and_output_as_text


class DFTWorkflowError(Exception):
    """Exception raised for errors in the DFT workflow."""
    pass


# Default environment variables to preserve for subprocess execution
_DEFAULT_ENV_VARS = ['PATH', 'LD_LIBRARY_PATH', 'SHELL', 'PWD', 'HOME', 'OMP_NUM_THREADS', 'OMPI_MCA_btl_vader_single_copy_mechanism']

@dataclass
class MPIHandler:
    """
    Handles MPI configuration and execution environment for parallel DFT calculations.

    Attributes
    ----------
    mpi_exec : str
        The MPI executor command (e.g., "mpirun -np 16").
    """

    mpi_exec : str = "mpirun"

    def __post_init__(self):
        if not self.mpi_exec.strip():
            raise ValueError("mpi_exec cannot be empty")

    def __repr__(self):
        return f"MPIHandler(mpi_exec={self.mpi_exec})"

    __str__ = __repr__

    def get_env_vars(self):
        """Retrieve essential environment variables for subprocess execution."""
        env_vars = {}
        for var_name in _DEFAULT_ENV_VARS:
            var = os.getenv(var_name)
            if var: env_vars[var_name] = var
        return env_vars

    def report(self, *args):
        return mpi.report(*args)

    def is_master_node(self) -> bool:
        return mpi.is_master_node()


class Driver:
    """
    Driver for orchestrating VASP DFT calculations with PLO projectors.

    VASP is started once as a persistent forked process and communicates
    with the Python driver via a lock file mechanism:
    - VASP creates vasp.lock when starting, deletes it when SCF is done
    - Python creates vasp.lock to signal VASP to resume with updated charge
    - VASP reads the charge correction, does another SCF step, deletes lock

    Attributes
    ----------
    seedname : str
        Base name for HDF5 files.
    plo_cfg : str
        Path to the PLO configuration file.
    vasp_command : str
        VASP executable name (e.g., "vasp_std").
    mpi_handler : MPIHandler
        MPI configuration handler.
    use_ibz : bool or None
        Passed to Converter.convert_dft_input(): run on the irreducible BZ
        with symmetrization. None (default) enables it whenever VASP used
        symmetry and vaspout.h5 holds the symmetry data.
    """

    def __init__(self, seedname: str, plo_cfg: str, mpi_handler: MPIHandler,
                 vasp_command: str = "vasp_std", use_ibz: bool | None = None) -> None:
        self.seedname = seedname
        self.plo_cfg = plo_cfg
        self.vasp_command = vasp_command
        self.mpi_handler = mpi_handler
        self.use_ibz = use_ibz
        self._vasp_process_id = None

    def __repr__(self):
        return (f"VaspDriver(seedname={self.seedname}, plo_cfg={self.plo_cfg}, "
                f"vasp_command={self.vasp_command}, use_ibz={self.use_ibz})")

    __str__ = __repr__

    def _fork_and_start_vasp(self):
        """
        Fork a child process that runs VASP via MPI.

        The child process running VASP never returns from this function.
        The parent process returns the child's PID and continues.
        Only called on the master MPI node.

        Returns
        -------
        int
            Process ID of the VASP child process.
        """
        env_vars = self.mpi_handler.get_env_vars()

        mpi_parts = shlex.split(self.mpi_handler.mpi_exec)
        mpi_exe = mpi_parts[0]

        # Resolve full path to MPI executable
        for path_dir in env_vars.get('PATH', '').split(':'):
            candidate = os.path.join(path_dir, mpi_exe)
            if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                mpi_exe = candidate
                break

        arguments = [mpi_exe] + mpi_parts[1:] + [self.vasp_command]

        self.mpi_handler.report(f"[{datetime.now()}] Starting VASP: {' '.join(arguments)}")

        vasp_process_id = os.fork()
        if vasp_process_id == 0:
            # Child process: close inherited file descriptors to avoid h5 handle leaks
            for fd in range(3, 256):
                try:
                    os.close(fd)
                except OSError:
                    pass
            print('\n Starting VASP now\n', flush=True)
            os.execve(mpi_exe, arguments, env_vars)
            print('\n VASP exec failed\n', flush=True)
            os._exit(127)

        return vasp_process_id

    def _is_lock_file_present(self):
        """Check if vasp.lock exists (VASP is still running its SCF step)."""
        res_bool = False
        if mpi.is_master_node():
            res_bool = os.path.isfile('./vasp.lock')
        res_bool = mpi.bcast(res_bool)
        return res_bool

    def _is_vasp_running(self):
        """Check if the VASP child process is still alive."""
        if self._vasp_process_id is None:
            return False
        try:
            os.kill(self._vasp_process_id, 0)
            return True
        except OSError:
            return False

    def _wait_for_lock_to_appear(self, timeout=600):
        """Wait for VASP to create the lock file (indicating it has started)."""
        start = time.time()
        while not self._is_lock_file_present():
            if time.time() - start > timeout:
                raise DFTWorkflowError(f"VASP did not create vasp.lock within {timeout}s")
            time.sleep(1)
        mpi.barrier(poll_msec=100)

    def _wait_for_vasp(self, timeout=3600):
        """Wait for VASP to finish (lock file disappears)."""
        start = time.time()
        while self._is_lock_file_present():
            if time.time() - start > timeout:
                raise DFTWorkflowError(f"VASP did not finish within {timeout}s")
            if mpi.is_master_node() and not self._is_vasp_running():
                raise DFTWorkflowError("VASP process has died unexpectedly")
            time.sleep(1)
        mpi.barrier(poll_msec=100)

    def _run_plo_converter(self):
        """Run the PLO converter to generate H(k) and projectors in the HDF5 file."""
        if mpi.is_master_node():
            if not os.path.exists(self.plo_cfg):
                raise FileNotFoundError(f"PLO config file not found: {self.plo_cfg}")

            self.mpi_handler.report(f"[{datetime.now()}] Running PLO converter")
            generate_and_output_as_text(self.plo_cfg, vasp_dir='./')

            self.mpi_handler.report(f"[{datetime.now()}] Running VASP HDF5 converter")
            Converter(filename=self.seedname).convert_dft_input(use_ibz=self.use_ibz)

        mpi.barrier(poll_msec=100)

    def run_initial_stage(self, **kwargs) -> int:
        """
        Start VASP, run initial SCF, and convert output to HDF5.

        This forks VASP as a persistent background process, waits for the
        initial SCF to complete, then runs the PLO converter.
        """
        # Remove stale control files left over from a previous (crashed) run.
        # Only done here, on the initial call: a leftover vasp.lock would make
        # the driver think VASP is already running, and a leftover STOPCAR would
        # abort VASP immediately.
        if mpi.is_master_node():
            for stale in ('STOPCAR', 'vasp.lock'):
                if os.path.isfile(stale):
                    os.remove(stale)
        mpi.barrier(poll_msec=100)

        # Fork VASP on master node
        vasp_process_id = 0
        if mpi.is_master_node():
            vasp_process_id = self._fork_and_start_vasp()
        vasp_process_id = mpi.bcast(vasp_process_id)
        self._vasp_process_id = vasp_process_id

        # Wait for VASP to start (lock appears) then finish SCF (lock disappears)
        self._wait_for_lock_to_appear()
        self._wait_for_vasp()

        # Convert VASP output to HDF5
        self._run_plo_converter()

        return 0

    def run_update_stage(self, N_k, Eint_m_dc, **kwargs):
        """
        Run a single VASP charge-update SCF step and reconvert.

        Steps:
        1. Read DFT energy and compute band energy correction
        2. Write the charge density correction to vaspgamma.h5
        3. Create vasp.lock to trigger the VASP charge update
        4. Wait for VASP to finish
        5. Re-run the PLO converter with the updated projectors

        Running several VASP steps per DMFT iteration (with the self-energy held
        fixed) is left to the caller: invoke this once per step and recompute the
        charge density correction from the updated one-body elements in between.
        """
        dft_energy = self.read_dft_energy()

        band_energy_correction = self.band_energy_and_write_charge_update(N_k)

        mpi.report(f"DFT + DMFT Total Energy: {dft_energy + band_energy_correction + Eint_m_dc}")

        # Trigger VASP to resume with updated charge
        if mpi.is_master_node():
            open('./vasp.lock', 'a').close()
        mpi.barrier(poll_msec=100)

        # Wait for VASP to finish charge update
        self._wait_for_vasp()

        # Re-convert with updated projectors
        self._run_plo_converter()

        return 0

    def band_energy_and_write_charge_update(self, N_k):
        """
        Compute band energy correction and write charge density update to HDF5.

        Parameters
        ----------
        N_k : np.ndarray
            DMFT density matrix in band basis, shape (n_k, n_spin, n_max_bands, n_max_bands).

        Returns
        -------
        float
            Band energy correction in eV.
        """
        n_k, n_spin, n_max_bands, _ = N_k.shape

        fermi_weights, n_bands_per_k, Hk, bz_weights = None, None, None, None
        spin_to_data_index, SO = [0, 0], False
        if mpi.is_master_node():
            with HDFArchive(f"{self.seedname}.h5", "r") as ar:
                fermi_weights = ar['dft_misc_input']['dft_fermi_weights']
                n_bands_per_k = ar['dft_input']['n_orbitals']
                Hk            = ar['dft_input']['hopping']
                bz_weights    = ar['dft_input']['bz_weights']
                SO            = bool(ar['dft_input']['SO'])
                # spin-orbit: one spinor channel; collinear spin-pol: up/down; else: index 0 twice
                spin_to_data_index = [0] if SO else ([0, 1] if ar['dft_input']['SP'] else [0, 0])

        fermi_weights = mpi.bcast(fermi_weights)
        n_bands_per_k = mpi.bcast(n_bands_per_k)
        Hk            = mpi.bcast(Hk)
        bz_weights    = mpi.bcast(bz_weights)
        spin_to_data_index = mpi.bcast(spin_to_data_index)
        SO            = mpi.bcast(SO)
        mpi.barrier(poll_msec=100)

        density_matrix_dft = [[fermi_weights[ik, idx, :].astype(complex) for ik in range(n_k)] for idx in spin_to_data_index]
        band_energy_correction = 0.0
        for ik in range(n_k):
            for spin, idx in enumerate(spin_to_data_index):
                nb = n_bands_per_k[ik, idx]
                diag_indices = np.diag_indices(nb)
                N_k[ik, spin][diag_indices] -= density_matrix_dft[spin][ik][:nb]
                band_energy_correction += np.dot(N_k[ik][spin][:nb, :nb], Hk[ik, idx, :nb, :nb]).trace().real * bz_weights[ik]

        mpi.report(f"\nBand Energy Correction= {band_energy_correction} eV\n")

        # Build the per-k correction block: a single spinor channel for spin-orbit,
        # otherwise the spin-average of the two DMFT spin channels (non-magnetic / SP=0).
        delta_Nk = np.zeros((n_k, n_max_bands, n_max_bands), dtype=complex)
        for ik in range(n_k):
            nb0 = n_bands_per_k[ik, 0]
            if SO:
                delta_Nk[ik, :nb0, :nb0] = N_k[ik, 0, :nb0, :nb0]
            else:
                delta_Nk[ik, :nb0, :nb0] = 0.5 * (N_k[ik, 0, :nb0, :nb0] + N_k[ik, 1, :nb0, :nb0])

        if mpi.is_master_node():
            # Write the charge density correction in the format VASP's ADD_GAMMA_FROM_FILE
            # reads (vaspgamma.h5): a top-level band_window plus a deltaN group with one
            # (nb x nb) block per IBZ k-point, nb = number of bands in that k's window.
            # Mirrors triqs_dft_tools SumkDFT.calc_density_correction(dm_type='vasp'):
            #   - spin-orbit (SO):  a single spinor-band channel  -> deltaN/ud
            #   - otherwise (non-magnetic / spin-averaged): up == down -> deltaN/up, deltaN/down
            with HDFArchive(f"{self.seedname}.h5", "r") as ar:
                band_window = ar['dft_misc_input']['band_window']
                n_k_ibz     = ar['dft_misc_input']['n_k_ibz']
            deltaN_ibz = [delta_Nk[ik, :n_bands_per_k[ik, 0], :n_bands_per_k[ik, 0]]
                          for ik in range(n_k_ibz)]
            with HDFArchive('vaspgamma.h5', 'w') as vasp_h5:
                vasp_h5['band_window'] = [band_window[0][:n_k_ibz, :]]
                vasp_h5.create_group('deltaN')
                if SO:
                    vasp_h5['deltaN']['ud'] = deltaN_ibz
                else:
                    vasp_h5['deltaN']['up']   = deltaN_ibz
                    vasp_h5['deltaN']['down'] = deltaN_ibz

            # also keep the raw (full-grid) correction in the seedname h5 for the record
            subgrp = 'dft_update'
            with HDFArchive(f"{self.seedname}.h5", "a") as ar:
                if subgrp not in ar: ar.create_group(subgrp)
                ar[subgrp]["delta_N"] = delta_Nk

        return band_energy_correction

    def read_dft_energy(self):
        """
        Read DFT total energy from VASP output.

        Tries vaspout.h5 first, falls back to OSZICAR.
        VASP energies are already in eV.

        Returns
        -------
        float
            DFT total energy in eV.
        """
        dft_energy = 0.0
        if mpi.is_master_node():
            if os.path.isfile('vaspout.h5'):
                with HDFArchive('vaspout.h5', 'r') as h5:
                    if 'oszicar' in h5['intermediate/ion_dynamics']:
                        dft_energy = float(np.asarray(h5['intermediate/ion_dynamics/oszicar'][-1, 1]).item())
                        dft_energy = mpi.bcast(dft_energy)
                        return dft_energy

            # Fallback: read from OSZICAR
            last_nonempty_line = None
            with open('OSZICAR', 'r') as f:
                for line in f:
                    if line.strip():
                        last_nonempty_line = line

            if last_nonempty_line is None:
                raise DFTWorkflowError('OSZICAR is empty (cannot read DFT energy)')

            parts = last_nonempty_line.split()
            try:
                dft_energy = float(parts[2])
            except (IndexError, ValueError) as err:
                raise DFTWorkflowError(f'Failed to parse DFT energy from OSZICAR: {last_nonempty_line!r}') from err

        dft_energy = mpi.bcast(dft_energy)
        return dft_energy

    def kill(self):
        """Terminate the VASP process cleanly."""
        if self._vasp_process_id is not None and mpi.is_master_node():
            self.mpi_handler.report(f"[{datetime.now()}] Stopping VASP (PID {self._vasp_process_id})")
            with open('STOPCAR', 'wt') as f:
                f.write('LABORT = .TRUE.\n')
            try:
                os.kill(self._vasp_process_id, signal.SIGTERM)
            except OSError:
                pass
            self._vasp_process_id = None
