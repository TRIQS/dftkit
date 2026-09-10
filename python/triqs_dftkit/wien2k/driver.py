"""
WIEN2k driver for TRIQS+DFT workflow automation.

This module provides a driver class for automating WIEN2k calculations in the
context of charge self-consistent (CSC) DMFT calculations with TRIQS/modest.

WIEN2k is a *chain* of Fortran programs -- lapw0, lapw1, lapw2, lcore,
mixer -- each launched through the ``x`` tcsh script, which writes the ``.def``
file of unit -> filename assignments that the program reads.  The SCF loop that
normally drives them lives in the ``run_lapw`` tcsh script.

This driver replaces ``run_lapw``: it owns the loop in Python so that modest's
DftDriver can inject a DMFT charge density correction between iterations.
WIEN2k's own QDMFT support cannot be reused because it inverts the control flow
-- ``run_lapw`` calls out to a user python script from inside its cycle, whereas
DftDriver requires python to be the caller.

Scope: serial, non-magnetic (SP=0, SO=0).  ``lapw2 -qdmft`` cannot be run in
parallel at all.
"""
import os, re, shutil, subprocess
from datetime import datetime

import numpy as np
import triqs.utility.mpi as mpi

from .converter import Converter
from ..converter_tools import ConverterTools


class DFTWorkflowError(Exception):
    """Exception raised for errors in the DFT workflow."""
    pass


# Default environment variables to preserve for subprocess execution.  WIENROOT
# and SCRATCH are WIEN2k specific: x interpolates $WIENROOT into every .def file
# it writes (e.g. the xc_funcs.h path for lapw0) and tcsh aborts outright on an
# undefined variable, so the environment cannot simply be stripped to the
# defaults the other drivers use.
_DEFAULT_ENV_VARS = ['PATH', 'LD_LIBRARY_PATH', 'SHELL', 'PWD', 'HOME', 'OMP_NUM_THREADS',
                     'OMPI_MCA_btl_vader_single_copy_mechanism', 'WIENROOT', 'SCRATCH']

# WIEN2k works in Rydberg, TRIQS in eV.
_RY_IN_EV = 13.605698

# Flags that make x rewrite the first five characters of case.in2 in place and
# restore the original from .oldin2 afterwards.
_IN2_MODE_FLAGS = ('-almd', '-qdmft', '-fermi', '-qtl', '-alm', '-efg')

# Programs that run_lapw calls with -c on a complex (no inversion symmetry) case,
# which makes x launch lapw1c/lapw2c and read case.in1c/case.in2c.  lapw0, lcore
# and mixer have no complex variant and reject the flag.
_CMPLX_PROGRAMS = ('lapw1', 'lapw2')

# Files saved to <name>_old before lapw0 and before mixer as run_lapw.
# case.clmsum_old is mixer's previous-iteration density
# on unit 10, so the second set is required for mixing to work at all.
_LAPW0_SAVE = ('vsp', 'vns', 'r2v')
_MIXER_SAVE = ('clmsum', 'vrespsum', 'tausum')

# Partial scf files concatenated into case.scf each cycle, in run_lapw's order
# for the non-spin-polarised non-HF branch.  Missing ones are skipped.
_SCF_PARTS = ('0', '1', 'so', '2', '1s', '2s', 'c')

# ':DIS  :  CHARGE DISTANCE      (<f11.7> for atom<i5> spin<i2>)<f15.7>'.  The
# value inside the parentheses is the per-atom maximum, which is what testconv
# tests against the -cc limit; the trailing number is the cell total.
_DIS_RE = re.compile(r'\(\s*([-+0-9.EDed]+)\s+for atom')

# DMFT label marker in the : LABEL style of run_lapw                     
_DMFT_MARKER = ':QDMFT:  CHARGE DENSITY UPDATED FROM DMFT OCCUPATIONS'


class Driver(object):
    """
    Driver orchestrating the WIEN2k program chain for CSC DFT+DMFT.

    Satisfies modest's DftDriver contract: it exposes ``seedname``,
    ``run_initial_stage`` and ``run_update_stage``, and leaves
    ``<seedname>.h5`` fully re-converted before either hook returns.

    Attributes
    ----------
    seedname : str
        WIEN2k case name.  All WIEN2k files are ``<seedname>.<ext>`` and the
        archive modest reads is ``<seedname>.h5``.
    wienroot : str
        WIEN2k installation directory; ``<wienroot>/x`` launches every program.
    dmftproj_exe : str
        The dmftproj executable.  Looked up on $PATH by default, matching how
        run_lapw invokes it.
    max_scf_iter : int
        Iteration cap for the initial SCF loop.
    ecut, ccut : float
        Energy (Ry) and charge convergence limits, as run_lapw's -ec and -cc.
    verbosity : int
        0 -- warnings only; 1 -- also per-cycle convergence and restart messages
        (default); 2 -- also one line per WIEN2k program launched.  Output is
        written on the master rank only.  Settable at any time as an attribute.

    Notes
    -----
    ``wienroot`` and ``dmftproj_exe`` are the only coupling to WIEN2k, so a fake
    WIEN2k can be substituted for testing by pointing them elsewhere.
    """

    def __init__(self, seedname, wienroot=None, dmftproj_exe='dmftproj',
                 max_scf_iter=100, ecut=1e-4, ccut=1e-4, verbosity=1):
        self.seedname = seedname
        self.wienroot = wienroot if wienroot is not None else os.getenv('WIENROOT')
        if not self.wienroot:
            raise DFTWorkflowError(
                "WIENROOT is not set and no wienroot= was given; cannot locate the WIEN2k 'x' script.")
        self.x_exe = os.path.join(self.wienroot, 'x')
        if not os.path.isfile(self.x_exe):
            raise DFTWorkflowError(f"No 'x' script at {self.x_exe} (is wienroot= correct?)")
        self.dmftproj_exe = dmftproj_exe
        self.max_scf_iter = max_scf_iter
        self.ecut = ecut
        self.ccut = ccut
        self.verbosity = verbosity
        self.fortran_to_replace = {'D': 'E'}
        # Set when this process inherits an SCF it did not run, and cleared by the
        # next mixer call.  See _flush_mixing_history.
        self._mixing_flush_pending = False

    def __repr__(self):
        return (f"Wien2kDriver(seedname={self.seedname}, wienroot={self.wienroot}, "
                f"dmftproj_exe={self.dmftproj_exe})")

    __str__ = __repr__

    # ------------------------------------------------------------------ files

    def _report(self, message, level=1):
        """
        Print progress on the master rank, if verbosity allows.

        Deliberately a per-driver setting rather than TRIQS's module-level
        ``mpi.Verbosity_Level_Report_Max``, since turning that down to quieten a
        few hundred lines of DFT bookkeeping would also quieten the solver.
        """
        if level <= self.verbosity:
            mpi.report(message)

    @staticmethod
    def _warn(message):
        """
        Report a problem that does not stop the run.

        Not gated on verbosity, and sent to stderr: these mean the calculation is
        probably wrong rather than merely noisy, so silencing progress output
        should not silence them too.
        """
        mpi.report(f"WARNING: {message}", stderr=True)

    def _f(self, ext):
        """Path of the WIEN2k file ``<seedname>.<ext>``."""
        return f"{self.seedname}.{ext}"

    @property
    def _cmplx(self):
        """'c' when this is a complex (no inversion symmetry) case, else ''."""
        in1c = self._f('in1c')
        return 'c' if os.path.isfile(in1c) and os.path.getsize(in1c) > 0 else ''

    def _in2_file(self):
        return self._f('in2' + self._cmplx)

    def _set_in2_mode(self, mode):
        """
        Overwrite the first five characters of case.in2 line 1 with ``mode``.

        x does this itself for -almd/-qdmft, but a crashed run can leave the
        wrong keyword behind, which would silently turn a plain density pass
        into a projector pass.
        """
        path = self._in2_file()
        with open(path) as fh:
            lines = fh.readlines()
        if not lines:
            raise DFTWorkflowError(f"{path} is empty")
        lines[0] = f"{mode:<5}" + lines[0].rstrip('\n')[5:] + '\n'
        with open(path, 'w') as fh:
            fh.writelines(lines)

    # -------------------------------------------------------------- execution

    def _env(self):
        """Environment for WIEN2k subprocesses, with WIENROOT forced to ours."""
        env = {}
        for name in _DEFAULT_ENV_VARS:
            value = os.getenv(name)
            if value:
                env[name] = value
        env['WIENROOT'] = self.wienroot
        return env

    def _run_x(self, program, *flags):
        """
        Run ``x <program> <flags>`` and verify it succeeded.  Master node only.

        Success is checked twice, because WIEN2k signals failure two ways: x
        exits 9 when the program returns non-zero, and every program writes a
        message into ``<program>.error`` on entry and truncates it again just
        before a successful exit.  The error file therefore reports failures
        that never reach the exit status, including untrapped runtime errors.
        """
        if not mpi.is_master_node():
            return 0

        # x does not detect a complex case on its own; run_lapw passes -c.
        if self._cmplx and program in _CMPLX_PROGRAMS:
            flags = ('-c',) + flags

        # x backs case.in2 up to .oldin2 only when that file is absent, and
        # -qdmft seds .oldin2 rather than case.in2 -- so a .oldin2 left behind by
        # an earlier crash is silently used as the source and then restored over
        # the real input.  x warns about this but does not clean it up.
        if any(flag in _IN2_MODE_FLAGS for flag in flags):
            for stale in ('.oldin2', '.oldin2a'):
                if os.path.isfile(stale):
                    os.remove(stale)

        return self._run_checked([self.x_exe, program, *flags], f'{program}.error',
                                 f"x {program} {' '.join(flags)}".rstrip())

    def _run_checked(self, command, error_file, label):
        """
        Run a WIEN2k command and verify it succeeded.

        Success is checked twice, because WIEN2k signals failure two ways: a
        non-zero exit status, and a message left in the .error file.  Every program
        writes that message on entry and truncates it again just before a
        successful exit, so it catches failures that never reach the exit status,
        including untrapped runtime errors.
        """
        self._report(f"[{datetime.now()}] running {' '.join(command)}", level=2)
        result = subprocess.run(command, capture_output=True, text=True, env=self._env())
        tail = f"{result.stdout[-2000:]}{result.stderr[-2000:]}"
        if result.returncode != 0:
            raise DFTWorkflowError(f"{label} failed with code {result.returncode}\n"
                                   f"{self._error_text(error_file)}{tail}")
        message = self._error_text(error_file)
        if message:
            raise DFTWorkflowError(f"{label} exited cleanly but left {error_file}:\n{message}")
        return 0

    @staticmethod
    def _error_text(error_file):
        """Contents of a WIEN2k .error file; empty string means success."""
        if not os.path.isfile(error_file) or os.path.getsize(error_file) == 0:
            return ''
        with open(error_file) as fh:
            return fh.read().strip()

    def _run_dmftproj(self):
        """Run dmftproj to build the projectors the converter reads."""
        if not mpi.is_master_node():
            return 0
        if not os.path.isfile(self._f('indmftpr')):
            raise DFTWorkflowError(
                f"{self._f('indmftpr')} not found; copy and edit "
                f"{os.path.join(self.wienroot, 'SRC_templates', 'case.indmftpr')}"
                "or use the command line tool init_dmftpr")
        return self._run_checked([self.dmftproj_exe], 'dmftproj.error', 'dmftproj')

    # ------------------------------------------------------------- scf output

    def _mark_dmft_cycle(self):
        """Record in case.scf that the cycle starting here is a DMFT charge update."""
        if not mpi.is_master_node():
            return
        with open(self._f('scf'), 'a') as out:
            out.write(f"\n{_DMFT_MARKER}\n\n")

    def _remove_broyden_files(self):
        """
        Delete case.broyd*, which is how run_lapw discards a mixing history.

        run_lapw removes them whenever it starts an SCF rather than continuing
        one, unless -NI is given.  Deleting the files is unambiguous: mixer
        rebuilds the history from scratch because there is nothing to read.
        """
        prefix = self._f('broyd')
        for name in os.listdir('.'):
            if name.startswith(prefix):
                os.remove(name)

    def _run_mixer(self):
        """
        Run mixer, first flushing the mixing history if this process inherited one.

        The flush is deferred to here rather than done where it is decided, so
        that case.broyd* is only ever removed immediately before the mixer call
        that would otherwise have read it, leaving the history in place for every
        other mixer call in the run.
        """
        if self._mixing_flush_pending:
            if mpi.is_master_node():
                self._remove_broyden_files()
            self._mixing_flush_pending = False
            mpi.barrier(poll_msec=100)
        self._run_x('mixer')

    def _flush_mixing_history(self, reason):
        """
        Arrange for the next mixer call to restart its mixing, and say so.

        Called when this process is about to mix against a case.broyd* history it
        did not build.  The Broyden vectors encode a sequence of density updates,
        so continuing against a history whose density had a different character --
        a DMFT-corrected one, or one from an interrupted run -- mixes the current
        residual into predictions derived from unrelated ones.
        """
        self._warn(f"{reason}; removing case.broyd* as run_lapw does, so the next "
                   "mixer call rebuilds its mixing history from scratch")
        self._mixing_flush_pending = True

    def _append_to_scf(self, parts):
        """
        Append partial scf files to case.scf, as run_lapw does.

        Called twice per cycle: once for the programs before mixer, once for
        case.scfm after it.
        """
        if not mpi.is_master_node():
            return
        with open(self._f('scf'), 'a') as out:
            for part in parts:
                path = self._f('scf' + part)
                if os.path.isfile(path):
                    with open(path) as fh:
                        out.write(fh.read())

    def _save_old(self, exts):
        """
        Copy case.<ext> to case.<ext>_old for each ext that exists.

        run_lapw does this at two points in every cycle and both are load
        bearing.  Most importantly mixer reads case.clmsum_old on unit 10 as the
        previous-iteration density: without the copy it mixes against a missing 
        or stale density, which emits no :DIS line and sends the SCF diverging until 
        the linearisation energies go bad.
        """
        if not mpi.is_master_node():
            return
        for ext in exts:
            path = self._f(ext)
            if os.path.isfile(path):
                shutil.copyfile(path, self._f(ext + '_old'))

    @staticmethod
    def _scf_tags(path, stop_at_dmft=False):
        """
        Parse a WIEN2k scf file into one ``(ene, dis)`` pair per mixer record.

        :ENE is written in three variants (**INFO****, *WARNING**, **********), so
        the energy is the last whitespace token rather than a fixed column.  :DIS
        carries two numbers; the one inside the parentheses is the per-atom
        maximum, which is what testconv compares against the -cc limit, and the
        trailing one is the cell total.

        Pairing is done by record rather than by zipping two independent lists.
        mixer writes :DIS only inside, i.e. only when it could read case.clmsum_old, 
        but writes :ENE unconditionally, so a cycle with no previous density contributes 
        an :ENE with no :DIS.  Zipping and padding the shorter list at the end would attribute
        every later :DIS to the cycle before its own and leave the newest cycle
        with ``dis=None``, which silently disables the charge convergence test.
        Within one record :DIS precedes :ENE, so :ENE closes the record.

        ``stop_at_dmft`` stops at the first _DMFT_MARKER, leaving the caller only
        the plain-DFT part of the history.
        """
        history, dis = [], None
        if not os.path.isfile(path):
            return history
        with open(path) as fh:
            for line in fh:
                if stop_at_dmft and line.startswith(_DMFT_MARKER):
                    break
                if line.startswith(':ENE'):
                    history.append((float(line.split()[-1]), dis))
                    dis = None
                elif line.startswith(':DIS'):
                    match = _DIS_RE.search(line)
                    dis = (float(match.group(1).replace('D', 'E').replace('d', 'e'))
                           if match else float(line.split()[-1]))
        return history

    def _read_scfm(self):
        """
        Return ``(ene, dis)`` from case.scfm: the total energy in Ry and the
        per-atom charge distance.

        :ENE is written in three variants (**INFO****, *WARNING**, **********),
        so the value is taken as the last whitespace token rather than by column.
        """
        path = self._f('scfm')
        if not os.path.isfile(path):
            raise DFTWorkflowError(f"{path} was not written; mixer did not run")
        history = self._scf_tags(path)
        if not history:
            raise DFTWorkflowError(f"no :ENE line in {path}")
        ene, dis = history[-1]
        if dis is None:
            # mixer omits :DIS when it has no previous density to compare
            # against, which means case.clmsum_old was missing -- the mixing is
            # then meaningless even though mixer exits cleanly.
            self._warn(f"no :DIS line in {path}; mixer had no previous density, "
                       "so the charge convergence test is being skipped")
        return ene, dis

    def read_dft_energy(self):
        """
        Total energy in eV, from the last :ENE line in case.scf.

        After a charge update this **already includes** the interaction energy:
        lapw2 -qdmft folds ``correner`` into :SUM (ETOT = ETOT + correner) and
        subtracts the DFT in-window band energy itself.  Callers must therefore
        not add ``Eint_m_dc`` again, and there is no separate band energy
        correction to compute.
        """
        energy = None
        if mpi.is_master_node():
            path = self._f('scf')
            if not os.path.isfile(path):
                raise DFTWorkflowError(f"{path} does not exist; no SCF has run")
            # Deliberately not stop_at_dmft: this must report the current total
            # energy, which after a charge update is the DMFT one.
            history = self._scf_tags(path)
            if not history:
                raise DFTWorkflowError(f"no :ENE line in {path}")
            energy = history[-1][0] * _RY_IN_EV
        return mpi.bcast(energy)

    # -------------------------------------------------------------- scf cycle

    def _check_inputs(self, need_projectors=True):
        """
        Verify the inputs are in place before launching anything.

        ``case.indmftpr`` is only consumed at the very end, by dmftproj, so
        without an up-front check a missing one would not surface until a full
        SCF had already run.
        """
        if not mpi.is_master_node():
            return
        # x takes the case name from the directory it runs in, while this driver
        # takes it from seedname.  If the two disagree, x reads and writes a
        # different set of case.* files than the driver looks at.
        cwd_case = os.path.basename(os.getcwd())
        if self.seedname != cwd_case:
            self._warn(f"seedname is '{self.seedname}' but the working directory is "
                       f"'{cwd_case}'; the x script derives the case name from the "
                       "directory, so these must normally match")

        required = ['struct', 'in0', 'in1' + self._cmplx, 'in2' + self._cmplx, 'inm', 'inc']
        if need_projectors:
            required.append('indmftpr')
        missing = [self._f(ext) for ext in required if not os.path.isfile(self._f(ext))]
        if missing:
            hint = 'run init_lapw'
            if self._f('indmftpr') in missing:
                hint += ("; case.indmftpr is not made by init_lapw -- generate it with "
                         "init_dmftpr, or copy and edit "
                         f"{os.path.join(self.wienroot, 'SRC_templates', 'case.indmftpr')}")
            raise DFTWorkflowError(f"missing WIEN2k input file(s): {', '.join(missing)}; {hint}")
        if not os.path.isfile(self._f('clmsum')):
            if os.path.isfile(self._f('clmsum_old')):
                self._warn(f"{self._f('clmsum')} missing, recovering from clmsum_old")
                with open(self._f('clmsum_old'), 'rb') as src, open(self._f('clmsum'), 'wb') as dst:
                    dst.write(src.read())
            else:
                raise DFTWorkflowError(
                    f"no {self._f('clmsum')} or {self._f('clmsum_old')}, which lapw0 needs; run dstart")

    def _scf_iteration(self):
        """
        One SCF cycle: lapw0 -> lapw1 -> lapw2 -> lcore -> mixer.

        Returns ``(ene, dis)`` read from case.scfm and broadcast, so that every
        rank reaches the same convergence verdict and stays in lockstep.
        """
        if mpi.is_master_node():
            self._set_in2_mode('TOT')

        self._save_old(_LAPW0_SAVE)
        for program in ('lapw0', 'lapw1', 'lapw2', 'lcore'):
            self._run_x(program)
        self._append_to_scf(_SCF_PARTS)
        self._save_old(_MIXER_SAVE)
        self._run_mixer()
        self._append_to_scf(('m',))

        values = self._read_scfm() if mpi.is_master_node() else None
        return mpi.bcast(values)

    def _converged(self, history):
        """
        Apply testconv's criterion to the (ene, dis) history.

        The energy test is the mean of the last two |dE| over three iterations,
        not a single difference, and needs three points before it can fire.
        """
        if len(history) < 3:
            return False
        e3, e2, e1 = (entry[0] for entry in history[-3:])
        ene_ok = 0.5 * (abs(e1 - e3) + abs(e1 - e2)) < self.ecut
        dis = history[-1][1]
        dis_ok = dis is None or dis < self.ccut
        return ene_ok and dis_ok

    def _prepare_fresh_scf(self):
        """
        Housekeeping before the first cycle of a new SCF.

        Only ever for a genuine fresh start.  Dropping case.broyd* part way
        through would discard the mixing history and stall convergence, which
        matters when the cycle is being stepped one iteration at a time.
        """
        if not mpi.is_master_node():
            return
        self._set_in2_mode('TOT')
        self._remove_broyden_files()
        # Truncate case.scf too.  Without this the abandoned cycles stay on disk
        # and _scf_history_on_disk splices them onto the new ones, so a later
        # restart applies the convergence test across the seam between two
        # unrelated runs -- and the caller was promised the history was discarded.
        if os.path.isfile(self._f('scf')):
            self._report(f"discarding the {len(self._scf_history_on_disk())} cycles "
                         f"already in {self._f('scf')}")
            open(self._f('scf'), 'w').close()

    def _run_scf(self, n_iter=None, history=None):
        """
        Iterate the SCF cycle.

        ``n_iter=None`` runs to convergence, up to ``max_scf_iter``, and raises if
        it is never reached.  ``n_iter=N`` runs exactly N cycles and returns
        whatever state they reached, without requiring convergence.

        ``history`` seeds the convergence test; when omitted it is recovered from
        case.scf.  That matters for a stepped run: the energy criterion needs
        three points, which a history rebuilt from scratch on every call would
        never accumulate.
        """
        if history is None:
            history = self._scf_history_on_disk() if mpi.is_master_node() else None
            history = mpi.bcast(history)
        history = list(history)
        done = len(history)
        limit = self.max_scf_iter if n_iter is None else n_iter

        for step in range(1, limit + 1):
            history.append(self._scf_iteration())
            self._report(f"    cycle {done + step}: :ENE = {history[-1][0]:.8f} Ry  "
                       f":DIS = {history[-1][1]}")
            if n_iter is None and self._converged(history):
                self._report(f"SCF converged after {done + step} cycles")
                return history

        if n_iter is None:
            raise DFTWorkflowError(
                f"SCF did not converge in {self.max_scf_iter} cycles "
                f"(ecut={self.ecut} Ry, ccut={self.ccut})")
        return history

    def _ensure_converged_scf(self, force_scf=False):
        """
        Bring the SCF to convergence, reusing whatever case.scf already holds.

        Already converged -> nothing to run.  Partially converged -> continue from
        there, keeping the Broyden history, since restarting would throw away both
        the completed cycles and the mixing state.  Cold start, or force_scf ->
        begin afresh.
        """
        history = self._scf_history_on_disk() if mpi.is_master_node() else None
        history = mpi.bcast(history)

        if not force_scf and self._converged(history):
            self._report(f"case.scf already holds a converged SCF ({len(history)} cycles); "
                       "reusing it (force_scf=True to redo)")
            resumed = self._has_dmft_cycles() if mpi.is_master_node() else None
            if mpi.bcast(resumed):
                # A CSC run being resumed: case.clmsum is the DMFT-corrected
                # density from the interrupted run, and the mixing history behind
                # it belongs to that run.  Nothing here reruns the SCF -- the
                # density is kept as it stands -- but the next charge update will
                # mix against that inherited history.
                self._flush_mixing_history(
                    "case.scf records DMFT charge updates from an earlier run")
            return history

        if force_scf or not history:
            self._prepare_fresh_scf()
            mpi.barrier(poll_msec=100)
            return self._run_scf(history=[])

        self._report(f"continuing the SCF from {len(history)} cycles already in case.scf")
        self._flush_mixing_history(
            f"resuming an SCF from {len(history)} cycles this process did not run")
        return self._run_scf(history=history)

    def _has_dmft_cycles(self):
        """Whether case.scf records a DMFT charge update from an earlier run."""
        path = self._f('scf')
        if not os.path.isfile(path):
            return False
        with open(path) as fh:
            return any(line.startswith(_DMFT_MARKER) for line in fh)

    def scf_converged(self):
        """
        Whether the SCF history in case.scf already meets the ecut/ccut criteria.

        Exposed so the cycle can be driven a step at a time::

            while not driver.scf_converged():
                driver.run_dft_only(n_iter=1)
        """
        history = self._scf_history_on_disk() if mpi.is_master_node() else None
        return self._converged(mpi.bcast(history))

    def run_dft_only(self, n_iter=None, force_scf=False):
        """
        Run the DFT SCF cycle and stop: no projectors, no dmftproj, no HDF5.

        Not part of the DftDriver contract.  This is for plain WIEN2k runs and for
        driving the cycle under external control; case.indmftpr is not needed.

        Parameters
        ----------
        n_iter : int, optional
            Run exactly this many cycles and return, converged or not.  Pass 1 to
            advance a single step.  The default, None, iterates to convergence and
            raises if max_scf_iter is exhausted.
        force_scf : bool, optional
            Start over from scratch, discarding the accumulated case.scf history
            and the Broyden files, instead of reusing a converged result or
            continuing a partial one.  Ignored when n_iter is given, since a fixed
            number of cycles was asked for explicitly.

        Returns
        -------
        list of (float, float)
            The (:ENE in Ry, :DIS) history, including cycles already on disk.
        """
        self._check_inputs(need_projectors=False)

        if n_iter is None:
            return self._ensure_converged_scf(force_scf)

        # Stepping: prepare only on a genuine cold start, so that repeated
        # single-step calls keep their mixing history and their :ENE record.
        if not self._scf_history_on_disk():
            self._prepare_fresh_scf()
            mpi.barrier(poll_msec=100)
        return self._run_scf(n_iter=n_iter)

    def _scf_history_on_disk(self):
        """
        (ene, dis) pairs recovered from an existing case.scf, for restart.

        Stops at the first DMFT marker: the cycles after it were produced by
        run_update_stage against a DMFT-corrected density, and the DFT ecut/ccut
        limits say nothing useful about them.
        """
        return self._scf_tags(self._f('scf'), stop_at_dmft=True)

    # --------------------------------------------------------- dmft interface

    def _read_oubwin(self):
        """
        Read case.oubwin, written by dmftproj.

        Returns ``(iso, windows)`` where windows is a list of
        ``(included, nb_bot, nb_top, weight)``, one per k-point, with 1-based
        inclusive band indices.  This file -- not dft_input/n_orbitals -- is what
        lapw2 -qdmft cross-checks case.qdmft against, so it is the authority on
        the per-k window.
        """
        path = self._f('oubwin')
        if not os.path.isfile(path):
            raise DFTWorkflowError(f"{path} not found; dmftproj must run before the charge update")
        reader = ConverterTools.read_fortran_file(self, path, self.fortran_to_replace)
        try:
            n_k = int(next(reader))
            iso = int(next(reader))
            windows = []
            for _ in range(n_k):
                included = int(next(reader))
                if included == 1:
                    nb_bot, nb_top = int(next(reader)), int(next(reader))
                    weight = next(reader)
                    windows.append((included, nb_bot, nb_top, weight))
                else:
                    windows.append((included, None, None, None))
        except StopIteration:
            raise DFTWorkflowError(f"wien2k: reading file {path} failed!")
        return iso, windows

    def _write_qdmft(self, N_k, Eint_m_dc, mu=0.0, beta=0.0):
        """
        Write the DMFT band occupation matrices to case.qdmft for lapw2 -qdmft.

        The layout is fixed by the reader in SRC_lapw2/qdmft.F::readdata_qdmft::

            mu                       in Ry  (read, then unused)
            beta                     in Ry^-1 (read, then unused)
            per k-point:
                nn                   must equal nb_top - nb_bot + 1
                nn records of 2*nn reals: Re Im, one record per matrix *row*
                one throwaway record
            correner                 in eV; lapw2 divides it by 13.605698

        Three conventions differ from the VASP and QE writers:

        * the **full** occupation matrix is written, not the deviation from the
          Kohn-Sham density -- lapw2 excises the in-window DFT bands itself;
        * the matrix must be **unweighted**, because lapw2 multiplies it by the
          k-point weight it reads from case.oubwin;
        * ``correner`` is written in eV with no conversion.

        The blank line after each matrix is mandatory: the reader issues a
        ``READ(32,*)`` with an empty io-list, which consumes one whole record.
        """
        if not mpi.is_master_node():
            return

        _, windows = self._read_oubwin()
        if len(windows) != N_k.shape[0]:
            raise DFTWorkflowError(
                f"case.oubwin has {len(windows)} k-points but N_k has {N_k.shape[0]}")

        excluded = [ik for ik, w in enumerate(windows) if w[0] != 1]
        if excluded:
            # lapw2 loops over every k-point when reading case.qdmft but only
            # sets nn for included ones, so an excluded k-point is compared
            # against uninitialised memory.
            raise DFTWorkflowError(
                f"case.oubwin marks k-point(s) {excluded} as not included; lapw2 -qdmft "
                "requires every k-point to be inside the correlated window "
                "(widen the energy window in case.indmftpr)")

        n_sigma = N_k.shape[1]
        if n_sigma != 2:
            raise DFTWorkflowError(
                f"expected 2 spin channels for a non-magnetic case, got {n_sigma}; "
                "spin-polarised and spin-orbit cases are not supported yet")
        # Both channels are computed independently even when SP=0, so average
        # them.  Build a new array: modest hands over the caller's N_k uncopied.
        Nk_avg = 0.5 * (N_k[:, 0] + N_k[:, 1])

        for ik, (_, nb_bot, nb_top, _w) in enumerate(windows):
            nn = nb_top - nb_bot + 1
            if nn > Nk_avg.shape[1]:
                raise DFTWorkflowError(
                    f"k-point {ik}: case.oubwin window is {nn} bands but N_k only "
                    f"has {Nk_avg.shape[1]}")

        with open(self._f('qdmft'), 'w') as fh:
            # Converted to Ry / Ry^-1 for the file even though lapw2 discards
            # both, so the values match what the format documents and what
            # legacy dft_tools wrote.
            fh.write("%.14f\n" % (mu / _RY_IN_EV))
            fh.write("%.14f\n" % (beta * _RY_IN_EV))
            for ik, (_, nb_bot, nb_top, _w) in enumerate(windows):
                nn = nb_top - nb_bot + 1
                fh.write("%s\n" % nn)
                block = Nk_avg[ik, :nn, :nn]
                for row in range(nn):
                    fh.write(''.join(f"{block[row, col].real:.14f} {block[row, col].imag:.14f} "
                                     for col in range(nn)))
                    fh.write("\n")
                fh.write("\n")                       # the mandatory throwaway record
            fh.write("%.16f\n" % np.real(Eint_m_dc))

    def _regenerate_projectors(self):
        """
        Rebuild the projectors and reconvert the archive.
        """
        self._run_x('lapw2', '-almd')
        self._run_dmftproj()
        if mpi.is_master_node():
            Converter(filename=self.seedname).convert_dft_input()
            # lapw2 -almd writes one fort.225 record per band, l and m for every
            # atom and k-point through a unit with no .def entry.
            if os.path.isfile('fort.225'):
                os.remove('fort.225')
        mpi.barrier(poll_msec=100)

    # ---------------------------------------------------------------- the API

    def run_initial_stage(self, force_scf=False, **kwargs):
        """
        Converge the DFT SCF cycle, then build projectors and convert to HDF5.

        Restart-aware: a converged case.scf is reused and a partial one continued
        rather than recomputed, because modest calls this hook unconditionally and
        keeps no DFT state of its own across a resumed run.  Pass
        ``force_scf=True`` to redo the SCF from scratch regardless.
        """
        self._check_inputs()
        self._ensure_converged_scf(force_scf)
        self._regenerate_projectors()
        return 0

    def run_update_stage(self, N_k, Eint_m_dc, mu=0.0, beta=0.0, **kwargs):
        """
        Apply the DMFT charge density correction and reconverge one cycle.

        Runs lapw2 -qdmft (which rebuilds the valence density from ``N_k``),
        lcore and mixer to obtain the new density, then lapw0 and lapw1 to get
        the potential and eigenvectors that go with it, and finally rebuilds the
        projectors so the archive is current when this returns.

        Safe to call repeatedly with a fresh ``N_k``, which modest's CSC loop
        does several times per DMFT iteration.

        ``mu`` (eV) and ``beta`` (1/eV) are accepted only to fill their slots in
        case.qdmft, converted to Ry there.  lapw2 reads both and uses neither, so
        omitting them changes nothing.
        """
        self._write_qdmft(N_k, Eint_m_dc, mu=mu, beta=beta)
        mpi.barrier(poll_msec=100)

        # Serial by necessity: the DMFT path is compiled out of lapw2_mpi.
        self._run_x('lapw2', '-qdmft')
        self._run_x('lcore')
        self._mark_dmft_cycle()
        self._append_to_scf(_SCF_PARTS)
        self._save_old(_MIXER_SAVE)
        self._run_mixer()
        self._append_to_scf(('m',))

        self._save_old(_LAPW0_SAVE)
        self._run_x('lapw0')
        self._run_x('lapw1')
        self._regenerate_projectors()

        if mpi.is_master_node() and os.path.isfile('fort.77'):
            os.remove('fort.77')
        # Four decimals, not the full float repr: a single charge update carries
        # the solver's statistical error through N_k, so consecutive updates
        # scatter by a few meV even after the density has settled.  Printing 17
        # significant digits advertises a precision the number does not have.
        self._report(f"DFT + DMFT Total Energy: {self.read_dft_energy():.4f} eV")
        mpi.barrier(poll_msec=100)
        return 0
