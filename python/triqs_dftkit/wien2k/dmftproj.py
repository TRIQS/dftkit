################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by L. Pourovskii, V. Vildosola, C. Martins, M. Aichhorn
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

"""Pure-Python driver for the dmftproj file interface.

`run_dmftproj(case)` reads the Wien2k projector inputs (case.almblm{up,dn},
case.indmftpr, case.struct, case.dmftsym) and writes every case.* file the
Wien2k converter consumes: case.ctqmcout, case.symqmc, case.sympar,
case.parproj and case.oubwin{,up,dn}. With `band=True` it writes case.outband
from the band k-path (case.klist_band) instead, the dmftproj -band mode that
feeds convert_bands_input. Each output is produced by the corresponding
generator, all validated to machine precision against native
dmftproj.

Run as a script (`python -m triqs_dftkit.wien2k.dmftproj -sp -so [-band]`) it
mimics the Fortran binary: the case name is the working-directory basename.
`run_dmftproj` infers spin polarization from the presence of
case.almblm{up,dn} and spin-orbit from case.indmftpr; the command checks that
-sp and -so agree with these inputs and exits with an error otherwise.

Unsupported inputs raise explicitly: selected irreducible subspaces, cubic f
harmonics, multiple local orbitals per angular momentum, spin-mixing bases
without SO, differing SO flags across correlated sorts, and non-SO band paths.
Band paths currently require projection mode 0. Two-spin projector generation
requires matching selected and below-window bands in both spin channels.
"""

import os
import sys

from . import ctqmcout, outband, oubwin, parproj, symqmc, sympar
from ._dmftproj import read_indmftpr


def run_dmftproj(case, band=False):
    """Write the dmftproj output files for `case` (a path prefix). The band
    mode (dmftproj -band) writes case.outband; otherwise the SCF projector
    files the converter's convert_dft_input / convert_parproj_input read."""
    if band:
        outband.write_outband(case)
        return
    symqmc.write_symqmc(case)
    ctqmcout.write_ctqmcout(case)
    sympar.write_sympar(case)
    parproj.write_parproj(case)
    oubwin.write_oubwin(case)


def check_flags(case, sp, so, band=False):
    """Raise ValueError unless the dmftproj -sp/-so flags describe the inputs
    of `case` the way run_dmftproj reads them: -sp <=> case.almblm{up,dn}
    both present, -so <=> the case.indmftpr spin-orbit flag, -so => -sp."""
    if so and not sp:
        raise ValueError('-so requires spin-polarized input files and -sp')
    spin_files = all(os.path.exists(f'{case}.almblm{s}') for s in ('up', 'dn'))
    if sp and not spin_files:
        raise ValueError(f'-sp given but {case}.almblmup / {case}.almblmdn not found')
    if not sp and spin_files:
        raise ValueError(f'{case}.almblmup and {case}.almblmdn present without -sp; '
                         f'pass -sp, or remove them to use {case}.almblm')
    if band:
        info = outband._read_indmftpr_band(case + '.indmftpr')[0]
    else:
        info = read_indmftpr(case + '.indmftpr')
    if bool(info['so']) != so:
        raise ValueError(f'{case}.indmftpr has spin-orbit flag {info["so"]}, '
                         f'but -so was {"given" if so else "not given"}')


def main(argv=None):
    # Same command line as the Fortran dmftproj: the case name is the working
    # directory basename, flags are -sp, -so and -band.
    argv = sys.argv[1:] if argv is None else argv
    for flag in argv:
        if flag not in ('-sp', '-so', '-band'):
            raise SystemExit(f'dmftproj: unknown flag {flag}')
    case = os.path.basename(os.getcwd())
    band = '-band' in argv
    try:
        check_flags(case, '-sp' in argv, '-so' in argv, band=band)
    except ValueError as error:
        raise SystemExit(f'dmftproj: {error}')
    run_dmftproj(case, band=band)


if __name__ == '__main__':
    main()
