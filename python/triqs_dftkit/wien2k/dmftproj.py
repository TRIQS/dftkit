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

"""Pure-Python drop-in for the dmftproj Fortran executable.

`run_dmftproj(case)` reads the Wien2k projector inputs (case.almblm{up,dn},
case.indmftpr, case.struct, case.dmftsym) and writes every case.* file the
Wien2k converter consumes: case.ctqmcout, case.symqmc, case.sympar,
case.parproj and case.oubwin{,up,dn}. With `band=True` it writes case.outband
from the band k-path (case.klist_band) instead, the dmftproj -band mode that
feeds convert_bands_input. Each output is produced by the corresponding
generator, all validated to machine precision against the precision-fixed
dmftproj.

Run as a script (`python -m triqs_dftkit.wien2k.dmftproj -sp -so [-band]`) it
mimics the Fortran binary: the case name is the working-directory basename.
"""

import os
import sys

from . import ctqmcout, outband, oubwin, parproj, symqmc, sympar


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


def main(argv=None):
    # The Fortran dmftproj takes the case name from the working directory; -sp
    # and -so are inferred from the input files, so they are accepted for
    # compatibility and ignored. -band selects the band-structure output.
    argv = sys.argv[1:] if argv is None else argv
    for flag in argv:
        if flag not in ('-sp', '-so', '-band'):
            raise SystemExit(f'dmftproj: unknown flag {flag}')
    case = os.path.basename(os.getcwd())
    run_dmftproj(case, band='-band' in argv)


if __name__ == '__main__':
    main()
