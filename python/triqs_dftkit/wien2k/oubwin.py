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

"""Pure-Python generation of the dmftproj band-window file (case.oubwin),
replacing the corresponding output of the dmftproj Fortran executable.

Given the alm/blm coefficient file (case.almblm, or the spin pair
case.almblmup / case.almblmdn) and the projector definition (case.indmftpr),
this selects, per k-point and spin, the contiguous band range whose Fermi-
shifted Kohn-Sham eigenvalues fall in the energy window, and writes the
result in the case.oubwin format read by the charge self-consistency.

It reproduces the Fortran path: read elecn/nk/nloat/eferm from the almblm
header (outbwin.f / dmftproj.f), parse the energy window (e_bot, e_top,
proj_mode) from the last line of case.indmftpr, shift every band eigenvalue
by eferm, and apply the set_projections.f selection. For proj_mode==0 this is
the energy-window rule (strict lower bound e_bot < E, inclusive upper bound
E <= e_top, yielding a single contiguous [nb_bot, nb_top] per k). For
proj_mode 1 and 2 the window is a pair of band indices (b_bot, b_top): every
k-point is included with nb_bot/nb_top those indices clamped to the local
nbmin/nbmax (set_projections.f:70-88). In mode 2 the indices come straight
from the window line; in mode 1 they are the global min/max band index over
all spins/k whose Fermi-shifted energy lies in (e_bot, e_top]
(dmftproj.f:704-722). The weight written per included k-point is the
tetrahedron weight of the lowest band nbmin.

Spin-polarization is detected from the input files: case.almblmup /
case.almblmdn present => two spin files (case.oubwinup / case.oubwindn),
otherwise the single case.oubwin. The spin-orbit flag (ifSO) comes from
case.indmftpr; with -so the up and dn windows must coincide and dmftproj
aborts otherwise, a check this generator also performs.
"""

import os

from ._dmftproj import (band_index_window, read_almblm, read_indmftpr,
                        select_band_window, select_window)


def _windows_from_spin(sp, info, band_window):
    """Per-k (included, nb_bot, nb_top, weight) for one already-read spin dict.
    proj_mode==0 uses the energy window; modes 1/2 use the band-index window
    (b_bot, b_top) precomputed in band_window."""
    out = []
    for kp in sp['kp']:
        if info['proj_mode'] == 0:
            incl, nb_bot, nb_top = select_window(
                kp['nbmin'], kp['nbmax'], kp['eband'],
                info['e_bot'], info['e_top'])
        else:
            incl, nb_bot, nb_top = select_band_window(
                kp['nbmin'], kp['nbmax'], band_window[0], band_window[1])
        out.append((incl, nb_bot, nb_top, kp['weight']))
    return out


def _write_oubwin_file(path, ifso, windows):
    """Write one case.oubwin spin file (outbwin.f). Integer fields are i6;
    the per-k weight is a list-directed real. Not-included k-points emit only
    the flag line."""
    with open(path, 'w') as f:
        f.write('%6d\n' % len(windows))
        f.write('%6d\n' % (1 if ifso else 0))
        for incl, nb_bot, nb_top, weight in windows:
            f.write('%6d\n' % (1 if incl else 0))
            if incl:
                f.write('%6d%6d\n' % (nb_bot, nb_top))
                f.write('   %.16f     \n' % weight)


def write_oubwin(case):
    """Read <case>.almblm{up,dn} (or <case>.almblm) and write the matching
    <case>.oubwin{up,dn} (or <case>.oubwin) in the dmftproj format.

    Spin-polarization is inferred from the input files. The spin-orbit flag
    written into every file comes from case.indmftpr; under -sp+-so the up/dn
    windows must coincide (dmftproj aborts otherwise)."""
    info = read_indmftpr(case + '.indmftpr')
    ifso = bool(info['so'])

    up, dn = case + '.almblmup', case + '.almblmdn'
    spin_polarized = os.path.exists(up) and os.path.exists(dn)

    if spin_polarized:
        sp_up = read_almblm(up, info)
        sp_dn = read_almblm(dn, info)
        # mode 1 scans both spins for the global band-index window.
        bw = (band_index_window(info, [sp_up, sp_dn])
              if info['proj_mode'] != 0 else None)
        win_up = _windows_from_spin(sp_up, info, bw)
        win_dn = _windows_from_spin(sp_dn, info, bw)
        if ifso:
            for u, d in zip(win_up, win_dn):
                if u[0] != d[0] or u[1] != d[1] or u[2] != d[2]:
                    raise ValueError(
                        'spin-orbit run requires identical up/dn band '
                        'windows at every k-point')
        _write_oubwin_file(case + '.oubwinup', ifso, win_up)
        _write_oubwin_file(case + '.oubwindn', ifso, win_dn)
    else:
        sp = read_almblm(case + '.almblm', info)
        bw = (band_index_window(info, [sp])
              if info['proj_mode'] != 0 else None)
        win = _windows_from_spin(sp, info, bw)
        _write_oubwin_file(case + '.oubwin', ifso, win)
