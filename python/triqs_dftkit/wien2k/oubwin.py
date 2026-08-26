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
by eferm, and apply the proj_mode==0 selection of set_projections.f
(strict lower bound e_bot < E, inclusive upper bound E <= e_top, yielding a
single contiguous [nb_bot, nb_top] per k). The weight written per included
k-point is the tetrahedron weight of the lowest band nbmin.

Spin-polarization is detected from the input files: case.almblmup /
case.almblmdn present => two spin files (case.oubwinup / case.oubwindn),
otherwise the single case.oubwin. The spin-orbit flag (ifSO) comes from
case.indmftpr; with -so the up and dn windows must coincide and dmftproj
aborts otherwise, a check this generator also performs.
"""

import os
import numpy as np


# --- list-directed token stream over a Fortran free-form text file -----------

class _Reader:
    """Record-oriented reader mirroring Fortran list-directed and READ('()')
    semantics on a free-form almblm file.

    `record()` returns the whitespace tokens of the next physical line (the
    common case: every almblm logical record fits on one line). `skip()`
    consumes one physical line unconditionally, matching READ('()')."""

    def __init__(self, path):
        with open(path) as fh:
            self._lines = fh.read().splitlines()
        self._i = 0

    def skip(self):
        self._i += 1

    def record(self):
        toks = self._lines[self._i].split()
        self._i += 1
        return toks


def _to_float(tok):
    """Parse a Fortran real token, accepting D/d exponents (1.0D-3) and the
    list-directed 4.57E-002 form. Python float already handles E/e and the
    embedded sign in the exponent; only D/d needs translation."""
    return float(tok.replace('D', 'E').replace('d', 'e'))


# --- inputs ------------------------------------------------------------------

def _read_window(indmftpr):
    """Energy/band window from the last non-empty line of case.indmftpr.

    Try list-directed `e_bot e_top proj_mode`; on failure fall back to
    `e_bot e_top` with proj_mode=0 (dmftproj.f:182-197)."""
    lines = [l for l in open(indmftpr).read().splitlines() if l.strip()]
    toks = lines[-1].split()
    e_bot, e_top = _to_float(toks[0]), _to_float(toks[1])
    proj_mode = int(toks[2]) if len(toks) >= 3 else 0
    return e_bot, e_top, proj_mode


def _read_indmftpr_geometry(indmftpr):
    """Return (nsort, lmax) from case.indmftpr. nsort is its first value and
    lmax its third (dmftproj.f:95). Both are needed to walk the per-sort
    almblm header up to the k-loop."""
    raw = [l.split('!')[0].strip() for l in open(indmftpr)]
    raw = [l for l in raw if l != '']
    nsort = int(raw[0].split()[0])
    lmax = int(raw[2].split()[0])
    return nsort, lmax


def _read_so_flag(indmftpr):
    """Spin-orbit flag from case.indmftpr: the single integer line that
    follows a correlated sort's irep specification (the value dmftproj stores
    as ifSOflag). Returns 0 or 1."""
    raw = [l.split('!')[0].strip() for l in open(indmftpr)]
    raw = [l for l in raw if l != '']
    nsort = int(raw[0].split()[0])
    i = 3
    so = 0
    for _ in range(nsort):
        basis = raw[i].split()[0]
        i += 1
        if basis == 'fromfile':
            i += 1
        l_inc = [int(x) for x in raw[i].split()]
        i += 1
        ireps = [int(x) for x in raw[i].split()]
        i += 1
        correlated = any(v == 2 for v in l_inc)
        if any(n > 0 for n in ireps):
            i += 1
        if correlated:
            so = int(raw[i].split()[0])
            i += 1
    return so


# --- almblm parsing ----------------------------------------------------------

def _read_almblm_windows(path, nsort, lmax):
    """Read case.almblm and return (eferm, [(nbmin, nbmax, eband, weight) per
    k-point]) for the FIRST atomic sort. nbmin/nbmax/eband/weight are identical
    across sorts at a given k, so the first sort's k-loop suffices for window
    selection.

    Follows the dmftproj.f read order field by field: header, then the first
    sort's per-l header (u_dot_norm, nLO, nLO overlap lines), then its k-loop
    (two skipped banner lines, the idum/nbmin/nbmax line, the per-band
    rtetr/eband lines)."""
    r = _Reader(path)
    r.record()                       # elecn
    nk = int(r.record()[0])
    r.record()                       # nloat
    eferm = _to_float(r.record()[0])  # non-band path: eferm is header line 4

    # First sort's per-l header: u_dot_norm, nLO, then nLO overlap lines.
    for _ in range(lmax + 1):
        r.record()                   # u_dot_norm(l)
        nlo = int(r.record()[0])
        for _ in range(nlo):
            r.record()               # ovl_LO_u, ovl_LO_udot

    kpoints = []
    for _ in range(nk):
        r.skip()                     # "IK = .. jatom=.." banner
        r.skip()                     # 3-int line
        idum_nb = r.record()
        nbmin, nbmax = int(idum_nb[1]), int(idum_nb[2])
        eband = np.empty(nbmax - nbmin + 1)
        weight = None
        for off in range(nbmax - nbmin + 1):
            toks = r.record()
            rtetr = _to_float(toks[0])
            eband[off] = _to_float(toks[1])
            if off == 0:             # tetrahedron weight of the lowest band
                weight = rtetr
        kpoints.append((nbmin, nbmax, eband, weight))
    return eferm, kpoints


# --- band-window selection (set_projections.f, proj_mode==0) -----------------

def _select_window(nbmin, nbmax, eband, eferm, e_bot, e_top):
    """proj_mode==0 contiguous band selection for one k-point.

    Scan absolute band index ib = nbmin..nbmax over the Fermi-shifted energy
    E = eband - eferm. The first ib with e_bot < E <= e_top opens the window;
    the first ib above it with E > e_top closes it at ib-1 (scan exits). A
    window reaching nbmax closes at nbmax. Returns (included, nb_bot, nb_top)
    with nb_bot=nb_top=0 when no band qualifies."""
    included = False
    nb_bot = nb_top = 0
    for off, ib in enumerate(range(nbmin, nbmax + 1)):
        e = eband[off] - eferm
        if not included and e_bot < e <= e_top:
            included = True
            nb_bot = ib
        elif included and e > e_top:
            nb_top = ib - 1
            break
        elif ib == nbmax and e_bot < e <= e_top:
            nb_top = ib
            included = True
    if not included:
        nb_bot = nb_top = 0
    return included, nb_bot, nb_top


def _windows_for_spin(almblm, indmftpr, nsort, lmax, e_bot, e_top, proj_mode):
    """Per-k (included, nb_bot, nb_top, weight) for one spin's almblm file."""
    if proj_mode != 0:
        raise NotImplementedError(
            'oubwin generation is implemented for the energy-window '
            'projection mode (proj_mode==0); band-index modes 1 and 2 are '
            'not yet covered by a test fixture')
    eferm, kpoints = _read_almblm_windows(almblm, nsort, lmax)
    out = []
    for nbmin, nbmax, eband, weight in kpoints:
        incl, nb_bot, nb_top = _select_window(
            nbmin, nbmax, eband, eferm, e_bot, e_top)
        out.append((incl, nb_bot, nb_top, weight))
    return out


# --- output ------------------------------------------------------------------

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
    indmftpr = case + '.indmftpr'
    nsort, lmax = _read_indmftpr_geometry(indmftpr)
    ifso = _read_so_flag(indmftpr)
    e_bot, e_top, proj_mode = _read_window(indmftpr)

    up, dn = case + '.almblmup', case + '.almblmdn'
    spin_polarized = os.path.exists(up) and os.path.exists(dn)

    if spin_polarized:
        win_up = _windows_for_spin(
            up, indmftpr, nsort, lmax, e_bot, e_top, proj_mode)
        win_dn = _windows_for_spin(
            dn, indmftpr, nsort, lmax, e_bot, e_top, proj_mode)
        if ifso:
            for u, d in zip(win_up, win_dn):
                if u[0] != d[0] or u[1] != d[1] or u[2] != d[2]:
                    raise ValueError(
                        'spin-orbit run requires identical up/dn band '
                        'windows at every k-point')
        _write_oubwin_file(case + '.oubwinup', ifso, win_up)
        _write_oubwin_file(case + '.oubwindn', ifso, win_dn)
    else:
        win = _windows_for_spin(
            case + '.almblm', indmftpr, nsort, lmax, e_bot, e_top, proj_mode)
        _write_oubwin_file(case + '.oubwin', ifso, win)
