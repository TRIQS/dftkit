################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola, C. Martins
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

"""Shared dmftproj machinery for the pure-Python case.* generators.

The five generators (symqmc, oubwin, ctqmcout, sympar, parproj) reproduce
different outputs of the dmftproj Fortran executable but draw on the same
underlying pieces: the list-directed almblm reader, the Fortran real/complex
token parsers, the Wigner-D rotation in the dmftproj convention, the angular
basis transforms (transmat), the case.indmftpr / case.dmftsym parsers, the
proj_mode==0 band-window selection and the gfortran-style real formatters.
This module owns each of those exactly once.

A deliberate split survives in the angular basis: symqmc builds its symmetry
matrices from the EXACT cubic harmonics, while ctqmcout/sympar/parproj feed
projector/representation numbers through the single-precision CMPLX cast that
dmftproj applies in set_ang_trans.f:146 (0.70710676908 for 1/sqrt2, the
dft_tools #148 noise). `reptrans` exposes both via the `cast` flag; nothing
here "fixes" that truncation.
"""

import math
import os
import numpy as np


# --- list-directed token stream over a Fortran free-form text file -----------

class Reader:
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


def to_float(tok):
    """Parse a Fortran real token, accepting D/d exponents (1.0D-3) and the
    list-directed 4.57E-002 form. Python float already handles E/e and the
    embedded sign in the exponent; only D/d needs translation."""
    return float(tok.replace('D', 'E').replace('d', 'e'))


def to_complex(s):
    s = s.strip().lstrip('(').rstrip(')')
    re_s, im_s = s.split(',')
    return complex(to_float(re_s), to_float(im_s))


def read_complex(r):
    """Read one list-directed complex `(re,im)` value, possibly split over
    whitespace tokens."""
    buf = list(r.record())
    while ')' not in ''.join(buf):
        buf += r.record()
    return to_complex(''.join(buf))


def read_two_complex(r):
    """Read a record holding two list-directed complex values (Alm, Blm)."""
    buf = list(r.record())
    while ''.join(buf).count(')') < 2:
        buf += r.record()
    s = ''.join(buf)
    cut = s.index(')') + 1
    return to_complex(s[:cut]), to_complex(s[cut:])


# --- angular bases: transmat = <new_i|m>, m = -l..l --------------------------
# Standard cubic harmonics, Wien2k convention. dmftproj stores cubic/fromfile
# coefficients via a single-precision CMPLX cast (set_ang_trans.f:146); the
# projector/representation generators reproduce that float32 truncation
# (0.70710676908 for 1/sqrt2), while the symqmc symmetry matrices use the
# exact harmonics. `reptrans(..., cast=...)` selects between the two.

_COMPLEX = {l: np.eye(2 * l + 1, dtype=complex) for l in range(4)}

_CUBIC = {
    1: np.array([[0, 1, 0], [-1j, 0, -1j], [1, 0, -1]], dtype=complex),
    2: np.array([
        [0, 0, 1, 0, 0],
        [2 ** -0.5, 0, 0, 0, 2 ** -0.5],
        [-(2 ** -0.5), 0, 0, 0, 2 ** -0.5],
        [0, 2 ** -0.5, 0, -(2 ** -0.5), 0],
        [0, 2 ** -0.5, 0, 2 ** -0.5, 0],
    ], dtype=complex),
}


def reptrans(basis, l, cast=True):
    """transmat = <new_i|m>. `cast=True` truncates the cubic harmonics through
    the single-precision CMPLX cast dmftproj applies (set_ang_trans.f:146);
    `cast=False` keeps them exact. A complex basis is the exact identity."""
    if basis == 'cubic' and l in _CUBIC:
        c = _CUBIC[l]
        return c.astype(np.complex64).astype(np.complex128) if cast else c
    return _COMPLEX[l]


def read_fromfile(path, l):
    """Parse a dmftproj fromfile basis (each line: the coefficients of a new
    basis vector in {|m,up>, |m,dn>}, m = -l..l, real/imag interleaved; '*'
    marks the end of an irep). Return (P, mixing) where P = <new|m> is the
    transform matrix (reptrans.transmat), full 2(2l+1) for a spin-mixing basis
    or the (2l+1) up/up block otherwise."""
    n = 2 * (2 * l + 1)
    rows = []
    for line in open(path):
        line = line.rstrip('\n')
        if not line.strip():
            continue
        body = line[1:]
        vals = [float(x) for x in body.split()][:2 * n]
        rows.append([vals[2 * k] + 1j * vals[2 * k + 1] for k in range(n)])
        if len(rows) == n:
            break
    R = np.array(rows)                       # R[i, :] = |new_i> in old basis
    d = 2 * l + 1
    up_up, up_dn = R[:d, :d], R[:d, d:]
    dn_up, dn_dn = R[d:, :d], R[d:, d:]
    mixing = not (np.allclose(dn_dn, up_up) and
                  np.allclose(up_dn, 0) and np.allclose(dn_up, 0))
    P = np.conj(R) if mixing else np.conj(up_up)   # <new|m> = conj(<m|new>)
    return P, mixing


# --- Wigner D matrix (dmftproj convention, setsym.f) -------------------------

def small_d(l, m, n, b):
    f1 = (math.factorial(l + m) * math.factorial(l - m)) / \
         (math.factorial(l + n) * math.factorial(l - n))
    s = 0.0
    for t in range(0, 2 * l + 1):
        if (l - m - t) >= 0 and (l - n - t) >= 0 and (t + n + m) >= 0:
            f2 = (math.factorial(l + n) * math.factorial(l - n)) / \
                 (math.factorial(l - m - t) * math.factorial(m + n + t) *
                  math.factorial(l - n - t) * math.factorial(t))
            f3 = 1.0 if (2 * l - m - n - 2 * t) == 0 \
                else math.sin(b / 2) ** (2 * l - m - n - 2 * t)
            f4 = 1.0 if (2 * t + n + m) == 0 \
                else math.cos(b / 2) ** (2 * t + n + m)
            s += (-1) ** (l - m - t) * f2 * f3 * f4
    return math.sqrt(f1) * s


def dmat(l, a, b, c, det):
    D = np.zeros((2 * l + 1, 2 * l + 1), dtype=complex)
    for m in range(-l, l + 1):
        for n in range(-l, l + 1):
            v = np.exp(1j * n * a) * np.exp(1j * m * c) * small_d(l, m, n, b)
            if det < -0.5:
                v *= (-1) ** l
            D[m + l, n + l] = v
    return D


def tmat(l):
    """Complex-conjugation operator in the spherical-harmonic basis,
    T[m, -m] = (-1)^m (timeinv.f)."""
    T = np.zeros((2 * l + 1, 2 * l + 1), dtype=complex)
    for m in range(-l, l + 1):
        T[-m + l, m + l] = (-1) ** m
    return T


def timeinv_orbital(l, mat):
    return tmat(l) @ np.conj(mat)


# --- case.indmftpr -----------------------------------------------------------

def read_indmftpr(indmftpr):
    """Full structured parse of case.indmftpr.

    Returns a dict with nsort, mult, lmax, the spin-orbit flag `so`, the energy
    window (e_bot, e_top, proj_mode) and one `sorts` entry per atomic sort with
    its basis name, optional fromfile sourcefile path, and the correlated /
    included l lists (l_inc==2 / l_inc in {1,2}). Every generator derives its
    shell or orbital list from this one parse."""
    raw = [l.split('!')[0].strip() for l in open(indmftpr)]
    raw = [l for l in raw if l != '']
    nsort = int(raw[0].split()[0])
    mult = [int(x) for x in raw[1].split()][:nsort]
    lmax = int(raw[2].split()[0])
    i = 3
    so = 0
    sorts = []
    for isort in range(nsort):
        basis = raw[i].split()[0]
        i += 1
        sourcefile = None
        if basis == 'fromfile':
            sourcefile = os.path.join(os.path.dirname(indmftpr), raw[i])
            i += 1
        l_inc = [int(x) for x in raw[i].split()]
        i += 1
        ireps = [int(x) for x in raw[i].split()]
        i += 1
        correlated_ls = [l for l in range(len(l_inc)) if l_inc[l] == 2]
        included_ls = [l for l in range(len(l_inc)) if l_inc[l] in (1, 2)]
        if any(n > 0 for n in ireps):
            i += 1                       # skip the correps line
        if correlated_ls:
            so = int(raw[i].split()[0])  # SO flag follows a correlated sort
            i += 1
        sorts.append(dict(basis=basis, sourcefile=sourcefile,
                          correlated_ls=correlated_ls, included_ls=included_ls))
    last = raw[-1].split()
    e_bot, e_top = to_float(last[0]), to_float(last[1])
    proj_mode = int(last[2]) if len(last) >= 3 else 0
    return dict(nsort=nsort, mult=mult, lmax=lmax, sorts=sorts, so=so,
                e_bot=e_bot, e_top=e_top, proj_mode=proj_mode)


# --- case.dmftsym ------------------------------------------------------------

def read_dmftsym(path, rotloc=False):
    """Symmetry operations from case.dmftsym. Each op carries (perm, a, b, g,
    iprop, krotm); the third Euler angle is named both `c` and `g` for the
    callers that use either name. With rotloc=True also parse the
    'Global->local' representative rotloc per sort and return it as a third
    value (setsym.f:437-455)."""
    lines = open(path).read().split('\n')
    nsym = int(lines[0].split()[0])
    perms = [[int(x) for x in lines[1 + i].split()] for i in range(nsym)]
    starts = [i for i, l in enumerate(lines) if 'Sym. op.' in l]
    ops = []
    for k, s in enumerate(starts[:nsym]):
        toks = lines[s + 1].split()
        a, b, g = (math.radians(float(x)) for x in toks[:3])
        iprop = int(toks[3])
        krotm = np.array([[to_float(x) for x in lines[s + 2 + r].split()]
                          for r in range(3)])
        ops.append(dict(perm=perms[k], a=a, b=b, c=g, g=g, iprop=iprop,
                        krotm=krotm))
    if not rotloc:
        return nsym, ops

    gl = next(i for i, l in enumerate(lines) if 'Global->local' in l)
    rotloc_ref = []
    j = gl + 1
    while len(rotloc_ref) < nsym and j < len(lines):
        if lines[j].strip() == '' or not lines[j].split()[0].lstrip('-').isdigit():
            j += 1
            continue
        # sort index line, then 3 krotm rows, then the Euler/iprop line.
        if len(lines[j].split()) == 1:
            krotm = np.array([[to_float(x) for x in lines[j + 1 + r].split()]
                              for r in range(3)])
            ang = lines[j + 4].split()
            a, b, g = (math.radians(float(x)) for x in ang[:3])
            iprop = int(ang[3])
            rotloc_ref.append(dict(krotm=krotm, a=a, b=b, g=g, iprop=iprop))
            j += 5
        else:
            j += 1
    return nsym, ops, rotloc_ref


# --- almblm parsing ----------------------------------------------------------

def read_almblm(path, info, projectors=False):
    """Read one spin's almblm file in the dmftproj.f read order.

    With projectors=False keep only the header (elecn, nk, nloat, eferm), the
    per-(sort,l) overlap block (u_dot_norm, nLO, ovl_LO_u, ovl_LO_udot) and the
    per-k band data (nbmin/nbmax, Fermi-shifted eband, tetrahedron weights); the
    per-(l,m) coefficient records are still consumed to keep the reader in sync
    but their values are discarded (oubwin's window-only path).

    With projectors=True also accumulate the per-(l,m) Alm/Blm/Clm coefficients
    in the canonical Wien2k packing lm = l*l + (m+l), for every atom of every
    sort (ctqmcout/parproj)."""
    nsort = info['nsort']
    lmax = info['lmax']
    mult = info['mult']
    nlm = (lmax + 1) ** 2
    natom = sum(mult)

    r = Reader(path)
    elecn = to_float(r.record()[0])
    nk = int(r.record()[0])
    r.record()                         # nloat
    eferm = to_float(r.record()[0])

    nLO = {}
    u_dot_norm = {}
    ovl_LO_u = {}
    ovl_LO_udot = {}
    kp = [None] * nk

    for isrt in range(1, nsort + 1):
        for l in range(lmax + 1):
            u_dot_norm[(l, isrt)] = to_float(r.record()[0])
            n = int(r.record()[0])
            nLO[(l, isrt)] = n
            for ilo in range(1, n + 1):
                toks = r.record()
                ovl_LO_u[(ilo, l, isrt)] = to_float(toks[0])
                ovl_LO_udot[(ilo, l, isrt)] = to_float(toks[1])
        for ik in range(nk):
            r.skip()                   # "IK = .." banner
            r.skip()                   # 3-int line
            head = r.record()
            nbmin, nbmax = int(head[1]), int(head[2])
            nb = nbmax - nbmin + 1
            if kp[ik] is None:
                kp[ik] = dict(nbmin=nbmin, nbmax=nbmax, eband=None,
                              weight=None, tetr=None)
                if projectors:
                    kp[ik].update(
                        Alm=np.zeros((nlm, natom + 1, nb), dtype=complex),
                        Blm=np.zeros((nlm, natom + 1, nb), dtype=complex),
                        Clm=np.zeros((4, nlm, natom + 1, nb), dtype=complex))
            eband = np.empty(nb)
            tetr = np.empty(nb)
            for off in range(nb):
                toks = r.record()
                tetr[off] = to_float(toks[0])
                eband[off] = to_float(toks[1])
            eband = eband - eferm
            if kp[ik]['eband'] is None:
                kp[ik]['eband'] = eband
                kp[ik]['weight'] = tetr[0]
                kp[ik]['tetr'] = tetr
            for imu in range(1, mult[isrt - 1] + 1):
                iatom = sum(mult[:isrt - 1]) + imu
                r.skip()               # banner
                r.record()             # idum
                for off in range(nb):
                    lm = 0
                    for l in range(lmax + 1):
                        for m in range(-l, l + 1):
                            alm, blm = read_two_complex(r)
                            for ilo in range(nLO[(l, isrt)]):
                                clm = read_complex(r)
                                if projectors:
                                    kp[ik]['Clm'][ilo, lm, iatom, off] = clm
                            if projectors:
                                kp[ik]['Alm'][lm, iatom, off] = alm
                                kp[ik]['Blm'][lm, iatom, off] = blm
                            lm += 1

    return dict(elecn=elecn, eferm=eferm, nk=nk, nlm=nlm, natom=natom,
                nLO=nLO, u_dot_norm=u_dot_norm, ovl_LO_u=ovl_LO_u,
                ovl_LO_udot=ovl_LO_udot, kp=kp)


# --- band-window selection (set_projections.f, proj_mode==0) -----------------

def select_window(nbmin, nbmax, eband, e1, e2):
    """proj_mode==0 contiguous band selection for one k-point over the energy
    array `eband` (already Fermi-shifted). The first band with e1 < E <= e2
    opens the window; the first band above it with E > e2 closes it at the
    preceding index; a window reaching nbmax closes at nbmax. Returns
    (included, nb_bot, nb_top) with nb_bot=nb_top=0 when no band qualifies."""
    included = False
    nb_bot = nb_top = 0
    for off, ib in enumerate(range(nbmin, nbmax + 1)):
        e = eband[off]
        if not included and e > e1 and e <= e2:
            included = True
            nb_bot = ib
        elif included and e > e2:
            nb_top = ib - 1
            break
        elif ib == nbmax and e > e1 and e <= e2:
            nb_top = ib
            included = True
    if not included:
        nb_bot = nb_top = 0
    return included, nb_bot, nb_top


# --- gfortran-style real formatters ------------------------------------------

def fmt(x):
    """One list-directed real, gfortran-style (leading sign space, ~17 sig)."""
    return '  %.16E' % float(x)


def write_row(f, arr):
    f.write(''.join(fmt(x) for x in arr) + '\n')


def fmt_scalar(x):
    x = float(x)
    return '%.16f' % x if abs(x) < 1e5 else '%.16E' % x
