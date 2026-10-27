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

All generators use the exact double-precision cubic harmonics; with
full-precision SRC_templates dmftproj reads the same coefficients, so the port
reproduces its output to machine precision.
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
# Standard cubic harmonics, Wien2k convention, exact double precision.

_COMPLEX = {l: np.eye(2 * l + 1, dtype=complex) for l in range(4)}

_CUBIC = {
    1: np.array([[2**-0.5, 0, -2**-0.5], [-1j*2**-0.5, 0, -1j*2**-0.5], [0, 1, 0]], dtype=complex),
    2: np.array([
        [0, 0, 1, 0, 0],
        [2 ** -0.5, 0, 0, 0, 2 ** -0.5],
        [-(2 ** -0.5), 0, 0, 0, 2 ** -0.5],
        [0, 2 ** -0.5, 0, -(2 ** -0.5), 0],
        [0, 2 ** -0.5, 0, 2 ** -0.5, 0],
    ], dtype=complex),
}


def reptrans(basis, l):
    """transmat = <new_i|m>: the exact double-precision cubic harmonics for a
    cubic basis, the identity for a complex basis."""
    if basis == 'cubic' and l in _CUBIC:
        return _CUBIC[l]
    if basis == 'cubic' and l != 0:
        raise NotImplementedError('cubic f harmonics are not implemented yet')
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
        vals = [to_float(x) for x in body.split()[:2 * n]]
        if len(vals) != 2 * n:
            raise ValueError(f'{path}: each basis row requires {2 * n} real coefficients')
        rows.append([vals[2 * k] + 1j * vals[2 * k + 1] for k in range(n)])
        if len(rows) == n:
            break
    R = np.array(rows)                       # R[i, :] = |new_i> in old basis
    if R.shape != (n, n) or not np.all(np.isfinite(R)):
        raise ValueError(f'{path}: expected a finite {n} by {n} basis matrix')
    d = 2 * l + 1
    up_up, up_dn = R[:d, :d], R[:d, d:]
    dn_up, dn_dn = R[d:, :d], R[d:, d:]
    mixing = not (np.array_equal(dn_dn, up_up) and
                  np.all(up_dn == 0) and np.all(dn_up == 0))
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


def mixing_timeinv_op(l, P):
    """The spinor time-reversal operator -i sigma_y (x) T in the mixing basis:
    tinv_{new} = P tinv_{lm} P^T (timeinv.f, ifmixing branch). Returned as the
    bare operator; callers apply it as tinv @ conj(mat)."""
    d = 2 * l + 1
    tm = tmat(l)
    tinv = np.zeros((2 * d, 2 * d), dtype=complex)
    tinv[:d, d:] = -tm
    tinv[d:, :d] = tm
    return P @ tinv @ P.T


def rotloc_rotl_so(l, ref, ops, iatom, iref):
    """The composed 2(2l+1) Rloc rotation rotloc(iatom)%rotl(l) under SP+SO
    (setsym.f:496-528 + set_rotloc.f): the representative spinor rotloc
    spmt (x) D(rotloc_ref) composed with the first symmetry op R[isym] mapping
    iref onto iatom, with the orbital time-reversal applied for the magnetic op.
    Returns (rotl, timeinv); callers apply their own basis transform (the full
    transmat for a mixing basis, blkdiag(transmat, transmat) otherwise)."""
    d = 2 * l + 1
    Dref = dmat(l, ref['a'], ref['b'], ref['g'], float(ref['iprop']))
    f = (ref['a'] + ref['g']) / 2.0
    spmt = np.zeros((2, 2), dtype=complex)
    spmt[0, 0] = np.exp(1j * f) * math.cos(ref['b'] / 2.0)
    spmt[1, 1] = np.conj(spmt[0, 0])
    f = -(ref['a'] - ref['g']) / 2.0
    spmt[0, 1] = np.exp(1j * f) * math.sin(ref['b'] / 2.0)
    spmt[1, 0] = -np.conj(spmt[0, 1])
    rotl = np.zeros((2 * d, 2 * d), dtype=complex)
    rotl[:d, :d] = spmt[0, 0] * Dref
    rotl[d:, d:] = spmt[1, 1] * Dref
    rotl[:d, d:] = spmt[0, 1] * Dref
    rotl[d:, :d] = spmt[1, 0] * Dref

    op = next(o for o in ops if o['perm'][iref - 1] == iatom)
    det2 = (op['krotm'][0, 0] * op['krotm'][1, 1]
            - op['krotm'][0, 1] * op['krotm'][1, 0])
    timeinv = det2 < 0.0
    srot_phase = (op['g'] - op['a']) if timeinv else (op['a'] + op['g'])
    rotl_sym = dmat(l, op['a'], op['b'], op['g'], float(op['iprop']))
    if timeinv:
        rotl_sym = tmat(l) @ np.conj(rotl_sym)
    ephase = np.exp(1j * srot_phase / 2.0)
    tmp = np.zeros((2 * d, 2 * d), dtype=complex)
    tmp[:d, :d] = ephase * rotl_sym
    tmp[d:, d:] = np.conj(ephase) * rotl_sym
    rotl = (tmp @ np.conj(rotl)) if timeinv else (tmp @ rotl)
    return rotl, timeinv


def mixing_rotrep(op, l, P, ti):
    """The full 2(2l+1) spinor representation D(R)_{new_i} = P spinrot P^dag of
    one symmetry operation in a spin-coupling (mixing) basis, with the spinor
    time-reversal operator applied for the magnetic (timeinv) operations
    (setsym.f spinrotmat + timeinv_op). This is srot%rotrep(l,isrt)%mat, shared
    by the symqmc and sympar shell matrices."""
    rotl = dmat(l, op['a'], op['b'], op['c'], float(op['iprop']))
    phase = (op['c'] - op['a']) if ti else (op['a'] + op['c'])
    e = np.exp(1j * phase / 2)
    d = 2 * l + 1
    spinrot = np.zeros((2 * d, 2 * d), dtype=complex)
    if ti:                                      # beta = pi, block-antidiagonal
        spinrot[:d, d:] = e * rotl
        spinrot[d:, :d] = -np.conj(e) * rotl
    else:                                       # beta = 0, block-diagonal
        spinrot[:d, :d] = e * rotl
        spinrot[d:, d:] = np.conj(e) * rotl
    rotrep = P @ spinrot @ np.conj(P.T)
    if ti:
        rotrep = mixing_timeinv_op(l, P) @ np.conj(rotrep)
    return rotrep


# --- shell symmetry matrices (shared by symqmc and sympar) -------------------

def _phase(op, ti):
    a, c = op['a'], op['c']
    return (c - a) if ti else (a + c)       # (g-a) on magnetic ops, else (a+g)


def _l0_matrix(op, ti):
    """s shell: 2x2 spin phase block diag(e, conj(e)) (outputqmc.f l==0)."""
    e = np.exp(1j * _phase(op, ti) / 2)
    return np.array([[e, 0], [0, np.conj(e)]], dtype=complex)


def _nonmixing_matrix(op, shell, ti):
    """Non-mixing SP+SO whole shell: the up/up block scaled by +-(a+g)/2,
    block-diagonal over spin, with the orbital time-reversal operator on the
    magnetic operations (setsym.f, outputqmc.f:1292-1339)."""
    l, P = shell['l'], shell['P']
    rotl = dmat(l, op['a'], op['b'], op['c'], float(op['iprop']))
    if ti:
        rotl = timeinv_orbital(l, rotl)
    rotrep = P @ rotl @ np.conj(P.T)
    e = np.exp(1j * _phase(op, ti) / 2)
    d = 2 * l + 1
    mat = np.zeros((2 * d, 2 * d), dtype=complex)
    mat[:d, :d] = e * rotrep
    mat[d:, d:] = np.conj(e) * rotrep
    return mat


def shell_matrix(op, shell, ti):
    """Spinor symmetry matrix for one shell under one operation."""
    l = shell['l']
    if l == 0:
        return _l0_matrix(op, ti)
    if shell['mixing']:
        return mixing_rotrep(op, l, shell['P'], bool(ti))
    return _nonmixing_matrix(op, shell, ti)


def build_orbs(info):
    """One entry per INCLUDED orbital (sort, l, atom order)."""
    orbs = []
    for isort in range(1, info['nsort'] + 1):
        sortinfo = info['sorts'][isort - 1]
        for l in sortinfo['included_ls']:
            for imu in range(1, info['mult'][isort - 1] + 1):
                atom = sum(info['mult'][:isort - 1]) + imu
                orbs.append(dict(atom=atom, sort=isort, l=l,
                                 basis=sortinfo['basis']))
    return orbs


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
        for l in correlated_ls:
            if ireps[l] > 1:
                raise NotImplementedError('selection of individual irreducible representations is not implemented')
            if ireps[l] > 0:
                if raw[i].strip() != '1':
                    raise NotImplementedError('selection of individual irreducible representations is not implemented')
                i += 1                   # each correlated l has its own correps line
        if correlated_ls:
            sort_so = int(raw[i].split()[0])
            if sort_so not in (0, 1):
                raise ValueError('spin-orbit flags must be 0 or 1')
            if any(s['correlated_ls'] for s in sorts) and sort_so != so:
                raise NotImplementedError('different spin-orbit flags across correlated sorts are not implemented')
            so = sort_so                 # SO flag follows a correlated sort
            i += 1
        sorts.append(dict(basis=basis, sourcefile=sourcefile,
                          correlated_ls=correlated_ls, included_ls=included_ls))
    last = raw[-1].split()
    e_bot, e_top = to_float(last[0]), to_float(last[1])
    proj_mode = int(last[2]) if len(last) >= 3 else 0
    if proj_mode not in (0, 1, 2) or not np.isfinite([e_bot, e_top]).all() or e_bot > e_top:
        raise ValueError('projection window requires finite ordered bounds and mode 0, 1 or 2')
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
    while j < len(lines):
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
            if n > 1:
                raise NotImplementedError('more than one local orbital per angular momentum is not implemented')
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
    opened = False
    included = False
    nb_bot = nb_top = 0
    for off, ib in enumerate(range(nbmin, nbmax + 1)):
        e = eband[off]
        if not opened and e > e1 and e <= e2:
            opened = True
            nb_bot = ib
        elif opened and e > e2:
            nb_top = ib - 1
            included = True
            break
        elif ib == nbmax and e > e1 and e <= e2:
            nb_top = ib
            included = True
    if not included:
        nb_bot = nb_top = 0
    return included, nb_bot, nb_top


def select_band_window(nbmin, nbmax, b_bot, b_top):
    """proj_mode 1/2 band-index selection for one k-point (set_projections.f
    70-88). e1/e2 are band indices, not energies: every k-point is included,
    nb_bot = b_bot clamped up to nbmin (strict INT(e1) > nbmin), nb_top = b_top
    clamped down to nbmax. Returns (included=True, nb_bot, nb_top)."""
    nb_bot = b_bot if b_bot > nbmin else nbmin
    nb_top = b_top if b_top < nbmax else nbmax
    return True, nb_bot, nb_top


def band_index_window(info, spins):
    """Resolve the (b_bot, b_top) band-index window for proj_mode 1 and 2.

    proj_mode 2 (dmftproj.f:233-237): b_bot=INT(e_bot), b_top=INT(e_top) taken
    directly from the indmftpr window line.

    proj_mode 1 (dmftproj.f:704-722): e_bot/e_top are Fermi-shifted energies;
    scan every spin and k-point for bands with e_bot < E <= e_top and take the
    global min/max band index, seeded with b_bot=1000, b_top=1 so an empty scan
    keeps that seed. `spins` is the list of read_almblm dicts (one per spin)."""
    if info['proj_mode'] == 2:
        return int(info['e_bot']), int(info['e_top'])
    if info['proj_mode'] != 1:
        raise ValueError('band_index_window is only valid for proj_mode 1 or 2')
    e_bot, e_top = info['e_bot'], info['e_top']
    b_bot, b_top = 1000, 1
    for sp in spins:
        for kp in sp['kp']:
            for off, ib in enumerate(range(kp['nbmin'], kp['nbmax'] + 1)):
                e = kp['eband'][off]
                if e > e_bot and e <= e_top:
                    if ib > b_top:
                        b_top = ib
                    if ib < b_bot:
                        b_bot = ib
    return b_bot, b_top


# --- gfortran-style real formatters ------------------------------------------

def fmt(x):
    """One list-directed real, gfortran-style (leading sign space, ~17 sig)."""
    return '  %.16E' % float(x)


def write_row(f, arr):
    f.write(''.join(fmt(x) for x in arr) + '\n')


def fmt_scalar(x):
    x = float(x)
    return '%.16f' % x if abs(x) < 1e5 else '%.16E' % x


def require_identity_struct(case):
    """Fail explicitly until representative struct rotations are ported."""
    with open(case + '.struct') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if 'LOCAL ROT MATRIX:' in line:
            rows = [line.split(':', 1)[1], *lines[i + 1:i + 3]]
            rotation = np.array([[to_float(x) for x in row.split()] for row in rows])
            if not np.allclose(rotation, np.eye(3), atol=1e-7, rtol=0):
                raise NotImplementedError('non-identity struct local rotations are not implemented yet')


def validate_so_windows(spins, info, windows, below=None):
    """Require valid selected windows and the same absolute bands in both spins."""
    if len(spins) != 2:
        raise ValueError('spin-orbit projectors require both up and down almblm files')
    bw = band_index_window(info, spins) if info['proj_mode'] else None
    for spin in spins[1:]:
        if spin['nk'] != spins[0]['nk']:
            raise ValueError('up/down k-point counts differ')
        for ik, kp in enumerate(spin['kp']):
            selected = (select_window(kp['nbmin'], kp['nbmax'], kp['eband'], info['e_bot'], info['e_top'])
                        if bw is None else select_band_window(kp['nbmin'], kp['nbmax'], *bw))
            if selected != windows[ik]:
                raise ValueError('spin-orbit run requires identical up/dn band windows')
            if below is not None:
                lower = (select_window(kp['nbmin'], kp['nbmax'], kp['eband'], -1e6, info['e_bot'])
                         if bw is None else select_band_window(kp['nbmin'], kp['nbmax'], 1, bw[0] - 1))
                if lower != below[ik]:
                    raise ValueError('spin-orbit run requires identical up/dn below-window bands')
    if not all(incl and lo <= hi for incl, lo, hi in windows):
        raise ValueError('projectors require a nonempty included window at every k-point')
