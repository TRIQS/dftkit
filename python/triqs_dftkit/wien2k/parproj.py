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

"""Pure-Python generation of the dmftproj partial-projector file
(case.parproj), replacing the parproj path of the dmftproj Fortran executable
(outputqmc.f).

This is the sibling of ctqmcout.py. Where ctqmcout covers the correlated
shells (crorb, l_inc==2) and runs them through the Loewdin orthonormalization,
parproj covers ALL included shells (orb, l_inc in {1,2}) and writes the RAW
radial-normalized Theta projector. The three structural differences from
ctqmcout are:

1. No Loewdin orthonormalization. pr_orb%matn_rep is never touched by
   orthogonal_wannier; the only normalization is the per-orbital radial
   transform s12 = O_radial^{+1/2} (orthogonal_r with inv=.FALSE.) applied as
   coeff.s12 per (m, ib) in set_projections.f:273-332.

2. The Theta projector keeps the full radial vector coeff = [Alm, Blm, Clm...]
   of length n = nLO+2; s12 is built from u_dot_norm, ovl_LO_u and the second
   LO-overlap token ovl_LO_udot (set_projections.f:163-197). All n radial
   channels are written separately (the ir loop).

3. Three extra per-orbital blocks: the BZ-symmetrized, tetrahedron-weighted,
   Rloc-rotated density matrix densmat (density.f + symmetrize_mat + rotdens_mat),
   the per-orbital Rloc spinor rotrep (set_rotloc.f), and a time-reversal flag.

The reference case CaOs2 is SP+SO (ifSP=ifSO=1, ns=2), so every per-orbital
matrix is the full 2*(2l+1)=10 spin+orbital block. The cubic transmat (cast to
complex64 for the Os shells, as dmftproj does in set_ang_trans.f:146) and the
rot_projectmat local rotation are applied to matn_rep before writing, exactly
as in ctqmcout.

outputqmc.f parproj writer: 897-1183. set_projections.f s12/projector:
163-197/273-622. density.f Theta path: 588-915. symmetrize_mat.f /
rot_dens.f / setsym.f / set_rotloc.f for the symmetrization and Rloc rotation.
"""

import math
import os
import numpy as np


# --- list-directed token stream over a Fortran free-form text file -----------

class _Reader:
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
    return float(tok.replace('D', 'E').replace('d', 'e'))


def _to_complex(s):
    s = s.strip().lstrip('(').rstrip(')')
    re_s, im_s = s.split(',')
    return complex(_to_float(re_s), _to_float(im_s))


def _read_complex(r):
    buf = list(r.record())
    while ')' not in ''.join(buf):
        buf += r.record()
    return _to_complex(''.join(buf))


def _read_two_complex(r):
    """Read a record holding two list-directed complex values (Alm, Blm)."""
    buf = list(r.record())
    while ''.join(buf).count(')') < 2:
        buf += r.record()
    s = ''.join(buf)
    cut = s.index(')') + 1
    return _to_complex(s[:cut]), _to_complex(s[cut:])


# --- angular basis: cubic d transform transmat = <new_i|lm> ------------------

_CUBIC = {
    2: np.array([
        [0, 0, 1, 0, 0],
        [2 ** -0.5, 0, 0, 0, 2 ** -0.5],
        [-(2 ** -0.5), 0, 0, 0, 2 ** -0.5],
        [0, 2 ** -0.5, 0, -(2 ** -0.5), 0],
        [0, 2 ** -0.5, 0, 2 ** -0.5, 0],
    ], dtype=complex),
}


def _reptrans(basis, l):
    """transmat = <new_i|lm>. dmftproj stores the cubic coefficients via a
    single-precision CMPLX cast (set_ang_trans.f:146); reproduce the float32
    truncation so the written numbers match bit-for-bit."""
    if basis == 'cubic' and l in _CUBIC:
        return _CUBIC[l].astype(np.complex64).astype(np.complex128)
    return np.eye(2 * l + 1, dtype=complex)


# --- Wigner D matrix (dmftproj convention, setsym.f dmat) --------------------

def _small_d(l, m, n, b):
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


def _dmat(l, a, b, c, det):
    D = np.zeros((2 * l + 1, 2 * l + 1), dtype=complex)
    for m in range(-l, l + 1):
        for n in range(-l, l + 1):
            v = np.exp(1j * n * a) * np.exp(1j * m * c) * _small_d(l, m, n, b)
            if det < -0.5:
                v *= (-1) ** l
            D[m + l, n + l] = v
    return D


def _tmat(l):
    """Complex-conjugation operator in the spherical-harmonic basis,
    T[m, -m] = (-1)^m (timeinv.f)."""
    T = np.zeros((2 * l + 1, 2 * l + 1), dtype=complex)
    for m in range(-l, l + 1):
        T[-m + l, m + l] = (-1) ** m
    return T


def _sqrtm_real_sym(O):
    """O^{+1/2} of a real symmetric matrix (orthogonal_r with inv=.FALSE.:
    Z diag(sqrt(w)) Z^T, real part kept). dmftproj evaluates sqrt of the
    eigenvalue as a complex sqrt (W_comp = CMPLX(W,0)); reproduce that with a
    complex power so a negative eigenvalue gives i*sqrt(|w|)."""
    w, Z = np.linalg.eigh(O)
    D1 = Z * np.sqrt(w.astype(complex))
    return (D1 @ Z.T).real


# --- inputs ------------------------------------------------------------------

def _read_indmftpr(indmftpr):
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
        if basis == 'fromfile':
            i += 1
        l_inc = [int(x) for x in raw[i].split()]
        i += 1
        ireps = [int(x) for x in raw[i].split()]
        i += 1
        correlated_ls = [l for l in range(len(l_inc)) if l_inc[l] == 2]
        included_ls = [l for l in range(len(l_inc)) if l_inc[l] in (1, 2)]
        if any(n > 0 for n in ireps):
            i += 1
        if correlated_ls:
            so = int(raw[i].split()[0])
            i += 1
        sorts.append(dict(basis=basis, correlated_ls=correlated_ls,
                          included_ls=included_ls))
    last = raw[-1].split()
    e_bot, e_top = _to_float(last[0]), _to_float(last[1])
    proj_mode = int(last[2]) if len(last) >= 3 else 0
    return dict(nsort=nsort, mult=mult, lmax=lmax, sorts=sorts, so=so,
                e_bot=e_bot, e_top=e_top, proj_mode=proj_mode)


def _read_dmftsym(path):
    """Symmetry operations and the per-sort representative Rloc rotation from
    case.dmftsym. The symop block carries (perm, a, b, g, iprop, krotm); the
    'Global->local' tail carries the representative rotloc krotm + Euler angles
    per sort (setsym.f:437-455)."""
    lines = open(path).read().split('\n')
    nsym = int(lines[0].split()[0])
    perms = [[int(x) for x in lines[1 + i].split()] for i in range(nsym)]
    starts = [i for i, l in enumerate(lines) if 'Sym. op.' in l]
    ops = []
    for k, s in enumerate(starts[:nsym]):
        toks = lines[s + 1].split()
        a, b, g = (math.radians(float(x)) for x in toks[:3])
        iprop = int(toks[3])
        krotm = np.array([[_to_float(x) for x in lines[s + 2 + r].split()]
                          for r in range(3)])
        ops.append(dict(perm=perms[k], a=a, b=b, g=g, iprop=iprop,
                        krotm=krotm))
    # Global->local representative rotloc per sort.
    gl = next(i for i, l in enumerate(lines) if 'Global->local' in l)
    rotloc_ref = []
    j = gl + 1
    while len(rotloc_ref) < nsym and j < len(lines):
        if lines[j].strip() == '' or not lines[j].split()[0].lstrip('-').isdigit():
            j += 1
            continue
        # sort index line, then 3 krotm rows, then the Euler/iprop line.
        if len(lines[j].split()) == 1:
            krotm = np.array([[_to_float(x) for x in lines[j + 1 + r].split()]
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

def _read_almblm(path, info):
    """Read one spin's almblm file. Keeps Alm, Blm and Clm (parproj needs the
    full radial vector), plus u_dot_norm, ovl_LO_u, ovl_LO_udot for the s12
    build and the per-band tetrahedron weights for the density matrix."""
    nsort = info['nsort']
    lmax = info['lmax']
    mult = info['mult']
    nlm = (lmax + 1) ** 2
    natom = sum(mult)
    nloat_max = 1  # max nLO over (l, sort); grown below

    r = _Reader(path)
    elecn = _to_float(r.record()[0])
    nk = int(r.record()[0])
    r.record()                         # nloat
    eferm = _to_float(r.record()[0])

    nLO = {}
    u_dot_norm = {}
    ovl_LO_u = {}
    ovl_LO_udot = {}
    kp = [None] * nk

    for isrt in range(1, nsort + 1):
        for l in range(lmax + 1):
            u_dot_norm[(l, isrt)] = _to_float(r.record()[0])
            n = int(r.record()[0])
            nLO[(l, isrt)] = n
            nloat_max = max(nloat_max, n)
            for ilo in range(1, n + 1):
                toks = r.record()
                ovl_LO_u[(ilo, l, isrt)] = _to_float(toks[0])
                ovl_LO_udot[(ilo, l, isrt)] = _to_float(toks[1])
        for ik in range(nk):
            r.skip()                   # "IK = .." banner
            r.skip()                   # 3-int line
            head = r.record()
            nbmin, nbmax = int(head[1]), int(head[2])
            nb = nbmax - nbmin + 1
            if kp[ik] is None:
                kp[ik] = dict(
                    nbmin=nbmin, nbmax=nbmax, eband=None, weight=None,
                    tetr=None,
                    Alm=np.zeros((nlm, natom + 1, nb), dtype=complex),
                    Blm=np.zeros((nlm, natom + 1, nb), dtype=complex),
                    Clm=np.zeros((4, nlm, natom + 1, nb), dtype=complex))
            eband = np.empty(nb)
            tetr = np.empty(nb)
            for off in range(nb):
                toks = r.record()
                tetr[off] = _to_float(toks[0])
                eband[off] = _to_float(toks[1])
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
                            alm, blm = _read_two_complex(r)
                            kp[ik]['Alm'][lm, iatom, off] = alm
                            kp[ik]['Blm'][lm, iatom, off] = blm
                            for ilo in range(nLO[(l, isrt)]):
                                clm = _read_complex(r)
                                kp[ik]['Clm'][ilo, lm, iatom, off] = clm
                            lm += 1

    return dict(elecn=elecn, eferm=eferm, nk=nk, nlm=nlm, natom=natom,
                nLO=nLO, u_dot_norm=u_dot_norm, ovl_LO_u=ovl_LO_u,
                ovl_LO_udot=ovl_LO_udot, kp=kp)


# --- band-window selection (set_projections.f, proj_mode==0) -----------------

def _select_window(nbmin, nbmax, eband, e1, e2):
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


# --- included orbital descriptors --------------------------------------------

def _build_orbs(info):
    orbs = []
    for isort in range(1, info['nsort'] + 1):
        sortinfo = info['sorts'][isort - 1]
        for l in sortinfo['included_ls']:
            for imu in range(1, info['mult'][isort - 1] + 1):
                atom = sum(info['mult'][:isort - 1]) + imu
                orbs.append(dict(atom=atom, sort=isort, l=l,
                                 basis=sortinfo['basis']))
    return orbs


# --- s12 radial transform (set_projections.f:163-197) ------------------------

def _build_s12(l, isrt, n, spin):
    """O_radial^{+1/2} for one (l, sort, spin): the n x n radial overlap built
    from u_dot_norm + ovl_LO_u + ovl_LO_udot, then its matrix square root."""
    s = np.zeros((n, n))
    s[0, 0] = 1.0
    s[1, 1] = spin['u_dot_norm'][(l, isrt)]
    for ilo in range(1, spin['nLO'][(l, isrt)] + 1):
        s[1 + ilo, 1 + ilo] = 1.0
        s[1 + ilo, 0] = spin['ovl_LO_u'][(ilo, l, isrt)]
        s[0, 1 + ilo] = spin['ovl_LO_u'][(ilo, l, isrt)]
        s[1 + ilo, 1] = spin['ovl_LO_udot'][(ilo, l, isrt)]
        s[1, 1 + ilo] = spin['ovl_LO_udot'][(ilo, l, isrt)]
    return _sqrtm_real_sym(s)


# --- srot rotrep (setsym.f) --------------------------------------------------

def _srot_rotrep_nonmixing(op, l, transmat, ifSP, ifSO):
    """srot(isym)%rotrep(l,isrt)%mat, the (2l+1) up/up block of the symmetry op
    in the new basis: transmat . D(R)_{lm} . transmat^H, with the orbital
    time-reversal operator applied for the magnetic (timeinv) operations.

    Returns (rotrep[2l+1,2l+1], timeinv, phase)."""
    a, b, g, iprop = op['a'], op['b'], op['g'], op['iprop']
    krotm = op['krotm']
    det2 = krotm[0, 0] * krotm[1, 1] - krotm[0, 1] * krotm[1, 0]
    if ifSP and ifSO:
        timeinv = det2 < 0.0
        phase = (g - a) if timeinv else (a + g)
    else:
        timeinv = False
        phase = 0.0
    rotl = _dmat(l, a, b, g, float(iprop))
    rotrep = transmat @ rotl @ np.conj(transmat.T)
    if timeinv:
        # timeinv_op in the new basis: reptrans T reptrans^T applied to conj.
        tmat = _tmat(l)
        tinv = transmat @ tmat @ transmat.T
        rotrep = tinv @ np.conj(rotrep)
    return rotrep, timeinv, phase


# --- rotloc rotrep (set_rotloc.f) under SP+SO, non-mixing --------------------

def _euler(Rot):
    """Euler angles (a, b, g) of a 3x3 proper rotation, dmftproj convention
    (set_rotloc.f euler, called on TRANSPOSE(krotm))."""
    def dot(u, v):
        return float(np.dot(u, v))

    def vecprod(u, v):
        return np.cross(u, v)
    R = Rot
    y = np.array([0.0, 1.0, 0.0])
    z = np.array([0.0, 0.0, 1.0])
    yyy = R[:, 1]
    zz = R[:, 2]
    yy = vecprod(z, zz)
    y_norm = math.sqrt(dot(yy, yy))
    pi = math.pi
    if y_norm < 1e-10:
        d = dot(y, yyy)
        a = math.acos(d / abs(d)) if abs(d) > 1.0 else math.acos(d)
        if dot(z, zz) > 0.0:
            c = 0.0
            b = 0.0
            if yyy[0] > 0.0:
                a = 2 * pi - a
        else:
            c = a
            a = 0.0
            b = pi
            if yyy[0] < 0.0:
                c = 2 * pi - c
    else:
        yy = yy / y_norm
        aa, bb, cc = dot(y, yy), dot(z, zz), dot(yy, yyy)
        if abs(aa) > 1.0:
            aa = aa / abs(aa)
        if abs(bb) > 1.0:
            bb = bb / abs(bb)
        if abs(cc) > 1.0:
            cc = cc / abs(cc)
        b = math.acos(bb)
        a = math.acos(aa)
        c = math.acos(cc)
        if yy[0] > 0.0:
            a = 2 * pi - a
        pom = vecprod(yy, yyy)
        if dot(pom, zz) < 0.0:
            c = 2 * pi - c
    return a, b, c


def _rotloc_rotrep_so(orb, ops, info, ref, transmat):
    """rotloc(iatom)%rotrep(l)%mat, the full 2*(2l+1) Rloc spinor rotation in
    the new basis under SP+SO, non-mixing (set_rotloc.f), plus timeinv flag.

    ref is the representative-sort rotloc (krotm, a, b, g, iprop). set_rotloc
    runs the same composition for EVERY atom of the sort, including the
    representative: it finds the FIRST symmetry op R[isym] with
    perm(iref)==iatom and composes srot%rotl with the representative spinor
    rotloc. For the representative this op need not be the identity."""
    l = orb['l']
    isrt = orb['sort']
    iatom = orb['atom']
    iref = sum(info['mult'][:isrt - 1]) + 1
    d = 2 * l + 1

    # rotloc%rotl(2d) initial value (setsym.f:496-528): spmt (x) D(rotloc_ref).
    Dref = _dmat(l, ref['a'], ref['b'], ref['g'], float(ref['iprop']))
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

    # Compose with the first symmetry op mapping iref onto iatom (set_rotloc.f:60).
    isym = next(k for k, o in enumerate(ops)
                if o['perm'][iref - 1] == iatom)
    op = ops[isym]
    det2 = op['krotm'][0, 0] * op['krotm'][1, 1] - \
        op['krotm'][0, 1] * op['krotm'][1, 0]
    timeinv = det2 < 0.0
    srot_phase = (op['g'] - op['a']) if timeinv else (op['a'] + op['g'])
    # srot%rotl D(R[isym])_{lm}, with the timeinv operator applied if magnetic
    # (setsym.f:318-331 applies timeinv_op to srot%rotl before set_rotloc).
    rotl_sym = _dmat(l, op['a'], op['b'], op['g'], float(op['iprop']))
    if timeinv:
        rotl_sym = _tmat(l) @ np.conj(rotl_sym)
    ephase = np.exp(1j * srot_phase / 2.0)
    tmp = np.zeros((2 * d, 2 * d), dtype=complex)
    tmp[:d, :d] = ephase * rotl_sym
    tmp[d:, d:] = np.conj(ephase) * rotl_sym
    if timeinv:
        rotl = tmp @ np.conj(rotl)
    else:
        rotl = tmp @ rotl

    # rotrep = S rotl S^(H or T), S = blkdiag(transmat, transmat).
    S = np.zeros((2 * d, 2 * d), dtype=complex)
    S[:d, :d] = transmat
    S[d:, d:] = transmat
    if timeinv:
        rotrep = S @ rotl @ S.T
    else:
        rotrep = S @ rotl @ np.conj(S.T)
    return rotrep, timeinv


# --- projector matn_rep (set_projections.f Theta path) -----------------------

def _build_matn_rep(orb, info, spins, ns, windows, transmat, rot):
    """matn_rep[(ik, is)] -> (2l+1, nbsel, n): the raw radial-normalized Theta
    projector with the local rotation and the basis transform applied per
    radial channel (non-mixing SP+SO path, set_projections.f:591-621)."""
    l = orb['l']
    atom = orb['atom']
    sort = orb['sort']
    out = {}
    for ik in range(info['nk']):
        incl, nb_bot, nb_top = windows[ik]
        if not incl:
            continue
        kp0 = spins[0]['kp'][ik]
        off_bot = nb_bot - kp0['nbmin']
        off_top = nb_top - kp0['nbmin']
        nbsel = off_top - off_bot + 1
        for ispin in range(ns):
            sp = spins[ispin]
            kp = sp['kp'][ik]
            n = sp['nLO'][(l, sort)] + 2
            s12 = _build_s12(l, sort, n, sp)
            matn = np.zeros((2 * l + 1, nbsel, n), dtype=complex)
            for mi, m in enumerate(range(-l, l + 1)):
                lm = l * l + (m + l)
                for j, off in enumerate(range(off_bot, off_top + 1)):
                    coeff = np.zeros(n, dtype=complex)
                    coeff[0] = kp['Alm'][lm, atom, off]
                    coeff[1] = kp['Blm'][lm, atom, off]
                    for ilo in range(n - 2):
                        coeff[2 + ilo] = kp['Clm'][ilo, lm, atom, off]
                    matn[mi, j, :] = coeff @ s12
            # rot_projectmat then the basis transform, per radial channel.
            for ir in range(n):
                matn[:, :, ir] = transmat @ (rot @ matn[:, :, ir])
            out[(ik, ispin)] = matn
    return out


# --- density matrix (density.f Theta path, non-mixing SP+SO) -----------------

def _orbital_densmat_blocks(orb, info, spins, ns, windows, matn_reps):
    """The nsp=4 raw density-matrix blocks (up/up, dn/dn, up/dn, dn/up) for one
    orbital, point-integrated with the geometric k-weight (density.f:889-907,
    tetr=.FALSE. path). The window is the bands below e_bot (dmftproj.f:798)."""
    l = orb['l']
    d = 2 * l + 1
    blocks = [np.zeros((d, d), dtype=complex) for _ in range(4)]
    for ik in range(info['nk']):
        incl, nb_bot, nb_top = windows[ik]
        if not incl:
            continue
        kp0 = spins[0]['kp'][ik]
        off_bot = nb_bot - kp0['nbmin']
        off_top = nb_top - kp0['nbmin']
        n = spins[0]['nLO'][(l, orb['sort'])] + 2
        for iss in range(1, 5):
            if iss <= 2:
                is_, is1 = iss, iss
            else:
                is_ = iss - 2
                is1 = 3 - is_
            isp, isp1 = is_ - 1, is1 - 1
            weight = spins[isp]['kp'][ik]['weight']
            acc = np.zeros((d, d), dtype=complex)
            for i in range(n):
                mat = matn_reps[orb['atom']][(ik, isp)][:, :, i]
                cmat = matn_reps[orb['atom']][(ik, isp1)][:, :, i]
                acc += mat @ np.conj(cmat.T)
            blocks[iss - 1] += acc * weight
    return blocks


def _symmetrize_densmat(orbs, info, ops, nsym, srot_rotreps, dens_raw, ns,
                        ifSP, ifSO):
    """symmetrize_mat for the non-mixing SP+SO path (symmetrize_mat.f:196-276).
    dens_raw[atom] -> [4 blocks]. Returns symmetrized blocks per atom.

    The blocks are summed over symmetry ops with srot%rotrep, the spin block
    phase ephase (up/dn, dn/up), and the orbital-time-reversal conjugation for
    magnetic ops; equivalent atoms of a sort scatter into one another via the
    permutation. Result divided by nsym."""
    out = {}
    # group orbitals by sort, in orb order (they are already sort-major).
    isort_groups = {}
    for o in orbs:
        isort_groups.setdefault(o['sort'], []).append(o)
    for isrt, group in isort_groups.items():
        l = group[0]['l']
        d = 2 * l + 1
        mult = len(group)
        sym = [[np.zeros((d, d), dtype=complex) for _ in range(4)]
               for _ in range(mult)]
        for imult in range(mult):
            iatom = group[imult]['atom']
            for isym in range(nsym):
                op = ops[isym]
                rotrep, timeinv, phase = srot_rotreps[(isrt, isym)]
                jorb = op['perm'][iatom - 1] - iatom + imult  # 0-based target
                for iss in range(1, 5):
                    is_ = iss if iss <= 2 else iss - 2
                    isp = is_ - 1
                    tmp = dens_raw[iatom][iss - 1].copy()
                    if ifSP and timeinv:
                        tmp = np.conj(tmp)
                    ephase = 1.0
                    if iss == 3:
                        ephase = np.exp(1j * phase)
                    elif iss == 4:
                        ephase = np.exp(-1j * phase)
                    val = rotrep @ (tmp @ np.conj(rotrep.T)) * ephase
                    sym[jorb][iss - 1] += val
        for imult in range(mult):
            iatom = group[imult]['atom']
            out[iatom] = [sym[imult][k] / nsym for k in range(4)]
    return out


def _rotdens_densmat(orb, blocks, rotrep_loc, timeinv):
    """rotdens_mat for non-mixing SP+SO (rot_dens.f:151-193): assemble the four
    blocks into a 2*(2l+1) matrix, apply inverse(Rloc) D Rloc, return it as the
    full 2*(2l+1) densprint matrix."""
    l = orb['l']
    d = 2 * l + 1
    rd = np.zeros((2 * d, 2 * d), dtype=complex)
    rd[:d, :d] = blocks[0]
    rd[d:, d:] = blocks[1]
    rd[:d, d:] = blocks[2]
    rd[d:, :d] = blocks[3]
    if timeinv:
        rd = np.conj(rd @ rotrep_loc)
        rd = rotrep_loc.T @ rd
    else:
        rd = rd @ rotrep_loc
        rd = np.conj(rotrep_loc.T) @ rd
    return rd


# --- parproj writer ----------------------------------------------------------

def _fmt(x):
    return '  %.16E' % float(x)


def _w(f, arr):
    f.write(''.join(_fmt(x) for x in arr) + '\n')


def write_parproj(case):
    """Read <case>.almblm{up,dn}, <case>.indmftpr, <case>.struct,
    <case>.dmftsym and write <case>.parproj in the dmftproj format."""
    info = _read_indmftpr(case + '.indmftpr')
    nsym, ops, rotloc_ref = _read_dmftsym(case + '.dmftsym')

    up, dn = case + '.almblmup', case + '.almblmdn'
    if os.path.exists(up) and os.path.exists(dn):
        spin_files = [up, dn]
    else:
        spin_files = [case + '.almblm']
    ifSP = len(spin_files) == 2
    ifSO = bool(info['so'])
    ns = 2 if ifSP else 1

    spins = [_read_almblm(p, info) for p in spin_files]
    nk = spins[0]['nk']
    info['nk'] = nk

    orbs = _build_orbs(info)
    if any(o['basis'] == 'fromfile' for o in orbs):
        raise NotImplementedError(
            'parproj for fromfile/mixing bases is not yet covered by a test '
            'fixture; cubic and complex bases are supported')
    norb = len(orbs)

    # Two band ranges (dmftproj.f:757,836): the energy window [e_bot, e_top]
    # for the written Theta projector (block A), and the full band range
    # (-Elarge, Elarge) over which the density matrix is integrated. dmftproj
    # computes the density matrix FIRST over the full range, then overwrites
    # the projectors with the energy window before outputqmc.
    # Two band ranges (dmftproj.f:798,836): the energy window [e_bot, e_top]
    # for the written Theta projector (block A), and the bands below e_bot
    # (-Elarge, e_bot) for the density matrix. dmftproj computes the density
    # matrix over the below-e_bot range with point integration LAST (the third
    # density call is correlated-only), so densmat in outputqmc is that one.
    windows = []
    below_windows = []
    for ik in range(nk):
        kp = spins[0]['kp'][ik]
        windows.append(_select_window(kp['nbmin'], kp['nbmax'], kp['eband'],
                                      info['e_bot'], info['e_top']))
        below_windows.append(_select_window(kp['nbmin'], kp['nbmax'],
                                            kp['eband'], -1e6, info['e_bot']))

    transmats = {(o['l'], o['sort']): _reptrans(o['basis'], o['l'])
                 for o in orbs}

    # rot_projectmat local rotation: the op mapping the representative atom of
    # the sort onto this atom (set_projections.f via rot_projectmat).
    rotloc_op = {}
    for o in orbs:
        iref = sum(info['mult'][:o['sort'] - 1]) + 1
        rotloc_op[o['atom']] = next(op for op in ops
                                    if op['perm'][iref - 1] == o['atom'])

    # ---- projectors matn_rep per orbital (energy window for block A, full
    # band range for the density matrix) ----
    matn_reps = {}
    matn_reps_full = {}
    for o in orbs:
        l = o['l']
        transmat = transmats[(l, o['sort'])]
        op = rotloc_op[o['atom']]
        rot = _dmat(l, op['a'], op['b'], op['g'], float(op['iprop']))
        matn_reps[o['atom']] = _build_matn_rep(
            o, info, spins, ns, windows, transmat, rot)
        matn_reps_full[o['atom']] = _build_matn_rep(
            o, info, spins, ns, below_windows, transmat, rot)

    # ---- srot rotrep per (sort, isym) for symmetrize_mat ----
    srot_rotreps = {}
    for isrt in range(1, info['nsort'] + 1):
        if not info['sorts'][isrt - 1]['included_ls']:
            continue
        l = info['sorts'][isrt - 1]['included_ls'][0]
        transmat = transmats[(l, isrt)]
        for isym in range(nsym):
            srot_rotreps[(isrt, isym)] = _srot_rotrep_nonmixing(
                ops[isym], l, transmat, ifSP, ifSO)

    # ---- density matrices: raw -> symmetrize -> rotdens ----
    dens_raw = {}
    for o in orbs:
        dens_raw[o['atom']] = _orbital_densmat_blocks(
            o, info, spins, ns, below_windows, matn_reps_full)
    dens_sym = _symmetrize_densmat(orbs, info, ops, nsym, srot_rotreps,
                                   dens_raw, ns, ifSP, ifSO)

    # ---- rotloc rotrep per orbital ----
    rotloc_rotrep = {}
    for o in orbs:
        transmat = transmats[(o['l'], o['sort'])]
        ref = rotloc_ref[o['sort'] - 1]
        rotloc_rotrep[o['atom']] = _rotloc_rotrep_so(
            o, ops, info, ref, transmat)

    densprint = {}
    for o in orbs:
        rotrep_loc, timeinv = rotloc_rotrep[o['atom']]
        densprint[o['atom']] = _rotdens_densmat(
            o, dens_sym[o['atom']], rotrep_loc, timeinv)

    # ---- write ----
    with open(case + '.parproj', 'w') as f:
        for o in orbs:
            n = spins[0]['nLO'][(o['l'], o['sort'])] + 2
            f.write('%6d\n' % n)

        for o in orbs:
            l = o['l']
            atom = o['atom']
            n = spins[0]['nLO'][(l, o['sort'])] + 2

            # (A) Theta projector, non-mixing SP+SO (outputqmc.f:955-973).
            for ik in range(nk):
                incl, nb_bot, nb_top = windows[ik]
                for ir in range(n):
                    for ispin in range(ns):
                        P = matn_reps[atom][(ik, ispin)][:, :, ir]
                        for mi in range(2 * l + 1):
                            _w(f, P[mi, :].real)
                    for ispin in range(ns):
                        P = matn_reps[atom][(ik, ispin)][:, :, ir]
                        for mi in range(2 * l + 1):
                            _w(f, P[mi, :].imag)

            # (B) density matrix, non-mixing SP+SO 2*(2l+1) (outputqmc.f:1033-1049).
            dp = densprint[atom]
            for m in range(2 * (2 * l + 1)):
                _w(f, dp[m, :].real)
            for m in range(2 * (2 * l + 1)):
                _w(f, dp[m, :].imag)

            # (C) Rloc rotrep, non-mixing SP+SO (outputqmc.f:1148-1158).
            rotrep_loc, timeinv = rotloc_rotrep[atom]
            for m in range(2 * (2 * l + 1)):
                _w(f, rotrep_loc[m, :].real)
            for m in range(2 * (2 * l + 1)):
                _w(f, rotrep_loc[m, :].imag)
            if ifSP:
                f.write('%6d\n' % (1 if timeinv else 0))
