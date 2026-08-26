################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
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

"""Pure-Python generation of the dmftproj included-shell symmetry file
(case.sympar), replacing the corresponding output of the dmftproj Fortran
executable (outputqmc.f, unit ousympar).

case.sympar is the sibling of case.symqmc. Both write the same per-operation
spinor symmetry matrix srot(isym)%rotrep(l,isrt)%mat (setsym.f), built
identically for every shell as transmat . D(R)_lm . transmat^dag. They differ
only in the shell list: symqmc covers the CORRELATED shells (l_inc==2), sympar
covers ALL INCLUDED shells (l_inc in {1,2}), iterating orb(1:norb) instead of
crorb(1:ncrorb). The orb ordering is (sort, l, atom) per dmftproj.f:357-387.

The per-shell basis transmat is per sort: a complex shell uses the exact
identity, a cubic shell uses the single-precision CMPLX cast of the cubic-d
harmonics (set_ang_trans.f:146), giving 0.70710676908 for 1/sqrt2 — the same
truncation ctqmcout carries. Mixed within one file (Ca complex exact, Os cubic
truncated for CaOs2).

sympar drops three things symqmc/ctqmcout carry: no orbital-description block
at the top, no ifsplit/irep sub-block selection (always the full matrix), and
no crorb%ifSOat/correp metadata. It is purely group-theoretic: header, perms,
optional timeflag line (ifSP), the representation matrices, and — only when
.not.ifSP — a paramagnetic time-reversal operator per orbital.

This reference case (CaOs2) is SP+SO: indmftpr SO flag=1 sets ifSO, and
ifSO=>ifSP, so it exercises the non-mixing 2*(2l+1)=10-wide block-diag spinor
path and writes the timeflag line; the paramagnetic tail is absent.
"""

import math
import os
import numpy as np


# --- angular bases -----------------------------------------------------------
# transpose(P) = <new|m>, m = -l..l. dmftproj stores cubic/fromfile coefficients
# via a single-precision CMPLX cast (set_ang_trans.f:146); reproduce the
# float32 truncation so the written numbers match (0.70710676908 for 1/sqrt2).
# A complex basis is the exact identity (no cast).

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


def _reptrans(basis, l):
    if basis == 'cubic' and l in _CUBIC:
        return _CUBIC[l].astype(np.complex64).astype(np.complex128)
    return np.eye(2 * l + 1, dtype=complex)


# --- Wigner D matrix (dmftproj convention, setsym.f) -------------------------

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
    """T[m, -m] = (-1)^m, complex conjugation in the |lm> basis (timeinv.f)."""
    T = np.zeros((2 * l + 1, 2 * l + 1), dtype=complex)
    for m in range(-l, l + 1):
        T[-m + l, m + l] = (-1) ** m
    return T


def _timeinv_orbital(l, mat):
    return _tmat(l) @ np.conj(mat)


# --- inputs ------------------------------------------------------------------

def _read_indmftpr(indmftpr):
    """Return one entry per INCLUDED shell (l_inc in {1,2}), one per atom of its
    sort, in orb order (sort, l, atom). Each shell carries l, basis name and the
    transform P. Also returns the SO flag."""
    raw = [l.split('!')[0].strip() for l in open(indmftpr)]
    raw = [l for l in raw if l != '']
    nsort = int(raw[0].split()[0])
    mult = [int(x) for x in raw[1].split()][:nsort]
    i = 3
    so = 0
    shells = []
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
            i += 1                          # skip the correps line
        if correlated_ls:
            so = int(raw[i].split()[0])     # SO flag follows a correlated sort
            i += 1
        for l in included_ls:
            P = _reptrans(basis, l)
            for imu in range(1, mult[isort] + 1):
                atom = sum(mult[:isort]) + imu
                shells.append(dict(l=l, sort=isort + 1, atom=atom,
                                   basis=basis, P=P))
    return shells, so


def _read_dmftsym(path):
    lines = open(path).read().split('\n')
    nsym = int(lines[0].split()[0])
    perms = [[int(x) for x in lines[1 + i].split()] for i in range(nsym)]
    starts = [i for i, l in enumerate(lines) if 'Sym. op.' in l]
    ops = []
    for k, s in enumerate(starts[:nsym]):
        toks = lines[s + 1].split()
        a, b, c = (math.radians(float(x)) for x in toks[:3])
        krotm = np.array([[float(x.replace('D', 'E').replace('d', 'e'))
                           for x in lines[s + 2 + r].split()] for r in range(3)])
        ops.append(dict(perm=perms[k], a=a, b=b, c=c, krotm=krotm))
    return nsym, ops


# --- matrix construction (shared with symqmc) --------------------------------

def _phase(op, ti):
    a, c = op['a'], op['c']
    return (c - a) if ti else (a + c)       # (g-a) on magnetic ops, else (a+g)


def _l0_matrix(op, ti):
    e = np.exp(1j * _phase(op, ti) / 2)
    return np.array([[e, 0], [0, np.conj(e)]], dtype=complex)


def _nonmixing_matrix(op, shell, ti):
    """Non-mixing SP+SO whole shell: the up/up block scaled by +-(a+g)/2,
    block-diagonal over spin, with the orbital time-reversal operator on the
    magnetic operations (setsym.f, outputqmc.f:1292-1339)."""
    l, P = shell['l'], shell['P']
    rotl = _dmat(l, op['a'], op['b'], op['c'], np.linalg.det(op['krotm']))
    if ti:
        rotl = _timeinv_orbital(l, rotl)
    rotrep = P @ rotl @ np.conj(P.T)
    e = np.exp(1j * _phase(op, ti) / 2)
    d = 2 * l + 1
    mat = np.zeros((2 * d, 2 * d), dtype=complex)
    mat[:d, :d] = e * rotrep
    mat[d:, d:] = np.conj(e) * rotrep
    return mat


def _shell_matrix(op, shell, ti):
    l = shell['l']
    if l == 0:
        return _l0_matrix(op, ti)
    return _nonmixing_matrix(op, shell, ti)


# --- output ------------------------------------------------------------------

def _fmt(x):
    """One gfortran list-directed real (leading sign space, ~17 sig figs)."""
    return '  %.16E' % float(x)


def _w(f, arr):
    f.write(''.join(_fmt(x) for x in arr) + '\n')


def _write_matrix(f, mat):
    for m in range(mat.shape[0]):
        _w(f, mat[m, :].real)
    for m in range(mat.shape[0]):
        _w(f, mat[m, :].imag)


def write_sympar(case):
    """Write <case>.sympar from <case>.dmftsym, <case>.indmftpr, <case>.struct.

    Detects ifSP/ifSO from the presence of <case>.almblm{up,dn} and the indmftpr
    SO flag, exactly as oubwin/symqmc/ctqmcout. For SP+SO the file carries the
    timeflag line and 2*(2l+1)-wide block-diag spinor matrices; the paramagnetic
    time-reversal tail is written only when .not.ifSP."""
    shells, so = _read_indmftpr(case + '.indmftpr')
    if any(sh['basis'] == 'fromfile' for sh in shells):
        raise NotImplementedError(
            'sympar for fromfile/mixing bases is not yet covered by a test '
            'fixture; cubic and complex bases are supported')
    nsym, ops = _read_dmftsym(case + '.dmftsym')
    natom = len(ops[0]['perm'])

    ifSP = os.path.exists(case + '.almblmup') and os.path.exists(case + '.almblmdn')
    ifSO = bool(so)
    ifSP = ifSP or ifSO                      # ifSO => ifSP

    timeflag = []
    for op in ops:
        det2 = (op['krotm'][0, 0] * op['krotm'][1, 1]
                - op['krotm'][0, 1] * op['krotm'][1, 0])
        timeflag.append(1 if (ifSO and det2 < 0.0) else 0)

    with open(case + '.sympar', 'w') as f:
        f.write('%6d %6d\n' % (nsym, natom))
        for op in ops:
            f.write(''.join('%6d ' % p for p in op['perm']) + '\n')
        if ifSP:
            f.write(''.join('%6d ' % t for t in timeflag) + '\n')

        for isym, op in enumerate(ops):
            for sh in shells:
                l = sh['l']
                ti = bool(timeflag[isym])
                if l == 0 and not (ifSP and ifSO):
                    f.write(_fmt(1.0) + '\n')
                    f.write(_fmt(0.0) + '\n')
                    continue
                _write_matrix(f, _shell_matrix(op, sh, ti))

        if not ifSP:
            for sh in shells:
                l = sh['l']
                if l == 0:
                    f.write(_fmt(1.0) + '\n')
                    f.write(_fmt(0.0) + '\n')
                    continue
                tm = _tmat(l)
                op = sh['P'] @ tm @ sh['P'].T
                ident = np.eye(2 * l + 1, dtype=complex)
                time_op = op @ np.conj(ident)
                _write_matrix(f, time_op)
