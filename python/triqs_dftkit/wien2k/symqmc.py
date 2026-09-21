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

"""Pure-Python generation of the dmftproj correlated-shell symmetry file
(case.symqmc), replacing the corresponding output of the dmftproj Fortran
executable.

Given the dmftproj symmetry input (case.dmftsym), the projector definition
(case.indmftpr) and the structure (case.struct), this builds the spinor symmetry
matrices and writes them in the case.symqmc format the converter reads. It
reproduces the Fortran construction: the Wigner D matrix (dmat), the orbital
time-reversal operator for the magnetic (SP+SO) operations, the basis transform
to the chosen angular harmonics, and the spin-1/2 phase blocks.

Covers all dmftproj spin-orbit cases:

- non-mixing spin-diagonal bases (complex, cubic, and a fromfile basis whose
  spin-up and spin-down blocks coincide): the spin-reduced up/up block scaled
  by the +-(a+g)/2 phase, with the orbital time-reversal operator on the
  magnetic operations;
- mixing bases (a fromfile basis that couples spin, e.g. the |j, m_j> basis):
  the full 2(2l+1) spinor representation P spinrot P^dag, with the spinor
  time-reversal operator -i sigma_y (x) T applied to the magnetic operations;
- l = 0 (s) shells: the 2x2 spin phase block.

The mixing fromfile path is the one dft_tools #148 singles out. dmftproj reads
that basis with a single-precision CMPLX cast and a 250-column line cap (Fortran
set_ang_trans.f), so its case.symqmc carries a ~1e-7 error there; this generator
is full double precision.
"""

import math
import os
import numpy as np

# --- angular bases (transpose(P) = <new|m>, m = -l..l) -----------------------

_COMPLEX = {l: np.eye(2 * l + 1, dtype=complex) for l in range(4)}

# standard cubic harmonics, Wien2k convention
_CUBIC = {
    1: np.array([[0, 1, 0], [-1j, 0, -1j], [1, 0, -1]], dtype=complex) / 1.0,
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
        return _CUBIC[l]
    return _COMPLEX[l]


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
            f3 = 1.0 if (2 * l - m - n - 2 * t) == 0 else math.sin(b / 2) ** (2 * l - m - n - 2 * t)
            f4 = 1.0 if (2 * t + n + m) == 0 else math.cos(b / 2) ** (2 * t + n + m)
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


def _timeinv_orbital(l, mat):
    return _tmat(l) @ np.conj(mat)


# --- fromfile angular basis ---------------------------------------------------

def _read_fromfile(path, l):
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


# --- input parsing -----------------------------------------------------------

def _read_dmftsym(path):
    lines = open(path).read().split('\n')
    nsym = int(lines[0].split()[0])
    perms = [[int(x) for x in lines[1 + i].split()] for i in range(nsym)]
    rest = lines[1 + nsym:]
    starts = [i for i, l in enumerate(rest) if 'Sym. op.' in l]
    ops = []
    for k, s in enumerate(starts):
        ang = rest[s + 1].split()
        a, b, c = (math.radians(float(x)) for x in ang[:3])
        krotm = np.array([[float(x) for x in rest[s + 2 + r].split()] for r in range(3)])
        ops.append(dict(perm=perms[k], a=a, b=b, c=c, krotm=krotm))
    return nsym, ops


def _read_correlated_shells(indmftpr, struct):
    """Return the list of correlated shells, one entry per correlated atom, plus
    the SO flag, from case.indmftpr and case.struct multiplicities. Each shell
    carries l, the basis name, its transform matrix P = <new|m> and a mixing
    flag (True for a spin-coupling fromfile basis)."""
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
        sourcefile = None
        if basis == 'fromfile':
            sourcefile = os.path.join(os.path.dirname(indmftpr), raw[i])
            i += 1
        l_inc = [int(x) for x in raw[i].split()]
        i += 1
        ireps = [int(x) for x in raw[i].split()]
        i += 1
        correlated_ls = [l for l in range(len(l_inc)) if l_inc[l] == 2]
        if any(n > 0 for n in ireps):
            i += 1                       # skip the correps line
        if correlated_ls:
            so = int(raw[i].split()[0])  # SO flag follows a correlated sort
            i += 1
            for l in correlated_ls:
                if basis == 'fromfile':
                    P, mixing = _read_fromfile(sourcefile, l)
                else:
                    P, mixing = _reptrans(basis, l), False
                for _ in range(mult[isort]):
                    shells.append(dict(l=l, basis=basis, P=P, mixing=mixing))
    return shells, so


def write_symqmc(case):
    """Write <case>.symqmc from <case>.dmftsym, <case>.indmftpr, <case>.struct."""
    nsym, ops = _read_dmftsym(case + '.dmftsym')
    shells, so = _read_correlated_shells(case + '.indmftpr', case + '.struct')
    natom = len(ops[0]['perm'])

    timeinv = []
    for op in ops:
        det2 = op['krotm'][0, 0] * op['krotm'][1, 1] - op['krotm'][0, 1] * op['krotm'][1, 0]
        timeinv.append(1 if (so and det2 < 0.0) else 0)

    with open(case + '.symqmc', 'w') as f:
        f.write('%6d %6d\n' % (nsym, natom))
        for op in ops:
            f.write(''.join('%6d ' % p for p in op['perm']) + '\n')
        if so:
            f.write(''.join('%6d ' % t for t in timeinv) + '\n')
        for isym, op in enumerate(ops):
            for sh in shells:
                f.write(_format_matrix(_shell_matrix(op, sh, timeinv[isym])))


def _shell_matrix(op, shell, ti):
    """Spinor symmetry matrix for one correlated shell under one operation."""
    l = shell['l']
    if l == 0:
        return _l0_matrix(op, ti)
    if shell['mixing']:
        return _mixing_matrix(op, shell, ti)
    return _nonmixing_matrix(op, shell, ti)


def _phase(op, ti):
    a, c = op['a'], op['c']
    return (c - a) if ti else (a + c)   # (g-a) on magnetic ops, else (a+g)


def _l0_matrix(op, ti):
    """s shell: 2x2 spin phase block diag(e, conj(e)) (outputqmc.f l==0)."""
    e = np.exp(1j * _phase(op, ti) / 2)
    return np.array([[e, 0], [0, np.conj(e)]], dtype=complex)


def _nonmixing_matrix(op, shell, ti):
    """Spin-diagonal basis: the up/up block scaled by +-(a+g)/2, with the
    orbital time-reversal operator on the magnetic operations."""
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


def _mixing_matrix(op, shell, ti):
    """Spin-coupling basis: full 2(2l+1) spinor representation. The orbital
    rotation enters block-diagonal (beta=0) or block-antidiagonal (beta=pi, the
    magnetic operations); the spinor time-reversal -i sigma_y (x) T is then
    applied to the magnetic operations (setsym.f, timeinv.f)."""
    l, P = shell['l'], shell['P']
    rotl = _dmat(l, op['a'], op['b'], op['c'], np.linalg.det(op['krotm']))
    e = np.exp(1j * _phase(op, ti) / 2)
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
        tm = _tmat(l)
        tinv = np.zeros((2 * d, 2 * d), dtype=complex)
        tinv[:d, d:] = -tm
        tinv[d:, :d] = tm
        rotrep = (P @ tinv @ P.T) @ np.conj(rotrep)
    return rotrep


def _format_matrix(mat):
    out = []
    for part in (mat.real, mat.imag):
        for row in part:
            out.append(''.join('  %.14E' % x for x in row) + '\n')
    return ''.join(out)
