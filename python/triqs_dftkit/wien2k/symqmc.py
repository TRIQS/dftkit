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

import os

import numpy as np

from ._dmftproj import (dmat, mixing_rotrep, read_dmftsym, read_fromfile,
                        read_indmftpr, reptrans, timeinv_orbital, tmat)


# --- correlated shells -------------------------------------------------------

def _correlated_shells(info):
    """One entry per correlated atom (l_inc==2), in sort order. Each shell
    carries l, basis name, its transform matrix P = <new|m> and a mixing flag
    (True for a spin-coupling fromfile basis). symqmc builds its symmetry
    matrices from the EXACT cubic harmonics (no single-precision cast)."""
    shells = []
    for isort in range(info['nsort']):
        s = info['sorts'][isort]
        basis = s['basis']
        for l in s['correlated_ls']:
            if basis == 'fromfile':
                P, mixing = read_fromfile(s['sourcefile'], l)
            else:
                P, mixing = reptrans(basis, l, cast=False), False
            for _ in range(info['mult'][isort]):
                shells.append(dict(l=l, basis=basis, P=P, mixing=mixing))
    return shells


# The symqmc test reaches into these two names directly; keep them as the
# module's parsing entry points.
_read_dmftsym = read_dmftsym


def _read_correlated_shells(indmftpr, struct):
    info = read_indmftpr(indmftpr)
    return _correlated_shells(info), info['so']


def write_symqmc(case):
    """Write <case>.symqmc from <case>.dmftsym, <case>.indmftpr, <case>.struct.

    ifSO comes from the indmftpr SO flag; ifSP from the presence of
    <case>.almblm{up,dn} (ifSO => ifSP). With SO the matrices are the 2*(2l+1)
    block-diag spinor representation and a timeflag line is written; without SO
    they are the bare (2l+1) representation, with a timeflag line (all zero)
    only when ifSP. The paramagnetic time-reversal tail follows when .not.ifSP."""
    nsym, ops = read_dmftsym(case + '.dmftsym')
    info = read_indmftpr(case + '.indmftpr')
    shells, so = _correlated_shells(info), info['so']
    natom = len(ops[0]['perm'])

    ifSO = bool(so)
    ifSP = (os.path.exists(case + '.almblmup')
            and os.path.exists(case + '.almblmdn')) or ifSO

    timeinv = []
    for op in ops:
        det2 = op['krotm'][0, 0] * op['krotm'][1, 1] - op['krotm'][0, 1] * op['krotm'][1, 0]
        timeinv.append(1 if (ifSO and det2 < 0.0) else 0)

    with open(case + '.symqmc', 'w') as f:
        f.write('%6d %6d\n' % (nsym, natom))
        for op in ops:
            f.write(''.join('%6d ' % p for p in op['perm']) + '\n')
        if ifSP:
            f.write(''.join('%6d ' % t for t in timeinv) + '\n')
        for isym, op in enumerate(ops):
            for sh in shells:
                f.write(_format_matrix(_shell_matrix(op, sh, timeinv[isym], ifSO)))
        if not ifSP:
            for sh in shells:
                f.write(_format_matrix(_time_op_matrix(sh)))


def _time_op_matrix(shell):
    """Paramagnetic time-reversal operator for one correlated shell (.not.ifSP,
    outputqmc.f:835-887): the (2l+1) operator P T P^T in the new basis; the s
    shell reduces to the 1x1 identity."""
    l = shell['l']
    if l == 0:
        return np.array([[1.0 + 0j]], dtype=complex)
    P = shell['P']
    return (P @ tmat(l) @ P.T) @ np.eye(2 * l + 1, dtype=complex)


def _shell_matrix(op, shell, ti, ifSO):
    """Symmetry matrix for one correlated shell under one operation: the
    2*(2l+1) block-diag spinor representation under SO, the bare (2l+1)
    representation otherwise."""
    l = shell['l']
    if l == 0:
        return _l0_matrix(op, ti) if ifSO else np.array([[1.0 + 0j]], dtype=complex)
    if shell['mixing']:
        return _mixing_matrix(op, shell, ti)
    if not ifSO:
        return _nonmixing_orbital(op, shell)
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
    rotl = dmat(l, op['a'], op['b'], op['c'], float(op['iprop']))
    if ti:
        rotl = timeinv_orbital(l, rotl)
        rotrep = P @ rotl @ P.T          # antiunitary op: transpose, not dagger
    else:
        rotrep = P @ rotl @ np.conj(P.T)
    e = np.exp(1j * _phase(op, ti) / 2)
    d = 2 * l + 1
    mat = np.zeros((2 * d, 2 * d), dtype=complex)
    mat[:d, :d] = e * rotrep
    mat[d:, d:] = np.conj(e) * rotrep
    return mat


def _nonmixing_orbital(op, shell):
    """Non-SO spin-diagonal basis: the bare (2l+1) representation
    P D(R)_{lm} P^H (setsym.f non-SO branch). Under non-SO srot%timeinv is
    always false, so no orbital time-reversal or spin phase enters."""
    l, P = shell['l'], shell['P']
    rotl = dmat(l, op['a'], op['b'], op['c'], float(op['iprop']))
    return P @ rotl @ np.conj(P.T)


def _mixing_matrix(op, shell, ti):
    """Spin-coupling basis: full 2(2l+1) spinor representation (setsym.f,
    timeinv.f), shared with sympar via _dmftproj.mixing_rotrep."""
    return mixing_rotrep(op, shell['l'], shell['P'], bool(ti))


def _format_matrix(mat):
    out = []
    for part in (mat.real, mat.imag):
        for row in part:
            out.append(''.join('  %.14E' % x for x in row) + '\n')
    return ''.join(out)
