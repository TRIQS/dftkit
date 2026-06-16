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

import os
import numpy as np

from ._dmftproj import (dmat, fmt, mixing_rotrep, read_dmftsym, read_fromfile,
                        read_indmftpr, reptrans, timeinv_orbital, tmat,
                        write_row)


# --- included shells ---------------------------------------------------------

def _included_shells(info):
    """One entry per INCLUDED shell (l_inc in {1,2}), one per atom of its sort,
    in orb order (sort, l, atom). Each shell carries l, basis name, a mixing
    flag and the transform P: a spin-coupling fromfile basis gives the full
    2(2l+1) P (P spinrot P^dag matrix), otherwise the (2l+1) up/up block (the
    single-precision-cast cubic harmonics dmftproj writes)."""
    shells = []
    for isort in range(info['nsort']):
        s = info['sorts'][isort]
        for l in s['included_ls']:
            if s['basis'] == 'fromfile':
                P, mixing = read_fromfile(s['sourcefile'], l)
            else:
                P, mixing = reptrans(s['basis'], l), False
            for imu in range(1, info['mult'][isort] + 1):
                atom = sum(info['mult'][:isort]) + imu
                shells.append(dict(l=l, sort=isort + 1, atom=atom,
                                   basis=s['basis'], P=P, mixing=mixing))
    return shells


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
    rotl = dmat(l, op['a'], op['b'], op['c'], np.linalg.det(op['krotm']))
    if ti:
        rotl = timeinv_orbital(l, rotl)
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
    if shell['mixing']:
        return mixing_rotrep(op, l, shell['P'], bool(ti))
    return _nonmixing_matrix(op, shell, ti)


# --- output ------------------------------------------------------------------

def _write_matrix(f, mat):
    for m in range(mat.shape[0]):
        write_row(f, mat[m, :].real)
    for m in range(mat.shape[0]):
        write_row(f, mat[m, :].imag)


def write_sympar(case):
    """Write <case>.sympar from <case>.dmftsym, <case>.indmftpr, <case>.struct.

    Detects ifSP/ifSO from the presence of <case>.almblm{up,dn} and the indmftpr
    SO flag, exactly as oubwin/symqmc/ctqmcout. For SP+SO the file carries the
    timeflag line and 2*(2l+1)-wide block-diag spinor matrices; the paramagnetic
    time-reversal tail is written only when .not.ifSP."""
    info = read_indmftpr(case + '.indmftpr')
    shells, so = _included_shells(info), info['so']
    nsym, ops = read_dmftsym(case + '.dmftsym')
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
                    f.write(fmt(1.0) + '\n')
                    f.write(fmt(0.0) + '\n')
                    continue
                _write_matrix(f, _shell_matrix(op, sh, ti))

        if not ifSP:
            for sh in shells:
                l = sh['l']
                if l == 0:
                    f.write(fmt(1.0) + '\n')
                    f.write(fmt(0.0) + '\n')
                    continue
                tm = tmat(l)
                op = sh['P'] @ tm @ sh['P'].T
                ident = np.eye(2 * l + 1, dtype=complex)
                time_op = op @ np.conj(ident)
                _write_matrix(f, time_op)
