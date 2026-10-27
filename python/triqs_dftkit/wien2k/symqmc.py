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

All matrices are built in double precision from the exact cubic harmonics.
A case.symqmc written by a dmftproj that casts the basis or the spin phases
through single-precision CMPLX differs from this output at about 1e-7.
"""

import os

from ._dmftproj import (paramag_time_op, read_dmftsym, read_fromfile, read_indmftpr,
                        reptrans, shell_matrix)

_shell_matrix = shell_matrix


# --- correlated shells -------------------------------------------------------

def _correlated_shells(info):
    """One entry per correlated atom (l_inc==2), in sort order. Each shell
    carries l, basis name, its transform matrix P = <new|m> and a mixing flag
    (True for a spin-coupling fromfile basis). symqmc builds its symmetry
    matrices from the exact cubic harmonics."""
    shells = []
    for isort in range(info['nsort']):
        s = info['sorts'][isort]
        basis = s['basis']
        for l in s['correlated_ls']:
            if basis == 'fromfile':
                P, mixing = read_fromfile(s['sourcefile'], l)
            else:
                P, mixing = reptrans(basis, l), False
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
                f.write(_format_matrix(paramag_time_op(sh)))


def _format_matrix(mat):
    out = []
    for part in (mat.real, mat.imag):
        for row in part:
            out.append(''.join('  %.14E' % x for x in row) + '\n')
    return ''.join(out)
