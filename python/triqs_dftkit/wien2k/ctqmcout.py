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

"""Pure-Python generation of the dmftproj correlated-shell projector file
(case.ctqmcout), replacing the projector path of the dmftproj Fortran
executable (outputqmc.f subroutine outqmc).

Given the alm/blm coefficient files (case.almblmup / case.almblmdn, one per
spin), the projector definition (case.indmftpr), the structure
(case.struct) and the symmetry input (case.dmftsym), this rebuilds the
correlated (Wannier) projector pr_crorb%mat_rep and writes case.ctqmcout in
the format the converter reads.

It follows the Fortran path field by field:

- almblm read order (dmftproj.f:559-689): header, per-sort per-l overlap
  block, then per-k banners + idum/nbmin/nbmax + rtetr/eband lines, then the
  per-(l,m) (Alm,Blm) / Clm complex coefficients in the canonical Wien2k lm
  packing lm = l*l + (m+l) + 1.
- band-window selection (set_projections.f:89-129, proj_mode==0).
- raw correlated projector P(m,ib) = Alm(lm) + sum_ilo Clm(ilo,lm)*ovl_LO_u
  (set_projections.f:243-267); Blm and the radial s12 do not enter pr_crorb.
- local rotation rot_projectmat (Wigner-D of rotloc, identity here) then the
  angular-basis transform mat_rep = transmat . (Rloc . P) where transmat =
  <new_i|lm> is the cubic-harmonics transform (set_projections.f:450-478).
- Loewdin orthonormalization of the stacked correlated projectors,
  P <- O^{-1/2} P with O = D D^H over all correlated rows
  (orthogonal_wannier_SO, dmftproj.f:847). This step is what the ctqmcout
  numbers carry; it is applied per k over the full ndim x nbands stack.
- the Rloc spinor block rotloc%rotrep (outputqmc.f:303-313), built from the
  symmetry operation mapping the representative atom onto each equivalent one
  (set_rotloc.f / setsym.f).

ctqmcout layout is outputqmc.f:66-613. The cubic transform coefficients are
stored in dmftproj via a single-precision CMPLX cast (set_ang_trans.f:146);
to match the committed file the transmat is cast to complex64 before use
(the dft_tools #148 noise, e.g. 0.70710676908 for 1/sqrt2).
"""

import os
import numpy as np

from ._dmftproj import (dmat, fmt_scalar, read_almblm, read_dmftsym,
                        read_indmftpr, reptrans, select_window, tmat, write_row)


def _sqrt_inv(O):
    """O^{-1/2} of a Hermitian matrix (orthogonal.f sqrtm with inv=.TRUE.:
    Z diag(w^{-1/2}) Z^H). dmftproj evaluates 1/sqrt of the eigenvalue as a
    *complex* sqrt (W_comp = CMPLX(W,0)), so a (numerically) negative
    eigenvalue contributes i/sqrt(|w|) rather than NaN; reproduce that with a
    complex power. The final result D1 @ conj(Z).T matches sqrtm exactly."""
    w, Z = np.linalg.eigh(O)
    D1 = Z * (w.astype(complex) ** -0.5)        # Z @ diag(w^{-1/2})
    return D1 @ np.conj(Z).T


# --- correlated / included orbital descriptors -------------------------------

def _build_crorbs(info):
    crorbs = []
    for isort in range(1, info['nsort'] + 1):
        sortinfo = info['sorts'][isort - 1]
        for l in sortinfo['correlated_ls']:
            for imu in range(1, info['mult'][isort - 1] + 1):
                atom = sum(info['mult'][:isort - 1]) + imu
                crorbs.append(dict(atom=atom, sort=isort, l=l,
                                  basis=sortinfo['basis'], first=(imu == 1)))
    return crorbs


def _build_orbs(info):
    orbs = []
    for isort in range(1, info['nsort'] + 1):
        sortinfo = info['sorts'][isort - 1]
        for l in sortinfo['included_ls']:
            for imu in range(1, info['mult'][isort - 1] + 1):
                atom = sum(info['mult'][:isort - 1]) + imu
                orbs.append(dict(atom=atom, sort=isort, l=l))
    return orbs


# --- Rloc spinor representation (set_rotloc.f / setsym.f) ---------------------

def _rloc_rotrep(op, l, transmat):
    """rotloc%rotrep(l)%mat, the 2(2l+1) spinor rotation in the new basis for a
    non-mixing SO shell, with identity struct local rotation (rotloc_ref
    Euler = 0). rotloc%rotl = blkdiag(ephase*D, conj(ephase)*D) with
    D = D(R[isym])_{lm}, ephase = exp(i*phase/2); rotrep = S rotl S^H,
    S = blkdiag(transmat, transmat). Returns (rotrep, timeinv)."""
    a, b, g, iprop = op['a'], op['b'], op['g'], op['iprop']
    krotm = op['krotm']
    det2 = krotm[0, 0] * krotm[1, 1] - krotm[0, 1] * krotm[1, 0]
    timeinv = det2 < 0.0
    phase = (g - a) if timeinv else (a + g)
    D = dmat(l, a, b, g, float(iprop))
    if timeinv:
        D = tmat(l) @ np.conj(D)       # setsym.f:320-326 orbital time reversal
    ephase = np.exp(1j * phase / 2)
    d = 2 * l + 1
    rotl = np.zeros((2 * d, 2 * d), dtype=complex)
    rotl[:d, :d] = ephase * D
    rotl[d:, d:] = np.conj(ephase) * D
    S = np.zeros((2 * d, 2 * d), dtype=complex)
    S[:d, :d] = transmat
    S[d:, d:] = transmat
    rotrep = S @ rotl @ np.conj(S.T)
    return rotrep, timeinv


# --- ctqmcout writer ---------------------------------------------------------

def write_ctqmcout(case):
    """Read <case>.almblm{up,dn}, <case>.indmftpr, <case>.struct,
    <case>.dmftsym and write <case>.ctqmcout in the dmftproj format."""
    info = read_indmftpr(case + '.indmftpr')
    nsym, ops = read_dmftsym(case + '.dmftsym')

    up, dn = case + '.almblmup', case + '.almblmdn'
    if os.path.exists(up) and os.path.exists(dn):
        spin_files = [up, dn]
    else:
        spin_files = [case + '.almblm']
    ifSP = len(spin_files) == 2
    ifSO = bool(info['so'])
    ns = 2 if ifSP else 1

    spins = [read_almblm(p, info, projectors=True) for p in spin_files]
    nk = spins[0]['nk']
    elecn = spins[0]['elecn']

    crorbs = _build_crorbs(info)
    if any(cr['basis'] == 'fromfile' for cr in crorbs):
        raise NotImplementedError(
            'ctqmcout for fromfile/mixing correlated bases is not yet covered '
            'by a test fixture; cubic and complex bases are supported')
    orbs = _build_orbs(info)
    ncrorb = len(crorbs)
    norb = len(orbs)

    # window in [e_bot, e_top]
    windows = []
    for ik in range(nk):
        kp = spins[0]['kp'][ik]
        windows.append(select_window(kp['nbmin'], kp['nbmax'], kp['eband'],
                                     info['e_bot'], info['e_top']))

    # window below e_bot (for qbbot)
    win_below = []
    for ik in range(nk):
        kp = spins[0]['kp'][ik]
        win_below.append(select_window(kp['nbmin'], kp['nbmax'], kp['eband'],
                                       -1e6, info['e_bot']))

    # qbbot: point integration over bands below e_bot, is=1 only under SO
    qbbot = 0.0
    for ispin in range(ns):
        for ik in range(nk):
            incl, nb_bot, nb_top = win_below[ik]
            if incl:
                qbbot += (nb_top - nb_bot + 1) * spins[ispin]['kp'][ik]['weight']
        if ifSO:
            break

    transmats = {(cr['l'], cr['sort']): reptrans(cr['basis'], cr['l'])
                 for cr in crorbs}

    # rotloc Euler angles per crorb: the symmetry op mapping the representative
    # atom of the sort onto this atom (set_rotloc.f). With identity struct
    # local rotation, rotloc%(a,b,g,iprop) are that op's Euler angles, and
    # rot_projectmat applies dmat(l, a, b, g, iprop) in the |lm> basis
    # (rot_projectmat.f:59-65). For atom 2 this is a non-trivial in-plane
    # rotation (a=270) that mixes the m-rows; omitting it flips the m=1,2 rows.
    rotloc_op = {}
    for icr, cr in enumerate(crorbs):
        iref = sum(info['mult'][:cr['sort'] - 1]) + 1
        rotloc_op[icr] = next(o for o in ops
                              if o['perm'][iref - 1] == cr['atom'])

    # ---- raw correlated projector mat_rep, then Loewdin orthonormalize ----
    # mat_rep[(icr, ik, is)] -> (2l+1, nbsel)
    mat_rep = {}
    for icr, cr in enumerate(crorbs):
        l = cr['l']
        atom = cr['atom']
        sort = cr['sort']
        transmat = transmats[(l, sort)]
        op = rotloc_op[icr]
        rot = dmat(l, op['a'], op['b'], op['g'], float(op['iprop']))
        for ik in range(nk):
            incl, nb_bot, nb_top = windows[ik]
            if not incl:
                continue
            kp0 = spins[0]['kp'][ik]
            off_bot = nb_bot - kp0['nbmin']
            off_top = nb_top - kp0['nbmin']
            for ispin in range(ns):
                sp = spins[ispin]
                kp = sp['kp'][ik]
                nlo = sp['nLO'][(l, sort)]
                nbsel = off_top - off_bot + 1
                P = np.zeros((2 * l + 1, nbsel), dtype=complex)
                for mi, m in enumerate(range(-l, l + 1)):
                    lm = l * l + (m + l)        # 0-based packed index
                    for j, off in enumerate(range(off_bot, off_top + 1)):
                        val = kp['Alm'][lm, atom, off]
                        for ilo in range(nlo):
                            val += kp['Clm'][ilo, lm, atom, off] * \
                                sp['ovl_LO_u'][(ilo + 1, l, sort)]
                        P[mi, j] = val
                # rot_projectmat (local rotation) then the basis transform
                P = transmat @ (rot @ P)
                mat_rep[(icr, ik, ispin)] = P

    # Loewdin: per k stack all crorb rows (is=1 block then is=2 block per crorb)
    for ik in range(nk):
        incl, nb_bot, nb_top = windows[ik]
        if not incl:
            continue
        blocks = []
        layout = []  # (icr, is, nrows)
        for icr, cr in enumerate(crorbs):
            l = cr['l']
            for ispin in range(ns):
                blocks.append(mat_rep[(icr, ik, ispin)])
                layout.append((icr, ispin, 2 * l + 1))
        D = np.vstack(blocks)                 # ndim x nbnd
        O = D @ np.conj(D.T)
        S = _sqrt_inv(O)
        D_orth = S @ D
        row = 0
        for icr, ispin, nrows in layout:
            mat_rep[(icr, ik, ispin)] = D_orth[row:row + nrows, :]
            row += nrows

    # ---- Rloc rotrep per crorb ----
    rloc_blocks = []
    for cr in crorbs:
        l = cr['l']
        transmat = transmats[(l, cr['sort'])]
        iref = sum(info['mult'][:cr['sort'] - 1]) + 1
        op = next(o for o in ops if o['perm'][iref - 1] == cr['atom'])
        rloc_blocks.append(_rloc_rotrep(op, l, transmat))

    # ---- write ----
    with open(case + '.ctqmcout', 'w') as f:
        f.write('13.605698\n')
        f.write('%6d\n' % nk)
        f.write('%6d\n' % (1 if ifSP else 0))
        f.write('%6d\n' % (1 if ifSO else 0))
        f.write(' %s\n' % fmt_scalar(qbbot))
        f.write(' %s\n' % fmt_scalar(elecn))

        f.write('%6d\n' % norb)
        for o in orbs:
            dim = 2 * (2 * o['l'] + 1) if ifSO else 2 * o['l'] + 1
            f.write('%6d %6d %6d %6d \n' % (o['atom'], o['sort'], o['l'], dim))

        f.write('%6d\n' % ncrorb)
        for cr in crorbs:
            l = cr['l']
            size = 2 * (2 * l + 1) if ifSO else 2 * l + 1
            f.write('%6d %6d %6d %6d %6d %6d \n' %
                    (cr['atom'], cr['sort'], l, size, 1, 1))

        # Rloc block per crorb (non-mixing SP+SO whole shell)
        for rotrep, timeinv in rloc_blocks:
            for m in range(rotrep.shape[0]):
                write_row(f, rotrep[m, :].real)
            for m in range(rotrep.shape[0]):
                write_row(f, rotrep[m, :].imag)
            f.write('%6d\n' % (1 if timeinv else 0))

        # complex-harmonics -> basis transform block (crorb%first only)
        for cr in crorbs:
            if not cr['first']:
                continue
            l = cr['l']
            transmat = transmats[(l, cr['sort'])]
            d = 2 * l + 1
            spinrot = np.zeros((2 * d, 2 * d), dtype=complex)
            spinrot[:d, :d] = transmat
            spinrot[d:, d:] = transmat
            f.write('%6d %6d \n' % (1, 2 * d))
            for m in range(2 * d):
                write_row(f, spinrot[m, :].real)
            for m in range(2 * d):
                write_row(f, spinrot[m, :].imag)

        # number of bands per k (skip is=2 under SO)
        for ispin in range(ns):
            if ifSP and ifSO and ispin == 1:
                continue
            for ik in range(nk):
                incl, nb_bot, nb_top = windows[ik]
                f.write('%6d\n' % abs(nb_top - nb_bot + 1))

        # projector block: DO ik, DO icrorb (non-mixing SO whole-shell)
        for ik in range(nk):
            for icr, cr in enumerate(crorbs):
                l = cr['l']
                for ispin in range(ns):
                    P = mat_rep[(icr, ik, ispin)]
                    for mi in range(2 * l + 1):
                        write_row(f, P[mi, :].real)
                for ispin in range(ns):
                    P = mat_rep[(icr, ik, ispin)]
                    for mi in range(2 * l + 1):
                        write_row(f, P[mi, :].imag)

        # k-weights
        for ik in range(nk):
            f.write(' %s\n' % fmt_scalar(spins[0]['kp'][ik]['weight']))

        # H(k) eigenvalues (skip is=2 under SO)
        for ispin in range(ns):
            if ifSP and ifSO and ispin == 1:
                continue
            for ik in range(nk):
                incl, nb_bot, nb_top = windows[ik]
                kp = spins[ispin]['kp'][ik]
                for ib in range(nb_bot, nb_top + 1):
                    f.write(' %s\n' % fmt_scalar(kp['eband'][ib - kp['nbmin']]))
