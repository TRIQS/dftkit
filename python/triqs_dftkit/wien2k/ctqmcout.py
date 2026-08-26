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

ctqmcout layout is outputqmc.f:66-613. The cubic transform uses the exact
analytic harmonics; with the precision-fixed dmftproj (full double precision
templates and KIND=8 casts) the port reproduces case.ctqmcout to machine
precision on a full-rank window.
"""

import os
import numpy as np
from scipy.linalg.lapack import zheev
from scipy.linalg.blas import zgemm

from ._dmftproj import (band_index_window, dmat, fmt_scalar, read_almblm,
                        read_dmftsym, read_fromfile, read_indmftpr, reptrans,
                        rotloc_rotl_so, select_band_window, select_window,
                        tmat, write_row)


def _sqrt_inv(O):
    """O^{-1/2} of a Hermitian matrix, reproducing orthogonal.f sqrtm
    (inv=.TRUE.: Z diag(w^{-1/2}) Z^H). Use the same LAPACK routine the Fortran
    calls, ZHEEV('V','U'), via scipy so the eigenvectors match rather than the
    divide-and-conquer ZHEEVD of numpy.linalg.eigh. dmftproj takes 1/sqrt of the
    eigenvalue as a *complex* sqrt (W_comp = CMPLX(W,0)); reproduce that with a
    complex power. The result D1 @ conj(Z).T matches the Fortran's ZGEMM('N','T').

    With a full-rank overlap (enough bands for the correlated spin-orbitals) the
    eigenvectors are unique and ctqmcout matches dmftproj to machine precision.
    A rank-deficient overlap (narrow window, more orbitals than bands) makes the
    near-null eigenvectors non-unique, so O^{-1/2} amplifies the last-ULP libm
    difference; that is a numerical property of the degenerate case, not the
    port (see the ctqmcout test, which uses a full-rank window)."""
    w, Z, info = zheev(O, compute_v=1, lower=0)   # 'V', 'U'
    D1 = Z * (w.astype(complex) ** -0.5)          # Z @ diag(w^{-1/2})
    return zgemm(1.0, D1, np.conj(Z), trans_b=1)  # ZGEMM('N','T'): D1 @ conj(Z)^T


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

def _rloc_rotrep(op, l, transmat, ifSO):
    """rotloc%rotrep(l)%mat, the Rloc rotation in the new basis for a non-mixing
    shell, with identity struct local rotation (rotloc_ref Euler = 0).

    Under SO it is the 2(2l+1) spinor rotation: rotloc%rotl =
    blkdiag(ephase*D, conj(ephase)*D), rotrep = S rotl S^H, S =
    blkdiag(transmat, transmat). Without SO srot%timeinv is always false, there
    is no spin phase, and it reduces to the bare (2l+1) transmat D transmat^H
    (set_rotloc.f non-SO branch). Returns (rotrep, timeinv)."""
    a, b, g, iprop = op['a'], op['b'], op['g'], op['iprop']
    D = dmat(l, a, b, g, float(iprop))
    if not ifSO:
        return transmat @ D @ np.conj(transmat.T), False
    krotm = op['krotm']
    det2 = krotm[0, 0] * krotm[1, 1] - krotm[0, 1] * krotm[1, 0]
    timeinv = det2 < 0.0
    phase = (g - a) if timeinv else (a + g)
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


def _rloc_rotrep_mixing(l, ref, ops, iatom, iref, P):
    """rotloc(iatom)%rotrep(l)%mat for a mixing SO shell: the composed Rloc
    spinor rotation (rotloc_rotl_so) put into the new basis with the full
    2(2l+1) transmat P (set_rotloc.f mixing branch). Returns (rotrep, timeinv)."""
    rotl, timeinv = rotloc_rotl_so(l, ref, ops, iatom, iref)
    rotrep = (P @ rotl @ P.T) if timeinv else (P @ rotl @ np.conj(P.T))
    return rotrep, timeinv


# --- ctqmcout writer ---------------------------------------------------------

def write_ctqmcout(case):
    """Read <case>.almblm{up,dn}, <case>.indmftpr, <case>.struct,
    <case>.dmftsym and write <case>.ctqmcout in the dmftproj format."""
    info = read_indmftpr(case + '.indmftpr')
    nsym, ops, rotloc_ref = read_dmftsym(case + '.dmftsym', rotloc=True)

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
    orbs = _build_orbs(info)
    ncrorb = len(crorbs)
    norb = len(orbs)

    # band-index window for proj_mode 1/2 (set_projections is called with band
    # indices, not energies). proj_mode 1 scans all spins for the global window.
    bw = band_index_window(info, spins) if info['proj_mode'] != 0 else None

    # window in [e_bot, e_top]  (proj_mode 0) or [b_bot, b_top]  (mode 1/2)
    windows = []
    for ik in range(nk):
        kp = spins[0]['kp'][ik]
        if info['proj_mode'] == 0:
            windows.append(select_window(
                kp['nbmin'], kp['nbmax'], kp['eband'],
                info['e_bot'], info['e_top']))
        else:
            windows.append(select_band_window(
                kp['nbmin'], kp['nbmax'], bw[0], bw[1]))

    # window below e_bot (mode 0) or below b_bot (mode 1/2), for qbbot.
    # Mode 1/2: set_projections(1, b_bot-1) -> select_band_window with top=b_bot-1.
    win_below = []
    for ik in range(nk):
        kp = spins[0]['kp'][ik]
        if info['proj_mode'] == 0:
            win_below.append(select_window(
                kp['nbmin'], kp['nbmax'], kp['eband'], -1e6, info['e_bot']))
        else:
            win_below.append(select_band_window(
                kp['nbmin'], kp['nbmax'], 1, bw[0] - 1))

    # qbbot: point integration over bands below e_bot, is=1 only under SO
    qbbot = 0.0
    for ispin in range(ns):
        for ik in range(nk):
            incl, nb_bot, nb_top = win_below[ik]
            if incl:
                qbbot += (nb_top - nb_bot + 1) * spins[ispin]['kp'][ik]['weight']
        if ifSO:
            break

    transmats = {}
    mixing = {}
    for cr in crorbs:
        key = (cr['l'], cr['sort'])
        if cr['basis'] == 'fromfile':
            transmats[key], mixing[key] = read_fromfile(
                info['sorts'][cr['sort'] - 1]['sourcefile'], cr['l'])
        else:
            transmats[key], mixing[key] = reptrans(cr['basis'], cr['l']), False

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
    # Non-mixing: mat_rep[(icr, ik, is)] -> (2l+1, nbsel) per spin block.
    # Mixing: mat_rep[(icr, ik)] -> (2*(2l+1), nbsel), the stacked spinor block.
    mat_rep = {}
    for icr, cr in enumerate(crorbs):
        l = cr['l']
        atom = cr['atom']
        sort = cr['sort']
        d = 2 * l + 1
        transmat = transmats[(l, sort)]
        ismix = mixing[(l, sort)]
        op = rotloc_op[icr]
        rot = dmat(l, op['a'], op['b'], op['g'], float(op['iprop']))
        for ik in range(nk):
            incl, nb_bot, nb_top = windows[ik]
            if not incl:
                continue
            kp0 = spins[0]['kp'][ik]
            off_bot = nb_bot - kp0['nbmin']
            off_top = nb_top - kp0['nbmin']
            nbsel = off_top - off_bot + 1
            spin_blocks = []
            for ispin in range(ns):
                sp = spins[ispin]
                kp = sp['kp'][ik]
                nlo = sp['nLO'][(l, sort)]
                P = np.zeros((d, nbsel), dtype=complex)
                for mi, m in enumerate(range(-l, l + 1)):
                    lm = l * l + (m + l)        # 0-based packed index
                    for j, off in enumerate(range(off_bot, off_top + 1)):
                        val = kp['Alm'][lm, atom, off]
                        for ilo in range(nlo):
                            val += kp['Clm'][ilo, lm, atom, off] * \
                                sp['ovl_LO_u'][(ilo + 1, l, sort)]
                        P[mi, j] = val
                spin_blocks.append(rot @ P)     # rot_projectmat per spin
            if ismix:
                stack = np.vstack(spin_blocks)  # (2*(2l+1), nbsel), up then dn
                mat_rep[(icr, ik)] = transmat @ stack
            else:
                for ispin in range(ns):
                    mat_rep[(icr, ik, ispin)] = transmat @ spin_blocks[ispin]

    # Loewdin orthonormalization. Under SO (orthogonal_wannier_SO) the up+dn
    # blocks of each crorb are stacked together and a single overlap is
    # orthonormalized per k. Without SO (orthogonal_wannier) spin is a good
    # quantum number, so each spin is orthonormalized independently over the
    # stack of (2l+1)-wide crorb blocks of that spin.
    if ifSO:
        ortho_groups = [[(ik, None)] for ik in range(nk)]
    else:
        ortho_groups = [[(ik, ispin)] for ik in range(nk) for ispin in range(ns)]
    for group in ortho_groups:
        (ik, gspin) = group[0]
        incl, nb_bot, nb_top = windows[ik]
        if not incl:
            continue
        blocks = []
        layout = []  # (key, nrows)
        for icr, cr in enumerate(crorbs):
            l = cr['l']
            if ifSO and mixing[(l, cr['sort'])]:
                blocks.append(mat_rep[(icr, ik)])
                layout.append(((icr, ik), 2 * (2 * l + 1)))
            elif ifSO:
                for ispin in range(ns):
                    blocks.append(mat_rep[(icr, ik, ispin)])
                    layout.append(((icr, ik, ispin), 2 * l + 1))
            else:
                blocks.append(mat_rep[(icr, ik, gspin)])
                layout.append(((icr, ik, gspin), 2 * l + 1))
        D = np.vstack(blocks)                 # ndim x nbnd
        # match the Fortran's exact BLAS calls (orthogonal_wannier[_SO]): the
        # near-singular O^{-1/2} amplifies any last-bit difference, so use the
        # identical ZGEMM trans flags rather than numpy's conj-transpose copies.
        O = zgemm(1.0, D, D, trans_b=2)        # ZGEMM('N','C'): D @ D^H
        S = _sqrt_inv(O)
        D_orth = zgemm(1.0, S, D)             # ZGEMM('N','N'): O^{-1/2} @ D
        row = 0
        for key, nrows in layout:
            mat_rep[key] = D_orth[row:row + nrows, :]
            row += nrows

    # ---- Rloc rotrep per crorb ----
    rloc_blocks = []
    for cr in crorbs:
        l = cr['l']
        transmat = transmats[(l, cr['sort'])]
        iref = sum(info['mult'][:cr['sort'] - 1]) + 1
        if mixing[(l, cr['sort'])]:
            rloc_blocks.append(_rloc_rotrep_mixing(
                l, rotloc_ref[cr['sort'] - 1], ops, cr['atom'], iref, transmat))
        else:
            op = next(o for o in ops if o['perm'][iref - 1] == cr['atom'])
            rloc_blocks.append(_rloc_rotrep(op, l, transmat, ifSO))

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
                    (cr['atom'], cr['sort'], l, size, 1 if ifSO else 0, 1))

        # Rloc block per crorb. SO writes the 2*(2l+1) spinor whole shell;
        # non-SO the bare (2l+1) shell. The time-reversal flag follows under
        # SO always, under non-SO only when ifSP (outputqmc.f:347-351).
        for rotrep, timeinv in rloc_blocks:
            for m in range(rotrep.shape[0]):
                write_row(f, rotrep[m, :].real)
            for m in range(rotrep.shape[0]):
                write_row(f, rotrep[m, :].imag)
            if ifSO or ifSP:
                f.write('%6d\n' % (1 if timeinv else 0))

        # complex-harmonics -> basis transform block (crorb%first only). Mixing
        # writes the full 2(2l+1) transmat; non-mixing the spin block-diagonal.
        for cr in crorbs:
            if not cr['first']:
                continue
            l = cr['l']
            transmat = transmats[(l, cr['sort'])]
            d = 2 * l + 1
            if mixing[(l, cr['sort'])]:
                spinrot = transmat
            elif ifSO:
                spinrot = np.zeros((2 * d, 2 * d), dtype=complex)
                spinrot[:d, :d] = transmat
                spinrot[d:, d:] = transmat
            else:
                spinrot = transmat       # bare (2l+1) transform, non-SO
            dim = spinrot.shape[0]
            f.write('%6d %6d \n' % (1, dim))
            for m in range(dim):
                write_row(f, spinrot[m, :].real)
            for m in range(dim):
                write_row(f, spinrot[m, :].imag)

        # number of bands per k (skip is=2 under SO)
        for ispin in range(ns):
            if ifSP and ifSO and ispin == 1:
                continue
            for ik in range(nk):
                incl, nb_bot, nb_top = windows[ik]
                f.write('%6d\n' % abs(nb_top - nb_bot + 1))

        # projector block: DO ik, DO icrorb. Mixing writes the full 2(2l+1)
        # block from is=1; non-mixing the two (2l+1) spin blocks.
        for ik in range(nk):
            for icr, cr in enumerate(crorbs):
                l = cr['l']
                if mixing[(l, cr['sort'])]:
                    P = mat_rep[(icr, ik)]
                    for m in range(2 * (2 * l + 1)):
                        write_row(f, P[m, :].real)
                    for m in range(2 * (2 * l + 1)):
                        write_row(f, P[m, :].imag)
                    continue
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
