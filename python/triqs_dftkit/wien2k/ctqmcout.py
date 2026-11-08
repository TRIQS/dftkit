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
- irrep selection (indmftpr irrep count and flag line): without spin-orbit
  coupling, and for a spin-mixing basis, only the selected irreps of a shell
  are orthonormalized and written, one correlated-shell entry per irrep; a
  spin-diagonal basis under spin-orbit coupling keeps the whole shell.
- the Rloc spinor block rotloc%rotrep (outputqmc.f:303-313), built from the
  symmetry operation mapping the representative atom onto each equivalent one
  (set_rotloc.f / setsym.f).

ctqmcout layout is outputqmc.f:66-613. The cubic transform uses the exact
analytic harmonics; against a dmftproj with double-precision casts and
full-precision templates the port reproduces case.ctqmcout to machine
precision on a full-rank window.
"""

import os
import numpy as np
from scipy.linalg.lapack import zheev
from scipy.linalg.blas import zgemm

from ._dmftproj import (build_orbs as _build_orbs, projector_windows, dmat, fmt_scalar, read_almblm,
                        basis_dreps, block_rows, irrep_blocks,
                        read_dmftsym, read_fromfile, read_indmftpr, reptrans,
                        rotloc_orbital, rotloc_rotl_so,
                        tmat, write_row)


def _sqrt_inv(O):
    """O^{-1/2} of a Hermitian matrix, reproducing orthogonal.f sqrtm
    (inv=.TRUE.: Z diag(w^{-1/2}) Z^H). Use the same LAPACK routine the Fortran
    calls, ZHEEV('V','U'), via scipy so the eigenvectors match rather than the
    divide-and-conquer ZHEEVD of numpy.linalg.eigh. dmftproj takes 1/sqrt of the
    eigenvalue as a *complex* sqrt (W_comp = CMPLX(W,0)); reproduce that with a
    complex power. The result D1 @ conj(Z).T matches the Fortran's ZGEMM('N','T').

    Full-rank overlaps give a well-defined inverse square root. Singular
    overlaps cannot define orthonormal correlated projectors and are rejected.
    """
    w, Z, info = zheev(O, compute_v=1, lower=0)   # 'V', 'U'
    if info != 0:
        raise np.linalg.LinAlgError(f'Loewdin eigensolver failed (info={info})')
    if w[0] <= np.finfo(float).eps * len(w) * abs(w[-1]):
        raise ValueError('Loewdin overlap is rank deficient; widen the projection window')
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


def _select_irreps(crorbs, info, mixing, ifSO):
    """Attach the correlated irreps to each crorb: `blocks` (offset, dim,
    irep) of the written rows, `rows` their indices, `split` whether the
    shell is written per irep, and `dreps` the irrep dimensions of the
    basis. A shell without selection is one block over all its rows."""
    for cr in crorbs:
        l, sort = cr['l'], cr['sort']
        sortinfo = info['sorts'][sort - 1]
        ismix = mixing[(l, sort)]
        dreps = None
        if l > 0 and l in sortinfo['correps']:
            dreps = basis_dreps(cr['basis'], l, sortinfo['sourcefile'])
        blocks = irrep_blocks(sortinfo, l, ismix, ifSO, dreps)
        cr['split'] = blocks is not None
        cr['dreps'] = dreps
        if blocks is None:
            blocks = [(0, 2 * (2 * l + 1) if ismix else 2 * l + 1, 1)]
        cr['blocks'] = blocks
        cr['rows'] = block_rows(blocks)


def write_projector_blocks(f, cr, mat_rep, icr, ik, ns, ismix):
    """Write the correlated projector of one crorb at one k: per irep block,
    the real rows then the imaginary rows, each over the spins of a
    spin-diagonal basis (outputqmc.f, outband.f)."""
    for off, size, _ in cr['blocks']:
        if ismix:
            P = mat_rep[(icr, ik)]
            for m in range(off, off + size):
                write_row(f, P[m, :].real)
            for m in range(off, off + size):
                write_row(f, P[m, :].imag)
            continue
        for ispin in range(ns):
            P = mat_rep[(icr, ik, ispin)]
            for m in range(off, off + size):
                write_row(f, P[m, :].real)
        for ispin in range(ns):
            P = mat_rep[(icr, ik, ispin)]
            for m in range(off, off + size):
                write_row(f, P[m, :].imag)


# --- Rloc spinor representation (set_rotloc.f / setsym.f) ---------------------

def _rloc_rotrep(l, ref, ops, iatom, iref, transmat, ifSO):
    """rotloc(iatom)%rotrep(l)%mat for a non-mixing shell: the composed Rloc
    spinor rotation rotloc_ref then the symmetry op mapping iref onto iatom
    (set_rotloc.f / setsym.f), put into the new basis with S = blkdiag(transmat,
    transmat). Shares rotloc_rotl_so with the mixing path so the non-trivial
    struct reference rotation rotloc_ref (e.g. Te, beta=pi/2) is not dropped.
    Without SO it reduces to the bare (2l+1) D(op)D(ref) transmat sandwich
    (set_rotloc.f non-SO branch). Returns (rotrep, timeinv)."""
    if not ifSO:
        D = rotloc_orbital(l, ref, ops, iatom, iref)
        return transmat @ D @ np.conj(transmat.T), False
    rotl, timeinv = rotloc_rotl_so(l, ref, ops, iatom, iref)
    d = 2 * l + 1
    S = np.zeros((2 * d, 2 * d), dtype=complex)
    S[:d, :d] = transmat
    S[d:, d:] = transmat
    rotrep = (S @ rotl @ S.T) if timeinv else (S @ rotl @ np.conj(S.T))
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

    windows, win_below = projector_windows(info, spins)

    # qbbot: point integration over bands below e_bot, is=1 only under SO
    qbbot = 0.0
    for ispin in range(ns):
        for ik in range(nk):
            incl, nb_bot, nb_top = win_below[ispin][ik]
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
    if not ifSO and any(mixing.values()):
        raise NotImplementedError('spin-mixing bases require spin-orbit coupling')
    _select_irreps(crorbs, info, mixing, ifSO)

    # rotloc Euler angles per crorb: the symmetry op mapping the representative
    # atom of the sort onto this atom (set_rotloc.f). With identity struct
    # local rotation, rotloc%(a,b,g,iprop) are that op's Euler angles, and
    # rot_projectmat applies dmat(l, a, b, g, iprop) in the |lm> basis
    # (rot_projectmat.f:59-65). For atom 2 this is a non-trivial in-plane
    # rotation (a=270) that mixes the m-rows; omitting it flips the m=1,2 rows.
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
        ref = rotloc_ref[sort - 1]
        rot = rotloc_orbital(l, ref, ops, atom, sum(info['mult'][:sort - 1]) + 1)
        for ik in range(nk):
            spin_blocks = []
            for ispin in range(ns):
                _, nb_bot, nb_top = windows[ispin][ik]
                nbsel = nb_top - nb_bot + 1
                sp = spins[ispin]
                kp = sp['kp'][ik]
                off_bot = nb_bot - kp['nbmin']
                off_top = nb_top - kp['nbmin']
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
        layout = []  # (key, rows): the rows of mat_rep[key] that are orthonormalized
        for icr, cr in enumerate(crorbs):
            if ifSO and mixing[(cr['l'], cr['sort'])]:
                layout.append(((icr, ik), cr['rows']))
            elif ifSO:
                for ispin in range(ns):
                    layout.append(((icr, ik, ispin), cr['rows']))
            else:
                layout.append(((icr, ik, gspin), cr['rows']))
        D = np.vstack([mat_rep[key][rows] for key, rows in layout])  # ndim x nbnd
        # match the Fortran's exact BLAS calls (orthogonal_wannier[_SO]): the
        # near-singular O^{-1/2} amplifies any last-bit difference, so use the
        # identical ZGEMM trans flags rather than numpy's conj-transpose copies.
        O = zgemm(1.0, D, D, trans_b=2)        # ZGEMM('N','C'): D @ D^H
        S = _sqrt_inv(O)
        D_orth = zgemm(1.0, S, D)             # ZGEMM('N','N'): O^{-1/2} @ D
        row = 0
        for key, rows in layout:
            mat_rep[key] = mat_rep[key].copy()
            mat_rep[key][rows] = D_orth[row:row + len(rows), :]
            row += len(rows)

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
            rloc_blocks.append(_rloc_rotrep(
                l, rotloc_ref[cr['sort'] - 1], ops, cr['atom'], iref, transmat, ifSO))

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
            if cr['split']:
                sizes = [(size, irep) for _, size, irep in cr['blocks']]
            else:
                sizes = [(2 * (2 * l + 1) if ifSO else 2 * l + 1, 1)]
            for size, irep in sizes:
                f.write('%6d %6d %6d %6d %6d %6d \n' %
                        (cr['atom'], cr['sort'], l, size, 1 if ifSO else 0, irep))

        # Rloc block per crorb. SO writes the 2*(2l+1) spinor whole shell;
        # non-SO the bare (2l+1) shell; a split shell the diagonal block of
        # each correlated irrep. The time-reversal flag follows under
        # SO always, under non-SO only when ifSP (outputqmc.f:347-351).
        for cr, (rotrep, timeinv) in zip(crorbs, rloc_blocks):
            blocks = cr['blocks'] if cr['split'] else [(0, rotrep.shape[0], 1)]
            for off, size, _ in blocks:
                sub = rotrep[off:off + size, off:off + size]
                for m in range(size):
                    write_row(f, sub[m, :].real)
                for m in range(size):
                    write_row(f, sub[m, :].imag)
            if ifSO or ifSP:
                f.write('%6d\n' % (1 if timeinv else 0))

        # complex-harmonics -> basis transform block (crorb%first only). Mixing
        # writes the full 2(2l+1) transmat; non-mixing the spin block-diagonal.
        # A split shell heads it with the number and dimensions of all irreps.
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
            if cr['split']:
                dreps = cr['dreps']
                f.write(''.join('%6d ' % x for x in [len(dreps)] + dreps) + '\n')
            else:
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
                incl, nb_bot, nb_top = windows[ispin][ik]
                f.write('%6d\n' % abs(nb_top - nb_bot + 1))

        # projector block: DO ik, DO icrorb, DO correlated irrep.
        for ik in range(nk):
            for icr, cr in enumerate(crorbs):
                write_projector_blocks(f, cr, mat_rep, icr, ik, ns,
                                       mixing[(cr['l'], cr['sort'])])

        # k-weights
        for ik in range(nk):
            f.write(' %s\n' % fmt_scalar(spins[0]['kp'][ik]['weight']))

        # H(k) eigenvalues (skip is=2 under SO)
        for ispin in range(ns):
            if ifSP and ifSO and ispin == 1:
                continue
            for ik in range(nk):
                incl, nb_bot, nb_top = windows[ispin][ik]
                kp = spins[ispin]['kp'][ik]
                for ib in range(nb_bot, nb_top + 1):
                    f.write(' %s\n' % fmt_scalar(kp['eband'][ib - kp['nbmin']]))
