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

from ._dmftproj import (dmat, mixing_rotrep, read_almblm, read_dmftsym,
                        read_fromfile, read_indmftpr, reptrans,
                        rotloc_rotl_so, select_window, tmat, write_row)


def _op_timeinv(op, ifSO):
    det2 = (op['krotm'][0, 0] * op['krotm'][1, 1]
            - op['krotm'][0, 1] * op['krotm'][1, 0])
    return bool(ifSO and det2 < 0.0)


def _sqrtm_real_sym(O):
    """O^{+1/2} of a real symmetric matrix (orthogonal_r with inv=.FALSE.:
    Z diag(sqrt(w)) Z^T, real part kept). dmftproj evaluates sqrt of the
    eigenvalue as a complex sqrt (W_comp = CMPLX(W,0)); reproduce that with a
    complex power so a negative eigenvalue gives i*sqrt(|w|)."""
    w, Z = np.linalg.eigh(O)
    D1 = Z * np.sqrt(w.astype(complex))
    return (D1 @ Z.T).real


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
    rotl = dmat(l, a, b, g, float(iprop))
    rotrep = transmat @ rotl @ np.conj(transmat.T)
    if timeinv:
        # timeinv_op in the new basis: reptrans T reptrans^T applied to conj.
        tinv = transmat @ tmat(l) @ transmat.T
        rotrep = tinv @ np.conj(rotrep)
    return rotrep, timeinv, phase


# --- rotloc rotrep (set_rotloc.f) under SP+SO, non-mixing --------------------

def _rotloc_rotrep_so(orb, ops, info, ref, transmat, mixing):
    """rotloc(iatom)%rotrep(l)%mat, the full 2*(2l+1) Rloc spinor rotation in
    the new basis under SP+SO (set_rotloc.f), plus timeinv flag.

    ref is the representative-sort rotloc; rotloc_rotl_so composes the
    representative spinor rotloc with the first symmetry op mapping the sort
    representative onto this atom. The new-basis transform is the full transmat
    for a mixing (spin-coupling) basis, and blkdiag(transmat, transmat) for a
    spin-diagonal one."""
    l = orb['l']
    iref = sum(info['mult'][:orb['sort'] - 1]) + 1
    d = 2 * l + 1
    rotl, timeinv = rotloc_rotl_so(l, ref, ops, orb['atom'], iref)
    if mixing:
        S = transmat
    else:
        S = np.zeros((2 * d, 2 * d), dtype=complex)
        S[:d, :d] = transmat
        S[d:, d:] = transmat
    rotrep = (S @ rotl @ S.T) if timeinv else (S @ rotl @ np.conj(S.T))
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


def _raw_matn(orb, sp, ik, nb_bot, nb_top, rot):
    """The rot_projectmat'd (2l+1, nbsel, n) Theta projector in the |lm> basis
    for one spin, before the angular transform (set_projections.f matn_rep)."""
    l, atom, sort = orb['l'], orb['atom'], orb['sort']
    kp = sp['kp'][ik]
    off_bot = nb_bot - kp['nbmin']
    off_top = nb_top - kp['nbmin']
    nbsel = off_top - off_bot + 1
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
    for ir in range(n):
        matn[:, :, ir] = rot @ matn[:, :, ir]
    return matn


def _build_matn_rep_mixing(orb, info, spins, windows, transmat, rot):
    """matn_rep[ik] -> (2*(2l+1), nbsel, n): the mixing Theta projector. Each
    spin's rot_projectmat'd |lm> block is stacked (up then dn) and multiplied by
    the full 2(2l+1) transmat per radial channel (set_projections.f:515-585)."""
    l = orb['l']
    d = 2 * l + 1
    out = {}
    for ik in range(info['nk']):
        incl, nb_bot, nb_top = windows[ik]
        if not incl:
            continue
        up = _raw_matn(orb, spins[0], ik, nb_bot, nb_top, rot)
        dn = _raw_matn(orb, spins[1], ik, nb_bot, nb_top, rot)
        nbsel, n = up.shape[1], up.shape[2]
        matn = np.zeros((2 * d, nbsel, n), dtype=complex)
        matn[:d] = up
        matn[d:] = dn
        for ir in range(n):
            matn[:, :, ir] = transmat @ matn[:, :, ir]
        out[ik] = matn
    return out


def _orbital_densmat_mixing(orb, info, spins, windows, matn_rep):
    """The single 2*(2l+1) raw density block for a mixing orbital, point-
    integrated with the geometric k-weight (density.f:786-812): D = sum over k
    and radial channels of matn_rep matn_rep^H, already in the new basis."""
    l = orb['l']
    d = 2 * l + 1
    n = spins[0]['nLO'][(l, orb['sort'])] + 2
    dens = np.zeros((2 * d, 2 * d), dtype=complex)
    for ik in range(info['nk']):
        incl, _, _ = windows[ik]
        if not incl:
            continue
        weight = spins[0]['kp'][ik]['weight']
        for i in range(n):
            mat = matn_rep[orb['atom']][ik][:, :, i]
            dens += (mat @ np.conj(mat.T)) * weight
    return dens


def _symmetrize_densmat_mixing(orb, ops, nsym, rotreps, dens_raw):
    """symmetrize_mat for a mixing orbital (symmetrize_mat.f:140-188). For the
    representative atom: sum over symmetry ops of rotrep (conj(D) if magnetic)
    rotrep^H, divided by nsym. rotreps[isym] is srot%rotrep (mixing_rotrep)."""
    d = 2 * (2 * orb['l'] + 1)
    sym = np.zeros((d, d), dtype=complex)
    for isym in range(nsym):
        rotrep, timeinv = rotreps[isym]
        tmp = np.conj(dens_raw) if timeinv else dens_raw
        sym += rotrep @ (tmp @ np.conj(rotrep.T))
    return sym / nsym


def _rotdens_densmat_mixing(blk, rotrep_loc, timeinv):
    """rotdens_mat for a mixing orbital (rot_dens.f:119-141): inverse(Rloc) D
    Rloc on the single 2*(2l+1) block."""
    if timeinv:
        return rotrep_loc.T @ np.conj(blk @ rotrep_loc)
    return np.conj(rotrep_loc.T) @ (blk @ rotrep_loc)


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

def write_parproj(case):
    """Read <case>.almblm{up,dn}, <case>.indmftpr, <case>.struct,
    <case>.dmftsym and write <case>.parproj in the dmftproj format."""
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
    info['nk'] = nk

    orbs = _build_orbs(info)
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
        windows.append(select_window(kp['nbmin'], kp['nbmax'], kp['eband'],
                                      info['e_bot'], info['e_top']))
        below_windows.append(select_window(kp['nbmin'], kp['nbmax'],
                                            kp['eband'], -1e6, info['e_bot']))

    transmats = {}
    mixing = {}
    for o in orbs:
        key = (o['l'], o['sort'])
        if o['basis'] == 'fromfile':
            transmats[key], mixing[key] = read_fromfile(
                info['sorts'][o['sort'] - 1]['sourcefile'], o['l'])
        else:
            transmats[key], mixing[key] = reptrans(o['basis'], o['l']), False

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
        rot = dmat(l, op['a'], op['b'], op['g'], float(op['iprop']))
        if mixing[(l, o['sort'])]:
            matn_reps[o['atom']] = _build_matn_rep_mixing(
                o, info, spins, windows, transmat, rot)
            matn_reps_full[o['atom']] = _build_matn_rep_mixing(
                o, info, spins, below_windows, transmat, rot)
        else:
            matn_reps[o['atom']] = _build_matn_rep(
                o, info, spins, ns, windows, transmat, rot)
            matn_reps_full[o['atom']] = _build_matn_rep(
                o, info, spins, ns, below_windows, transmat, rot)

    # ---- rotloc rotrep per orbital ----
    rotloc_rotrep = {}
    for o in orbs:
        transmat = transmats[(o['l'], o['sort'])]
        ref = rotloc_ref[o['sort'] - 1]
        rotloc_rotrep[o['atom']] = _rotloc_rotrep_so(
            o, ops, info, ref, transmat, mixing[(o['l'], o['sort'])])

    # ---- density matrices: raw -> symmetrize -> rotdens ----
    # Non-mixing sorts use the 4-block path; mixing sorts the single-block path.
    nonmix_orbs = [o for o in orbs if not mixing[(o['l'], o['sort'])]]
    densprint = {}

    if nonmix_orbs:
        srot_rotreps = {}
        for isrt in range(1, info['nsort'] + 1):
            ls = info['sorts'][isrt - 1]['included_ls']
            if not ls or mixing[(ls[0], isrt)]:
                continue
            l = ls[0]
            for isym in range(nsym):
                srot_rotreps[(isrt, isym)] = _srot_rotrep_nonmixing(
                    ops[isym], l, transmats[(l, isrt)], ifSP, ifSO)
        dens_raw = {o['atom']: _orbital_densmat_blocks(
            o, info, spins, ns, below_windows, matn_reps_full)
            for o in nonmix_orbs}
        dens_sym = _symmetrize_densmat(nonmix_orbs, info, ops, nsym,
                                       srot_rotreps, dens_raw, ns, ifSP, ifSO)
        for o in nonmix_orbs:
            rotrep_loc, timeinv = rotloc_rotrep[o['atom']]
            densprint[o['atom']] = _rotdens_densmat(
                o, dens_sym[o['atom']], rotrep_loc, timeinv)

    for o in orbs:
        l, sort = o['l'], o['sort']
        if not mixing[(l, sort)]:
            continue
        rotreps = []
        for isym in range(nsym):
            ti = _op_timeinv(ops[isym], ifSO)
            rotreps.append(
                (mixing_rotrep(ops[isym], l, transmats[(l, sort)], ti), ti))
        raw = _orbital_densmat_mixing(o, info, spins, below_windows,
                                      matn_reps_full)
        sym = _symmetrize_densmat_mixing(o, ops, nsym, rotreps, raw)
        rotrep_loc, timeinv = rotloc_rotrep[o['atom']]
        densprint[o['atom']] = _rotdens_densmat_mixing(sym, rotrep_loc, timeinv)

    # ---- write ----
    with open(case + '.parproj', 'w') as f:
        for o in orbs:
            n = spins[0]['nLO'][(o['l'], o['sort'])] + 2
            f.write('%6d\n' % n)

        for o in orbs:
            l = o['l']
            atom = o['atom']
            n = spins[0]['nLO'][(l, o['sort'])] + 2
            ismix = mixing[(l, o['sort'])]

            # (A) Theta projector (outputqmc.f:935-973). Mixing writes the full
            # 2(2l+1) block from is=1 only; non-mixing the two (2l+1) spin blocks.
            for ik in range(nk):
                incl, nb_bot, nb_top = windows[ik]
                for ir in range(n):
                    if ismix:
                        P = matn_reps[atom][ik][:, :, ir]
                        for m in range(2 * (2 * l + 1)):
                            write_row(f, P[m, :].real)
                        for m in range(2 * (2 * l + 1)):
                            write_row(f, P[m, :].imag)
                        continue
                    for ispin in range(ns):
                        P = matn_reps[atom][(ik, ispin)][:, :, ir]
                        for mi in range(2 * l + 1):
                            write_row(f, P[mi, :].real)
                    for ispin in range(ns):
                        P = matn_reps[atom][(ik, ispin)][:, :, ir]
                        for mi in range(2 * l + 1):
                            write_row(f, P[mi, :].imag)

            # (B) density matrix, SP+SO 2*(2l+1) (outputqmc.f:1012-1066).
            dp = densprint[atom]
            for m in range(2 * (2 * l + 1)):
                write_row(f, dp[m, :].real)
            for m in range(2 * (2 * l + 1)):
                write_row(f, dp[m, :].imag)

            # (C) Rloc rotrep, SP+SO 2*(2l+1) (outputqmc.f:1130-1177).
            rotrep_loc, timeinv = rotloc_rotrep[atom]
            for m in range(2 * (2 * l + 1)):
                write_row(f, rotrep_loc[m, :].real)
            for m in range(2 * (2 * l + 1)):
                write_row(f, rotrep_loc[m, :].imag)
            if ifSP:
                f.write('%6d\n' % (1 if timeinv else 0))
