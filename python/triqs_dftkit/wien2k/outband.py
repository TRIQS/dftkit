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

"""Pure-Python generation of the dmftproj band-structure projector file
(case.outband), replacing the -band path of the dmftproj Fortran executable
(outband.f).

case.outband feeds the converter method convert_bands_input: it is the
band-structure analog of case.ctqmcout. The correlated-shell projectors are
built exactly as in ctqmcout.py (raw Alm/Clm projector -> rot_projectmat local
rotation -> cubic/fromfile transmat -> Loewdin orthonormalization
orthogonal_wannier[_SO]), but evaluated on a band k-path and written in the
outband layout. case.outband also carries the raw radial-normalized Theta
projector pr_orb%matn_rep for every included shell (the parproj.py path, no
Loewdin) and the k-path labels.

Two band-mode differences from the ctqmcout path (dmftproj.f:557-581):

- nkband, the number of k-points along the plotted k-path, is read from
  <case>.klist_band (read_k_list.f). It is written in the header and drives the
  k-label block; the per-k projector/band data loop over the almblm nk, which
  must match nkband.
- the Fermi energy is read from the LAST line of <case>.indmftpr, not from the
  almblm header; the band almblm keeps the eferm record as a placeholder that
  the reader skips (READ(iualmblm,*) with no target).

outband layout is outband.f:56-278: nkband, the per-spin per-k band counts, the
Loewdin correlated projectors (DO ik, DO icrorb), H(k) eigenvalues, the
included-orbital radial sizes, the raw Theta projectors, and the k-labels.
"""

import os
import numpy as np
from scipy.linalg.blas import zgemm

from ._dmftproj import (band_index_window, select_band_window, validate_spin_windows,
                        dmat, read_almblm, read_dmftsym, read_fromfile,
                        read_indmftpr, reptrans, rotloc_orbital, select_window, to_float,
                        write_row)
from .ctqmcout import _build_crorbs, _select_irreps, _sqrt_inv, write_projector_blocks
from .parproj import (_build_matn_rep, _build_matn_rep_mixing, _build_orbs)


def _read_indmftpr_band(indmftpr):
    """case.indmftpr in band mode: the standard parse, plus the Fermi energy
    appended as the last physical line after the (e_bot e_top proj_mode) window
    line. The trailing scalar would be mistaken for the window line, so feed
    read_indmftpr the lines without that record and return eferm separately."""
    with open(indmftpr) as fh:
        raw = [l.split('!')[0].rstrip('\n') for l in fh]
    nonblank = [i for i, l in enumerate(raw) if l.strip() != '']
    eferm = to_float(raw[nonblank[-1]].split()[0])
    info = read_indmftpr(indmftpr, lines=raw[:nonblank[-1]])
    return info, eferm


def _read_almblm_band(path, info, eferm):
    """Band-mode almblm read: identical to read_almblm but the header eferm
    record is a placeholder the band path ignores; the Fermi shift uses the
    eferm passed in (from case.indmftpr, dmftproj.f:576-581). read_almblm reads
    the 4th record as eferm, so overwrite that record with the indmftpr value."""
    with open(path) as fh:
        text = fh.read().splitlines()
    text[3] = ' %.16E' % eferm          # Reader record 4 (physical line index 3)
    return read_almblm(path, info, projectors=True, lines=text)


def read_k_list(path):
    """nkband and the k-path labels from <case>.klist_band (read_k_list.f).

    Every physical line up to the END marker is a k-point; a line whose first
    character is not a blank carries a label, whose position is the k-point
    index and whose name is the leading non-blank token. Returns
    (nkband, [(pos, name), ...])."""
    lines = open(path).read().splitlines()
    nkband = 0
    labels = []
    for line in lines:
        if line[:3] == 'END':
            break
        nkband += 1
        if line[:1] != ' ':
            labels.append((nkband, line.split()[0]))
    return nkband, labels


def write_outband(case):
    """Read <case>.almblm{up,dn}, <case>.indmftpr, <case>.struct,
    <case>.dmftsym and <case>.klist_band and write <case>.outband in the
    dmftproj band format."""
    info, eferm = _read_indmftpr_band(case + '.indmftpr')
    nsym, ops, rotloc_ref = read_dmftsym(case + '.dmftsym', rotloc=True)

    up, dn = case + '.almblmup', case + '.almblmdn'
    if os.path.exists(up) and os.path.exists(dn):
        spin_files = [up, dn]
    else:
        spin_files = [case + '.almblm']
    ifSP = len(spin_files) == 2
    ifSO = bool(info['so'])
    ns = 2 if ifSP else 1

    spins = [_read_almblm_band(p, info, eferm) for p in spin_files]
    nk = spins[0]['nk']
    info['nk'] = nk

    nkband, labels = read_k_list(case + '.klist_band')
    if nkband != nk:
        raise ValueError('klist_band and almblm k-point counts differ')

    crorbs = _build_crorbs(info)
    orbs = _build_orbs(info)
    ncrorb = len(crorbs)
    norb = len(orbs)

    if any(sp['nk'] != nk for sp in spins):
        raise ValueError('up/down k-point counts differ')
    bw = band_index_window(info, spins) if info['proj_mode'] else None
    windows = [[select_window(kp['nbmin'], kp['nbmax'], kp['eband'],
                              info['e_bot'], info['e_top'])
                if bw is None else select_band_window(kp['nbmin'], kp['nbmax'], *bw)
                for kp in sp['kp']] for sp in spins]
    if ifSO:
        validate_spin_windows(spins, info, windows[0])
    elif not all(incl and lo <= hi for spin_windows in windows
                 for incl, lo, hi in spin_windows):
        raise ValueError('projectors require a nonempty included window at every k-point')

    transmats = {}
    mixing = {}
    for o in orbs:
        key = (o['l'], o['sort'])
        if key in transmats:
            continue
        if o['basis'] == 'fromfile':
            transmats[key], mixing[key] = read_fromfile(
                info['sorts'][o['sort'] - 1]['sourcefile'], o['l'])
        else:
            transmats[key], mixing[key] = reptrans(o['basis'], o['l']), False
    if not ifSO and any(mixing.values()):
        raise NotImplementedError('spin-mixing bases require spin-orbit coupling')
    _select_irreps(crorbs, info, mixing, ifSO, ops, transmats)

    # rot_projectmat local rotation: the op mapping the sort representative onto
    # this atom (set_projections.f via rot_projectmat). Same as ctqmcout/parproj.
    # ---- correlated projector mat_rep, then Loewdin orthonormalize ----
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
                    lm = l * l + (m + l)
                    for j, off in enumerate(range(off_bot, off_top + 1)):
                        val = kp['Alm'][lm, atom, off]
                        for ilo in range(nlo):
                            val += kp['Clm'][ilo, lm, atom, off] * \
                                sp['ovl_LO_u'][(ilo + 1, l, sort)]
                        P[mi, j] = val
                spin_blocks.append(rot @ P)
            if ismix:
                stack = np.vstack(spin_blocks)
                mat_rep[(icr, ik)] = transmat @ stack
            else:
                for ispin in range(ns):
                    mat_rep[(icr, ik, ispin)] = transmat @ spin_blocks[ispin]

    # Without SO the two spins are independent Hilbert spaces.
    groups = [(ik, None) for ik in range(nk)] if ifSO else [
        (ik, ispin) for ik in range(nk) for ispin in range(ns)]
    for ik, gspin in groups:
        layout = []  # (key, rows): the rows of mat_rep[key] that are orthonormalized
        for icr, cr in enumerate(crorbs):
            if ifSO and mixing[(cr['l'], cr['sort'])]:
                layout.append(((icr, ik), cr['rows']))
            elif ifSO:
                for ispin in range(ns):
                    layout.append(((icr, ik, ispin), cr['rows']))
            else:
                layout.append(((icr, ik, gspin), cr['rows']))
        D = np.vstack([mat_rep[key][rows] for key, rows in layout])
        # match the Fortran's exact BLAS calls (orthogonal_wannier_SO): the
        # near-singular O^{-1/2} amplifies any last-bit difference, so use the
        # identical ZGEMM trans flags rather than numpy's conj-transpose copies.
        O = zgemm(1.0, D, D, trans_b=2)       # ZGEMM('N','C'): D @ D^H
        S = _sqrt_inv(O)
        D_orth = zgemm(1.0, S, D)             # ZGEMM('N','N'): O^{-1/2} @ D
        row = 0
        for key, rows in layout:
            mat_rep[key] = mat_rep[key].copy()
            mat_rep[key][rows] = D_orth[row:row + len(rows), :]
            row += len(rows)

    # ---- raw Theta projectors matn_rep per included orbital (parproj path) ----
    matn_reps = {}
    for o in orbs:
        l = o['l']
        transmat = transmats[(l, o['sort'])]
        ref = rotloc_ref[o['sort'] - 1]
        iref = sum(info['mult'][:o['sort'] - 1]) + 1
        rot = rotloc_orbital(l, ref, ops, o['atom'], iref)
        if mixing[(l, o['sort'])]:
            matn_reps[(o['atom'], o['l'])] = _build_matn_rep_mixing(
                o, info, spins, windows[0], transmat, rot)
        else:
            blocks = {}
            for ispin, sp in enumerate(spins):
                blocks.update({(ik, ispin): block for (ik, _), block in
                               _build_matn_rep(o, info, [sp], windows[ispin],
                                               transmat, rot).items()})
            matn_reps[(o['atom'], o['l'])] = blocks

    # ---- write ----
    with open(case + '.outband', 'w') as f:
        f.write('%6d\n' % nkband)

        # number of bands per k (skip is=2 under SP+SO)
        for ispin in range(ns):
            if ifSP and ifSO and ispin == 1:
                continue
            for ik in range(nk):
                incl, nb_bot, nb_top = windows[ispin][ik]
                f.write('%6d\n' % abs(nb_top - nb_bot + 1))

        # correlated projector block: DO ik, DO icrorb (Loewdin orthonormalized).
        for ik in range(nk):
            for icr, cr in enumerate(crorbs):
                write_projector_blocks(f, cr, mat_rep, icr, ik, ns,
                                       mixing[(cr['l'], cr['sort'])])

        # H(k) eigenvalues (skip is=2 under SO)
        for ispin in range(ns):
            if ifSO and ispin == 1:
                continue
            for ik in range(nk):
                incl, nb_bot, nb_top = windows[ispin][ik]
                kp = spins[ispin]['kp'][ik]
                for ib in range(nb_bot, nb_top + 1):
                    f.write(' %.16E\n' % kp['eband'][ib - kp['nbmin']])

        # included-orbital radial sizes norm_radf%n
        for o in orbs:
            n = spins[0]['nLO'][(o['l'], o['sort'])] + 2
            f.write('%6d\n' % n)

        # Theta projector block: DO iorb, DO ik, DO ir. Mixing writes the full
        # 2(2l+1) block from is=1; non-mixing the two (2l+1) spin blocks.
        for o in orbs:
            l = o['l']
            atom = o['atom']
            n = spins[0]['nLO'][(l, o['sort'])] + 2
            ismix = mixing[(l, o['sort'])]
            for ik in range(nk):
                for ir in range(n):
                    if ismix:
                        P = matn_reps[(atom, l)][ik][:, :, ir]
                        for m in range(2 * (2 * l + 1)):
                            write_row(f, P[m, :].real)
                        for m in range(2 * (2 * l + 1)):
                            write_row(f, P[m, :].imag)
                        continue
                    for ispin in range(ns):
                        P = matn_reps[(atom, l)][(ik, ispin)][:, :, ir]
                        for mi in range(2 * l + 1):
                            write_row(f, P[mi, :].real)
                    for ispin in range(ns):
                        P = matn_reps[(atom, l)][(ik, ispin)][:, :, ir]
                        for mi in range(2 * l + 1):
                            write_row(f, P[mi, :].imag)

        # k-labels
        for i, (pos, name) in enumerate(labels, start=1):
            f.write('%6d%6d%s\n' % (i, pos, name))
