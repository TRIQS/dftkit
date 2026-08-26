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
orthogonal_wannier_SO), but evaluated on a band k-path and written in the
outband layout. case.outband also carries the raw radial-normalized Theta
projector pr_orb%matn_rep for every included shell (the parproj.py path, no
Loewdin) and the k-path labels.

Two band-mode differences from the ctqmcout path (dmftproj.f:557-581):

- nkband, the number of k-points along the plotted k-path, is read from
  <case>.klist_band (read_k_list.f). It is written in the header and drives the
  k-label block; the per-k projector/band data loop over the almblm nk, which is
  the single representative k-point lapw2 -almd writes in band mode.
- the Fermi energy is read from the LAST line of <case>.indmftpr, not from the
  almblm header; the band almblm keeps the eferm record as a placeholder that
  the reader skips (READ(iualmblm,*) with no target).

outband layout is outband.f:56-278: nkband, the per-spin per-k band counts, the
Loewdin correlated projectors (DO ik, DO icrorb), H(k) eigenvalues, the
included-orbital radial sizes, the raw Theta projectors, and the k-labels.
"""

import os
import tempfile
import numpy as np
from scipy.linalg.blas import zgemm

from ._dmftproj import (dmat, read_almblm, read_dmftsym, read_fromfile,
                        read_indmftpr, reptrans, select_window, to_float,
                        write_row)
from .ctqmcout import _build_crorbs, _sqrt_inv
from .parproj import (_build_matn_rep, _build_matn_rep_mixing, _build_orbs)


def _read_indmftpr_band(indmftpr):
    """case.indmftpr in band mode: the standard parse, plus the Fermi energy
    appended as the last physical line after the (e_bot e_top proj_mode) window
    line. The trailing scalar would be mistaken for the window line, so feed
    read_indmftpr the file without that record and return eferm separately."""
    raw = [l.split('!')[0].rstrip('\n') for l in open(indmftpr)]
    nonblank = [i for i, l in enumerate(raw) if l.strip() != '']
    eferm = to_float(raw[nonblank[-1]].split()[0])
    tmp = tempfile.NamedTemporaryFile('w', suffix='.indmftpr', delete=False,
                                      dir=os.path.dirname(os.path.abspath(indmftpr)))
    tmp.write('\n'.join(raw[:nonblank[-1]]) + '\n')
    tmp.close()
    try:
        info = read_indmftpr(tmp.name)
    finally:
        os.unlink(tmp.name)
    return info, eferm


def _read_almblm_band(path, info, eferm):
    """Band-mode almblm read: identical to read_almblm but the header eferm
    record is a placeholder the band path ignores; the Fermi shift uses the
    eferm passed in (from case.indmftpr, dmftproj.f:576-581). read_almblm reads
    the 4th record as eferm, so overwrite that record with the indmftpr value."""
    text = open(path).read().splitlines()
    text[3] = ' %.16E' % eferm          # Reader record 4 (physical line index 3)
    tmp = tempfile.NamedTemporaryFile('w', suffix='.almblm', delete=False)
    tmp.write('\n'.join(text) + '\n')
    tmp.close()
    try:
        return read_almblm(tmp.name, info, projectors=True)
    finally:
        os.unlink(tmp.name)


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

    crorbs = _build_crorbs(info)
    orbs = _build_orbs(info)
    ncrorb = len(crorbs)
    norb = len(orbs)

    # proj_mode 0 energy window [e_bot, e_top] per k (band path: proj_mode 0).
    windows = []
    for ik in range(nk):
        kp = spins[0]['kp'][ik]
        windows.append(select_window(kp['nbmin'], kp['nbmax'], kp['eband'],
                                     info['e_bot'], info['e_top']))

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

    # rot_projectmat local rotation: the op mapping the sort representative onto
    # this atom (set_projections.f via rot_projectmat). Same as ctqmcout/parproj.
    rotloc_op = {}
    for o in orbs:
        iref = sum(info['mult'][:o['sort'] - 1]) + 1
        rotloc_op[o['atom']] = next(op for op in ops
                                    if op['perm'][iref - 1] == o['atom'])

    # ---- correlated projector mat_rep, then Loewdin orthonormalize ----
    mat_rep = {}
    for icr, cr in enumerate(crorbs):
        l = cr['l']
        atom = cr['atom']
        sort = cr['sort']
        d = 2 * l + 1
        transmat = transmats[(l, sort)]
        ismix = mixing[(l, sort)]
        op = rotloc_op[atom]
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

    for ik in range(nk):
        incl, nb_bot, nb_top = windows[ik]
        if not incl:
            continue
        blocks = []
        layout = []
        for icr, cr in enumerate(crorbs):
            l = cr['l']
            if mixing[(l, cr['sort'])]:
                blocks.append(mat_rep[(icr, ik)])
                layout.append(((icr, ik), 2 * (2 * l + 1)))
            else:
                for ispin in range(ns):
                    blocks.append(mat_rep[(icr, ik, ispin)])
                    layout.append(((icr, ik, ispin), 2 * l + 1))
        D = np.vstack(blocks)
        # match the Fortran's exact BLAS calls (orthogonal_wannier_SO): the
        # near-singular O^{-1/2} amplifies any last-bit difference, so use the
        # identical ZGEMM trans flags rather than numpy's conj-transpose copies.
        O = zgemm(1.0, D, D, trans_b=2)       # ZGEMM('N','C'): D @ D^H
        S = _sqrt_inv(O)
        D_orth = zgemm(1.0, S, D)             # ZGEMM('N','N'): O^{-1/2} @ D
        row = 0
        for key, nrows in layout:
            mat_rep[key] = D_orth[row:row + nrows, :]
            row += nrows

    # ---- raw Theta projectors matn_rep per included orbital (parproj path) ----
    matn_reps = {}
    for o in orbs:
        l = o['l']
        transmat = transmats[(l, o['sort'])]
        op = rotloc_op[o['atom']]
        rot = dmat(l, op['a'], op['b'], op['g'], float(op['iprop']))
        if mixing[(l, o['sort'])]:
            matn_reps[o['atom']] = _build_matn_rep_mixing(
                o, info, spins, windows, transmat, rot)
        else:
            matn_reps[o['atom']] = _build_matn_rep(
                o, info, spins, ns, windows, transmat, rot)

    # ---- write ----
    with open(case + '.outband', 'w') as f:
        f.write('%6d\n' % nkband)

        # number of bands per k (skip is=2 under SP+SO)
        for ispin in range(ns):
            if ifSP and ifSO and ispin == 1:
                continue
            for ik in range(nk):
                incl, nb_bot, nb_top = windows[ik]
                f.write('%6d\n' % abs(nb_top - nb_bot + 1))

        # correlated projector block: DO ik, DO icrorb (Loewdin orthonormalized).
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

        # H(k) eigenvalues (skip is=2 under SO)
        for ispin in range(ns):
            if ifSO and ispin == 1:
                continue
            for ik in range(nk):
                incl, nb_bot, nb_top = windows[ik]
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
                incl, nb_bot, nb_top = windows[ik]
                if not incl:
                    continue
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

        # k-labels
        for i, (pos, name) in enumerate(labels, start=1):
            f.write('%6d%6d%s\n' % (i, pos, name))
