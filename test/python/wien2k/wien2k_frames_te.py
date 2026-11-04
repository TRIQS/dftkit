################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Frame convention of the Wien2k converter output on full-BZ tellurium (P3_121,
# three equivalent atoms with different local frames, SO + SP, cubic p).
#
# proj_mat and the dft_symmcorr_input matrices are in the global frame; rot_mat
# maps global to local. SumkDFT symmetrizes in the global frame and rotates to
# the local frame afterwards. On a full-BZ k-sum the local Hamiltonian
#     H_i = sum_k w_k P_i(k) H(k) P_i(k)^dag
# is already symmetric, so symmetrizing it in that order must change nothing,
# while applying the global-frame matrices to local-frame data must not.
#
# Input: Te_fbz.ctqmcout.gz and Te_fbz.oubwinup.gz (3x3x2 full-BZ mesh, see
# Te.README.md) with the shared Te.struct and Te.symqmc.gz. Energies in Ry.

import gzip
import os
import shutil
import tempfile

import numpy as np
from h5 import HDFArchive

from triqs_dftkit.wien2k import Converter

HERE = os.path.dirname(os.path.abspath(__file__))
TOL_SYMMETRIC = 1e-6
MIN_MISMATCH = 1e-3


def _gunzip(src, dst):
    with gzip.open(src, 'rb') as fi, open(dst, 'wb') as fo:
        shutil.copyfileobj(fi, fo)


def symmetrize(mats, symm, n_orbits):
    """Same operation as triqs_dft_tools.symmetry.Symmetry.symmetrize on matrices."""
    out = [np.zeros_like(m) for m in mats]
    for isym in range(symm['n_symm']):
        for iorb in range(n_orbits):
            jorb = symm['orb_map'][isym][iorb]
            m = mats[iorb].conj() if symm['time_inv'][isym] else mats[iorb]
            r = symm['mat'][isym][iorb]
            out[jorb] += r @ m @ r.conj().T / symm['n_symm']
    return out


def to_local(mats, rot_mat, rot_mat_time_inv, SO):
    """SumkDFT.rotloc(..., 'toLocal') on Hermitian matrices."""
    return [(r.conj().T @ m @ r).conj() if (tinv and SO) else r.conj().T @ m @ r
            for m, r, tinv in zip(mats, rot_mat, rot_mat_time_inv)]


def eig_drift(a, b):
    return max(np.max(np.abs(np.linalg.eigvalsh(x) - np.linalg.eigvalsh(y))) for x, y in zip(a, b))


tmp = tempfile.mkdtemp()
try:
    case = os.path.join(tmp, 'Te_fbz')
    shutil.copy(os.path.join(HERE, 'Te.struct'), case + '.struct')
    _gunzip(os.path.join(HERE, 'Te_fbz.ctqmcout.gz'), case + '.ctqmcout')
    _gunzip(os.path.join(HERE, 'Te.symqmc.gz'), case + '.symqmc')
    _gunzip(os.path.join(HERE, 'Te_fbz.oubwinup.gz'), case + '.oubwinup')
    conv = Converter(filename=case, repacking=False)
    conv.convert_dft_input()
    with HDFArchive(conv.hdf_file, 'r') as ar:
        dft = {k: ar['dft_input'][k] for k in ar['dft_input'].keys()}
        symm = {k: ar['dft_symmcorr_input'][k] for k in ar['dft_symmcorr_input'].keys()}
finally:
    shutil.rmtree(tmp)

assert dft['symm_op'] == 1 and dft['SO'] == 1 and dft['SP'] == 1
assert dft['n_k'] == 18 and len(dft['corr_shells']) == 3
shells = dft['corr_shells']
n_orbits = len(shells)

orbits = symm['orbits']
symm['orb_map'] = [[orbits.index({**orbits[i], 'atom': symm['perm'][isym][orbits[i]['atom'] - 1]})
                    for i in range(n_orbits)] for isym in range(symm['n_symm'])]

# The three atoms have pairwise different local frames, so a frame mix-up is visible.
rot = dft['rot_mat']
assert all(np.max(np.abs(rot[i] - rot[j])) > 0.1 for i in range(3) for j in range(i))

hglob = [np.zeros((s['dim'], s['dim']), complex) for s in shells]
for ik in range(dft['n_k']):
    n = dft['n_orbitals'][ik, 0]
    h = dft['hopping'][ik, 0, :n, :n] / dft['energy_unit']
    for i, s in enumerate(shells):
        p = dft['proj_mat'][ik, 0, i, :s['dim'], :n]
        hglob[i] += dft['bz_weights'][ik] * p @ h @ p.conj().T
hloc = to_local(hglob, rot, dft['rot_mat_time_inv'], dft['SO'])

# Equivalent atoms have the same local levels.
for i in range(1, n_orbits):
    np.testing.assert_allclose(np.linalg.eigvalsh(hglob[i]), np.linalg.eigvalsh(hglob[0]),
                               atol=TOL_SYMMETRIC, rtol=0)

# SumkDFT order: symmetrize the global-frame H_i, then rotate. Nothing moves.
sym_glob = symmetrize(hglob, symm, n_orbits)
residual = max(np.max(np.abs(a - b)) for a, b in zip(sym_glob, hglob))
assert residual < TOL_SYMMETRIC, residual
assert eig_drift(to_local(sym_glob, rot, dft['rot_mat_time_inv'], dft['SO']), hloc) < TOL_SYMMETRIC

# Global-frame matrices applied to local-frame data: wrong frame, detected.
drift = eig_drift(symmetrize(hloc, symm, n_orbits), hloc)
assert drift > MIN_MISMATCH, drift

# Matrices rotated into the local frames applied to global-frame data: also wrong.
symm_local = dict(symm)
symm_local['mat'] = [[rot[symm['orb_map'][isym][i]].conj().T @ symm['mat'][isym][i]
                      @ (rot[i].conj() if symm['time_inv'][isym] else rot[i])
                      for i in range(n_orbits)] for isym in range(symm['n_symm'])]
drift = eig_drift(symmetrize(hglob, symm_local, n_orbits), hglob)
assert drift > MIN_MISMATCH, drift

print('wien2k_frames_te: ok')
