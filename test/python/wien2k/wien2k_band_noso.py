################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# convert_bands_input on spin-polarized band projectors without spin-orbit
# coupling. Native dmftproj -band writes each projector block as the real parts
# of all spins followed by the imaginary parts of all spins (outband.f). The
# archive must hold the same complex projectors, decoded here independently.
#
# Inputs (native dmftproj -sp, CaOs2 without SO, Os d correlated, Ca d
# included): CaOs2_noso_partial.{ctqmcout,symqmc} for the shell metadata and
# CaOs2_band_noso{,_unequal}.outband along a one-point path. In the _unequal
# file the two spins keep different numbers of bands.

import gzip
import os
import shutil
import tempfile

import numpy as np
from h5 import HDFArchive

from triqs_dftkit.wien2k import Converter

HERE = os.path.dirname(os.path.abspath(__file__))


def _gunzip(src, dst):
    with gzip.open(src, 'rb') as fi, open(dst, 'wb') as fo:
        shutil.copyfileobj(fi, fo)


def _native_band_data(path, conv):
    with gzip.open(path, 'rt') as fh:
        tokens = iter(fh.read().split())

    def take(n):
        return np.fromiter((float(next(tokens)) for _ in range(n)), float, count=n)

    nk = int(next(tokens))
    ns = conv.n_spin_blocs
    counts = take(nk * ns).astype(int).reshape(ns, nk).T
    nb = counts.max()

    def blocks(dim, ik):
        sizes = dim * counts[ik]
        data = take(sum(sizes)) + 1j * take(sum(sizes))
        return [a.reshape(dim, n) for a, n in
                zip(np.split(data, np.cumsum(sizes)[:-1]), counts[ik])]

    corr = np.zeros((nk, ns, conv.n_corr_shells,
                     max(s['dim'] for s in conv.corr_shells), nb), complex)
    for ik in range(nk):
        for ish, shell in enumerate(conv.corr_shells):
            for isp, block in enumerate(blocks(shell['dim'], ik)):
                corr[ik, isp, ish, :shell['dim'], :counts[ik, isp]] = block
    hopping = np.zeros((nk, ns, nb, nb), complex)
    for isp in range(ns):
        for ik in range(nk):
            n = counts[ik, isp]
            hopping[ik, isp, :n, :n] = np.diag(take(n) * conv.energy_unit)
    nr = take(conv.n_shells).astype(int)
    partial = np.zeros((nk, ns, conv.n_shells, max(nr),
                        max(s['dim'] for s in conv.shells), nb), complex)
    for ish, shell in enumerate(conv.shells):
        for ik in range(nk):
            for ir in range(nr[ish]):
                for isp, block in enumerate(blocks(shell['dim'], ik)):
                    partial[ik, isp, ish, ir, :shell['dim'], :counts[ik, isp]] = block
    return dict(n_k=nk, n_orbitals=counts, proj_mat=corr, hopping=hopping,
                n_parproj=nr, proj_mat_all=partial)


def _check(band):
    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, 'case')
        for ext in ('ctqmcout', 'symqmc'):
            _gunzip(os.path.join(HERE, f'CaOs2_noso_partial.{ext}.gz'), f'{prefix}.{ext}')
        reference = os.path.join(HERE, f'{band}.outband.gz')
        _gunzip(reference, prefix + '.outband')
        conv = Converter(filename=prefix)
        conv.convert_dft_input()
        conv.convert_bands_input()
        assert conv.SO == 0 and conv.n_spin_blocs == 2
        expected = _native_band_data(reference, conv)
        with HDFArchive(conv.hdf_file, 'r') as archive:
            data = archive[conv.bands_subgrp]
            for key, value in expected.items():
                np.testing.assert_allclose(data[key], value, atol=1e-12, rtol=0,
                                           err_msg=f'{band}: {key}')
        return expected


_check('CaOs2_band_noso')
counts = _check('CaOs2_band_noso_unequal')['n_orbitals']
assert counts[0, 0] != counts[0, 1]

print('wien2k_band_noso: ok')
