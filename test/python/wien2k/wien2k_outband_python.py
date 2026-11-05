################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

# Verify the pure-Python case.outband generator (triqs_dftkit.wien2k.outband)
# reproduces the dmftproj Fortran band-structure projectors to machine
# precision: regenerate case.outband from the band-mode almblm/struct/dmftsym,
# the k-path (klist_band) and the band indmftpr, and compare to the native
# dmftproj -band reference (full-precision cubic-d template).
#
# case.outband feeds convert_bands_input. It is the band analog of
# case.ctqmcout: the same correlated-shell projector construction (raw Alm/Clm
# projector -> rot_projectmat -> fromfile/cubic transmat -> Loewdin
# orthonormalization) plus the raw Theta projectors of the parproj path, written
# in the outband layout. The CaOs2_band fixture is SP+SO with the spin-mixing
# fromfile (|j, m_j>) Os d shell (jbasis_d.dat); the wide window (-2.0, 3.0)
# keeps 50 bands > 20 correlated spin-orbitals, so the Loewdin overlap is full
# rank and the projectors match to machine precision.
#
# Band-mode specifics (dmftproj.f:557-581): nkband, the number of k-points along
# the plotted path, is read from CaOs2_band.klist_band; the Fermi energy is the
# last line of CaOs2_band.indmftpr (the band almblm keeps a placeholder eferm
# record that the reader skips).
#
# CaOs2_band_cubic runs the same band inputs with a cubic (spin-diagonal) Os d
# shell and an uncorrelated included Ca d shell, covering the non-mixing
# projector path; its reference is native dmftproj -band output.

import gzip
import os
import shutil
import tempfile
from unittest.mock import patch
import numpy as np
from h5 import HDFArchive

from triqs_dftkit.wien2k import Converter, ctqmcout, oubwin, outband, symqmc

HERE = os.path.dirname(os.path.abspath(__file__))


def _gunzip(src, dst):
    with gzip.open(src, 'rb') as fi, open(dst, 'wb') as fo:
        shutil.copyfileobj(fi, fo)


def _compare(got_path, ref_lines):
    got = open(got_path).read().split()
    ref = ''.join(ref_lines).split()
    assert len(got) == len(ref), (len(got), len(ref))
    floats = []
    for a, b in zip(got, ref):
        if a.lstrip('-').isdigit() and b.lstrip('-').isdigit():
            assert int(a) == int(b), (a, b)
        elif _isfloat(a) and _isfloat(b):
            floats.append((float(a), float(b)))
        else:
            assert a == b, (a, b)        # k-label tokens
    a, b = np.array(floats).T
    np.testing.assert_allclose(a, b, atol=1e-11, rtol=0, equal_nan=False)


def _isfloat(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def _check(case, inputs=None):
    inputs = inputs or case
    tmp = tempfile.mkdtemp()
    try:
        shutil.copy(os.path.join(HERE, f'{case}.indmftpr'),
                    os.path.join(tmp, f'{case}.indmftpr'))
        shutil.copy(os.path.join(HERE, f'{inputs}.klist_band'),
                    os.path.join(tmp, f'{case}.klist_band'))
        shutil.copy(os.path.join(HERE, 'jbasis_d.dat'),
                    os.path.join(tmp, 'jbasis_d.dat'))
        for ext in ('struct', 'dmftsym'):
            shutil.copy(os.path.join(HERE, f'CaOs2.{ext}'),
                        os.path.join(tmp, f'{case}.{ext}'))
        for ext in ('almblmup', 'almblmdn'):
            _gunzip(os.path.join(HERE, f'{inputs}.{ext}.gz'),
                    os.path.join(tmp, f'{case}.{ext}'))
        # Band parsing must also work without writable temporary storage.
        with patch('tempfile.NamedTemporaryFile', side_effect=AssertionError('temporary file created')):
            outband.write_outband(os.path.join(tmp, case))
        with gzip.open(os.path.join(HERE, f'{case}.outband.gz'), 'rt') as fh:
            _compare(os.path.join(tmp, f'{case}.outband'), fh.readlines())
        with open(os.path.join(tmp, f'{case}.klist_band'), 'w') as fh:
            fh.write('GAMMA            0    0    0    8  2.0\nX                4    0    0    8  2.0\nEND\n')
        assert outband.read_k_list(os.path.join(tmp, f'{case}.klist_band')) == (2, [(1, 'GAMMA'), (2, 'X')])
        try:
            outband.write_outband(os.path.join(tmp, case))
        except ValueError as exc:
            assert 'k-point counts differ' in str(exc)
        else:
            raise AssertionError('inconsistent band-path dimensions must be rejected')
    finally:
        shutil.rmtree(tmp)


_check('CaOs2_band')
_check('CaOs2_band_cubic', inputs='CaOs2_band')


def _native_band_data(path, converter):
    """Decode native blocks as packed real and imaginary vectors, across spins."""
    with gzip.open(path, 'rt') as fh:
        tokens = iter(fh.read().split())

    def take(n):
        return np.fromiter((float(next(tokens)) for _ in range(n)), float, count=n)

    nk = int(next(tokens))
    ns = converter.n_spin_blocs
    counts = take(nk * ns).astype(int).reshape(ns, nk).T
    nb = counts.max()
    corr = np.zeros((nk, ns, converter.n_corr_shells,
                     max(s['dim'] for s in converter.corr_shells), nb), complex)

    def blocks(dim, ik):
        sizes = dim * counts[ik]
        data = take(sum(sizes)) + 1j * take(sum(sizes))
        return [a.reshape(dim, n) for a, n in
                zip(np.split(data, np.cumsum(sizes)[:-1]), counts[ik])]

    for ik in range(nk):
        for ish, shell in enumerate(converter.corr_shells):
            for isp, block in enumerate(blocks(shell['dim'], ik)):
                corr[ik, isp, ish, :shell['dim'], :counts[ik, isp]] = block
    hopping = np.zeros((nk, ns, nb, nb), complex)
    for isp in range(ns):
        for ik in range(nk):
            n = counts[ik, isp]
            hopping[ik, isp, :n, :n] = np.diag(take(n) * converter.energy_unit)
    nr = take(converter.n_shells).astype(int)
    partial = np.zeros((nk, ns, converter.n_shells, max(nr),
                        max(s['dim'] for s in converter.shells), nb), complex)
    for ish, shell in enumerate(converter.shells):
        for ik in range(nk):
            for ir in range(nr[ish]):
                for isp, block in enumerate(blocks(shell['dim'], ik)):
                    partial[ik, isp, ish, ir, :shell['dim'], :counts[ik, isp]] = block
    assert list(tokens) == ['1', '1GAMMA']
    return dict(n_k=nk, n_orbitals=counts, proj_mat=corr, hopping=hopping,
                n_parproj=nr, proj_mat_all=partial)


def _check_nonso(case, main_input, paramagnetic=False):
    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, case)
        for ext in ('struct', 'dmftsym', 'outputs'):
            shutil.copy(os.path.join(HERE, f'CaOs2.{ext}'), prefix + '.' + ext)
        for spin in ('up',) if paramagnetic else ('up', 'dn'):
            _gunzip(os.path.join(HERE, f'CaOs2_noso.almblm{spin}.gz'),
                    prefix + ('.almblm' if paramagnetic else '.almblm' + spin))
        shutil.copy(os.path.join(HERE, 'CaOs2_band.klist_band'), prefix + '.klist_band')

        # The ordinary DFT archive supplies the shell metadata for band conversion.
        with open(os.path.join(HERE, main_input + '.indmftpr')) as fh:
            main = fh.read()
        if case.endswith('_complex'):
            main = main.replace('cubic', 'complex')
        with open(prefix + '.indmftpr', 'w') as fh:
            fh.write(main)
        ctqmcout.write_ctqmcout(prefix)
        symqmc.write_symqmc(prefix)
        oubwin.write_oubwin(prefix)
        converter = Converter(filename=prefix)
        converter.convert_dft_input()

        shutil.copy(os.path.join(HERE, case + '.indmftpr'), prefix + '.indmftpr')
        outband.write_outband(prefix)
        reference = os.path.join(HERE, case + '.outband.gz')
        with gzip.open(reference, 'rt') as fh:
            _compare(prefix + '.outband', fh.readlines())
        converter.convert_bands_input()
        expected = _native_band_data(reference, converter)
        assert converter.n_spin_blocs == (1 if paramagnetic else 2)
        if case.endswith('_unequal'):
            assert expected['n_orbitals'][0, 0] != expected['n_orbitals'][0, 1]
        with HDFArchive(converter.hdf_file, 'r') as archive:
            band = archive[converter.bands_subgrp]
            for key, value in expected.items():
                assert np.isfinite(band[key]).all(), key
                np.testing.assert_allclose(band[key], value, atol=1e-11, rtol=0,
                                           equal_nan=False, err_msg=key)
            assert band['kpts_labels'] == ['GAMMA']
            np.testing.assert_array_equal(band['kpts_labels_idx'], [0])
            for isp in range(converter.n_spin_blocs):
                n = band['n_orbitals'][0, isp]
                projectors = np.vstack([band['proj_mat'][0, isp, ish, :s['dim'], :n]
                                         for ish, s in enumerate(converter.corr_shells)])
                np.testing.assert_allclose(projectors @ projectors.conj().T,
                                           np.eye(len(projectors)), atol=1e-12, rtol=0)


_check_nonso('CaOs2_band_noso', 'CaOs2_noso_partial')
_check_nonso('CaOs2_band_noso_pm', 'CaOs2_noso_partial', paramagnetic=True)
_check_nonso('CaOs2_band_noso_t2g', 'CaOs2_irrep_t2g')
_check_nonso('CaOs2_band_noso_eg', 'CaOs2_irrep_eg')
_check_nonso('CaOs2_band_noso_complex', 'CaOs2_noso_partial', paramagnetic=True)
_check_nonso('CaOs2_band_noso_unequal', 'CaOs2_noso_partial')

print('wien2k_outband_python: ok')
