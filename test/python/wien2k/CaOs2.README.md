# CaOs2 SOC converter fixture

This fixture covers spin-polarized spin-orbit coupling (`dmftproj -so -sp`),
two equivalent correlated Os d shells, and sixteen magnetic symmetry
operations, eight flagged antiunitary. The per-shell time-reversal matrices
are identity matrices and the local rotation time-reversal flags are zero;
this fixture does not cover nontrivial local time-reversal rotations.
The retained input contains one
k-point with unit normalized weight. It does not test variation over k-points.

The original narrow window (recorded in the stack inputs at `5a805f7`),
-0.15 to 0.30 Ry relative to the Fermi energy,
selected bands 60–76 and could not accommodate the twenty correlated spin
orbitals. The shipped `CaOs2.indmftpr` uses -2.0 to 3.0 Ry and selects bands
43–92 (50 bands). The native producer gives a full twenty-dimensional
correlated subspace: the maximum stacked projector overlap error is about
3.1e-8. The test checks this numerical invariant independently of the golden
archive, including overlaps between different shells. It also checks the
eight operations that exchange Os atoms and the unitarity of the shell
symmetry matrices.

## Producer provenance

The compressed `CaOs2.almblmup.gz` and `CaOs2.almblmdn.gz` are the retained
Wien2k projector inputs already used by the dmftproj port stack. Their
decompressed content matches commit
`5a805f74e37854ea0f4576e39516ce47fafa25e7`, which introduced those inputs.
The original Wien2k version, initialization commands, and parent k-mesh were
not retained; reproducing the underlying DFT calculation
is therefore outside this fixture's documented provenance. No new k-points
or projector coefficients have been constructed for this correction.

The files include the complete inputs needed to regenerate the converter
fixture with the native Fortran producer: the two `almblm` files, `struct`,
`dmftsym`, `indmftpr`, and retained `case.cf_d_eg_t2g` basis template. The
original cubic basis is preserved. The fixture was regenerated with the native
CMake producer (GNU Fortran 15.3.0) from upstream `43435dc`; the converter
reference uses TRIQS 4.0.2. Use the producer built from this repository:

```sh
mkdir -p wienroot/SRC_templates CaOs2
cp /path/to/fixtures/case.cf_d_eg_t2g wienroot/SRC_templates/
export WIENROOT="$PWD/wienroot"
cd CaOs2
# Copy CaOs2.struct, CaOs2.dmftsym, CaOs2.indmftpr and CaOs2.outputs here.
gzip -dc /path/to/fixtures/CaOs2.almblmup.gz > CaOs2.almblmup
gzip -dc /path/to/fixtures/CaOs2.almblmdn.gz > CaOs2.almblmdn
/path/to/build/fortran/dmftproj/dmftproj -so -sp
```

The directory must be named `CaOs2`, since dmftproj derives the case name from
the working directory. Retain the resulting `.ctqmcout`, `.symqmc`, and
`.oubwinup`; `.oubwindn` is identical and the SOC converter reads `.oubwinup`.
Only trailing whitespace is removed from native text output.
The original `.struct` and `.outputs` supply the converter's miscellaneous
symmetry data. Regenerate the HDF reference with the same installed converter:

```python
from triqs_dftkit.wien2k import Converter
converter = Converter(filename='CaOs2')
converter.hdf_file = 'wien2k_soc_convert.ref.h5'
converter.convert_dft_input()
```

The archive comparison checks every entry at an absolute tolerance of `1e-12`
for conversion of the same fixture bytes. It does not use
`triqs.utility.h5diff`, which compares arrays at `1e-6` regardless of its
`precision` argument. The physical overlap check uses `1e-6` because the
retained Wien2k inputs and native text output have finite precision.
