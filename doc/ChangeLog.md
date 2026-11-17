(changelog)=

# Changelog

## Unstable

Changes on the `unstable` branch that are not part of a release yet. Rename
this section to the version number when the next release is cut.

### VASP
* Add `KPOINTS_OPT` band conversion from `vaspout.h5`: when `LOCPROJ_OPT` data are available, the converter writes `dft_bands_input` for band/spectral workflows, applying the same PLO config settings (`EWINDOW`, `TRANSFORM`, `NORMALIZE`, and optional `EFERMI`) as the regular VASP conversion, and stores the high-symmetry k-path labels
* Run DMFT on the irreducible Brillouin zone: `convert_dft_input(use_ibz=...)` reduces the data to the irreducible k-points, sets `symm_op = 1` and writes the correlated-shell symmetry operations (real-harmonic rotations for `l = 0, 1, 2` and atom permutations, built from the space group in `vaspout.h5`) to `dft_symmcorr_input`, so that DFTTools and ModEST restore the full-BZ average by symmetrization. Enabled by default when VASP used symmetry and `vaspout.h5` holds the symmetry data, and only when it is known to work: correlated shells with `l <= 2`, a symmetry-closed orbital set with an orthonormal `TRANSFORM`, all equivalent atoms projected, no spin-orbit coupling, and the symmetrized IBZ occupation and local levels reproduce the full grid within `ibz_tol` (1e-3). Otherwise the converter falls back to the full grid with a note; `use_ibz=True` raises instead. The VASP `Driver` passes `use_ibz` through for CSC runs
* PLOVASP projected DOS (`DOSMESH`): fix the orbital-resolved DOS. With symmetry (`vaspout.h5` input) the tetrahedra refer to the IBZ k-points while the projectors are on the full grid; the DOS used the projector weights of one k-point per star only, which broke the degeneracy of symmetry-equivalent orbitals (e.g. t2g / eg of NiO). The weights are now averaged over the star. Also pair the projector weights with the right band when the band window starts at different bands at different k-points (energy window without `BANDS`). The total DOS of a shell and the PLOs written for DFTTools are unchanged
* PLOVASP: pair the projectors with the right bands in the density matrix and local Hamiltonian printed during PLO generation when the band window starts at different bands at different k-points (energy window without `BANDS`). The printed occupations were wrong in that case (NiO with `EWINDOW = -7 10`: d-shell trace 9.20 instead of 8.50); the PLOs written for DFTTools were not affected
* PLOVASP projected DOS: evaluate the tetrahedron weights for many energies at once (`dos_tetra_weights_3d` accepts an array of energies), about 5x faster PLO generation with `DOSMESH`

### Wien2k
* Support non-spin-orbit band projectors in the Python `case.outband` generator, including independent spin windows; correct the spin ordering when converting their partial projectors and read wrapped atom permutations from Wien2k symmetry files
* Support band-index projection modes 1 and 2 in the Python `case.outband` generator
* Use independent selected and below-window bands for ordinary non-spin-orbit projector generation in each spin channel
* Support the selection of individual irreducible representations of a correlated shell in `case.indmftpr` (e.g. only t2g or eg of a cubic d shell, or j=5/2 of a `fromfile` |j, m_j> basis) in the Python `case.ctqmcout`, `case.symqmc` and `case.outband` generators, matching native dmftproj. As in dmftproj, a spin-diagonal basis under spin-orbit coupling keeps the whole shell. Selecting more than one irrep of a shell is refused, since the Wien2k converter cannot read the resulting files
* Add pure-Python band-mode projectors, projection modes 1/2, and custom `fromfile` bases, including spin-mixing bases, to the Wien2k projector generators
* Read the high-symmetry k-path labels from the end of `case.outband` and store them as `kpts_labels` / `kpts_labels_idx` in `dft_bands_input`, matching the VASP band conversion

### test
* Check the Wien2k frame convention on full-BZ tellurium with spin-orbit coupling: the symmetry matrices from `case.symqmc` leave the global-frame local Hamiltonian invariant, as `SumkDFT` uses them, and applying them in the local frames is detected. Document the convention for `rot_mat` and the symmetry matrices in the HDF5 reference
* Check that the Wannier and Bloch basis conversions of the same wannier90 run describe the same operator, by downfolding the Bloch Hamiltonian with the projectors and comparing with the Wannier one, `P(k) H(k) P(k)^dag = H_W(k)`. `h5diff` only compares each archive with its own frozen reference and cannot catch an error present in both modes

## Version 4.0.1

dftkit version 4.0.1 is a patch release. It fixes the conversion of Wien2k
band partial projectors for spin-polarized runs without spin-orbit coupling,
and is the first dftkit release published on PyPI (`pip install triqs_dftkit`).

We thank all contributors: Christopher Albert, Thomas Hahn, Alexander Hampel, Nils Wentzell

Find below an itemized list of changes in this release.

### Wien2k
* Read the band partial projectors of spin-polarized runs without spin-orbit coupling in the order native `dmftproj -band` writes them (real parts of both spins, then imaginary parts). `proj_mat_all` in `dft_bands_input` was wrong before for such runs

### VASP
* Write the SP flag, not the number of spin channels, to `dft_symmcorr_input/SP`. No calculation used the old value

### ghaction
* Publish to PyPI on release tags, with the workflow of DFTTools

### doc
* Document the PyPI package in the install instructions; fix the clone URL
* Add a version switcher flyout and a canonical URL to the sphinx docs
* Fix the documentation website URLs and add the Zenodo DOI badge to the README


## Version 4.0.0

dftkit version 4.0.0 is the first release of the toolkit and is compatible with
TRIQS version 4.0.0. dftkit bundles the DFT-code interfaces that turn the output
of various DFT codes into the TRIQS-compatible HDF5 format consumed by
[DFTTools](https://triqs.github.io/dft_tools), providing a modern, modular
successor to the converters historically shipped with DFTTools. For each
supported code it provides a converter and, where available, a driver for
charge self-consistent DFT+DMFT calculations.

We thank all contributors: Markus Aichhorn, Jennifer Coulter, Michel Ferrero, Olivier Gingras, Thomas Hahn, Alexander Hampel, Alyn James, Harrison LaBollita, Henri Menke, Olivier Parcollet, Oleg E. Peil, Malte Schüler, Priyanka Seth, Dylan Simon, Nils Wentzell, Manuel Zingl

Find below an itemized list of changes in this release.

### General
* Port the DFT converters (Elk, Wien2k, Wannier90, generic H(k), VASP/PLOVasp) from the `triqs_dft_tools` repository
* Port the `dmftproj` Fortran component used by the Wien2k interface from `triqs_dft_tools`

### Drivers
* Add drivers that wrap the DFT executable for charge self-consistent DFT+DMFT runs, covering Quantum Espresso, Wien2k, Abinit and VASP

### VASP
* Add a VASP driver for charge self-consistent DFT+DMFT calculations
* Warn on misplaced or unknown tags in the PLOVASP configuration
* Read VASP `ICHARG=5` miscellaneous input from `vaspout.h5`

### Wannier90
* Add ABINIT support to `Wannier90Converter` for charge self-consistent calculations

### Fix
* Fix a bug in the `deltaN` write for the Quantum Espresso and Abinit interfaces

### doc
* Document the converter + driver architecture and the per-code interface table
* Rework the landing page and documentation index

### build
* Compatible with TRIQS 4.0.0 and the latest app4triqs/4.0.x skeleton
