(changelog)=

# Changelog

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

### Wannier90
* Add ABINIT support to `Wannier90Converter` for charge self-consistent calculations
* Read VASP `ICHARG=5` miscellaneous input from `vaspout.h5`

### Fix
* Fix a bug in the `deltaN` write for the Quantum Espresso and Abinit interfaces

### doc
* Document the converter + driver architecture and the per-code interface table
* Rework the landing page and documentation index

### build
* Compatible with TRIQS 4.0.0 and the latest app4triqs/4.0.x skeleton
