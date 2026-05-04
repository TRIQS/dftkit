(changelog)=

# Changelog

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
* Add `KPOINTS_OPT` band conversion from `vaspout.h5`: when `LOCPROJ_OPT` data are available, the converter writes `dft_bands_input` for band/spectral workflows, applying the same PLO config settings (`EWINDOW`, `TRANSFORM`, `NORMALIZE`, and optional `EFERMI`) as the regular VASP conversion, and stores the high-symmetry k-path labels
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
