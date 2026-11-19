[![build](https://github.com/TRIQS/dftkit/workflows/build/badge.svg)](https://github.com/TRIQS/dftkit/actions?query=workflow%3Abuild)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21691784-blue.svg)](https://doi.org/10.5281/zenodo.21691784)

# triqs_dftkit

DFT converters for [TRIQS](https://triqs.github.io), enabling DFT+DMFT calculations by converting output from various DFT codes to TRIQS-compatible HDF5 format.

## Supported DFT Codes

| Converter | Description |
|-----------|-------------|
| **Elk** | Elk all-electron DFT code |
| **VASP** | Vienna Ab initio Simulation Package (includes PLOVasp tools) |
| **Wien2k** | Wien2k LAPW code (includes the dmftproj projector generator, in Python) |
| **Wannier90** | Wannier tight-binding models |
| **Quantum Espresso** | QE plane-wave DFT code |
| **H(k)** | General tight-binding Hamiltonians |


## Documentation

Documentation is available at https://triqs.github.io/dftkit

## License

triqs_dftkit is released under the GNU General Public License v3.
