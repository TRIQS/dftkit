# Wien2k Converter Tests

# Copy input files
file(COPY SrVO3.ctqmcout SrVO3.symqmc SrVO3.parproj SrVO3.struct
     SrVO3.outputs SrVO3.oubwin SrVO3.pmat SrVO3.sympar
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Add test
add_test(NAME Py_wien2k_convert
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_convert.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_convert APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Spin-orbit + spin-polarized converter test (CaOs2: cubic, two equivalent
# correlated atoms, 8 time-reversal symmetry operations).
file(COPY CaOs2.ctqmcout CaOs2.symqmc CaOs2.struct CaOs2.outputs
          CaOs2.oubwinup CaOs2.oubwindn wien2k_soc_convert.ref.h5
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME Py_wien2k_soc_convert
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_soc_convert.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_soc_convert APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Pure-Python case.symqmc generator vs the dmftproj Fortran output (reuses the
# SOC reference h5; adds only the small text inputs the generator reads).
file(COPY CaOs2.dmftsym CaOs2.indmftpr DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME Py_wien2k_symqmc_python
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_symqmc_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_symqmc_python APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Generator on the two paths the cubic test does not reach: the spin-mixing
# fromfile (|j,m_j>) basis and an l=0 shell. Reuses CaOs2.dmftsym/struct; the
# references are the dmftproj matrix data (single precision, ~26 kB total).
file(COPY CaOs2_jbasis.indmftpr CaOs2_l0.indmftpr jbasis_d.dat
          wien2k_symqmc_jbasis.ref.npy wien2k_symqmc_l0.ref.npy
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME Py_wien2k_symqmc_mixing
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_symqmc_mixing.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_symqmc_mixing APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Pure-Python case.oubwin band-window generator vs the dmftproj Fortran output.
# Reuses the committed CaOs2.oubwin{up,dn} references and the gzipped CaOs2
# almblm inputs (shared with the ctqmcout test). CaOs2_mode1 adds the band-index
# projection path (proj_mode 1, window 60..76) against its own dmftproj refs.
file(COPY CaOs2.almblmup.gz CaOs2.almblmdn.gz
          CaOs2_mode1.indmftpr CaOs2_mode1.oubwinup CaOs2_mode1.oubwindn
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME Py_wien2k_oubwin_python
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_oubwin_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_oubwin_python APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Wide-window fromfile (|j, m_j>) inputs and references shared by the ctqmcout,
# sympar and parproj generator tests (the spin-mixing dft_tools #148 path).
file(COPY CaOs2_jbasis_full.indmftpr CaOs2_jbasis_full.ctqmcout.gz
          CaOs2_jbasis_full.sympar.gz CaOs2_jbasis_full.parproj.gz
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Pure-Python case.ctqmcout correlated-shell projector generator vs the dmftproj
# Fortran output. Reuses CaOs2.ctqmcout / indmftpr / struct / dmftsym and the
# gzipped almblm; floats compared at 1e-11 (full-rank wide window). The
# CaOs2_mode{1,2}_full fixtures exercise the band-index projection modes
# (proj_mode 1 scan, proj_mode 2 explicit indices), both resolving to bands
# 43..92, against their own precision-fixed dmftproj references.
file(COPY CaOs2_mode1_full.indmftpr CaOs2_mode1_full.ctqmcout.gz
          CaOs2_mode2_full.indmftpr CaOs2_mode2_full.ctqmcout.gz
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME Py_wien2k_ctqmcout_python
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_ctqmcout_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_ctqmcout_python APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Pure-Python case.sympar partial-shell symmetry generator vs the dmftproj
# Fortran output. The CaOs2_partial fixture adds an uncorrelated Ca d-shell to
# the correlated Os d-shells; inputs reuse the CaOs2 struct/dmftsym/almblm.
file(COPY CaOs2_partial.indmftpr CaOs2_partial.sympar.gz
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME Py_wien2k_sympar_python
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_sympar_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_sympar_python APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Pure-Python case.parproj partial projectors + density matrices vs the dmftproj
# Fortran output. Same CaOs2_partial fixture as the sympar test.
file(COPY CaOs2_partial.parproj.gz DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME Py_wien2k_parproj_python
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_parproj_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_parproj_python APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Unit tests for the shared dmftproj machinery (_dmftproj), the module the five
# case.* generators are built on. Reuses CaOs2 / CaOs2_partial indmftpr + dmftsym.
add_test(NAME Py_wien2k_dmftproj_common
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_dmftproj_common_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_dmftproj_common APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})
# Pure-Python case.outband band-structure projector generator (the -band path
# that feeds convert_bands_input) vs the precision-fixed dmftproj -band output.
# The CaOs2_band fixture is the spin-mixing fromfile (|j, m_j>) Os d shell on a
# 4-point k-path; the band almblm carries a single representative k-point and the
# Fermi energy comes from the last line of the band indmftpr. Reuses
# CaOs2.struct/dmftsym and jbasis_d.dat; floats compared at 1e-11 (full-rank
# wide window).
file(COPY CaOs2_band.indmftpr CaOs2_band.klist_band
          CaOs2_band.almblmup.gz CaOs2_band.almblmdn.gz CaOs2_band.outband.gz
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME Py_wien2k_outband_python
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_outband_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_outband_python APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})
