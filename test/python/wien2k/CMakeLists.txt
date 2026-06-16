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

# Spin-orbit + spin-polarized converter test (two equivalent Os atoms).
file(COPY CaOs2.ctqmcout CaOs2.symqmc CaOs2.struct CaOs2.outputs
          CaOs2.oubwinup CaOs2.parproj
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# dmftproj writes CaOs2.sympar byte-identical to CaOs2.symqmc for this fixture.
configure_file(CaOs2.symqmc ${CMAKE_CURRENT_BINARY_DIR}/CaOs2.sympar COPYONLY)

add_test(NAME Py_wien2k_soc_convert
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_soc_convert.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_soc_convert APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Pure-Python case.symqmc generator vs the shipped dmftproj Fortran output.
# The generator tests read source fixtures and write only to scratch directories.
add_test(NAME Py_wien2k_symqmc_python
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_symqmc_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_symqmc_python APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Generator on the two paths the cubic test does not reach: the spin-mixing
# fromfile (|j,m_j>) basis and an l=0 shell. Reuses CaOs2.dmftsym/struct; the
# references are the dmftproj matrix data (single precision, ~26 kB total).
add_test(NAME Py_wien2k_symqmc_mixing
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_symqmc_mixing.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_symqmc_mixing APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Generator on a non-centrosymmetric spin-orbit cell (tellurium, p shell), the
# dft_tools #148 regime the centrosymmetric d-shell fixtures never reached. The
# reference is the dmftproj matrix data (single precision, ~10 kB).
file(COPY Te.dmftsym Te.indmftpr Te.struct wien2k_symqmc_te.ref.npy
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME Py_wien2k_symqmc_noncentrosym
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_symqmc_noncentrosym.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_symqmc_noncentrosym APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Pure-Python case.oubwin band-window generator vs the dmftproj Fortran output.
# Reads the existing compressed almblm inputs from the source fixtures and
# compares both spin outputs against CaOs2.oubwinup in a scratch directory.
add_test(NAME Py_wien2k_oubwin_python
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_oubwin_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_oubwin_python APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Pure-Python correlated-shell projector generator vs the native fixture.
add_test(NAME Py_wien2k_ctqmcout_python
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_ctqmcout_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_ctqmcout_python APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Included-shell symmetry generator, including an uncorrelated Ca d shell.
add_test(NAME Py_wien2k_sympar_python
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_sympar_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_sympar_python APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

add_test(NAME Py_wien2k_parproj_python
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_parproj_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_parproj_python APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

add_test(NAME Py_wien2k_kpath_labels
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_kpath_labels.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_kpath_labels APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

add_test(NAME Py_wien2k_dmftproj_common_python
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_dmftproj_common_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_dmftproj_common_python APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

add_test(NAME Py_wien2k_outband_python
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_outband_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_outband_python APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Pure-Python generators on the non-spin-orbit path (ifSO=0, the common DMFT
# case): dmftproj writes (2l+1) matrices, per-spin orthonormalization. Validated
# against the precision-fixed dmftproj run without -so. Reuses CaOs2.struct/dmftsym.
file(COPY CaOs2_noso.indmftpr CaOs2_noso.almblmup.gz CaOs2_noso.almblmdn.gz
          CaOs2_noso.ctqmcout.gz CaOs2_noso.symqmc.gz
          CaOs2_noso_partial.indmftpr CaOs2_noso_partial.sympar.gz CaOs2_noso_partial.parproj.gz
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME Py_wien2k_noso_python
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_noso_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_noso_python APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Odd-l shells (p l=1, f l=3): the improper-rotation parity (-1)^l is invisible
# for the even-l (d) fixtures, so these pin dmat's use of iprop. p uses the
# standard almblm (wide window, full rank); f uses a high-band-cutoff almblm so
# the high-energy window carries real l=3 character and the Loewdin overlap is
# full rank. Both reproduce the precision-fixed dmftproj to machine precision.
file(COPY CaOs2_pshell.indmftpr CaOs2_pshell.ctqmcout.gz CaOs2_pshell.symqmc.gz
          CaOs2_fshell.indmftpr CaOs2_fshell.almblmup.gz CaOs2_fshell.almblmdn.gz
          CaOs2_fshell.ctqmcout.gz CaOs2_fshell.symqmc.gz
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME Py_wien2k_oddl_python
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_oddl_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_oddl_python APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})

# Capstone: the pure-Python run_dmftproj driver replaces the Fortran executable.
# Both tests run the driver through the converter and compare every HDF value
# at machine precision (1e-12). The Loewdin overlap must be full rank: the physical
# narrow window is rank-deficient (more correlated spin-orbitals than bands, so
# O^{-1/2} has non-unique null eigenvectors), which is why the wide window is
# used here. The narrow physical case is still locked exactly by the
# Fortran-output converter test wien2k_soc_convert.
#   Non-SO partial: convert_dft_input on the SP/non-SO oubwin path + convert_parproj_input.
#   SO full-rank:   convert_dft_input on a wide-window SO cell.
file(COPY CaOs2_noso_partial.indmftpr CaOs2_noso.almblmup.gz CaOs2_noso.almblmdn.gz
          CaOs2.struct CaOs2.dmftsym CaOs2.outputs wien2k_noso_convert.ref.h5
          CaOs2_socfull.indmftpr CaOs2.almblmup.gz CaOs2.almblmdn.gz
          wien2k_socfull_convert.ref.h5
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME Py_wien2k_noso_convert
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_noso_convert.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_noso_convert APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})
add_test(NAME Py_wien2k_socfull_convert
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_socfull_convert.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_socfull_convert APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})
