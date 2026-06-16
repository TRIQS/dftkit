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
# Reuses the committed CaOs2.oubwin{up,dn} references; adds the truncated almblm
# inputs (header + first-sort band block, ~5 kB each) the selection reads.
file(COPY CaOs2.almblmup CaOs2.almblmdn DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME Py_wien2k_oubwin_python
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_oubwin_python.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_oubwin_python APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})
