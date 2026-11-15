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

# Spin-polarized band projectors without spin-orbit coupling (native dmftproj)
file(COPY CaOs2_noso_partial.ctqmcout.gz CaOs2_noso_partial.symqmc.gz
     CaOs2_band_noso.outband.gz CaOs2_band_noso_unequal.outband.gz
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME Py_wien2k_band_noso
         COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wien2k_band_noso.py
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_property(TEST Py_wien2k_band_noso APPEND PROPERTY ENVIRONMENT
             PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})
