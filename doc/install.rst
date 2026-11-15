.. highlight:: bash

.. _install:

Install triqs_dftkit
********************

.. _pypi:

PyPI
====

We provide a Python package through the `Python Package Index <https://pypi.org/>`_ (PyPI). After :ref:`installing TRIQS <triqslibs:triqs_install>` you can install dftkit with::

        pip install triqs_dftkit

The PyPI package requires a TRIQS installation with a matching major.minor version number. It contains the Python converters and drivers only. The Fortran ``dmftproj`` program of the Wien2k interface is not included; to use it, compile dftkit from source as described below.

See also `pypi.org/project/triqs_dftkit <https://pypi.org/project/triqs_dftkit/>`_.

Compiling triqs_dftkit from source
==================================

.. note:: To guarantee reproducibility in scientific calculations we strongly recommend the use of a stable `release <https://github.com/TRIQS/triqs/releases>`_ of both TRIQS and its applications.

Prerequisites
-------------

#. The :ref:`TRIQS <triqslibs:welcome>` library, see :ref:`TRIQS installation instruction <triqslibs:triqs_install>`.
   In the following, we assume that TRIQS is installed in the directory ``path_to_triqs``.

Installation steps
------------------

#. Download the source code of the latest stable version by cloning the ``TRIQS/dftkit`` repository from GitHub::

     $ git clone https://github.com/TRIQS/dftkit triqs_dftkit.src

#. Create and move to a new directory where you will compile the code::

     $ mkdir triqs_dftkit.build && cd triqs_dftkit.build

#. Ensure that your shell contains the TRIQS environment variables by sourcing the ``triqsvars.sh`` file from your TRIQS installation::

     $ source path_to_triqs/share/triqs/triqsvars.sh

#. In the build directory call cmake, including any additional custom CMake options, see below::

     $ cmake ../triqs_dftkit.src

#. Compile the code, run the tests and install the application::

     $ make test
     $ make install

Version compatibility
---------------------

Keep in mind that the version of ``triqs_dftkit`` must be compatible with your TRIQS library version,
see :ref:`TRIQS website <triqslibs:versions>`.
In particular the Major and Minor Version numbers have to be the same.
To use a particular version, go into the directory with the sources, and look at all available versions::

     $ cd triqs_dftkit.src && git tag

Checkout the version of the code that you want::

     $ git checkout 2.1.0

and follow steps 2 to 4 above to compile the code.

Custom CMake options
--------------------

The compilation of ``triqs_dftkit`` can be configured using CMake-options::

    cmake ../triqs_dftkit.src -DOPTION1=value1 -DOPTION2=value2 ...

+-----------------------------------------------------------------+-----------------------------------------------+
| Options                                                         | Syntax                                        |
+=================================================================+===============================================+
| Specify an installation path other than path_to_triqs           | -DCMAKE_INSTALL_PREFIX=path_to_triqs_dftkit      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build in Debugging Mode                                         | -DCMAKE_BUILD_TYPE=Debug                      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Disable testing (not recommended)                               | -DBuild_Tests=OFF                             |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build the documentation                                         | -DBuild_Documentation=ON                      |
+-----------------------------------------------------------------+-----------------------------------------------+
