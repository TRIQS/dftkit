.. _welcome:

dftkit - DFT converters for TRIQS
*********************************

.. sidebar:: dftkit |PROJECT_VERSION|

   This is the homepage of dftkit |PROJECT_VERSION|.
   For changes see the :ref:`changelog page <changelog>`.

      .. image:: _static/logo_github.png
         :width: 75%
         :align: center
         :target: https://github.com/triqs/dftkit


dftkit is an official :ref:`TRIQS <triqslibs:welcome>` application that provides
converters turning the output of various DFT codes into the TRIQS-compatible HDF5
format required for DFT+DMFT calculations. Supported codes include
`Elk <https://elk.sourceforge.io/>`_,
`VASP <https://www.vasp.at>`_ (with the PLOVasp projection tools),
`Wien2k <http://www.wien2k.at>`_ (with the ``dmftproj`` projector generator),
`Quantum Espresso <https://www.quantum-espresso.org/>`_,
`Wannier90 <http://www.wannier.org/>`_ and generic tight-binding Hamiltonians.
It provides a modern, modular successor to the converters historically shipped with
`DFTTools <https://triqs.github.io/dft_tools>`_.

Learn how to use dftkit in the :ref:`documentation`.


.. image:: _static/logo_cea.png
   :width: 14%
   :target: http://ipht.cea.fr

.. image:: _static/logo_x.png
   :width: 14%
   :target: "https://www.cpht.polytechnique.fr

.. image:: _static/logo_cnrs.png
   :width: 14%
   :target: https://www.cnrs.fr

.. image:: _static/logo_erc.jpg
   :width: 14%

.. image:: _static/logo_flatiron.png
   :width: 20%
   :target: https://www.simonsfoundation.org/flatiron

.. image:: _static/logo_simons.jpg
   :width: 20%
   :target: https://www.simonsfoundation.org


.. toctree::
   :maxdepth: 2
   :hidden:

   install
   documentation
   issues
   ChangeLog.md
   about
