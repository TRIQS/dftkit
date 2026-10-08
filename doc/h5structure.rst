.. _hdfstructure:

standardized hdf5 structure
===========================

All DFT input data produced by the dftkit converters is stored in a single hdf5
archive, using the conventions of the TRIQS ``HDFArchive``. This page documents
the layout
of that archive: **which top-level group** each piece of data goes into, what it
means, and **which converter writes it**.

Every converter writes into the same archive file (``seedname.h5``), one
top-level group per kind of data. The group names are not hardcoded: on the
writing side they are constructor arguments of the converters (``dft_subgrp``,
``misc_subgrp``, ``symmcorr_subgrp``, …), and on the reading side of
:py:class:`SumkDFT <triqs_dft_tools.sumk_dft.SumkDFT>` and ``SumkDFTTools``
(``dft_data``, ``misc_data``, ``symmcorr_data``, …). Both default to the names
used below, and this is what all of
`DFTTools <https://triqs.github.io/dft_tools>`_ assumes.
:py:class:`SumkDFT <triqs_dft_tools.sumk_dft.SumkDFT>` only ever needs
``dft_input`` (plus ``dft_symmcorr_input`` if ``symm_op = 1``); all other groups
serve specific post-processing or charge self-consistency tasks.

.. _hdf_groups_overview:

Overview of the groups
----------------------

.. list-table::
   :header-rows: 1
   :widths: 22 12 30 36

   * - Group
     - Written by
     - Contents
     - Consumed by
   * - :ref:`dft_input <hdf_dft_input>`
     - all
     - The main data: shells, projectors, Hamiltonian, k-mesh. Always required.
     - :py:class:`SumkDFT <triqs_dft_tools.sumk_dft.SumkDFT>`
   * - :ref:`dft_misc_input <hdf_dft_misc_input>`
     - all but H(k)
     - Everything that is code specific and not part of the standard: Fermi
       weights, band windows, lattice information, k-paths.
     - charge self-consistency, transport, band plotting
   * - :ref:`dft_symmcorr_input <hdf_dft_symm_input>`
     - Wien2k, VASP, Elk
     - Symmetry operations acting on the **correlated** shells.
     - :py:class:`Symmetry <triqs_dft_tools.symmetry.Symmetry>` via
       ``SumkDFT.symmcorr``
   * - :ref:`dft_symmpar_input <hdf_dft_symm_input>`
     - Wien2k
     - Symmetry operations acting on **all** shells (partial projectors).
     - ``SumkDFTTools.symmpar``
   * - :ref:`dft_parproj_input <hdf_dft_parproj_input>`
     - Wien2k
     - Partial projectors onto all (also uncorrelated) shells.
     - ``SumkDFTTools`` DOS / partial charges
   * - :ref:`dft_bands_input <hdf_dft_bands_input>`
     - Wien2k, Elk
     - Projectors and eigenvalues on a k-path, for spectral functions.
     - ``SumkDFTTools.spaghettis``
   * - :ref:`dft_transp_input <hdf_dft_transp_input>`
     - Wien2k, Elk
     - Velocity (momentum) matrix elements for optics / transport.
     - ``transport_distribution()``
   * - :ref:`dft_contours_input <hdf_dft_contours_input>`
     - Elk
     - Projectors and eigenvalues on a k-plane, for Fermi surface contours.
     - ``SumkDFTTools.spectral_contours``
   * - :ref:`dft_update <hdf_dft_update>`
     - QE, and as an unused copy VASP
     - The DMFT charge density correction handed back to Quantum Espresso.
     - Quantum Espresso, in charge self-consistent runs. For VASP the correction
       goes to :file:`vaspgamma.h5` instead and this group is ignored.

Throughout this page the converters are abbreviated as follows:

.. list-table::
   :header-rows: 1
   :widths: 12 44 14

   * - Short name
     - Converter
     - ``dft_code``
   * - Wien2k
     - :py:class:`triqs_dftkit.wien2k.Converter`
     - ``'wien2k'``
   * - VASP
     - :py:class:`triqs_dftkit.vasp.Converter` (PLOVasp)
     - ``'vasp'``
   * - W90
     - :py:class:`triqs_dftkit.wannier90.Converter`
     - ``'w90'``
   * - Elk
     - :py:class:`triqs_dftkit.elk.Converter`
     - ``'elk'``
   * - H(k)
     - :py:class:`triqs_dftkit.hk.Converter`
     - ``'hk'``

.. _hdf_conventions:

Storage conventions
-------------------

The archive is written through ``HDFArchive``, which means the raw hdf5 file
looks slightly different from the Python objects:

* A **Python list** becomes a subgroup whose members are named ``'0'``, ``'1'``,
  … (with a ``Format = List`` attribute). So ``rot_mat`` is not one array but a
  group with one array per correlated shell, and ``dim_reps`` is a group of
  groups.
* A **dict** becomes a subgroup with one member per key (attribute
  ``Format = Dict``). Each entry of ``shells`` / ``corr_shells`` is such a group.
* A **complex** array is stored as a real array with one extra trailing axis of
  length 2 (real, imaginary part), flagged by the ``__complex__`` attribute.
  Reading through ``HDFArchive`` hides this; inspecting with ``h5py`` does not.
  A ``hopping`` of logical shape ``[n_k, n_spin_blocs, nb, nb]`` therefore shows
  up as ``[n_k, n_spin_blocs, nb, nb, 2]``.

Further conventions:

* The **spin index** runs over ``n_spin_blocs = SP + 1 - SO`` blocks: one block
  for a paramagnetic or a spin-orbit coupled calculation, two blocks for a
  collinear spin-polarised one. Below this is written as ``SP+1-SO``.
  (Exception: the Elk ``band_window``, see :ref:`there <hdf_dft_misc_input>`.)
* Arrays over bands are allocated with the **maximum** number of bands
  ``max(n_orbitals)`` occurring on the k-mesh and zero-padded at k-points with
  fewer bands, so that all matrices have the same size.
* ``band_window`` entries are **1-based** Kohn-Sham band indices, as used by the
  DFT codes, not 0-based Python indices.
* All energies are converted to **eV** by the converters, so ``energy_unit`` is
  the factor that was applied: ``13.605698`` for Wien2k (Ry),
  ``27.2113850560`` for Elk (Ha), and ``1.0`` for VASP, W90 and H(k), which
  already work in eV.

.. _hdf_dft_input:

Main data: ``dft_input``
------------------------

The main group of the archive, written by every converter (default name
``dft_input``, argument ``dft_subgrp``). It holds a lot of information, so that
DMFT calculations for all kinds of situations are possible, e.g. d-p
Hamiltonians, more than one correlated atomic shell, or symmetry operations for
the k-summation.

Required entries
^^^^^^^^^^^^^^^^

These are written by **all** converters and are what
:py:class:`SumkDFT <triqs_dft_tools.sumk_dft.SumkDFT>` reads.

.. list-table::
   :header-rows: 1
   :widths: 17 28 55

   * - Name
     - Type
     - Meaning
   * - ``energy_unit``
     - float
     - Unit of energy that was used to convert the DFT output to eV, see
       :ref:`conventions <hdf_conventions>`.
   * - ``dft_code``
     - string
     - DFT code the input data comes from, one of ``'wien2k'``, ``'vasp'``,
       ``'w90'``, ``'elk'``, ``'hk'``. Used by DFTTools where the codes have to
       be treated differently.
   * - ``n_k``
     - int
     - Number of k-points used for the BZ integration.
   * - ``k_dep_projection``
     - int
     - 1 if the dimension of the projection operators depends on the k-point, 0
       otherwise. 1 for Wien2k, VASP and Elk (projection onto a band window);
       0 for W90 and H(k) (Wannier basis, same size at every k).
   * - ``SP``
     - int
     - 1 for a spin-polarised Hamiltonian, 0 for a paramagnetic one.
   * - ``SO``
     - int
     - 1 if spin-orbit interaction is included, 0 otherwise.
   * - ``charge_below``
     - float
     - Number of electrons in the crystal below the correlated orbitals. Set to
       0 where the concept does not apply (VASP, W90), since the projection
       window is defined by band indices there.
   * - ``density_required``
     - float
     - Required total electron density. Needed to determine the chemical
       potential. The density in the projection window is then
       ``density_required - charge_below``.
   * - ``symm_op``
     - int
     - 1 if symmetry operations are used for the BZ sums, 0 if all k-points are
       directly included in the input. 1 for Wien2k and Elk, 0 for VASP, W90 and
       H(k). If 1, the group ``dft_symmcorr_input`` must be present.
   * - ``n_shells``
     - int
     - Number of atomic shells for which post-processing is possible. Note: this
       is *not* the number of correlated orbitals! If there are two equivalent
       atoms in the unit cell, ``n_shells`` is 2.
   * - ``shells``
     - list of dict {string:int}, dim ``n_shells``
     - Atomic shell information. For each shell a dict with at least the keys
       ``['atom', 'sort', 'l', 'dim']``. ``'atom'`` is the atom index, ``'sort'``
       defines the equivalency of the atoms, ``'l'`` is the angular quantum
       number, ``'dim'`` is the dimension of the atomic shell. E.g. for two
       equivalent atoms in the unit cell, ``atom`` runs from 0 to 1, but ``sort``
       can take only one value 0. Some converters add extra keys, see
       :ref:`below <hdf_shell_keys>`.
   * - ``n_corr_shells``
     - int
     - Number of correlated atomic shells. If there are two correlated
       equivalent atoms in the unit cell, ``n_corr_shells`` is 2.
   * - ``corr_shells``
     - list of dict {string:int}, dim ``n_corr_shells``
     - Correlated orbital information. For each correlated shell a dict with at
       least the keys ``['atom', 'sort', 'l', 'dim', 'SO', 'irep']``. The first
       four are as for ``shells``; ``'SO'`` is 1 if spin-orbit is included, 0
       otherwise; ``'irep'`` is the index of the irreducible representation
       (a dummy 0 for all converters except Wien2k). See
       :ref:`below <hdf_shell_keys>`.
   * - ``n_inequiv_shells``
     - int
     - Number of inequivalent atomic shells. Needs to be smaller than or equal
       to ``n_corr_shells``. The up / downfolding routines mediate between all
       correlated shells and the actual inequivalent shells, by using the
       self-energy etc. for all equal shells belonging to the same class of
       inequivalent shells. The mapping is performed with the information stored
       in ``corr_to_inequiv`` and ``inequiv_to_corr``.
   * - ``corr_to_inequiv``
     - list of int, dim ``n_corr_shells``
     - Mapping from correlated shells to inequivalent correlated shells. A list
       of length ``n_corr_shells`` containing integers, where same numbers mark
       equivalent sites.
   * - ``inequiv_to_corr``
     - list of int, dim ``n_inequiv_shells``
     - A list of length ``n_inequiv_shells`` containing list indices as integers
       pointing to the corresponding sites in ``corr_to_inequiv``.
   * - ``use_rotations``
     - int
     - 1 if local and global coordinate systems are used, 0 otherwise.
   * - ``rot_mat``
     - list of numpy.array.complex, dim ``n_corr_shells`` x
       [``corr_shells['dim']``, ``corr_shells['dim']``]
     - Rotation matrices for correlated shells, if ``use_rotations``. These
       rotations are automatically applied for up / downfolding. Set to the
       unity matrix if no rotations are used.
   * - ``rot_mat_time_inv``
     - list of int, dim ``n_corr_shells``
     - If ``SP`` is 1, 1 if the coordinate transformation contains inversion, 0
       otherwise. If ``use_rotations`` or ``SP`` is 0, a list of zeros.
   * - ``n_reps``
     - list of int, dim ``n_inequiv_shells``
     - Number of irreducible representations of each inequivalent correlated
       shell, e.g. 2 if eg/t2g splitting is used. Only Wien2k and H(k) use a
       value other than 1.
   * - ``dim_reps``
     - list of list of int, dim ``n_inequiv_shells`` x ``n_reps[ish]``
     - Dimension of the representations of each inequivalent correlated shell,
       e.g. ``[[2, 3]]`` for one shell split into eg/t2g. VASP and Elk write a
       single entry per shell (the full shell dimension); W90 writes a dummy
       scalar ``0`` per shell instead of a list.
   * - ``T``
     - list of numpy.array.complex, dim ``n_inequiv_shells`` x
       [``(2*l+1)*(SO+1)``, ``(2*l+1)*(SO+1)``]
     - Transformation matrix from the complex spherical harmonics to the
       impurity problem basis (normally the real cubic harmonics). Note that the
       dimension is the *full* :math:`2l+1` shell size (doubled with spin-orbit),
       independent of how many orbitals ``corr_shells['dim']`` actually holds:
       a t2g-only shell still gets a 5x5 matrix. This matrix can be used to
       calculate the 4-index U matrix, which is not done automatically. Only
       Wien2k and Elk fill it with real data; VASP writes the identity, H(k) a
       hardcoded d-shell matrix, and W90 writes zeros.
   * - ``n_orbitals``
     - numpy.array.int, dim [``n_k``, ``SP+1-SO``]
     - Number of Bloch bands included in the projection window for each
       k-point. If ``SP+1-SO = 2``, the number of included bands may depend on
       the spin projection up/down. For W90 and H(k) this is k-independent and
       equal to the number of Wannier functions.
   * - ``proj_mat``
     - numpy.array.complex, dim [``n_k``, ``SP+1-SO``, ``n_corr_shells``,
       ``max(corr_shells['dim'])``, ``max(n_orbitals)``]
     - Projection matrices from Bloch bands to Wannier orbitals. For efficient
       storage reasons, all matrices must be of the same size (given by the last
       two indices). For k-points with fewer bands, only the first entries are
       used, the rest are zero, e.g. if the number of Bloch bands ranges from
       4-6, all matrices are of size 6.
   * - ``bz_weights``
     - numpy.array.float, dim [``n_k``]
     - Weights of the k-points for the k summation. Deprecated duplicate of
       ``kpt_weights``, kept for backwards compatibility; where both are
       present they hold identical data.
   * - ``hopping``
     - numpy.array.complex, dim [``n_k``, ``SP+1-SO``, ``max(n_orbitals)``,
       ``max(n_orbitals)``]
     - Non-interacting Hamiltonian matrix for each k-point. As for ``proj_mat``,
       all matrices have to be of the same size. In the band-basis case
       (Wien2k, VASP in ``proj`` mode, Elk, W90 in Bloch basis) this is
       diagonal, holding the Kohn-Sham eigenvalues; in the Wannier-basis case
       (VASP in ``hk`` mode, W90 in Wannier basis, H(k)) it is the full
       Hamiltonian in orbital basis.

.. _hdf_shell_keys:

.. note::

   The exact set of keys in the ``shells`` / ``corr_shells`` dicts is not
   uniform across converters. Wien2k, W90 and H(k) write the four documented
   keys in ``shells``; VASP additionally writes ``'SO'`` and ``'irep'`` there
   (the ``shells`` and ``corr_shells`` entries are the same dicts); Elk adds a
   ``'natom'`` key to both. Consumers should therefore access keys by name and
   not rely on the number of entries.

Optional and converter specific entries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These live in ``dft_input`` as well, but only some converters write them.

.. list-table::
   :header-rows: 1
   :widths: 15 24 12 49

   * - Name
     - Type
     - Written by
     - Meaning
   * - ``kpts``
     - numpy.array.float, dim [``n_k``, 3]
     - VASP, W90
     - k-points given in fractional (reciprocal) coordinates.
   * - ``kpt_weights``
     - numpy.array.float, dim [``n_k``]
     - VASP, W90
     - Weights of the k-points for the k summation. The successor of
       ``bz_weights``.
   * - ``kpt_basis``
     - numpy.array.float, dim [3, 3]
     - VASP, W90 (Bloch basis)
     - Basis of the k-point mesh, i.e. the reciprocal lattice vectors as
       columns.
   * - ``proj_or_hk``
     - string
     - VASP
     - Switch determining whether the VASP converter ran in projection mode
       ``'proj'`` or in Hamiltonian mode ``'hk'``. In Hamiltonian mode
       ``hopping`` is written in orbital basis, in projection mode in band
       basis.
   * - ``proj_mat_csc``
     - numpy.array.complex, dim [``n_k``, ``SP+1-SO``,
       ``sum(shells['dim'])``, ``max(n_orbitals)``]
     - VASP (``hk`` mode)
     - Projection matrices from Bloch bands to Wannier orbitals for the
       Hamiltonian based ``hk`` approach. No site index is given, since h(k) is
       written in orbital basis; the orbital index runs over *all* shells
       (correlated and uncorrelated) of the projector group.
   * - ``wan_centres``
     - numpy.array.float, dim [``SP+1-SO``, ``n_wannier``, 3]
     - W90
     - Wannier centres in cartesian coordinates (Angstrom), read from
       :file:`seedname_centres.xyz`. Only written if that file exists, which
       requires ``write_xyz = true`` and ``translate_home_cell = false`` in the
       :program:`wannier90` input.

.. _hdf_dft_misc_input:

Code specific data: ``dft_misc_input``
--------------------------------------

Everything that a given DFT code provides on top of the standard goes into
``dft_misc_input`` (argument ``misc_subgrp``). None of it is needed for a
one-shot DMFT calculation, but charge self-consistency, transport and band
plotting rely on parts of it. The H(k) converter writes no such group.

.. list-table::
   :header-rows: 1
   :widths: 17 25 12 46

   * - Name
     - Type
     - Written by
     - Meaning
   * - ``band_window``
     - list of numpy.array.int, dim ``SP+1-SO`` x [``n_k``, 2]
     - Wien2k, VASP, W90 (Bloch basis), Elk
     - Lowest and highest Kohn-Sham band index (1-based) inside the projection
       window, per spin channel and k-point. Needed to write the charge density
       correction back to the DFT code. W90 stores this as a single
       numpy.array.int of dim [``SP+1-SO``, ``n_k``, 2] rather than a list, and
       Elk writes ``SP+1`` instead of ``SP+1-SO`` channels (it works with
       spinors, so the spin-orbit case still has two entries).
   * - ``dft_fermi_weights``
     - numpy.array, dim [``n_k``, ``SP+1-SO``, ``max(n_orbitals)``]
     - VASP, W90 (Bloch basis)
     - DFT Fermi weights (occupations) of the Kohn-Sham states for each
       k-point, used to compute the density matrix correction. Complex for
       VASP, real for W90.
   * - ``dft_fermi_energy``
     - float
     - W90
     - Kohn-Sham Fermi energy in eV, as read from the DFT output.
   * - ``kpts_cart``
     - numpy.array.float, dim [``n_k``, 3]
     - VASP, W90 (Bloch basis)
     - The k-points in cartesian coordinates.
   * - ``n_k_ibz``
     - int
     - VASP, W90 (Bloch basis)
     - Number of k-points in the irreducible wedge, when the DFT mesh has a
       reducible/irreducible split. The first ``n_k_ibz`` k-points of the mesh
       are the irreducible ones, and the density correction is written for
       those only. W90 writes it only if ``n_k_ibz < n_k``.
   * - ``symm_kpath``
     - subgroup with ``kpts`` (numpy.array.float, dim [``n_kpath``, 3]),
       ``labels`` (string), ``label_idx`` (numpy.array.int)
     - W90
     - High-symmetry k-path read from :file:`seedname_band.kpt`, for band
       plotting. ``labels`` is one character per high-symmetry point (e.g.
       ``'GXMG'``) and ``label_idx`` gives their positions along ``kpts``.
   * - ``lattice_type``
     - string
     - Wien2k
     - Bravais lattice type as defined by Wien2k.
   * - ``lattice_constants``
     - numpy.array.float, dim [3]
     - Wien2k
     - Unit cell parameters in atomic units.
   * - ``lattice_angles``
     - numpy.array.float, dim [3]
     - Wien2k
     - Unit cell angles in rad.
   * - ``n_symmetries``
     - int
     - Wien2k, Elk (transport)
     - Number of space group symmetry operations.
   * - ``rot_symmetries``
     - list of numpy.array.float, dim ``n_symmetries`` x [3, 3]
     - Wien2k, Elk (transport)
     - Matrix representation of all space group symmetry operations, used to
       symmetrise the transport quantities.
   * - ``cell_vol``
     - float
     - Elk (transport)
     - Unit cell volume, needed to normalise the transport distribution.
   * - ``nstsv``
     - int
     - Elk
     - Total number of (second-variational) states in the Elk calculation, i.e.
       the number of bands before the projection window is applied.
   * - ``vkl``
     - numpy.array.float, dim [``n_k``, 3]
     - Elk
     - The k-points in fractional coordinates (Elk's naming).

.. _hdf_dft_symm_input:

Symmetry operations: ``dft_symmcorr_input`` and ``dft_symmpar_input``
---------------------------------------------------------------------

These groups hold the data needed to apply symmetry operations in the DMFT loop
and are read by :py:class:`Symmetry <triqs_dft_tools.symmetry.Symmetry>`.
``dft_symmcorr_input`` (argument ``symmcorr_subgrp``) contains the symmetries
acting on the **correlated** shells and is required whenever ``symm_op = 1``;
``dft_symmpar_input`` (argument ``symmpar_subgrp``) contains the same data for
**all** shells and is only used together with the partial projectors. Both
groups have identical structure, they only differ in what ``orbits`` runs over.

.. list-table::
   :header-rows: 1
   :widths: 15 32 53

   * - Name
     - Type
     - Meaning
   * - ``n_symm``
     - int
     - Number of symmetry operations.
   * - ``n_atoms``
     - int
     - Number of atoms involved in the permutations.
   * - ``orbits``
     - list of dict {string:int}, dim ``n_orbits``
     - The shells the symmetry operations act on: a copy of ``corr_shells``
       (``dft_symmcorr_input``) or of ``shells`` (``dft_symmpar_input``).
   * - ``perm``
     - list of list of int, dim ``n_symm`` x ``n_atoms``
     - For each symmetry operation, the permutation of the atoms it induces
       (1-based atom indices).
   * - ``time_inv``
     - list of int, dim ``n_symm``
     - 1 if the symmetry operation has to be combined with time inversion, 0
       otherwise. All zeros unless ``SP = 1``.
   * - ``mat``
     - list of list of numpy.array.complex, dim ``n_symm`` x ``n_orbits`` x
       [``orbits['dim']``, ``orbits['dim']``]
     - Representation of each symmetry operation in the orbital basis of each
       shell.
   * - ``mat_tinv``
     - list of numpy.array.complex, dim ``n_orbits`` x [``orbits['dim']``,
       ``orbits['dim']``]
     - Additional matrix for the time inversion operation, per shell. Only read
       from file for ``SO = 0`` and ``SP = 0``, the identity otherwise.
   * - ``SP``
     - int
     - Spin-polarisation flag, mirroring ``dft_input/SP``.
   * - ``SO``
     - int
     - Spin-orbit flag, mirroring ``dft_input/SO``.

.. note::

   The VASP converter writes this group as a **placeholder** only: the
   symmetries are not analysed, so ``n_symm = 0``, ``n_atoms = 1``, and ``mat``
   / ``mat_tinv`` are 1x1 identities. This is consistent with VASP setting
   ``symm_op = 0``, i.e. the full BZ is written out and no symmetrisation is
   applied in the k-sum. The Elk converter writes real symmetry data but no
   time inversion (``time_inv`` is all zeros). The H(k) converter writes no
   symmetry group at all.

.. _hdf_dft_parproj_input:

Partial projectors: ``dft_parproj_input``
-----------------------------------------

Written by :py:meth:`Converter.convert_parproj_input()
<triqs_dftkit.wien2k.converter.Converter.convert_parproj_input>` of the Wien2k
converter (argument ``parproj_subgrp``), together with ``dft_symmpar_input``. It provides
projectors onto **all** shells, correlated and uncorrelated, which is what
``SumkDFTTools`` needs for partial densities of states and partial charges.
Only Wien2k supports this.

.. list-table::
   :header-rows: 1
   :widths: 22 34 44

   * - Name
     - Type
     - Meaning
   * - ``n_parproj``
     - numpy.array.int, dim [``n_shells``]
     - Number of partial projectors (theta projectors) per shell.
   * - ``proj_mat_all``
     - numpy.array.complex, dim [``n_k``, ``SP+1-SO``, ``n_shells``,
       ``max(n_parproj)``, ``max(shells['dim'])``, ``max(n_orbitals)``]
     - Projection matrices from Bloch bands onto all shells.
   * - ``dens_mat_below``
     - list of list of numpy.array.complex, dim ``SP+1-SO`` x ``n_shells`` x
       [``shells['dim']``, ``shells['dim']``]
     - Density matrix of each shell from the bands *below* the projection
       window. Divided by 2 if ``SP = 0``.
   * - ``rot_mat_all``
     - list of numpy.array.complex, dim ``n_shells`` x [``shells['dim']``,
       ``shells['dim']``]
     - Global-to-local rotation matrices for all shells, the analogue of
       ``rot_mat``.
   * - ``rot_mat_all_time_inv``
     - list of int, dim ``n_shells``
     - Whether the coordinate transformation of each shell contains inversion.
       Only read for ``SP = 1``, zeros otherwise.

.. _hdf_dft_bands_input:

Data on a k-path: ``dft_bands_input``
-------------------------------------

Written by ``convert_bands_input()`` of the Wien2k and Elk converters (argument
``bands_subgrp``). It repeats the projectors and eigenvalues, but on a
high-symmetry k-path instead of the BZ mesh, so that ``SumkDFTTools`` can
compute correlated band structures. The quantities have the same meaning as
their counterparts in ``dft_input`` / ``dft_parproj_input``, with ``n_k`` now
counting the k-points along the path.

.. list-table::
   :header-rows: 1
   :widths: 22 34 44

   * - Name
     - Type
     - Meaning
   * - ``n_k``
     - int
     - Number of k-points along the path.
   * - ``n_orbitals``
     - numpy.array.int, dim [``n_k``, ``SP+1-SO``]
     - Number of bands in the projection window at each k-point of the path.
   * - ``proj_mat``
     - numpy.array.complex, dim [``n_k``, ``SP+1-SO``, ``n_corr_shells``,
       ``max(corr_shells['dim'])``, ``max(n_orbitals)``]
     - Projectors onto the correlated shells along the path.
   * - ``hopping``
     - numpy.array.complex, dim [``n_k``, ``SP+1-SO``, ``max(n_orbitals)``,
       ``max(n_orbitals)``]
     - Kohn-Sham eigenvalues along the path, on the diagonal.
   * - ``n_parproj``
     - numpy.array.int, dim [``n_shells``]
     - As in ``dft_parproj_input``. Elk does not generate partial projectors
       and writes a dummy ``array([0])``.
   * - ``proj_mat_all``
     - numpy.array.complex, dim [``n_k``, ``SP+1-SO``, ``n_shells``,
       ``max(n_parproj)``, ``max(shells['dim'])``, ``max(n_orbitals)``]
     - As in ``dft_parproj_input``, along the path. Elk writes a dummy
       ``array([0])``.
   * - ``kpts_labels``
     - list of string
     - Names of the high-symmetry points of the path (e.g. ``'GAMMA'``, ``'X'``),
       for labelling the ticks of a band plot. Only written by Wien2k and VASP,
       and only when the underlying DFT output provides the labels.
   * - ``kpts_labels_idx``
     - numpy.array.int, dim [``len(kpts_labels)``]
     - Position of each entry of ``kpts_labels`` along the path, as a 0-based
       index into the ``n_k`` k-points. Each converter validates these indices
       as far as its own DFT output allows, so consumers should not rely on
       them being sorted or in range without checking.

.. note::

   The VASP converter does not currently write a ``dft_bands_input`` group; band
   structures from VASP projectors are obtained by running the converter on a
   separate band-path calculation.

.. _hdf_dft_transp_input:

Transport data: ``dft_transp_input``
------------------------------------

Written by ``convert_transport_input()`` of the Wien2k and Elk converters
(argument ``transp_subgrp``), and read by DFTTools'
``transport_distribution()``. Note that the accompanying
``n_symmetries``, ``rot_symmetries`` and (Elk) ``cell_vol`` go to
``dft_misc_input``, not here.

.. list-table::
   :header-rows: 1
   :widths: 22 34 44

   * - Name
     - Type
     - Meaning
   * - ``band_window_optics``
     - list of numpy.array.int, dim ``SP+1-SO`` x [``n_k``, 2]
     - Band window for which the velocity matrix elements are available. Same
       idea as ``band_window``, but rearranged to match the optics output of
       the DFT code, and in general a larger window.
   * - ``velocities_k``
     - list of list of numpy.array.complex, dim ``SP+1-SO`` x ``n_k`` x
       [``nb_k``, ``nb_k``, 3]
     - Velocity (momentum) matrix elements between all bands inside
       ``band_window_optics``, for the three cartesian directions. ``nb_k``
       varies from k-point to k-point, hence the nested-list storage.

.. _hdf_dft_contours_input:

Fermi surface contours: ``dft_contours_input``
----------------------------------------------

Written by :py:meth:`Converter.convert_contours_input()
<triqs_dftkit.elk.converter.Converter.convert_contours_input>` of the Elk
converter (argument ``cont_subgrp``), for spectral functions on a plane or a
k-grid in the BZ (Fermi surfaces).

.. list-table::
   :header-rows: 1
   :widths: 22 34 44

   * - Name
     - Type
     - Meaning
   * - ``n_k``
     - int
     - Number of irreducible k-points on the plane / grid.
   * - ``n_orbitals``
     - numpy.array.int, dim [``n_k``, ``SP+1-SO``]
     - Number of bands in the projection window at each of these k-points.
   * - ``proj_mat``
     - numpy.array.complex, dim [``n_k``, ``SP+1-SO``, ``n_corr_shells``,
       ``max(corr_shells['dim'])``, ``max(n_orbitals)``]
     - Projectors onto the correlated shells at these k-points.
   * - ``hopping``
     - numpy.array.complex, dim [``n_k``, ``SP+1-SO``, ``max(n_orbitals)``,
       ``max(n_orbitals)``]
     - Kohn-Sham eigenvalues at these k-points, on the diagonal.
   * - ``bmat``
     - numpy.array.float, dim [3, 3]
     - Reciprocal lattice vectors, used to convert the k-points to cartesian
       coordinates for plotting.
   * - ``BZ_n_k``
     - int
     - Number of k-points in the full (reducible) plane / grid.
   * - ``BZ_vkl``
     - numpy.array.float, dim [``BZ_n_k``, 3]
     - The full set of k-points in fractional coordinates.
   * - ``BZ_iknr``
     - numpy.array.int, dim [``BZ_n_k``]
     - For each of these k-points, the index of the irreducible k-point it maps
       onto, i.e. the mapping from the full plane into the ``n_k`` entries
       above.

.. _hdf_dft_update:

Charge density correction: ``dft_update``
-----------------------------------------

This is how the DMFT charge density correction is handed back to **Quantum
Espresso** during a charge self-consistent DFT+DMFT run: written either by
DFTTools' ``SumkDFT.calc_density_correction(dm_type='qe')`` or, equivalently, by
the ``calc_density_correction()`` step of the dftkit
:py:class:`Quantum Espresso driver <triqs_dftkit.qe.driver.Driver>`. It also
turns up in VASP archives, but is not used there, see the
:ref:`note below <hdf_dft_update_vasp>`.

.. list-table::
   :header-rows: 1
   :widths: 22 34 44

   * - Name
     - Type
     - Meaning
   * - ``delta_N``
     - numpy.array.complex, dim [``n_k``, ``max(n_orbitals)``,
       ``max(n_orbitals)``]
     - The DMFT correction to the Kohn-Sham density matrix in band basis, on
       the full k-mesh. For spin-orbit calculations this is the single spinor
       channel, otherwise the average of the two spin channels.

The other codes do **not** use this group: their corrections are written to
separate files that are not part of this archive.

* **VASP** reads :file:`vaspgamma.h5`, a standalone archive with a top-level
  ``band_window`` and a ``deltaN`` group holding one block per irreducible
  k-point (keyed ``up`` / ``down``, or ``ud`` for spin-orbit). This is what
  ``SumkDFT.calc_density_correction(dm_type='vasp')`` writes.
* **Wien2k** reads :file:`dens_mat.dat`, **Elk** reads :file:`DMATDMFT.OUT`
  (both plain text). For VASP a legacy text :file:`GAMMA` file can still be
  requested via ``dm_type='vasp', filename='GAMMA'``.
* **Abinit** reads a :file:`seedname.deltaN` text file, written by the dftkit
  Abinit driver.

.. _hdf_dft_update_vasp:

.. note::

   The dftkit VASP driver additionally leaves a copy of the correction in
   ``dft_update`` after writing :file:`vaspgamma.h5`, so the group can show up in
   a VASP archive. That copy is incidental: nothing reads it back, and it is not
   part of the interface to VASP.

General and simple H(k) converter
---------------------------------

The converters described above read the full output of a DFT code, including
information about symmetry operations, band windows and partial projectors.
Sometimes, however, we want to use a light implementation where the input
consists basically only of the Hamiltonian matrix in Wannier basis, given at a
grid of k-points in the first Brillouin zone. For this purpose the simple
:py:class:`triqs_dftkit.hk.Converter` is included in dftkit, implemented for
the simplest case of paramagnetic DFT calculations without spin-orbit coupling.
It reads a simple, easy to construct text file, and produces an archive that
can be used for the DMFT calculations. An example input file for a structure
with one correlated site with 3 t2g orbitals in the unit cell contains the
following:

  10               <- n_k

  1.0              <- density_required

  1                <- n_shells

  1 1 2 3          <- shells, as above: atom, sort, l, dim

  1                <- n_corr_shells

  1 1 2 3 0 0      <- corr_shells, as above: atom, sort, l, dim, SO, dummy

  2 2 3            <- n_reps, dim_reps (length 2, because eg/t2g splitting) for each inequivalent correlated shell

After this header, we give the Hamiltonian matrices for all the k-points. For
each k-point we give first the matrix of the real part, then the matrix of the
imaginary part. The projection matrices are set automatically to unity
matrices, no rotations, no symmetry operations are used. That means that the
symmetry sub group in the hdf5 archive needs not be set, since it is not used.
It is furthermore assumed that all k-points have equal weight in the k-sum.
Note that the input file should contain only the numbers, not the comments
given in above example.

The Hamiltonian matrices can be taken, e.g., from Wannier90, which constructs
the Hamiltonian in a maximally localized Wannier basis.

Note that with this simplified converter, no full charge self consistent
calculations are possible!
