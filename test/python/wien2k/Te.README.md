# Te projector fixture

`Te.almblmup.gz` and `Te.almblmdn.gz` are the Wien2k `lapw2 -almd` projector
inputs for trigonal tellurium at a single k-point. Te is non-centrosymmetric
(P3_121) with three equivalent atoms, so these files exercise the composed
local frame in `ctqmcout`, `parproj` and `outband`, which the centrosymmetric
fixtures never reach. They pair with the shipped `Te.struct`, `Te.dmftsym` and
`Te.indmftpr`, which are the files of the same calculation with trailing
whitespace removed.

## Producer provenance

- Wien2k 19.1 (gfortran build), serial, no MPI.
- Structure: `Te.struct`, P3_121, a = 9.41 bohr, c = 11.199 bohr, one
  inequivalent Te with RMT = 2.10 bohr; spin-orbit with M || z, which leaves
  the six operations of the space group; the three two-fold axes
  perpendicular to z reverse M and enter combined with time reversal.
- Density: PBE, spin-polarized scalar-relativistic self-consistent run on
  12 irreducible k-points, converged to 1e-4 Ry in energy and 1e-3 e in
  charge. The solution is non-magnetic (cell moment -0.001 μB).
- Projector step on that density, spin-polarized with spin-orbit coupling, at
  the single general k-point (3, 7, 11)/35 in reciprocal-lattice units with
  weight 1. `Te.klist`:

  ```
           1         3         7        11        35  1.0 -7.0  1.5         1 k, div: (  1  1  1)
  END
  ```

- Band range: 34 spin-orbit bands, 31 to 64 (Te 5s, 5p and the lowest
  conduction bands). The Te 4d semicore bands are left out with `EMIN = -1.0`
  in `Te.in2c` (NE stays 48, since the Fermi search counts all bands). The
  spin-orbit basis is limited to scalar-relativistic states below 0.85 Ry
  (`Te.inso` Emax; 1.9 Ry in the original projector run), which moves the
  energies of bands 31 to 64 by at most 6.3e-5 Ry.
- Fermi level: 0.26554 Ry from Fermi smearing (`TEMP 0.002` in `Te.in2c`),
  because the tetrahedron method needs a mesh. The full-mesh run gave
  0.24443 Ry. With the shipped `Te.indmftpr` (cubic p, -0.50 to 0.50 Ry) the
  window holds bands 37 to 62, as in the full-mesh run.

Commands, in a case directory named `Te`, after writing the klist and the
changed `Te.in2c` and `Te.inso` lines above:

```sh
x lapw1 -up; x lapw1 -dn; x lapwso -up
x lapw2 -up -c -so -almd; x lapw2 -dn -c -so -almd
gzip -9 -n Te.almblmup Te.almblmdn
```

## Native references

dmftproj derives the case name from the working directory, so run it in a
directory named `Te` that holds `Te.struct`, `Te.dmftsym`, `Te.indmftpr` and the
decompressed `Te.almblmup` and `Te.almblmdn`, with `$WIENROOT/SRC_templates`
containing `case.cf_p_cubic` from `fortran/dmftproj/SRC_templates`:

```sh
/path/to/build/fortran/dmftproj/dmftproj -sp -so
```

Without spin-orbit (indmftpr SO flag 0) run `dmftproj -sp`; for the
paramagnetic path, use `Te.almblmup` as `Te.almblm` and run `dmftproj` without
flags, as for the `CaOs2_noso` paramagnetic references.
