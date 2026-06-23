################################################################################
#
# TRIQS / dftkit
#
# Construction of correlated-shell symmetry operations for VASP, enabling DMFT
# to run on the irreducible Brillouin zone (IBZ) with symmetrization instead of
# unfolding to the full k-grid.
#
# VASP exposes, in vaspout.h5, the space-group operations that relate every
# full-grid k-point to its IBZ representative (as integer matrices in the basis
# of the reciprocal lattice). From these we build the rotation matrices acting
# on the correlated local orbitals (the `mat` entries of `dft_symmcorr_input`)
# that triqs_modest / dft_tools use to restore the full-BZ average of any
# k-summed local quantity.
#
################################################################################
"""
VASP symmetry tools for IBZ-based DMFT.
"""
import numpy as np


# ---------------------------------------------------------------------------
# Real cubic-harmonic rotation matrices D^l(R)
#
# The local orbitals are real (cubic) harmonics in the VASP/PLO ordering. Under
# a Cartesian point-group operation R (proper or improper, det = +-1) they
# transform with a (2l+1)x(2l+1) orthogonal matrix D^l(R), built here without
# any external dependency:
#
#   l = 1 : real harmonics are the Cartesian components; D^1 is R itself,
#           re-ordered to the VASP basis [y, z, x].
#   l = 2 : each real d harmonic maps to a symmetric traceless 3x3 tensor T_i;
#           under R the tensor transforms as T -> R T R^T, and
#           D^2[i, j] = <T_i, R T_j R^T> (Frobenius product on the orthonormal
#           tensor basis). Correct for improper R (d harmonics are parity-even).
#
# Basis orderings match modest's dft_tools/vasp/spherical_rotation.cpp:
#   l = 1 : "y", "z", "x"            (m = -1, 0, +1)
#   l = 2 : "xy","yz","z2","xz","x2-y2"
# ---------------------------------------------------------------------------

# l = 1: real-harmonic basis vectors in VASP order [y, z, x]
_P_BASIS = [np.array([0.0, 1.0, 0.0]),   # y
            np.array([0.0, 0.0, 1.0]),   # z
            np.array([1.0, 0.0, 0.0])]   # x


def _d_tensor_basis():
    """Symmetric traceless 3x3 tensors for the 5 real d harmonics (VASP order),
    Frobenius-orthonormal: <T_i, T_j> = delta_ij."""
    s2 = np.sqrt(2.0)
    T = []
    m = np.zeros((3, 3)); m[0, 1] = m[1, 0] = 1 / s2; T.append(m)          # xy
    m = np.zeros((3, 3)); m[1, 2] = m[2, 1] = 1 / s2; T.append(m)          # yz
    T.append(np.diag([-1.0, -1.0, 2.0]) / np.sqrt(6.0))                    # z^2
    m = np.zeros((3, 3)); m[0, 2] = m[2, 0] = 1 / s2; T.append(m)          # xz
    T.append(np.diag([1.0, -1.0, 0.0]) / s2)                               # x^2-y^2
    return T


_D_BASIS = _d_tensor_basis()


def real_harmonic_rotation(R, l):
    """
    Rotation matrix of the real (cubic) harmonics of angular momentum ``l``
    under the Cartesian orthogonal operation ``R`` (det = +-1).

    Parameters
    ----------
    R : (3, 3) array
        Cartesian point-group operation (proper or improper).
    l : int
        Angular momentum (1 or 2 are implemented).

    Returns
    -------
    (2l+1, 2l+1) ndarray
        Orthogonal rotation matrix in the VASP real-harmonic basis.
    """
    if l == 1:
        D = np.zeros((3, 3))
        for j, bj in enumerate(_P_BASIS):
            Rbj = R @ bj
            for i, bi in enumerate(_P_BASIS):
                D[i, j] = bi @ Rbj
        return D
    if l == 2:
        D = np.zeros((5, 5))
        for j, Tj in enumerate(_D_BASIS):
            RTjRt = R @ Tj @ R.T
            for i, Ti in enumerate(_D_BASIS):
                D[i, j] = np.tensordot(Ti, RTjRt)
        return D
    raise NotImplementedError(
        f"real_harmonic_rotation is only implemented for l = 1, 2 (got l = {l}). "
        "Higher l (f electrons) needs the Wigner-D / rank-l tensor extension.")


def _subspace_invariance(T, R_cart, l):
    """
    Maximum 'leakage' of the correlated subspace out of itself under the real-
    harmonic representation D(R).

    The correlated orbitals span the row space of ``T`` (shape ``(dim, 2l+1)``)
    inside the (2l+1)-dimensional real-harmonic space. The returned value is
    ``max_R || (1 - P) D(R) P ||_2`` with ``P`` the orthogonal projector onto
    that subspace. It is zero iff the subspace is invariant under every ``R``
    (the condition for ``Q = T D T^dag`` to be a group representation). The
    measure is gauge-free: it depends only on the span of ``T``, not on its
    particular rows, so a non-orthonormal or rotated transform gives the same
    answer.
    """
    nm = 2 * l + 1
    T = np.asarray(T, dtype=complex)
    # Orthonormal basis of the row span via SVD (rows of Vh with non-zero sv).
    _, s, Vh = np.linalg.svd(T, full_matrices=False)
    Q = Vh[s > 1e-10]                 # (r, nm), orthonormal rows
    P = Q.conj().T @ Q                # (nm, nm) Hermitian projector
    eye = np.eye(nm)
    viol = 0.0
    for R in R_cart:
        D = real_harmonic_rotation(R, l).astype(complex)
        viol = max(viol, np.linalg.norm((eye - P) @ D @ P, ord=2))
    return viol


# ---------------------------------------------------------------------------
# Reading the VASP symmetry data from vaspout.h5
# ---------------------------------------------------------------------------
def _close_group(ops, max_size=200):
    """
    Close a set of integer matrices under multiplication to recover the full
    group. ``ops`` is an (n, 3, 3) integer array; returns an (m, 3, 3) array
    with the complete group. Exact (integer arithmetic).
    """
    group = [np.asarray(o, dtype=int) for o in ops]
    seen = {o.tobytes() for o in group}
    changed = True
    while changed:
        changed = False
        for a in list(group):
            for b in list(group):
                c = a @ b
                key = c.tobytes()
                if key not in seen:
                    seen.add(key)
                    group.append(c)
                    changed = True
        if len(group) > max_size:
            raise RuntimeError("symmetry group closure exceeded %d operations; "
                               "stored symmetry data may be inconsistent" % max_size)
    return np.array(group)


def read_vasp_symmetry(vasp_h5):
    """
    Read the symmetry-relevant data from a vaspout.h5 file.

    Returns a dict with:
        n_k_ibz       : number of irreducible k-points
        ibz_weights   : (n_k_ibz,) IBZ k-point weights (sum to 1)
        R_cart        : list of unique Cartesian point-group operations (3x3)
        positions     : (natom, 3) fractional ion positions
        lattice       : (3, 3) Cartesian lattice vectors (rows)
        type_of_ion   : (natom,) integer ion-type index
    or ``None`` if the file does not contain k-point symmetry data (no symmetry
    was used / not a symmetry-reduced run).
    """
    from h5 import HDFArchive
    with HDFArchive(vasp_h5, 'r') as ar:
        if 'results' not in ar or 'electron_eigenvalues' not in ar['results']:
            return None
        ee = ar['results']['electron_eigenvalues']
        if 'kpoints_symmetry_symop' not in ee:
            return None
        symop = np.asarray(ee['kpoints_symmetry_symop'])      # (nktot, 3, 3) int, reciprocal-fractional
        ibz_weights = np.asarray(ee['kpoints_symmetry_weight'])  # (nkibz,) sum = 1
        n_k_ibz = int(ee['kpoints']) if 'kpoints' in ee else len(ibz_weights)
        positions_grp = ar['results']['positions']
        lattice = np.asarray(positions_grp['lattice_vectors'])    # rows = Cartesian lattice vectors
        positions = np.asarray(positions_grp['position_ions'])    # (natom, 3) fractional
        number_ion_types = np.asarray(positions_grp['number_ion_types'])

    # type index per ion
    type_of_ion = np.concatenate([[it] * n for it, n in enumerate(number_ion_types)])

    # reciprocal basis (columns = reciprocal vectors): B = 2pi (A^-1)^T, A rows = lattice vectors
    B = (2 * np.pi * np.linalg.inv(lattice).T).T
    Binv = np.linalg.inv(B)
    # VASP stores, per full-grid k-point, the single operation mapping the IBZ
    # representative to that point (a coset representative). The set of stored
    # operations is therefore complete only if some k-point has a full-size star.
    # On coarse / special meshes it is a proper subset, so we close it under
    # multiplication to recover the full point group (exact: integer matrices).
    uniq = _close_group(np.unique(symop.reshape(-1, 9), axis=0).reshape(-1, 3, 3))
    R_cart = [B @ s.astype(float) @ Binv for s in uniq]

    return dict(n_k_ibz=n_k_ibz, ibz_weights=ibz_weights, R_cart=R_cart,
                positions=positions, lattice=lattice, type_of_ion=type_of_ion)


# ---------------------------------------------------------------------------
# Atom permutations under the space-group operations
# ---------------------------------------------------------------------------
def _frac_mod(x):
    """Fold fractional coordinates into [0, 1)."""
    return x - np.floor(x + 1e-8)


def atom_permutation(R_cart, positions, lattice, type_of_ion, tol=1e-4):
    """
    Permutation of the crystal atoms under the Cartesian operation ``R_cart``.

    A global fractional translation ``t`` (space-group translation) is detected
    so the routine works for symmorphic operations regardless of cell origin.
    Returns a 1-indexed permutation array ``perm`` such that atom ``j`` maps to
    atom ``perm[j]`` (1-indexed), or ``None`` if no consistent permutation is
    found (e.g. a genuinely non-symmorphic operation needing its own glide/screw
    translation, not handled here).
    """
    natom = len(positions)
    pos_cart = positions @ lattice                       # (natom, 3) Cartesian
    rot_cart = pos_cart @ R_cart.T                        # apply R to each position
    Linv = np.linalg.inv(lattice)
    rot_frac = rot_cart @ Linv                            # back to fractional

    # candidate global translations from mapping atom 0 -> any same-type atom
    for i0 in range(natom):
        if type_of_ion[i0] != type_of_ion[0]:
            continue
        t = positions[i0] - rot_frac[0]
        perm = np.full(natom, -1, dtype=int)
        ok = True
        for j in range(natom):
            img = _frac_mod(rot_frac[j] + t)
            # find matching atom of the same type
            match = -1
            for i in range(natom):
                if type_of_ion[i] != type_of_ion[j]:
                    continue
                d = _frac_mod(positions[i] - img)
                d = np.minimum(d, 1.0 - d)
                if np.linalg.norm(d) < tol:
                    match = i
                    break
            if match < 0:
                ok = False
                break
            perm[j] = match + 1   # 1-indexed
        if ok and len(set(perm)) == natom:
            return perm
    return None


# ---------------------------------------------------------------------------
# Top-level: build the dft_symmcorr_input payload
# ---------------------------------------------------------------------------
def build_symmcorr(vasp_h5, corr_shells, transforms, SP, SO):
    """
    Build the ``dft_symmcorr_input`` data for the correlated shells.

    Parameters
    ----------
    vasp_h5 : str
        Path to vaspout.h5.
    corr_shells : list of dict
        Correlated shells (as built by the converter); each has 'atom', 'l',
        'dim'.
    transforms : list of (dim, 2l+1) arrays
        Per-correlated-shell transformation matrix T mapping the real harmonics
        onto the correlated orbitals (the PLO TRANSFORM). One per corr shell.
    SP, SO : int
        Spin-polarization and spin-orbit flags.

    Returns
    -------
    dict with keys n_symm, n_atoms, perm, orbits, SO, SP, time_inv, mat, mat_tinv
    and additionally n_k_ibz, ibz_weights for the IBZ slicing.
    Returns ``None`` if no symmetry data is available.
    """
    sym = read_vasp_symmetry(vasp_h5)
    if sym is None:
        return None

    R_cart = sym['R_cart']
    n_symm = len(R_cart)
    n_corr = len(corr_shells)

    # The projected operations Q = T D(R) T^dag form a representation of the
    # point group on the correlated orbitals only if the orbital subspace (the
    # row span of T) is invariant under every D(R). If it is not (e.g. an
    # arbitrary LOCPROJ subset that splits a degenerate multiplet such as t2g),
    # the symmetrization of any k-summed quantity on the IBZ is invalid: the
    # missing partners would be needed to close the star. Refuse rather than
    # silently produce a wrong full-BZ average.
    for ish, csh in enumerate(corr_shells):
        viol = _subspace_invariance(transforms[ish], R_cart, csh['l'])
        if viol > 1e-6:
            raise RuntimeError(
                f"IBZ symmetrization: the correlated subspace of shell {ish} "
                f"(atom {csh.get('atom')}, l={csh['l']}, dim={csh['dim']}) is not "
                f"invariant under the crystal point group (leakage {viol:.2e}). The "
                "chosen orbital set is not symmetry-closed, so k-summed quantities "
                "cannot be symmetrized on the irreducible BZ. Use a symmetry-adapted "
                "orbital set (the full 2l+1 shell, or a complete irrep such as t2g or "
                "eg in cubic symmetry), or run on the full grid (use_ibz=False).")

    # symmetrization matrices: mat[g] = list over corr shells of (dim x dim)
    mat = []
    for R in R_cart:
        per_shell = []
        for ish, csh in enumerate(corr_shells):
            l = csh['l']
            T = np.asarray(transforms[ish], dtype=complex)   # (dim, 2l+1)
            D = real_harmonic_rotation(R, l).astype(complex)  # (2l+1, 2l+1)
            Q = T @ D @ T.conj().T                            # (dim, dim) in the correlated basis
            per_shell.append(Q)
        mat.append(per_shell)

    # atom permutations (over all crystal atoms, 1-indexed)
    natom = len(sym['positions'])
    perm = []
    for R in R_cart:
        p = atom_permutation(R, sym['positions'], sym['lattice'], sym['type_of_ion'])
        if p is None:
            # fall back to identity; correct whenever the correlated atoms are
            # each left invariant by the point operations (warn and continue).
            import triqs.utility.mpi as mpi
            mpi.report("  WARNING: could not determine atom permutation for a symmetry op; "
                       "assuming identity (valid only if correlated atoms are fixed).")
            p = np.arange(1, natom + 1)
        perm.append(list(map(int, p)))

    time_inv = [0] * n_symm
    mat_tinv = [np.identity(csh['dim'], complex) for csh in corr_shells]

    return dict(n_symm=n_symm, n_atoms=natom, perm=perm, orbits=corr_shells,
                SO=SO, SP=SP, time_inv=time_inv, mat=mat, mat_tinv=mat_tinv,
                n_k_ibz=sym['n_k_ibz'], ibz_weights=sym['ibz_weights'])
