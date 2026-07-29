"""
Pure mathematics for NICS probe point calculation.
No Qt, RDKit, or PyVista imports — safe to test headlessly.
"""

import math

import numpy as np

NICS1_HEIGHT = 1.0  # Angstroms above/below ring plane

#: A ring whose atoms sit this far (RMS, Å) off their own best-fit plane has no
#: well-defined "1 Å above the ring". The probe is still placed — cyclohexane
#: and tub-shaped COT are legitimate NICS subjects — but the value is reported
#: so callers can flag it. 0.1 Å passes every planar arene and catches the
#: cyclohexane chair (0.23) and the COT tub (0.38).
PLANARITY_TOLERANCE = 0.1


def ring_centroid(positions: np.ndarray) -> np.ndarray:
    """Geometric centroid of ring atom positions (shape N×3)."""
    return positions.mean(axis=0)


def orient_normal(normal: np.ndarray, reference: np.ndarray = None) -> np.ndarray:
    """Resolve the arbitrary sign of a plane normal.

    SVD returns an eigenvector, so its sign is whatever LAPACK produced — for
    naphthalene the two ring normals come back antiparallel, which would put
    ring 1's "above" probe on the opposite molecular face from ring 0's. Every
    ring must agree on which face is "above" or cross-ring comparisons of
    NICS(1) are meaningless on any molecule whose two faces differ.

    Aligns to *reference* when the two are not near-perpendicular; otherwise
    (rings at right angles, e.g. a spiro junction, where no shared "up" exists)
    falls back to a deterministic rule so repeat runs at least agree with
    each other.
    """
    if reference is not None:
        dot = float(normal @ reference)
        if abs(dot) > 1e-6:
            return normal if dot > 0 else -normal
    k = int(np.argmax(np.abs(normal)))
    return normal if normal[k] > 0 else -normal


def ring_normal(positions: np.ndarray, reference: np.ndarray = None) -> np.ndarray:
    """
    Best-fit plane normal via SVD of centred positions.
    Returns unit normal vector (shape (3,)), sign-resolved by `orient_normal`.
    """
    centered = positions - positions.mean(axis=0)
    _, _, vh = np.linalg.svd(centered)
    n = vh[-1]
    norm = np.linalg.norm(n)
    if norm <= 1e-10:
        return np.array([0.0, 0.0, 1.0])
    return orient_normal(n / norm, reference)


def planarity_rms(positions: np.ndarray, normal: np.ndarray = None) -> float:
    """RMS displacement of the ring atoms from their own best-fit plane, in Å."""
    if normal is None:
        normal = ring_normal(positions)
    oop = (positions - positions.mean(axis=0)) @ normal
    return float(np.sqrt((oop**2).mean()))


def compute_nics_points(
    positions: np.ndarray,
    height: float = NICS1_HEIGHT,
    reference: np.ndarray = None,
) -> dict:
    """
    Compute NICS probe positions for a ring.

    Parameters
    ----------
    positions : np.ndarray, shape (N, 3)
        3D coordinates of ring atoms (N >= 3).
    height : float
        Distance in Å for NICS(1) probes above/below ring plane.
    reference : np.ndarray, optional
        Molecule-level "up" direction used to make `nics1_above` mean the same
        face for every ring — see `orient_normal`.

    Returns
    -------
    dict with keys:
        'nics0'       – centroid (in-plane)
        'nics1_above' – centroid + height * normal
        'nics1_below' – centroid - height * normal
        'normal'      – the sign-resolved unit normal used
        'planarity'   – RMS out-of-plane deviation of the ring atoms (Å)
    """
    if len(positions) < 3:
        raise ValueError("Need at least 3 atoms to define a ring plane")
    centroid = ring_centroid(positions)
    normal = ring_normal(positions, reference)
    return {
        "nics0": centroid.copy(),
        "nics1_above": centroid + height * normal,
        "nics1_below": centroid - height * normal,
        "normal": normal,
        "planarity": planarity_rms(positions, normal),
    }


def molecular_reference_normal(mol) -> np.ndarray:
    """Best-fit plane normal over every ring atom in *mol*.

    Serves as the shared "up" direction for `orient_normal`. Using only ring
    atoms keeps a flexible side chain from tilting the reference away from the
    ring system it is supposed to describe.
    """
    ring_atoms = sorted({i for r in mol.GetRingInfo().AtomRings() for i in r})
    if len(ring_atoms) < 3:
        return np.array([0.0, 0.0, 1.0])
    pts = get_ring_positions(mol, tuple(ring_atoms))
    return ring_normal(pts)


def get_ring_positions(mol, ring_atoms: tuple) -> np.ndarray:
    """
    Extract 3D positions for a list of atom indices from an RDKit mol.
    Returns np.ndarray shape (N, 3).
    """
    conf = mol.GetConformer()
    return np.array(
        [[*conf.GetAtomPosition(i)] for i in ring_atoms],
        dtype=float,
    )


def get_rings(mol) -> list:
    """
    Return all rings in *mol* that have 3D coordinates.

    Each entry is a dict:
        {'atoms': tuple[int, ...], 'is_aromatic': bool}

    Returns [] if mol has no conformer.
    """
    try:
        mol.GetConformer()
    except Exception:
        return []

    arom_mol = _reperceived_aromaticity(mol)
    ring_info = mol.GetRingInfo()
    rings = []
    for atom_ring in ring_info.AtomRings():
        is_aromatic = all(arom_mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ring)
        rings.append({"atoms": tuple(atom_ring), "is_aromatic": is_aromatic})
    return rings


def _reperceived_aromaticity(mol):
    """Return a copy of *mol* with aromaticity freshly perceived.

    The flags already on the molecule reflect however it was last processed,
    not what it is. MMFF optimisation in particular kekulises in place, and
    RDKit does not restore the flags afterwards: a molecule that arrives here
    straight from the force field reports azulene as non-aromatic. Perceiving
    on a copy keeps the caller's molecule untouched.

    Falls back to the molecule as-is when perception fails (a partially built
    or deliberately odd valence state is not a reason to lose the ring table).
    """
    try:
        from rdkit import Chem
    except ImportError:  # pragma: no cover - rdkit absent only in unit stubs
        return mol
    try:
        copy = Chem.Mol(mol)
        Chem.SetAromaticity(copy, Chem.AromaticityModel.AROMATICITY_RDKIT)
        return copy
    except Exception:
        return mol


# ---------------------------------------------------------------------------
# 2D probe grids (NICS scan / isochemical-shielding-surface style sampling)
# ---------------------------------------------------------------------------

#: Ring-frame planes (u, v in the ring plane; n along the ring normal) followed
#: by fixed laboratory planes. Ring-frame grids follow the ring wherever it sits
#: and are what you want for a single ring; lab planes cut the whole molecule on
#: a fixed axis, which is easier to compare between molecules and to line up
#: with a crystallographic frame.
RING_GRID_PLANES = ("parallel", "perpendicular_u", "perpendicular_v")
LAB_GRID_PLANES = ("xy", "xz", "yz")
GRID_PLANES = RING_GRID_PLANES + LAB_GRID_PLANES

GRID_PLANE_LABELS = {
    "parallel": "Ring frame — parallel to ring plane (u–v)",
    "perpendicular_u": "Ring frame — perpendicular, along u (u–n)",
    "perpendicular_v": "Ring frame — perpendicular, along v (v–n)",
    "xy": "Lab frame — XY plane",
    "xz": "Lab frame — XZ plane",
    "yz": "Lab frame — YZ plane",
}

_LAB_AXES = {
    "xy": (np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])),
    "xz": (np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0])),
    "yz": (np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0])),
}


def heavy_atom_positions(mol) -> np.ndarray:
    """Coordinates of every non-hydrogen atom, shape (N, 3).

    Ghost atoms (atomic number 0) are excluded as well — sizing a grid from
    probes placed by an earlier run would make the grid grow every time it is
    regenerated.
    """
    conf = mol.GetConformer()
    pts = [
        [*conf.GetAtomPosition(a.GetIdx())]
        for a in mol.GetAtoms()
        if a.GetAtomicNum() > 1
    ]
    return np.array(pts, dtype=float) if pts else np.zeros((0, 3))


def molecule_bounds(mol) -> tuple:
    """Return (center, half_extents) of the heavy-atom bounding box, in Å."""
    pts = heavy_atom_positions(mol)
    if len(pts) == 0:
        return np.zeros(3), np.zeros(3)
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    return (lo + hi) / 2.0, (hi - lo) / 2.0


def axis_extents(
    mol,
    plane: str = "xy",
    positions: np.ndarray = None,
    reference: np.ndarray = None,
    center: np.ndarray = None,
    margin: float = 2.0,
    minimum: float = 2.0,
) -> tuple:
    """Half-widths (e1, e2, e3) that cover the molecule along the grid's own axes.

    Every heavy atom is projected onto each of the three grid axes and the
    largest |projection| taken, so the answer is exact for a tilted ring frame
    as well as for the lab planes — a bounding box in lab coordinates would
    otherwise either clip a tilted grid or, if padded by its diagonal,
    massively oversize it.

    *margin* defaults to 2 Å — roughly a van der Waals radius — so the
    ring-current field is sampled where it has decayed rather than clipped at
    the last atom.
    """
    ax1, ax2, ax3 = grid_axes(plane, positions, reference)
    pts = heavy_atom_positions(mol)
    if center is None:
        center = (
            molecule_bounds(mol)[0] if plane in _LAB_AXES else ring_centroid(positions)
        )
    if len(pts) == 0:
        return (minimum, minimum, minimum)
    rel = pts - np.asarray(center, dtype=float)
    return tuple(
        max(float(np.abs(rel @ axis).max()) + margin, minimum)
        for axis in (ax1, ax2, ax3)
    )


def counts_for_spacing(extents, spacing: float) -> tuple:
    """Point counts giving (at most) *spacing* Å between neighbours on each axis.

    Used by the uniform-spacing option: with independent half-widths, matching
    the counts is not the same as matching the step, and it is the step that
    decides whether the sampled field is smooth. A zero-width axis collapses to
    a single point rather than to a division by zero.
    """
    if spacing <= 0:
        raise ValueError("spacing must be positive")
    counts = []
    for e in np.atleast_1d(np.asarray(extents, dtype=float)):
        counts.append(1 if e <= 0 else int(math.ceil(2.0 * e / spacing)) + 1)
    return tuple(counts)


def ring_frame(positions: np.ndarray, reference: np.ndarray = None) -> tuple:
    """Orthonormal frame (centroid, u, v, n) for a ring.

    *n* is the sign-resolved ring normal; *u* points from the centroid toward
    the first ring atom (projected into the plane), so the grid keeps a fixed
    orientation relative to the ring's own atoms rather than to the lab axes;
    *v* completes a right-handed set.
    """
    centroid = ring_centroid(positions)
    n = ring_normal(positions, reference)
    u = positions[0] - centroid
    u = u - (u @ n) * n
    norm = np.linalg.norm(u)
    if norm < 1e-8:
        # First atom sits on the axis — impossible for a real ring, but pick
        # any in-plane direction rather than emitting NaNs.
        seed = np.array([1.0, 0.0, 0.0])
        if abs(n @ seed) > 0.9:
            seed = np.array([0.0, 1.0, 0.0])
        u = seed - (seed @ n) * n
        norm = np.linalg.norm(u)
    u = u / norm
    v = np.cross(n, u)
    return centroid, u, v, n


def grid_axes(
    plane: str, positions: np.ndarray = None, reference: np.ndarray = None
) -> tuple:
    """Return (ax1, ax2, ax3) unit vectors spanning *plane*, ax3 its normal.

    *positions* (ring atoms) is required only for the ring-frame planes.
    """
    if plane not in GRID_PLANES:
        raise ValueError(f"Unknown plane {plane!r}; expected one of {GRID_PLANES}")
    if plane in _LAB_AXES:
        a1, a2 = _LAB_AXES[plane]
        return a1, a2, np.cross(a1, a2)
    if positions is None:
        raise ValueError(f"plane {plane!r} needs ring positions")
    _c, u, v, n = ring_frame(positions, reference)
    if plane == "parallel":
        return u, v, n
    if plane == "perpendicular_u":
        return u, n, v
    return v, n, u


def _per_axis(value, count, name):
    """Broadcast a scalar to *count* axes, or validate a per-axis sequence.

    Lets every caller pass either one number for a square/cubic grid or one
    per axis for a rectangular one, without two sets of parameters.
    """
    if np.isscalar(value):
        return [value] * count
    values = list(value)
    if len(values) != count:
        raise ValueError(f"{name} needs {count} values, got {len(values)}")
    return values


def _axis_steps(n_points, extent, name):
    """Sample points along one axis, validating the pair."""
    if n_points < 1:
        raise ValueError(f"{name}: need at least 1 point")
    if extent < 0:
        raise ValueError(f"{name}: extent must not be negative")
    # One point means "just this plane/line"; it belongs at the centre, not at
    # -extent, or a single-layer grid would sit off the feature it describes.
    if n_points == 1:
        return np.array([0.0])
    return np.linspace(-extent, extent, n_points)


def compute_nics_grid(
    positions: np.ndarray = None,
    plane: str = "parallel",
    n_points=9,
    extent=3.0,
    offset: float = 0.0,
    reference: np.ndarray = None,
    center: np.ndarray = None,
) -> list:
    """Probe positions on a rectangular grid.

    Parameters
    ----------
    positions : np.ndarray, shape (N, 3), optional
        Ring atom coordinates. Required for the ring-frame planes, ignored for
        the lab planes.
    plane : one of GRID_PLANES
        'parallel'        – spans the ring plane itself, displaced along the
                            normal by *offset* (offset=1.0 is the NICS(1)
                            face map).
        'perpendicular_u' – contains the ring normal; the classic side-on
                            scan that shows the ring-current cone.
        'perpendicular_v' – as above, rotated 90° about the normal.
        'xy' / 'xz' / 'yz' – fixed laboratory planes.
    n_points : int or (int, int)
        Points along each in-plane axis. A scalar gives a square grid.
    extent : float or (float, float)
        Half-width in Å along each in-plane axis. A scalar gives a square grid.
        Molecules are rarely square, so the two are independent.
    offset : float
        Displacement of the whole plane along its own normal.
    center : np.ndarray, optional
        Grid centre. Defaults to the ring centroid for ring-frame planes and
        to the origin for lab planes — callers wanting a molecule-centred lab
        grid should pass `molecule_bounds(mol)[0]`.

    Returns
    -------
    list of dicts: {'pos': np.ndarray(3,), 'i': int, 'j': int, 'a': float, 'b': float}
    where (a, b) are the in-plane coordinates in Å and (i, j) the grid indices.
    """
    n1, n2 = _per_axis(n_points, 2, "n_points")
    e1, e2 = _per_axis(extent, 2, "extent")

    ax1, ax2, ax3 = grid_axes(plane, positions, reference)
    if center is None:
        center = np.zeros(3) if plane in _LAB_AXES else ring_centroid(positions)
    origin = np.asarray(center, dtype=float) + offset * ax3

    steps1 = _axis_steps(n1, e1, "axis 1")
    steps2 = _axis_steps(n2, e2, "axis 2")
    return [
        {
            "pos": origin + a * ax1 + b * ax2,
            "i": i,
            "j": j,
            "a": float(a),
            "b": float(b),
        }
        for i, a in enumerate(steps1)
        for j, b in enumerate(steps2)
    ]


def compute_nics_volume(
    positions: np.ndarray = None,
    plane: str = "parallel",
    n_points=9,
    extent=3.0,
    n_normal: int = None,
    normal_extent: float = None,
    offset: float = 0.0,
    reference: np.ndarray = None,
    center: np.ndarray = None,
) -> list:
    """Probe positions filling a 3D box — the sampling an ICSS plot needs.

    The box is spanned by the same frame `compute_nics_grid` uses: *plane*
    fixes the two in-plane axes and the third axis is that plane's normal.
    So a 'parallel' volume is a slab stacked along the ring normal, and an
    'xy' volume is an axis-aligned lab box.

    All three axes are sized independently, because the useful sampling is
    rarely cubic: the ring-current field decays over ~3-4 Å along the normal
    while the in-plane window usually has to span the whole ring system. Pass
    scalars for a cube, or a value per axis. *n_normal* / *normal_extent*
    override the third axis and default to the in-plane values.

    Point count is n1 * n2 * n_normal and grows fast — 15 x 15 x 15 is 3375
    ghost centres. Callers are expected to guard it.

    Returns
    -------
    list of dicts: {'pos', 'i', 'j', 'k', 'a', 'b', 'c'} where (a, b, c) are
    coordinates in Å along (ax1, ax2, normal).
    """
    n1, n2 = _per_axis(n_points, 2, "n_points")
    e1, e2 = _per_axis(extent, 2, "extent")
    if n_normal is None:
        n_normal = n2
    if normal_extent is None:
        normal_extent = e2

    ax1, ax2, ax3 = grid_axes(plane, positions, reference)
    if center is None:
        center = np.zeros(3) if plane in _LAB_AXES else ring_centroid(positions)
    origin = np.asarray(center, dtype=float) + offset * ax3

    steps1 = _axis_steps(n1, e1, "axis 1")
    steps2 = _axis_steps(n2, e2, "axis 2")
    normal_steps = _axis_steps(n_normal, normal_extent, "normal axis")

    return [
        {
            "pos": origin + a * ax1 + b * ax2 + c * ax3,
            "i": i,
            "j": j,
            "k": k,
            "a": float(a),
            "b": float(b),
            "c": float(c),
        }
        for i, a in enumerate(steps1)
        for j, b in enumerate(steps2)
        for k, c in enumerate(normal_steps)
    ]
