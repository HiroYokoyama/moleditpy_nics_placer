"""
Pure mathematics for NICS probe point calculation.
No Qt, RDKit, or PyVista imports — safe to test headlessly.
"""
import numpy as np

NICS1_HEIGHT = 1.0  # Angstroms above/below ring plane


def ring_centroid(positions: np.ndarray) -> np.ndarray:
    """Geometric centroid of ring atom positions (shape N×3)."""
    return positions.mean(axis=0)


def ring_normal(positions: np.ndarray) -> np.ndarray:
    """
    Best-fit plane normal via SVD of centred positions.
    Returns unit normal vector (shape (3,)).
    """
    centered = positions - positions.mean(axis=0)
    _, _, vh = np.linalg.svd(centered)
    n = vh[-1]
    norm = np.linalg.norm(n)
    return n / norm if norm > 1e-10 else np.array([0.0, 0.0, 1.0])


def compute_nics_points(positions: np.ndarray, height: float = NICS1_HEIGHT) -> dict:
    """
    Compute NICS probe positions for a ring.

    Parameters
    ----------
    positions : np.ndarray, shape (N, 3)
        3D coordinates of ring atoms (N >= 3).
    height : float
        Distance in Å for NICS(1) probes above/below ring plane.

    Returns
    -------
    dict with keys:
        'nics0'       – centroid (in-plane)
        'nics1_above' – centroid + height * normal
        'nics1_below' – centroid - height * normal
    """
    if len(positions) < 3:
        raise ValueError("Need at least 3 atoms to define a ring plane")
    centroid = ring_centroid(positions)
    normal = ring_normal(positions)
    return {
        "nics0": centroid.copy(),
        "nics1_above": centroid + height * normal,
        "nics1_below": centroid - height * normal,
    }


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

    ring_info = mol.GetRingInfo()
    rings = []
    for atom_ring in ring_info.AtomRings():
        is_aromatic = all(
            mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ring
        )
        rings.append({"atoms": tuple(atom_ring), "is_aromatic": is_aromatic})
    return rings
