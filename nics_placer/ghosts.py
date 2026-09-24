"""
Ghost-atom handling shared by the NICS Placer and NICS Grid windows.

Both windows add, remove, relabel and look for the same ghost atoms, so the
rules live here once: which labels count as a ghost, when a new probe would
sit on top of an existing one, and how to tell that the molecule has changed
under an open window.
"""

import logging

import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D

GHOST_SYMBOLS = ("Bq", "H:")  # all recognised ghost atom labels

#: Two ghosts closer than this (Å) are the same probe. Coincident ghosts add
#: nothing but a second copy of the same shielding value, and an ORCA "H:"
#: ghost carries basis functions, so a stacked pair makes the basis linearly
#: dependent.
DUPLICATE_TOL = 0.01

#: Window keys the plugin registers its two dialogs under.
WINDOW_KEYS = ("main_panel", "grid_panel")


def is_ghost(atom) -> bool:
    return atom.HasProp("custom_symbol") and atom.GetProp("custom_symbol") in (
        GHOST_SYMBOLS
    )


def ghost_positions(mol) -> np.ndarray:
    """Coordinates of every ghost atom in *mol*, shape (N, 3)."""
    if mol is None or not mol.GetNumConformers():
        return np.zeros((0, 3))
    conf = mol.GetConformer()
    pts = [[*conf.GetAtomPosition(a.GetIdx())] for a in mol.GetAtoms() if is_ghost(a)]
    return np.array(pts, dtype=float) if pts else np.zeros((0, 3))


def _cell(p):
    return tuple(int(c) for c in np.floor(np.asarray(p, dtype=float) / DUPLICATE_TOL))


def new_probe_positions(mol, positions) -> tuple:
    """Split *positions* into (to_place, n_skipped).

    A position is skipped when a ghost already sits within DUPLICATE_TOL of it,
    or when an earlier position in the same batch does. Uses a spatial hash so
    a grid of thousands of probes is checked in linear time.
    """
    buckets = {}

    def _near(p):
        cx, cy, cz = _cell(p)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for q in buckets.get((cx + dx, cy + dy, cz + dz), ()):
                        if float(((q - p) ** 2).sum()) < DUPLICATE_TOL**2:
                            return True
        return False

    def _add(p):
        buckets.setdefault(_cell(p), []).append(p)

    for q in ghost_positions(mol):
        _add(q)
    kept, skipped = [], 0
    for p in positions:
        p = np.asarray(p, dtype=float)
        if _near(p):
            skipped += 1
            continue
        _add(p)
        kept.append(p)
    return kept, skipped


def _finish(rw):
    """Refresh derived state after adding or removing isolated atoms.

    Deliberately not SanitizeMol: that re-perceives aromaticity over the
    whole molecule and turns a Kekulé structure's bonds into AROMATIC, which
    rewrites the user's molecule just because a probe was added. An unbonded
    dummy needs only its property cache and the ring table refreshed.
    """
    rw.UpdatePropertyCache(strict=False)
    try:
        Chem.GetSymmSSSR(rw)
    except Exception as _e:  # pragma: no cover - ring perception rarely fails
        logging.warning("[ghosts.py] ring perception: %s", _e)
    return rw.GetMol()


def add_ghost_atoms(mol, positions, symbol: str = "Bq"):
    """Return a new Mol with a ghost dummy atom appended at each of *positions*.

    Batched deliberately: rebuilding the conformer once per atom makes placing
    a grid quadratic in the number of probes, and a 3D volume can run to
    thousands.
    """
    positions = list(positions)
    rw = Chem.RWMol(mol)
    new_idx = []
    for _p in positions:
        atom = Chem.Atom(0)
        atom.SetProp("custom_symbol", symbol)
        new_idx.append(rw.AddAtom(atom))

    old_conf = mol.GetConformer()
    new_conf = Chem.Conformer(rw.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        p = old_conf.GetAtomPosition(i)
        new_conf.SetAtomPosition(i, Point3D(p.x, p.y, p.z))
    for idx, pos in zip(new_idx, positions):
        new_conf.SetAtomPosition(
            idx, Point3D(float(pos[0]), float(pos[1]), float(pos[2]))
        )
    rw.RemoveAllConformers()
    rw.AddConformer(new_conf)
    return _finish(rw)


def remove_ghost_atoms(mol):
    """Return a new Mol with every Bq / H: ghost atom removed."""
    rw = Chem.RWMol(mol)
    for idx in sorted((a.GetIdx() for a in rw.GetAtoms() if is_ghost(a)), reverse=True):
        rw.RemoveAtom(idx)
    return _finish(rw)


def retag_dummy_atoms(mol, symbol: str) -> bool:
    """Relabel every atomic-number-0 atom in *mol* to *symbol*, in place.

    Returns True when anything changed. Keeping every probe on one label
    matters: ORCA reads "H:" and does not know "Bq", so a molecule carrying
    both cannot be run in either program.
    """
    changed = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue
        old = atom.GetProp("custom_symbol") if atom.HasProp("custom_symbol") else None
        if old != symbol:
            atom.SetProp("custom_symbol", symbol)
            changed = True
    return changed


def molecule_signature(mol):
    """Something that changes whenever the geometry an open window uses does.

    Object identity alone misses the host's in-place edits (dragging atoms,
    alignment, geometry optimisation all move coordinates inside the same Mol),
    so the coordinates are part of it.
    """
    if mol is None:
        return None
    if not mol.GetNumConformers():
        return (mol.GetNumAtoms(), None)
    return (mol.GetNumAtoms(), mol.GetConformer().GetPositions().tobytes())


class MoleculeWatcher:
    """Tells a window when the host's current molecule has changed.

    Holds a reference to the last molecule rather than its id(): a freed
    molecule's address can be reused by the next one, which an id() check
    would read as "unchanged".
    """

    def __init__(self, mol=None):
        self.reset(mol)

    def reset(self, mol):
        self._mol = mol
        self._sig = molecule_signature(mol)

    def changed(self, mol) -> bool:
        sig = molecule_signature(mol)
        if mol is self._mol and sig == self._sig:
            return False
        self._mol, self._sig = mol, sig
        return True


def sync_other_windows(context, source) -> None:
    """Push the shared ghost label into every other open plugin window."""
    for key in WINDOW_KEYS:
        try:
            win = context.get_window(key)
        except Exception:
            win = None
        if win is None or win is source or not hasattr(
            win, "sync_symbol_from_settings"
        ):
            continue
        try:
            win.sync_symbol_from_settings()
        except Exception as _e:
            logging.warning("[ghosts.py] sync %s: %s", key, _e)
