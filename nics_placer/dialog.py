"""
NicsPlacerDialog — ring detection, NICS probe visualisation, and Bq placement.

Sphere colours:
  Yellow (semi-transparent) = unselected probe position
  Red    (semi-transparent) = staged for placement (click to toggle)
  Green                     = already placed (ghost atom exists)

Workflow:
  1. Click a yellow sphere  → turns red  (staged)
  2. Click a red sphere     → turns yellow (unstaged)
  3. Press "Apply"          → places ghost atom for all red spheres
     OR press "Place All"  → stages + places everything in one shot

Ghost atom symbol options:
  Bq  — Gaussian convention; also valid in ORCA (recognised by ORCA Input Generator Pro)
  H:  — ORCA native ghost atom notation

The 'custom_symbol' atom property is shared with XYZ Editor and ORCA Input Generator Pro.
"""
import logging

import numpy as np
from PyQt6.QtCore import Qt, QTimer, QObject, QEvent
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)
from rdkit import Chem
from rdkit.Geometry import Point3D
import pyvista as pv

from . import PLUGIN_NAME, PLUGIN_VERSION
from .nics_math import (
    NICS1_HEIGHT,
    compute_nics_points,
    get_ring_positions,
    get_rings,
)

_GHOST_SYMBOLS = ("Bq", "H:")   # all recognised ghost atom labels
_PICK_DIST_SQ = 1.0             # Å² snap threshold
_SPHERE_RADIUS = 0.25


# ---------------------------------------------------------------------------
# Qt event filter for non-drag click detection
# ---------------------------------------------------------------------------

class _ClickFilter(QObject):
    def __init__(self, callback, parent=None):
        super().__init__(parent)
        self._callback = callback
        self._press_pos = None

    def eventFilter(self, obj, event):
        t = event.type()
        if t == QEvent.Type.MouseButtonPress:
            if event.button() == Qt.MouseButton.LeftButton:
                self._press_pos = event.position().toPoint()
        elif t == QEvent.Type.MouseButtonRelease:
            if event.button() == Qt.MouseButton.LeftButton and self._press_pos is not None:
                rel = event.position().toPoint()
                dx = rel.x() - self._press_pos.x()
                dy = rel.y() - self._press_pos.y()
                if dx * dx + dy * dy <= 25:
                    self._callback(rel.x(), rel.y(), obj)
                self._press_pos = None
        return False


# ---------------------------------------------------------------------------
# RDKit molecule helpers
# ---------------------------------------------------------------------------

def _add_bq_atom(mol, position: np.ndarray, symbol: str = "Bq"):
    """Return a new Mol with a ghost dummy atom appended at *position*."""
    rw = Chem.RWMol(mol)
    atom = Chem.Atom(0)
    atom.SetProp("custom_symbol", symbol)
    new_idx = rw.AddAtom(atom)

    old_conf = mol.GetConformer()
    new_conf = Chem.Conformer(rw.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        p = old_conf.GetAtomPosition(i)
        new_conf.SetAtomPosition(i, Point3D(p.x, p.y, p.z))
    new_conf.SetAtomPosition(
        new_idx,
        Point3D(float(position[0]), float(position[1]), float(position[2])),
    )
    rw.RemoveAllConformers()
    rw.AddConformer(new_conf)
    try:
        Chem.SanitizeMol(rw)
    except Exception:
        rw.UpdatePropertyCache(strict=False)
    return rw.GetMol()


def _remove_all_bq(mol):
    """Return a new Mol with every ghost atom (Bq or H:) removed."""
    rw = Chem.RWMol(mol)
    to_remove = [
        a.GetIdx()
        for a in rw.GetAtoms()
        if a.HasProp("custom_symbol") and a.GetProp("custom_symbol") in _GHOST_SYMBOLS
    ]
    for idx in sorted(to_remove, reverse=True):
        rw.RemoveAtom(idx)
    try:
        Chem.SanitizeMol(rw)
    except Exception:
        rw.UpdatePropertyCache(strict=False)
    return rw.GetMol()


# ---------------------------------------------------------------------------
# Main dialog
# ---------------------------------------------------------------------------

_STATE_UNSET   = "unset"    # yellow — not staged
_STATE_STAGED  = "staged"   # red    — will be placed on Apply
_STATE_PLACED  = "placed"   # green  — Bq atom exists in molecule


class NicsPlacerDialog(QDialog):
    """
    Dialog for visualising and placing Bq atoms at NICS probe positions.

    Each NICS point has one of three states:
      unset   → yellow sphere (click to stage)
      staged  → red sphere    (click to un-stage; Apply to commit)
      placed  → green sphere  (already in molecule; Clear All to remove)
    """

    def __init__(self, context, parent=None):
        super().__init__(parent)
        self._context = context
        self.setWindowTitle(f"NICS Placer  —  {PLUGIN_NAME} v{PLUGIN_VERSION}")
        self.resize(620, 450)
        self.setWindowFlags(
            self.windowFlags() | Qt.WindowType.WindowMinMaxButtonsHint
        )

        # List of dicts: {'ring':int, 'type':str, 'pos':np.ndarray, 'state':str}
        self._nics_points: list = []
        self._ghost_symbol: str = "Bq"   # "Bq" or "H:"
        self._click_filter = None
        # Actor references for pick-detection
        self._actor_yellow = None
        self._actor_red = None
        self._actor_green = None

        self._setup_ui()
        self._load_rings()
        self._enable_picking()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _setup_ui(self):
        layout = QVBoxLayout(self)

        info = QLabel(
            "<b>Yellow</b> = available &nbsp;|&nbsp; "
            "<b style='color:red'>Red</b> = staged (click Apply to place) &nbsp;|&nbsp; "
            "<b style='color:green'>Green</b> = placed.  "
            "Click a sphere to toggle staging."
        )
        info.setWordWrap(True)
        layout.addWidget(info)

        sym_row = QHBoxLayout()
        sym_row.addWidget(QLabel("Ghost atom label:"))
        self._sym_combo = QComboBox()
        self._sym_combo.addItem("Bq  (Gaussian / ORCA)", "Bq")
        self._sym_combo.addItem("H:  (ORCA native)",     "H:")
        self._sym_combo.currentIndexChanged.connect(self._on_symbol_changed)
        sym_row.addWidget(self._sym_combo)
        sym_row.addStretch()
        layout.addLayout(sym_row)

        self._table = QTableWidget()
        self._table.setColumnCount(4)
        self._table.setHorizontalHeaderLabels(["Ring", "Size", "Aromatic", "Status"])
        self._table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        layout.addWidget(self._table)

        # Row 1: select-and-stage helpers
        row1 = QHBoxLayout()
        btn_stage_nics0 = QPushButton("Stage NICS(0) for selected")
        btn_stage_nics0.clicked.connect(self._stage_selected_nics0)
        row1.addWidget(btn_stage_nics0)

        btn_stage_nics1 = QPushButton("Stage NICS(1)± for selected")
        btn_stage_nics1.clicked.connect(self._stage_selected_nics1)
        row1.addWidget(btn_stage_nics1)
        layout.addLayout(row1)

        # Row 2: Apply / bulk / clear
        row2 = QHBoxLayout()

        self._btn_apply = QPushButton("Apply  (place red Bq)")
        self._btn_apply.setStyleSheet("font-weight:bold; padding:6px;")
        self._btn_apply.clicked.connect(self._apply_staged)
        row2.addWidget(self._btn_apply)

        btn_place_all = QPushButton("Place All")
        btn_place_all.clicked.connect(self._place_all)
        row2.addWidget(btn_place_all)

        btn_clear = QPushButton("Clear All Bq")
        btn_clear.clicked.connect(self._clear_all_bq)
        row2.addWidget(btn_clear)

        btn_refresh = QPushButton("Refresh Rings")
        btn_refresh.clicked.connect(self._load_rings)
        row2.addWidget(btn_refresh)

        row2.addStretch()
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        row2.addWidget(btn_close)
        layout.addLayout(row2)

    # ------------------------------------------------------------------
    # Ring detection & table
    # ------------------------------------------------------------------

    def _load_rings(self):
        mol = self._context.current_molecule
        self._nics_points = []
        self._table.setRowCount(0)

        if not mol or not mol.GetNumConformers():
            self._render_spheres()
            return

        rings = get_rings(mol)
        self._table.setRowCount(len(rings))

        for i, ring in enumerate(rings):
            try:
                pts = get_ring_positions(mol, ring["atoms"])
                nics = compute_nics_points(pts)
            except Exception as _e:
                logging.warning("[dialog.py] ring %d NICS calc: %s", i, _e)
                for col, val in enumerate([
                    str(i + 1), str(len(ring["atoms"])),
                    "Yes" if ring["is_aromatic"] else "No", "error",
                ]):
                    self._table.setItem(i, col, QTableWidgetItem(val))
                continue

            for ntype, pos in (
                ("nics0",       nics["nics0"]),
                ("nics1_above", nics["nics1_above"]),
                ("nics1_below", nics["nics1_below"]),
            ):
                self._nics_points.append(
                    {"ring": i, "type": ntype, "pos": pos, "state": _STATE_UNSET}
                )

            for col, val in enumerate([
                str(i + 1),
                str(len(ring["atoms"])),
                "Yes" if ring["is_aromatic"] else "No",
                "0/3",
            ]):
                self._table.setItem(i, col, QTableWidgetItem(val))

        self._sync_placed_status()
        self._render_spheres()

    def _sync_placed_status(self):
        """Mark points whose Bq atom already exists in the molecule as placed."""
        mol = self._context.current_molecule
        if not mol or not mol.GetNumConformers():
            return
        conf = mol.GetConformer()
        bq_pos = []
        for atom in mol.GetAtoms():
            if atom.HasProp("custom_symbol") and atom.GetProp("custom_symbol") in _GHOST_SYMBOLS:
                p = conf.GetAtomPosition(atom.GetIdx())
                bq_pos.append(np.array([p.x, p.y, p.z]))

        for pt in self._nics_points:
            for bp in bq_pos:
                if np.sum((pt["pos"] - bp) ** 2) < 0.01:
                    pt["state"] = _STATE_PLACED
                    break

        self._update_table_status()

    def _update_table_status(self):
        for ring_idx in range(self._table.rowCount()):
            pts = [p for p in self._nics_points if p["ring"] == ring_idx]
            n_placed  = sum(1 for p in pts if p["state"] == _STATE_PLACED)
            n_staged  = sum(1 for p in pts if p["state"] == _STATE_STAGED)
            n_total   = len(pts)
            if n_total == 0:
                label = "—"
            elif n_placed == n_total:
                label = f"✓ all placed ({n_total})"
            elif n_staged:
                label = f"{n_placed} placed, {n_staged} staged"
            else:
                label = f"{n_placed}/{n_total} placed"
            self._table.setItem(ring_idx, 3, QTableWidgetItem(label))

    # ------------------------------------------------------------------
    # 3D sphere rendering
    # ------------------------------------------------------------------

    def _render_spheres(self):
        plotter = self._context.plotter
        if plotter is None:
            return
        try:
            plotter.remove_actor("nics_yellow")
            plotter.remove_actor("nics_red")
            plotter.remove_actor("nics_green")
            self._actor_yellow = None
            self._actor_red    = None
            self._actor_green  = None

            def _glyph(positions, radius=_SPHERE_RADIUS):
                poly = pv.PolyData(np.array(positions, dtype=float))
                poly["r"] = [radius] * len(positions)
                return poly.glyph(geom=pv.Sphere(radius=1.0), scale="r", orient=False)

            yellow_pts = [p["pos"] for p in self._nics_points if p["state"] == _STATE_UNSET]
            red_pts    = [p["pos"] for p in self._nics_points if p["state"] == _STATE_STAGED]
            green_pts  = [p["pos"] for p in self._nics_points if p["state"] == _STATE_PLACED]

            if yellow_pts:
                self._actor_yellow = plotter.add_mesh(
                    _glyph(yellow_pts), name="nics_yellow",
                    color="yellow", opacity=0.55, pickable=True,
                )
            if red_pts:
                self._actor_red = plotter.add_mesh(
                    _glyph(red_pts), name="nics_red",
                    color="red", opacity=0.60, pickable=True,
                )
            if green_pts:
                self._actor_green = plotter.add_mesh(
                    _glyph(green_pts), name="nics_green",
                    color="green", opacity=0.55, pickable=False,
                )

            plotter.render()
        except Exception as _e:
            logging.warning("[dialog.py] _render_spheres: %s", _e)

    def _clear_actors(self):
        try:
            plotter = self._context.plotter
            if plotter:
                for name in ("nics_yellow", "nics_red", "nics_green"):
                    plotter.remove_actor(name)
                plotter.render()
        except Exception as _e:
            logging.warning("[dialog.py] _clear_actors: %s", _e)
        self._actor_yellow = self._actor_red = self._actor_green = None

    # ------------------------------------------------------------------
    # 3D click picking → toggle stage
    # ------------------------------------------------------------------

    def _enable_picking(self):
        try:
            plotter = self._context.plotter
            if plotter is None:
                return
            self._click_filter = _ClickFilter(self._on_plotter_click, parent=self)
            plotter.installEventFilter(self._click_filter)
        except Exception as _e:
            logging.warning("[dialog.py] _enable_picking: %s", _e)

    def _disable_picking(self):
        try:
            plotter = self._context.plotter
            if plotter and self._click_filter:
                plotter.removeEventFilter(self._click_filter)
        except Exception as _e:
            logging.warning("[dialog.py] _disable_picking: %s", _e)
        self._click_filter = None

    def _on_plotter_click(self, x, y, widget):
        try:
            import vtk
            plotter = self._context.plotter
            if plotter is None:
                return

            vtk_y = widget.height() - y
            picker = vtk.vtkCellPicker()
            picker.SetTolerance(0.005)
            picker.Pick(x, vtk_y, 0, plotter.renderer)
            picked = picker.GetActor()

            # Only react to yellow or red sphere actors
            if picked not in (self._actor_yellow, self._actor_red):
                return

            pick_pos = np.array(picker.GetPickPosition(), dtype=float)
            # Find nearest toggleable point
            best_idx = -1
            best_sq = float("inf")
            for i, pt in enumerate(self._nics_points):
                if pt["state"] == _STATE_PLACED:
                    continue
                diff = pt["pos"] - pick_pos
                sq = float(np.dot(diff, diff))
                if sq < best_sq:
                    best_sq = sq
                    best_idx = i

            if best_idx >= 0 and best_sq < _PICK_DIST_SQ:
                pt = self._nics_points[best_idx]
                pt["state"] = (
                    _STATE_UNSET if pt["state"] == _STATE_STAGED else _STATE_STAGED
                )
                self._update_table_status()
                self._render_spheres()

        except Exception as _e:
            logging.warning("[dialog.py] _on_plotter_click: %s", _e)

    # ------------------------------------------------------------------
    # Stage helpers
    # ------------------------------------------------------------------

    def _on_symbol_changed(self, _index):
        self._ghost_symbol = self._sym_combo.currentData()

    def _selected_ring_indices(self) -> set:
        selected = {idx.row() for idx in self._table.selectedIndexes()}
        return selected if selected else set(range(self._table.rowCount()))

    def _stage_rings(self, ring_indices: set, types: set):
        changed = False
        for pt in self._nics_points:
            if pt["ring"] in ring_indices and pt["type"] in types and pt["state"] == _STATE_UNSET:
                pt["state"] = _STATE_STAGED
                changed = True
        if changed:
            self._update_table_status()
            self._render_spheres()

    def _stage_selected_nics0(self):
        self._stage_rings(self._selected_ring_indices(), {"nics0"})

    def _stage_selected_nics1(self):
        self._stage_rings(self._selected_ring_indices(), {"nics1_above", "nics1_below"})

    # ------------------------------------------------------------------
    # Placement
    # ------------------------------------------------------------------

    def _apply_staged(self):
        """Place Bq atoms for all red (staged) points."""
        staged = [pt for pt in self._nics_points if pt["state"] == _STATE_STAGED]
        if not staged:
            return
        mol = self._context.current_molecule
        if not mol:
            return
        for pt in staged:
            mol = _add_bq_atom(mol, pt["pos"], symbol=self._ghost_symbol)
            pt["state"] = _STATE_PLACED
        self._context.current_molecule = mol
        self._context.push_undo_checkpoint()
        self._update_table_status()
        QTimer.singleShot(150, self._render_spheres)

    def _place_all(self):
        """Stage every unset point and immediately apply."""
        for pt in self._nics_points:
            if pt["state"] == _STATE_UNSET:
                pt["state"] = _STATE_STAGED
        self._apply_staged()

    def _clear_all_bq(self):
        mol = self._context.current_molecule
        if not mol:
            return
        new_mol = _remove_all_bq(mol)
        self._context.current_molecule = new_mol
        self._context.push_undo_checkpoint()
        for pt in self._nics_points:
            if pt["state"] == _STATE_PLACED:
                pt["state"] = _STATE_UNSET
        self._update_table_status()
        QTimer.singleShot(150, self._render_spheres)

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def closeEvent(self, event):
        self._disable_picking()
        self._clear_actors()
        super().closeEvent(event)
