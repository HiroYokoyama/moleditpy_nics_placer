"""
NicsGridDialog — place a sheet or a volume of NICS probes for a NICS scan.

Where the main NICS Placer dialog puts three probes per ring, this one lays
down an N x N plane (2D) or an N x N x M box (3D) so the shielding can be
plotted as a surface or contoured as an isochemical-shielding surface (ICSS)
rather than read off as single numbers. Two families of plane are offered —
in 3D they fix the orientation of the box:

  Ring frame  – follows the selected ring. 'parallel' at offset 1.0 is the
                face map (a NICS(1) surface); the two perpendicular planes
                contain the ring normal and give the side-on scan that shows
                the ring-current cone.
  Lab frame   – fixed XY / XZ / YZ cuts through the whole molecule.

By default the grid is centred on the molecular centre of mass; unchecking
that falls back to the ring centroid (ring frame) or the heavy-atom bounding
box (lab frame).

Every probe is an ordinary ghost atom (Bq for Gaussian, H: for ORCA), so the
result drops straight into ORCA Input Generator Pro (with H:) or Gaussian Input
Generator Neo (with Bq) exactly like the single-point probes do. Grids get
large fast: N=21 is 441 atoms, and NMR cost grows with the number of centres,
so the point count is shown live and confirmed above a threshold.
"""

import logging

import numpy as np
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)
import pyvista as pv

from . import PLUGIN_NAME, PLUGIN_VERSION, _plugin_settings, _save_plugin_settings
from .dialog import _GHOST_SYMBOLS, _add_bq_atoms, _remove_all_bq
from .nics_math import (
    GRID_PLANE_LABELS,
    GRID_PLANES,
    LAB_GRID_PLANES,
    axis_extents,
    center_of_mass,
    compute_nics_grid,
    compute_nics_volume,
    counts_for_spacing,
    get_ring_positions,
    get_rings,
    molecular_reference_normal,
    molecule_bounds,
    ring_centroid,
    snap_extents_to_spacing,
)

#: Above this many probes the dialog asks before committing. An NMR job scales
#: badly with ghost centres, and 441 (N=21) is already a long calculation --
#: this is a "did you mean it", not a hard limit.
_CONFIRM_ABOVE = 400

_GRID_SPHERE_RADIUS = 0.12  # smaller than the single-probe spheres: grids are dense

#: Refuse to build a preview larger than this. Not a chemistry limit — an NMR
#: job with this many ghost centres is not a real calculation anyway — but the
#: grid is rebuilt on every spinbox tick, and 101^3 takes ~3 s, which reads as
#: a hang. 200k builds in well under a second.
_MAX_BUILD = 200_000

#: Marks which ring a ring-frame grid is anchored to. Bigger than a probe
#: sphere so it reads as a landmark rather than as one more grid point.
_RING_MARKER_ACTOR = "nics_grid_ring_center"
_RING_MARKER_RADIUS = 0.22

#: Preview spheres are decimated above this count — see _render_spheres.
_PREVIEW_MAX = 2000

#: Control defaults, shared by the initial build and the Reset button so the
#: two cannot drift apart.
_DEFAULT_POINTS = 9
_DEFAULT_EXTENT = 3.0
_DEFAULT_SPACING = 1.0
_DEFAULT_MARGIN = 2.0
_DEFAULT_CONFIRM_ABOVE = _CONFIRM_ABOVE


class NicsGridDialog(QDialog):
    """Configure and place a 2D grid of NICS probes."""

    def __init__(self, context, parent=None):
        super().__init__(parent)
        self._context = context
        self.setWindowTitle(f"NICS Grid (2D, 3D)  —  {PLUGIN_NAME} v{PLUGIN_VERSION}")
        self.resize(600, 820)
        self.setWindowFlags(self.windowFlags() | Qt.WindowType.WindowMinMaxButtonsHint)
        self._rings: list = []
        self._grid_points: list = []
        self._ghost_symbol: str = _plugin_settings.get("ghost_symbol", "Bq")
        self._setup_ui()

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------

    def _setup_ui(self):
        layout = QVBoxLayout(self)

        # Mode and plane come first: they decide what the rest of the dialog
        # even means -- which axes exist, and whether a ring anchors the grid.
        # Asking for them first stops the settings below from silently
        # changing meaning underneath the user.
        top = QFormLayout()

        self._mode_combo = QComboBox()
        self._mode_combo.addItem("2D - plane of probes", "2d")
        self._mode_combo.addItem("3D - volume of probes (ICSS)", "3d")
        self._mode_combo.currentIndexChanged.connect(self._on_mode_changed)
        top.addRow("Mode:", self._mode_combo)

        self._plane_combo = QComboBox()
        for key in GRID_PLANES:
            self._plane_combo.addItem(GRID_PLANE_LABELS[key], key)
        self._plane_combo.currentIndexChanged.connect(self._on_plane_changed)
        top.addRow("Plane:", self._plane_combo)
        layout.addLayout(top)

        self._hint_label = QLabel()
        self._hint_label.setWordWrap(True)
        layout.addWidget(self._hint_label)

        # The ring table lives in its own box so it can be hidden wholesale for
        # the lab planes, which are molecule-wide and have no anchoring ring.
        # A disabled-but-visible table just invites clicking at it.
        self._ring_group = QGroupBox("Anchor ring")
        ring_layout = QVBoxLayout(self._ring_group)
        self._table = QTableWidget()
        self._table.setColumnCount(4)
        self._table.setHorizontalHeaderLabels(["Ring", "Size", "Aromatic", "Planarity"])
        self._table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        # Enough for four or five rings without scrolling — fused systems
        # routinely have that many, and the table is how you pick between
        # them. It still grows if the dialog is resized.
        self._table.setMinimumHeight(170)
        self._table.setMaximumHeight(260)
        self._table.itemSelectionChanged.connect(self._on_params_changed)
        ring_layout.addWidget(self._table)
        layout.addWidget(self._ring_group)

        # ---- Extent -------------------------------------------------------
        extent_group = QGroupBox("Extent")
        extent_form = QFormLayout(extent_group)

        # Centre of mass is the frame most NICS work is quoted in, so it is the
        # default. Note what unchecking does: a ring-frame grid falls back to
        # its own ring's centroid (which is where a single-ring face map wants
        # to be), a lab grid to the heavy-atom bounding-box centre.
        self._use_com = QCheckBox("Molecular centre of mass")
        self._use_com.setChecked(True)
        self._use_com.setToolTip(
            "Unchecked: ring-frame grids centre on their own ring, lab grids "
            "on the heavy-atom bounding box. The plane still sets the "
            "orientation either way."
        )
        self._use_com.toggled.connect(self._on_params_changed)
        extent_form.addRow("Centre:", self._use_com)

        margin_row = QHBoxLayout()
        self._auto_extent = QCheckBox("Fit each axis to the molecule, margin")
        self._auto_extent.setChecked(True)
        self._auto_extent.toggled.connect(self._on_auto_toggled)
        margin_row.addWidget(self._auto_extent)
        self._margin_spin = QDoubleSpinBox()
        self._margin_spin.setRange(0.0, 20.0)
        self._margin_spin.setSingleStep(0.5)
        self._margin_spin.setDecimals(2)
        self._margin_spin.setValue(_DEFAULT_MARGIN)
        self._margin_spin.setSuffix(" A")
        self._margin_spin.setToolTip(
            "Padding beyond the outermost atom on each axis, so the "
            "ring-current field is sampled where it has decayed rather than "
            "clipped at the last atom. ~2 A is about one van der Waals radius."
        )
        self._margin_spin.valueChanged.connect(self._on_params_changed)
        margin_row.addWidget(self._margin_spin)
        margin_row.addStretch()
        extent_form.addRow("Auto:", margin_row)

        self._offset_spin = QDoubleSpinBox()
        self._offset_spin.setRange(-50.0, 50.0)
        self._offset_spin.setSingleStep(0.5)
        self._offset_spin.setDecimals(2)
        self._offset_spin.setValue(0.0)
        self._offset_spin.setSuffix(" A")
        self._offset_spin.valueChanged.connect(self._on_params_changed)
        self._offset_row_label = QLabel("Offset along normal:")
        extent_form.addRow(self._offset_row_label, self._offset_spin)
        layout.addWidget(extent_group)

        # ---- Sampling -----------------------------------------------------
        sampling_group = QGroupBox("Sampling")
        sampling_form = QFormLayout(sampling_group)

        # Uniform spacing derives the counts from a step size instead. With
        # independent half-widths, equal counts do NOT mean equal spacing, and
        # it is the spacing that decides whether the sampled field is smooth
        # enough to contour.
        spacing_row = QHBoxLayout()
        self._uniform_spacing = QCheckBox("Uniform spacing (cubic cells)")
        self._uniform_spacing.toggled.connect(self._on_uniform_toggled)
        spacing_row.addWidget(self._uniform_spacing)
        self._spacing_spin = QDoubleSpinBox()
        self._spacing_spin.setRange(0.05, 10.0)
        self._spacing_spin.setSingleStep(0.1)
        self._spacing_spin.setDecimals(2)
        self._spacing_spin.setValue(_DEFAULT_SPACING)
        self._spacing_spin.setSuffix(" A")
        self._spacing_spin.setEnabled(False)
        self._spacing_spin.valueChanged.connect(self._on_params_changed)
        spacing_row.addWidget(self._spacing_spin)
        spacing_row.addStretch()
        sampling_form.addRow("Step:", spacing_row)

        # One row per axis: count + half-width, so a rectangular grid is as
        # easy to ask for as a square one. Molecules are rarely square.
        self._axis_rows = []
        for slot in range(3):
            row = QHBoxLayout()
            n_spin = QSpinBox()
            n_spin.setRange(1, 101)
            n_spin.setValue(_DEFAULT_POINTS)
            n_spin.valueChanged.connect(self._on_params_changed)
            row.addWidget(n_spin)
            row.addWidget(QLabel("points over +/-"))
            e_spin = QDoubleSpinBox()
            e_spin.setRange(0.0, 100.0)
            e_spin.setSingleStep(0.5)
            e_spin.setDecimals(2)
            e_spin.setValue(_DEFAULT_EXTENT)
            e_spin.setSuffix(" A")
            e_spin.setEnabled(False)
            e_spin.valueChanged.connect(self._on_params_changed)
            row.addWidget(e_spin)
            row.addStretch()
            label = QLabel(f"Axis {slot + 1}:")
            sampling_form.addRow(label, row)
            self._axis_rows.append({"label": label, "n": n_spin, "e": e_spin})
        layout.addWidget(sampling_group)

        # ---- Output -------------------------------------------------------
        output_group = QGroupBox("Output")
        output_form = QFormLayout(output_group)

        self._sym_combo = QComboBox()
        self._sym_combo.addItem("Bq  (Gaussian)", "Bq")
        self._sym_combo.addItem("H:  (ORCA native)", "H:")
        self._sym_combo.setToolTip(
            "Select ghost atom symbol: 'Bq' for Gaussian, 'H:' for ORCA."
        )
        self._sym_combo.currentIndexChanged.connect(self._on_symbol_changed)
        output_form.addRow("Ghost atom label:", self._sym_combo)

        confirm_row = QHBoxLayout()
        self._confirm_spin = QSpinBox()
        self._confirm_spin.setRange(0, 100000)
        self._confirm_spin.setValue(_DEFAULT_CONFIRM_ABOVE)
        self._confirm_spin.setToolTip(
            "Ask before placing more than this many probes. NMR cost grows "
            "with the number of ghost centres. Set to 0 to be asked every time."
        )
        self._confirm_spin.valueChanged.connect(self._on_params_changed)
        confirm_row.addWidget(self._confirm_spin)
        confirm_row.addWidget(QLabel("probes"))
        confirm_row.addStretch()
        output_form.addRow("Confirm above:", confirm_row)
        layout.addWidget(output_group)

        self._count_label = QLabel()
        self._count_label.setWordWrap(True)
        layout.addWidget(self._count_label)

        row = QHBoxLayout()
        self._btn_place = QPushButton("Place Grid")
        self._btn_place.setStyleSheet("font-weight:bold; padding:6px;")
        self._btn_place.clicked.connect(self._place_grid)
        row.addWidget(self._btn_place)

        btn_clear = QPushButton("Clear All Probes")
        btn_clear.clicked.connect(self._clear_all_bq)
        row.addWidget(btn_clear)

        btn_refresh = QPushButton("Refresh")
        btn_refresh.clicked.connect(self._reload)
        row.addWidget(btn_refresh)

        btn_reset = QPushButton("Reset Settings")
        btn_reset.setToolTip("Restore the grid controls to their defaults.")
        btn_reset.clicked.connect(self._reset_settings)
        row.addWidget(btn_reset)

        row.addStretch()
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        row.addWidget(btn_close)
        layout.addLayout(row)

    # ------------------------------------------------------------------
    # Parameter handling
    # ------------------------------------------------------------------

    def _current_plane(self) -> str:
        data = self._plane_combo.currentData()
        return data if data in GRID_PLANES else "parallel"

    def _selected_ring(self) -> int:
        rows = {idx.row() for idx in self._table.selectedIndexes()}
        return min(rows) if rows else 0

    def _is_3d(self) -> bool:
        return self._mode_combo.currentData() == "3d"

    def _axis_labels(self) -> tuple:
        """Names for the three grid axes under the current plane."""
        plane = self._current_plane()
        if plane in LAB_GRID_PLANES:
            first, second = plane[0].upper(), plane[1].upper()
            third = next(c for c in "XYZ" if c not in (first, second))
            return (first, second, third)
        if plane == "parallel":
            return ("u (in plane)", "v (in plane)", "n (normal)")
        if plane == "perpendicular_u":
            return ("u (in plane)", "n (normal)", "v (in plane)")
        return ("v (in plane)", "n (normal)", "u (in plane)")

    def _active_axis_count(self) -> int:
        return 3 if self._is_3d() else 2

    def _refresh_axis_rows(self):
        """Show one row per active axis, named for the current plane."""
        names = self._axis_labels()
        active = self._active_axis_count()
        for slot, row in enumerate(self._axis_rows):
            visible = slot < active
            row["label"].setVisible(visible)
            row["n"].setVisible(visible)
            row["e"].setVisible(visible)
            row["label"].setText(f"Axis {names[slot]}:")

    def _apply_plane_state(self):
        """Widget enabling that follows from the plane. Does not rebuild.

        Split out from the signal handler so `_reset_settings`, which changes
        the plane with signals blocked, can reapply it without triggering a
        second grid build.
        """
        # A ring-frame grid is anchored to one ring; lab planes are
        # molecule-wide, so the table does not apply. Greyed out rather than
        # hidden: hiding it reflows everything below on every plane change,
        # and the rings are still worth reading even when they are not
        # steering the grid.
        ring_frame = self._current_plane() not in LAB_GRID_PLANES
        self._ring_group.setEnabled(ring_frame)
        self._table.setEnabled(ring_frame)
        self._hint_label.setText(
            "Follows the ring selected below; its centre is marked green in "
            "the 3D view."
            if ring_frame
            else "Fixed laboratory axes across the whole molecule; no anchor "
            "ring applies."
        )

    def _on_plane_changed(self, _index):
        self._apply_plane_state()
        self._refresh_axis_rows()
        self._on_params_changed()

    def _apply_mode_state(self):
        """Widget enabling that follows from the mode. Does not rebuild."""
        is3d = self._is_3d()
        # In 2D the offset picks which slice you are sampling -- 1 A is the
        # NICS(1) face map. In 3D the box already spans the normal direction
        # symmetrically, so an offset only slides it off the centre you just
        # chose: redundant at best, misleading at worst. Force it to zero.
        if is3d:
            self._offset_spin.blockSignals(True)
            self._offset_spin.setValue(0.0)
            self._offset_spin.blockSignals(False)
        self._offset_spin.setEnabled(not is3d)
        self._offset_row_label.setEnabled(not is3d)
        self._offset_row_label.setText(
            "Offset along normal (2D only):" if is3d else "Offset along normal:"
        )
        self._plane_combo.setToolTip(
            "In 3D mode the plane only fixes the box orientation." if is3d else ""
        )

    def _on_mode_changed(self, _index):
        self._apply_mode_state()
        self._refresh_axis_rows()
        self._on_params_changed()

    def _on_auto_toggled(self, checked):
        for row in self._axis_rows:
            row["e"].setEnabled(not checked)
        self._on_params_changed()

    def _on_uniform_toggled(self, checked):
        # Counts become a function of the step size, so editing them directly
        # would be a lie -- lock them and let the label report what came out.
        self._spacing_spin.setEnabled(checked)
        for row in self._axis_rows:
            row["n"].setEnabled(not checked)
        self._on_params_changed()

    def _reset_settings(self):
        """Restore every grid control to its default.

        Signals are blocked throughout so the grid is rebuilt once at the end
        rather than once per widget — a half-reset intermediate state can be an
        expensive grid (small step, large extent) that nobody asked for.
        The ghost-atom label is deliberately left alone: it is a persisted
        user/project preference, not a grid parameter.
        """
        widgets = [
            self._mode_combo,
            self._plane_combo,
            self._use_com,
            self._auto_extent,
            self._uniform_spacing,
            self._spacing_spin,
            self._offset_spin,
        ]
        widgets += [row["n"] for row in self._axis_rows]
        widgets += [row["e"] for row in self._axis_rows]
        for w in widgets:
            w.blockSignals(True)
        try:
            self._mode_combo.setCurrentIndex(0)
            self._plane_combo.setCurrentIndex(0)
            self._use_com.setChecked(True)
            self._auto_extent.setChecked(True)
            self._uniform_spacing.setChecked(False)
            self._spacing_spin.setValue(_DEFAULT_SPACING)
            self._spacing_spin.setEnabled(False)
            self._offset_spin.setValue(0.0)
            self._margin_spin.setValue(_DEFAULT_MARGIN)
            self._confirm_spin.setValue(_DEFAULT_CONFIRM_ABOVE)
            for row in self._axis_rows:
                row["n"].setValue(_DEFAULT_POINTS)
                row["n"].setEnabled(True)
                row["e"].setValue(_DEFAULT_EXTENT)
                row["e"].setEnabled(False)
        finally:
            for w in widgets:
                w.blockSignals(False)
        # Reapply every state that normally rides on a signal. The plane was
        # changed with signals blocked, so without this a reset from a lab
        # plane leaves the grid ring-anchored again while the ring table stays
        # greyed out and unpickable.
        self._apply_mode_state()
        self._apply_plane_state()
        self._refresh_axis_rows()
        self._on_params_changed()

    def _on_symbol_changed(self, _index):
        self._ghost_symbol = self._sym_combo.currentData()
        _plugin_settings["ghost_symbol"] = self._ghost_symbol
        _save_plugin_settings(_plugin_settings)

    def _on_params_changed(self, *_args):
        self._rebuild_grid()
        self._render_spheres()

    def _rebuild_grid(self):
        """Recompute probe positions from the current controls."""
        self._grid_points = []
        mol = self._context.current_molecule
        if not mol or not mol.GetNumConformers():
            self._count_label.setText("<i>No 3D molecule loaded.</i>")
            return

        plane = self._current_plane()
        is3d = self._is_3d()
        n_axes = self._active_axis_count()
        try:
            if plane in LAB_GRID_PLANES:
                positions = None
                center = molecule_bounds(mol)[0]
            else:
                if not self._rings:
                    self._count_label.setText("<i>No rings detected.</i>")
                    return
                ring = self._rings[min(self._selected_ring(), len(self._rings) - 1)]
                positions = get_ring_positions(mol, ring["atoms"])
                # None lets the grid builders fall back to the ring centroid.
                center = None
            if self._use_com.isChecked():
                center = center_of_mass(mol)

            reference = molecular_reference_normal(mol)
            if self._auto_extent.isChecked():
                extents = list(
                    axis_extents(
                        mol,
                        plane,
                        positions=positions,
                        reference=reference,
                        center=center,
                        margin=self._margin_spin.value(),
                    )
                )
            else:
                extents = [row["e"].value() for row in self._axis_rows]

            if self._uniform_spacing.isChecked():
                # Grow each half-width to a whole number of steps FIRST.
                # Rounding only the count leaves each axis with a different
                # actual step -- half-widths of 4.4/4.3/4.2 at a 1.0 A request
                # all give 10 points but step 0.978/0.956/0.933.
                extents = list(
                    snap_extents_to_spacing(extents, self._spacing_spin.value())
                )

            for row, value in zip(self._axis_rows, extents):
                row["e"].blockSignals(True)
                row["e"].setValue(value)
                row["e"].blockSignals(False)
            extents = [row["e"].value() for row in self._axis_rows]

            if self._uniform_spacing.isChecked():
                counts = list(counts_for_spacing(extents, self._spacing_spin.value()))
                for row, value in zip(self._axis_rows, counts):
                    row["n"].blockSignals(True)
                    row["n"].setValue(value)
                    row["n"].blockSignals(False)
                    # setValue clamps to the spinbox range; read back what the
                    # widget actually holds so the label cannot claim a count
                    # the grid does not have.
                counts = [row["n"].value() for row in self._axis_rows]
            else:
                counts = [row["n"].value() for row in self._axis_rows]

            # Check the size before building it. The grid is rebuilt on every
            # spinbox tick, and at the maximum 101 points on all three axes
            # that is over a million probes and ~3 s per keystroke — the
            # dialog would appear to hang while you were still typing.
            requested = counts[0] * counts[1] * (counts[2] if is3d else 1)
            if requested > _MAX_BUILD:
                self._grid_points = []
                self._count_label.setText(
                    f"<b style='color:red'>{requested:,} probes is too many to "
                    f"preview</b> (limit {_MAX_BUILD:,}). Increase the step or "
                    "reduce the point counts."
                )
                return

            common = dict(
                plane=plane,
                n_points=counts[:2],
                extent=extents[:2],
                offset=self._offset_spin.value(),
                reference=reference,
                center=center,
            )
            if is3d:
                self._grid_points = compute_nics_volume(
                    positions,
                    n_normal=counts[2],
                    normal_extent=extents[2],
                    **common,
                )
            else:
                self._grid_points = compute_nics_grid(positions, **common)
        except Exception as _e:
            logging.warning("[grid_dialog.py] _rebuild_grid: %s", _e)
            self._count_label.setText(f"<b style='color:red'>Grid error:</b> {_e}")
            return

        total = len(self._grid_points)
        warn = (
            "  <b style='color:#b26b00'>&mdash; large; NMR cost grows with "
            "ghost centres.</b>"
            if total > self._confirm_spin.value()
            else ""
        )
        names = self._axis_labels()
        shape = " x ".join(str(c) for c in counts[:n_axes])
        detail = "; ".join(
            f"{names[k]} {self._step(counts[k], extents[k]):.2f} A "
            f"over +/-{extents[k]:.2f} A"
            for k in range(n_axes)
        )
        self._count_label.setText(
            f"{shape} = <b>{total}</b> probes &mdash; {detail}.{warn}"
        )

    @staticmethod
    def _step(count, extent):
        return (2.0 * extent) / (count - 1) if count > 1 else 0.0

    # ------------------------------------------------------------------
    # Rings
    # ------------------------------------------------------------------

    def _reload(self):
        mol = self._context.current_molecule
        self._rings = get_rings(mol) if mol and mol.GetNumConformers() else []
        self._table.setRowCount(len(self._rings))
        for i, ring in enumerate(self._rings):
            planarity = ""
            try:
                pts = get_ring_positions(mol, ring["atoms"])
                from .nics_math import planarity_rms

                planarity = f"{planarity_rms(pts):.3f} A"
            except Exception as _e:
                logging.warning("[grid_dialog.py] ring %d planarity: %s", i, _e)
            for col, val in enumerate(
                [
                    str(i + 1),
                    str(len(ring["atoms"])),
                    "Yes" if ring["is_aromatic"] else "No",
                    planarity,
                ]
            ):
                self._table.setItem(i, col, QTableWidgetItem(val))
        if self._rings and not self._table.selectedIndexes():
            self._table.selectRow(0)
        self._on_params_changed()

    # ------------------------------------------------------------------
    # Rendering
    # ------------------------------------------------------------------

    def _render_spheres(self):
        plotter = self._context.plotter
        if plotter is None:
            return
        try:
            plotter.remove_actor("nics_grid")
            if self._grid_points:
                pts = np.array([p["pos"] for p in self._grid_points], dtype=float)
                # A 3D box runs to thousands of spheres; glyphing all of them
                # stalls the viewer on every spinbox tick. The preview only has
                # to show where the box sits, so thin it out — placement still
                # uses every point.
                preview_max = _PREVIEW_MAX
                if len(pts) > preview_max:
                    stride = int(np.ceil(len(pts) / preview_max))
                    pts = pts[::stride]
                poly = pv.PolyData(pts)
                poly["r"] = [_GRID_SPHERE_RADIUS] * len(pts)
                mesh = poly.glyph(geom=pv.Sphere(radius=1.0), scale="r", orient=False)
                plotter.add_mesh(
                    mesh,
                    name="nics_grid",
                    color="deepskyblue",
                    opacity=0.5,
                    pickable=False,
                )
            self._render_ring_marker(plotter)
            plotter.render()
        except Exception as _e:
            logging.warning("[grid_dialog.py] _render_spheres: %s", _e)

    def _render_ring_marker(self, plotter):
        """Green sphere on the centroid of the ring the grid is anchored to.

        Only for ring-frame planes: a lab grid is molecule-wide and the ring
        table is disabled there, so a highlight would point at a ring that has
        no bearing on the grid.
        """
        plotter.remove_actor(_RING_MARKER_ACTOR)
        if self._current_plane() in LAB_GRID_PLANES or not self._rings:
            return
        mol = self._context.current_molecule
        if not mol or not mol.GetNumConformers():
            return
        ring = self._rings[min(self._selected_ring(), len(self._rings) - 1)]
        centroid = ring_centroid(get_ring_positions(mol, ring["atoms"]))
        plotter.add_mesh(
            pv.Sphere(
                radius=_RING_MARKER_RADIUS, center=tuple(float(c) for c in centroid)
            ),
            name=_RING_MARKER_ACTOR,
            color="limegreen",
            opacity=0.75,
            pickable=False,
        )

    def _clear_actors(self):
        try:
            plotter = self._context.plotter
            if plotter:
                plotter.remove_actor("nics_grid")
                plotter.remove_actor(_RING_MARKER_ACTOR)
                plotter.render()
        except Exception as _e:
            logging.warning("[grid_dialog.py] _clear_actors: %s", _e)

    # ------------------------------------------------------------------
    # Placement
    # ------------------------------------------------------------------

    def _place_grid(self):
        if not self._grid_points:
            return
        mol = self._context.current_molecule
        if not mol:
            return
        total = len(self._grid_points)
        if total > self._confirm_spin.value():
            reply = QMessageBox.question(
                self,
                "Large grid",
                f"This will add {total} ghost atoms to the molecule.\n\n"
                "NMR calculation cost grows with the number of ghost centres, "
                "and the 3D view will slow down.\n\nPlace them anyway?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if reply != QMessageBox.StandardButton.Yes:
                return
        positions = [p["pos"] for p in self._grid_points]
        self._context.current_molecule = _add_bq_atoms(
            mol, positions, symbol=self._ghost_symbol
        )
        self._context.push_undo_checkpoint()
        self._context.show_status_message(f"Placed {total} NICS grid probes.", 3000)
        QTimer.singleShot(150, self._render_spheres)

    def _clear_all_bq(self):
        mol = self._context.current_molecule
        if not mol:
            return
        self._context.current_molecule = _remove_all_bq(mol)
        self._context.push_undo_checkpoint()
        QTimer.singleShot(150, self._render_spheres)

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def sync_symbol_from_settings(self):
        sym = _plugin_settings.get("ghost_symbol", "Bq")
        idx = self._sym_combo.findData(sym)
        if idx >= 0:
            self._sym_combo.blockSignals(True)
            self._sym_combo.setCurrentIndex(idx)
            self._sym_combo.blockSignals(False)
        self._ghost_symbol = sym if sym in _GHOST_SYMBOLS else "Bq"

    def showEvent(self, event):
        super().showEvent(event)
        self.sync_symbol_from_settings()
        self._on_mode_changed(self._mode_combo.currentIndex())
        self._on_plane_changed(self._plane_combo.currentIndex())
        self._reload()

    def closeEvent(self, event):
        self._clear_actors()
        super().closeEvent(event)
