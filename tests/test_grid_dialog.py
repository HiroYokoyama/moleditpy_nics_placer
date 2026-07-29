"""
Tests for NicsGridDialog — driven headlessly against the rich PyQt6 / pyvista
stubs installed by tests/conftest.py.
"""

import logging
import os
import sys
from unittest.mock import MagicMock

import numpy as np
import pytest

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), "..")))

_AVAILABLE = False
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    import nics_placer.grid_dialog as grid_mod
    from nics_placer.grid_dialog import NicsGridDialog, _CONFIRM_ABOVE, _PREVIEW_MAX

    _AVAILABLE = Chem is not None
except Exception:
    pass

needs_dialog = pytest.mark.skipif(
    not _AVAILABLE, reason="RDKit or nics_placer.grid_dialog not importable"
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _mol_3d(smiles="c1ccccc1", seed=42):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    return mol


class _StubContext:
    def __init__(self, mol=None):
        self.current_molecule = mol
        self.plotter = MagicMock(name="plotter")
        self.plotter.add_mesh.side_effect = lambda mesh, name=None, **kw: MagicMock(
            name=f"actor_{name}"
        )
        self._windows = {}
        self.push_undo_checkpoint = MagicMock()
        self.show_status_message = MagicMock()

    def get_window(self, key):
        return self._windows.get(key)

    def register_window(self, key, win):
        self._windows[key] = win

    def get_main_window(self):
        return None


def _dialog(mol=None, smiles="c1ccccc1"):
    ctx = _StubContext(mol if mol is not None else _mol_3d(smiles))
    dlg = NicsGridDialog(ctx)
    dlg.showEvent(MagicMock())
    return dlg


def _select_plane(dlg, plane):
    from nics_placer.nics_math import GRID_PLANES

    dlg._plane_combo.setCurrentIndex(GRID_PLANES.index(plane))
    dlg._on_plane_changed(dlg._plane_combo.currentIndex())


def _set_counts(dlg, *counts):
    for row, n in zip(dlg._axis_rows, counts):
        row["n"].setValue(n)
    dlg._on_params_changed()


def _set_extents(dlg, *extents):
    dlg._auto_extent.setChecked(False)
    dlg._on_auto_toggled(False)
    for row, e in zip(dlg._axis_rows, extents):
        row["e"].setValue(e)
    dlg._on_params_changed()


def _counts(dlg):
    return [row["n"].value() for row in dlg._axis_rows]


def _extents(dlg):
    return [row["e"].value() for row in dlg._axis_rows]


def _select_mode(dlg, mode):
    dlg._mode_combo.setCurrentIndex(0 if mode == "2d" else 1)
    dlg._on_mode_changed(dlg._mode_combo.currentIndex())


@pytest.fixture(autouse=True)
def _reset_msgbox():
    if _AVAILABLE:
        grid_mod.QMessageBox._calls = []
        grid_mod.QMessageBox._answer = grid_mod.QMessageBox.StandardButton.No
    yield


# ---------------------------------------------------------------------------
# Construction and ring table
# ---------------------------------------------------------------------------


@needs_dialog
def test_dialog_lists_the_rings_with_planarity():
    dlg = _dialog(smiles="c1ccc2ccccc2c1")
    assert dlg._table.rowCount() == 2
    assert dlg._table.item(0, 2).text() == "Yes"
    assert dlg._table.item(0, 3).text().endswith("A")


@needs_dialog
def test_puckered_ring_reports_a_nonzero_planarity():
    dlg = _dialog(smiles="C1CCCCC1")
    value = float(dlg._table.item(0, 3).text().split()[0])
    assert value > 0.1, "a cyclohexane chair is not planar"


@needs_dialog
def test_first_ring_is_selected_by_default():
    dlg = _dialog()
    assert dlg._selected_ring() == 0


@needs_dialog
def test_no_molecule_reports_instead_of_raising():
    dlg = NicsGridDialog(_StubContext(None))
    dlg.showEvent(MagicMock())
    assert dlg._grid_points == []
    assert "No 3D molecule" in dlg._count_label.setText.call_args[0][0]


@needs_dialog
def test_ring_frame_plane_without_rings_reports_no_rings():
    dlg = _dialog(smiles="CCCC")
    _select_plane(dlg, "parallel")
    assert dlg._grid_points == []
    assert "No rings" in dlg._count_label.setText.call_args[0][0]


@needs_dialog
def test_lab_plane_works_without_any_ring():
    dlg = _dialog(smiles="CCCC")
    _select_plane(dlg, "xy")
    n1, n2, _n3 = _counts(dlg)
    assert len(dlg._grid_points) == n1 * n2


# ---------------------------------------------------------------------------
# Grid construction
# ---------------------------------------------------------------------------


@needs_dialog
def test_default_grid_is_2d():
    dlg = _dialog()
    n1, n2, _n3 = _counts(dlg)
    assert len(dlg._grid_points) == n1 * n2


@needs_dialog
def test_changing_the_point_count_rebuilds_the_grid():
    dlg = _dialog()
    _set_counts(dlg, 5, 5)
    assert len(dlg._grid_points) == 25


@needs_dialog
def test_switching_to_3d_produces_a_volume():
    dlg = _dialog()
    _select_mode(dlg, "3d")
    _set_counts(dlg, 4, 4, 3)
    assert len(dlg._grid_points) == 4 * 4 * 3


@needs_dialog
def test_the_third_axis_row_is_hidden_in_2d_mode():
    dlg = _dialog()
    third = dlg._axis_rows[2]
    _select_mode(dlg, "2d")
    third["n"].setVisible.assert_called_with(False)
    third["label"].setVisible.assert_called_with(False)
    _select_mode(dlg, "3d")
    third["n"].setVisible.assert_called_with(True)


@needs_dialog
def test_axis_rows_are_named_for_the_current_plane():
    dlg = _dialog()
    _select_plane(dlg, "xz")
    assert dlg._axis_labels() == ("X", "Z", "Y")
    _select_plane(dlg, "parallel")
    assert dlg._axis_labels() == ("u (in plane)", "v (in plane)", "n (normal)")
    _select_plane(dlg, "perpendicular_u")
    assert dlg._axis_labels() == ("u (in plane)", "n (normal)", "v (in plane)")
    _select_plane(dlg, "perpendicular_v")
    assert dlg._axis_labels() == ("v (in plane)", "n (normal)", "u (in plane)")


@needs_dialog
def test_ring_table_is_disabled_for_lab_planes():
    """A lab-frame grid is molecule-wide, so picking a ring would be
    misleading."""
    dlg = _dialog()
    _select_plane(dlg, "xy")
    dlg._table.setEnabled.assert_called_with(False)
    _select_plane(dlg, "parallel")
    dlg._table.setEnabled.assert_called_with(True)


@needs_dialog
def test_lab_grid_defaults_to_the_centre_of_mass():
    from nics_placer.nics_math import center_of_mass

    dlg = _dialog(smiles="c1ccc2ccccc2c1")
    _select_plane(dlg, "xy")
    _set_counts(dlg, 3, 3)
    centre = next(p["pos"] for p in dlg._grid_points if p["i"] == 1 and p["j"] == 1)
    expected = center_of_mass(dlg._context.current_molecule)
    expected = expected + dlg._offset_spin.value() * np.array([0.0, 0.0, 1.0])
    np.testing.assert_allclose(centre, expected, atol=1e-9)


@needs_dialog
def test_unchecking_com_falls_back_to_the_bounding_box_for_a_lab_grid():
    from nics_placer.nics_math import molecule_bounds

    dlg = _dialog(smiles="c1ccc2ccccc2c1")
    _select_plane(dlg, "xy")
    dlg._use_com.setChecked(False)
    _set_counts(dlg, 3, 3)
    centre = next(p["pos"] for p in dlg._grid_points if p["i"] == 1 and p["j"] == 1)
    expected = molecule_bounds(dlg._context.current_molecule)[0]
    np.testing.assert_allclose(centre, expected, atol=1e-9)


@needs_dialog
def test_unchecking_com_puts_a_ring_grid_back_on_its_ring():
    """A single-ring face map wants to sit over the ring, not over the whole
    molecule -- that is what unchecking is for."""
    from nics_placer.nics_math import get_ring_positions, ring_centroid

    dlg = _dialog(smiles="c1ccc2ccccc2c1")
    _select_plane(dlg, "parallel")
    dlg._use_com.setChecked(False)
    _set_counts(dlg, 3, 3)
    centre = next(p["pos"] for p in dlg._grid_points if p["i"] == 1 and p["j"] == 1)
    ring = dlg._rings[0]
    expected = ring_centroid(
        get_ring_positions(dlg._context.current_molecule, ring["atoms"])
    )
    np.testing.assert_allclose(centre, expected, atol=1e-9)


@needs_dialog
def test_com_and_ring_centroid_differ_for_an_off_centre_ring():
    """Guard against the two options being silently identical -- on
    naphthalene each ring centroid is well away from the molecular centre."""
    from nics_placer.nics_math import center_of_mass, get_ring_positions, ring_centroid

    mol = _mol_3d("c1ccc2ccccc2c1")
    dlg = _dialog(mol=mol)
    ring_c = ring_centroid(get_ring_positions(mol, dlg._rings[0]["atoms"]))
    assert float(np.linalg.norm(center_of_mass(mol) - ring_c)) > 1.0


@needs_dialog
def test_auto_extent_grows_with_the_molecule():
    small = _dialog(smiles="c1ccccc1")
    big = _dialog(smiles="c1ccc2cc3cc4ccccc4cc3cc2c1")
    assert max(_extents(big)) > max(_extents(small))


@needs_dialog
def test_auto_extent_covers_every_heavy_atom_in_a_lab_plane():
    from nics_placer.nics_math import heavy_atom_positions

    from nics_placer.nics_math import molecule_bounds

    dlg = _dialog(smiles="c1ccc2ccccc2c1")
    _select_plane(dlg, "xy")
    e1, e2, e3 = _extents(dlg)
    pts = heavy_atom_positions(dlg._context.current_molecule)
    rel = pts - molecule_bounds(dlg._context.current_molecule)[0]
    # Each axis is sized on its own, so a long flat molecule gets a rectangle
    # rather than a square padded out to its longest dimension.
    assert np.abs(rel[:, 0]).max() <= e1
    assert np.abs(rel[:, 1]).max() <= e2
    assert np.abs(rel[:, 2]).max() <= e3


@needs_dialog
def test_unchecking_auto_keeps_the_manual_extent():
    dlg = _dialog()
    _set_extents(dlg, 7.5, 2.5)
    assert _extents(dlg)[:2] == pytest.approx([7.5, 2.5])
    assert max(abs(p["a"]) for p in dlg._grid_points) == pytest.approx(7.5)
    assert max(abs(p["b"]) for p in dlg._grid_points) == pytest.approx(2.5)


@needs_dialog
def test_auto_toggle_enables_the_manual_spinbox_only_when_off():
    dlg = _dialog()
    dlg._on_auto_toggled(True)
    assert all(not row["e"].isEnabled() for row in dlg._axis_rows)
    dlg._on_auto_toggled(False)
    assert all(row["e"].isEnabled() for row in dlg._axis_rows)


@needs_dialog
def test_count_label_reports_the_2d_shape_and_spacing():
    dlg = _dialog()
    _set_counts(dlg, 5, 5)
    text = dlg._count_label.setText.call_args[0][0]
    assert "5 x 5" in text and "25" in text


@needs_dialog
def test_count_label_reports_the_3d_shape():
    dlg = _dialog()
    _select_mode(dlg, "3d")
    _set_counts(dlg, 4, 4, 3)
    text = dlg._count_label.setText.call_args[0][0]
    assert "4 x 4 x 3" in text and "48" in text


@needs_dialog
def test_large_grid_is_flagged_in_the_label():
    dlg = _dialog()
    _set_counts(dlg, 31, 31)
    assert len(dlg._grid_points) > _CONFIRM_ABOVE
    assert "large" in dlg._count_label.setText.call_args[0][0]


@needs_dialog
def test_grid_error_is_surfaced_not_swallowed(monkeypatch):
    monkeypatch.setattr(
        grid_mod,
        "compute_nics_grid",
        MagicMock(side_effect=ValueError("boom")),
    )
    dlg = _dialog()
    assert dlg._grid_points == []
    assert "boom" in dlg._count_label.setText.call_args[0][0]


@needs_dialog
def test_planarity_failure_leaves_the_row_blank(monkeypatch, caplog):
    import nics_placer.nics_math as math_mod

    monkeypatch.setattr(
        math_mod, "planarity_rms", MagicMock(side_effect=RuntimeError("nope"))
    )
    with caplog.at_level(logging.WARNING):
        dlg = _dialog()
    assert dlg._table.item(0, 3).text() == ""


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------


@needs_dialog
def test_preview_renders_one_glyph_actor():
    dlg = _dialog()
    names = [c.kwargs.get("name") for c in dlg._context.plotter.add_mesh.call_args_list]
    assert "nics_grid" in names


@needs_dialog
def test_preview_is_decimated_for_a_large_volume(monkeypatch):
    captured = {}

    class _Poly:
        def __init__(self, pts):
            captured["n"] = len(pts)

        def __setitem__(self, k, v):
            pass

        def glyph(self, **kw):
            return MagicMock()

    monkeypatch.setattr(grid_mod.pv, "PolyData", _Poly)
    dlg = _dialog()
    _select_mode(dlg, "3d")
    _set_counts(dlg, 21, 21, 21)
    assert len(dlg._grid_points) == 21**3
    assert captured["n"] <= _PREVIEW_MAX


@needs_dialog
def test_render_survives_a_missing_plotter():
    dlg = _dialog()
    dlg._context.plotter = None
    dlg._render_spheres()  # must not raise


@needs_dialog
def test_render_failure_is_logged_not_raised(caplog):
    dlg = _dialog()
    dlg._context.plotter.add_mesh.side_effect = RuntimeError("vtk sulked")
    with caplog.at_level(logging.WARNING):
        dlg._render_spheres()
    assert "vtk sulked" in caplog.text


@needs_dialog
def test_close_removes_the_preview_actor():
    dlg = _dialog()
    dlg.closeEvent(MagicMock())
    dlg._context.plotter.remove_actor.assert_any_call("nics_grid")


@needs_dialog
def test_clear_actors_failure_is_logged(caplog):
    dlg = _dialog()
    dlg._context.plotter.remove_actor.side_effect = RuntimeError("gone")
    with caplog.at_level(logging.WARNING):
        dlg._clear_actors()
    assert "gone" in caplog.text


# ---------------------------------------------------------------------------
# Placement
# ---------------------------------------------------------------------------


@needs_dialog
def test_placing_a_small_grid_adds_one_ghost_atom_per_probe():
    dlg = _dialog()
    _set_counts(dlg, 3, 3)
    before = dlg._context.current_molecule.GetNumAtoms()
    dlg._place_grid()
    assert dlg._context.current_molecule.GetNumAtoms() == before + 9
    dlg._context.push_undo_checkpoint.assert_called_once()


@needs_dialog
def test_placed_atoms_carry_the_ghost_label_and_the_probe_coordinates():
    dlg = _dialog()
    _set_counts(dlg, 3, 3)
    expected = [p["pos"] for p in dlg._grid_points]
    dlg._place_grid()
    mol = dlg._context.current_molecule
    conf = mol.GetConformer()
    ghosts = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    assert len(ghosts) == 9
    assert all(a.GetProp("custom_symbol") == "Bq" for a in ghosts)
    got = np.array([list(conf.GetAtomPosition(a.GetIdx())) for a in ghosts])
    np.testing.assert_allclose(got, np.array(expected), atol=1e-6)


@needs_dialog
def test_a_small_grid_is_placed_without_asking():
    dlg = _dialog()
    _set_counts(dlg, 3, 3)
    dlg._place_grid()
    assert grid_mod.QMessageBox._calls == []


@needs_dialog
def test_a_large_grid_asks_first_and_declining_places_nothing():
    dlg = _dialog()
    _set_counts(dlg, 31, 31)
    before = dlg._context.current_molecule.GetNumAtoms()
    dlg._place_grid()
    assert len(grid_mod.QMessageBox._calls) == 1
    assert dlg._context.current_molecule.GetNumAtoms() == before
    dlg._context.push_undo_checkpoint.assert_not_called()


@needs_dialog
def test_confirming_a_large_grid_places_it():
    grid_mod.QMessageBox._answer = grid_mod.QMessageBox.StandardButton.Yes
    dlg = _dialog()
    _set_counts(dlg, 21, 21)
    before = dlg._context.current_molecule.GetNumAtoms()
    dlg._place_grid()
    assert dlg._context.current_molecule.GetNumAtoms() == before + 441


@needs_dialog
def test_placing_an_empty_grid_is_a_no_op():
    dlg = _dialog()
    dlg._grid_points = []
    dlg._place_grid()
    dlg._context.push_undo_checkpoint.assert_not_called()


@needs_dialog
def test_placing_without_a_molecule_is_a_no_op():
    dlg = _dialog()
    dlg._context.current_molecule = None
    dlg._place_grid()
    dlg._context.push_undo_checkpoint.assert_not_called()


@needs_dialog
def test_clear_removes_every_ghost_atom():
    dlg = _dialog()
    _set_counts(dlg, 3, 3)
    dlg._place_grid()
    dlg._clear_all_bq()
    mol = dlg._context.current_molecule
    assert not any(a.GetAtomicNum() == 0 for a in mol.GetAtoms())


@needs_dialog
def test_clear_without_a_molecule_is_a_no_op():
    dlg = _dialog()
    dlg._context.current_molecule = None
    dlg._clear_all_bq()
    dlg._context.push_undo_checkpoint.assert_not_called()


@needs_dialog
def test_batch_placement_leaves_real_atoms_untouched():
    """The batch adder rebuilds the conformer once; the original coordinates
    must survive that rebuild intact."""
    dlg = _dialog()
    mol = dlg._context.current_molecule
    before = np.array(mol.GetConformer().GetPositions())
    _set_counts(dlg, 3, 3)
    dlg._place_grid()
    after = np.array(dlg._context.current_molecule.GetConformer().GetPositions())
    np.testing.assert_allclose(after[: len(before)], before, atol=1e-9)


# ---------------------------------------------------------------------------
# Ghost symbol
# ---------------------------------------------------------------------------


@needs_dialog
def test_symbol_change_is_persisted_and_used(monkeypatch):
    saved = {}
    monkeypatch.setattr(grid_mod, "_save_plugin_settings", lambda s: saved.update(s))
    dlg = _dialog()
    dlg._sym_combo.setCurrentIndex(1)
    dlg._on_symbol_changed(1)
    assert dlg._ghost_symbol == "H:"
    assert saved["ghost_symbol"] == "H:"
    _set_counts(dlg, 3, 3)
    dlg._place_grid()
    ghosts = [
        a for a in dlg._context.current_molecule.GetAtoms() if a.GetAtomicNum() == 0
    ]
    assert all(a.GetProp("custom_symbol") == "H:" for a in ghosts)


@needs_dialog
def test_symbol_syncs_from_settings_on_show(monkeypatch):
    monkeypatch.setitem(grid_mod._plugin_settings, "ghost_symbol", "H:")
    dlg = _dialog()
    assert dlg._ghost_symbol == "H:"


@needs_dialog
def test_an_unrecognised_stored_symbol_falls_back_to_bq(monkeypatch):
    monkeypatch.setitem(grid_mod._plugin_settings, "ghost_symbol", "Xx")
    dlg = _dialog()
    assert dlg._ghost_symbol == "Bq"


@needs_dialog
def test_unknown_plane_data_falls_back_to_parallel():
    dlg = _dialog()
    dlg._plane_combo.setCurrentIndex(99)
    assert dlg._current_plane() == "parallel"


# ---------------------------------------------------------------------------
# Rectangular grids, uniform spacing, and the offset on every plane
# ---------------------------------------------------------------------------


@needs_dialog
def test_a_rectangular_grid_is_allowed():
    """Molecules are rarely square; the two in-plane axes are independent."""
    dlg = _dialog()
    _set_extents(dlg, 6.0, 2.0)
    _set_counts(dlg, 7, 3)
    assert len(dlg._grid_points) == 21
    assert max(abs(p["a"]) for p in dlg._grid_points) == pytest.approx(6.0)
    assert max(abs(p["b"]) for p in dlg._grid_points) == pytest.approx(2.0)


@needs_dialog
def test_a_3d_box_can_differ_on_all_three_axes():
    dlg = _dialog()
    _select_mode(dlg, "3d")
    _set_extents(dlg, 6.0, 4.0, 2.0)
    _set_counts(dlg, 4, 3, 2)
    assert len(dlg._grid_points) == 24
    assert max(abs(p["c"]) for p in dlg._grid_points) == pytest.approx(2.0)


@needs_dialog
def test_auto_sizing_gives_a_rectangle_for_a_long_molecule():
    """Anthracene is far longer than it is wide, so the auto extents must not
    come out equal -- that was the whole point of per-axis sizing."""
    dlg = _dialog(smiles="c1ccc2cc3ccccc3cc2c1")
    _select_plane(dlg, "xy")
    e1, e2, _e3 = _extents(dlg)
    assert abs(e1 - e2) > 1.0


@needs_dialog
def test_uniform_spacing_matches_the_step_on_every_axis():
    dlg = _dialog()
    _select_mode(dlg, "3d")
    _set_extents(dlg, 6.0, 4.0, 2.0)
    dlg._uniform_spacing.setChecked(True)
    dlg._spacing_spin.setValue(1.0)
    dlg._on_uniform_toggled(True)
    n1, n2, n3 = _counts(dlg)
    assert (n1, n2, n3) == (13, 9, 5)
    for count, extent in zip((n1, n2, n3), (6.0, 4.0, 2.0)):
        assert (2.0 * extent) / (count - 1) == pytest.approx(1.0)


@needs_dialog
def test_uniform_spacing_locks_the_count_spinboxes():
    dlg = _dialog()
    dlg._uniform_spacing.setChecked(True)
    dlg._on_uniform_toggled(True)
    assert all(not row["n"].isEnabled() for row in dlg._axis_rows)
    assert dlg._spacing_spin.isEnabled()
    dlg._uniform_spacing.setChecked(False)
    dlg._on_uniform_toggled(False)
    assert all(row["n"].isEnabled() for row in dlg._axis_rows)


@needs_dialog
def test_a_finer_step_produces_more_probes():
    dlg = _dialog()
    _set_extents(dlg, 4.0, 4.0)
    dlg._uniform_spacing.setChecked(True)
    dlg._on_uniform_toggled(True)
    dlg._spacing_spin.setValue(1.0)
    dlg._on_params_changed()
    coarse = len(dlg._grid_points)
    dlg._spacing_spin.setValue(0.5)
    dlg._on_params_changed()
    assert len(dlg._grid_points) > coarse


@needs_dialog
def test_offset_defaults_to_zero():
    """NICS(0) -- the ring plane itself -- is the neutral starting point; the
    grid dialog should not silently start you 1 A off it."""
    dlg = _dialog()
    assert dlg._offset_spin.value() == pytest.approx(0.0)


@needs_dialog
@pytest.mark.parametrize(
    "plane,axis",
    [
        ("xy", np.array([0.0, 0.0, 1.0])),
        ("xz", np.array([0.0, -1.0, 0.0])),
        ("yz", np.array([1.0, 0.0, 0.0])),
    ],
)
def test_offset_shifts_a_lab_plane_along_its_own_normal(plane, axis):
    from nics_placer.nics_math import center_of_mass

    dlg = _dialog()
    _select_plane(dlg, plane)
    centre = center_of_mass(dlg._context.current_molecule)
    dlg._offset_spin.setValue(2.0)
    dlg._on_params_changed()
    for p in dlg._grid_points:
        assert float((p["pos"] - centre) @ axis) == pytest.approx(2.0)


@needs_dialog
def test_a_single_point_axis_is_accepted():
    dlg = _dialog()
    _set_counts(dlg, 1, 5)
    assert len(dlg._grid_points) == 5
    assert all(p["a"] == 0.0 for p in dlg._grid_points)


@needs_dialog
def test_count_label_names_the_axes_of_the_current_plane():
    dlg = _dialog()
    _select_plane(dlg, "xz")
    text = dlg._count_label.setText.call_args[0][0]
    assert "X" in text and "Z" in text


# ---------------------------------------------------------------------------
# Reset
# ---------------------------------------------------------------------------


@needs_dialog
def test_reset_restores_every_grid_control():
    dlg = _dialog()
    _select_mode(dlg, "3d")
    _select_plane(dlg, "yz")
    dlg._use_com.setChecked(False)
    dlg._uniform_spacing.setChecked(True)
    dlg._on_uniform_toggled(True)
    dlg._spacing_spin.setValue(0.2)
    dlg._offset_spin.setValue(4.0)
    _set_extents(dlg, 9.0, 8.0, 7.0)
    dlg._reset_settings()

    assert dlg._mode_combo.currentData() == "2d"
    assert dlg._current_plane() == "parallel"
    assert dlg._use_com.isChecked()
    assert dlg._auto_extent.isChecked()
    assert not dlg._uniform_spacing.isChecked()
    assert dlg._spacing_spin.value() == pytest.approx(grid_mod._DEFAULT_SPACING)
    assert dlg._offset_spin.value() == pytest.approx(0.0)
    assert _counts(dlg) == [grid_mod._DEFAULT_POINTS] * 3


@needs_dialog
def test_reset_leaves_the_ghost_label_alone():
    """The ghost symbol is a persisted user/project preference, not a grid
    parameter -- resetting the grid must not silently retag it."""
    dlg = _dialog()
    dlg._sym_combo.setCurrentIndex(1)
    dlg._on_symbol_changed(1)
    dlg._reset_settings()
    assert dlg._ghost_symbol == "H:"


@needs_dialog
def test_reset_re_enables_the_controls_it_may_have_locked():
    dlg = _dialog()
    dlg._uniform_spacing.setChecked(True)
    dlg._on_uniform_toggled(True)
    assert all(not row["n"].isEnabled() for row in dlg._axis_rows)
    dlg._reset_settings()
    assert all(row["n"].isEnabled() for row in dlg._axis_rows)
    assert all(not row["e"].isEnabled() for row in dlg._axis_rows)
    assert not dlg._spacing_spin.isEnabled()


@needs_dialog
def test_reset_rebuilds_the_grid_exactly_once():
    """Half-reset intermediate states can be very expensive grids (small step,
    large extent); the rebuild must happen after every widget has settled."""
    dlg = _dialog()
    _select_mode(dlg, "3d")
    dlg._context.plotter.add_mesh.reset_mock()
    dlg._reset_settings()
    assert dlg._context.plotter.add_mesh.call_count == 1
    assert len(dlg._grid_points) == grid_mod._DEFAULT_POINTS**2


@needs_dialog
def test_reset_leaves_a_usable_grid():
    dlg = _dialog()
    _set_counts(dlg, 2, 2)
    dlg._reset_settings()
    assert len(dlg._grid_points) > 0
    dlg._place_grid()
    assert dlg._context.push_undo_checkpoint.called
