"""
Tests for NicsPlacerDialog itself — instantiated headlessly against the
rich PyQt6 / pyvista stubs installed by tests/conftest.py.
"""
import sys
import os
import types
import logging
from unittest.mock import MagicMock

import pytest
import numpy as np

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), "..")))

_DIALOG_AVAILABLE = False
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import nics_placer.dialog as dialog_mod
    from nics_placer.dialog import (
        NicsPlacerDialog,
        _STATE_UNSET,
        _STATE_STAGED,
        _STATE_PLACED,
    )
    _DIALOG_AVAILABLE = Chem is not None
except Exception:
    pass

needs_dialog = pytest.mark.skipif(
    not _DIALOG_AVAILABLE, reason="RDKit or nics_placer.dialog not importable"
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _benzene_3d():
    mol = Chem.MolFromSmiles("c1ccccc1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


class _StubContext:
    def __init__(self, mol=None):
        self.current_molecule = mol
        self.plotter = MagicMock(name="plotter")
        self.plotter.add_mesh.side_effect = (
            lambda mesh, name=None, **kw: MagicMock(name=f"actor_{name}")
        )
        self._windows = {}
        self.push_undo_checkpoint = MagicMock()
        self.show_status_message = MagicMock()

    def get_window(self, key):
        return self._windows.get(key)

    def register_window(self, key, win):
        self._windows[key] = win

    def get_main_window(self):
        return MagicMock()


class _FakeVtkModule:
    """Swapped into sys.modules['vtk'] for _on_plotter_click tests."""
    picked_actor = None
    pick_position = (0.0, 0.0, 0.0)

    class vtkCellPicker:
        def __init__(self):
            pass

        def SetTolerance(self, t):
            pass

        def Pick(self, x, y, z, renderer):
            pass

        def GetActor(self):
            return _FakeVtkModule.picked_actor

        def GetPickPosition(self):
            return _FakeVtkModule.pick_position


@pytest.fixture
def fake_vtk():
    saved = sys.modules.get("vtk")
    sys.modules["vtk"] = _FakeVtkModule
    yield _FakeVtkModule
    _FakeVtkModule.picked_actor = None
    _FakeVtkModule.pick_position = (0.0, 0.0, 0.0)
    if saved is not None:
        sys.modules["vtk"] = saved
    else:
        del sys.modules["vtk"]


def _make_widget(ratio=1.0, height=100):
    w = MagicMock()
    w.devicePixelRatioF.return_value = ratio
    w.height.return_value = height
    return w


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

@needs_dialog
def test_dialog_construction_does_not_raise():
    ctx = _StubContext()
    dlg = NicsPlacerDialog(ctx)
    assert dlg._context is ctx


@needs_dialog
def test_dialog_construction_initial_state():
    ctx = _StubContext()
    dlg = NicsPlacerDialog(ctx)
    assert dlg._nics_points == []
    assert dlg._ghost_symbol == "Bq"
    assert dlg._click_filter is None


# ---------------------------------------------------------------------------
# _load_rings
# ---------------------------------------------------------------------------

@needs_dialog
def test_load_rings_no_molecule():
    ctx = _StubContext(mol=None)
    dlg = NicsPlacerDialog(ctx)
    dlg._load_rings()
    assert dlg._nics_points == []
    assert dlg._table.rowCount() == 0


@needs_dialog
def test_load_rings_no_conformer():
    ctx = _StubContext(mol=Chem.MolFromSmiles("c1ccccc1"))
    dlg = NicsPlacerDialog(ctx)
    dlg._load_rings()
    assert dlg._nics_points == []


@needs_dialog
def test_load_rings_benzene_populates_table_and_points():
    mol = _benzene_3d()
    ctx = _StubContext(mol=mol)
    dlg = NicsPlacerDialog(ctx)
    dlg._load_rings()
    assert dlg._table.rowCount() == 1
    assert len(dlg._nics_points) == 3
    types_present = {p["type"] for p in dlg._nics_points}
    assert types_present == {"nics0", "nics1_above", "nics1_below"}
    status_item = dlg._table.item(0, 3)
    assert status_item.text() == "0/3 placed"


@needs_dialog
def test_load_rings_ring_calc_exception_sets_error_status(monkeypatch):
    mol = _benzene_3d()
    ctx = _StubContext(mol=mol)
    dlg = NicsPlacerDialog(ctx)

    def _boom(*a, **kw):
        raise ValueError("boom")

    monkeypatch.setattr(dialog_mod, "get_ring_positions", _boom)
    dlg._load_rings()
    # regression: _sync_placed_status()/_update_table_status() must not clobber
    # the "error" cell just because the errored ring has zero staged points.
    assert dlg._table.item(0, 3).text() == "error"
    assert dlg._nics_points == []


# ---------------------------------------------------------------------------
# _sync_placed_status / _update_table_status
# ---------------------------------------------------------------------------

@needs_dialog
def test_sync_placed_status_marks_matching_bq_as_placed():
    mol = _benzene_3d()
    ctx = _StubContext(mol=mol)
    dlg = NicsPlacerDialog(ctx)
    dlg._load_rings()
    nics0_pt = next(p for p in dlg._nics_points if p["type"] == "nics0")

    mol2 = dialog_mod._add_bq_atom(mol, nics0_pt["pos"])
    ctx.current_molecule = mol2
    dlg._sync_placed_status()
    assert nics0_pt["state"] == _STATE_PLACED


@needs_dialog
def test_sync_placed_status_no_molecule_is_noop():
    ctx = _StubContext(mol=None)
    dlg = NicsPlacerDialog(ctx)
    dlg._sync_placed_status()  # should not raise


@needs_dialog
def test_update_table_status_no_points_shows_dash():
    mol = _benzene_3d()
    ctx = _StubContext(mol=mol)
    dlg = NicsPlacerDialog(ctx)
    dlg._table.setRowCount(1)
    dlg._nics_points = []
    dlg._update_table_status()
    assert dlg._table.item(0, 3).text() == "—"


@needs_dialog
def test_update_table_status_all_placed():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._table.setRowCount(1)
    dlg._nics_points = [
        {"ring": 0, "type": "nics0", "pos": np.zeros(3), "state": _STATE_PLACED}
    ]
    dlg._update_table_status()
    assert "all placed" in dlg._table.item(0, 3).text()


@needs_dialog
def test_update_table_status_staged_mix():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._table.setRowCount(1)
    dlg._nics_points = [
        {"ring": 0, "type": "nics0", "pos": np.zeros(3), "state": _STATE_PLACED},
        {"ring": 0, "type": "nics1_above", "pos": np.zeros(3), "state": _STATE_STAGED},
    ]
    dlg._update_table_status()
    assert "staged" in dlg._table.item(0, 3).text()


@needs_dialog
def test_update_table_status_plain_fraction():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._table.setRowCount(1)
    dlg._nics_points = [
        {"ring": 0, "type": "nics0", "pos": np.zeros(3), "state": _STATE_UNSET},
    ]
    dlg._update_table_status()
    assert dlg._table.item(0, 3).text() == "0/1 placed"


# ---------------------------------------------------------------------------
# _render_spheres / _clear_actors
# ---------------------------------------------------------------------------

@needs_dialog
def test_render_spheres_no_plotter_is_noop():
    ctx = _StubContext(mol=_benzene_3d())
    ctx.plotter = None
    dlg = NicsPlacerDialog(ctx)
    dlg._render_spheres()  # should not raise


@needs_dialog
def test_render_spheres_creates_actors_for_each_state():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._nics_points = [
        {"ring": 0, "type": "nics0", "pos": np.array([0.0, 0.0, 0.0]), "state": _STATE_UNSET},
        {"ring": 0, "type": "nics1_above", "pos": np.array([0.0, 0.0, 1.0]), "state": _STATE_STAGED},
        {"ring": 0, "type": "nics1_below", "pos": np.array([0.0, 0.0, -1.0]), "state": _STATE_PLACED},
    ]
    dlg._render_spheres()
    assert dlg._actor_yellow is not None
    assert dlg._actor_red is not None
    assert dlg._actor_green is not None
    assert ctx.plotter.render.called


@needs_dialog
def test_render_spheres_swallows_exception(monkeypatch):
    ctx = _StubContext(mol=_benzene_3d())
    ctx.plotter.remove_actor.side_effect = RuntimeError("boom")
    dlg = NicsPlacerDialog(ctx)
    dlg._render_spheres()  # should not raise


@needs_dialog
def test_clear_actors_resets_refs_and_removes():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._actor_yellow = MagicMock()
    dlg._clear_actors()
    assert dlg._actor_yellow is None
    assert dlg._actor_red is None
    assert dlg._actor_green is None
    assert ctx.plotter.render.called


@needs_dialog
def test_clear_actors_swallows_exception():
    ctx = _StubContext(mol=_benzene_3d())
    ctx.plotter.remove_actor.side_effect = RuntimeError("boom")
    dlg = NicsPlacerDialog(ctx)
    dlg._clear_actors()  # should not raise
    assert dlg._actor_yellow is None


# ---------------------------------------------------------------------------
# _enable_picking / _disable_picking
# ---------------------------------------------------------------------------

@needs_dialog
def test_enable_picking_installs_filter():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._enable_picking()
    assert dlg._click_filter is not None
    ctx.plotter.installEventFilter.assert_called_once()


@needs_dialog
def test_enable_picking_idempotent():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._enable_picking()
    first = dlg._click_filter
    dlg._enable_picking()
    assert dlg._click_filter is first
    ctx.plotter.installEventFilter.assert_called_once()


@needs_dialog
def test_enable_picking_no_plotter():
    ctx = _StubContext(mol=_benzene_3d())
    ctx.plotter = None
    dlg = NicsPlacerDialog(ctx)
    dlg._enable_picking()
    assert dlg._click_filter is None


@needs_dialog
def test_disable_picking_removes_filter():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._enable_picking()
    dlg._disable_picking()
    assert dlg._click_filter is None
    ctx.plotter.removeEventFilter.assert_called_once()


@needs_dialog
def test_disable_picking_without_filter_is_safe():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._disable_picking()  # should not raise
    assert dlg._click_filter is None


@needs_dialog
def test_disable_picking_swallows_exception():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._enable_picking()
    ctx.plotter.removeEventFilter.side_effect = RuntimeError("boom")
    dlg._disable_picking()  # should not raise
    assert dlg._click_filter is None


# ---------------------------------------------------------------------------
# _on_plotter_click
# ---------------------------------------------------------------------------

@needs_dialog
def test_on_plotter_click_no_plotter_is_noop(fake_vtk):
    ctx = _StubContext(mol=_benzene_3d())
    ctx.plotter = None
    dlg = NicsPlacerDialog(ctx)
    dlg._on_plotter_click(1, 1, _make_widget())  # should not raise


@needs_dialog
def test_on_plotter_click_missing_pick_ignored_even_with_no_red_actor(fake_vtk):
    """Regression: a miss (picked is None) must not be treated as a hit on
    the (also None) red actor just because there are no staged points yet."""
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._load_rings()
    dlg._render_spheres()
    assert dlg._actor_red is None  # nothing staged -> no red actor rendered
    fake_vtk.picked_actor = None
    fake_vtk.pick_position = tuple(
        float(v) for v in next(p for p in dlg._nics_points if p["type"] == "nics0")["pos"]
    )
    dlg._on_plotter_click(1, 1, _make_widget())
    assert all(p["state"] == _STATE_UNSET for p in dlg._nics_points)


@needs_dialog
def test_on_plotter_click_toggles_nearby_point(fake_vtk):
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._load_rings()
    dlg._render_spheres()
    nics0_pt = next(p for p in dlg._nics_points if p["type"] == "nics0")

    fake_vtk.picked_actor = dlg._actor_yellow
    fake_vtk.pick_position = tuple(float(v) for v in nics0_pt["pos"])
    dlg._on_plotter_click(5, 5, _make_widget())
    assert nics0_pt["state"] == _STATE_STAGED

    # click again -> toggles back to unset
    dlg._render_spheres()
    fake_vtk.picked_actor = dlg._actor_red
    dlg._on_plotter_click(5, 5, _make_widget())
    assert nics0_pt["state"] == _STATE_UNSET


@needs_dialog
def test_on_plotter_click_too_far_does_not_toggle(fake_vtk):
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._load_rings()
    dlg._render_spheres()
    fake_vtk.picked_actor = dlg._actor_yellow
    fake_vtk.pick_position = (999.0, 999.0, 999.0)
    dlg._on_plotter_click(5, 5, _make_widget())
    assert all(p["state"] == _STATE_UNSET for p in dlg._nics_points)


@needs_dialog
def test_on_plotter_click_skips_placed_points(fake_vtk):
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._load_rings()
    for p in dlg._nics_points:
        p["state"] = _STATE_PLACED
    dlg._render_spheres()
    fake_vtk.picked_actor = MagicMock()  # not yellow/red -> ignored anyway
    dlg._on_plotter_click(5, 5, _make_widget())
    assert all(p["state"] == _STATE_PLACED for p in dlg._nics_points)


@needs_dialog
def test_on_plotter_click_swallows_exception(fake_vtk):
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._load_rings()
    dlg._render_spheres()
    widget = _make_widget()
    widget.devicePixelRatioF.side_effect = RuntimeError("boom")
    dlg._on_plotter_click(1, 1, widget)  # should not raise


# ---------------------------------------------------------------------------
# Symbol combo sync
# ---------------------------------------------------------------------------

@needs_dialog
def test_sync_symbol_from_settings_default_bq():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dialog_mod._plugin_settings["ghost_symbol"] = "Bq"
    dlg.sync_symbol_from_settings()
    assert dlg._ghost_symbol == "Bq"
    assert dlg._sym_combo.currentData() == "Bq"


@needs_dialog
def test_sync_symbol_from_settings_h_colon():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dialog_mod._plugin_settings["ghost_symbol"] = "H:"
    dlg.sync_symbol_from_settings()
    assert dlg._ghost_symbol == "H:"
    assert dlg._sym_combo.currentData() == "H:"
    dialog_mod._plugin_settings["ghost_symbol"] = "Bq"  # restore


@needs_dialog
def test_on_symbol_changed_saves_settings_and_retags(monkeypatch):
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._sym_combo.setCurrentIndex(1)  # "H:"

    saved = {}
    monkeypatch.setattr(dialog_mod, "_save_plugin_settings", lambda s: saved.update(s))
    called = {}
    monkeypatch.setattr(dlg, "_retag_placed_atoms", lambda sym: called.setdefault("sym", sym))

    dlg._on_symbol_changed(1)
    assert dlg._ghost_symbol == "H:"
    assert saved.get("ghost_symbol") == "H:"
    assert called.get("sym") == "H:"
    dialog_mod._plugin_settings["ghost_symbol"] = "Bq"  # restore module-level state


# ---------------------------------------------------------------------------
# _retag_placed_atoms / _retag_bare_dummy_atoms
# ---------------------------------------------------------------------------

@needs_dialog
def test_retag_placed_atoms_no_molecule_is_noop():
    ctx = _StubContext(mol=None)
    dlg = NicsPlacerDialog(ctx)
    dlg._retag_placed_atoms("H:")  # should not raise
    assert not ctx.push_undo_checkpoint.called


@needs_dialog
def test_retag_placed_atoms_relabels_dummy_atoms():
    mol = _benzene_3d()
    mol = dialog_mod._add_bq_atom(mol, np.array([0.0, 0.0, 1.5]), symbol="Bq")
    ctx = _StubContext(mol=mol)
    dlg = NicsPlacerDialog(ctx)
    dlg._retag_placed_atoms("H:")
    dummy = ctx.current_molecule.GetAtomWithIdx(mol.GetNumAtoms() - 1)
    assert dummy.GetProp("custom_symbol") == "H:"
    assert ctx.push_undo_checkpoint.called


@needs_dialog
def test_retag_placed_atoms_no_change_skips_checkpoint():
    mol = _benzene_3d()
    mol = dialog_mod._add_bq_atom(mol, np.array([0.0, 0.0, 1.5]), symbol="Bq")
    ctx = _StubContext(mol=mol)
    dlg = NicsPlacerDialog(ctx)
    dlg._retag_placed_atoms("Bq")  # already Bq -> no change
    assert not ctx.push_undo_checkpoint.called


@needs_dialog
def test_retag_bare_dummy_atoms_no_molecule_is_noop():
    ctx = _StubContext(mol=None)
    dlg = NicsPlacerDialog(ctx)
    dlg._retag_bare_dummy_atoms()  # should not raise


@needs_dialog
def test_retag_bare_dummy_atoms_stamps_bare_dummies():
    mol = _benzene_3d()
    rw = Chem.RWMol(mol)
    atom = Chem.Atom(0)  # bare dummy, no custom_symbol
    idx = rw.AddAtom(atom)
    conf = rw.GetConformer()
    from rdkit.Geometry import Point3D
    conf.SetAtomPosition(idx, Point3D(0.0, 0.0, 2.0))
    mol2 = rw.GetMol()

    ctx = _StubContext(mol=mol2)
    dlg = NicsPlacerDialog(ctx)
    dlg._ghost_symbol = "Bq"
    dlg._retag_bare_dummy_atoms()
    dummy = ctx.current_molecule.GetAtomWithIdx(idx)
    assert dummy.GetProp("custom_symbol") == "Bq"


# ---------------------------------------------------------------------------
# _selected_ring_indices / _stage_rings
# ---------------------------------------------------------------------------

@needs_dialog
def test_selected_ring_indices_none_selected_returns_all_rows():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._table.setRowCount(3)
    assert dlg._selected_ring_indices() == {0, 1, 2}


@needs_dialog
def test_selected_ring_indices_returns_selected_rows():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._table.setRowCount(3)
    dlg._table._set_selected_rows([1])
    assert dlg._selected_ring_indices() == {1}


@needs_dialog
def test_stage_selected_nics0():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._load_rings()
    dlg._stage_selected_nics0()
    nics0_states = {p["state"] for p in dlg._nics_points if p["type"] == "nics0"}
    other_states = {p["state"] for p in dlg._nics_points if p["type"] != "nics0"}
    assert nics0_states == {_STATE_STAGED}
    assert other_states == {_STATE_UNSET}


@needs_dialog
def test_stage_selected_nics1():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._load_rings()
    dlg._stage_selected_nics1()
    for p in dlg._nics_points:
        if p["type"] in ("nics1_above", "nics1_below"):
            assert p["state"] == _STATE_STAGED
        else:
            assert p["state"] == _STATE_UNSET


@needs_dialog
def test_stage_rings_noop_when_nothing_changes():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._load_rings()
    for p in dlg._nics_points:
        p["state"] = _STATE_PLACED
    dlg._stage_rings({0}, {"nics0"})  # nothing UNSET -> no-op branch
    assert all(p["state"] == _STATE_PLACED for p in dlg._nics_points)


# ---------------------------------------------------------------------------
# Placement: _apply_staged / _place_all / _clear_all_bq
# ---------------------------------------------------------------------------

@needs_dialog
def test_apply_staged_noop_when_nothing_staged():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._apply_staged()
    assert not ctx.push_undo_checkpoint.called


@needs_dialog
def test_apply_staged_noop_when_no_molecule():
    ctx = _StubContext(mol=None)
    dlg = NicsPlacerDialog(ctx)
    dlg._nics_points = [
        {"ring": 0, "type": "nics0", "pos": np.zeros(3), "state": _STATE_STAGED}
    ]
    dlg._apply_staged()
    assert not ctx.push_undo_checkpoint.called


@needs_dialog
def test_apply_staged_places_bq_atoms():
    mol = _benzene_3d()
    n0 = mol.GetNumAtoms()
    ctx = _StubContext(mol=mol)
    dlg = NicsPlacerDialog(ctx)
    dlg._load_rings()
    dlg._stage_selected_nics0()
    dlg._apply_staged()

    assert ctx.push_undo_checkpoint.called
    assert ctx.current_molecule.GetNumAtoms() == n0 + 1
    nics0_pt = next(p for p in dlg._nics_points if p["type"] == "nics0")
    assert nics0_pt["state"] == _STATE_PLACED


@needs_dialog
def test_place_all_stages_and_applies_everything():
    mol = _benzene_3d()
    n0 = mol.GetNumAtoms()
    ctx = _StubContext(mol=mol)
    dlg = NicsPlacerDialog(ctx)
    dlg._load_rings()
    dlg._place_all()
    assert ctx.current_molecule.GetNumAtoms() == n0 + 3
    assert all(p["state"] == _STATE_PLACED for p in dlg._nics_points)


@needs_dialog
def test_clear_all_bq_noop_when_no_molecule():
    ctx = _StubContext(mol=None)
    dlg = NicsPlacerDialog(ctx)
    dlg._clear_all_bq()  # should not raise
    assert not ctx.push_undo_checkpoint.called


@needs_dialog
def test_clear_all_bq_removes_placed_atoms():
    mol = _benzene_3d()
    n0 = mol.GetNumAtoms()
    ctx = _StubContext(mol=mol)
    dlg = NicsPlacerDialog(ctx)
    dlg._load_rings()
    dlg._place_all()
    dlg._clear_all_bq()
    assert ctx.current_molecule.GetNumAtoms() == n0
    assert all(p["state"] == _STATE_UNSET for p in dlg._nics_points)
    assert ctx.push_undo_checkpoint.call_count == 2  # once for place, once for clear


# ---------------------------------------------------------------------------
# _check_molecule_changed
# ---------------------------------------------------------------------------

@needs_dialog
def test_check_molecule_changed_same_id_no_reload(monkeypatch):
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._load_rings()
    called = {"n": 0}
    monkeypatch.setattr(dlg, "_load_rings", lambda: called.__setitem__("n", called["n"] + 1))
    dlg._check_molecule_changed()
    assert called["n"] == 0


@needs_dialog
def test_check_molecule_changed_new_mol_triggers_reload(monkeypatch):
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    called = {"n": 0}
    monkeypatch.setattr(dlg, "_load_rings", lambda: called.__setitem__("n", called["n"] + 1))
    ctx.current_molecule = _benzene_3d()  # new object -> different id()
    dlg._check_molecule_changed()
    assert called["n"] == 1


class _ToggleRaiseContext(_StubContext):
    """current_molecule reads normally until .boom is armed, then raises."""

    def __init__(self, mol=None):
        self._mol = mol
        self.boom = False
        super().__init__(mol=mol)

    @property
    def current_molecule(self):
        if self.boom:
            raise RuntimeError("boom")
        return self._mol

    @current_molecule.setter
    def current_molecule(self, value):
        self._mol = value


@needs_dialog
def test_check_molecule_changed_swallows_exception():
    ctx = _ToggleRaiseContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    ctx.boom = True
    dlg._check_molecule_changed()  # should not raise


# ---------------------------------------------------------------------------
# showEvent / closeEvent
# ---------------------------------------------------------------------------

@needs_dialog
def test_show_event_full_flow():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg.showEvent(MagicMock())
    assert dlg._table.rowCount() == 1
    assert dlg._click_filter is not None
    assert dlg._poll_timer.isActive() is True


@needs_dialog
def test_show_event_does_not_restart_active_timer():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg._poll_timer.start()
    dlg.showEvent(MagicMock())
    assert dlg._poll_timer.isActive() is True


@needs_dialog
def test_close_event_stops_timer_and_disables_picking():
    ctx = _StubContext(mol=_benzene_3d())
    dlg = NicsPlacerDialog(ctx)
    dlg.showEvent(MagicMock())
    assert dlg._poll_timer.isActive() is True
    dlg.closeEvent(MagicMock())
    assert dlg._poll_timer.isActive() is False
    assert dlg._click_filter is None
    assert dlg._actor_yellow is None
