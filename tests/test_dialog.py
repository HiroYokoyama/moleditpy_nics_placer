"""
Tests for nics_placer/dialog.py — Bq atom helpers and click-filter logic.
PyQt6 and pyvista are stubbed so tests run headlessly.
"""
import sys
import os
import types
from unittest.mock import MagicMock

import pytest
import numpy as np

# ---------------------------------------------------------------------------
# Qt and pyvista stubs — must install BEFORE importing nics_placer.dialog
# ---------------------------------------------------------------------------

def _install_stubs():
    class _QBase:
        def __init__(self, *args, **kwargs):
            pass

    class _Qt:
        class MouseButton:
            LeftButton = 1

        class WindowType:
            WindowMinMaxButtonsHint = 0
            WindowCloseButtonHint = 0

    _QEvent = MagicMock()
    _QEvent.Type.MouseButtonPress = 2
    _QEvent.Type.MouseButtonRelease = 3

    qt_core = types.ModuleType("PyQt6.QtCore")
    qt_core.Qt = _Qt
    qt_core.QTimer = MagicMock()
    qt_core.QObject = _QBase
    qt_core.QEvent = _QEvent

    qt_widgets = types.ModuleType("PyQt6.QtWidgets")
    qt_widgets.QDialog = _QBase
    qt_widgets.QAbstractItemView = MagicMock()
    qt_widgets.QComboBox = MagicMock()
    qt_widgets.QHBoxLayout = MagicMock()
    qt_widgets.QHeaderView = MagicMock()
    qt_widgets.QLabel = MagicMock()
    qt_widgets.QPushButton = MagicMock()
    qt_widgets.QTableWidget = MagicMock()
    qt_widgets.QTableWidgetItem = MagicMock()
    qt_widgets.QVBoxLayout = MagicMock()

    pyqt6 = types.ModuleType("PyQt6")
    pyqt6.QtCore = qt_core
    pyqt6.QtWidgets = qt_widgets
    for name, mod in [
        ("PyQt6", pyqt6),
        ("PyQt6.QtCore", qt_core),
        ("PyQt6.QtWidgets", qt_widgets),
    ]:
        sys.modules.setdefault(name, mod)
    sys.modules.setdefault("pyvista", MagicMock())

    return _Qt, _QEvent


_Qt_cls, _QEvent_cls = _install_stubs()

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), "..")))

_DIALOG_AVAILABLE = False
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from nics_placer.dialog import (
        _add_bq_atom,
        _remove_all_bq,
        _ClickFilter,
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


def _make_event(event_type, button_val, x, y):
    event = MagicMock()
    event.type.return_value = event_type
    event.button.return_value = button_val
    pt = MagicMock()
    pt.x.return_value = x
    pt.y.return_value = y
    event.position.return_value.toPoint.return_value = pt
    return event


# ---------------------------------------------------------------------------
# _add_bq_atom
# ---------------------------------------------------------------------------

@needs_dialog
def test_add_bq_atom_increases_atom_count():
    mol = _benzene_3d()
    new_mol = _add_bq_atom(mol, np.array([0.0, 0.0, 1.5]))
    assert new_mol.GetNumAtoms() == mol.GetNumAtoms() + 1


@needs_dialog
def test_add_bq_atom_sets_custom_symbol_bq():
    mol = _benzene_3d()
    new_mol = _add_bq_atom(mol, np.array([0.0, 0.0, 1.5]), symbol="Bq")
    added = new_mol.GetAtomWithIdx(new_mol.GetNumAtoms() - 1)
    assert added.HasProp("custom_symbol")
    assert added.GetProp("custom_symbol") == "Bq"


@needs_dialog
def test_add_bq_atom_h_colon_symbol():
    mol = _benzene_3d()
    new_mol = _add_bq_atom(mol, np.array([0.0, 0.0, 1.5]), symbol="H:")
    added = new_mol.GetAtomWithIdx(new_mol.GetNumAtoms() - 1)
    assert added.GetProp("custom_symbol") == "H:"


@needs_dialog
def test_add_bq_atom_is_dummy_atom():
    mol = _benzene_3d()
    new_mol = _add_bq_atom(mol, np.array([0.0, 0.0, 2.0]))
    added = new_mol.GetAtomWithIdx(new_mol.GetNumAtoms() - 1)
    assert added.GetAtomicNum() == 0


@needs_dialog
def test_add_bq_atom_position_stored_correctly():
    mol = _benzene_3d()
    pos = np.array([1.23, 4.56, 7.89])
    new_mol = _add_bq_atom(mol, pos)
    last_idx = new_mol.GetNumAtoms() - 1
    p = new_mol.GetConformer().GetAtomPosition(last_idx)
    np.testing.assert_allclose([p.x, p.y, p.z], pos, atol=1e-5)


@needs_dialog
def test_add_bq_atom_does_not_mutate_original():
    mol = _benzene_3d()
    n = mol.GetNumAtoms()
    _add_bq_atom(mol, np.array([0.0, 0.0, 1.5]))
    assert mol.GetNumAtoms() == n


@needs_dialog
def test_add_two_bq_atoms_accumulate():
    mol = _benzene_3d()
    n = mol.GetNumAtoms()
    mol = _add_bq_atom(mol, np.array([0.0, 0.0, 1.0]))
    mol = _add_bq_atom(mol, np.array([0.0, 0.0, -1.0]))
    assert mol.GetNumAtoms() == n + 2


# ---------------------------------------------------------------------------
# _remove_all_bq
# ---------------------------------------------------------------------------

@needs_dialog
def test_remove_all_bq_removes_ghost_atoms():
    mol = _benzene_3d()
    n = mol.GetNumAtoms()
    mol = _add_bq_atom(mol, np.array([0.0, 0.0, 1.5]))
    mol = _add_bq_atom(mol, np.array([0.0, 0.0, -1.5]))
    cleaned = _remove_all_bq(mol)
    assert cleaned.GetNumAtoms() == n


@needs_dialog
def test_remove_all_bq_no_ghosts_is_idempotent():
    mol = _benzene_3d()
    n = mol.GetNumAtoms()
    cleaned = _remove_all_bq(mol)
    assert cleaned.GetNumAtoms() == n


@needs_dialog
def test_remove_all_bq_removes_h_colon_atoms():
    mol = _benzene_3d()
    n = mol.GetNumAtoms()
    mol = _add_bq_atom(mol, np.array([0.0, 0.0, 1.5]), symbol="H:")
    cleaned = _remove_all_bq(mol)
    assert cleaned.GetNumAtoms() == n


@needs_dialog
def test_remove_then_add_roundtrip():
    mol = _benzene_3d()
    n_orig = mol.GetNumAtoms()
    mol = _add_bq_atom(mol, np.array([0.0, 0.0, 1.5]))
    mol = _remove_all_bq(mol)
    assert mol.GetNumAtoms() == n_orig


@needs_dialog
def test_remove_all_bq_no_ghost_custom_symbol_remains():
    """Atoms without custom_symbol property are left untouched by _remove_all_bq."""
    mol = _benzene_3d()
    n = mol.GetNumAtoms()
    cleaned = _remove_all_bq(mol)
    assert cleaned.GetNumAtoms() == n


# ---------------------------------------------------------------------------
# _ClickFilter
# ---------------------------------------------------------------------------

@needs_dialog
def test_click_filter_press_position_starts_none():
    f = _ClickFilter(MagicMock())
    assert f._press_pos is None


@needs_dialog
def test_click_filter_records_press_on_left_button():
    callback = MagicMock()
    f = _ClickFilter(callback)
    event = _make_event(_QEvent_cls.Type.MouseButtonPress, _Qt_cls.MouseButton.LeftButton, 50, 80)
    f.eventFilter(MagicMock(), event)
    assert f._press_pos is not None


@needs_dialog
def test_click_filter_eventfilter_returns_false():
    f = _ClickFilter(MagicMock())
    event = _make_event(_QEvent_cls.Type.MouseButtonPress, _Qt_cls.MouseButton.LeftButton, 10, 10)
    assert f.eventFilter(MagicMock(), event) is False


@needs_dialog
def test_click_filter_fires_callback_on_small_drag():
    callback = MagicMock()
    f = _ClickFilter(callback)
    obj = MagicMock()
    # Press at (100, 100)
    f.eventFilter(obj, _make_event(_QEvent_cls.Type.MouseButtonPress, _Qt_cls.MouseButton.LeftButton, 100, 100))
    # Release at (101, 100) — dx=1, dy=0 → distance² = 1 ≤ 25
    f.eventFilter(obj, _make_event(_QEvent_cls.Type.MouseButtonRelease, _Qt_cls.MouseButton.LeftButton, 101, 100))
    callback.assert_called_once_with(101, 100, obj)


@needs_dialog
def test_click_filter_does_not_fire_on_large_drag():
    callback = MagicMock()
    f = _ClickFilter(callback)
    obj = MagicMock()
    # Press at (0, 0), release at (20, 20) → distance² = 800 > 25
    f.eventFilter(obj, _make_event(_QEvent_cls.Type.MouseButtonPress, _Qt_cls.MouseButton.LeftButton, 0, 0))
    f.eventFilter(obj, _make_event(_QEvent_cls.Type.MouseButtonRelease, _Qt_cls.MouseButton.LeftButton, 20, 20))
    callback.assert_not_called()


@needs_dialog
def test_click_filter_resets_press_pos_after_release():
    f = _ClickFilter(MagicMock())
    obj = MagicMock()
    f.eventFilter(obj, _make_event(_QEvent_cls.Type.MouseButtonPress, _Qt_cls.MouseButton.LeftButton, 10, 10))
    f.eventFilter(obj, _make_event(_QEvent_cls.Type.MouseButtonRelease, _Qt_cls.MouseButton.LeftButton, 10, 10))
    assert f._press_pos is None


# ---------------------------------------------------------------------------
# State constants
# ---------------------------------------------------------------------------

@needs_dialog
def test_state_constants_are_distinct():
    assert _STATE_UNSET != _STATE_STAGED
    assert _STATE_STAGED != _STATE_PLACED
    assert _STATE_UNSET != _STATE_PLACED
