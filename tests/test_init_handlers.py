"""
Tests for nics_placer/__init__.py — show_dialog, settings load/save error
branches, and the deferred vs. immediate _apply_labels path in on_load.
Relies on the rich PyQt6/pyvista stubs installed by tests/conftest.py so
that `from .dialog import NicsPlacerDialog` succeeds headlessly.
"""
import builtins
import importlib
import json
import os
import sys
from unittest.mock import MagicMock, mock_open, patch

import pytest

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), "..")))

import nics_placer as pkg
from nics_placer import initialize

try:
    import rdkit  # noqa: F401
    _RDKIT_AVAILABLE = True
except Exception:
    _RDKIT_AVAILABLE = False

needs_rdkit = pytest.mark.skipif(
    not _RDKIT_AVAILABLE, reason="RDKit not importable (dialog import requires it)"
)


class _StubContext:
    def __init__(self):
        self._windows = {}
        self._load_handler = None
        self.current_molecule = None
        self.plotter = None
        self.add_menu_action = MagicMock()
        self.show_status_message = MagicMock()

    def get_main_window(self):
        return MagicMock()

    def get_window(self, key):
        return self._windows.get(key)

    def register_window(self, key, win):
        self._windows[key] = win

    def register_save_handler(self, fn):
        pass

    def register_load_handler(self, fn):
        self._load_handler = fn

    def register_document_reset_handler(self, fn):
        pass

    def load(self, data):
        self._load_handler(data)


def _get_show_dialog(ctx):
    initialize(ctx)
    return ctx.add_menu_action.call_args[0][1]


@pytest.fixture(autouse=True)
def _reset_dialog_opened():
    pkg._dialog_opened = False
    yield
    pkg._dialog_opened = False


# ---------------------------------------------------------------------------
# show_dialog
# ---------------------------------------------------------------------------

def test_show_dialog_reuses_existing_window():
    ctx = _StubContext()
    win = MagicMock()
    ctx._windows["main_panel"] = win
    show_dialog = _get_show_dialog(ctx)
    show_dialog()
    win.show.assert_called_once()
    win.raise_.assert_called_once()
    win.activateWindow.assert_called_once()


def test_show_dialog_no_molecule_shows_status():
    ctx = _StubContext()
    ctx.current_molecule = None
    show_dialog = _get_show_dialog(ctx)
    show_dialog()
    ctx.show_status_message.assert_called_once()
    assert "No molecule" in ctx.show_status_message.call_args[0][0]
    assert ctx.get_window("main_panel") is None


def test_show_dialog_no_conformers_shows_status():
    ctx = _StubContext()
    mol = MagicMock()
    mol.GetNumConformers.return_value = 0
    ctx.current_molecule = mol
    show_dialog = _get_show_dialog(ctx)
    show_dialog()
    ctx.show_status_message.assert_called_once()
    assert "3D conversion" in ctx.show_status_message.call_args[0][0]


@needs_rdkit
def test_show_dialog_creates_and_registers_dialog():
    ctx = _StubContext()
    mol = MagicMock()
    mol.GetNumConformers.return_value = 1
    ctx.current_molecule = mol
    show_dialog = _get_show_dialog(ctx)
    show_dialog()
    assert ctx.get_window("main_panel") is not None
    assert pkg._dialog_opened is True


# ---------------------------------------------------------------------------
# _load_plugin_settings / _save_plugin_settings error branches
# ---------------------------------------------------------------------------

def test_load_plugin_settings_file_not_found(tmp_path, monkeypatch):
    missing = tmp_path / "does_not_exist.json"
    monkeypatch.setattr(pkg, "_SETTINGS_FILE", str(missing))
    result = pkg._load_plugin_settings()
    assert result == {"ghost_symbol": "Bq"}


def test_load_plugin_settings_invalid_json(tmp_path, monkeypatch):
    bad = tmp_path / "settings.json"
    bad.write_text("{not valid json", encoding="utf-8")
    monkeypatch.setattr(pkg, "_SETTINGS_FILE", str(bad))
    result = pkg._load_plugin_settings()
    assert result == {"ghost_symbol": "Bq"}


def test_load_plugin_settings_non_dict_json(tmp_path, monkeypatch):
    f = tmp_path / "settings.json"
    f.write_text("[1, 2, 3]", encoding="utf-8")
    monkeypatch.setattr(pkg, "_SETTINGS_FILE", str(f))
    result = pkg._load_plugin_settings()
    assert result == {"ghost_symbol": "Bq"}


def test_load_plugin_settings_invalid_ghost_symbol_defaults_to_bq(tmp_path, monkeypatch):
    f = tmp_path / "settings.json"
    f.write_text(json.dumps({"ghost_symbol": "nonsense"}), encoding="utf-8")
    monkeypatch.setattr(pkg, "_SETTINGS_FILE", str(f))
    result = pkg._load_plugin_settings()
    assert result["ghost_symbol"] == "Bq"


def test_load_plugin_settings_valid_h_colon(tmp_path, monkeypatch):
    f = tmp_path / "settings.json"
    f.write_text(json.dumps({"ghost_symbol": "H:"}), encoding="utf-8")
    monkeypatch.setattr(pkg, "_SETTINGS_FILE", str(f))
    result = pkg._load_plugin_settings()
    assert result["ghost_symbol"] == "H:"


def test_save_plugin_settings_writes_file(tmp_path, monkeypatch):
    f = tmp_path / "settings.json"
    monkeypatch.setattr(pkg, "_SETTINGS_FILE", str(f))
    pkg._save_plugin_settings({"ghost_symbol": "H:"})
    assert json.loads(f.read_text(encoding="utf-8")) == {"ghost_symbol": "H:"}


def test_save_plugin_settings_swallows_write_error(monkeypatch):
    monkeypatch.setattr(pkg, "_SETTINGS_FILE", "Z:\\this\\path\\does\\not\\exist\\settings.json")
    pkg._save_plugin_settings({"ghost_symbol": "Bq"})  # should not raise


# ---------------------------------------------------------------------------
# on_load: sync_symbol_from_settings callback + exception branch
# ---------------------------------------------------------------------------

def test_on_load_calls_sync_symbol_from_settings_on_window():
    ctx = _StubContext()
    initialize(ctx)
    win = MagicMock()
    ctx._windows["main_panel"] = win
    ctx.load({"ghost_symbol": "H:"})
    win.sync_symbol_from_settings.assert_called_once()


def test_on_load_sync_symbol_exception_is_swallowed():
    ctx = _StubContext()
    initialize(ctx)
    win = MagicMock()
    win.sync_symbol_from_settings.side_effect = RuntimeError("boom")
    ctx._windows["main_panel"] = win
    ctx.load({"ghost_symbol": "H:"})  # should not raise


def test_on_load_window_without_sync_method_is_skipped():
    ctx = _StubContext()
    initialize(ctx)

    class _NoSyncWindow:
        pass

    ctx._windows["main_panel"] = _NoSyncWindow()
    ctx.load({"ghost_symbol": "H:"})  # should not raise (hasattr check fails)


def test_on_load_applies_labels_immediately_when_headless():
    """No Qt event loop in tests -> _apply_labels runs synchronously."""
    ctx = _StubContext()
    initialize(ctx)
    mol = MagicMock()
    atom = MagicMock()
    mol.GetAtomWithIdx.return_value = atom
    ctx.current_molecule = mol
    ctx.load({"bq_labels": {"0": "Bq"}})
    atom.SetProp.assert_called_once_with("custom_symbol", "Bq")


def test_on_load_apply_labels_warns_when_molecule_still_none():
    ctx = _StubContext()
    initialize(ctx)
    ctx.current_molecule = None
    ctx.load({"bq_labels": {"0": "Bq"}})  # should not raise, just logs a warning


def test_on_load_defers_via_qtimer_when_event_loop_present(monkeypatch):
    """When a real Qt event loop is detected, on_load must defer the label
    restore via QTimer.singleShot instead of applying immediately."""
    ctx = _StubContext()
    initialize(ctx)

    fake_qtimer = MagicMock()
    fake_app = MagicMock()
    fake_app.instance.return_value = MagicMock()  # truthy -> "running" app

    monkeypatch.setattr(pkg, "_QTimer", fake_qtimer)
    monkeypatch.setattr(pkg, "_QCoreApplication", fake_app)

    mol = MagicMock()
    ctx.current_molecule = mol
    ctx.load({"bq_labels": {"0": "Bq"}})

    fake_qtimer.singleShot.assert_called_once()
    args = fake_qtimer.singleShot.call_args[0]
    assert args[0] == 0
    # Immediate apply must NOT have happened yet (only scheduled)
    mol.GetAtomWithIdx.assert_not_called()
    # Invoking the deferred callback now performs the actual restore
    args[1]()
    mol.GetAtomWithIdx.assert_called_once_with(0)
