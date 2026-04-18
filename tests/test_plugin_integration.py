"""
Integration tests for nics_placer/__init__.py — verifies plugin contract
(save / load / reset handlers) without Qt, RDKit, or PyVista.
"""
import sys
import os
import unittest
from unittest.mock import MagicMock, call

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), "..")))

import nics_placer as _pkg
from nics_placer import initialize, PLUGIN_NAME, PLUGIN_VERSION


# ---------------------------------------------------------------------------
# Stub PluginContext
# ---------------------------------------------------------------------------

class _StubContext:
    def __init__(self):
        self._windows = {}
        self._save_handler = None
        self._load_handler = None
        self._reset_handler = None
        self.current_molecule = None
        self.plotter = None
        self.add_analysis_tool = MagicMock()
        self.add_menu_action = MagicMock()
        self.show_status_message = MagicMock()

    def get_main_window(self):
        return MagicMock()

    def get_window(self, key):
        return self._windows.get(key)

    def register_window(self, key, win):
        self._windows[key] = win

    def register_save_handler(self, fn):
        self._save_handler = fn

    def register_load_handler(self, fn):
        self._load_handler = fn

    def register_document_reset_handler(self, fn):
        self._reset_handler = fn

    # convenience
    def save(self):
        return self._save_handler() if self._save_handler else None

    def load(self, data):
        if self._load_handler:
            self._load_handler(data)

    def reset(self):
        if self._reset_handler:
            self._reset_handler()


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _make_mol_with_bq(bq_indices):
    """Return a mock mol that has Bq custom_symbol on given atom indices."""
    mol = MagicMock()
    atoms = []
    for i in range(max(bq_indices) + 1 if bq_indices else 1):
        a = MagicMock()
        a.GetIdx.return_value = i
        if i in bq_indices:
            a.HasProp.return_value = True
            a.GetProp.return_value = "Bq"
        else:
            a.HasProp.return_value = False
            a.GetProp.return_value = ""
        atoms.append(a)
    mol.GetAtoms.return_value = atoms
    return mol


# ---------------------------------------------------------------------------
# Tests: metadata
# ---------------------------------------------------------------------------

class TestMetadata(unittest.TestCase):
    def test_plugin_name(self):
        self.assertEqual(PLUGIN_NAME, "NICS Placer")

    def test_plugin_version(self):
        self.assertIsInstance(PLUGIN_VERSION, str)
        parts = PLUGIN_VERSION.split(".")
        self.assertEqual(len(parts), 3)

    def test_initialize_registers_menu_action(self):
        ctx = _StubContext()
        initialize(ctx)
        ctx.add_menu_action.assert_called_once()
        path_arg = ctx.add_menu_action.call_args[0][0]
        self.assertIn("3D Edit", path_arg)
        self.assertIn("NICS", path_arg)


# ---------------------------------------------------------------------------
# Tests: save handler
# ---------------------------------------------------------------------------

class TestSaveHandler(unittest.TestCase):
    def setUp(self):
        self.ctx = _StubContext()
        initialize(self.ctx)
        import nics_placer as pkg
        pkg._dialog_opened = False  # start clean

    def test_save_returns_none_before_dialog_opened(self):
        self.ctx.current_molecule = None
        self.assertIsNone(self.ctx.save())

    def test_save_includes_ghost_symbol_after_dialog_opened(self):
        import nics_placer as pkg
        pkg._dialog_opened = True
        self.ctx.current_molecule = None
        result = self.ctx.save()
        self.assertIsNotNone(result)
        self.assertIn("ghost_symbol", result)

    def test_save_no_bq_labels_key_when_no_ghost_atoms(self):
        import nics_placer as pkg
        pkg._dialog_opened = True
        mol = MagicMock()
        atom = MagicMock()
        atom.HasProp.return_value = False
        mol.GetAtoms.return_value = [atom]
        self.ctx.current_molecule = mol
        result = self.ctx.save()
        self.assertNotIn("bq_labels", result)

    def test_save_returns_bq_labels(self):
        import nics_placer as pkg
        pkg._dialog_opened = True
        mol = _make_mol_with_bq({2, 5})
        self.ctx.current_molecule = mol
        result = self.ctx.save()
        self.assertIsNotNone(result)
        self.assertIn("bq_labels", result)
        labels = result["bq_labels"]
        self.assertIn("2", labels)
        self.assertIn("5", labels)
        self.assertEqual(labels["2"], "Bq")
        self.assertEqual(labels["5"], "Bq")

    def test_save_excludes_non_bq_atoms(self):
        import nics_placer as pkg
        pkg._dialog_opened = True
        mol = MagicMock()
        atoms = []
        for i in range(3):
            a = MagicMock()
            a.GetIdx.return_value = i
            if i == 1:
                a.HasProp.return_value = True
                a.GetProp.return_value = "Bq"
            else:
                a.HasProp.return_value = False
                a.GetProp.return_value = ""
            atoms.append(a)
        mol.GetAtoms.return_value = atoms
        self.ctx.current_molecule = mol
        result = self.ctx.save()
        labels = result["bq_labels"]
        self.assertEqual(list(labels.keys()), ["1"])


# ---------------------------------------------------------------------------
# Tests: load handler
# ---------------------------------------------------------------------------

class TestLoadHandler(unittest.TestCase):
    def setUp(self):
        self.ctx = _StubContext()
        initialize(self.ctx)

    def test_load_none_is_safe(self):
        self.ctx.load(None)  # should not raise

    def test_load_empty_dict_is_safe(self):
        self.ctx.load({})  # should not raise

    def test_load_sets_dialog_opened_flag(self):
        import nics_placer as pkg
        pkg._dialog_opened = False
        self.ctx.load({"ghost_symbol": "Bq"})
        self.assertTrue(pkg._dialog_opened)

    def test_load_restores_ghost_symbol(self):
        import nics_placer as pkg
        pkg._plugin_settings["ghost_symbol"] = "Bq"
        self.ctx.load({"ghost_symbol": "H:"})
        self.assertEqual(pkg._plugin_settings["ghost_symbol"], "H:")

    def test_load_ignores_invalid_ghost_symbol(self):
        import nics_placer as pkg
        pkg._plugin_settings["ghost_symbol"] = "Bq"
        self.ctx.load({"ghost_symbol": "INVALID"})
        self.assertEqual(pkg._plugin_settings["ghost_symbol"], "Bq")

    def test_load_stamps_bq_labels(self):
        mol = MagicMock()
        atom = MagicMock()
        mol.GetAtomWithIdx.return_value = atom
        self.ctx.current_molecule = mol
        self.ctx.load({"bq_labels": {"3": "Bq", "7": "Bq"}})
        calls = atom.SetProp.call_args_list
        self.assertEqual(len(calls), 2)
        for c in calls:
            self.assertEqual(c[0][0], "custom_symbol")
            self.assertEqual(c[0][1], "Bq")

    def test_load_ignores_invalid_index(self):
        mol = MagicMock()
        mol.GetAtomWithIdx.side_effect = Exception("bad index")
        self.ctx.current_molecule = mol
        # Should not propagate exception
        self.ctx.load({"bq_labels": {"999": "Bq"}})


# ---------------------------------------------------------------------------
# Tests: reset handler
# ---------------------------------------------------------------------------

class TestResetHandler(unittest.TestCase):
    def setUp(self):
        self.ctx = _StubContext()
        initialize(self.ctx)

    def test_reset_restores_default_ghost_symbol(self):
        import nics_placer as pkg
        pkg._plugin_settings["ghost_symbol"] = "H:"
        self.ctx.reset()
        self.assertEqual(pkg._plugin_settings["ghost_symbol"], "Bq")

    def test_reset_clears_dialog_opened_flag(self):
        import nics_placer as pkg
        pkg._dialog_opened = True
        self.ctx.reset()
        self.assertFalse(pkg._dialog_opened)
        self.assertIsNone(self.ctx.save())

    def test_reset_when_no_window_is_safe(self):
        self.ctx.reset()  # should not raise

    def test_reset_closes_open_window(self):
        win = MagicMock()
        self.ctx._windows["main_panel"] = win
        self.ctx.reset()
        win.close.assert_called_once()

    def test_reset_handles_close_exception(self):
        win = MagicMock()
        win.close.side_effect = Exception("already closed")
        self.ctx._windows["main_panel"] = win
        self.ctx.reset()  # should not raise


# ---------------------------------------------------------------------------
# Tests: roundtrip save → load
# ---------------------------------------------------------------------------

class TestRoundtrip(unittest.TestCase):
    def test_roundtrip_preserves_bq_labels(self):
        import nics_placer as pkg
        ctx1 = _StubContext()
        initialize(ctx1)
        pkg._dialog_opened = True
        mol = _make_mol_with_bq({1, 4})
        ctx1.current_molecule = mol
        saved = ctx1.save()
        self.assertIsNotNone(saved)

        ctx2 = _StubContext()
        initialize(ctx2)
        target_mol = MagicMock()
        target_atom = MagicMock()
        target_mol.GetAtomWithIdx.return_value = target_atom
        ctx2.current_molecule = target_mol
        ctx2.load(saved)
        # SetProp should have been called for each bq label
        self.assertEqual(target_atom.SetProp.call_count, 2)


if __name__ == "__main__":
    unittest.main()
