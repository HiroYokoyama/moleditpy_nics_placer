"""
MoleditPy NICS Placer Plugin
=============================
Detects rings in the loaded molecule and places Bq ghost atoms at
NICS(0) and NICS(1) probe positions.

Uses the same ``custom_symbol`` atom property as the XYZ Editor, so ORCA
Input Generator Pro automatically renders the Bq labels in the coordinate
block without any additional configuration.
"""
import json
import logging
import os

PLUGIN_NAME = "NICS Placer"
PLUGIN_VERSION = "1.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = (
    "Detect rings and place Bq ghost atoms at NICS(0)/NICS(1) probe positions. "
    "Compatible with ORCA Input Generator Pro via the custom_symbol property."
)
PLUGIN_CATEGORY = "3D Edit"

_SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "settings.json")
_GHOST_SYMBOLS = {"Bq", "H:"}
_dialog_opened: bool = False   # guard: don't write to project file until plugin is used


def _default_settings() -> dict:
    return {"ghost_symbol": "Bq"}


def _load_plugin_settings() -> dict:
    """Load settings.json; fall back to defaults on any error."""
    try:
        with open(_SETTINGS_FILE, encoding="utf-8") as f:
            data = json.load(f)
        if not isinstance(data, dict):
            return _default_settings()
        # Validate known keys
        if data.get("ghost_symbol") not in _GHOST_SYMBOLS:
            data["ghost_symbol"] = "Bq"
        return data
    except FileNotFoundError:
        return _default_settings()
    except Exception as _e:
        logging.warning("[nics_placer/__init__.py] load settings.json: %s", _e)
        return _default_settings()


def _save_plugin_settings(settings: dict) -> None:
    """Write settings.json; log on failure."""
    try:
        with open(_SETTINGS_FILE, "w", encoding="utf-8") as f:
            json.dump(settings, f, indent=2)
    except Exception as _e:
        logging.warning("[nics_placer/__init__.py] save settings.json: %s", _e)


# Loaded once at import time — persists across documents within a session
_plugin_settings: dict = _load_plugin_settings()


def initialize(context):
    """MoleditPy V3 plugin entry point."""

    # ---------------------------------------------------------------
    # Menu action
    # ---------------------------------------------------------------
    def show_dialog():
        win = context.get_window("main_panel")
        if win:
            win.show()
            win.raise_()
            win.activateWindow()
            return

        mol = context.current_molecule
        if not mol:
            context.show_status_message("No molecule loaded.", 3000)
            return
        if not mol.GetNumConformers():
            context.show_status_message(
                "Molecule has no 3D coordinates — run 3D conversion first.", 4000
            )
            return

        from .dialog import NicsPlacerDialog
        mw = context.get_main_window()
        dlg = NicsPlacerDialog(context, parent=mw)
        context.register_window("main_panel", dlg)
        global _dialog_opened
        _dialog_opened = True
        dlg.show()

    context.add_menu_action("3D Edit/NICS Placer...", show_dialog)

    # ---------------------------------------------------------------
    # Project file persistence
    # ghost_symbol: per-project override (on_load overwrites plugin setting
    #               for this session; on_reset restores from settings.json)
    # bq_labels:    molecule-specific ghost atom indices
    # ---------------------------------------------------------------
    def on_save():
        if not _dialog_opened:
            return None
        result = {"ghost_symbol": _plugin_settings["ghost_symbol"]}
        mol = context.current_molecule
        if mol:
            labels = {
                str(a.GetIdx()): a.GetProp("custom_symbol")
                for a in mol.GetAtoms()
                if a.HasProp("custom_symbol")
                and a.GetProp("custom_symbol") in _GHOST_SYMBOLS
            }
            if labels:
                result["bq_labels"] = labels
        return result

    def on_load(data):
        if not isinstance(data, dict):
            return
        global _dialog_opened
        _dialog_opened = True
        # Per-project ghost_symbol overrides the plugin default for this session
        sym = data.get("ghost_symbol")
        if sym in _GHOST_SYMBOLS:
            _plugin_settings["ghost_symbol"] = sym
            win = context.get_window("main_panel")
            if win and hasattr(win, "sync_symbol_from_settings"):
                try:
                    win.sync_symbol_from_settings()
                except Exception as _e:
                    logging.warning("[nics_placer/__init__.py] sync combo: %s", _e)
        # Restore ghost atom labels onto molecule.
        # IMPORTANT: the app fires load handlers BEFORE restoring the 3D molecule,
        # so context.current_molecule is None at this point.  We defer the actual
        # label restore to the next Qt event-loop tick (after the project loader
        # has set view_3d_manager.current_mol) using QTimer.singleShot.
        labels = data.get("bq_labels", {})
        if not labels:
            return

        def _apply_labels():
            mw = context.get_main_window()
            mol = (
                mw.view_3d_manager.current_mol
                if mw and hasattr(mw, "view_3d_manager")
                else None
            )
            if not mol:
                logging.warning(
                    "[nics_placer/__init__.py] on_load: molecule still None after deferred restore"
                )
                return
            for idx_str, lbl in labels.items():
                try:
                    mol.GetAtomWithIdx(int(idx_str)).SetProp("custom_symbol", lbl)
                except Exception as _e:
                    logging.warning(
                        "[nics_placer/__init__.py] on_load atom %s: %s", idx_str, _e
                    )
            # Labels are applied in-place; no redraw needed here — the app's own
            # draw_molecule_3d call that follows the load will pick them up.

        try:
            from PyQt6.QtCore import QTimer
            QTimer.singleShot(0, _apply_labels)
        except ImportError:
            # Fallback for test environments without PyQt6
            _apply_labels()

    def on_reset():
        global _dialog_opened
        _dialog_opened = False
        # Restore plugin-level default from settings.json
        _plugin_settings.clear()
        _plugin_settings.update(_load_plugin_settings())
        win = context.get_window("main_panel")
        if win:
            try:
                win.close()
            except Exception as _e:
                logging.warning("[nics_placer/__init__.py] on_reset close: %s", _e)

    context.register_save_handler(on_save)
    context.register_load_handler(on_load)
    context.register_document_reset_handler(on_reset)
