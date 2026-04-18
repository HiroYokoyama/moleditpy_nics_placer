"""
MoleditPy NICS Placer Plugin
=============================
Detects rings in the loaded molecule and places Bq ghost atoms at
NICS(0) and NICS(1) probe positions.

Uses the same ``custom_symbol`` atom property as the XYZ Editor, so ORCA
Input Generator Pro automatically renders the Bq labels in the coordinate
block without any additional configuration.
"""
import logging

PLUGIN_NAME = "NICS Placer"
PLUGIN_VERSION = "1.0.0"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = (
    "Detect rings and place Bq ghost atoms at NICS(0)/NICS(1) probe positions. "
    "Compatible with ORCA Input Generator Pro via the custom_symbol property."
)
PLUGIN_CATEGORY = "Analysis"


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
        dlg.show()

    context.add_analysis_tool("NICS Placer...", show_dialog)

    # ---------------------------------------------------------------
    # Project file persistence
    # Save only the Bq atom positions (index → custom_symbol).
    # On load, re-stamp custom_symbol onto the matching atom indices.
    # ---------------------------------------------------------------
    def on_save():
        mol = context.current_molecule
        if not mol:
            return None
        labels = {
            str(a.GetIdx()): a.GetProp("custom_symbol")
            for a in mol.GetAtoms()
            if a.HasProp("custom_symbol") and a.GetProp("custom_symbol") == "Bq"
        }
        return {"bq_labels": labels} if labels else None

    def on_load(data):
        if not isinstance(data, dict):
            return
        labels = data.get("bq_labels", {})
        if not labels:
            return
        mol = context.current_molecule
        if not mol:
            return
        for idx_str, lbl in labels.items():
            try:
                mol.GetAtomWithIdx(int(idx_str)).SetProp("custom_symbol", lbl)
            except Exception as _e:
                logging.warning("[nics_placer/__init__.py] on_load atom %s: %s", idx_str, _e)
        context.current_molecule = mol

    def on_reset():
        win = context.get_window("main_panel")
        if win:
            try:
                win.close()
            except Exception as _e:
                logging.warning("[nics_placer/__init__.py] on_reset close: %s", _e)

    context.register_save_handler(on_save)
    context.register_load_handler(on_load)
    context.register_document_reset_handler(on_reset)
