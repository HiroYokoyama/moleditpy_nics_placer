"""Real Qt lifecycle checks; run separately from the mocked tests/ suite."""

import os
from unittest.mock import MagicMock

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("PYTEST_QT_API", "pyqt6")

import pytest
from PyQt6.QtCore import Qt
from PyQt6.QtTest import QTest
from PyQt6.QtWidgets import QApplication
from rdkit import Chem
from rdkit.Chem import AllChem

from nics_placer import dialog, grid_dialog


@pytest.fixture(scope="session")
def qapp():
    app = QApplication.instance() or QApplication([])
    app.setQuitOnLastWindowClosed(False)
    yield app


@pytest.fixture
def context(monkeypatch):
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    ctx = MagicMock()
    ctx.current_molecule = mol
    ctx.plotter = None
    ctx.get_main_window.return_value = None
    monkeypatch.setattr(dialog, "_save_plugin_settings", lambda *_: None)
    monkeypatch.setattr(grid_dialog, "_save_plugin_settings", lambda *_: None)
    return ctx


@pytest.mark.parametrize("cls", [dialog.NicsPlacerDialog, grid_dialog.NicsGridDialog])
@pytest.mark.parametrize("action", ["escape", "accept", "close"])
def test_native_close_and_reopen_stops_polling(qapp, context, cls, action):
    window = cls(context)
    window.show()
    qapp.processEvents()
    assert window._poll_timer.isActive()
    if action == "escape":
        QTest.keyClick(window, Qt.Key.Key_Escape)
    else:
        getattr(window, action)()
    assert not window.isVisible()
    assert not window._poll_timer.isActive()
    context.plotter = MagicMock()
    window._render_spheres()
    context.plotter.add_mesh.assert_not_called()
    context.plotter = None
    window.show()
    qapp.processEvents()
    assert window._poll_timer.isActive()
    window.close()


def test_late_preview_does_not_redraw_after_escape(qapp, context):
    window = dialog.NicsPlacerDialog(context)
    window.show()
    window._stage_selected_nics0()
    window._apply_staged()
    QTest.keyClick(window, Qt.Key.Key_Escape)
    context.plotter = MagicMock()
    QTest.qWait(200)
    context.plotter.add_mesh.assert_not_called()
    window.close()
