"""
Regression tests for the 2.4.0 fixes: every test here fails on 2.3.2.

Each one pins a way the plugin used to put ghost atoms in the wrong place, or
put the wrong ones there, without any error.
"""

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
    from rdkit.Geometry import Point3D

    import nics_placer as pkg
    import nics_placer.dialog as dialog_mod
    import nics_placer.grid_dialog as grid_mod
    from nics_placer import ghosts
    from nics_placer.nics_math import get_rings

    _AVAILABLE = True
except Exception:
    pass

pytestmark = pytest.mark.skipif(not _AVAILABLE, reason="RDKit not importable")


def _mol_3d(smiles="c1ccccc1", seed=42):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    return mol


class _Ctx:
    def __init__(self, mol):
        self.current_molecule = mol
        self.plotter = MagicMock(name="plotter")
        self.plotter.add_mesh.side_effect = lambda mesh, name=None, **kw: MagicMock()
        self._windows = {}
        self.push_undo_checkpoint = MagicMock()
        self.show_status_message = MagicMock()

    def get_window(self, key):
        return self._windows.get(key)

    def register_window(self, key, win):
        self._windows[key] = win

    def get_main_window(self):
        return None


def _placer(ctx):
    dlg = dialog_mod.NicsPlacerDialog(ctx)
    ctx.register_window("main_panel", dlg)
    dlg.showEvent(MagicMock())
    return dlg


def _grid(ctx):
    dlg = grid_mod.NicsGridDialog(ctx)
    ctx.register_window("grid_panel", dlg)
    dlg.showEvent(MagicMock())
    return dlg


def _n_ghosts(mol):
    return sum(1 for a in mol.GetAtoms() if ghosts.is_ghost(a))


@pytest.fixture(autouse=True)
def _keep_settings(monkeypatch):
    if not _AVAILABLE:
        yield
        return
    saved = dict(pkg._plugin_settings)
    # Never write the real settings.json from a test.
    monkeypatch.setattr(dialog_mod, "_save_plugin_settings", lambda s: None)
    monkeypatch.setattr(grid_mod, "_save_plugin_settings", lambda s: None)
    pkg._plugin_settings["ghost_symbol"] = "Bq"
    try:
        yield
    finally:
        pkg._plugin_settings.clear()
        pkg._plugin_settings.update(saved)


# ---------------------------------------------------------------------------
# Probes follow the geometry, not the object
# ---------------------------------------------------------------------------


class TestProbesFollowInPlaceEdits:
    def test_moving_atoms_in_place_moves_the_probes(self):
        """The host drags, aligns and optimises inside the same Mol object."""
        mol = _mol_3d()
        ctx = _Ctx(mol)
        dlg = _placer(ctx)
        before = next(p["pos"] for p in dlg._nics_points if p["type"] == "nics0")

        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, Point3D(p.x + 5.0, p.y, p.z))
        dlg._check_molecule_changed()

        after = next(p["pos"] for p in dlg._nics_points if p["type"] == "nics0")
        assert after == pytest.approx(before + np.array([5.0, 0.0, 0.0]), abs=1e-6)

    def test_watcher_holds_the_molecule_not_its_id(self):
        watcher = ghosts.MoleculeWatcher(_mol_3d())
        assert watcher.changed(_mol_3d()) is True


class TestGridFollowsMoleculeChanges:
    def test_new_molecule_reloads_the_ring_table(self):
        ctx = _Ctx(_mol_3d("c1ccccc1"))
        dlg = _grid(ctx)
        assert len(dlg._rings) == 1
        ctx.current_molecule = _mol_3d("c1ccc2ccccc2c1")
        dlg._check_molecule_changed()
        assert len(dlg._rings) == 2
        assert dlg._table.rowCount() == 2

    def test_grid_is_anchored_to_the_new_molecule(self):
        ctx = _Ctx(_mol_3d("c1ccccc1"))
        dlg = _grid(ctx)
        dlg._use_com.setChecked(False)
        mol2 = _mol_3d("c1ccccc1", seed=7)
        conf = mol2.GetConformer()
        for i in range(mol2.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, Point3D(p.x + 10.0, p.y, p.z))
        ctx.current_molecule = mol2
        dlg._check_molecule_changed()
        centre = np.mean([p["pos"] for p in dlg._grid_points], axis=0)
        ring = get_rings(mol2)[0]["atoms"]
        want = np.mean([list(conf.GetAtomPosition(i)) for i in ring], axis=0)
        assert centre == pytest.approx(want, abs=1e-6)

    def test_close_stops_the_poll(self):
        dlg = _grid(_Ctx(_mol_3d()))
        assert dlg._poll_timer.isActive()
        dlg.closeEvent(MagicMock())
        assert not dlg._poll_timer.isActive()


# ---------------------------------------------------------------------------
# No stacked ghosts
# ---------------------------------------------------------------------------


class TestNoDuplicateGhosts:
    def test_place_grid_twice_adds_nothing_the_second_time(self):
        ctx = _Ctx(_mol_3d())
        dlg = _grid(ctx)
        dlg._place_grid()
        first = _n_ghosts(ctx.current_molecule)
        assert first == len(dlg._grid_points)
        dlg._rebuild_grid()
        dlg._place_grid()
        assert _n_ghosts(ctx.current_molecule) == first

    def test_place_all_skips_a_nics0_already_holding_a_grid_probe(self):
        mol = _mol_3d()
        ctx = _Ctx(mol)
        placer = _placer(ctx)
        nics0 = next(p["pos"] for p in placer._nics_points if p["type"] == "nics0")
        ctx.current_molecule = dialog_mod._add_bq_atom(mol, nics0)
        placer._check_molecule_changed()
        placer._place_all()
        assert _n_ghosts(ctx.current_molecule) == 3  # not 4

    def test_filter_is_linear_and_exact(self):
        mol = dialog_mod._add_bq_atom(_mol_3d(), np.array([0.0, 0.0, 1.0]))
        pts = [np.array([0.0, 0.0, 1.0 + 1e-9]), np.array([0.0, 0.0, 2.0])] * 2
        kept, skipped = ghosts.new_probe_positions(mol, pts)
        assert len(kept) == 1 and skipped == 3


# ---------------------------------------------------------------------------
# One ghost label per molecule
# ---------------------------------------------------------------------------


class TestOneGhostLabel:
    def test_grid_label_change_relabels_and_reaches_the_placer(self):
        ctx = _Ctx(_mol_3d())
        placer = _placer(ctx)
        grid = _grid(ctx)
        placer._place_all()

        grid._sym_combo.setCurrentIndex(grid._sym_combo.findData("H:"))
        grid._on_symbol_changed(0)

        labels = {
            a.GetProp("custom_symbol")
            for a in ctx.current_molecule.GetAtoms()
            if a.GetAtomicNum() == 0
        }
        assert labels == {"H:"}
        assert placer._ghost_symbol == "H:"

    def test_placer_label_change_reaches_the_grid(self):
        ctx = _Ctx(_mol_3d())
        placer = _placer(ctx)
        grid = _grid(ctx)
        placer._sym_combo.setCurrentIndex(placer._sym_combo.findData("H:"))
        placer._on_symbol_changed(0)
        assert grid._ghost_symbol == "H:"

    def test_project_load_syncs_the_grid_window(self):
        class _InitCtx(_Ctx):
            def __init__(self):
                super().__init__(None)
                self.add_menu_action = MagicMock()

            def register_save_handler(self, fn):
                pass

            def register_load_handler(self, fn):
                self.load = fn

            def register_document_reset_handler(self, fn):
                pass

        ctx = _InitCtx()
        pkg.initialize(ctx)
        grid = MagicMock()
        ctx._windows["grid_panel"] = grid
        ctx.load({"ghost_symbol": "H:"})
        grid.sync_symbol_from_settings.assert_called_once()


# ---------------------------------------------------------------------------
# The molecule is not rewritten
# ---------------------------------------------------------------------------


class TestMoleculeLeftAlone:
    def _kekule(self):
        mol = _mol_3d()
        Chem.Kekulize(mol, clearAromaticFlags=True)
        return mol

    def test_adding_a_probe_keeps_kekule_bonds(self):
        out = dialog_mod._add_bq_atom(self._kekule(), np.zeros(3))
        assert {str(b.GetBondType()) for b in out.GetBonds()} == {"SINGLE", "DOUBLE"}

    def test_clearing_probes_keeps_kekule_bonds(self):
        mol = dialog_mod._add_bq_atom(self._kekule(), np.zeros(3))
        out = dialog_mod._remove_all_bq(mol)
        assert {str(b.GetBondType()) for b in out.GetBonds()} == {"SINGLE", "DOUBLE"}

    def test_ring_table_is_ready_after_placement(self):
        out = dialog_mod._add_bq_atom(_mol_3d(), np.zeros(3))
        assert out.GetRingInfo().NumRings() == 1


class TestRingsOnUnsanitisableMolecule:
    def test_rings_found_when_the_ring_table_is_empty(self):
        # Built atom by atom, the way an editor does, so no sanitisation ever
        # fills in the ring table.
        rw = Chem.RWMol()
        for _ in range(6):
            rw.AddAtom(Chem.Atom(6))
        for i in range(6):
            rw.AddBond(i, (i + 1) % 6, Chem.BondType.SINGLE)
        rw.AddBond(0, 3, Chem.BondType.SINGLE)
        rw.AddBond(0, 2, Chem.BondType.DOUBLE)  # C0 now over-valent
        rw.UpdatePropertyCache(strict=False)
        conf = Chem.Conformer(rw.GetNumAtoms())
        for i in range(rw.GetNumAtoms()):
            a = 2 * np.pi * i / 6
            conf.SetAtomPosition(i, Point3D(1.4 * np.cos(a), 1.4 * np.sin(a), 0.0))
        rw.AddConformer(conf)
        mol = rw.GetMol()
        assert mol.GetRingInfo().AtomRings() == ()
        assert len(get_rings(mol)) == 3


# ---------------------------------------------------------------------------
# Visible leftovers and honest spacing
# ---------------------------------------------------------------------------


class TestOtherGhostsAreReported:
    def test_probes_from_an_old_height_are_reported(self):
        ctx = _Ctx(_mol_3d())
        dlg = _placer(ctx)
        dlg._place_all()
        dlg._check_molecule_changed()
        dlg._other_ghosts_label.setVisible.assert_called_with(False)

        dlg._height_spin.setValue(2.0)
        dlg._on_height_changed(2.0)
        dlg._other_ghosts_label.setVisible.assert_called_with(True)
        assert "2 other ghost atoms" in dlg._other_ghosts_label.setText.call_args[0][0]


class TestUniformSpacingCap:
    def test_capped_count_is_flagged(self):
        dlg = _grid(_Ctx(_mol_3d()))
        dlg._uniform_spacing.setChecked(True)
        dlg._on_uniform_toggled(True)
        dlg._spacing_spin.setValue(0.05)
        dlg._on_params_changed()
        text = dlg._count_label.setText.call_args[0][0]
        assert "not uniform" in text

    def test_uncapped_grid_is_not_flagged(self):
        dlg = _grid(_Ctx(_mol_3d()))
        dlg._uniform_spacing.setChecked(True)
        dlg._on_uniform_toggled(True)
        text = dlg._count_label.setText.call_args[0][0]
        assert "not uniform" not in text
