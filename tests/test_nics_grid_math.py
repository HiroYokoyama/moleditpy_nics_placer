"""
Geometry of the NICS probe grids, and the two audit fixes that preceded them.

Most of this is pure numpy and runs everywhere; the handful of cases that need
real molecules are guarded, since CI installs numpy and pytest only.
"""

import math
import sys
import unittest

import numpy as np

sys.path.insert(0, __file__.rsplit("tests", 1)[0])

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:  # pragma: no cover - CI runs without rdkit
    Chem = None
    AllChem = None

from nics_placer.nics_math import (  # noqa: E402
    GRID_PLANES,
    GRID_PLANE_LABELS,
    LAB_GRID_PLANES,
    PLANARITY_TOLERANCE,
    RING_GRID_PLANES,
    _reperceived_aromaticity,
    axis_extents,
    compute_nics_grid,
    compute_nics_points,
    compute_nics_volume,
    counts_for_spacing,
    grid_axes,
    heavy_atom_positions,
    molecular_reference_normal,
    molecule_bounds,
    orient_normal,
    planarity_rms,
    ring_centroid,
    ring_frame,
    ring_normal,
)


def _hexagon(z=0.0, radius=1.4):
    """Planar six-ring in the z = *z* plane."""
    return np.array(
        [
            [radius * math.cos(t), radius * math.sin(t), z]
            for t in np.linspace(0, 2 * math.pi, 6, endpoint=False)
        ]
    )


def _puckered_hexagon(amplitude=0.25, radius=1.4):
    """Chair-like alternating pucker — a cyclohexane stand-in."""
    pts = _hexagon(radius=radius)
    pts[:, 2] = [amplitude if i % 2 == 0 else -amplitude for i in range(6)]
    return pts


def _build(smiles, seed=7):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    return mol


# ---------------------------------------------------------------------------
# Fix 1: normal sign must be resolved, not left to LAPACK
# ---------------------------------------------------------------------------


class TestOrientNormal(unittest.TestCase):
    def test_keeps_a_normal_already_aligned_with_the_reference(self):
        n = np.array([0.0, 0.0, 1.0])
        out = orient_normal(n, reference=np.array([0.0, 0.0, 1.0]))
        np.testing.assert_allclose(out, n)

    def test_flips_a_normal_opposed_to_the_reference(self):
        n = np.array([0.0, 0.0, -1.0])
        out = orient_normal(n, reference=np.array([0.0, 0.0, 1.0]))
        np.testing.assert_allclose(out, np.array([0.0, 0.0, 1.0]))

    def test_perpendicular_reference_falls_back_to_the_deterministic_rule(self):
        """Two rings at right angles share no "up"; the answer must at least be
        reproducible rather than whatever the eigensolver returned."""
        n = np.array([0.0, 0.0, -1.0])
        ref = np.array([1.0, 0.0, 0.0])
        first = orient_normal(n, reference=ref)
        second = orient_normal(n.copy(), reference=ref)
        np.testing.assert_allclose(first, second)
        np.testing.assert_allclose(first, np.array([0.0, 0.0, 1.0]))

    def test_no_reference_makes_the_largest_component_positive(self):
        np.testing.assert_allclose(
            orient_normal(np.array([0.0, -0.8, 0.6])),
            np.array([0.0, 0.8, -0.6]),
        )

    def test_ring_normal_honours_the_reference(self):
        pts = _hexagon()
        up = ring_normal(pts, reference=np.array([0.0, 0.0, 1.0]))
        down = ring_normal(pts, reference=np.array([0.0, 0.0, -1.0]))
        np.testing.assert_allclose(up, -down, atol=1e-12)
        self.assertGreater(up[2], 0.0)

    def test_degenerate_ring_returns_a_usable_axis(self):
        """Collinear points have no plane; the caller still needs a unit vector
        rather than a NaN that silently poisons every probe position."""
        pts = np.zeros((4, 3))
        n = ring_normal(pts)
        self.assertAlmostEqual(float(np.linalg.norm(n)), 1.0)
        self.assertFalse(np.isnan(n).any())


@unittest.skipIf(Chem is None, "rdkit not installed")
class TestFaceConsistencyAcrossRings(unittest.TestCase):
    """The defect this fixes: naphthalene's two ring normals came back
    antiparallel, so 'nics1_above' meant opposite molecular faces for ring 0
    and ring 1. Any cross-ring comparison of NICS(1) on a molecule whose faces
    differ was then comparing two different things."""

    def _above_sides(self, smiles):
        from nics_placer.nics_math import get_ring_positions, get_rings

        mol = _build(smiles)
        ref = molecular_reference_normal(mol)
        conf = mol.GetConformer()
        heavy = np.array(
            [
                list(conf.GetAtomPosition(a.GetIdx()))
                for a in mol.GetAtoms()
                if a.GetAtomicNum() > 1
            ]
        )
        center = heavy.mean(axis=0)
        sides = []
        for ring in get_rings(mol):
            pts = get_ring_positions(mol, ring["atoms"])
            nics = compute_nics_points(pts, reference=ref)
            sides.append(float((nics["nics1_above"] - center) @ ref))
        return sides

    def test_naphthalene_rings_agree_on_which_face_is_above(self):
        sides = self._above_sides("c1ccc2ccccc2c1")
        self.assertEqual(len(sides), 2)
        self.assertTrue(all(s > 0 for s in sides), sides)

    def test_anthracene_rings_agree(self):
        sides = self._above_sides("c1ccc2cc3ccccc3cc2c1")
        self.assertEqual(len(sides), 3)
        self.assertTrue(all(s > 0 for s in sides), sides)

    def test_pyrene_rings_agree(self):
        sides = self._above_sides("c1cc2ccc3cccc4ccc(c1)c2c34")
        self.assertTrue(all(s > 0 for s in sides), sides)

    def test_above_and_below_stay_on_opposite_faces(self):
        from nics_placer.nics_math import get_ring_positions, get_rings

        mol = _build("c1ccccc1")
        pts = get_ring_positions(mol, get_rings(mol)[0]["atoms"])
        nics = compute_nics_points(pts, reference=molecular_reference_normal(mol))
        span = nics["nics1_above"] - nics["nics1_below"]
        self.assertAlmostEqual(float(np.linalg.norm(span)), 2.0, places=6)

    def test_reference_normal_ignores_a_molecule_without_rings(self):
        mol = _build("CCCC")
        np.testing.assert_allclose(
            molecular_reference_normal(mol), np.array([0.0, 0.0, 1.0])
        )


# ---------------------------------------------------------------------------
# Fix 2: aromaticity must be perceived, not read from stale flags
# ---------------------------------------------------------------------------


@unittest.skipIf(Chem is None, "rdkit not installed")
class TestAromaticityReperception(unittest.TestCase):
    """MMFF optimisation kekulises in place and RDKit does not restore the
    aromatic flags, so a molecule arriving straight from the force field
    reported azulene as non-aromatic."""

    def test_azulene_is_still_aromatic_after_mmff(self):
        from nics_placer.nics_math import get_rings

        mol = _build("c1ccc2cccc2cc1")
        self.assertFalse(
            any(a.GetIsAromatic() for a in mol.GetAtoms()),
            "precondition: MMFF should have cleared the cached flags",
        )
        self.assertTrue(all(r["is_aromatic"] for r in get_rings(mol)))

    def test_benzene_naphthalene_pyrrole_furan_are_aromatic(self):
        from nics_placer.nics_math import get_rings

        for smiles in ("c1ccccc1", "c1ccc2ccccc2c1", "c1cc[nH]c1", "c1ccoc1"):
            with self.subTest(smiles=smiles):
                self.assertTrue(
                    all(r["is_aromatic"] for r in get_rings(_build(smiles)))
                )

    def test_saturated_and_isolated_alkene_rings_are_not_aromatic(self):
        from nics_placer.nics_math import get_rings

        for smiles in ("C1CCCCC1", "C1CCC=CC1", "C1CCCC1"):
            with self.subTest(smiles=smiles):
                self.assertFalse(
                    any(r["is_aromatic"] for r in get_rings(_build(smiles)))
                )

    def test_perception_does_not_mutate_the_caller_molecule(self):
        mol = _build("c1ccc2cccc2cc1")
        before = [a.GetIsAromatic() for a in mol.GetAtoms()]
        _reperceived_aromaticity(mol)
        self.assertEqual([a.GetIsAromatic() for a in mol.GetAtoms()], before)

    def test_unperceivable_molecule_falls_back_to_the_original(self):
        """A partially built or odd-valence molecule is not a reason to lose
        the ring table entirely."""
        sentinel = object()
        self.assertIs(_reperceived_aromaticity(sentinel), sentinel)


# ---------------------------------------------------------------------------
# Planarity reporting
# ---------------------------------------------------------------------------


class TestPlanarity(unittest.TestCase):
    def test_flat_ring_has_zero_deviation(self):
        self.assertAlmostEqual(planarity_rms(_hexagon()), 0.0, places=12)

    def test_puckered_ring_reports_its_amplitude(self):
        # Every atom sits exactly +/-0.25 off the mean plane, so the RMS is 0.25.
        self.assertAlmostEqual(planarity_rms(_puckered_hexagon(0.25)), 0.25, places=6)

    def test_tolerance_separates_arenes_from_puckered_rings(self):
        self.assertLessEqual(planarity_rms(_hexagon()), PLANARITY_TOLERANCE)
        self.assertGreater(planarity_rms(_puckered_hexagon()), PLANARITY_TOLERANCE)

    def test_compute_nics_points_reports_normal_and_planarity(self):
        result = compute_nics_points(_puckered_hexagon(0.25))
        self.assertAlmostEqual(result["planarity"], 0.25, places=6)
        self.assertAlmostEqual(float(np.linalg.norm(result["normal"])), 1.0)


# ---------------------------------------------------------------------------
# Ring frame
# ---------------------------------------------------------------------------


class TestRingFrame(unittest.TestCase):
    def test_frame_is_orthonormal_and_right_handed(self):
        c, u, v, n = ring_frame(_hexagon())
        for vec in (u, v, n):
            self.assertAlmostEqual(float(np.linalg.norm(vec)), 1.0, places=12)
        self.assertAlmostEqual(float(u @ v), 0.0, places=12)
        self.assertAlmostEqual(float(u @ n), 0.0, places=12)
        self.assertAlmostEqual(float(v @ n), 0.0, places=12)
        self.assertAlmostEqual(
            float(np.linalg.det(np.array([u, v, n]))), 1.0, places=12
        )
        np.testing.assert_allclose(c, ring_centroid(_hexagon()), atol=1e-12)

    def test_u_points_at_the_first_ring_atom(self):
        """Anchoring u to an atom keeps the grid oriented to the ring rather
        than to the lab axes, so rotating the molecule rotates the grid with it."""
        pts = _hexagon()
        c, u, _v, _n = ring_frame(pts)
        direction = pts[0] - c
        direction /= np.linalg.norm(direction)
        np.testing.assert_allclose(u, direction, atol=1e-12)

    def test_frame_survives_a_first_atom_on_the_axis(self):
        """Geometrically impossible for a real ring, but must not emit NaNs.

        The remaining three points sum to zero, so the first — at the origin —
        is exactly the centroid and the usual "point u at atom 0" rule gives a
        zero-length vector.
        """
        pts = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-0.5, 0.87, 0.0], [-0.5, -0.87, 0.0]]
        )
        np.testing.assert_allclose(ring_centroid(pts), pts[0], atol=1e-12)
        _c, u, v, n = ring_frame(pts)
        for vec in (u, v, n):
            self.assertFalse(np.isnan(vec).any())
            self.assertAlmostEqual(float(np.linalg.norm(vec)), 1.0, places=10)
        self.assertAlmostEqual(float(u @ n), 0.0, places=10)

    def test_axis_fallback_avoids_a_seed_parallel_to_the_normal(self):
        """Same degeneracy, but with the ring normal along x — where the x seed
        is itself unusable and the y seed has to take over, or u would come out
        zero-length all over again."""
        pts = np.array(
            [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -0.5, 0.87], [0.0, -0.5, -0.87]]
        )
        np.testing.assert_allclose(ring_centroid(pts), pts[0], atol=1e-12)
        _c, u, _v, n = ring_frame(pts)
        self.assertAlmostEqual(abs(float(n[0])), 1.0, places=10)
        self.assertAlmostEqual(float(np.linalg.norm(u)), 1.0, places=10)
        self.assertAlmostEqual(float(u @ n), 0.0, places=10)

    def test_degenerate_svd_still_yields_a_unit_normal(self):
        """`np.linalg.svd` always returns unit rows, so the zero-length guard in
        `ring_normal` cannot be reached through real coordinates — it exists for
        a degenerate solver result. Force that result to prove the guard works
        rather than leaving it as untested defence."""
        import nics_placer.nics_math as math_mod

        real_svd = np.linalg.svd
        try:
            np.linalg.svd = lambda *a, **kw: (None, None, np.zeros((3, 3)))
            n = math_mod.ring_normal(_hexagon())
        finally:
            np.linalg.svd = real_svd
        np.testing.assert_allclose(n, np.array([0.0, 0.0, 1.0]))


# ---------------------------------------------------------------------------
# grid_axes
# ---------------------------------------------------------------------------


class TestGridAxes(unittest.TestCase):
    def test_every_plane_has_a_label(self):
        self.assertEqual(set(GRID_PLANES), set(GRID_PLANE_LABELS))
        self.assertEqual(GRID_PLANES, RING_GRID_PLANES + LAB_GRID_PLANES)

    def test_lab_planes_use_the_lab_axes(self):
        expected = {
            "xy": ([1, 0, 0], [0, 1, 0], [0, 0, 1]),
            "xz": ([1, 0, 0], [0, 0, 1], [0, -1, 0]),
            "yz": ([0, 1, 0], [0, 0, 1], [1, 0, 0]),
        }
        for plane, axes in expected.items():
            with self.subTest(plane=plane):
                got = grid_axes(plane)
                for g, e in zip(got, axes):
                    np.testing.assert_allclose(g, np.array(e, dtype=float), atol=1e-12)

    def test_ring_planes_permute_the_ring_frame(self):
        pts = _hexagon()
        _c, u, v, n = ring_frame(pts)
        for plane, expected in (
            ("parallel", (u, v, n)),
            ("perpendicular_u", (u, n, v)),
            ("perpendicular_v", (v, n, u)),
        ):
            with self.subTest(plane=plane):
                for got, exp in zip(grid_axes(plane, pts), expected):
                    np.testing.assert_allclose(got, exp, atol=1e-12)

    def test_all_axis_triples_are_orthonormal(self):
        for plane in GRID_PLANES:
            with self.subTest(plane=plane):
                a1, a2, a3 = grid_axes(plane, _hexagon())
                self.assertAlmostEqual(float(a1 @ a2), 0.0, places=12)
                self.assertAlmostEqual(float(a1 @ a3), 0.0, places=12)
                self.assertAlmostEqual(float(a2 @ a3), 0.0, places=12)

    def test_unknown_plane_is_rejected(self):
        with self.assertRaises(ValueError):
            grid_axes("diagonal")

    def test_ring_plane_without_positions_is_rejected(self):
        with self.assertRaises(ValueError):
            grid_axes("parallel", None)


# ---------------------------------------------------------------------------
# 2D grid
# ---------------------------------------------------------------------------


class TestComputeNicsGrid(unittest.TestCase):
    def test_point_count_is_n_squared(self):
        for n in (2, 5, 9, 21):
            with self.subTest(n=n):
                self.assertEqual(len(compute_nics_grid(_hexagon(), n_points=n)), n * n)

    def test_indices_cover_the_full_lattice(self):
        grid = compute_nics_grid(_hexagon(), n_points=4)
        self.assertEqual(
            {(p["i"], p["j"]) for p in grid},
            {(i, j) for i in range(4) for j in range(4)},
        )

    def test_spacing_is_uniform_and_spans_the_requested_extent(self):
        grid = compute_nics_grid(_hexagon(), n_points=5, extent=3.0)
        a = sorted({round(p["a"], 9) for p in grid})
        self.assertEqual(a, [-3.0, -1.5, 0.0, 1.5, 3.0])

    def test_parallel_grid_lies_in_a_plane_at_the_requested_offset(self):
        pts = _hexagon()
        _c, _u, _v, n = ring_frame(pts)
        grid = compute_nics_grid(pts, plane="parallel", n_points=5, offset=1.0)
        heights = [(p["pos"] - ring_centroid(pts)) @ n for p in grid]
        np.testing.assert_allclose(heights, np.ones(25), atol=1e-12)

    def test_parallel_grid_centre_is_exactly_the_nics1_probe(self):
        """The grid must agree with the single-probe dialog, or the two menu
        entries would disagree about where NICS(1) is."""
        pts = _hexagon()
        ref = np.array([0.0, 0.0, 1.0])
        grid = compute_nics_grid(
            pts, plane="parallel", n_points=5, offset=1.0, reference=ref
        )
        centre = next(p["pos"] for p in grid if p["i"] == 2 and p["j"] == 2)
        nics = compute_nics_points(pts, reference=ref)
        np.testing.assert_allclose(centre, nics["nics1_above"], atol=1e-12)

    def test_perpendicular_grid_contains_the_ring_normal(self):
        pts = _hexagon()
        _c, _u, _v, n = ring_frame(pts)
        grid = compute_nics_grid(pts, plane="perpendicular_u", n_points=5, extent=3.0)
        heights = [(p["pos"] - ring_centroid(pts)) @ n for p in grid]
        self.assertAlmostEqual(min(heights), -3.0, places=9)
        self.assertAlmostEqual(max(heights), 3.0, places=9)

    def test_lab_grid_defaults_to_the_origin_and_needs_no_ring(self):
        grid = compute_nics_grid(plane="xy", n_points=3, extent=1.0)
        np.testing.assert_allclose(
            next(p["pos"] for p in grid if p["i"] == 1 and p["j"] == 1),
            np.zeros(3),
            atol=1e-12,
        )

    def test_lab_grid_honours_an_explicit_centre(self):
        centre = np.array([1.0, 2.0, 3.0])
        grid = compute_nics_grid(plane="xy", n_points=3, extent=1.0, center=centre)
        np.testing.assert_allclose(
            next(p["pos"] for p in grid if p["i"] == 1 and p["j"] == 1),
            centre,
            atol=1e-12,
        )

    def test_offset_moves_a_lab_grid_along_its_own_normal(self):
        grid = compute_nics_grid(plane="xy", n_points=3, extent=1.0, offset=2.0)
        self.assertTrue(all(abs(p["pos"][2] - 2.0) < 1e-12 for p in grid))

    def test_a_rectangular_grid_uses_independent_counts(self):
        grid = compute_nics_grid(_hexagon(), n_points=(3, 7))
        self.assertEqual(len(grid), 21)
        self.assertEqual(max(p["i"] for p in grid), 2)
        self.assertEqual(max(p["j"] for p in grid), 6)

    def test_a_rectangular_grid_uses_independent_extents(self):
        grid = compute_nics_grid(_hexagon(), n_points=3, extent=(1.0, 5.0))
        self.assertAlmostEqual(max(p["a"] for p in grid), 1.0)
        self.assertAlmostEqual(max(p["b"] for p in grid), 5.0)

    def test_a_single_point_axis_sits_at_the_centre(self):
        """n=1 means "just this line"; landing it at -extent would shift the
        whole row off the feature it describes."""
        grid = compute_nics_grid(_hexagon(), n_points=(1, 5), extent=4.0)
        self.assertEqual(len(grid), 5)
        self.assertTrue(all(p["a"] == 0.0 for p in grid))

    def test_rejects_a_degenerate_point_count(self):
        with self.assertRaises(ValueError):
            compute_nics_grid(_hexagon(), n_points=0)

    def test_rejects_a_negative_extent(self):
        with self.assertRaises(ValueError):
            compute_nics_grid(_hexagon(), extent=-1.0)

    def test_rejects_a_wrong_length_per_axis_sequence(self):
        with self.assertRaises(ValueError):
            compute_nics_grid(_hexagon(), n_points=(3, 3, 3))
        with self.assertRaises(ValueError):
            compute_nics_grid(_hexagon(), extent=(1.0,))


# ---------------------------------------------------------------------------
# 3D volume
# ---------------------------------------------------------------------------


class TestComputeNicsVolume(unittest.TestCase):
    def test_point_count_is_n_squared_times_the_normal_layers(self):
        vol = compute_nics_volume(_hexagon(), n_points=5, n_normal=3)
        self.assertEqual(len(vol), 5 * 5 * 3)

    def test_defaults_to_a_cube(self):
        self.assertEqual(len(compute_nics_volume(_hexagon(), n_points=4)), 64)

    def test_indices_cover_the_full_lattice(self):
        vol = compute_nics_volume(_hexagon(), n_points=3, n_normal=2)
        self.assertEqual(
            {(p["i"], p["j"], p["k"]) for p in vol},
            {(i, j, k) for i in range(3) for j in range(3) for k in range(2)},
        )

    def test_a_single_normal_layer_sits_on_the_plane_itself(self):
        """n_normal=1 means "just this plane"; landing it at -normal_extent
        would silently shift the whole sheet off the ring."""
        pts = _hexagon()
        _c, _u, _v, n = ring_frame(pts)
        vol = compute_nics_volume(pts, n_points=3, n_normal=1, normal_extent=4.0)
        heights = [(p["pos"] - ring_centroid(pts)) @ n for p in vol]
        np.testing.assert_allclose(heights, np.zeros(9), atol=1e-12)

    def test_the_middle_layer_reproduces_the_2d_grid(self):
        pts = _hexagon()
        grid = compute_nics_grid(pts, n_points=5, extent=3.0, offset=1.0)
        vol = compute_nics_volume(
            pts, n_points=5, extent=3.0, offset=1.0, n_normal=3, normal_extent=2.0
        )
        middle = sorted((p for p in vol if p["k"] == 1), key=lambda p: (p["i"], p["j"]))
        for got, exp in zip(middle, sorted(grid, key=lambda p: (p["i"], p["j"]))):
            np.testing.assert_allclose(got["pos"], exp["pos"], atol=1e-12)

    def test_normal_extent_is_independent_of_the_in_plane_extent(self):
        pts = _hexagon()
        _c, _u, _v, n = ring_frame(pts)
        vol = compute_nics_volume(
            pts, n_points=3, extent=6.0, n_normal=3, normal_extent=1.5
        )
        heights = [(p["pos"] - ring_centroid(pts)) @ n for p in vol]
        self.assertAlmostEqual(min(heights), -1.5, places=9)
        self.assertAlmostEqual(max(heights), 1.5, places=9)

    def test_lab_volume_is_axis_aligned(self):
        vol = compute_nics_volume(
            plane="xy", n_points=3, extent=1.0, n_normal=3, normal_extent=1.0
        )
        zs = sorted({round(float(p["pos"][2]), 9) for p in vol})
        self.assertEqual(zs, [-1.0, 0.0, 1.0])

    def test_all_three_axes_can_differ(self):
        vol = compute_nics_volume(
            _hexagon(),
            n_points=(2, 3),
            extent=(1.0, 2.0),
            n_normal=4,
            normal_extent=3.0,
        )
        self.assertEqual(len(vol), 24)
        self.assertAlmostEqual(max(p["a"] for p in vol), 1.0)
        self.assertAlmostEqual(max(p["b"] for p in vol), 2.0)
        self.assertAlmostEqual(max(p["c"] for p in vol), 3.0)

    def test_the_third_axis_defaults_to_the_second(self):
        vol = compute_nics_volume(_hexagon(), n_points=(2, 3), extent=(1.0, 2.0))
        self.assertEqual(len(vol), 2 * 3 * 3)
        self.assertAlmostEqual(max(p["c"] for p in vol), 2.0)

    def test_rejects_degenerate_parameters(self):
        with self.assertRaises(ValueError):
            compute_nics_volume(_hexagon(), n_points=0)
        with self.assertRaises(ValueError):
            compute_nics_volume(_hexagon(), extent=-1.0)
        with self.assertRaises(ValueError):
            compute_nics_volume(_hexagon(), n_normal=0)
        with self.assertRaises(ValueError):
            compute_nics_volume(_hexagon(), normal_extent=-1.0)


# ---------------------------------------------------------------------------
# Molecule sizing
# ---------------------------------------------------------------------------


@unittest.skipIf(Chem is None, "rdkit not installed")
class TestMoleculeSizing(unittest.TestCase):
    def test_heavy_atom_positions_excludes_hydrogens(self):
        mol = _build("c1ccccc1")
        self.assertEqual(len(heavy_atom_positions(mol)), 6)

    def test_heavy_atom_positions_excludes_ghost_atoms(self):
        """Sizing from probes placed by an earlier run would make the grid grow
        every time it was regenerated."""
        from nics_placer.dialog import _add_bq_atoms

        mol = _build("c1ccccc1")
        grown = _add_bq_atoms(mol, [np.array([50.0, 50.0, 50.0])])
        self.assertEqual(len(heavy_atom_positions(grown)), 6)
        np.testing.assert_allclose(
            molecule_bounds(grown)[1], molecule_bounds(mol)[1], atol=1e-9
        )

    def test_bounds_centre_and_half_extents_cover_the_heavy_atoms(self):
        mol = _build("c1ccccc1")
        center, half = molecule_bounds(mol)
        pts = heavy_atom_positions(mol)
        self.assertTrue(np.all(pts >= center - half - 1e-9))
        self.assertTrue(np.all(pts <= center + half + 1e-9))

    def test_benzene_ring_radius_is_about_1_4_angstrom(self):
        """Sanity-check the sizing against a known geometry rather than only
        against itself."""
        mol = _build("c1ccccc1")
        _center, half = molecule_bounds(mol)
        self.assertAlmostEqual(float(half.max()), 1.4, delta=0.15)

    def test_lab_extents_match_the_bounding_box_plus_margin(self):
        mol = _build("c1ccc2ccccc2c1")
        _center, half = molecule_bounds(mol)
        e1, e2, e3 = axis_extents(mol, "xy", margin=1.0, minimum=0.0)
        for got, exp in zip((e1, e2, e3), half):
            self.assertAlmostEqual(got, float(exp) + 1.0, places=9)

    def test_every_heavy_atom_falls_inside_the_suggested_box(self):
        """The whole point of auto-sizing: nothing may be clipped."""
        mol = _build("c1ccc2ccccc2c1")
        pts = heavy_atom_positions(mol)
        for plane in GRID_PLANES:
            with self.subTest(plane=plane):
                positions = (
                    None
                    if plane in LAB_GRID_PLANES
                    else pts[list(mol.GetRingInfo().AtomRings()[0])]
                )
                center = (
                    molecule_bounds(mol)[0]
                    if plane in LAB_GRID_PLANES
                    else ring_centroid(positions)
                )
                extents = axis_extents(
                    mol,
                    plane,
                    positions=positions,
                    center=center,
                    margin=0.0,
                    minimum=0.0,
                )
                axes = grid_axes(plane, positions)
                rel = pts - center
                for axis, extent in zip(axes, extents):
                    self.assertLessEqual(float(np.abs(rel @ axis).max()), extent + 1e-9)

    def test_ring_frame_extents_are_measured_from_the_ring_centroid(self):
        """A ring-frame grid is centred on its ring, not on the molecule, so an
        off-centre ring needs a longer reach in one direction to still cover
        the far end — larger than the molecule's own half-width, and correctly
        so."""
        mol = _build("c1ccc2ccccc2c1")
        ring = list(mol.GetRingInfo().AtomRings()[0])
        positions = heavy_atom_positions(mol)[ring]
        centroid = ring_centroid(positions)
        extents = axis_extents(
            mol, "parallel", positions=positions, margin=0.0, minimum=0.0
        )
        axes = grid_axes("parallel", positions)
        rel = heavy_atom_positions(mol) - centroid
        for axis, extent in zip(axes, extents):
            self.assertAlmostEqual(float(np.abs(rel @ axis).max()), extent, places=9)

    def test_a_flat_ring_needs_almost_no_room_along_its_own_normal(self):
        """Benzene is planar, so the normal axis of its own ring frame collapses
        to the margin — a useful check that the axes really are the ring's."""
        mol = _build("c1ccccc1")
        positions = heavy_atom_positions(mol)
        _e1, _e2, e3 = axis_extents(
            mol, "parallel", positions=positions, margin=0.0, minimum=0.0
        )
        self.assertLess(e3, 0.05)

    def test_extents_respect_the_minimum(self):
        mol = _build("C")  # methane: one heavy atom, zero-size box
        self.assertTrue(all(e >= 2.0 for e in axis_extents(mol, "xy", margin=0.0)))


class TestCountsForSpacing(unittest.TestCase):
    def test_counts_deliver_at_most_the_requested_step(self):
        for extent in (1.0, 3.0, 4.7):
            for spacing in (0.25, 0.5, 1.0):
                with self.subTest(extent=extent, spacing=spacing):
                    (n,) = counts_for_spacing([extent], spacing)
                    self.assertGreater(n, 1)
                    self.assertLessEqual((2.0 * extent) / (n - 1), spacing + 1e-12)

    def test_an_exact_division_is_not_padded(self):
        self.assertEqual(counts_for_spacing([3.0], 0.5), (13,))

    def test_independent_extents_give_independent_counts(self):
        counts = counts_for_spacing([1.0, 3.0, 5.0], 1.0)
        self.assertEqual(counts, (3, 7, 11))

    def test_a_zero_width_axis_collapses_to_one_point(self):
        self.assertEqual(counts_for_spacing([0.0, 2.0], 1.0), (1, 5))

    def test_rejects_a_non_positive_spacing(self):
        for bad in (0.0, -0.5):
            with self.subTest(spacing=bad):
                with self.assertRaises(ValueError):
                    counts_for_spacing([1.0], bad)


class TestSizingWithoutHeavyAtoms(unittest.TestCase):
    class _EmptyMol:
        """A molecule of nothing but hydrogens — degenerate, but the sizing
        helpers must return numbers rather than raise."""

        class _Conf:
            def GetAtomPosition(self, _i):
                return [0.0, 0.0, 0.0]

        def GetConformer(self):
            return self._Conf()

        def GetAtoms(self):
            return []

    def test_bounds_of_an_empty_selection_are_zero(self):
        center, half = molecule_bounds(self._EmptyMol())
        np.testing.assert_allclose(center, np.zeros(3))
        np.testing.assert_allclose(half, np.zeros(3))

    def test_extents_fall_back_to_the_minimum(self):
        self.assertEqual(axis_extents(self._EmptyMol(), "xy"), (2.0, 2.0, 2.0))


if __name__ == "__main__":
    unittest.main()
