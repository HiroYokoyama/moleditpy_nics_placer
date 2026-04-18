"""
Unit tests for nics_math.py — pure numpy, no Qt/RDKit required.
"""
import math
import sys
import os
import unittest
from unittest.mock import MagicMock

import numpy as np

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), "..")))

from nics_placer.nics_math import (
    NICS1_HEIGHT,
    compute_nics_points,
    get_ring_positions,
    get_rings,
    ring_centroid,
    ring_normal,
)


def _hexagon(z=0.0):
    """Regular hexagon in the z=const plane, unit circumradius."""
    return np.array(
        [[math.cos(i * math.pi / 3), math.sin(i * math.pi / 3), z] for i in range(6)],
        dtype=float,
    )


# ---------------------------------------------------------------------------
# ring_centroid
# ---------------------------------------------------------------------------

class TestRingCentroid(unittest.TestCase):
    def test_origin(self):
        pts = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]], dtype=float)
        np.testing.assert_allclose(ring_centroid(pts), [0, 0, 0], atol=1e-10)

    def test_offset(self):
        pts = np.array([[2, 3, 4], [4, 3, 4], [3, 5, 4]], dtype=float)
        expected = np.array([3, 11 / 3, 4])
        np.testing.assert_allclose(ring_centroid(pts), expected, atol=1e-10)

    def test_benzene(self):
        np.testing.assert_allclose(ring_centroid(_hexagon()), [0, 0, 0], atol=1e-10)


# ---------------------------------------------------------------------------
# ring_normal
# ---------------------------------------------------------------------------

class TestRingNormal(unittest.TestCase):
    def test_xy_plane(self):
        pts = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]], dtype=float)
        n = ring_normal(pts)
        self.assertAlmostEqual(abs(n[2]), 1.0, places=5)

    def test_xz_plane(self):
        pts = np.array([[1, 0, 0], [-1, 0, 0], [0, 0, 1], [0, 0, -1]], dtype=float)
        n = ring_normal(pts)
        self.assertAlmostEqual(abs(n[1]), 1.0, places=5)

    def test_unit_length(self):
        n = ring_normal(_hexagon(z=0.1))
        self.assertAlmostEqual(np.linalg.norm(n), 1.0, places=8)

    def test_benzene_normal_is_z(self):
        n = ring_normal(_hexagon())
        self.assertAlmostEqual(abs(n[2]), 1.0, places=5)


# ---------------------------------------------------------------------------
# compute_nics_points
# ---------------------------------------------------------------------------

class TestComputeNicsPoints(unittest.TestCase):
    def test_nics0_at_centroid(self):
        result = compute_nics_points(_hexagon())
        np.testing.assert_allclose(result["nics0"], [0, 0, 0], atol=1e-5)

    def test_default_height(self):
        result = compute_nics_points(_hexagon())
        d_above = np.linalg.norm(result["nics1_above"] - result["nics0"])
        d_below = np.linalg.norm(result["nics1_below"] - result["nics0"])
        self.assertAlmostEqual(d_above, NICS1_HEIGHT, places=5)
        self.assertAlmostEqual(d_below, NICS1_HEIGHT, places=5)

    def test_custom_height(self):
        result = compute_nics_points(_hexagon(), height=2.0)
        d = np.linalg.norm(result["nics1_above"] - result["nics0"])
        self.assertAlmostEqual(d, 2.0, places=5)

    def test_zero_height_collapses_to_nics0(self):
        result = compute_nics_points(_hexagon(), height=0.0)
        np.testing.assert_allclose(result["nics1_above"], result["nics0"], atol=1e-10)
        np.testing.assert_allclose(result["nics1_below"], result["nics0"], atol=1e-10)

    def test_nics1_symmetric_about_nics0(self):
        result = compute_nics_points(_hexagon())
        mid = (result["nics1_above"] + result["nics1_below"]) / 2
        np.testing.assert_allclose(mid, result["nics0"], atol=1e-10)

    def test_raises_on_fewer_than_3_atoms(self):
        with self.assertRaises(ValueError):
            compute_nics_points(np.array([[0, 0, 0], [1, 0, 0]], dtype=float))

    def test_nics0_is_copy(self):
        result = compute_nics_points(_hexagon())
        result["nics0"][0] = 999.0
        result2 = compute_nics_points(_hexagon())
        self.assertNotAlmostEqual(result2["nics0"][0], 999.0)

    def test_nics1_above_z_positive_for_xy_ring(self):
        result = compute_nics_points(_hexagon())
        # Above should differ from below in one direction
        self.assertFalse(
            np.allclose(result["nics1_above"], result["nics1_below"])
        )

    def test_keys_present(self):
        result = compute_nics_points(_hexagon())
        self.assertIn("nics0", result)
        self.assertIn("nics1_above", result)
        self.assertIn("nics1_below", result)


# ---------------------------------------------------------------------------
# get_ring_positions
# ---------------------------------------------------------------------------

class TestGetRingPositions(unittest.TestCase):
    def _make_mol(self, positions):
        mol = MagicMock()
        conf = MagicMock()
        def get_pos(i):
            p = positions[i]
            # Return a list so [*pos] unpacking works (mirrors real RDKit Point3D iteration)
            return [p[0], p[1], p[2]]
        conf.GetAtomPosition.side_effect = get_pos
        mol.GetConformer.return_value = conf
        return mol

    def test_returns_ndarray(self):
        mol = self._make_mol({0: (0, 0, 0), 1: (1, 0, 0), 2: (0, 1, 0)})
        pts = get_ring_positions(mol, (0, 1, 2))
        self.assertIsInstance(pts, np.ndarray)
        self.assertEqual(pts.shape, (3, 3))

    def test_correct_values(self):
        mol = self._make_mol({0: (1.0, 2.0, 3.0), 1: (4.0, 5.0, 6.0)})
        pts = get_ring_positions(mol, (0, 1))
        np.testing.assert_allclose(pts[0], [1.0, 2.0, 3.0])
        np.testing.assert_allclose(pts[1], [4.0, 5.0, 6.0])


# ---------------------------------------------------------------------------
# get_rings
# ---------------------------------------------------------------------------

class TestGetRings(unittest.TestCase):
    def _make_mol(self, ring_atoms, aromatic_set, has_conformer=True):
        mol = MagicMock()
        if not has_conformer:
            mol.GetConformer.side_effect = Exception("no conformer")
        else:
            mol.GetConformer.return_value = MagicMock()

        ring_info = MagicMock()
        ring_info.AtomRings.return_value = [tuple(ring_atoms)]
        mol.GetRingInfo.return_value = ring_info

        def get_atom(i):
            a = MagicMock()
            a.GetIsAromatic.return_value = i in aromatic_set
            return a
        mol.GetAtomWithIdx.side_effect = get_atom
        return mol

    def test_no_conformer_returns_empty(self):
        mol = self._make_mol([0, 1, 2, 3, 4, 5], {0, 1, 2, 3, 4, 5}, has_conformer=False)
        self.assertEqual(get_rings(mol), [])

    def test_aromatic_ring(self):
        mol = self._make_mol([0, 1, 2, 3, 4, 5], {0, 1, 2, 3, 4, 5})
        rings = get_rings(mol)
        self.assertEqual(len(rings), 1)
        self.assertTrue(rings[0]["is_aromatic"])
        self.assertEqual(rings[0]["atoms"], (0, 1, 2, 3, 4, 5))

    def test_non_aromatic_ring(self):
        mol = self._make_mol([0, 1, 2, 3], set())
        rings = get_rings(mol)
        self.assertEqual(len(rings), 1)
        self.assertFalse(rings[0]["is_aromatic"])

    def test_partially_aromatic_is_not_aromatic(self):
        # Only some atoms aromatic → ring not fully aromatic
        mol = self._make_mol([0, 1, 2, 3, 4, 5], {0, 1, 2})
        rings = get_rings(mol)
        self.assertFalse(rings[0]["is_aromatic"])


if __name__ == "__main__":
    unittest.main()
