# NICS Placer — Test Suite

43 tests across 2 files. Core mathematical functions are tested without any
stubs. Plugin integration tests use a stub `PluginContext` and optionally the
real `PluginContext` when the main app is present as a sibling repo.

---

## Running the tests

```bash
# Full suite
python -m pytest tests/ -v

# Single file
python -m pytest tests/test_nics_math.py -v
```

---

## Test files

| File | Tests | Area |
|---|---|---|
| `test_nics_math.py` | ~35 | Ring centroid, ring normal, NICS point placement, ring position utilities |
| `test_plugin_integration.py` | ~8 | PluginContext registration contract |

---

## Test files — detailed

### `test_nics_math.py` — Core NICS mathematics

Pure Python calculations — no Qt, no RDKit stubs required.

| Class | What is tested |
|---|---|
| `TestRingCentroid` | Centroid of triangle, square, irregular polygon; single-atom edge case |
| `TestRingNormal` | Normal vector of planar ring; vector is unit length; perpendicular to ring plane; sign consistent with right-hand rule |
| `TestComputeNicsPoints` | NICS(0) point at centroid; NICS(1) offset 1 Å above and below plane; custom distance; coplanar rings share plane normal |
| `TestGetRingPositions` | Atom positions extracted correctly from molecule; ordering matches ring atom list |
| `TestGetRings` | Detects 6-membered rings; detects 5-membered rings; fused ring systems return each ring separately; single ring; no rings → empty list |

---

### `test_plugin_integration.py` — PluginContext contract

| Class | What is tested |
|---|---|
| `TestMetadata` | `PLUGIN_NAME`, `PLUGIN_VERSION`, `PLUGIN_AUTHOR`, `PLUGIN_DESCRIPTION` present and non-empty |
| `TestSaveHandler` | Save handler returns a dict with expected keys; round-trip preserves data |
| `TestLoadHandler` | Load handler restores state from saved dict; empty/None data does not raise |
| `TestResetHandler` | Reset clears NICS point list without raising |
| `TestRoundtrip` | Save → load produces identical state |

---

## CI

| Job | Python matrix | Tests run |
|---|---|---|
| `test` | 3.11, 3.12 | Full suite |
