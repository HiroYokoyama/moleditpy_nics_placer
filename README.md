# MoleditPy NICS Placer Plugin

![NICS Placer](img/main.png)

A [MoleditPy](https://github.com/HiroYokoyama/python_molecular_editor) plugin that detects rings in the loaded molecule and places ghost atoms at NICS(0) and NICS(1) probe positions for use in NICS (Nucleus-Independent Chemical Shift) calculations with ORCA or Gaussian.

## Features

- Detects **all rings** (aromatic and non-aromatic) from the 3D structure
- Computes NICS probe positions:
  - **NICS(0)** — ring centroid (in-plane)
  - **NICS(1)±** — ±1 Å above and below the ring plane (SVD best-fit plane)
- Interactive **3-state sphere preview** in the 3D viewport:

  | Colour | Meaning | Interaction |
  |--------|---------|-------------|
  | Yellow (semi-transparent) | Available, not staged | Click → turns red |
  | Red (semi-transparent) | Staged for placement | Click → turns yellow; Apply → places atom |
  | Green | Already placed in molecule | Clear All to remove |

- **Ghost atom symbol selector** — choose per project or set a persistent default:

  | Symbol | Software |
  |--------|----------|
  | `Bq` | Gaussian; also valid in ORCA (default) |
  | `H:` | ORCA native ghost atom notation |

- Table view with per-ring status (size, aromaticity, placed/staged count)
- Helper buttons: Stage NICS(0) / Stage NICS(1)± for selected rings, Place All, Clear All, Refresh
- **Auto-refresh** — when the molecule changes (load, undo/redo), rings and spheres update automatically
- Compatible with **ORCA Input Generator Pro** and **Gaussian Input Generator Neo** via the shared `custom_symbol` atom property

## Workflow

1. Load a molecule with 3D coordinates in MoleditPy.
2. Open **3D Edit → NICS Placer…**
3. Yellow spheres appear at all NICS(0), NICS(1)+, and NICS(1)− positions for every ring.
4. Click yellow spheres to stage them (turns red), or use the **Stage** buttons for bulk selection.
5. Select ghost atom symbol (`Bq` or `H:`) from the combo box.
6. Press **Apply (place red Bq)** to insert ghost atoms at all staged positions.
7. Open **ORCA Input Generator Pro** or **Gaussian Input Generator Neo** — ghost labels appear automatically in the coordinate block.
8. Run NICS calculation.

## Settings & Persistence

### Plugin setting (`settings.json`)
The selected ghost atom symbol is saved to `nics_placer/settings.json` whenever it is changed. This is the **user default** — it persists across all sessions and documents.

### Project setting (`.moleditpy` project file)
When a project is saved, the current ghost symbol and all placed ghost atom indices are stored in the project file. Loading a project restores both, overriding the plugin default for that session. Closing the project (File → New) reverts to the plugin default from `settings.json`.

## Installation

Copy the `nics_placer/` folder into your MoleditPy plugins directory:

```
moleditpy_nics_placer/
    nics_placer/
        __init__.py      ← plugin entry point
        dialog.py        ← NicsPlacerDialog (PyQt6 + PyVista)
        nics_math.py     ← pure numpy ring geometry
        settings.json    ← auto-created on first symbol change (gitignored)
```

## Requirements

- MoleditPy ≥ v3 (V3 plugin API)
- PyQt6
- RDKit
- PyVista
- numpy

## Running Tests

The test suite runs fully headless (no Qt, RDKit, or PyVista required):

```bash
cd moleditpy_nics_placer
python -m pytest tests/ -v
```

| Test file | Coverage |
|-----------|----------|
| `test_nics_math.py` | Pure geometry: centroid, SVD normal, NICS point computation, ring extraction |
| `test_plugin_integration.py` | Plugin contract: save/load/reset handlers, ghost label persistence, symbol persistence |

## Implementation Notes

### Ring plane (best-fit plane)

The ring normal is computed via SVD of the mean-centred atom positions:

```
centered = positions - mean(positions)
U, S, Vt = SVD(centered)
normal = Vt[-1]       # last singular vector = direction of minimum variance
normal /= ‖normal‖
```

This gives the least-squares best-fit plane normal, robust for all planar and near-planar rings.

### Ghost atom convention

Ghost atoms are `rdkit.Chem.Atom(0)` (atomic number 0, dummy atom) with `SetProp("custom_symbol", symbol)`. The `custom_symbol` property is the shared convention used by MoleditPy's XYZ Editor, ORCA Input Generator Pro, and Gaussian Input Generator Neo.

## Version

1.0.0 — HiroYokoyama
