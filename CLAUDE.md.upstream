# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

maskLib is a Python library for drawing superconducting qubit/circuit mask designs and writing them out as DXF files. All dimensions are in microns. This folder is a downloaded snapshot, **not a git repository**.

## Environment

Use the local virtualenv at `masklib/` (Python 3.14, maskLib already installed in editable mode):

```
masklib\Scripts\python.exe <script.py>
```

Do not install packages into the global Python — use `masklib\Scripts\pip.exe`. Core dependencies: `dxfwrite` (the main drawing engine), `numpy`, `matplotlib`, `ezdxf`, `klayout` (DXF→GDS conversion in `maskLib.gdsExport`); optional: `gdspy`, `opencv-python`. `setup_env.sh` documents an alternative conda-based setup.

## Running and verifying

There is no test suite or linter. Verify changes by running an example script and checking that the DXF output is generated without errors:

```
masklib\Scripts\python.exe example\CPWResonatorExample.py
```

Scripts write `.dxf` files to a path given when constructing the `Wafer` (usually `DXF/`). Root-level scripts (`TmonDoseArrayPrathu.py`, `junction array.py`, etc.) are personal design scripts; `example/` holds the canonical usage examples (duplicated under `src/maskLib/example/`).

## Architecture

The package lives in `src/maskLib/` (setuptools src layout, `pip install -e .`).

Core object model in `MaskLib.py`:
- **`Wafer`** — top-level drawing: owns the dxfwrite drawing, layer setup (`SetupLayers`), chip grid layout, and dicing borders. Typical script flow: create `Wafer` → `SetupLayers([...])` → `w.init()` → `w.DicingBorder()` → define/draw chips → `w.populate()` → `w.save()`.
- **`Chip`** — one chip design; subclasses define standard sizes/launcher positions (`Chip7mm`, `Chip10mm`, `ChipLL_*`). A chip carries a `defaults` dict (e.g. `{'w','s','radius','r_out','r_ins','curve_pts'}`) that component functions fall back on for CPW width, gap, bend radius, etc. `chip.save(wafer, ...)` registers it with the wafer.
- **`Structure`** — a drawing cursor (position + direction). Component functions take `(chip, structure, ...)` and advance the structure as they draw, so calls chain: `CPW_launcher(...)` → `CPW_straight(...)` → `CPW_bend(...)` → `CPW_wiggles(...)`. Chips pre-define structures at launcher positions in `self.structures`.

Component libraries build on this: `microwaveLib.py` (CPW transmission lines — the largest), `qubitLib.py`, `junctionLib.py`, `resonatorLib.py`, `fluxoniumLib.py`, `dcLib.py`, `markerLib.py`, `mmWaveLib.py`. Low-level shape entities (e.g. `SolidPline`, `SkewRect`) are in `Entities.py`; geometry helpers (`curveAB`, `cornerRound`, `doMirrored`, `transformedQuadrants`) in `utilities.py`. `arrayLib.py` provides `Sweep3D`, the 3D parameter-sweep engine for dose/geometry arrays (tiled field grids, auto-generated per-dose layers, field labels like `A0103`, xlsx export with per-parameter gradient maps) — see `junction array.py` for usage.

Drawing engine notes:
- Primary engine is `dxfwrite`; modules set `const.POLYLINE_3D_POLYLINE = 0` at import to force 2D polylines — keep that when adding modules.
- `ezdxf` is used secondarily (text-to-path, DXF reading).
- Shapes are drawn twice conceptually: outline (frame) and solid fill, controlled by the wafer's `frame`/`solid`/`multiLayer` flags and per-call `bgcolor`.

## Conventions and pitfalls

- Deprecated code is kept around: `*_old.py` modules, `fluxoniumLib_newandbroken.py`, and deprecated top-level functions in `MaskLib.py` (marked with "Deprecated - use X instead" comments). Don't extend these — use the replacements in `markerLib`, `utilities`, and `Entities`.
- `MaskLib.py` has duplicated import blocks at the top (historical accident); harmless but don't replicate the pattern.
- Layers are referenced by name strings (e.g. `layer='MARKERS'`) that must match names passed to `SetupLayers`.
- `Chip.add()` shifts objects that have a `.points` attribute (e.g. `Entities.SolidPline`) by `chip.origin_offset` (−chipsize/2 each axis) and grid-snaps them; plain `dxf.polyline`/`dxf.rectangle` objects are added as-is. Scripts that position `SolidPline`-based components therefore pre-compensate with `+chipsize/2` — when mixing entity types at the same location, apply `chip.origin_offset` manually to the plain entities (see `JunctionWithLeads` in `junction array.py`).
- The layer name `703/0` is invalid per the DXF spec (`/` not allowed); `dxfwrite` writes it but strict readers like `ezdxf` refuse to open the file.
