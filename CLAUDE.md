# CLAUDE.md — maskLib fork (Eddie)

## What this repo is

A fork of `agrm930/masklib`, a Python toolkit for superconducting device chip layout
and EBL mask design. Output is primarily DXF (via `dxfwrite`-style entity calls), with
downstream conversion to GDSII for a Heidelberg DWL 66+ laser writer. Used for CPS/CPW
transmission line structures, double-Y baluns, resonators, and related microwave
geometry.

- `origin` = my fork; `upstream` = `agrm930/masklib`. Keep `main` a clean mirror of
  upstream; all custom work goes on `eddie-dev` (or feature branches off it).
- Custom extensions and scripts live in this fork. Do not assume file layout matches
  upstream exactly.

## Environment rules (important — history of split-brain issues)

- This project uses a local venv at `.venv/`. **Before running any `pip install`,
  verify both `which python` and `which pip` resolve inside `.venv/`.** pyenv shims
  have previously hijacked `python` after venv activation, causing pip to write
  metadata into the venv while module files landed in pyenv's site-packages. If the
  two `which` results disagree, stop and fix the environment before installing.
- maskLib must be installed in **editable mode** (`pip install -e .`) from this clone.
  Never install maskLib from any other path.
- There are other copies of maskLib on the network (NAS-hosted) that have caused
  import shadowing before. After any environment change, sanity-check with:
  `python -c "import maskLib; print(maskLib.__file__)"` — the printed path must be
  inside this repo. If it isn't, do not proceed; diagnose sys.path/shadowing first.
- Do not add `sys.path` hacks to make imports work. Fix the install instead.
- Keep `.venv/` in `.gitignore`. Never commit environment directories or generated
  DXF/GDS output unless explicitly asked.

## Version pinning

- Declare loose dependencies in `pyproject.toml`; keep `requirements.txt` as a frozen
  lockfile (`pip freeze`). Update the lockfile deliberately, not as a side effect.
- Related projects in my stack are sensitive to version drift (e.g., pyEPR work pins
  numpy 1.26.4 and pandas 2.0.3). Before bumping numpy/pandas/shapely here, check it
  won't break interop with scripts shared across projects.

## maskLib API patterns

- Entities are added to chips via `self.add(dxf.entity(...))` — follow existing
  patterns in the codebase rather than inventing new drawing paths.
- Layer selection goes through `wafer.lyr()`; respect the existing layer-management
  conventions.
- CPW gap geometry uses an **XOR layer strategy**: gaps are drawn on a dedicated layer
  and boolean-subtracted in KLayout. Don't "simplify" this into direct subtraction in
  Python — the XOR-in-KLayout step is intentional.
- Positive-draw CPW functions take explicit ground plane widths; keep that convention
  for new primitives.

## Output / fabrication constraints

- The Heidelberg DWL 66+ GDSII converter **rejects single-vertex PATH elements**.
  Any code path that could emit degenerate (zero-length / single-vertex) paths must
  guard against it. When generating or cleaning GDS files, check for and remove these.
- Preserve exact dimensions. Never silently round, rescale, or "clean up" coordinates
  in geometry code — micron-level and sub-micron features are load-bearing.
- When modifying drawing functions, prefer adding a test/demo script that renders the
  primitive to DXF over asserting correctness from code reading alone.

## Style and workflow preferences

- Ask before large refactors; prefer plan mode for multi-file changes.
- Small, focused commits with descriptive messages.
- When a geometry function changes, note in the commit message whether output DXF for
  existing designs is expected to change.

## Future direction

- gdsfactory may be added to this workflow (as a sibling project, not inside this
  repo). If asked to integrate: it works natively in GDSII, so it may eventually
  replace the DXF→converter step. Check numpy pin compatibility before installing it
  into this venv.
