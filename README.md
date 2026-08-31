# 2D Multi-Element Wing Optimisation Tool

## Contents

1. [Features / Capabilities](#1-features--capabilities)
2. [File / Folder Structure and Workflow](#2-file--folder-structure-and-workflow)
3. [Pre-requisites](#3-pre-requisites)
4. [Limitations](#4-limitations)
5. [Examples](#5-examples)
6. [Future Work](#6-future-work)

A Python-based mesh generator that uses [Gmsh](https://gmsh.info/) to build 2D CFD grids for multi-element wings (e.g. main element + flap/slat combinations), driven entirely by simple TOML config files. Meshes are exported in a UNV2 format specifically post-processed for direct compatibility with the [XCALibre.jl](https://github.com/mberto79/XCALibre.jl) finite volume CFD solver, and the repo includes example Julia scripts for running RANS cases on the generated grids.

## 1. Features / Capabilities

### Mesher (`02_Tool_Scripts/02_Meshing_Script`)

Currently only supports modified .unv outputs intended for use with the XCALibre.jl CFD solver. Generic .unv and OpenFOAM coming soon.

- **Multi-element geometry** — any number of aerofoil elements per case, each with its own Selig-format `.dat` coordinate file, chord, angle of attack and (optional) blunt trailing-edge thickness.
- **Trailing-edge-relative positioning** — each element after the first is placed relative to the *previous* element's trailing edge via `OVERLAP` (how far its leading edge sits upstream of that TE; positive = overlapping) and `GAP_HEIGHT` (the vertical offset from that TE). The first element has no predecessor and is always placed at the world origin. Moving an upstream element automatically carries every element downstream of it along with it.
- **Config-driven** — every case is defined in a single TOML file (see [`01_Input_Files`](01_Input_Files)); no code edits needed for standard runs.
- **Boundary layer (inflation) meshing** — configurable first-layer height, growth rate and total thickness, with automatic fanning around blunt trailing-edge corners.
- **Wake refinement** — an AOA-scaled refinement box generated per element, sized and oriented to capture each element's own shed wake.
- **Curvature-based refinement** — cell size adapts to local surface curvature (e.g. leading edges).
- **XCALibre-ready export** — writes a `.unv` mesh with the Fortran exponent, boundary-group dataset, and node/element numbering fixes required for [XCALibre.jl](https://github.com/mberto79/XCALibre.jl) to read it correctly.
- **Interactive preview** — the generated mesh opens automatically in the Gmsh GUI for inspection.

### CFD Scripts (`02_Tool_Scripts/03_CFD_Scripts`)

Example [XCALibre.jl](https://github.com/mberto79/XCALibre.jl) run scripts for meshes produced by this repo:

- Steady, incompressible RANS via XCALibre's `Physics`/`Configuration` model.
- Two turbulence closures demonstrated: **k-ω** and **k-ω SST**, with wall functions on the airfoil boundary and a symmetry condition on the farfield.
- Boundary conditions matched directly to the physical groups written by the mesher (`inlet`, `outlet`, `airfoil`, `farfield`).
- Per-field solver configuration (Bicgstab/CG, Jacobi preconditioning, relaxation factors, convergence tolerances).
- CPU multithreaded execution via XCALibre's hardware backend (the same scripts can be retargeted to GPU backends supported by XCALibre).
- `Cl`/`Cd` are normalised against the **sum of every element's chord** (not `flow.REYNOLDS_LENGTH`, which is a separate mesh-sizing/Reynolds-number reference length) — so retuning the mesh doesn't silently change reported force coefficients.

See the [XCALibre.jl repository](https://github.com/mberto79/XCALibre.jl) for full solver documentation.

### Batch Runner (`01_Input_Files`, `02_Tool_Scripts/01_Batch_Run_Script`)

Sweeps one or more fields of a case TOML (AOA, velocity, chord, ...) across a
cartesian grid of values and runs the full mesh → solve pipeline for every
combination:

- **Config-driven grid sweeps** — set `[batch].sweep_keys` (one or more dotted
  field paths, e.g. `foils.AOA.0`) and a matching `sweep_values` list per key;
  every combination across all keys is run (a 2-key, 4-value-each sweep is 16
  points). It derives one case file per grid point automatically.
- **Persistent CFD worker** — meshing runs per point as usual, but all CFD
  solves for the sweep are handed to a single long-lived Julia process, so
  XCALibre's JIT/precompile cost is paid once instead of once per point.
- **Resilient to per-point failure** — a mesh or solve failure for one point
  is logged and the sweep continues; only a hard crash of the Julia worker
  itself stops remaining points.
- **CSV + plot output** — every point's inputs and `Cl`/`Cd`/lift/drag are
  written to a CSV (one column per swept key); a single-key sweep also gets a
  two-panel PNG plot of `Cl` and `Cd` against the swept value.
- **Kept separate from hand-authored cases** — each sweep's derived case
  TOMLs, meshes and results all land under their own `Batch`/`<batch name>`
  subfolders (see the table below), so batch runs never clutter or collide
  with cases you've built by hand.

## 2. File / Folder Structure and Workflow

| Folder | Contents |
|---|---|
| [`01_Input_Files`](01_Input_Files) | `RUN_CASE.py`: the single entry point — reads a unified case TOML's `[case].TYPE` (`mesh`/`cfd`/`batch`) and routes it to the right tool script. Hand-authored case TOMLs live here too, alongside `Best_Practice_Parameters/` (shared mesh-defaults TOML), `Cases/` (resolved single-run case files), and `Batch/<batch name>/` (per-point case files a batch sweep derives). |
| [`02_Tool_Scripts`](02_Tool_Scripts) | The three pipeline stages `RUN_CASE.py` drives, each in its own numbered subfolder: `01_Batch_Run_Script`, `02_Meshing_Script`, `03_CFD_Scripts`. Each `RUN_SCRIPT.py` (or `cfd_worker.jl`) can still be run directly for dev/debugging. |
| [`03_Foils`](03_Foils) | Raw aerofoil coordinates (Selig format, one `.dat` per aerofoil) used to build each element. |
| [`04_Meshes`](04_Meshes) | Output `.unv` mesh files, named after each TOML file's `title`. Batch- and single-run meshes are kept separate under `Batch/`/`Cases/`. |
| [`05_CFD_Results`](05_CFD_Results) | `results.csv` + plot for a batch sweep, or `result.csv` for a single `TYPE="cfd"` run — one subfolder per batch name / case title. |
| [`06_PostProc`](06_PostProc) | Post-processing notebooks/scripts for plotting results. |

### Workflow

1. Add/confirm your aerofoil coordinate file(s) in `03_Foils` (Selig `.dat` format).
2. Create or edit a unified case file in `01_Input_Files` (copy one of the `EXAMPLE_*.toml` files as a starting point) — set `[case].TYPE` (`mesh`, `cfd`, or `batch`), each element's chord, overlap/gap height, AOA and TE thickness, the farfield box, and any mesh-parameter overrides (everything not overridden is inherited from `01_Input_Files/Best_Practice_Parameters/`). Remember `OVERLAP`/`GAP_HEIGHT` are trailing-edge-relative for every element after the first (see [Mesher](#mesher-02_tool_scripts02_meshing_script) above) — one entry per downstream element, offset from the previous element's TE, not a world coordinate.
3. Run it:
   ```
   python RUN_CASE.py your_case.toml
   ```
   (run from inside `01_Input_Files`). `RUN_CASE.py` meshes (and, for `cfd`/`batch`, solves) automatically based on `[case].TYPE`:
   - `TYPE="mesh"` — writes the mesh to `04_Meshes/Cases/<title>.unv`.
   - `TYPE="cfd"` — meshes, solves one point via XCALibre.jl, and writes `05_CFD_Results/Cases/<title>/result.csv`.
   - `TYPE="batch"` — cartesian-expands `[batch].sweep_keys`/`sweep_values` (one or more dotted fields, each with its own list of values) into every combination, meshing and solving each point, and writes `05_CFD_Results/<case file name>/results.csv` (+ a plot, for a single-key sweep).
4. For direct control over one pipeline stage (e.g. iterating on mesh settings without re-solving), the underlying scripts in `02_Tool_Scripts/*/RUN_SCRIPT.py` can still be run by hand — see the comments at the top of each for usage.

## 3. Pre-requisites

- **[Python](https://www.python.org/)** 3.9+, with:
  - [`gmsh`](https://pypi.org/project/gmsh/)
  - [`toml`](https://pypi.org/project/toml/)
  - [`matplotlib`](https://pypi.org/project/matplotlib/) (batch runner result plots only)
  ```
  pip install gmsh toml matplotlib
  ```
- **[Gmsh](https://gmsh.info/)** — installed automatically via the `gmsh` Python package above (no separate install needed).
- **[Julia](https://julialang.org/)** 1.10+, with:
  - [`XCALibre.jl`](https://github.com/mberto79/XCALibre.jl)
  - [`Plots.jl`](https://github.com/JuliaPlots/Plots.jl)
  ```julia
  using Pkg
  Pkg.add(["XCALibre", "Plots"])
  ```

## 4. Limitations

- **2D only** — no spanwise/3D effects, so finite-wing behaviour (tip vortices, spanwise flow) isn't captured.
- **XCALibre.jl-specific export** — the mesher currently only writes the modified `.unv` format XCALibre.jl expects; generic `.unv` and OpenFOAM export aren't supported yet.
- **Global, not local, boundary-layer sizing** — the first inflation-layer cell height comes from one flat-plate correlation evaluated at freestream velocity and a single reference length, applied uniformly around every element. It isn't locally resolved, so regions of locally accelerated flow (e.g. a leading-edge suction peak, or a slot between elements at high incidence) can need a finer target `Y_PLUS` than the reported value guarantees — check y+ contours after solving rather than trusting the single number blindly.
- **Batch sweep plots are 1-D only** — `RUN_CASE.py` supports a multi-key cartesian grid sweep, but only a single-key sweep gets a `Cl`/`Cd` plot; higher-dimensional sweeps produce `results.csv` only.
- **No built-in mesh-quality or convergence reporting** — mesh quality is inspected visually in the Gmsh preview, and solver convergence/y+ are checked manually from XCALibre's own output; the pipeline doesn't automate either. `[cfd].LIFT_CONVERGENCE`/`DRAG_CONVERGENCE` are accepted in a case file but not yet enforced by the solver.
- **Script-driven, no GUI** — every step (meshing, solving, batch sweeps) is run from the command line against TOML/Julia files; there's no packaged app or graphical front end.
- **Julia + XCALibre.jl are a separate install** — this repo doesn't bundle or auto-install them; see [Pre-requisites](#3-pre-requisites).

## 5. Examples

### 2-element wing

![2-element mesh](docs/2_el_wing_mesh.png)
Example **[gmsh](https://pypi.org/project/gmsh/)** preview window.

### 3-element wing

![3-element mesh](docs/3_el_wing_mesh.png)
Example **[gmsh](https://pypi.org/project/gmsh/)** preview window.

### CFD results

![3-element wing CFD](docs/3_el_wing_CFD.png)
Contour made with **[ParaView](https://www.paraview.org/)**.

## 6. Future Work

- **Generic `.unv` and OpenFOAM export** — broaden the mesher beyond the XCALibre.jl-specific output it currently writes.
- **Locally-resolved boundary-layer sizing** — replace the single global first-layer-height correlation with a per-point value driven by an actual local velocity estimate (e.g. a lightweight panel method), instead of one freestream-based number applied uniformly around every element.
- **Multi-dimensional sweep plotting** — extend batch-sweep plotting beyond the single-key case to visualise a multi-key cartesian grid.
- **Automated mesh-quality and convergence reporting** — surface mesh quality metrics and solver convergence/y+ checks directly from the pipeline (wiring up the `[cfd].LIFT_CONVERGENCE`/`DRAG_CONVERGENCE` fields already accepted in a case file) instead of relying on manual inspection.
- **Gradient based optimisation** - Allow batch runs to efficiently find optimal solutions for multi-variable design spaces.
