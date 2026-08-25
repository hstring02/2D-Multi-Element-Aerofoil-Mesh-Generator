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

### Mesher (`05_Meshing_Script`)

Currently only supports modified .unv outputs intended for use with the XCALibre.jl CFD solver. Generic .unv and OpenFOAM coming soon.

- **Multi-element geometry** — any number of aerofoil elements per case, each with its own Selig-format `.dat` coordinate file, chord, angle of attack and (optional) blunt trailing-edge thickness.
- **Trailing-edge-relative positioning** — each element's `POSITION` is an `(x, y)` offset from the *previous* element's trailing edge, not an absolute world coordinate (the first element is the exception — its `POSITION` is absolute). Moving an upstream element automatically carries every element downstream of it along with it.
- **Config-driven** — every case is defined in a single TOML file (see [`04_Mesh_Input_File`](04_Mesh_Input_File)); no code edits needed for standard runs.
- **Boundary layer (inflation) meshing** — configurable first-layer height, growth rate and total thickness, with automatic fanning around blunt trailing-edge corners.
- **Wake refinement** — an AOA-scaled refinement box generated per element, sized and oriented to capture each element's own shed wake.
- **Curvature-based refinement** — cell size adapts to local surface curvature (e.g. leading edges).
- **XCALibre-ready export** — writes a `.unv` mesh with the Fortran exponent, boundary-group dataset, and node/element numbering fixes required for [XCALibre.jl](https://github.com/mberto79/XCALibre.jl) to read it correctly.
- **Interactive preview** — the generated mesh opens automatically in the Gmsh GUI for inspection.

### CFD Scripts (`07_CFD_Scripts`)

Example [XCALibre.jl](https://github.com/mberto79/XCALibre.jl) run scripts for meshes produced by this repo:

- Steady, incompressible RANS via XCALibre's `Physics`/`Configuration` model.
- Two turbulence closures demonstrated: **k-ω** and **k-ω SST**, with wall functions on the airfoil boundary and a symmetry condition on the farfield.
- Boundary conditions matched directly to the physical groups written by the mesher (`inlet`, `outlet`, `airfoil`, `farfield`).
- Per-field solver configuration (Bicgstab/CG, Jacobi preconditioning, relaxation factors, convergence tolerances).
- CPU multithreaded execution via XCALibre's hardware backend (the same scripts can be retargeted to GPU backends supported by XCALibre).
- `Cl`/`Cd` are normalised against the **sum of every element's chord** (not `flow.REYNOLDS_LENGTH`, which is a separate mesh-sizing/Reynolds-number reference length) — so retuning the mesh doesn't silently change reported force coefficients.

See the [XCALibre.jl repository](https://github.com/mberto79/XCALibre.jl) for full solver documentation.

### Batch Runner (`01_Batch_Run_Input_File`, `02_Batch_Run_Script`)

Sweeps any single field of a case TOML (AOA, velocity, chord, ...) across a
list of values and runs the full mesh → solve pipeline for each point:

- **Config-driven sweeps** — point it at a base case TOML plus a dotted field
  path and list of values (e.g. `foils.AOA.0` over `[-4, -2, 0, 2, 4]`); it
  derives one case file per value automatically.
- **Persistent CFD worker** — meshing runs per point as usual, but all CFD
  solves for the sweep are handed to a single long-lived Julia process, so
  XCALibre's JIT/precompile cost is paid once instead of once per point.
- **Resilient to per-point failure** — a mesh or solve failure for one value
  is logged and the sweep continues; only a hard crash of the Julia worker
  itself stops remaining points.
- **CSV + plot output** — every point's inputs and `Cl`/`Cd`/lift/drag are
  written to a CSV, plus a two-panel PNG plot of `Cl` and `Cd` against the
  swept value.
- **Kept separate from hand-authored cases** — each sweep's derived case
  TOMLs, meshes and results all land under their own `Batch`/`<batch name>`
  subfolders (see the table below), so batch runs never clutter or collide
  with cases you've built by hand.

## 2. File / Folder Structure and Workflow

| Folder | Contents |
|---|---|
| [`01_Batch_Run_Input_File`](01_Batch_Run_Input_File) | `batch_config.toml`: the base case, the dotted field to sweep, and the list of values. |
| [`02_Batch_Run_Script`](02_Batch_Run_Script) | `RUN_SCRIPT.py` sweeps a case TOML field over a list of values, meshing and solving each point and collecting the results into a CSV + plot. |
| [`03_Foils`](03_Foils) | Raw aerofoil coordinates (Selig format, one `.dat` per aerofoil) used to build each element. |
| [`04_Mesh_Input_File`](04_Mesh_Input_File) | Per-case TOML files: element geometry (chord, position, AOA, TE thickness), farfield extents, and all mesh refinement settings. Also holds `Batch/<batch name>/`, the per-point case files a batch sweep derives from the base case. |
| [`05_Meshing_Script`](05_Meshing_Script) | `RUN_SCRIPT.py` is used to make the mesh plus various modules containing meshing functions. |
| [`06_Meshes`](06_Meshes) | Output `.unv` mesh files, named after each TOML file's `title`. Batch-run meshes are kept separate under `Batch/`. |
| [`07_CFD_Scripts`](07_CFD_Scripts) | Example Julia scripts that load a mesh from `06_Meshes` and run a case in XCALibre.jl, plus the batch runner's persistent solve worker (`07_CFD_Scripts/Batch`). |
| [`08_Batch_Run_Results`](08_Batch_Run_Results) | Per-batch `results.csv` + plot, one subfolder per batch name. |

### Workflow

1. Add/confirm your aerofoil coordinate file(s) in `03_Foils` (Selig `.dat` format).
2. Create or edit a case file in `04_Mesh_Input_File` (copy `2_el_wing.toml` or `3_el_wing.toml` as a starting point) — set each element's chord, position, AOA and TE thickness, the farfield box, and the refinement/boundary-layer/wake settings. Remember `POSITION` is trailing-edge-relative for every element after the first (see [Mesher](#mesher-05_meshing_script) above) — it's an offset from the previous element's TE, not a world coordinate.
3. Run the mesher, pointing it at your config file:
   ```
   python RUN_SCRIPT.py ../04_Mesh_Input_File/your_case.toml
   ```
   (run from inside `05_Meshing_Script`; omitting the argument defaults to `2_el_wing.toml`). The mesh is written to `06_Meshes/<title>.unv` and the Gmsh GUI opens automatically so you can inspect it.
4. Point a Julia script in `07_CFD_Scripts` at your new `.unv` file and run it to solve the flow with XCALibre.jl.
5. To sweep a parameter across many points instead of running one case by
   hand, edit `01_Batch_Run_Input_File/batch_config.toml` (base case, field
   to sweep, values) and run:
   ```
   python RUN_SCRIPT.py ../01_Batch_Run_Input_File/batch_config.toml
   ```
   (run from inside `02_Batch_Run_Script`; omitting the argument defaults to
   `01_Batch_Run_Input_File/batch_config.toml`). Results land in
   `08_Batch_Run_Results/<batch name>/`.

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
- **Single-parameter sweeps** — the batch runner varies one dotted TOML field at a time; there's no built-in multi-dimensional (grid) sweep across several fields in one run.
- **No built-in mesh-quality or convergence reporting** — mesh quality is inspected visually in the Gmsh preview, and solver convergence/y+ are checked manually from XCALibre's own output; the pipeline doesn't automate either.
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
- **Multi-parameter batch sweeps** — extend the batch runner beyond a single dotted field to multi-dimensional (grid) sweeps in one run.
- **Automated mesh-quality and convergence reporting** — surface mesh quality metrics and solver convergence/y+ checks directly from the pipeline instead of relying on manual inspection.
- **Gradient based optimisation** - Allow batch runs to efficiently find optimal solutions for multi-variable design spaces.
