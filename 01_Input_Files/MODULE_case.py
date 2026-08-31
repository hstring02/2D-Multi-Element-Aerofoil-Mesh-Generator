import itertools
import sys
from pathlib import Path

import toml

REPO_ROOT = Path(__file__).resolve().parent.parent
INPUT_FILES_DIR = REPO_ROOT / "01_Input_Files"
DEFAULTS_DIR = INPUT_FILES_DIR / "Best_Practice_Parameters"
MESHES_DIR = REPO_ROOT / "04_Meshes"
RESULTS_DIR = REPO_ROOT / "05_CFD_Results"

sys.path.insert(0, str(REPO_ROOT / "02_Tool_Scripts" / "01_Batch_Run_Script"))
from MODULE_batch import (  # noqa: E402
    run_mesher, CFDWorker,
    build_case_toml_multi, write_results_csv, write_results_csv_multi, plot_results,
)

# Sections merged key-by-key against [mesh_settings].DEFAULTS (a case file only
# needs to set the keys it wants to override, e.g. ENABLED = false).
MERGED_SECTIONS = ("mesh_settings", "global_refinement", "boundary_layer", "wake_refinement", "te_refinement")
# Sections that are case-file-only — copied through unchanged, never merged.
PASSTHROUGH_SECTIONS = ("case", "batch", "cfd", "flow", "foils", "farfield")


# ============================================================
# LOAD + MERGE
# ============================================================

def load_and_merge_case(case_path):
    """
    Loads a unified case TOML plus its [mesh_settings].DEFAULTS file and
    merges them per-section, per-key (case file wins). Returns a dict still
    shaped like the unified schema (not yet flattened for the mesher).
    """
    case_data = toml.load(case_path)

    mesh_settings = case_data.get("mesh_settings", {})
    if "DEFAULTS" not in mesh_settings:
        raise ValueError(f"{case_path}: [mesh_settings] is missing DEFAULTS")

    defaults_arg = Path(mesh_settings["DEFAULTS"])
    defaults_path = defaults_arg if defaults_arg.parent != Path(".") else DEFAULTS_DIR / defaults_arg
    if not defaults_path.is_file():
        raise FileNotFoundError(f"mesh_settings.DEFAULTS file not found: {defaults_path}")
    defaults_data = toml.load(defaults_path)

    resolved = {}
    for section in MERGED_SECTIONS:
        resolved[section] = {**defaults_data.get(section, {}), **case_data.get(section, {})}
    for section in PASSTHROUGH_SECTIONS:
        if section in case_data:
            resolved[section] = case_data[section]

    return resolved


def to_flat_mesh_format(resolved):
    """
    Translates a resolved (merged) unified-schema case into the flat schema
    02_Tool_Scripts/02_Meshing_Script/RUN_SCRIPT.py actually reads (top-level
    'title', [visualisation].ENABLED, etc).
    """
    return {
        "title": resolved["case"]["TITLE"],
        "flow": resolved["flow"],
        "foils": resolved["foils"],
        "farfield": resolved["farfield"],
        "visualisation": {"ENABLED": resolved["mesh_settings"]["VISUALISATION"]},
        "global_refinement": resolved["global_refinement"],
        "boundary_layer": resolved["boundary_layer"],
        "wake_refinement": resolved["wake_refinement"],
        "te_refinement": resolved["te_refinement"],
    }


def warn_if_output_exists(*paths):
    for path in paths:
        if path.exists():
            print(f"Warning: overwriting existing output at {path}")


def _write_case_toml(flat, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    title = flat["title"]
    case_toml_path = out_dir / f"{title}.toml"
    warn_if_output_exists(case_toml_path)
    with open(case_toml_path, "w") as f:
        toml.dump(flat, f)
    return case_toml_path, title


# ============================================================
# TYPE = "mesh"
# ============================================================

def run_mesh_case(case_path):
    resolved = load_and_merge_case(case_path)
    flat = to_flat_mesh_format(resolved)

    case_toml_path, title = _write_case_toml(flat, INPUT_FILES_DIR / "Cases")

    warn_if_output_exists(MESHES_DIR / "Cases" / f"{title}.unv")
    mesh_path = run_mesher(case_toml_path, title, dest_subdir="Cases")

    if resolved["mesh_settings"].get("SAVE_MESH", True):
        print(f"Mesh saved to {mesh_path}")
    else:
        mesh_path.unlink()
        print(f"SAVE_MESH is false — deleted {mesh_path} after generation (dry run).")

    return mesh_path


# ============================================================
# TYPE = "cfd"
# ============================================================

def run_cfd_case(case_path):
    resolved = load_and_merge_case(case_path)
    flat = to_flat_mesh_format(resolved)

    case_toml_path, title = _write_case_toml(flat, INPUT_FILES_DIR / "Cases")

    warn_if_output_exists(MESHES_DIR / "Cases" / f"{title}.unv")
    mesh_path = run_mesher(case_toml_path, title, dest_subdir="Cases")

    worker = CFDWorker()
    try:
        result = worker.solve(str(mesh_path), str(case_toml_path))
    finally:
        worker.quit()

    print(f"Cl={result['Cl']:.4g}  Cd={result['Cd']:.4g}  "
          f"lift={result['lift']:.4g}N  drag={result['drag']:.4g}N")

    results_path = RESULTS_DIR / "Cases" / title / "result.csv"
    row = {
        "index": 0, "sweep_key": "", "value": "", "status": "ok",
        "Cl": result["Cl"], "Cd": result["Cd"], "lift": result["lift"], "drag": result["drag"],
        "mesh": str(mesh_path), "error": "",
    }
    write_results_csv([row], results_path)
    print(f"Wrote {results_path}")

    if not resolved["mesh_settings"].get("SAVE_MESH", True):
        mesh_path.unlink()
        print(f"SAVE_MESH is false — deleted {mesh_path} after solving.")

    return result


# ============================================================
# TYPE = "batch"
# ============================================================

def run_batch_case(case_path):
    resolved = load_and_merge_case(case_path)
    flat = to_flat_mesh_format(resolved)

    sweep_keys = resolved["batch"]["sweep_keys"]
    sweep_values = resolved["batch"]["sweep_values"]

    for key in sweep_keys:
        if key.startswith("case.") or key.startswith("mesh_settings."):
            raise ValueError(
                f"Cannot sweep '{key}': [case] and [mesh_settings] fields are consumed "
                "before the sweep is applied (e.g. case.TITLE -> title, "
                "mesh_settings.VISUALISATION -> visualisation.ENABLED), so overriding "
                "them per sweep point would have no effect on the mesh/CFD run."
            )

    grid_points = list(itertools.product(*sweep_values))

    batch_name = Path(case_path).stem
    runs_dir = INPUT_FILES_DIR / "Batch" / batch_name
    results_dir = RESULTS_DIR / batch_name
    save_mesh = resolved["mesh_settings"].get("SAVE_MESH", True)

    print(f"Batch '{batch_name}': grid sweep over {sweep_keys} -> {len(grid_points)} point(s)")

    worker = CFDWorker()
    rows = []

    try:
        for i, combo in enumerate(grid_points):
            overrides = dict(zip(sweep_keys, combo))
            row = {
                "index": i, **overrides, "status": "ok",
                "Cl": "", "Cd": "", "lift": "", "drag": "", "mesh": "", "error": "",
            }
            try:
                case_toml_path, title = build_case_toml_multi(flat, overrides, runs_dir, i)
                mesh_path = run_mesher(case_toml_path, title, dest_subdir="Batch")
                row["mesh"] = str(mesh_path)

                result = worker.solve(str(mesh_path), str(case_toml_path))
                row.update(result)

                if not save_mesh:
                    mesh_path.unlink()

                print(f"[{i + 1}/{len(grid_points)}] {overrides} "
                      f"-> Cl={result['Cl']:.4g} Cd={result['Cd']:.4g}")
            except Exception as e:
                row["status"] = "error"
                row["error"] = str(e)
                print(f"[{i + 1}/{len(grid_points)}] {overrides} -> FAILED: {e}")

            rows.append(row)

            # See 02_Tool_Scripts/01_Batch_Run_Script/RUN_SCRIPT.py for why a dead worker
            # stops the sweep instead of failing through remaining points.
            if worker.proc.poll() is not None:
                remaining = grid_points[i + 1:]
                for j, rc in enumerate(remaining, start=i + 1):
                    ov = dict(zip(sweep_keys, rc))
                    rows.append({
                        "index": j, **ov, "status": "skipped",
                        "Cl": "", "Cd": "", "lift": "", "drag": "", "mesh": "",
                        "error": "CFD worker process had already exited",
                    })
                print(f"Julia CFD worker exited — skipping {len(remaining)} remaining point(s).")
                break
    finally:
        worker.quit()

    csv_path = results_dir / "results.csv"
    write_results_csv_multi(rows, sweep_keys, csv_path)
    print(f"\nWrote {csv_path}")

    n_ok = sum(1 for r in rows if r["status"] == "ok")
    if n_ok > 0 and len(sweep_keys) == 1:
        plot_rows = [{**r, "value": r[sweep_keys[0]]} for r in rows]
        plot_path = results_dir / "results.png"
        plot_results(plot_rows, sweep_keys[0], plot_path)
        print(f"Wrote {plot_path}")
    elif n_ok > 0:
        print("N-dimensional sweep — skipping plot (only 1-key sweeps are plotted).")
    else:
        print("No successful points — skipped plot.")

    n_bad = len(rows) - n_ok
    if n_bad:
        print(f"{n_bad}/{len(rows)} point(s) failed or were skipped — "
              f"see the 'status'/'error' columns in {csv_path.name}")

    return rows


# ============================================================
# DISPATCH
# ============================================================

def run(case_path):
    case_path = Path(case_path)
    resolved = load_and_merge_case(case_path)
    case_type = resolved["case"]["TYPE"]

    if case_type == "mesh":
        return run_mesh_case(case_path)
    elif case_type == "cfd":
        return run_cfd_case(case_path)
    elif case_type == "batch":
        return run_batch_case(case_path)
    else:
        raise ValueError(f"Unknown [case].TYPE: {case_type!r} (expected 'mesh', 'cfd', or 'batch')")
