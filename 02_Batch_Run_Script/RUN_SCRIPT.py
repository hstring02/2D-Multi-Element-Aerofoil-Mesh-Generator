# ============================================================
# LOAD FUNCTIONS
# ============================================================

import sys
from pathlib import Path

import toml

from MODULE_batch import (
    REPO_ROOT,
    build_case_toml, run_mesher, CFDWorker,
    write_results_csv, plot_results,
)

# ============================================================
# LOAD BATCH CONFIG
# ============================================================

BATCH_INPUT_DIR = REPO_ROOT / "01_Batch_Run_Input_File"
DEFAULT_CONFIG = "batch_config.toml"    # Defaults to this if no command-line argument is provided.

if len(sys.argv) > 1:
    arg_path = Path(sys.argv[1])
    CONFIG_PATH = arg_path if arg_path.parent != Path(".") else BATCH_INPUT_DIR / arg_path
else:
    CONFIG_PATH = BATCH_INPUT_DIR / DEFAULT_CONFIG

if not CONFIG_PATH.is_file():
    raise FileNotFoundError(f"Batch config file not found: {CONFIG_PATH}")

batch_config = toml.load(CONFIG_PATH)
batch_name = CONFIG_PATH.stem

base_case_path = REPO_ROOT / batch_config["base_case"]
if not base_case_path.is_file():
    raise FileNotFoundError(f"base_case file not found: {base_case_path}")
base_data = toml.load(base_case_path)

sweep_key = batch_config["sweep_key"]
sweep_values = batch_config["sweep_values"]

runs_dir = REPO_ROOT / "04_Mesh_Input_File" / "Batch" / batch_name
results_dir = REPO_ROOT / "08_Batch_Run_Results" / batch_name

# ============================================================
# RUN SWEEP
# ============================================================

print(f"Batch '{batch_name}': sweeping {sweep_key} over {len(sweep_values)} values")

worker = CFDWorker()
rows = []

try:
    for i, value in enumerate(sweep_values):
        row = {
            "index": i, "sweep_key": sweep_key, "value": value, "status": "ok",
            "Cl": "", "Cd": "", "lift": "", "drag": "", "mesh": "", "error": "",
        }
        try:
            case_toml_path, title = build_case_toml(base_data, sweep_key, value, runs_dir, i)
            mesh_path = run_mesher(case_toml_path, title)
            row["mesh"] = str(mesh_path)

            result = worker.solve(str(mesh_path), str(case_toml_path))
            row.update(result)

            print(f"[{i + 1}/{len(sweep_values)}] {sweep_key}={value} "
                  f"-> Cl={result['Cl']:.4g} Cd={result['Cd']:.4g}")
        except Exception as e:
            row["status"] = "error"
            row["error"] = str(e)
            print(f"[{i + 1}/{len(sweep_values)}] {sweep_key}={value} -> FAILED: {e}")

        rows.append(row)

        # The CFD worker is a single long-lived process for the whole sweep
        # (see 07_CFD_Scripts/Batch/cfd_worker.jl). If it has died, further
        # points can't be solved — record the rest as skipped and stop
        # instead of failing through them one at a time.
        if worker.proc.poll() is not None:
            remaining = sweep_values[i + 1:]
            for j, v in enumerate(remaining, start=i + 1):
                rows.append({
                    "index": j, "sweep_key": sweep_key, "value": v, "status": "skipped",
                    "Cl": "", "Cd": "", "lift": "", "drag": "", "mesh": "",
                    "error": "CFD worker process had already exited",
                })
            print(f"Julia CFD worker exited — skipping {len(remaining)} remaining point(s).")
            break
finally:
    worker.quit()

# ============================================================
# RESULTS
# ============================================================

csv_path = results_dir / "results.csv"
plot_path = results_dir / "results.png"

write_results_csv(rows, csv_path)
print(f"\nWrote {csv_path}")

n_ok = sum(1 for r in rows if r["status"] == "ok")
if n_ok > 0:
    plot_results(rows, sweep_key, plot_path)
    print(f"Wrote {plot_path}")
else:
    print("No successful points — skipped plot.")

n_bad = len(rows) - n_ok
if n_bad:
    print(f"{n_bad}/{len(rows)} point(s) failed or were skipped — see the 'status'/'error' columns in {csv_path.name}")
