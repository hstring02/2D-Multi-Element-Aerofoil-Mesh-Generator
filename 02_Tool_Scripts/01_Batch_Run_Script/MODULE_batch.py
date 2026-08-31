import copy
import csv
import os
import subprocess
import sys
from pathlib import Path

import toml

# Forces UTF-8 for the meshing subprocess's stdout/stderr. When Python's
# stdout is piped (as it is under subprocess.run/capture_output) rather than
# attached to a real console, Windows falls back to the system codepage
# (e.g. cp1252) instead of UTF-8, which crashes on any non-ASCII character a
# print() emits (this repo's mesher does print one, an arrow in a status
# line). Without this, that crash looks like a meshing failure.
SUBPROCESS_ENV = {**os.environ, "PYTHONIOENCODING": "utf-8"}

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
MESHING_DIR = REPO_ROOT / "02_Tool_Scripts" / "02_Meshing_Script"
MESHES_DIR = REPO_ROOT / "04_Meshes"
CFD_WORKER_SCRIPT = REPO_ROOT / "02_Tool_Scripts" / "03_CFD_Scripts" / "cfd_worker.jl"


# ============================================================
# CASE TOML GENERATION
# ============================================================

def set_by_path(data, path, value):
    """
    Sets a value in a nested dict/list structure from a dotted path, e.g.
    "flow.VELOCITY" or "foils.AOA.0" (numeric segments index into lists).
    """
    segments = path.split(".")
    node = data
    for seg in segments[:-1]:
        node = node[int(seg)] if seg.isdigit() else node[seg]

    last = segments[-1]
    key = int(last) if last.isdigit() else last
    node[key] = value


def slugify_value(value):
    return str(value).replace(".", "p").replace("-", "neg")


def build_case_toml(base_data, sweep_key, value, out_dir, index):
    """
    Deep-copies the base case, applies the sweep override, disables the
    interactive Gmsh preview (batch runs must not block on a GUI window),
    and writes the derived case out under a unique title. Returns
    (case_toml_path, title).
    """
    data = copy.deepcopy(base_data)
    set_by_path(data, sweep_key, value)

    data.setdefault("visualisation", {})["ENABLED"] = False

    base_title = base_data.get("title", "case")
    key_slug = sweep_key.replace(".", "_")
    title = f"{base_title}__{key_slug}_{slugify_value(value)}__{index:03d}"
    data["title"] = title

    out_dir.mkdir(parents=True, exist_ok=True)
    case_toml_path = out_dir / f"{title}.toml"
    with open(case_toml_path, "w") as f:
        toml.dump(data, f)

    return case_toml_path, title


def build_case_toml_multi(base_data, overrides, out_dir, index):
    """
    Like build_case_toml, but applies any number of dotted-path overrides at
    once instead of a single (sweep_key, value) pair — for cartesian-grid
    batch sweeps over multiple keys. `overrides` is a dict of
    {dotted_path: value}, applied via set_by_path in insertion order.
    """
    data = copy.deepcopy(base_data)
    for key, value in overrides.items():
        set_by_path(data, key, value)

    data.setdefault("visualisation", {})["ENABLED"] = False

    base_title = base_data.get("title", "case")
    override_slug = "_".join(
        f"{key.replace('.', '_')}_{slugify_value(value)}"
        for key, value in overrides.items()
    )
    title = f"{base_title}__{override_slug}__{index:03d}"
    data["title"] = title

    out_dir.mkdir(parents=True, exist_ok=True)
    case_toml_path = out_dir / f"{title}.toml"
    with open(case_toml_path, "w") as f:
        toml.dump(data, f)

    return case_toml_path, title


# ============================================================
# MESHING
# ============================================================

def run_mesher(case_toml_path, title, dest_subdir="Batch"):
    """
    Runs 02_Tool_Scripts/02_Meshing_Script/RUN_SCRIPT.py on the given case
    file, then moves the resulting .unv into 04_Meshes/<dest_subdir>
    (keeping derived-run meshes out of the top-level 04_Meshes folder) and
    returns its new path.
    Raises RuntimeError (with the mesher's stderr) on failure.
    """
    result = subprocess.run(
        [sys.executable, "RUN_SCRIPT.py", str(case_toml_path)],
        cwd=MESHING_DIR,
        capture_output=True,
        text=True,
        encoding="utf-8",
        env=SUBPROCESS_ENV,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Meshing failed for {case_toml_path.name}:\n{result.stderr.strip()}"
        )

    generated_path = MESHES_DIR / f"{title}.unv"
    if not generated_path.is_file():
        raise RuntimeError(
            f"Mesher reported success but no mesh was found at {generated_path}"
        )

    dest_dir = MESHES_DIR / dest_subdir
    dest_dir.mkdir(parents=True, exist_ok=True)
    mesh_path = dest_dir / f"{title}.unv"
    generated_path.replace(mesh_path)
    return mesh_path


# ============================================================
# PERSISTENT CFD WORKER
# ============================================================

class CFDWorker:
    """
    Wraps a single long-lived `julia cfd_worker.jl` process so XCALibre's
    JIT/precompile cost is paid once for the whole sweep instead of once
    per mesh. See 02_Tool_Scripts/03_CFD_Scripts/cfd_worker.jl for the line protocol.
    """

    def __init__(self, julia_exe="julia"):
        self.proc = subprocess.Popen(
            [julia_exe, str(CFD_WORKER_SCRIPT)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            bufsize=1,
            env=SUBPROCESS_ENV,
        )
        self._wait_ready()

    def _drain_output(self):
        try:
            return self.proc.stdout.read()
        except Exception:
            return ""

    def _check_alive(self, context):
        if self.proc.poll() is not None:
            raise RuntimeError(
                f"Julia CFD worker exited unexpectedly ({context}).\n"
                f"{self._drain_output()}"
            )

    def _wait_ready(self):
        while True:
            self._check_alive("while starting up")
            line = self.proc.stdout.readline()
            if line.strip() == "##READY##":
                return

    def solve(self, mesh_path, toml_path):
        self._check_alive(f"before solving {mesh_path}")
        self.proc.stdin.write(f"{mesh_path}|{toml_path}\n")
        self.proc.stdin.flush()

        while True:
            self._check_alive(f"while solving {mesh_path}")
            line = self.proc.stdout.readline()
            if not line:
                continue
            line = line.strip()

            if line.startswith("##RESULT## "):
                _, cl, cd, lift, drag = line[len("##RESULT## "):].split("|")
                return {"Cl": float(cl), "Cd": float(cd), "lift": float(lift), "drag": float(drag)}

            if line.startswith("##ERROR## "):
                _, message = line[len("##ERROR## "):].split("|", 1)
                raise RuntimeError(f"CFD solve failed for {mesh_path}: {message}")

            # Anything else (solver residual logging, warnings, etc.) is
            # just worker chatter — ignore it and keep reading.

    def quit(self):
        if self.proc.poll() is None:
            try:
                self.proc.stdin.write("##QUIT##\n")
                self.proc.stdin.flush()
                self.proc.stdin.close()
            except (BrokenPipeError, OSError):
                pass
            self.proc.wait(timeout=30)


# ============================================================
# RESULTS
# ============================================================

FIELDNAMES = ["index", "sweep_key", "value", "status", "Cl", "Cd", "lift", "drag", "mesh", "error"]


def write_results_csv(rows, out_path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in FIELDNAMES})


def write_results_csv_multi(rows, sweep_keys, out_path):
    """
    Like write_results_csv, but for a cartesian-grid sweep over several
    sweep_keys at once: one column per swept key instead of a single
    'sweep_key'/'value' pair.
    """
    fieldnames = ["index", *sweep_keys, "status", "Cl", "Cd", "lift", "drag", "mesh", "error"]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def plot_results(rows, sweep_key, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ok_rows = [r for r in rows if r["status"] == "ok"]
    ok_rows.sort(key=lambda r: r["value"])
    values = [r["value"] for r in ok_rows]
    cls = [r["Cl"] for r in ok_rows]
    cds = [r["Cd"] for r in ok_rows]

    fig, (ax_cl, ax_cd) = plt.subplots(1, 2, figsize=(10, 4))

    ax_cl.plot(values, cls, "o-")
    ax_cl.set_xlabel(sweep_key)
    ax_cl.set_ylabel("Cl")
    ax_cl.grid(True)

    ax_cd.plot(values, cds, "o-")
    ax_cd.set_xlabel(sweep_key)
    ax_cd.set_ylabel("Cd")
    ax_cd.grid(True)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
