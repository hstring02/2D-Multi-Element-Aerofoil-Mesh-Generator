# ============================================================
# PLOT A BATCH RUN'S results.csv
# ============================================================
#
# Usage:
#   python plot_results.py [path/to/results.csv]
#
# Defaults to batch_config/results.csv (next to this script) if no path is
# given. Saves a PNG next to the CSV and opens an interactive window.

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CSV = SCRIPT_DIR / "batch_config" / "results.csv"


def load_rows(csv_path):
    with open(csv_path, newline="") as f:
        return list(csv.DictReader(f))


def main():
    csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_CSV
    if not csv_path.is_file():
        raise FileNotFoundError(f"results CSV not found: {csv_path}")

    rows = load_rows(csv_path)
    if not rows:
        raise ValueError(f"{csv_path} has no data rows")

    sweep_key = rows[0]["sweep_key"]

    ok_rows = [r for r in rows if r["status"] == "ok"]
    ok_rows.sort(key=lambda r: float(r["value"]))
    values = [float(r["value"]) for r in ok_rows]
    cls = [float(r["Cl"]) for r in ok_rows]
    cds = [float(r["Cd"]) for r in ok_rows]

    fig, (ax_cl, ax_cd) = plt.subplots(1, 2, figsize=(10, 4.5))
    fig.suptitle(f"Batch sweep: {sweep_key}")

    ax_cl.plot(values, cls, "o-")
    ax_cl.set_xlabel(sweep_key)
    ax_cl.set_ylabel("Cl")
    ax_cl.grid(True, alpha=0.3)

    ax_cd.plot(values, cds, "o-")
    ax_cd.set_xlabel(sweep_key)
    ax_cd.set_ylabel("Cd")
    ax_cd.grid(True, alpha=0.3)

    fig.tight_layout()

    out_path = csv_path.with_suffix(".png")
    fig.savefig(out_path, dpi=150)
    print(f"Wrote {out_path}")

    bad_rows = [r for r in rows if r["status"] != "ok"]
    if bad_rows:
        print(f"{len(bad_rows)} point(s) not plotted (status != ok):")
        for r in bad_rows:
            print(f"  value={r['value']}: {r['status']} - {r['error']}")

    plt.show()


if __name__ == "__main__":
    main()
