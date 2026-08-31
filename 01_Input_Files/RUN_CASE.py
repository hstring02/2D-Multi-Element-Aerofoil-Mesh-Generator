# ============================================================
# CASE MANAGER
# ============================================================
#
# Single entry point for a unified case TOML (see
# INWORK_new_input_template.toml): reads [case].TYPE and routes the case to
# meshing only, a single mesh+solve, or a batch grid sweep. The underlying
# 02_Tool_Scripts/02_Meshing_Script/RUN_SCRIPT.py and
# 02_Tool_Scripts/01_Batch_Run_Script/RUN_SCRIPT.py
# entry points are unchanged and can still be run directly.

import sys
from pathlib import Path

from MODULE_case import run

INPUT_DIR = Path(__file__).resolve().parent

if len(sys.argv) > 1:
    arg_path = Path(sys.argv[1])
    CASE_PATH = arg_path if arg_path.parent != Path(".") else INPUT_DIR / arg_path
else:
    raise FileNotFoundError(
        "No case file given. Usage: python RUN_CASE.py <path-to-case-toml>"
    )

if not CASE_PATH.is_file():
    raise FileNotFoundError(f"Case file not found: {CASE_PATH}")

run(CASE_PATH)
