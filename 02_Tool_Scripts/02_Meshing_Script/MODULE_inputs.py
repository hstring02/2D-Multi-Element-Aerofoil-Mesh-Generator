# ============================================================
# INPUT LOADING
# ============================================================

from pathlib import Path

import toml

from MODULE_log import success
from MODULE_mesh import compute_first_layer_height

# Sections a defaults file may supply, merged key-by-key (case file wins) so
# a case only needs to set the keys it wants to override, e.g. just
# ENABLED = false to turn a feature off.
MERGED_SECTIONS = ("mesh_settings", "global_refinement", "boundary_layer", "wake_refinement", "te_refinement")


def load_case_with_defaults(config_path, defaults_dir):
    """
    Loads a case TOML and, if its [mesh_settings].DEFAULTS names a defaults
    file, merges that file in per-section, per-key underneath it (the case
    file's own values always win). Returns the merged data in the same
    nested shape as the case file itself, so callers can keep unpacking it
    exactly as before regardless of how much was actually filled in by the
    case file versus inherited from the defaults.

    If boundary_layer.ENABLED, also resolves the target y+ into an actual
    first-layer height and stores it as boundary_layer.FIRST_LAYER_HEIGHT,
    since that conversion depends on [flow] and only needs doing once here
    rather than in every caller.
    """
    config_path = Path(config_path)
    data = toml.load(config_path)

    defaults_name = data.get("mesh_settings", {}).get("DEFAULTS")
    if defaults_name:
        defaults_arg = Path(defaults_name)
        defaults_path = defaults_arg if defaults_arg.parent != Path(".") else Path(defaults_dir) / defaults_arg
        if not defaults_path.is_file():
            raise FileNotFoundError(f"mesh_settings.DEFAULTS file not found: {defaults_path}")
        defaults_data = toml.load(defaults_path)

        for section in MERGED_SECTIONS:
            if section in defaults_data or section in data:
                data[section] = {**defaults_data.get(section, {}), **data.get(section, {})}

    bl_config = data.get("boundary_layer", {})
    if bl_config.get("ENABLED", False):
        flow_config = data.get("flow", {})
        missing_flow_keys = [
            key for key in ("DENSITY", "VELOCITY", "VISCOSITY", "REYNOLDS_LENGTH", "Y_PLUS")
            if key not in flow_config
        ]
        if missing_flow_keys:
            raise ValueError(
                f"boundary_layer is ENABLED but [flow] is missing: {', '.join(missing_flow_keys)}"
            )

        first_layer_height = compute_first_layer_height(
            density=flow_config["DENSITY"],
            velocity=flow_config["VELOCITY"],
            viscosity=flow_config["VISCOSITY"],
            ref_length=flow_config["REYNOLDS_LENGTH"],
            y_plus=flow_config["Y_PLUS"],
        )
        success(
            f"[boundary_layer] target y+={flow_config['Y_PLUS']:.6g} -> "
            f"computed FIRST_LAYER_HEIGHT={first_layer_height:.6g} m"
        )
        bl_config["FIRST_LAYER_HEIGHT"] = first_layer_height

    return data
