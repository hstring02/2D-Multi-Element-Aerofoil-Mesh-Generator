# ============================================================
# LOAD FUNCTIONS
# ============================================================

import re
import sys
from pathlib import Path
import toml
import gmsh

# CUSTOM FUNCTIONS
from MODULE_airfoil import read_transform_airfoil, build_airfoil, airfoils
from MODULE_geometry import build_farfield, subtract_airfoil
from MODULE_output import output_unv_xcalibre
from MODULE_mesh import (
    get_airfoil_curves, set_background_field, add_boundary_layer_field,
    set_boundary_layer_field, get_straight_line_curves, enforce_cells_on_curves,
    get_curve_endpoints, add_min_field, compute_first_layer_height,
    classify_curves_by_points, clamp_boundary_layer_thickness,
    add_near_surface_fields, add_wake_refinement_fields,
    add_te_corner_refinement_field
)
from MODULE_log import info, success, subtitle

# ============================================================
# LOAD DATA
# ============================================================

INPUT_FILE = "2_el_wing.toml"    # Defaults to this if no command-line argument is provided. Can be overridden by passing a path to a different TOML file as the first argument to this script.

SCRIPT_DIR = Path(__file__).resolve().parent.parent.parent
MESH_INPUT_DIR = SCRIPT_DIR / "01_Input_Files"

if len(sys.argv) > 1:
    arg_path = Path(sys.argv[1])
    CONFIG_PATH = arg_path if arg_path.parent != Path(".") else MESH_INPUT_DIR / arg_path
else:
    CONFIG_PATH = MESH_INPUT_DIR / INPUT_FILE

if not CONFIG_PATH.is_file():
    raise FileNotFoundError(f"Mesh input file not found: {CONFIG_PATH}")

data = toml.load(CONFIG_PATH)

FOILS_DIR = SCRIPT_DIR / "03_Foils"

XMIN = data["farfield"]["XMIN"]
XMAX = data["farfield"]["XMAX"]
YMIN = data["farfield"]["YMIN"]
YMAX = data["farfield"]["YMAX"]
MESH_MAX = data["global_refinement"]["MESH_MAX"]
GLOBAL_GROWTH_RATE = data["global_refinement"]["GROWTH_RATE"]
TE_CELLS = data["global_refinement"].get("TE_CELLS", 3)
MAX_CURVATURE_ANGLE = data["global_refinement"].get("MAX_CURVATURE_ANGLE", 360.0)
TE_FAN_ELEMENTS = data["global_refinement"].get("TE_FAN_ELEMENTS", 5)
BL_CONFIG = data.get("boundary_layer", {})
BL_ENABLED = BL_CONFIG.get("ENABLED", False)
BL_GROWTH_RATE = BL_CONFIG.get("GROWTH_RATE")
BL_THICKNESS = BL_CONFIG.get("THICKNESS")

DIST_MIN = BL_THICKNESS if BL_ENABLED else 0.0

if BL_ENABLED:
    FLOW_CONFIG = data.get("flow", {})
    missing_flow_keys = [
        key for key in ("DENSITY", "VELOCITY", "VISCOSITY", "REYNOLDS_LENGTH", "Y_PLUS")
        if key not in FLOW_CONFIG
    ]
    if missing_flow_keys:
        raise ValueError(
            f"boundary_layer is ENABLED but [flow] is missing: {', '.join(missing_flow_keys)}"
        )

    BL_FIRST_LAYER_HEIGHT = compute_first_layer_height(
        density=FLOW_CONFIG["DENSITY"],
        velocity=FLOW_CONFIG["VELOCITY"],
        viscosity=FLOW_CONFIG["VISCOSITY"],
        ref_length=FLOW_CONFIG["REYNOLDS_LENGTH"],
        y_plus=FLOW_CONFIG["Y_PLUS"],
    )
    success(
        f"[boundary_layer] target y+={FLOW_CONFIG['Y_PLUS']:.6g} -> "
        f"computed FIRST_LAYER_HEIGHT={BL_FIRST_LAYER_HEIGHT:.6g} m"
    )
else:
    BL_FIRST_LAYER_HEIGHT = None

WAKE_CONFIG = data.get("wake_refinement", {})
WAKE_ENABLED = WAKE_CONFIG.get("ENABLED", False)
WAKE_MESH_SIZE = WAKE_CONFIG.get("MESH_SIZE")
WAKE_X_START_CHORDS = WAKE_CONFIG.get("X_START_CHORDS", 0.1)
WAKE_Y_HALF_HEIGHT_CHORDS = WAKE_CONFIG.get("Y_HALF_HEIGHT_CHORDS", 0.5)

TE_CONFIG = data.get("te_refinement", {})
TE_REFINEMENT_ENABLED = TE_CONFIG.get("ENABLED", False)
TE_REFINEMENT_MESH_SIZE = TE_CONFIG.get("MESH_SIZE")
TE_REFINEMENT_DIST_MAX = TE_CONFIG.get("DIST_MAX")

# ============================================================
# BUILD AIRFOIL
# ============================================================

gmsh.initialize()
gmsh.model.add("FSAE_Airfoil")
occ = gmsh.model.occ

airfoil_surfaces, element_info = airfoils(data, occ, FOILS_DIR)    # Reads airfoil data points, builds curves and returns surfaces + per-element chord/TE info

if BL_ENABLED:
    BL_THICKNESS = clamp_boundary_layer_thickness(element_info, BL_THICKNESS, BL_FIRST_LAYER_HEIGHT)
    DIST_MIN = BL_THICKNESS

# ============================================================
# BUILD FARFIELD
# ============================================================

farfield_surface, farfield_lines = build_farfield(
    occ,
    XMIN,
    XMAX,
    YMIN,
    YMAX,
    MESH_MAX  
)

farfield_corners = [        # Snap farfield corner coordinates so we can re-identify them post-cut
    (XMIN, YMIN),
    (XMAX, YMIN),
    (XMAX, YMAX),
    (XMIN, YMAX),
]

for index, surfaces in enumerate(airfoil_surfaces):
    
    airfoil_surface = surfaces

    fluid_dim_tag = subtract_airfoil(
        occ,
        farfield_surface,
        airfoil_surface
    )

occ.synchronize()


fluid_surface_tag = fluid_dim_tag[1]    # Extract the surface tag integer from the tuple
all_fluid_curves = get_airfoil_curves(fluid_surface_tag)

CORNER_TOL = 1e-6

def point_is_farfield_corner(x, y):
    for cx, cy in farfield_corners:
        if abs(x - cx) < CORNER_TOL and abs(y - cy) < CORNER_TOL:
            return True
    return False

farfield_boundary_curves = []
airfoil_boundary_curves  = []

for curve_tag in all_fluid_curves:
    # Get the two end-points of this curve
    bnd = gmsh.model.getBoundary([(1, curve_tag)], oriented=False)
    end_pts = [tag for _, tag in bnd]

    # Check if BOTH endpoints sit on farfield corners
    pts_on_ff = 0
    for pt_tag in end_pts:
        coords = gmsh.model.getValue(0, pt_tag, [])  # returns [x, y, z]
        if point_is_farfield_corner(coords[0], coords[1]):
            pts_on_ff += 1

    if pts_on_ff == len(end_pts):
        farfield_boundary_curves.append(curve_tag)
    else:
        airfoil_boundary_curves.append(curve_tag)

if not airfoil_boundary_curves:
    raise RuntimeError(
        "No airfoil boundary curves recovered after boolean cut. "
        "Check that the airfoil lies fully inside the farfield box."
    )

if len(farfield_boundary_curves) != 4:
    raise RuntimeError(
        f"Expected 4 farfield curves but found {len(farfield_boundary_curves)}. "
        "The farfield rectangle may have been fragmented by the boolean cut."
    )


# IDENTIFY INDIVIDUAL FARFIELD WALLS
def curve_midpoint(curve_tag):
    bbox = gmsh.model.getBoundingBox(1, curve_tag)
    # bbox = [xmin, ymin, zmin, xmax, ymax, zmax]
    mx = (bbox[0] + bbox[3]) / 2.0
    my = (bbox[1] + bbox[4]) / 2.0
    return mx, my

inlet_curve   = None
outlet_curve  = None
top_curve     = None
bottom_curve  = None

for c in farfield_boundary_curves:
    mx, my = curve_midpoint(c)
    if   abs(mx - XMIN) < CORNER_TOL * 1000:
        inlet_curve  = c
    elif abs(mx - XMAX) < CORNER_TOL * 1000:
        outlet_curve = c
    elif abs(my - YMAX) < CORNER_TOL * 1000:
        top_curve    = c
    elif abs(my - YMIN) < CORNER_TOL * 1000:
        bottom_curve = c

for name, val in [("inlet", inlet_curve), ("outlet", outlet_curve),
                    ("top", top_curve), ("bottom", bottom_curve)]:
    if val is None:
        raise RuntimeError(
            f"Could not identify the '{name}' farfield wall. "
            "Check XMIN/XMAX/YMIN/YMAX in config.py match the geometry exactly."
        )

# ============================================================
# CLASSIFY AIRFOIL CURVES BY ELEMENT
# ============================================================

element_indices = classify_curves_by_points(
    airfoil_boundary_curves,
    [info["points"] for info in element_info]
)

element_curves = [[] for _ in element_info]
for curve_tag, elem_idx in zip(airfoil_boundary_curves, element_indices):
    element_curves[elem_idx].append(curve_tag)

# ============================================================
# TRAILING-EDGE CELL COUNT
# ============================================================

te_line_curves = get_straight_line_curves(airfoil_boundary_curves)
enforce_cells_on_curves(te_line_curves, TE_CELLS)
te_fan_points = get_curve_endpoints(te_line_curves)

# ============================================================
# MESHING CONTROL
# ============================================================

_next_field_id = 0

def new_field_id():
    global _next_field_id
    _next_field_id += 1
    return _next_field_id

combined_field_ids = add_near_surface_fields(
    element_info, element_curves, MESH_MAX, GLOBAL_GROWTH_RATE, DIST_MIN, new_field_id
)

if WAKE_ENABLED:
    combined_field_ids += add_wake_refinement_fields(
        element_info, WAKE_MESH_SIZE, MESH_MAX, GLOBAL_GROWTH_RATE,
        WAKE_X_START_CHORDS, WAKE_Y_HALF_HEIGHT_CHORDS, new_field_id
    )

if TE_REFINEMENT_ENABLED:
    te_size_field_id = add_te_corner_refinement_field(
        te_fan_points, TE_REFINEMENT_MESH_SIZE, TE_REFINEMENT_DIST_MAX, MESH_MAX, new_field_id
    )
    if te_size_field_id is not None:
        combined_field_ids.append(te_size_field_id)

# ============================================================
# COMBINE SIZING FIELDS
# ============================================================

if len(combined_field_ids) > 1:
    min_field_id = new_field_id()
    add_min_field(field_id=min_field_id, in_fields=combined_field_ids)
    background_field_id = min_field_id
else:
    background_field_id = combined_field_ids[0]

set_background_field(background_field_id)

# INFLATION LAYER (BOUNDARY LAYER FIELD)
if BL_ENABLED:
    if te_fan_points:
        gmsh.option.setNumber("Mesh.BoundaryLayerFanElements", TE_FAN_ELEMENTS)
    bl_field_id = new_field_id()
    add_boundary_layer_field(
        field_id=bl_field_id,
        curves=airfoil_boundary_curves,
        first_layer_height=BL_FIRST_LAYER_HEIGHT,
        growth_rate=BL_GROWTH_RATE,
        thickness=BL_THICKNESS,
        fan_points=te_fan_points
    )
    set_boundary_layer_field(bl_field_id)

# ============================================================
# MESH STABILISATION OPTIONS
# ============================================================

gmsh.option.setNumber("Mesh.Algorithm", 5)           # Delaunay (Frontal-Delaunay/6 is more prone to intermittent "Edge not recovered" failures on tight multi-element geometry)
gmsh.option.setNumber("Mesh.Optimize", 1)
gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)    # Turning this one (to 1) would flood the convace part of the wing with cells. Not really needed.
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 360.0 / MAX_CURVATURE_ANGLE)

# ============================================================
# PHYSICAL GROUPS (Required for XCALibre UNV boundary groups)
# ============================================================

# 1D boundaries only — XCALibre cannot label 2D regions, only 1D boundaries
gmsh.model.addPhysicalGroup(1, airfoil_boundary_curves, name="airfoil")
gmsh.model.addPhysicalGroup(1, [inlet_curve],           name="inlet")
gmsh.model.addPhysicalGroup(1, [outlet_curve],          name="outlet")
gmsh.model.addPhysicalGroup(1, [top_curve, bottom_curve], name="farfield")

# ============================================================
# GENERATE MESH
# ============================================================

gmsh.model.mesh.generate(1)     # Generate 1D mesh for the boundary curves
gmsh.model.mesh.generate(2)     # Generate 2D mesh for the fluid domain

# ============================================================
# EXPORT OPTIONS
# ============================================================

output_filename = data["title"]
if not output_filename.lower().endswith(".unv"):
    output_filename += ".unv"
output_unv_xcalibre(output_filename)

node_count = len(gmsh.model.mesh.getNodes()[0])
element_count = sum(len(tags) for tags in gmsh.model.mesh.getElements()[1])

success(f"Mesh generated and saved to: {SCRIPT_DIR / '04_Meshes' / output_filename}")

subtitle("Summary")
info(f"  nodes    : {node_count}")
info(f"  elements : {element_count}")
info(f"  airfoil  : curves {airfoil_boundary_curves}")
info(f"  inlet    : curve  {inlet_curve}")
info(f"  outlet   : curve  {outlet_curve}")
info(f"  farfield : curves {[top_curve, bottom_curve]}")

# ============================================================
# VISUALISATION
# ============================================================

if data.get("visualisation", {}).get("ENABLED", False):
    gmsh.fltk.run()
    gmsh.finalize()