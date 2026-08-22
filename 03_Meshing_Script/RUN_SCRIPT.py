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
    get_airfoil_curves, add_distance_field, add_threshold_field,
    set_background_field, add_boundary_layer_field, set_boundary_layer_field,
    get_straight_line_curves, enforce_cells_on_curves, get_curve_endpoints,
    add_box_field, add_min_field, compute_growth_dist_max,
    classify_curves_by_points
)

# ============================================================
# LOAD DATA
# ============================================================

INPUT_FILE = "2_el_wing.toml"

SCRIPT_DIR = Path(__file__).resolve().parent.parent
CONFIG_PATH = Path(sys.argv[1]) if len(sys.argv) > 1 else SCRIPT_DIR / "02_Mesh_Input_File" / INPUT_FILE
FOILS_DIR = SCRIPT_DIR / "01_Foils"

data = toml.load(CONFIG_PATH)

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
BL_FIRST_LAYER_HEIGHT = BL_CONFIG.get("FIRST_LAYER_HEIGHT")
BL_GROWTH_RATE = BL_CONFIG.get("GROWTH_RATE")
BL_THICKNESS = BL_CONFIG.get("THICKNESS")

DIST_MIN = BL_THICKNESS if BL_ENABLED else 0.0

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

print(f"Fluid Surface Tag:       {fluid_surface_tag}")
print(f"Farfield Curves:         {farfield_boundary_curves}")
print(f"Airfoil Boundary Curves: {airfoil_boundary_curves}")

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
# The boolean cut above merges every element into one fluid surface and
# renumbers every curve, so which element originally produced a given
# curve isn't tracked directly any more. Recovered here by nearest-point
# comparison against each element's own (already known) point cloud —
# see classify_curves_by_points for why bounding boxes alone aren't
# reliable enough for closely-packed multi-element geometry.

element_indices = classify_curves_by_points(
    airfoil_boundary_curves,
    [info["points"] for info in element_info]
)

element_curves = [[] for _ in element_info]
for curve_tag, elem_idx in zip(airfoil_boundary_curves, element_indices):
    element_curves[elem_idx].append(curve_tag)

for info, curves in zip(element_info, element_curves):
    print(f"  {info['element']}: boundary curves {curves}")

# ============================================================
# TRAILING-EDGE CELL COUNT
# ============================================================

te_line_curves = get_straight_line_curves(airfoil_boundary_curves)
enforce_cells_on_curves(te_line_curves, TE_CELLS)

# Corner points where each blunt TE base line meets the upper/lower airfoil
# surface. These are convex ~90 degree corners, so the boundary layer's
# default miter join stretches a single quad column across them instead of
# filling the corner — this is what produces the high aspect ratio cells at
# the outer edge of the layer. Passing them as FanPointsList makes Gmsh fan
# several triangular columns across the corner instead. Also used below to
# anchor the TE-corner refinement field, since flow accelerates sharply
# around the same corners.
te_fan_points = get_curve_endpoints(te_line_curves)

# ============================================================
# MESHING CONTROL
# ============================================================

_next_field_id = 0

def new_field_id():
    global _next_field_id
    _next_field_id += 1
    return _next_field_id

# NEAR-SURFACE FIELD (ONE PER ELEMENT)
# ============================================================
# Each element grows its own near-surface mesh size from its own
# POINT_SIZE (rather than a single shared MESH_MIN floor) out to
# MESH_MAX at the shared GROWTH_RATE. A finely-pointed flap and a
# coarsely-pointed main element each get a DIST_MAX consistent with
# their own starting resolution, instead of both being forced through
# the same growth region regardless of how fine their own surface mesh
# actually is.
combined_field_ids = []

for info, curves in zip(element_info, element_curves):
    point_size = info["point_size"]
    dist_max = compute_growth_dist_max(point_size, MESH_MAX, GLOBAL_GROWTH_RATE, DIST_MIN)
    print(
        f"[{info['element']}] POINT_SIZE={point_size} GROWTH_RATE={GLOBAL_GROWTH_RATE} -> "
        f"computed DIST_MAX={dist_max:.6g} m"
    )

    dist_field_id = new_field_id()
    add_distance_field(field_id=dist_field_id, curves=curves)

    size_field_id = new_field_id()
    add_threshold_field(
        field_id=size_field_id,
        in_field=dist_field_id,
        size_min=point_size,
        size_max=MESH_MAX,
        dist_min=DIST_MIN,
        dist_max=dist_max
    )
    combined_field_ids.append(size_field_id)

# ============================================================
# WAKE REFINEMENT (BOX FIELD, ONE PER ELEMENT)
# ============================================================
# Each element gets its own wake box, anchored on that element's own TE
# point in x (world coordinates, so it already accounts for AOA/position)
# and sized as a multiple of that element's own chord. The box stays
# aligned with the global X axis (not the local chord line) since the
# wake convects downstream with the freestream, not along whichever way
# an individual flap element happens to be pitched.
#
# Box height spans the element's own Y-extent (its footprint projected
# onto the y-axis) plus a chord-scaled pad on each side, rather than a
# fixed height about the TE point — so an element pitched to a high AOA,
# whose shed wake spans a taller Y range, automatically gets a taller box.
#
# Box downstream length also scales with AOA: 1 chord at AOA = 0 degrees,
# ramping linearly up to 3 chords at AOA = 90 degrees (a more sharply
# pitched element sheds a larger, slower-resolving wake structure that
# needs a longer fine region to capture). AOA beyond 90 degrees clamps at
# the 3-chord ramp ceiling rather than continuing to grow.
WAKE_X_LENGTH_CHORDS_AOA_0  = 1.0
WAKE_X_LENGTH_CHORDS_AOA_90 = 3.0

if WAKE_ENABLED:
    # Same growth-rate logic as each element's own near-surface field:
    # the ramp-back distance outside the box is however far a geometric
    # series starting at the box's own cell size needs to grow, at the
    # shared GROWTH_RATE, to reach MESH_MAX. One shared value, since
    # WAKE_MESH_SIZE/MESH_MAX/GROWTH_RATE are all global (not per-element).
    wake_transition = compute_growth_dist_max(WAKE_MESH_SIZE, MESH_MAX, GLOBAL_GROWTH_RATE)
    print(
        f"[wake_refinement] MESH_SIZE={WAKE_MESH_SIZE} GROWTH_RATE={GLOBAL_GROWTH_RATE} -> "
        f"computed transition={wake_transition:.6g} m"
    )

    for info in element_info:
        chord = info["chord"]
        te_x, _ = info["te_point"]

        aoa_frac = min(abs(info["aoa"]), 90.0) / 90.0
        x_length_chords = WAKE_X_LENGTH_CHORDS_AOA_0 + aoa_frac * (
            WAKE_X_LENGTH_CHORDS_AOA_90 - WAKE_X_LENGTH_CHORDS_AOA_0
        )

        x_min = te_x + WAKE_X_START_CHORDS * chord
        x_max = x_min + x_length_chords * chord
        y_pad = WAKE_Y_HALF_HEIGHT_CHORDS * chord

        field_id = new_field_id()
        add_box_field(
            field_id=field_id,
            size_in=WAKE_MESH_SIZE,
            size_out=MESH_MAX,
            x_min=x_min,
            x_max=x_max,
            y_min=info["y_min"] - y_pad,
            y_max=info["y_max"] + y_pad,
            transition=wake_transition
        )
        combined_field_ids.append(field_id)

# ============================================================
# TRAILING-EDGE CORNER REFINEMENT
# ============================================================
# Flow accelerates sharply around each TE corner, so cell size needs to
# taper in well below each element's own near-surface size (POINT_SIZE)
# as the corner is approached, not just hold constant like the rest of
# the surface. Anchored on the TE corner points (te_fan_points, only
# populated for blunt TEs — a sharp TE has no separate corner to single
# out) so only that region gets the extra refinement.

if TE_REFINEMENT_ENABLED and te_fan_points:
    te_dist_field_id = new_field_id()
    add_distance_field(
        field_id=te_dist_field_id,
        points=te_fan_points
    )
    # size_max is MESH_MAX (not any element's POINT_SIZE) deliberately:
    # this field is combined with the others via Min, so its "far" value
    # only needs to stop competing with them once past dist_max — it
    # must never floor the WHOLE domain (including the farfield
    # boundary) at a size only meant to apply right at the airfoil
    # surface.
    te_size_field_id = new_field_id()
    add_threshold_field(
        field_id=te_size_field_id,
        in_field=te_dist_field_id,
        size_min=TE_REFINEMENT_MESH_SIZE,
        size_max=MESH_MAX,
        dist_min=0.0,
        dist_max=TE_REFINEMENT_DIST_MAX
    )
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

gmsh.option.setNumber("Mesh.Algorithm", 6)           # Frontal-Delaunay
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

print(f"\nMesh successfully generated and saved to: {SCRIPT_DIR / '04_Meshes' / output_filename}")
print("\nBoundary summary:")
print(f"  airfoil  : curves {airfoil_boundary_curves}")
print(f"  inlet    : curve  {inlet_curve}")
print(f"  outlet   : curve  {outlet_curve}")
print(f"  farfield : curves {[top_curve, bottom_curve]}")

# ============================================================
# VISUALISATION
# ============================================================

gmsh.fltk.run()
gmsh.finalize()