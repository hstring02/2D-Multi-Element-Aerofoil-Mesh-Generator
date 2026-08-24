import math

import gmsh


# ============================================================
# GEOMETRIC GROWTH -> DIST_MAX
# ============================================================

def compute_growth_dist_max(mesh_min, mesh_max, growth_rate, dist_min=0.0):
    """
    Converts a target geometric growth rate into the DistMax a Threshold
    field needs so its background cell size grows from `mesh_min` to
    `mesh_max` at that rate away from the surface.

    A Threshold field only interpolates cell size LINEARLY between
    DistMin and DistMax — it has no direct growth-rate control. This
    instead sizes DistMax to match the distance a geometric series of
    rings would cover, where each ring's thickness equals its own
    (growing) cell size: mesh_min, mesh_min*r, mesh_min*r^2, ... — the
    same growth-series logic used for the boundary layer field, just
    applied to coarser background rings instead of wall-normal layers.

    Parameters
    ----------
    mesh_min    : float – cell size at dist_min (MESH_MIN)
    mesh_max    : float – target cell size to grow out to (MESH_MAX)
    growth_rate : float – ring-to-ring growth ratio (> 1.0), e.g. 1.2
    dist_min    : float – distance at which growth starts (DIST_MIN)
    """
    if growth_rate <= 1.0:
        raise ValueError(f"GROWTH_RATE must be > 1.0 (got {growth_rate})")
    if mesh_max <= mesh_min:
        raise ValueError(
            f"MESH_MAX ({mesh_max}) must be greater than MESH_MIN ({mesh_min}) "
            "to compute a growth-rate-based DIST_MAX"
        )

    n_rings = math.ceil(math.log(mesh_max / mesh_min) / math.log(growth_rate))
    growth_distance = mesh_min * (growth_rate ** n_rings - 1) / (growth_rate - 1)

    return dist_min + growth_distance


# ============================================================
# FLOW CONDITIONS -> FIRST LAYER HEIGHT
# ============================================================

def compute_first_layer_height(density, velocity, viscosity, ref_length, y_plus):
    """
    Computes the wall-normal height of the boundary layer's first cell
    needed to hit a target y+ at the given freestream flow conditions.

    Wall shear stress is estimated from the Schlichting flat-plate
    turbulent skin-friction correlation (valid for Re < 10^9):

        Cf = (2 * log10(Re_L) - 0.65)^-2.3

    which gives the friction velocity u_tau = sqrt(tau_wall / rho) and,
    from the definition of y+, the required first-layer height:

        y = y_plus * viscosity / (density * u_tau)

    Parameters
    ----------
    density    : float – freestream fluid density, kg/m^3
    velocity   : float – freestream velocity, m/s
    viscosity  : float – dynamic viscosity, Pa.s
    ref_length : float – reference length for the Reynolds number
                 (e.g. the chord of the largest element), m
    y_plus     : float – target non-dimensional wall distance
    """
    reynolds = density * velocity * ref_length / viscosity
    if reynolds <= 1.0:
        raise ValueError(
            f"Reynolds number must be > 1 to estimate skin friction (got {reynolds:.6g})"
        )

    skin_friction = (2.0 * math.log10(reynolds) - 0.65) ** -2.3
    wall_shear_stress = 0.5 * skin_friction * density * velocity ** 2
    friction_velocity = math.sqrt(wall_shear_stress / density)

    return y_plus * viscosity / (density * friction_velocity)


# ============================================================
# CLASSIFY AIRFOIL CURVES BY ELEMENT
# ============================================================

def classify_curves_by_points(curves, element_points_list, samples=5):
    """
    Assigns each curve in `curves` (1D curve tags) to the index of the
    element whose own (transformed) airfoil points it lies closest to.

    The boolean cut that subtracts every airfoil from the farfield
    merges them into one fluid surface and renumbers every curve, so
    there's no direct record of which element originally produced a
    given curve. Each curve is sampled at a few points along its own
    parametrisation; the element whose stored point cloud is nearest
    those samples, on average, is the one that curve actually belongs
    to. Comparing against the real point cloud (rather than each
    element's bounding box) stays correct even when elements are packed
    closely enough for their bounding boxes to overlap.

    Parameters
    ----------
    curves              : list of int         – curve tags to classify
    element_points_list : list of list[(x,y)] – each element's own
                                                 transformed airfoil points,
                                                 in the same order the
                                                 returned indices refer to
    samples             : int                 – parametric samples per curve
    """
    assignments = []

    for curve_tag in curves:
        bounds = gmsh.model.getParametrizationBounds(1, curve_tag)
        t_min, t_max = bounds[0][0], bounds[1][0]

        sample_points = []
        for i in range(samples):
            t = t_min + (t_max - t_min) * i / (samples - 1) if samples > 1 else t_min
            x, y, _ = gmsh.model.getValue(1, curve_tag, [t])
            sample_points.append((x, y))

        best_index = None
        best_dist = None
        for elem_idx, points in enumerate(element_points_list):
            total = 0.0
            for sx, sy in sample_points:
                total += min((sx - px) ** 2 + (sy - py) ** 2 for px, py in points)
            if best_dist is None or total < best_dist:
                best_dist = total
                best_index = elem_idx

        assignments.append(best_index)

    return assignments


# ============================================================
# EXTRACT AIRFOIL BOUNDARY
# ============================================================

def get_airfoil_curves(fluid_surface):
    """
    Returns all 1D boundary curves of the given fluid surface tag (int).
    """
    model = gmsh.model

    # Normalise input to the list-of-dimtags format getBoundary expects
    if isinstance(fluid_surface, int):
        dim_tags = [(2, fluid_surface)]
    elif isinstance(fluid_surface, tuple):
        dim_tags = [fluid_surface]
    else:
        dim_tags = list(fluid_surface)

    boundary = model.getBoundary(
        dim_tags,
        oriented=False,
        recursive=False
    )

    return [tag for dim, tag in boundary if dim == 1]


# ============================================================
# DISTANCE FIELD
# ============================================================

def add_distance_field(field_id, curves=None, points=None, sampling=100):
    """
    Creates a Distance field from the given curve and/or point tags.

    Parameters
    ----------
    field_id  : int         – Gmsh field ID to assign
    curves    : list|None   – curve tags to measure distance from
    points    : list|None   – point tags to measure distance from (e.g.
                              trailing-edge corners, for a refinement
                              region anchored on a single point rather
                              than a whole curve)
    sampling  : int         – number of sample points per curve
                              (required in Gmsh >= 4.10; ignored in older
                              versions)
    """
    field = gmsh.model.mesh.field

    field.add("Distance", field_id)
    if curves:
        field.setNumbers(field_id, "CurvesList", curves)
    if points:
        field.setNumbers(field_id, "PointsList", points)

    # Sampling is mandatory in Gmsh >= 4.10; safe to set in all versions.
    field.setNumber(field_id, "Sampling", sampling)


# ============================================================
# THRESHOLD FIELD
# ============================================================

def add_threshold_field(
    field_id,
    in_field,
    size_min,
    size_max,
    dist_min,
    dist_max
):
    """
    Converts a Distance field into a mesh-size field via a linear ramp.

    size_min is used within dist_min of the curves;
    size_max is used beyond dist_max; linear interpolation in between.
    """
    field = gmsh.model.mesh.field

    field.add("Threshold", field_id)
    field.setNumber(field_id, "InField",  in_field)
    field.setNumber(field_id, "SizeMin",  size_min)
    field.setNumber(field_id, "SizeMax",  size_max)
    field.setNumber(field_id, "DistMin",  dist_min)
    field.setNumber(field_id, "DistMax",  dist_max)


# ============================================================
# BOX FIELD (WAKE / REGIONAL REFINEMENT)
# ============================================================

def add_box_field(
    field_id,
    size_in,
    size_out,
    x_min,
    x_max,
    y_min,
    y_max,
    transition=0.0
):
    """
    Creates a Box field: mesh size is `size_in` inside the given
    rectangle and relaxes out to `size_out` beyond it. Used for the wake
    refinement region behind the airfoils, but general-purpose for any
    rectangular region.

    Parameters
    ----------
    field_id    : int    – Gmsh field ID to assign
    size_in     : float  – mesh size inside the box
    size_out    : float  – mesh size far outside the box
    x_min/x_max : float  – box extents in x
    y_min/y_max : float  – box extents in y
    transition  : float  – distance outside the box over which the size
                            ramps from size_in to size_out (0 = sharp
                            cutoff at the box edge)
    """
    field = gmsh.model.mesh.field

    field.add("Box", field_id)
    field.setNumber(field_id, "VIn",  size_in)
    field.setNumber(field_id, "VOut", size_out)
    field.setNumber(field_id, "XMin", x_min)
    field.setNumber(field_id, "XMax", x_max)
    field.setNumber(field_id, "YMin", y_min)
    field.setNumber(field_id, "YMax", y_max)
    field.setNumber(field_id, "Thickness", transition)


# ============================================================
# MIN FIELD (COMBINE MULTIPLE SIZING FIELDS)
# ============================================================

def add_min_field(field_id, in_fields):
    """
    Creates a Min field that takes the smallest mesh size returned by
    `in_fields` at every point, so several sizing fields (e.g. the
    airfoil distance field and a wake box field) can be combined into a
    single background field.
    """
    field = gmsh.model.mesh.field

    field.add("Min", field_id)
    field.setNumbers(field_id, "FieldsList", in_fields)


# ============================================================
# BOUNDARY LAYER (INFLATION) FIELD
# ============================================================

def add_boundary_layer_field(
    field_id,
    curves,
    first_layer_height,
    growth_rate,
    thickness,
    fan_points=None,
    quads=True
):
    """
    Creates a structured inflation (boundary) layer field that grows
    quad cells off the given curves, from `first_layer_height` at the
    wall out to `thickness`, expanding by `growth_rate` each layer.

    Parameters
    ----------
    field_id            : int        – Gmsh field ID to assign
    curves               : list       – curve tags to grow the layer from
    first_layer_height   : float      – wall-normal height of the first cell
    growth_rate          : float      – layer-to-layer growth ratio (> 1.0)
    thickness             : float      – total boundary layer thickness
    fan_points           : list|None  – point tags (e.g. sharp trailing
                                        edges) where the layer should fan
                                        out instead of pinching
    quads                : bool       – generate quad elements in the layer
    """
    field = gmsh.model.mesh.field

    field.add("BoundaryLayer", field_id)
    field.setNumbers(field_id, "EdgesList", curves)
    field.setNumber(field_id, "hwall_n", first_layer_height)
    field.setNumber(field_id, "ratio", growth_rate)
    field.setNumber(field_id, "Thickness", thickness)
    field.setNumber(field_id, "Quads", 1 if quads else 0)
    field.setNumber(field_id, "IntersectMetrics", 1)

    if fan_points:
        field.setNumbers(field_id, "FanPointsList", fan_points)


def set_boundary_layer_field(field_id):
    """
    Activates the given field as THE boundary layer field.
    (Separate from the background mesh-size field — both can be active
    at once.)
    """
    gmsh.model.mesh.field.setAsBoundaryLayer(field_id)


# ============================================================
# TRAILING-EDGE CELL COUNT
# ============================================================

def get_straight_line_curves(curves):
    """
    Filters `curves` down to the ones that are straight OCC lines (as
    opposed to BSplines). The blunt-TE closing edge is built with
    occ.addLine, while the airfoil surfaces are BSplines, so this reliably
    picks out the TE edges from a curve list even after tags have been
    renumbered by a boolean cut.
    """
    return [c for c in curves if gmsh.model.getType(1, c) == "Line"]


def enforce_cells_on_curves(curves, num_cells):
    """
    Forces each curve in `curves` to mesh with exactly `num_cells` cells,
    evenly spaced, overriding whatever the background field would give it.
    """
    for curve_tag in curves:
        gmsh.model.mesh.setTransfiniteCurve(curve_tag, num_cells + 1)


def get_curve_endpoints(curves):
    """
    Returns the unique point tags at the endpoints of `curves`.

    Used to find the blunt-TE corner points (where the flat TE base line
    meets the upper/lower airfoil surfaces) so they can be passed as
    `fan_points` to the boundary layer field. Without fanning, the layer's
    default miter join at a convex corner like this stretches a single
    quad column across the corner instead of filling it, producing very
    high aspect ratio cells at the outer edge of the layer.
    """
    points = set()
    for curve_tag in curves:
        bnd = gmsh.model.getBoundary([(1, curve_tag)], oriented=False)
        for dim, tag in bnd:
            if dim == 0:
                points.add(tag)
    return list(points)


# ============================================================
# APPLY FIELD AS BACKGROUND MESH
# ============================================================

def set_background_field(field_id):
    """
    Activates the given field as the global mesh sizing function.
    """
    gmsh.model.mesh.field.setAsBackgroundMesh(field_id)


# ============================================================
# GENERATE MESH (utility wrapper)
# ============================================================

def generate_mesh(dim=2):
    """
    Generates mesh up to the given dimension (2 for 2-D CFD).
    """
    gmsh.model.mesh.generate(dim)
