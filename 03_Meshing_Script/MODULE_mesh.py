import math

import gmsh

from MODULE_log import success, warn


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
# ELEMENT-TO-ELEMENT CLEARANCE
# ============================================================

def _point_segment_distance(p, a, b):
    ax, ay = a
    bx, by = b
    px, py = p
    dx, dy = bx - ax, by - ay
    length_sq = dx * dx + dy * dy
    if length_sq == 0.0:
        return math.hypot(px - ax, py - ay)
    t = max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / length_sq))
    cx, cy = ax + t * dx, ay + t * dy
    return math.hypot(px - cx, py - cy)


def _min_distance_to_boundary(points, boundary):
    n = len(boundary)
    return min(
        _point_segment_distance(p, boundary[i], boundary[(i + 1) % n])
        for p in points
        for i in range(n)
    )


def min_element_clearance(element_points_list):
    """
    Smallest surface-to-surface gap between any two elements' (already
    scaled/rotated/translated) airfoil point clouds.

    Checked in both directions per pair — each element's own vertices
    against the other's edges — since the true closest approach can fall
    on either element's edges depending on which one samples more densely
    near the gap; checking only one direction can overstate the real
    clearance.

    Returns math.inf if there are fewer than two elements.
    """
    n_elements = len(element_points_list)
    if n_elements < 2:
        return math.inf

    overall_min = math.inf
    for i in range(n_elements):
        for j in range(i + 1, n_elements):
            a, b = element_points_list[i], element_points_list[j]
            overall_min = min(
                overall_min,
                _min_distance_to_boundary(a, b),
                _min_distance_to_boundary(b, a),
            )
    return overall_min


def clamp_boundary_layer_thickness(element_info, thickness, first_layer_height, safety_margin=0.95):
    """
    Caps `thickness` so no two elements' boundary layers can grow into
    each other. A single global THICKNESS grows outward from every
    element's own surface with no awareness of how close another element
    sits — on tight multi-element geometry (e.g. a flap tucked under a
    main element), two elements' layers can cross, which gmsh reports as
    a self-intersecting 1D mesh ("Edge not recovered") rather than a
    clean taper. Capping thickness so each element's layer can't reach
    past the midpoint of the tightest element-to-element gap rules that
    out regardless of how close together the elements are.

    Prints a warning if it had to clamp. Raises ValueError if the
    elements are too close together to fit even one boundary-layer cell.

    Returns the (possibly unchanged) thickness to use.
    """
    if len(element_info) < 2:
        return thickness

    min_clearance = min_element_clearance([info["points"] for info in element_info])
    max_safe_thickness = safety_margin * min_clearance / 2.0

    if thickness <= max_safe_thickness:
        return thickness

    if max_safe_thickness < first_layer_height:
        raise ValueError(
            f"Elements are only {min_clearance:.6g} m apart at their closest "
            f"approach — too tight to fit even one boundary_layer cell "
            f"(FIRST_LAYER_HEIGHT={first_layer_height:.6g} m). Move the "
            "elements further apart, reduce Y_PLUS, or disable boundary_layer "
            "for this case."
        )

    warn(
        f"[boundary_layer] THICKNESS {thickness:.6g} m would overlap a "
        f"neighbouring element (closest approach {min_clearance:.6g} m) "
        f"-> clamped to {max_safe_thickness:.6g} m"
    )
    return max_safe_thickness


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
# NEAR-SURFACE FIELD (ONE PER ELEMENT)
# ============================================================

def add_near_surface_fields(element_info, element_curves, mesh_max, growth_rate, dist_min, next_field_id):
    """
    Builds one background sizing field per element, growing from that
    element's own POINT_SIZE out to mesh_max at growth_rate — rather than
    a single shared MESH_MIN floor, so a finely-pointed flap and a
    coarsely-pointed main element each get a DIST_MAX consistent with
    their own starting resolution.

    Returns the list of Threshold field IDs (one per element), ready to
    be combined into a single background field via add_min_field.
    """
    field_ids = []
    for info, curves in zip(element_info, element_curves):
        point_size = info["point_size"]
        dist_max = compute_growth_dist_max(point_size, mesh_max, growth_rate, dist_min)
        success(
            f"[{info['element']}] POINT_SIZE={point_size:.6g} GROWTH_RATE={growth_rate:.6g} -> "
            f"computed DIST_MAX={dist_max:.6g} m"
        )

        dist_field_id = next_field_id()
        add_distance_field(field_id=dist_field_id, curves=curves)

        size_field_id = next_field_id()
        add_threshold_field(
            field_id=size_field_id,
            in_field=dist_field_id,
            size_min=point_size,
            size_max=mesh_max,
            dist_min=dist_min,
            dist_max=dist_max
        )
        field_ids.append(size_field_id)
    return field_ids


# ============================================================
# WAKE REFINEMENT (BOX FIELD, ONE PER ELEMENT)
# ============================================================

def add_wake_refinement_fields(
    element_info, wake_mesh_size, mesh_max, growth_rate,
    x_start_chords, y_half_height_chords, next_field_id
):
    """
    Builds one wake-refinement Box field per element, anchored on that
    element's own TE point (world coordinates, so it already accounts for
    AOA/position) and sized as a multiple of its own chord. The box stays
    aligned with the global X axis, not the local chord line, since the
    wake convects downstream with the freestream rather than along
    whichever way an individual flap element happens to be pitched.

    Box height spans the element's own Y-extent (its footprint projected
    onto the y-axis) plus a chord-scaled pad on each side, so an element
    pitched to a high AOA — whose shed wake spans a taller Y range —
    automatically gets a taller box.

    Box downstream length ramps from 1 chord at AOA=0 up to 3 chords at
    AOA=90 (a more sharply pitched element sheds a larger, slower-
    resolving wake that needs a longer fine region to capture); AOA
    beyond 90 clamps at the 3-chord ceiling rather than continuing to grow.

    Returns the list of Box field IDs, one per element.
    """
    AOA_0_LENGTH_CHORDS = 1.0
    AOA_90_LENGTH_CHORDS = 3.0

    wake_transition = compute_growth_dist_max(wake_mesh_size, mesh_max, growth_rate)
    success(
        f"[wake_refinement] MESH_SIZE={wake_mesh_size:.6g} GROWTH_RATE={growth_rate:.6g} -> "
        f"computed transition={wake_transition:.6g} m"
    )

    field_ids = []
    for info in element_info:
        chord = info["chord"]
        te_x, _ = info["te_point"]

        aoa_frac = min(abs(info["aoa"]), 90.0) / 90.0
        x_length_chords = AOA_0_LENGTH_CHORDS + aoa_frac * (AOA_90_LENGTH_CHORDS - AOA_0_LENGTH_CHORDS)

        x_min = te_x + x_start_chords * chord
        x_max = x_min + x_length_chords * chord
        y_pad = y_half_height_chords * chord

        field_id = next_field_id()
        add_box_field(
            field_id=field_id,
            size_in=wake_mesh_size,
            size_out=mesh_max,
            x_min=x_min,
            x_max=x_max,
            y_min=info["y_min"] - y_pad,
            y_max=info["y_max"] + y_pad,
            transition=wake_transition
        )
        field_ids.append(field_id)
    return field_ids


# ============================================================
# TRAILING-EDGE CORNER REFINEMENT
# ============================================================

def add_te_corner_refinement_field(te_fan_points, mesh_size, dist_max, mesh_max, next_field_id):
    """
    Builds a Distance+Threshold field pair that tapers cell size down to
    mesh_size right at each blunt-TE corner (te_fan_points), relaxing
    back out to the surrounding background size over dist_max — flow
    accelerates sharply around these corners, so they need finer
    resolution than the rest of the surface (each element's POINT_SIZE).

    size_max is mesh_max, not any element's own POINT_SIZE, deliberately:
    this field is combined with the others via Min, so its "far" value
    only needs to stop competing with them past dist_max — it must never
    floor the whole domain, including the farfield boundary, at a size
    only meant to apply right at the airfoil surface.

    Returns the Threshold field ID, or None if there are no fan points
    (a sharp TE has no separate corner to single out).
    """
    if not te_fan_points:
        return None

    dist_field_id = next_field_id()
    add_distance_field(field_id=dist_field_id, points=te_fan_points)

    size_field_id = next_field_id()
    add_threshold_field(
        field_id=size_field_id,
        in_field=dist_field_id,
        size_min=mesh_size,
        size_max=mesh_max,
        dist_min=0.0,
        dist_max=dist_max
    )
    return size_field_id


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
