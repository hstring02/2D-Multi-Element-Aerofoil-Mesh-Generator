import math
from pathlib import Path

def read_transform_airfoil(filename, chord, position, aoa_deg, te_thickness=0.0):

    points = []
    with open(filename, "r") as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            fields = line.split()

            if len(fields) != 2:
                continue

            try:
                x = float(fields[0])
                y = float(fields[1])
            except ValueError:
                continue

            points.append((x, y))

    if len(points) < 20:
        raise RuntimeError("Not enough airfoil points.")

    dx, dy = position
    rad = math.radians(aoa_deg)
    cos_a, sin_a = math.cos(rad), math.sin(rad)

    # Scale first — trailing-edge blunting needs true chordwise thickness,
    # which is only meaningful before the airfoil is rotated to its AoA.
    scaled_points = [(x * chord, y * chord) for x, y in points]

    if te_thickness > 0:
        scaled_points = blunt_trailing_edge(scaled_points, te_thickness)

    transformed_points = []

    for x_scaled, y_scaled in scaled_points:
        # Rotate around origin (0,0)
        x_rot = x_scaled * cos_a - y_scaled * sin_a
        y_rot = x_scaled * sin_a + y_scaled * cos_a

        # Translate
        x_trans = x_rot + dx
        y_trans = y_rot + dy

        transformed_points.append((x_trans, y_trans))

    return transformed_points


def split_airfoil(points):
    # Split Selig-format airfoil into upper and lower surfaces

    le_index = min(range(len(points)), key=lambda i: points[i][0])

    upper = points[:le_index + 1]
    lower = points[le_index:]

    return upper, lower


def _interp_y(curve_ascending_x, x):
    # Linear interpolation of y(x) over a curve sorted ascending by x.
    # Clamps to the nearest endpoint if x falls outside the curve's range.

    if x <= curve_ascending_x[0][0]:
        return curve_ascending_x[0][1]
    if x >= curve_ascending_x[-1][0]:
        return curve_ascending_x[-1][1]

    for (x0, y0), (x1, y1) in zip(curve_ascending_x, curve_ascending_x[1:]):
        if x0 <= x <= x1:
            if x1 == x0:
                return y0
            t = (x - x0) / (x1 - x0)
            return y0 + t * (y1 - y0)


def _segment_crossing(a, b, origin, plane_normal):
    """
    Finds where segment a->b crosses the line through `origin` whose
    normal is `plane_normal` (i.e. the line runs perpendicular to
    `plane_normal`). Returns the intersection point, or None if the
    segment doesn't cross it.
    """
    fa = (a[0] - origin[0]) * plane_normal[0] + (a[1] - origin[1]) * plane_normal[1]
    fb = (b[0] - origin[0]) * plane_normal[0] + (b[1] - origin[1]) * plane_normal[1]

    if (fa > 0) == (fb > 0):
        return None

    t = fa / (fa - fb)
    return (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]))


def _first_crossing(surface_te_to_le, origin, plane_normal):
    # Scans from the TE end inward and returns the crossing nearest the tip.
    for a, b in zip(surface_te_to_le, surface_te_to_le[1:]):
        hit = _segment_crossing(a, b, origin, plane_normal)
        if hit is not None:
            return hit
    return None


def blunt_trailing_edge(points, thickness):
    """
    Truncates the upper/lower surfaces near the trailing edge and closes
    them with a flat base of length `thickness`, oriented perpendicular
    to the local mean-camber-line direction at the tip — not the
    chordline. Upper and lower surfaces generally aren't mirror images of
    each other, so no single line is perpendicular to both; the camber
    line's local direction is the standard practical compromise.

    `points` must already be scaled (chord-length units); `thickness` is
    in the same units. Returns the truncated point list, still in Selig
    order (TE -> upper -> LE -> lower -> TE), now with an open TE gap
    equal to `thickness`.
    """

    upper, lower = split_airfoil(points)   # upper: TE->LE, lower: LE->TE

    upper_asc = list(reversed(upper))      # LE->TE, ascending x
    lower_asc = lower                      # already LE->TE, ascending x

    te_tip = (
        (points[0][0] + points[-1][0]) / 2.0,
        (points[0][1] + points[-1][1]) / 2.0,
    )

    # Local camber-line tangent at the tip, estimated from a small window
    # (1% of chord) just inboard of the TE. This direction doubles as the
    # cutting line's normal, since the cut runs perpendicular to it.
    max_x = max(upper_asc[-1][0], lower_asc[-1][0])
    min_x = min(upper_asc[0][0], lower_asc[0][0])
    x_near = max_x - 0.01 * (max_x - min_x)

    camber_tip = (_interp_y(upper_asc, max_x) + _interp_y(lower_asc, max_x)) / 2.0
    camber_near = (_interp_y(upper_asc, x_near) + _interp_y(lower_asc, x_near)) / 2.0

    dx, dy = x_near - max_x, camber_near - camber_tip
    mag = math.hypot(dx, dy)
    camber_dir = (dx / mag, dy / mag)   # unit vector, points TE -> LE

    lower_te_to_le = list(reversed(lower))

    def gap_at(s):
        origin = (te_tip[0] + s * camber_dir[0], te_tip[1] + s * camber_dir[1])
        p_upper = _first_crossing(upper, origin, camber_dir)
        p_lower = _first_crossing(lower_te_to_le, origin, camber_dir)
        if p_upper is None or p_lower is None:
            return None, None, None
        gap = math.hypot(p_upper[0] - p_lower[0], p_upper[1] - p_lower[1])
        return gap, p_upper, p_lower

    # Bracket the target thickness by growing s (distance inboard from the
    # tip, along the camber direction) until the gap catches up.
    chord_span = max_x - min_x
    s_lo, s_hi = 0.0, 0.01 * chord_span
    while True:
        gap_hi, _, _ = gap_at(s_hi)
        if gap_hi is not None and gap_hi >= thickness:
            break
        s_hi *= 2.0
        if s_hi > chord_span:
            raise RuntimeError(
                f"TE_THICKNESS {thickness} exceeds the airfoil's maximum "
                "thickness — cannot blunt the trailing edge."
            )

    for _ in range(60):
        s_mid = (s_lo + s_hi) / 2.0
        gap_mid, _, _ = gap_at(s_mid)
        if gap_mid is None or gap_mid < thickness:
            s_lo = s_mid
        else:
            s_hi = s_mid

    _, p_upper_cut, p_lower_cut = gap_at(s_hi)

    # Keep only the portion of each surface inboard of the cut (further
    # from the tip, along the camber direction, than the cut line).
    def inboard(p):
        proj = (p[0] - te_tip[0]) * camber_dir[0] + (p[1] - te_tip[1]) * camber_dir[1]
        return proj > s_hi

    new_upper = [p_upper_cut] + [p for p in upper if inboard(p)]
    new_lower = [p for p in lower if inboard(p)] + [p_lower_cut]

    return new_upper + new_lower[1:]


def add_points(occ, points, mesh_size):

    tags = []

    for x, y in points:

        tag = occ.addPoint(x, y, 0.0, mesh_size)
        tags.append(tag)

    return tags


def build_airfoil(occ, points, mesh_size):
    """
    Builds a closed airfoil surface using upper + lower splines.
    Returns:
        airfoil_surface, airfoil_curves
    """

    # Split geometry
    upper, lower = split_airfoil(points)

    # Create CAD points
    all_tags = add_points(occ, points, mesh_size)

    le_index = min(range(len(points)), key=lambda i: points[i][0])

    upper_tags = all_tags[:le_index + 1]
    lower_tags = all_tags[le_index:]

    # Create splines
    upper_curve = occ.addSpline(upper_tags)
    lower_curve = occ.addSpline(lower_tags)

    curves = [upper_curve, lower_curve]

    # Check trailing edge closure
    te_gap = math.sqrt(
        (points[0][0] - points[-1][0]) ** 2 +
        (points[0][1] - points[-1][1]) ** 2
    )

    if te_gap > 1e-8:
        te_line = occ.addLine(lower_tags[-1], upper_tags[0])
        curves.append(te_line)

    # Create surface
    loop = occ.addCurveLoop(curves)
    surface = occ.addPlaneSurface([loop])

    return surface, curves


def airfoils(data, occ, foils_dir="."):
    # GEOMETRY
    airfoil_surfaces = []
    element_info = []
    CELLS_PER_CHORD = data["global_refinement"]["CELLS_PER_CHORD"]
    for index, elements in enumerate(data["foils"]["ELEMENT"]):

        AIRFOIL_FILE = Path(foils_dir) / data["foils"]["AIRFOIL_FILE"][index]
        CHORD = data["foils"]["CHORD"][index]
        POSITION = data["foils"]["POSITION"][index]
        AOA = data["foils"]["AOA"][index]
        # Surface point spacing scales with each element's own chord, so
        # a small flap gets proportionally finer points than a large main
        # element rather than sharing one absolute size across both.
        POINT_SIZE = CHORD / CELLS_PER_CHORD
        TE_THICKNESS = data["foils"].get("TE_THICKNESS", [0.0] * len(data["foils"]["ELEMENT"]))[index]

        foil = read_transform_airfoil(AIRFOIL_FILE, CHORD, POSITION, AOA, TE_THICKNESS)

        airfoil_surface, airfoil_curves = build_airfoil(
            occ,
            foil,
            POINT_SIZE
        )
        airfoil_surfaces.append(airfoil_surface)

        # TE point in world coordinates (post-rotation, post-translation),
        # taken directly from the transformed point list rather than
        # re-derived from CHORD/POSITION/AOA, so it stays correct however
        # the airfoil was scaled/blunted/rotated above.
        te_point = (
            (foil[0][0] + foil[-1][0]) / 2.0,
            (foil[0][1] + foil[-1][1]) / 2.0,
        )

        # Y-extent of the transformed airfoil (its footprint projected onto
        # the y-axis). At high AOA this span approaches the chord length
        # rather than the airfoil's physical thickness, since the whole
        # foil is pitched over — used to size the wake box height so it
        # scales with how "tall" the airfoil actually sits at its AOA.
        foil_ys = [y for _, y in foil]
        y_min = min(foil_ys)
        y_max = max(foil_ys)

        element_info.append({
            "element": elements,
            "chord": CHORD,
            "aoa": AOA,
            "point_size": POINT_SIZE,
            "te_point": te_point,
            "y_min": y_min,
            "y_max": y_max,
            "points": foil,
        })

    return airfoil_surfaces, element_info