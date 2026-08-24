import gmsh



# ============================================================
# FARFIELD GEOMETRY
# ============================================================

def build_farfield(occ, xmin, xmax, ymin, ymax, mesh_size):
    """
    Builds rectangular farfield.
    Returns:
        surface, lines
    """

    p1 = occ.addPoint(xmin, ymin, 0, mesh_size)
    p2 = occ.addPoint(xmax, ymin, 0, mesh_size)
    p3 = occ.addPoint(xmax, ymax, 0, mesh_size)
    p4 = occ.addPoint(xmin, ymax, 0, mesh_size)

    l1 = occ.addLine(p1, p2)
    l2 = occ.addLine(p2, p3)
    l3 = occ.addLine(p3, p4)
    l4 = occ.addLine(p4, p1)

    loop = occ.addCurveLoop([l1, l2, l3, l4])
    surface = occ.addPlaneSurface([loop])

    return surface, [l1, l2, l3, l4]


# ============================================================
# BOOLEAN OPERATIONS
# ============================================================

def subtract_airfoil(occ, farfield_surface, airfoil_surface):
    """
    Returns fluid domain = farfield - airfoil
    """

    fluid = occ.cut(
        [(2, farfield_surface)],
        [(2, airfoil_surface)],
        removeObject=True,
        removeTool=True
    )

    # fluid[0] contains resulting surfaces
    return fluid[0][0]