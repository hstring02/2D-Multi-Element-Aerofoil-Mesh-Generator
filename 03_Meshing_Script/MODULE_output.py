import gmsh
import re
from pathlib import Path

# Gmsh's UNV writer leaves every node's coordinate-system reference (the
# "1 1" fields in each 2411 node record) pointing at coordinate system 1,
# but never defines it. Most UNV readers default to global Cartesian when
# a referenced coordinate system is undefined; XCALibre doesn't, and fails
# outright. This is the exact boilerplate SALOME inserts when a mesh is
# round-tripped through it (units dataset 164 + coordinate system dataset
# 2420, defining system 1 as identity/global Cartesian) — prepending it
# here removes the need for that manual import/export step.
UNV_UNITS_AND_COORD_SYS = """    -1
   164
         1  SI: Meter (newton)         2
    1.0000000000000000E+0    1.0000000000000000E+0    1.0000000000000000E+0
    2.7314999999999998E+2
    -1
    -1
  2420
         1
SMESH_Mesh
         1         0         0
Global Cartesian Coordinate System
    1.0000000000000000E+0    0.0000000000000000E+0    0.0000000000000000E+0
    0.0000000000000000E+0    1.0000000000000000E+0    0.0000000000000000E+0
    0.0000000000000000E+0    0.0000000000000000E+0    1.0000000000000000E+0
    0.0000000000000000E+0    0.0000000000000000E+0    0.0000000000000000E+0
    -1
"""


def _fmt_ints(values, width=10):
    return "".join(f"{v:{width}d}" for v in values)


def _renumber_unv(content):
    """
    Renumbers UNV dataset 2411 (nodes) and 2412 (elements) to contiguous
    1..N / 1..M labels in file order, rewriting every downstream reference
    (element vertex lists in 2412, entity tags in 2467) to match.

    Required because XCALibre's UNV2 reader ignores the label field on
    nodes entirely and treats every node/element reference elsewhere in
    the file as a direct 1-based positional index into the node/element
    arrays it built in file order. Gmsh writes each node/element's own
    internal (arbitrary, non-contiguous) tag as its label, so without this
    renumbering, references silently resolve to the wrong point/element —
    the mesh "loads" without error but its connectivity is scrambled.
    """
    lines = content.split("\n")
    n = len(lines)
    out = []
    i = 0

    node_map = {}
    elem_map = {}

    def read_ints(line):
        return [int(x) for x in line.split()]

    while i < n:
        line = lines[i]

        if line.strip() == "-1" and i + 1 < n and lines[i + 1].strip() in ("2411", "2412", "2467"):
            tag = lines[i + 1].strip()
            out.append(line)
            out.append(lines[i + 1])
            i += 2

            if tag == "2411":
                next_label = 1
                while lines[i].strip() != "-1":
                    header = read_ints(lines[i])
                    old_label = header[0]
                    node_map[old_label] = next_label
                    out.append(_fmt_ints([next_label] + header[1:]))
                    out.append(lines[i + 1])  # coordinates, unchanged
                    next_label += 1
                    i += 2
                continue

            if tag == "2412":
                next_label = 1
                while lines[i].strip() != "-1":
                    header = read_ints(lines[i])
                    old_label = header[0]
                    num_nodes = header[5]
                    elem_map[old_label] = next_label
                    out.append(_fmt_ints([next_label] + header[1:]))
                    i += 1

                    if num_nodes == 2:
                        # beam/line elements carry an extra orientation
                        # record ("0 0 0") between header and vertices
                        out.append(lines[i])
                        i += 1

                    vertex_tokens = read_ints(lines[i])
                    out.append(_fmt_ints([node_map[v] for v in vertex_tokens]))
                    i += 1

                    next_label += 1
                continue

            if tag == "2467":
                while lines[i].strip() != "-1":
                    header = read_ints(lines[i])
                    num_entities = header[-1]
                    out.append(_fmt_ints(header))  # group id/name, unchanged
                    i += 1
                    out.append(lines[i])  # group name string, unchanged
                    i += 1

                    remaining = num_entities
                    while remaining > 0:
                        tokens = read_ints(lines[i])
                        new_tokens = []
                        idx = 0
                        while idx < len(tokens) and remaining > 0:
                            etype, etag, eleaf, ecomp = tokens[idx:idx + 4]
                            new_tokens.extend([etype, elem_map[etag], eleaf, ecomp])
                            idx += 4
                            remaining -= 1
                        out.append(_fmt_ints(new_tokens))
                        i += 1
                continue

        out.append(line)
        i += 1

    return "\n".join(out)


def output_unv_xcalibre(output_filename):
    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.option.setNumber("Mesh.SaveElementTagType", 2)
    path = Path(__file__).resolve().parent.parent / "04_Meshes" / output_filename
    gmsh.write(str(path))


    # POST-PROCESS UNV FILE
    with open(path, "r") as f:
        unv_content = f.read()

    # Fix 1: Fortran exponent notation
    fixed_content = re.sub(r"(\d)[Dd]([+-]\d+)", r"\1E\2", unv_content)

    # Fix 2: Rename dataset 2477 → 2467 so XCALibre's reader finds the boundary groups.
    fixed_content = re.sub(r"(?m)^(  )2477$", r"\g<1>2467", fixed_content)

    # Fix 3: Renumber nodes/elements to contiguous labels (see _renumber_unv docstring).
    fixed_content = _renumber_unv(fixed_content)

    # Fix 4: Prepend the units + coordinate system datasets (see module docstring above).
    fixed_content = UNV_UNITS_AND_COORD_SYS + fixed_content

    with open(path, "w") as f:
        f.write(fixed_content)
